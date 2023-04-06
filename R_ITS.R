###############################################
####FUNGAL ITS SEQUENCE ANALYSIS WITH DADA2####
###############################################

# Load required libraries
library(pacman)
p_load(Rcpp, usethis, devtools, dada2, ShortRead, Biostrings, phyloseq, ggplot2,
       dplyr, magrittr, readr, stringr, tibble)

# PATH TO SAMPLES : 
  samplesPath <- "MiSeq_4Maladies"
  
# PATH TO TAXONOMY
  unite.ref <- paste0(samplesPath,"/sh_general_release_dynamic_10.05.2021.fasta")  # CHANGE ME to location on your machine
  
# Can you MULTITHREAD?
mthreads <- TRUE # TRUE if yes
  
list.files(samplesPath)

#######################################################################
# Primers identification ##############################################
#######################################################################

# Sort Rev and Fwd fastq files + add primers
fnFs <- sort(list.files(samplesPath, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(samplesPath, pattern = "_R2.fastq.gz", full.names = TRUE))

REV <- "GGAAGTAAAAGTCGTAACAAGG"
FWD <- "TCCTCCGCTTATTGATATGC"

# Check primer orientation
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients # it's all right

# Put N-filterd files in filtN/ subdirectory :
fnFs.filtN <- file.path(samplesPath, "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(samplesPath, "filtN", basename(fnRs))

# N-filtering
filterAndTrim(fnFs, fnFs.filtN, # input file, then output path
              fnRs, fnRs.filtN, 
              maxN = 0, multithread = mthreads)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
   REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
   REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#######################################################################
# Primers removal #####################################################
#######################################################################

system2("/Users/jorondo/miniconda3/bin/cutadapt", args = "--version")

# output dir and paths for primer-free sequences
path.cut <- paste0(getwd(), "/cut")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# cutadapt commands for forward and reverse removal
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for(i in seq_along(fnFs)) {
  c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
    "-o", fnFs.cut[i], 
    "-p", fnRs.cut[i], # output files
    fnFs.filtN[i], fnRs.filtN[i]) %>% 
      purrr::when(mthreads==TRUE ~ c(., "-j", 16)) %>% # input files
      system2(cutadapt, args = .)
}

# sanity check :
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#######################################################################
# Sample names formatting #############################################
#######################################################################

cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
list(cutFs)
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

# Extract sample names
get.sample.name <- function(fname) {
  paste0(strsplit(basename(fname), "_")[[1]][1],
         "_", strsplit(basename(fname), "_")[[1]][2]) %>% return
}
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)


#######################################################################
# Quality filtering ###################################################
#######################################################################

# Pre-filtering visualization:
plotQualityProfile(cutFs[1:2])

# Paths to output
filtFs <- file.path("filtered", basename(cutFs))
filtRs <- file.path("filtered", basename(cutRs))

# Quality filtering
out <- filterAndTrim(cutFs, filtFs, 
                     cutRs, filtRs,
                     maxN = 0, maxEE = c(2,2), 
                     truncQ = 2, minLen = 20, rm.phix = TRUE, 
                     compress = TRUE, multithread = mthreads)
head(out)
# Post-filtering visualization
plotQualityProfile(filtRs[1:2])

#######################################################################
# Error model #########################################################
#######################################################################

errF <- learnErrors(filtFs, multithread = mthreads)
errR <- learnErrors(filtRs, multithread = mthreads)
plotErrors(errF, nominalQ = TRUE)

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = mthreads)
dadaRs <- dada(derepRs, err = errR, multithread = mthreads)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                      verbose=TRUE, justConcatenate = TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#######################################################################
# Bimera removal #####################################################
#######################################################################

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread = mthreads, verbose=TRUE)

table(nchar(getSequences(seqtab.nochim)))
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", 
                     "nonchim")
rownames(track) <- sample.names
head(track) 

##assign taxa
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = mthreads, tryRC = TRUE)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
samples.out <- rownames(seqtab.nochim)
print(samples.out)

saveRDS(taxa,"taxa.RDS")
saveRDS(seqtab.nochim, "asv_table.RDS")

#######################################################################
### Phyloseq Object ###################################################
#######################################################################

asv <- readRDS("asv_table.RDS")

# Generate ASV ID
taxonomy <- 
  readRDS("taxa.RDS"); asvNames <- 
  data.frame(
    ASV = seq(1,nrow(taxonomy),1) %>% # sequence for unique ids
      str_pad(4,"left","0") %>%       # leading zeros
      paste0("ASV_",.),               # 
  ITS = rownames(taxonomy))

rownames(taxonomy) <- asvNames[,"ASV"]

# to make sure the ITS names are in the right order:
colnames(asv) %<>% 
  data.frame(ITS = .) %>% 
  full_join(asvNames, by="ITS") %$% ASV # expose the ASV IDs, assign them to colnames

# Samples data
sampData <- data.frame(sampID = rownames(asv),
             Disease = read_csv("meta.csv") %$% Disease,
             Status = read_csv("meta.csv") %$% Status) %>% 
  column_to_rownames("sampID")

ps <- phyloseq(otu_table(asv, taxa_are_rows = F),
         sample_data(sampData),
         tax_table(as.matrix(taxonomy)))

saveRDS(ps,"ps.RDS")
