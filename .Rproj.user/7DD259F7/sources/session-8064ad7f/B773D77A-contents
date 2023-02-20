###############################################
####FUNGAL ITS SEQUENCE ANALYSIS WITH DADA2####
###############################################


# Load required libraries
library(pacman)
p_load(Rcpp, usethis, devtools, dada2, ShortRead, Biostrings, phyloseq, ggplot2)

# PATH TO SAMPLES : 
  samplesPath <- "MiSeq_4Maladies"
  
# PATH TO TAXONOMY
  unite.ref <- paste0(samplesPath,"/sh_general_release_dynamic_10.05.2021.fasta")  # CHANGE ME to location on your machine
  
# Can you MULTITHREAD?
  mthreads <- TRUE # TRUE if yes
  
  
list.files(samplesPath)

# Sort Rev and Fwd fastq files + add primers
fnFs <- sort(list.files(samplesPath, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(samplesPath, pattern = "_R2.fastq.gz", full.names = TRUE))

REV <- "GGAAGTAAAAGTCGTAACAAGG"
FWD <- "TCCTCCGCTTATTGATATGC"

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
FWD.orients


fnFs.filtN <- file.path(samplesPath, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(samplesPath, "filtN", basename(fnRs))
# filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = mthreads)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
#    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
#   REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
#   REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#Cutting primers
#cutadapt <- "/Users/gaelelajeunesse/opt/miniconda3/bin/cutadapt"
#system2(cutadapt, args = "--version")

#path.cut <- getwd()

#if(!dir.exists(path.cut)) dir.create(path.cut)
#fnFs.cut <- file.path(path.cut, basename(fnFs))
#fnRs.cut <- file.path(path.cut, basename(fnRs))

#FWD.RC <- dada2:::rc(FWD)
#REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
#R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
#R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
#for(i in seq_along(fnFs)) {
#system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
#                           "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
#                          fnFs.filtN[i], fnRs.filtN[i])) # input files
#}

#rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
#    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
#    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
#    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(samplesPath, pattern = "_R1.fastq.gz", full.names = TRUE))
list(cutFs)


#cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#plotQualityProfile(cutFs[1:2])

filtFs <- file.path(samplesPath, "filtered", basename(cutFs))
#filtRs <- file.path(samplesPath, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2, 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = mthreads)  # on windows, set multithread = FALSE
head(out)

errF <- learnErrors(filtFs, multithread = mthreads)
#errR <- learnErrors(filtRs, multithread = mthreads)
plotErrors(errF, nominalQ = TRUE)

derepFs <- derepFastq(filtFs, verbose = TRUE)
#derepRs <- derepFastq(filtRs, verbose = TRUE)


names(derepFs) <- sample.names
#names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = mthreads)
#dadaRs <- dada(derepRs, err = errR, multithread = mthreads)

#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = 5)

seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread = mthreads, verbose=TRUE)
table(nchar(getSequences(seqtab.nochim)))

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
