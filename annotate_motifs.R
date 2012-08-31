#! /jgi/tools/bin/Rscript

instructions <- "
Creates a table with motif locations
 Parameters:
   reference_file: FASTA reference file 
   motif_file:     motifs.csv output file from MotifMaker
   outfile:        name of output file

  Creates a table that looks like this:
  strand  start  motif  onTarget  seqid
  +       1626   GATC   On        1
  +       1774   GATC   On        1
  +       1924   GATC   On        1
  +       2094   GATC   On        1
  +       2117   GATC   On        1
"


makeAnnotationTable <- function(reference_file,motif_file,outfile) {
  mt = read.csv(motif_file)
  motifs    <- as.character(mt$motifString)
  # centerPos in motifs.csv is zero-indexed, genomeAnnotation function expects 1-indexed
  positions <- mt$centerPos + 1
  dna_seq <- read.DNAStringSet(reference_file)
  genomeAnnotations <- genomeAnnotation(dna_seq,motifs,positions)
  table(genomeAnnotations$motif)
  write.table(genomeAnnotations,file=outfile,sep="\t",row.names=FALSE,quote=FALSE)
}

args <- commandArgs(TRUE)
if(length(args) == 3) {
  source('scripts.R')

  reference_file <- args[1]
  motif_file     <- args[2]
  outfile        <- args[3]
  
  makeAnnotationTable(reference_file,motif_file,outfile)
  
  q(save="no",status=0)
}

write("USAGE: Rscript annotate_motifs.R <reference_file> <motif_file> <outfile>",stdout())
write(instructions,stdout())