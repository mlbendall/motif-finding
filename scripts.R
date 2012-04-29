library(Biostrings)

# Parse an attirubte field from the GFF attributes field by name
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}

# Load a GFF file into a data.frame.  The attributes column is left alone
gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",  
                                "integer",
                                "integer", "character", "character", "character"))
  
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  return(gff)
}


# Helper to load a PacBio modifications.gff file.  Parse some custom attribute fields relevant
# for modification detection.
readModificationsGff <- function(gffFile, nrows=-1)
{
  gffData <- gffRead(gffFile, nrows)
  gffData$coverage <- as.numeric(getAttributeField(gffData$attributes, 'coverage'))
  gffData$context <- as.character(getAttributeField(gffData$attributes, 'context'))
  gffData$IPDRatio <- as.numeric(getAttributeField(gffData$attributes, 'IPDRatio'))  
  gffData$attributes <- NULL
  
  gffData
}

# example set of motifs to label

# Sequence of motifs
# motifs <- c('CGACCAG', 'gantc', 'CGACNNNNNNNTRGG', 'CCYANNNNNNNGTCG')

# Modified position in motif
# motifSites <- c(6,2,3,4)

# Label each element of contextString as containing a motifs, or 'None'
# motifs is a list of motif strings, possibly including ambiguous bases
# motifSites is the index of the modified position in the motif
# labels is an optional label to apply - by default the string of the matching 
# motif will be used
labelContexts <- function(contextStrings, motifs, motifSites, labels=motifs)
{
  labelVect <- character(length(contextStrings))
  labelVect[1:length(contextStrings)] <- 'None'

  for (m in 1:length(motifs))
  {
    isHit <- as.vector(isMatchingAt(motifs[m], 
                                    DNAStringSet(as.character(contextStrings)), 
                                    at=22-motifSites[m],  
                                    fixed=F))
    labelVect[isHit] <- motifs[m]
  }
  
  factor(labelVect)
}



# Generate a data.frame with one row for each example of motif in the reference
# meth pos is the position of the m6A in the motif, offset=0 means to expect a 
# peak at the metylated position, offset=-5 is the normal position of the 
# secondary peak for m6A.  If there are Ns in the motif tb.end must be set 
# to the end of the first block of non-Ns
genomeAnnotation <- function(dna_seq, motifs, positions, offsets=c(0))
{
  df <- data.frame()
  dna_seq_rc <- reverseComplement(dna_seq)
  
  for (ref in 1:length(dna_seq))
  {
    for(m in 1:length(motifs))
    {  
      motif <- motifs[m]
      methPos <- positions[m]
      
      maskFwd <- mask(dna_seq[[ref]], 'N')      
      maskRev <- mask(dna_seq_rc[[ref]], 'N')
      
      match_fwd <- matchPattern(motif, maskFwd, fixed=F)
      match_rev <- matchPattern(motif, maskRev, fixed=F)  
      
      for (o in offsets)
      {
        modName <- motif
        if(o > 0)
          modName = paste(motif, "p", o, sep='')
        if(o < 0)
          modName = paste(motif, "n", o, sep='')
        
        onTarget <- 'On'
        if(abs(o) > 0)
          onTarget <- 'Off'
        
        fwdDf <- data.frame(start=start(match_fwd))
        if(nrow(fwdDf) > 0)
        {
          fwdDf$seqid <- ref
          #fwdDf$contig <- names(dna_seq)[[ref]]
          
          fwdDf$strand <- '+'
          fwdDf$start <- fwdDf$start + methPos - 1 + o
          fwdDf$motif <- modName
          fwdDf$onTarget <- onTarget  
          df <- rbind(df, fwdDf)
        }
        
       revDf <- data.frame(start=start(match_rev))        
       if(nrow(revDf) > 0)
        {
          revDf$seqid <- ref
          #revDf$contig <- names(dna_seq)[[ref]]
          
          revDf$strand <- '-'
          revDf$start <- length(dna_seq[[1]]) + 1 - (revDf$start + methPos - 1 + o)
          revDf$motif <- modName
          revDf$onTarget <- onTarget
          df <- rbind(df, revDf)
        }
      }
    }
  }
  
  df <- df[,c('strand', 'start', 'motif',  'onTarget', 'seqid')]
  return(df)
}

