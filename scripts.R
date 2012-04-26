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
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",  
                                "integer",
                                "integer", "character", "character", "character"))
  
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  
  cat("found", nrow(gff), "rows with classes:", paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}


# Helper to load a PacBio modifications.gff file.  Parse some custom attribute fields relevant
# for modification detection.
readModsGff <- function(gffFile, nrows=-1)
{
  gffData <- gffRead(gffFile, nrows)
  gffData$coverage <- as.numeric(getAttributeField(gffData$attributes, 'coverage'))
  gffData$context <- as.character(getAttributeField(gffData$attributes, 'context'))
  gffData$IPDRatio <- as.numeric(getAttributeField(gffData$attributes, 'IPDRatio'))
  gffData
}





# example set of motifs to label

# Sequence of motifs
motifs <- c('CGACCAG', 'GANTC', 'CGACNNNNNNNTRGG', 'CCYANNNNNNNGTCG')

# Modified position in motif
motifSites <- c(6,2,3,4)


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

labels <- labelContexts(hits$context, motifs, motifSites)