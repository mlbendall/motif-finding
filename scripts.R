library(Biostrings)
library(grid)

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

filterGff <- function(gff,idx,labels,fOut)
{
  gffDat <- gffRead(gff)
  gffDat <- gffDat[idx,]
  gffDat$motif <- labels
  write.table(file=fOut,gffDat,row.names=F,col.names=T,sep=' ')
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
labelContexts <- function(contextStrings, motifs, motifSites, labels=motifs, contextStringCenter=21)
{
  labelVect <- character(length(contextStrings))
  labelVect[1:length(contextStrings)] <- 'None'

  for (m in 1:length(motifs))
  {
    isHit <- as.vector(isMatchingAt(motifs[m], 
                                    DNAStringSet(as.character(contextStrings)), 
                                    at=contextStringCenter + 1 -motifSites[m],  
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
          revDf$start <- length(dna_seq[[ref]]) + 1 - (revDf$start + methPos - 1 + o)
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

# Write the contextStrings to a FASTA file
writeContextFasta <- function(dat, filename)
{
  stringSet <- DNAStringSet(dat$context)
  names(stringSet) <- paste(as.character(dat$start),dat$strand,sep='')
  write.XStringSet(stringSet, filename)
}

# Plot kinetics = c('ipdRatio','score','coverage')
plotKinetics <- function(dat,ref,statMode,plotRange,rangePerPlot)
{
  dat <- subset(dat,seqid=ref)
  tableNames <- names(dat)
  idx <- match(statMode,tableNames)
  posIdx <- match('tpl',tableNames)
  baseIdx <- match('base',tableNames)
  plotTitle <- statMode
  fontSize <- floor(800/rangePerPlot)
  baseline <- median(dat[,idx],na.rm=TRUE)
  if (statMode =='coverage') baseline <- 0
  yFwdStr <- 'Forward Strand'
  yRevStr <- 'Reverse Strand'
 
  posRange <- plotRange
  nPlots <- ceiling((posRange[2]-posRange[1])/rangePerPlot)
  
  for (n in 1:nPlots)
  {
    sWin <- (n-1)*rangePerPlot+posRange[1]
    eWin <- sWin+rangePerPlot
    xbin <- 10
    subDat <- subset(dat,tpl >= sWin & tpl <= eWin)
    fwdDat <- subDat[subDat$strand=='+',c(posIdx,idx,baseIdx)]
    revDat <- subDat[subDat$strand=='-',c(posIdx,idx,baseIdx)]
    yMin <- -0.5
    yMax <- max(subDat[,idx],na.rm=TRUE)*1.05
    if (yMax < 5) yMax = 5
    plotDat <- data.frame(tpl = sWin:eWin, statFwd = baseline, statRev = baseline)
    plotDat$baseRev = "x"
    plotDat$baseFwd = "x"
    
    posMatch <- match(plotDat$tpl,fwdDat$tpl)
    fwdIdx <- posMatch[which(is.na(posMatch)==FALSE)]
    plotDat$statFwd[!is.na(posMatch)] = fwdDat[fwdIdx,2]
    plotDat$baseFwd[!is.na(posMatch)] = as.character(fwdDat$base[fwdIdx])
    
    posMatch <- match(plotDat$tpl,revDat$tpl)
    revIdx <- posMatch[which(is.na(posMatch)==FALSE)]
    plotDat$statRev[!is.na(posMatch)] = revDat[revIdx,2]
    plotDat$baseRev[!is.na(posMatch)] = as.character(revDat$base[revIdx])
    plotDat$baseline <- baseline
    
    p1 <- ggplot(data=plotDat, aes(xmin=tpl-0.4,xmax=tpl+0.4,ymin=baseline,ymax=statFwd))+geom_rect(fill='salmon',label='')+ylab(yFwdStr)+ylim(c(yMin,yMax))+ scale_x_continuous(breaks=seq(sWin,eWin,xbin))+opts(axis.text.x=theme_blank(),legend.position = "none",title=plotTitle,plot.margin= unit(c(0, 2, -.5, 0), "lines"))+coord_cartesian(xlim = c((sWin-1),(eWin+1)))
    pp1 <- p1 + geom_text(aes(x=tpl, y=0, vjust=1, label=toupper(baseFwd)), size=fontSize)+opts(axis.text.x = theme_blank(),axis.label.x = theme_blank())
    
    p2 <- ggplot(data=plotDat, aes(xmin=tpl-0.4,xmax=tpl+0.4,ymin=baseline,ymax=statRev))+geom_rect(fill='salmon',label='')+xlab('')+ylab(yRevStr)+ylim(c(yMax,yMin))+ scale_x_continuous(breaks=seq(sWin,eWin,xbin))+opts(legend.position = "none",plot.margin= unit(c(-0.5, 2, 0, 0), "lines"))+xlab('Ref Position')+coord_cartesian(xlim = c((sWin-1),(eWin+1)))
    pp2 <- p2 + geom_text(aes(x=tpl, y=0, vjust=0, label=toupper(baseRev)), size=fontSize)
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))
    print(pp1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(pp2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  }
}

readModCsv <- function(csvFile,genomeSize,chunkSize)
{
  fileSize <- genomeSize*2
  nchunks <- ceiling(fileSize/chunkSize)
  tmp <- read.table(file=csvFile,header=TRUE,sep=',',skip=0,nrow=5)
  colClasses <- sapply(tmp, class)
  csvHeader <- names(tmp)
  keepCols <- c('refId','tpl','strand','ipdRatio','pvalue','coverage')
  theTable <- NULL
  for (n in 1:nchunks)
  {   
    sk = (n-1)*chunkSize+1
    last = sk+chunkSize
    if (last > fileSize) chunkSize <- fileSize-sk
    tmp <- read.table(file=csvFile,header=FALSE,sep=',',skip=sk,nrow=chunkSize,colClasses=colClasses)
    names(tmp) <- csvHeader
    theTable <- rbind(theTable,tmp[,keepCols])
  }
  return(theTable)  
}
