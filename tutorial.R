#R can be downloaded for Linux, Mac and Windows from http://r-project.org We also recommend the graphical front-end RStudio, available from rstudio.org. The following R packages need to be installed to run through this demo. ggplot2, plyr are available through the built-in R package manager. Instructions for installing Bioconductor packages is available here: www.bioconductor.org/install. Biostrings are required from Bioconductor. The commands are:
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#scripts.R is included in this tutorial.
library(ggplot2)
source('~/motif-finding/scripts.R')

gff <-'/mnt/secondary/Smrtanalysis/userdata/jobs/040/040041/data/modifications.gff.gz' 

hits <-readModificationsGff(gff) 
hits$CognateBase <-substr(hits$context, 21, 21) 
head(hits) 

summary(hits$coverage) 


qplot(score, colour=CognateBase, geom='freqpoly', data=hits, binwidth=1,xlim=c(0,200))+scale_y_continuous(breaks=seq(0,30000,1000))+ylim(c(0,2000))

qplot(coverage, score, colour = CognateBase, alpha = I(0.2), data = subset(hits, CognateBase %in% c('A')))


#Figure 2: Score vs. Coverage 

qplot(coverage, score, colour = CognateBase, alpha = I(0.2), data = subset(hits, CognateBase == "A"), xlim=c(0,250), ylim=c(0,250))+geom_abline(slope=1,intercept=0)+geom_vline(x=25)

#Figure 3. Illustration of cutoffs
#We select these hits, then sort in decreasing order of score, so we consider the strongest signal ﬁrst.

goodHits <-subset(hits, coverage>25)
goodHits <-subset(goodHits, score > coverage)
goodHits <-goodHits[order(goodHits$score, decreasing=T),]
workHits <-goodHits

# MEMEChip Motif Finding 
#We use the online motif ﬁnding server ’MEMEChip’ (http://meme.sdsc.edu/meme/cgi-bin/meme-chip.cgi). The MEMEChip server requires the source sequences to uploaded in FASTA format. We write the context string of the top 1000 hits to ’contexts.fasta’. 

writeContextFasta(workHits[1:1000,], 'contexts.fasta') 

#We can count the number of hits matching to this motif:

motif <-c('GATC','NNNNGATC')
position <-c(2,1)
motifLabels <-labelContexts(workHits$context, motif, position) 
table(motifLabels) 
workHits$label = motifLabels

remain <-workHits[which(motifLabels == 'None'),] 
nrow(remain)
  
#We then iterate by submission to MEMEChip until all motifs have been identified. In our example, the 3 motifs are ’GATC’, ’GCACNNNNNNGTT’, and its complement ’AACNNNNNNGTGC’ as methylated sequence motifs. 
writeContextFasta(remain[1:1000,], 'remaining.fasta') 




# Genome Annotation 
#We can load the genome sequence and annotate its instances of each motif, to determine the genome-wide methylation fraction of our motifs. The genomeAnnotation function returns a data.frame containing one row for each match in the supplied genome of each motif. Genome positions can match multiple motifs, which gives multiple rows at the same genome position. 
#Load the genome sequence

seq_path <-'/mnt/secondary/Smrtanalysis/userdata/references/ecoli/sequence/ecoli.fasta' 

dna_seq <-read.DNAStringSet(seq_path)
motifs = c('GATC', 'GCACNNNNNNGTT', 'AACNNNNNNGTGC')
positions = c(2, 3, 2) #position of the modified base in the motifs string
genomeAnnotations <-genomeAnnotation(dna_seq, motifs, positions)
head(genomeAnnotations) 

table(genomeAnnotations$motif) 

#We merge the genomeAnnotation output with our GFF data.frame to count the number of motif instances that exist in the genome and the number that were detected as methylated. We merge the genomeAnnotation and workHits tables by genome position and strand, with all=T to include GFF hits that are not annotated with a motif, and genome motif instances that do not have a GFF hit. We adjust the merged data.frame to indicate these cases. 

workHits$seqid <-as.integer(substr(workHits$seqname, 4,11))
mm <-merge(workHits, genomeAnnotations, all=T) 
mm$motif[is.na(mm$motif)] <-'NoMotif'
mm$feature[is.na(mm$feature)] <-'not_detected' 
table(mm$feature, mm$motif) 

#Score distributions for identified motifs and remaining hits.

qplot(score, ..density.., colour=motif, geom='freqpoly', data=subset(mm, feature=='modified_base'), binwidth=3, xlim=c(0,200))





#Unmethylated Sites 
#Our comparison of the genome annotation and the detected methylations reveals that there are a small number of GATC sites in the genome that were not detected as methylated at our selected cutoﬀs. We may be interested in looking in more detail at the kinetic evidence of the undetected GATC sites to determine whether these sites are false negatives (through insufficient coverage), or are truly unmethylated in the organism. The GFF ﬁle only contains genome positions with a score > 20. The modifications.csv.gz ﬁle contains the complete statistics for all sites in the genome. 
#Here we load the CSV ﬁle, and merge it with the genomeAnnotations table, keeping only rows containing an identiﬁed motif. We can recreate our detection table by cutting on the selected cutoff criteria. Figure 4 shows the low end of the score distributions of the matching motif sites. It appears that there is a small population of GATC sites with very low score. We look at the table of hits with score < 10. These sites appear to have good coverage, and small ipd ratios. Interestingly, 12 of 16 of these sites are paired, showing no methylation on both strands, which gives us more conﬁdence in their unmethylated status. 

# if in desktop environment

#genomeSize <- sum(width(dna_seq))
#csv <-'modifications.csv.gz'
#rawKin <-readModCsv(csv,genomeSize,500000) 

# if in server environment

csv <-'/mnt/secondary/Smrtanalysis/userdata/jobs/040/040041/data/modifications.csv.gz'
rawKin <-read.csv(csv) 


rawKin$seqid <-as.integer(substr(rawKin$refId,4,11))
rawKin$strand[rawKin$strand == 1] <-'-'
rawKin$strand[rawKin$strand == 0] <-'+'
mga <-merge(genomeAnnotations, rawKin, by.x=c('start', 'strand'), by.y=c('tpl', 'strand'))

#Score distribution of motif sites 
qplot(score, ..density.., log='y',colour=motif, data=mga, geom='freqpoly', xlim=c(0,100), binwidth=1)

table(mga$motif, mga$score > mga$coverage) 


colsToShow <-c('start', 'strand', 'motif', 'score', 'coverage', 'ipdRatio')
subset(mga, score < 10 & motif=='GATC')[,colsToShow] 

#Plot ipdRatio, score or coverage by position

comp_dnaSeq <- complement(dna_seq)
pos <- 1:width(dna_seq)
tplBase <- data.frame(tpl=pos,strand='+',base=unlist(strsplit(as.character(dna_seq),split='')),row.names=NULL)
rTplBase <- data.frame(tpl=pos,strand='-',base=unlist(strsplit(as.character(comp_dnaSeq),split='')),row.names=NULL)
tplBase <- rbind(tplBase,rTplBase)
newRaw <- merge(tplBase,rawKin,by.x=c('tpl','strand'),by.y=c('tpl','strand'))
head(newRaw)

refId <- 1
statMode <- 'ipdRatio'
plotRange <- c(1000,2000)
rangePerPlot <- 200
pdfName <- 'test.pdf'
pdf(file=pdfName,height=6,width=25)
plotKinetics(newRaw,refId,statMode,plotRange,rangePerPlot)
dev.off()

### plot random selection of hits ###

targetRef <- 1
subHits <- workHits[which(workHits$label=='GATC'),]
sampleSize <- 5
subDat <- subHits[which(subHits$seqid==targetRef),]
if (nrow(subDat) < sampleSize) sampleSize <- nrow(subDat)
nsamples=sample(subDat$start,sampleSize,replace=F)
pdf(file='GATCExamples.pdf',height=6,width=25)
for (s in 1:length(nsamples))
{
  plotRange[1] <- nsamples[s]-round(rangePerPlot/2)
  plotRange[2] <- nsamples[s]+round(rangePerPlot/2)
  plotKinetics(newRaw,targetRef,statMode,plotRange,rangePerPlot)
}
dev.off()



### write out filtered gff

goodHitsIndx <-which(hits$coverage>25 & hits$score > hits$coverage)
filterGff(gff,goodHitsIndx,motifLabels[goodHitsIndx],'filteredGff.gff')






#A.1 Session Info 
#Here we give information about the version of R and installed packages used to generate this document. 
#> sessionInfo() 
#R version 2.14.2 Patched (2012-02-29 r59005) Platform: x86_64-unknown-linux-gnu (64-bit) 
#locale: 
#  [1] C 
#attached base packages: 
#  [1] grid stats graphics grDevices utils datasets methods #

#[8] base 
#other attached packages: 
#  [1] cosmo_1.18.0 seqLogo_1.18.0 Biostrings_2.22.0 IRanges_1.12.6 
#[5] plyr_1.7.1 ggplot2_0.9.0 
#loaded via a namespace (and not attached): 
#  [1] MASS_7.3-17 RColorBrewer_1.0-5 colorspace_1.1-1 dichromat_1.2-4 
#[5] digest_0.5.2 memoise_0.1 munsell_0.3 proto_0.3-9.2 
#[9] reshape2_1.2.1 scales_0.2.0 stringr_0.6 tools_2.14.2 
#For Research Use Only. Not for use in diagnostic procedures. Copyright 2011, Paciﬁc Biosciences of California, Inc. All rights reserved. Information in this document is subject to change without notice. Paciﬁc Biosciences assumes no responsibility for any errors or omissions in this document. Certain notices, terms, conditions and/or use restrictions may pertain to your use of Paciﬁc Biosciences products and/or third party products. Please refer to the applicable Paciﬁc Biosciences Terms and Conditions of Sale and to the applicable license terms at http://www.paciﬁcbiosciences.com/licenses.html. 
#Paciﬁc Biosciences, the Paciﬁc Biosciences logo, PacBio, SMRT and SMRTbell are trademarks of Paciﬁc Biosciences in the United States and/or certain other countries. All other trademarks are the sole property of their respective owners. 
