# R code from Venkat Addala, https://www.researchgate.net/post/Analyzing_Human_Genome_using_R-Scripts_Here_is_a_bunch_of_R_commands_you_can_use_if_you_are_one_of_those_who_does_daily_genomic_analyses_using_R_These_commands_are_basic_you_can_push_forward_the_anal

###### 1. explore sequence composition of human genome
library(BSgenome)
available.genomes()
library(BSgenome.Hsapiens.UCSC.hg18)
 
### get sequence for chromosome 1
Seq=Hsapiens[["chr1"]]
Seq # shows some summaries
Seq=unmasked(Seq)  ## remove the mask
 
## get GC content
bases=alphabetFrequency(Seq,baseOnly=TRUE)
bases[1:4]
ntotBases=sum(bases[1:4])
baseFreq=bases[1:4]/ntotBases
GCcontent=baseFreq["C"]+baseFreq["G"]
 
## look at CG dinucleotide content
cg=matchPattern("CG", Seq)
ncg=length(cg)
## compute the observed to expected ratio
ncg/(baseFreq["C"]*baseFreq["G"]*ntotBases) ## this shows CG rarely stay together.
 
## compare to TG
tg=matchPattern("TG", Seq)
ntg=length(tg)
ntg/(baseFreq["T"]*baseFreq["G"]*ntotBases) ## this shows TG presented more than expected.
 
######### look at GC content and CG dinucleotide distribution in 1000 bp windows in whole genome.
## Note I'm only doing it in chr1 now.
ss=seq(1, length(Seq), by=1000)
ss=ss[-length(ss)] ## remove the last one
Seq.set=DNAStringSet(Seq, start=ss, end=ss+999)
ff=alphabetFrequency(Seq.set, baseOnly=TRUE)
pCG=(ff[,"C"]+ff[,"G"])/rowSums(ff)
hist(pCG[pCG>0],100)
 
## CG occurance
Seq.set=DNAStringSet(Seq, start=ss, end=ss+999)
nCG=vcountPattern("CG", Seq.set)
obsExp=nCG*1000/(ff[,"C"]*ff[,"G"])
mean(obsExp,na.rm=TRUE)
hist(obsExp,100)
 
###### 2. Explore GC content at gene TSS (transcriptional starting site),
## and their overlaps with CpG island.
## first obtain genes from UCSC using GenomicFeatures functions
library(GenomicFeatures)
txdb=makeTranscriptDbFromUCSC(genom="hg18",tablename="refGene") ## this is slow, take a few minutes.
genes=as.data.frame(transcripts(txdb))
colnames(genes)[1:3]=c("chrom","txStart","txEnd")
## Alternatively, you can read in genes from the file downloaded at last lab
genes=read.table("hg18genes.txt", comment="", header=TRUE, stringsAsFactors=FALSE)
 
## get transcriptional start site.
## Note that txEnd is the start site for genes on - strand
tss=genes$txStart
idx=genes$strand=="-"
tss[idx]=genes$txEnd[idx]
## create GRanges object for the TSS, remove random and hap chromosomes
idx=c(grep("random", genes$chrom),grep("hap", genes$chrom))
## create ranges with TSS +/- 500 bp
TSS=GRanges(seqnames=Rle(genes$chrom[-idx]), ranges=IRanges(tss[-idx]-500,tss[-idx]+500))
 
## get GC contents of TSS for genes on chromosome 1
idx.chr1=seqnames(TSS)=="chr1"
Seq.set=DNAStringSet(Seq, start=start(TSS[idx.chr1]), end=end(TSS[idx.chr1]))
ff=alphabetFrequency(Seq.set, baseOnly=TRUE)
pCG.TSS=(ff[,"C"]+ff[,"G"])/rowSums(ff)
hist(pCG.TSS,100)
 
## compare with genome wide distribution of GC content
d1=density(pCG[pCG>0])
d2=density(pCG.TSS)
plot(d1, lwd=2)
lines(d2, col="red",lwd=2)
legend("topright", legend=c("Genome", "TSS"), lwd=2, col=c("black", "red"))
 
## CG dinucleotide around TSS, compute obsExp
nCG=vcountPattern("CG", Seq.set)
obsExp.TSS=nCG*1000/(ff[,"C"]*ff[,"G"])
hist(obsExp.TSS,100)
 
## compare with the genome wide one:
d1=density(obsExp[!is.na(obsExp)])
d2=density(obsExp.TSS)
plot(d1, lwd=2, xlim=c(0,2), main="Observed to expected CG ratio")
lines(d2, col="red",lwd=2)
legend("topright", legend=c("Genome", "TSS"), lwd=2, col=c("black", "red"))
 
 
##########
### now read in CpG island file and compare it to TSS
cgi=read.table("hg18CGI.txt", comment="", header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(cgi)[1]="chrom" ## fix column name
idx=c(grep("random", cgi$chrom),grep("hap", cgi$chrom))
CGI=GRanges(seqnames=Rle(cgi$chrom[-idx]), ranges=IRanges(cgi$chromStart[-idx],cgi$chromEnd[-idx]))
 
## look at % of TSS overlapping CGI
overlaps=countOverlaps(TSS, CGI)
mean(overlaps>0)
## or do
mean(TSS %in% CGI)
 
#####################################################################
## Time allows, perform following to explore worm gene annotations.
## otherwise students can keep and study the scripts.
#####################################################################
library(GenomicFeatures)
txdb=makeTranscriptDbFromUCSC(genom="ce2",tablename="refGene")
genes=transcriptsBy(txdb)
length(genes)
## number of +/- genes
strands=as.character(unlist(strand(genes)))
table(strands)
## gene lengths
geneLength=unlist(end(genes))-unlist(start(genes))
## distribution of gene lengths versus strand
boxplot(log(geneLength)~strands)
t.test(log(geneLength)~strands)
 
## number of genes on each chromosome
chrs=as.character(seqnames(genes))
ngenes.chr=table(chrs)
ngenes.chr
 
## number of genes on each chromosome versus chromosome lengths
## pay attention to the way to plot
chrlen=seqlengths(genes)
## there's no genes on chrM, remove it in plotting
chrlen2=chrlen[names(ngenes.chr)]/1e6
plot(chrlen2, as.numeric(ngenes.chr), type="n", xlab="Chromosome length (Mbp)", ylab="Number of genes", main="C. elegans")
text(chrlen2, ngenes.chr, names(ngenes.chr))
