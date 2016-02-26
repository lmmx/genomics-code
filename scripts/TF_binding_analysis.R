library('ChIPseeker') # also imports plyr, magrittr, GenomicRanges, TxDb.Hsapiens.UCSC.hg19.knownGene
library('tools')
library('biomaRt')

### Function Declarations

# iterate over guts of analysis in this function when finished drafting
AtacOverlap <- function(peak_name) {
  atac.peaks <<- get(peak_name)
  qc.atac.peaks <<- atac.peaks * -0.4 # scale to 40% width
  
  # Number of hits, including duplicated TF peaks (1-to-many TFBS-to-ATAC peak mapping)
  
  # length(which(duplicated(ranges(ref.peaks[subjectHits(findOverlaps(atac.peaks, ref.peaks))]))))
  
  qc.tfbs.at.atac <<- subsetByOverlaps(tf.peaks, qc.atac.peaks)
  n.untrimmed.overlaps <<- length(findOverlaps(tf.peaks, atac.peaks))
  n.unique.untrimmed.overlaps <<- length(subsetByOverlaps(tf.peaks, atac.peaks))
  n.qc.overlaps <<- length(findOverlaps(tf.peaks, qc.atac.peaks))
  n.unique.qc.overlaps <<- length(qc.tfbs.at.atac)
  
  paste0(peak_name %>%
           gsub(".peaks", "", .),
         " has ",
         n.unique.qc.overlaps,
         " (",
         round ( 100 * ( n.unique.qc.overlaps / n.qc.overlaps ), digits = 1 ),
         "% unique of ",
         n.qc.overlaps,
         " QC-trimmed) TFBS in ATAC-defined 'accessible regions', or ",
         n.unique.untrimmed.overlaps,
         " (",
         round ( 100 * ( n.unique.untrimmed.overlaps / n.untrimmed.overlaps ), digits = 1 ),
         "% unique of ",
         n.untrimmed.overlaps,
         ") without QC trimming.") %>%
    print()
}

mapSitesToTSS <- function() {
  tfbs.at.tss <- subsetByOverlaps(tss.us, qc.tfbs.at.atac)
  qc.tfbs.at.atac[qc.tfbs.at.atac %over% tss.us]
}

### Initial setup

# Use the reference list of peaks Yaoyong provided:
tf.peakfile <- '../../peaks/H1_FOXK2_strongest_peaks.narrowPeak'
tf.peaks <- readPeakFile(file.path(tf.peakfile))
seqlevels(tf.peaks) <- paste0('chr', seqlevels(tf.peaks))

# Make global Biomart
if (!exists('ensembl')) {
  ensembl <- useMart(mart = "ENSEMBL_MART_ENSEMBL",
                     dataset = "hsapiens_gene_ensembl",
                     host = "www.ensembl.org")
}

entrezToENSG <- function(entrez.gene.list) {
  getBM(attributes = c("entrezgene", "entrezgene_transcript_name", "ensembl_gene_id"),
        filters = 'entrezgene',
        mart = ensembl,
        values = entrez.gene.list)
}

# Set up a reference to the transcripts in the hg19 reference genome (for TSS coordinates)
if (!exists('txdb_hg19')) {
  txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genedf <- as.data.frame(transcripts(txdb_hg19))
  colnames(genedf)[1:3] <- c("chrom","txStart","txEnd")
  genedf$geneID <- select(TxDb.Hsapiens.UCSC.hg19.knownGene, as.character(genedf$tx_id), "GENEID", keytype = "TXID")[['GENEID']]
  
  id.vals <- genedf$tx_name[is.na(genedf$geneID)]
  entrezless.ensg.fetch <- getBM(attributes = c("ucsc","ensembl_gene_id"),
                                 filters = 'ucsc',
                                 mart = ensembl,
                                 values = id.vals)
  still.unknown.ucsc <- setdiff(id.vals, entrezless.ensg.fetch$ucsc)
  # First get all the gene-less UCSC ID's but not found in Biomart from UCSC ID,
  # failing that get ENST IDs from the UCSC hg19.knownGenes-linked knownToEnsembl table (downloaded from table browser)
  ucsc.ensembl.table <- read.table('~/Downloads/hg19_missing_known_genes_enst.txt',
                                   sep='\t', header = TRUE, na.strings = 'n/a',
                                   stringsAsFactors = FALSE)[,c(1,4)]
  ucsc.enst.vals <- ucsc.ensembl.table[!is.na(ucsc.ensembl.table$hg19.knownToEnsembl.value),]
  still.unknown.table <- ucsc.ensembl.table[ucsc.ensembl.table$hg19.knownGene.name %in% still.unknown.ucsc,]
  # nrow(still.unknown.table) = 8550
  # then get union of rows with still.unknown.ucscs in the ucsc.enst.vals
  known.enst.table <- still.unknown.table[still.unknown.table$hg19.knownGene.name %in% ucsc.enst.vals$hg19.knownGene.name,]
  # multiple duplicated ENST transcripts (280 of 2643)
  # 8550 - 2643 = 5907 still unknown
  final.unknown.table still.unknown.table[known.enst.table$]
  
  tss <- genedf$txStart
  n.str <- genedf$strand == '-'
  tss[n.str] <- genedf$txEnd[n.str]
  tss.kb.us <- tss - 1000
  tss.kb.us[n.str] <- tss[n.str] + 1000
  tes <- genedf$txEnd
  tes[n.str] <- genedf$txStart[n.str]
  genedf$txStart <- tss
  genedf$txKbUs <- tss.kb.us
  genedf$txEnd <- tes
  
  us.asc.coords <- transform(genedf[c("txStart","txKbUs")],
                             min = ifelse(n.str, txStart, txKbUs),
                             max = ifelse(n.str, txKbUs, txStart))
  min.coords <- us.asc.coords$min
  max.coords <- us.asc.coords$max
  
  # Create a GRanges object for upstream of the TSSs
  tss.us <- GRanges(seqnames = Rle(genedf$chrom),
                    strand = Rle(genedf$strand),
                    ranges = IRanges(min.coords, max.coords),
                    TxID = genedf$tx_id)
}

### Get stats against ATAC seq experimental/control samples

for (peakfile in list.files(peak_filepath <- '../../peaks/atac', '*.narrowPeak')) {
  peak_name <- file_path_sans_ext(peakfile) %>%
    gsub("_peaks", ".peaks", .)
  if (!exists(peak_name)) assign(peak_name, readPeakFile(file.path(peak_filepath, peakfile)))
  AtacOverlap(peak_name)
  mapSitesToTSS()
  rm(list = as.character(peak_name)) # avoid large files building up in the workspace
}
gc() # force garbage collection to ensure no system crashes
