library('GenomicRanges')
library('ChIPseeker')
library('tools')
library('magrittr')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')

### Initial setup

# Use the reference list of peaks Yaoyong provided:
ref.peakfile <- '../../peaks/H1_FOXK2_strongest_peaks.narrowPeak'
ref.peaks <- readPeakFile(file.path(ref.peakfile))
seqlevels(ref.peaks) <- paste0('chr', seqlevels(ref.peaks))

# Set up a reference to the transcripts in the hg19 reference genome (for TSS coordinates)
if (!exists txdb_hg19) {
  txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes <- as.data.frame(transcripts(txdb_hg19))
  colnames(genes)[1:3]=c("chrom","txStart","txEnd")
  tss = genes$txStart
  idx = genes$strand == '-'
  tss[idx] = genes$txEnd[idx]
  tes = genes$txEnd
  tes[idx] = genes$txStart[idx]
  genes$txStart <- tss
  genes$txEnd <- tes
  
  # Create a GRanges object for upstream of the TSSs
  tss.us <- GRanges(seqnames=Rle(genes$chrom), ranges=IRanges(tss-1000, tss))
}

### Loop through peak files

for (peakfile in list.files(peak_filepath <- '../../peaks/atac', '*.narrowPeak')) {
  peak_name <- file_path_sans_ext(peakfile) %>%
    gsub("_peaks", ".peaks", .)
  if (!exists(peak_name)) assign(peak_name, readPeakFile(file.path(peak_filepath, peakfile)))
  processPeakFile(peak_name)
  paste0(peak_name, " has ", qc.overlaps, " (of ", untrimmed.overlaps,  " untrimmed) overlapping peaks with the reference (comparing central 40% of each)") %>%
    print()
  mapPeaks()
}

# iterate over guts of analysis in this function when finished drafting
processPeakFile <- function(peak_name) {
  peaks <<- get(peak_name)
  qc.peaks <<- peaks * -0.4 # scale the coordinates to 40% width
  untrimmed.overlaps <<- sum(countOverlaps(ref.peaks, peaks))
  # equivalent to: length(findOverlaps(ref.peaks, peaks))
  qc.overlaps <<- sum(countOverlaps(ref.peaks, qc.peaks))
  
}

mapPeaks <- function() {
  qc.peaks[qc.peaks %over% tss.us]
}