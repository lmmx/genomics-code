library('magrittr')
library('ChIPseeker') # also imports plyr, GenomicRanges
library('tools')
library('plyr')
library('ggplot2')
library('RColorBrewer')
library('dplyr')
library('reshape2')
library('seqinr')
library('stringr')
library('biomaRt')
library('data.table')
library('tidyr')
library('VennDiagram')

### Source files
if (!exists('assembly.report')) {
  assembly.of.interest <- 'GRCh37.p13' # matches assembly version of Ensembl genome database
  assembly.version.map <- read.table('/gits/grc-issues/assembly_reports/assembly_version_table.tsv',
                                     sep = '\t')
  # I'm not coding for this but check if 'RefSeq Assembly and GenBank Assemblies Identical' in metadata = yes
  # ...in fact it's only MalDomGD1.0 (Malus_x_domestica_GoldenDelicious1.0) which is not identical
  # assume not using apple genome and just use the genbank assembly report rather than both
  
  genbank.assembly.report.name <- assembly.version.map[assembly.version.map[,2] == assembly.of.interest,1] %>%
    as.character %>%
    extract2(1)
  
  #   genbank.assembly.report <- read.table(paste0('/gits/grc-issues/assembly_reports/',genbank.assembly.report.name),
  #                                         sep = '\t', stringsAsFactors = FALSE)
  
  # the TSV files use `#` as comment character and in text... pseudoautosomal regions e.g. `PAR#1` break the TSV, it's fine:
  assembly.report.con <- file(paste0('/gits/grc-issues/assembly_reports/', genbank.assembly.report.name))
  assembly.report.lines <- readLines(assembly.report.con)
  comment.lines <- assembly.report.lines %>% substr(1,1) == '#'
  non.comment.lines <- assembly.report.lines[!comment.lines]
  report.table.headers <- which(comment.lines) %>%
    rev %>%
    extract2(1) %>%
    assembly.report.lines[.] %>%
    substring(3)
  assembly.report.table <- read.table(text=c(report.table.headers, non.comment.lines),
                                      sep="\t", comment.char = '', header = TRUE) 
  close(assembly.report.con)
  
  patch.genbank.ids <- assembly.report.table$Scaffold.GenBank.Accn %>% sort %>% unique
}
if (!exists('contigs.to.genbank')) { contigs.to.genbank <- read.table('/gits/ensembl-assembly-exceptions/tables/GenBank_IDs/GRCh37_ensembl2gencode.txt',
                                                                      sep = '\t', stringsAsFactors = FALSE)
}
if (!exists('grc.issues')) { grc.issues <- read.table('/gits/grc-issues/tables/human_summary.tsv',
                                                      sep = '\t', stringsAsFactors = FALSE, header = TRUE)
                             colnames(grc.issues) <- c('grc', 'chr', 'genbank')
}
if (!exists('grc.locations')) {
  grc.locations <- read.table('/gits/grc-issues/tables/human_locations.tsv',
                              sep = '\t', stringsAsFactors = FALSE, header = TRUE)
}
if (!exists('tf.peaks')) {
  # Use the reference list of peaks provided:
  tf.peakfile <- '../../peaks/H1_FOXK2_strongest_peaks.narrowPeak'
  tf.peaks <- readPeakFile(file.path(tf.peakfile))
  seqlevels(tf.peaks) <- paste0('chr', seqlevels(tf.peaks))
  qc.tf.peaks <- tf.peaks * -0.4 # scale to 40% width, as in Davie 2015 (10.1371/journal.pgen.1004994)
}
if (!exists('genome.contigs')) { 
  genome.contigs <- c(1:22, "X", "Y", "MT") # aim to determine which of these contig is in
}

### Function Declarations

# iterate over guts of analysis in this function when finished drafting
AtacOverlap <- function(peak_name) {
  atac.peaks <<- get(peak_name)
  
  # Number of hits, including duplicated TF peaks (1-to-many TFBS-to-ATAC peak mapping)
  
  # length(which(duplicated(ranges(ref.peaks[subjectHits(findOverlaps(atac.peaks, ref.peaks))]))))
  
  qc.tfbs.at.atac <<- subsetByOverlaps(qc.tf.peaks, atac.peaks)
  n.untrimmed.overlaps <<- length(findOverlaps(tf.peaks, atac.peaks))
  n.unique.untrimmed.overlaps <<- length(subsetByOverlaps(tf.peaks, atac.peaks))
  n.qc.overlaps <<- length(findOverlaps(qc.tf.peaks, atac.peaks))
  n.unique.qc.overlaps <<- length(qc.tfbs.at.atac)
  
  qc.tfbs.stats <<- rbind(qc.tfbs.stats,
    data.frame(uq.qc.over = n.unique.qc.overlaps,
       qc.over = n.qc.overlaps,
       uq.over = n.unique.untrimmed.overlaps,
       over = n.untrimmed.overlaps,
       stringsAsFactors = FALSE))
  
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
         ") without QC trimming.\n") %>%
    cat()
  
  sample.name <- peak_name %>% gsub(".peaks", "", .)
  
  sample.atac.lists[[sample.name]] <<- list(atac.peaks %>% ranges)
  
  mapSitesToTSS(sample.name)
}

mapSitesToTSS <- function(sample.name) {
  seqlevels(qc.tfbs.at.atac) <- seqlevels(qc.tfbs.at.atac) %>% str_split(pattern = 'chr') %>% vpluck(2)
  
  tss.with.tfbs <- subsetByOverlaps(tss.us.ranges, qc.tfbs.at.atac)
  # this object is useless in terms of biological meaning, the range needs to be reverted to a TSS coordinate
  
  tfbs.at.tss <- subsetByOverlaps(qc.tfbs.at.atac, tss.us.ranges)
  # can match these (3010) back to the TSS if you keep this object along with the mapping (Hits) object
  
  peak.tss.mapping <- findOverlaps(qc.tfbs.at.atac, tss.us.ranges)
  tss.sites <- tss.us.ranges
  ranges(tss.sites) <- ( ifelse(strand(tss.sites) == '+', start(tss.sites), end(tss.sites)) %>% IRanges(start = ., end = .) )
  
    #   system.time(( if(strand(tss.sites) == '+') {
    #     min(ranges(tss.sites)) %>% IRanges(start = ., end = .) 
    #   }
    #   else {
    #     max(ranges(tss.sites)) %>% IRanges(start = ., end = .)
    #   }))
  
  #  peak.tss.mapping %>% queryHits %>% unique %>% length
  #  [1] 3010
  #  matches the number of tfbs.at.tss
  
  #  > peak.tss.mapping %>% subjectHits %>% unique %>% length
  #  [1] 19095
  
  #   tss.sites[peak.tss.mapping %>% queryHits %>% unique]$genebiotype %>% sort %>% table
  #   
  #   IG_C_gene        IG_D_gene        IG_J_gene       IG_V_gene polymorphic_pseudogene   protein_coding    TR_C_gene     TR_D_gene 
  #   24               53               20              153       172                      2383              6             3 
  #   TR_J_gene        TR_V_gene 
  #   68               128
  
  peak.tss.mapping <- findOverlaps(qc.tfbs.at.atac, tss.us.ranges)
  pc.peak.tss.map <- peak.tss.mapping[ (tss.us.ranges[peak.tss.mapping %>% subjectHits]$genebiotype == 'protein_coding') ]
  cat(paste0("Only using protein coding genes: discarding ", (peak.tss.mapping %>% length) - (pc.peak.tss.map %>% length), ' of ', peak.tss.mapping %>% length,
      ' (', ((1 - ( (pc.peak.tss.map %>% length) / (peak.tss.mapping %>% length))) * 100) %>% round(digits = 2), '%)'))
  
  distance.table <- tss.sites[pc.peak.tss.map %>% subjectHits]
  distance.table$peak_index <- Rle(pc.peak.tss.map %>% queryHits)
  distance.table$peak_range <- IRanges(qc.tfbs.at.atac[pc.peak.tss.map %>% queryHits] %>% ranges)
  distance.table$distance <- distance(distance.table %>% ranges, distance.table$peak_range)
  
  melted.gene.distances <- melt(as.data.frame(distance.table)[,c('seqnames','gene','distance')], id.vars = c('seqnames',"gene","distance"))
  grouped.gene.distances <- group_by(melted.gene.distances, gene)
  # get the minimum distance for each gene and store it as an Rle back in the distance table
  min.distances <- summarise(grouped.gene.distances, min=min(distance)) %>%
    slice(match(unique(distance.table$gene), gene))
  
  distance.table$min_dist <- join(distance.table[,'gene'] %>% as.data.frame, min.distances, by = 'gene')$min
  
  min.distance.table <- distance.table[distance.table$distance == distance.table$min_dist][,c('gene','min_dist')]
  gene.min.distance.table <- min.distance.table[!duplicated(min.distance.table$gene),]
  
  print(paste0(sample.name," has ",gene.min.distance.table$gene %>% length," ENSG-identified genes with a peak within 1kb u/s"))
  sample.gene.lists[[sample.name]] <<- list(gene.min.distance.table$gene)
}

vpluck <- function(x, i) vapply(x, "[[", i, FUN.VALUE = x[[1]][[i]])
ContigToGenBank <- function(contig.name, accession.map = contigs.to.genbank) {
  accession.map[which(accession.map[,1] == contig.name),2]
}
IsUnplacedScaffold <- function(gb.id) {
  # technically unplaced or 'unlocalised' but either way cannot obtain TSS reference coordinates if so
  if (substr(gb.id,1,5) == 'GL000') {
    # range (for GRCh37!) GL000[191.1 - 210.1] is unlocalised, [211.1 - 249.1] is unplaced
    return ((substr(gb.id,6,8) %>% as.integer) %>% between(191, 249)) # TRUE if in the above range
  } else { return(FALSE) }
}

PatchIdToChr <- function(genbank.id, accession.map = assembly.report.table) {
  grc.chr <- accession.map[which(accession.map[, 'Scaffold.GenBank.Accn'] == genbank.id), 'Chromosome'] %>%
    sort %>% unique %>% intersect(genome.contigs)
  if (length(grc.chr) > 1) {
    warning("Unable to map via GRC issues or locations tables, multiple chromosome names ( ", paste0(grc.chr, collapse = ', ')," ) returned for GenBank ID: '",
            genbank.id, "' (via contig name: '", get('contig', envir = parent.frame(2)), "').")
    return(NA) # looked at what's there, and 11 mostly unlocalised transcripts will not be missed in 180,000 rows
  } else if (length(grc.chr) < 1) {
    warning("Unable to map, no chromosome name returned for GenBank ID: '",
            genbank.id, "' (via contig name: '", get('contig', envir = parent.frame(2)), "')",
            " (this should never happen since the script checked it was in the list of GenBank IDs - check assembly report.). ",
            "Note: the mappings at github.com/dpryan79/ChromosomeMappings/blob/master/GRCh37_gencode2ensembl.txt",
            " were 'curated manually', see github.com/dpryan79/ChromosomeMappings/issues/7 for more information.")
    return(NA)
  } else { return(grc.chr) }
}

# PatchIdToChr <- function(genbank.id, accession.map = grc.locations) {
#   grc.chr <- accession.map[which(accession.map[,'MappedSeqInfo.GenBankID'] == genbank.id), 'chr'] %>%
#     sort %>% unique %>% intersect(genome.contigs)
#   if (length(grc.chr) > 1) {
#     warning("Unable to map via GRC issues or locations tables, multiple chromosome names ( ", paste0(grc.chr, collapse = ', ')," ) returned for GenBank ID: '",
#          genbank.id, "' (via contig name: '", get('contig', envir = parent.frame(2)), "')",
#          ". Correct the GRC issues/locations table(s), or remove the transcript in question before using this function.")
#     return(NA) # looked at what's there, and 11 mostly unlocalised transcripts will not be missed in 180,000 rows
#   } else if (length(grc.chr) < 1) {
#        warning("Unable to map via GRC locations table, no chromosome name returned for GenBank ID: '",
#                genbank.id, "' (via contig name: '", get('contig', envir = parent.frame(2)), "')",
#                ". Correct the GRC issues/locations table(s), or remove the transcript in question before using this function. ",
#                "Note: the mappings at github.com/dpryan79/ChromosomeMappings/blob/master/GRCh37_gencode2ensembl.txt",
#                " were 'curated manually', see github.com/dpryan79/ChromosomeMappings/issues/7 for more information.")
#        return(NA)
#   } else { return(grc.chr) }
# }
GenBankToChr <- function(genbank.id, accession.map = grc.issues) {
  # must return single value not a vector, via unique
  # will error out if two issues for same transcript point to different chromosomes
  if (genbank.id %>% IsUnplacedScaffold) { return(NA) }
  if (genbank.id %in% patch.genbank.ids) { return(PatchIdToChr(genbank.id)) }
  grc.chr <- accession.map[which(accession.map[,'genbank'] == genbank.id), 'chr'] %>%
    sort %>% unique
  if (length(grc.chr) == 1) { return(grc.chr) }
  if (length(grc.chr) > 1) {
    stop("Oops! Duplicate GenBank ID in the GRC issues table, multiple chromosome names ( ", paste0(grc.chr, collapse = ', '),
         " ) returned for GenBank ID: '", genbank.id, "' (via contig name: '", get('contig', envir = parent.frame()), "')",
         ". Correct the GRC issues table, or remove the transcript in question before using this function.",
         " (I checked and this should not happen for GRCh37.p13... debug using the commented out code below this warning in source.)")
    return(NA) # looked at what's there, and 11 mostly unlocalised transcripts will not be missed in 180,000 rows
    #     test.var <- grc.issues$genbank[!grc.issues$genbank %in% c('na', 'None')] %>% sort
    #     dup.test <- test.var[grc.issues$genbank[!grc.issues$genbank %in% c('na', 'None')] %>% sort %>% duplicated]
    #     in.dups <- grc.issues[grc.issues$genbank %in% dup.test,]
    #     in.dups[with(in.dups, order('genbank')),]
  } else {
    # grc.chr is empty
    
    #     if (genbank.id %in% grc.locations$MappedSeqInfo.GenBankID %>% sort %>% unique) {
    #       return(PatchIdToChr(genbank.id))
    #     }
    return(NA) # looked at what's there, and 11 mostly unlocalised transcripts will not be missed in 180,000 rows
    #     warning("Unable to map via GRC issues table, no chromosome name returned for GenBank ID: '",
    #          genbank.id, "' (via contig name: '", get('contig', envir = parent.frame()), "')",
    #          ". Correct the GRC issues table, or remove the transcript in question before using this function. ",
    #          "Note: the mappings at github.com/dpryan79/ChromosomeMappings/blob/master/GRCh37_gencode2ensembl.txt",
    #          " were 'curated manually', see github.com/dpryan79/ChromosomeMappings/issues/7 for more information.")
  }
}
MapContig <- function(contig) {
  if (! contig %in% genome.contigs) {
    mapped.to.chr <- contig %>% ContigToGenBank %>% GenBankToChr # not piped for debugging purposes
    return (mapped.to.chr)
  } else {
    return (contig)
  }
}
ParseHeaderField <- function(header_vector, col_sep_i, space_sep_j) {
  return (header_vector[col_sep_i] %>%
            str_split(pattern = ":") %>%
            vpluck(space_sep_j))
}
ParseFastaHeader <- function(header) {
  spl.annot <- str_split(header, " ")[[1]]
  return (list(
    # status is always cDNA since source is Ensembl cDNA dataset [release 74, GRCh37.p13]
    location = ParseHeaderField(spl.annot, 3, 1),
    contig = ParseHeaderField(spl.annot, 3, 3) %>% MapContig, # NB if not found returns NA, will not become a level
    # store temporary variables within list assignment so txStart and txEnd can be set conditionally by strand
    min = current.min <- as.integer(ParseHeaderField(spl.annot, 3, 4)),
    max = current.max <- as.integer(ParseHeaderField(spl.annot, 3, 5)),
    strand = current.strand <- ifelse(is.p.str <- (ParseHeaderField(spl.annot, 3, 6) == '1'), '+', '-'),
    ensg = ParseHeaderField(spl.annot, 4, 2),
    gene.biotype = ParseHeaderField(spl.annot, 5, 2),
    transcript.biotype = ParseHeaderField(spl.annot, 6, 2),
    patch.id = ifelse((current.patch.candidate <- ParseHeaderField(spl.annot, 3, 3)) %in% genome.contigs, NA, current.patch.candidate),
    tx.start = current.tx.start <- ifelse(is.p.str, current.min, current.max),
    tx.end = ifelse(is.p.str, current.max, current.min),
    tx.kb.us = current.tx.start + ifelse(is.p.str, -1000, +1000)
  ))
}
ParseTranscript <- function(transcript) {
  # Increment counter of progress through transcript set [when used via ParseTranscriptSet]
  vlen <- get('transcripts.vlen', envir = globalenv())
  if (transcript.counter == vlen) {
    cat('100% !')
  } else if (transcript.counter %% (round(vlen/100) * 5) == 0) {
    cat(paste0(round(transcript.counter*100/vlen, digits = 0), '%')) # 5%
  } else if (transcript.counter %% round(vlen/100) == 0) {
    cat('.') # each of single % up to 5-multiple
  }
  transcript.counter <<- get('transcript.counter', envir = globalenv()) + 1
  # each %
  transcript %>%
    getAnnot %>%
    ParseFastaHeader %>%
    c(transcript = transcript %>% getName, .) %>%
    as.data.frame(stringsAsFactors = FALSE)
}
ParseTranscriptSet <- function(transcripts) {
  transcript.counter <<- 1
  transcripts.vlen <<- length(transcripts)
  parsed.transcripts <- lapply(transcripts, ParseTranscript) %>%
    unname %>%
    rbindlist %>%
    as.data.frame
  # only use the rows on which contig !is.na
  # (the unlocalised/unplaced ones in GL000 ranges, deliberately returned as such)
  # those are also the only rows on which location is 'supercontig' rather than 'chromosome'
  return( parsed.transcripts[!parsed.transcripts$contig %>% is.na, ] )
}

if (!exists('contig.map')) {
  odd.contigs <- genome.cdna[!(getAnnot(genome.cdna) %>% str_split(':') %>% vpluck(4) %in% genome.contigs)] %>% getAnnot %>% str_split(':') %>% vpluck(4)
  odd.contig.list <- odd.contigs[!duplicated(odd.contigs)]
  contig.map <- data.frame(contig = odd.contig.list[! (odd.contig.list %>% lapply(IsUnplacedScaffold) %>% unlist)],
                           stringsAsFactors = FALSE)
  contig.map$genbank <- contig.map$contig %>% lapply(ContigToGenBank) %>% unlist
}

# Set up a reference to the transcripts in the hg19 reference genome (for TSS coordinates)
if (!exists('genedf')) {
  
  ## COMMENT OUT THIS BLOCK TO OVERWRITE THE READ GENEDF FROM FILE
  
  if (file.exists('genedf_5kb.rds')) {
    genedf <- readRDS('genedf_5kb.rds')
#   if (file.exists('genedf_1kb.rds')) {
#     genedf <- readRDS('genedf_1kb.rds')
  } else {
  
  ## COMMENT OUT THIS BLOCK TO OVERWRITE THE READ GENEDF FROM FILE
    
    if (!exists('genome.cdna')) genome.cdna <- read.fasta('../../hg19/ensembl/Homo_sapiens.GRCh37.74.cdna.all.fa')
    genedf <- ParseTranscriptSet(genome.cdna)
   }
  is.p.str <- genedf$strand == '+'
  ascending.coords <- transform(genedf[c("tx.start","tx.kb.us")],
                                range.min = ifelse(is.p.str, tx.kb.us, tx.start),
                                range.max = ifelse(is.p.str, tx.start, tx.kb.us))
  
  # Create a GRanges object for upstream of the TSSs
  tss.us.ranges <- GRanges(seqnames = Rle(genedf$contig),
                           ranges = IRanges(ascending.coords$range.min, ascending.coords$range.max),
                           strand = Rle(genedf$strand),
                           transcript = genedf$transcript, # ENST transcript ID
                           gene = genedf$ensg,
                           genebiotype = genedf$gene.biotype,
                           txbiotype = genedf$transcript.biotype)
  
  saveRDS(genedf, 'genedf.rds') # save having to regenerate from source each time you reload project
}

### Get stats against ATAC seq experimental/control samples

sample.gene.lists <- list() # sample-to-gene list mapping
sample.atac.lists <- list() # sample-to-ATAC-open-region-count list mapping
qc.tfbs.stats <- data.frame() # sample-to-QC-trimming-count-statistics mapping

for (peakfile in list.files(peak_filepath <- '../../peaks/atac', '*.narrowPeak')) {
  peak_name <- file_path_sans_ext(peakfile) %>%
    gsub("_peaks", ".peaks", .)
  if (!exists(peak_name)) assign(peak_name, readPeakFile(file.path(peak_filepath, peakfile)))
  AtacOverlap(peak_name)
  rm(list = as.character(peak_name)) # avoid large files building up in the workspace
}
gc() # force garbage collection to ensure no system crashes

# OPTIONAL MEMORY SAVING RM STATEMENT

# rm('ascending.coords', 'assembly.report.table', 'assembly.version.map',
#    'genedf', 'atac.peaks', 'tss.us.ranges', 'tf.peaks', 'qc.tf.peaks')

paste0("Mean values were: ",
       round(qc.tfbs.stats$uq.qc.over %>% mean),
       " (",
       round ( 100 * ( (qc.tfbs.stats$uq.qc.over / qc.tfbs.stats$qc.over) %>% mean ), digits = 1 ),
       "% unique of ",
       round(qc.tfbs.stats$qc.over %>% mean),
       " QC-trimmed) TFBS in ATAC-defined 'accessible regions', or ",
       round(qc.tfbs.stats$uq.over %>% mean),
       " (",
       round ( 100 * ( (qc.tfbs.stats$uq.over / qc.tfbs.stats$over) %>% mean) , digits = 1 ),
       "% unique of ",
       round(qc.tfbs.stats$over %>% mean),
       ") without QC trimming.\n") %>%
  cat()

gc.length.list <- list()
for (list.item in names(sample.gene.lists)) {
  gc.length.list[[list.item]] <- length(unlist(sample.gene.lists[[list.item]]))
}

atac.length.list <- list()
for (list.item in names(sample.atac.lists)) {
  atac.length.list[[list.item]] <- unlist(sample.atac.lists[[list.item]])[[1]] %>% length
}

std.error <- function(x) sd(x)/sqrt(length(x))

atac.width.list <- list()
atac.se.list <- list()
for (list.item in names(sample.atac.lists)) {
  atac.width.list[[list.item]] <- unlist(sample.atac.lists[[list.item]])[[1]] %>% width %>% mean
  atac.se.list[[list.item]] <- unlist(sample.atac.lists[[list.item]])[[1]] %>% width %>% std.error
}

sample.names <- names(gc.length.list %>% unlist)
sample.names[1:9] <- c('1 (T)', '2 (N)', '3 (T)', '4 (N)',
                       '5 (T)', '6 (N)', '7 (T)', '8 (T)', '9 (T)')

sample.gene.list.length.df <- data.frame(sample = sample.names,
                                         gene.count = gc.length.list %>% unlist,
                                         group = sample.names %>%
                                           gsub(pattern = '_.*', '', x = .) %>%
                                           factor,
                                         row.names = NULL,
                                         stringsAsFactors = FALSE)

sample.atac.list.length.df <- data.frame(sample = sample.names,
                                         atac.count = atac.length.list %>% unlist,
                                         atac.width = atac.width.list %>% unlist,
                                         atac.width.se = atac.se.list %>% unlist,
                                         group = sample.names %>%
                                           gsub(pattern = '_.*', '', x = .) %>%
                                           factor,
                                         row.names = NULL,
                                         stringsAsFactors = FALSE)

sample.types <- data.frame(group = sample.gene.list.length.df$group %>% levels,
                           type = c('t','n','t','n','t','n',rep('t',3),rep('n',2),rep('t',3) ) )

sample.gene.list.length.df <- join(sample.gene.list.length.df, sample.types, by='group')
sample.atac.list.length.df <- join(sample.atac.list.length.df, sample.types, by='group')

# sample.gene.list.length.df <- data.frame(gene.count = gc.length.list %>% unlist %>% sort,
#                                          stringsAsFactors = FALSE)

# colourCount <- sample.gene.list.length.df$sample.group %>% levels %>% length
# getPalette = colorRampPalette(brewer.pal(9, cbPalette))
modAccentPalette <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99",
                      "#386CB0", "#dd3c3c", "#5af2a3", "#666666")
cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#b87ab8", "#D55E00", "#CC79A7")
getPalette <- c(modAccentPalette, cbPalette)

# Possible figure: A bar plot of the samples grouped by cell type
#   - no real differences shown between samples (though OE33_R3 looks out?)
#   - no real difference observable between tumour and normal (could be quantified/tested?)

plot.theme <- theme_bw() + theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

sample_gene_count_bar <- ggplot(sample.gene.list.length.df,
                                aes(x=sample,
                                    y=gene.count,
                                    fill=group)) + 
              geom_bar(stat="identity") +
              labs(fill='Cell line \n(if specified)\n\nT = tumour,\nN = normal\n',
                   x='Sample', y='Gene count') +
              scale_fill_manual(values = getPalette) + plot.theme +
              expand_limits(y=0) + scale_y_continuous(expand = c(0, 0))

# Change colour scheme for second version, comparing tumour vs. normal samples
sample_type_gene_count_bar <- ggplot(sample.gene.list.length.df,
                                     aes(x=sample,
                                         y=gene.count,
                                         fill=factor(type, labels=c('Normal','Tumour')))) +
                              geom_bar(stat="identity") +
                              labs(fill='Sample type',
                                   x='Sample', y='Gene count') + plot.theme +
                              expand_limits(y=0) + scale_y_continuous(expand = c(0, 0)) +
                              scale_fill_manual(values = brewer.pal(3, 'Set2'))

# TASK: Test whether numbers of genes between tumour and normal differ:
melt(sample.gene.list.length.df) %>%
  group_by(type) %>%
  summarise(n=length(value), mean=mean(value), sd(value))
# Using sample, group, type as id variables
# Source: local data frame [2 x 4]
# 
#   type       n        mean       sd(value)
#  (fctr)    (int)      (dbl)       (dbl)
#    n         7      8169.571     1332.282
#    t        10      7973.800     1824.760
t.test(sample.gene.list.length.df$gene.count~sample.gene.list.length.df$type)
# data:  sample.gene.list.length.df$gene.count by sample.gene.list.length.df$type
# t = 0.25562, df = 14.935, p-value = 0.8017
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1437.241  1828.784
# sample estimates:
#   mean in group n mean in group t 
# 8169.571        7973.800 
# 
# > t.test(sample.gene.list.length.df$gene.count~sample.gene.list.length.df$type)$p.value
# [1] 0.8017292

#TASK: test commonality of regulated genes between/within cancer and regular samples
all.in.all.cancer.genes <- sample.gene.lists[(sample.gene.list.length.df$type == 't') %>% which] %>% 
  unlist %>% unique
oes.cancer.genes <- sample.gene.lists[(sample.gene.list.length.df$group %in% c('OE19','OE33')) %>% which] %>%
  unlist %>% unique
all.in.all.normal.genes <- sample.gene.lists[(sample.gene.list.length.df$type == 'n') %>% which] %>% 
  unlist %>% unique
common.genes <- (all.in.all.cancer.genes)[(all.in.all.cancer.genes) %in% (all.in.all.normal.genes) %>% which]
common.oes.cancer.genes <- (oes.cancer.genes)[(oes.cancer.genes) %in% (all.in.all.cancer.genes) %>% which]
pc.unique.cancer <- (1 - ( common.genes %>% length / all.in.all.cancer.genes %>% length) ) * 100
pc.unique.oes <- (1 - ( common.oes.cancer.genes %>% length / all.in.all.cancer.genes %>% length) ) * 100
pc.unique.normal <- (1 - ( common.genes %>% length / all.in.all.normal.genes %>% length) ) * 100
cat(paste0('Percent of genes unique to cancer samples: ', pc.unique.cancer %>% round(digits = 1), '% (', 
             (all.in.all.cancer.genes %>% length) - (common.genes %>% length), ' of ', (all.in.all.cancer.genes %>% length), ')\n',
             'Percent of genes unique to normal samples: ', pc.unique.normal %>% round(digits = 1), '% (', 
             (all.in.all.normal.genes %>% length) - (common.genes %>% length), ' of ', (all.in.all.normal.genes %>% length), ')'))

# With 40% ATAC peak width QC trimming at 5kb upstream range:
# Percent of genes unique to cancer samples: 5.6% (594 of 10585)
# Percent of genes unique to normal samples: 3.2% (333 of 10324)

# With 40% ATAC peak width QC trimming at 1kb upstream range:
# Percent of genes unique to cancer samples: 3.7% (380 of 10331)
# Percent of genes unique to normal samples: 2.4% (247 of 10198)

sample_open_region_count_bar_by_type <- ggplot(sample.atac.list.length.df,
                                     aes(x=sample,
                                         y=atac.count,
                                         fill=factor(type, labels=c('Normal','Tumour')))) +
  geom_bar(stat="identity") +
  labs(fill='Sample type',
       x='Sample', y='ATAC accessible chromatin regions count') + plot.theme +
  expand_limits(y=0) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = brewer.pal(3, 'Set2'))

# TASK: Test whether numbers of ATAC open regions between tumour and normal differ:
melt(sample.atac.list.length.df) %>%
  group_by(type) %>%
  summarise(n=length(value), mean=mean(value), sd(value))
# Using sample, group, type as id variables
# Source: local data frame [2 x 4]
# 
# type     n     mean sd(value)
# (fctr) (int)    (dbl)     (dbl)
# 1      n     7 57232.57  42682.45
# 2      t    10 57640.50  40659.57
t.test(sample.atac.list.length.df[10:17,]$atac.count~sample.atac.list.length.df[10:17,]$type)
# Welch Two Sample t-test
# 
# data:  sample.atac.list.length.df$atac.count by sample.atac.list.length.df$type
# t = -0.019774, df = 12.643, p-value = 0.9845
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -45103.60  44287.74
# sample estimates:
#   mean in group n mean in group t 
# 57232.57        57640.50 

#TASK: above but with ATAC peak widths
sample_open_region_width_bar_by_type <- ggplot(sample.atac.list.length.df,
                                               aes(x=sample,
                                                   y=atac.width,
                                                   fill=factor(type, labels=c('Normal','Tumour')))) +
  geom_bar(stat="identity") +
  geom_errorbar(stat="identity", aes(ymin=atac.width-atac.width.se,
                                     ymax=atac.width+atac.width.se,
                                     width=.3,)) +
  labs(fill='Sample type',
       x='Sample', y='ATAC accessible chromatin regions, mean width') + plot.theme +
  expand_limits(y=0) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = brewer.pal(3, 'Set2'))
