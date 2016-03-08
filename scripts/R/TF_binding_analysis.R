library('magrittr')
library('ChIPseeker') # also imports plyr, GenomicRanges
library('tools')
library('dplyr')
library('seqinr')
library('stringr')
library('biomaRt')
library('data.table')
library('tidyr')

### Source files
if (!exists('genome.cdna')) { 
  genome.cdna <- read.fasta('../../hg19/ensembl/Homo_sapiens.GRCh37.74.cdna.all.fa')
}
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
  # Use the reference list of peaks Yaoyong provided:
  tf.peakfile <- '../../peaks/H1_FOXK2_strongest_peaks.narrowPeak'
  tf.peaks <- readPeakFile(file.path(tf.peakfile))
  seqlevels(tf.peaks) <- paste0('chr', seqlevels(tf.peaks))
}
if (!exists('genome.contigs') { 
  genome.contigs <- c(1:22, "X", "Y", "MT") # aim to determine which of these contig is in
}

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
  lapply(transcripts, ParseTranscript) %>%
    unname %>%
    rbindlist %>%
    as.data.frame
}

if (!exists('contig.map')) {
  odd.contigs <- genome.cdna[!(getAnnot(genome.cdna) %>% str_split(':') %>% vpluck(4) %in% genome.contigs)] %>% getAnnot %>% str_split(':') %>% vpluck(4)
  odd.contig.list <- odd.contigs[!duplicated(odd.contigs)]
  contig.map <- data.frame(contig = odd.contig.list[! (odd.contig.list %>% lapply(IsUnplacedScaffold) %>% unlist)],
                           stringsAsFactors = FALSE)
  contig.map$genbank <- contig.map$contig %>% lapply(ContigToGenBank) %>% unlist
}

if (!exists('transcriptdf')) {
  transcriptdf <- InitialParseTranscriptSet(genome.cdna[1:1000])
}

InitialParseTranscriptSet <- function(input.transcripts) {
  data.frame(transcripts = names(input.transcripts),
             contig = input.transcripts %>% getAnnot %>% str_split(":") %>% vpluck(4),
  )
  ContigToGenBank(x)
}

# Set up a reference to the transcripts in the hg19 reference genome (for TSS coordinates)
if (!exists('genedf')) {
  genedf <- ParseTranscriptSet(genome.cdna)
  # NB only use the rows on which contig !is.na (the unlocalised/unplaced ones in GL000 ranges, deliberately returned as such)
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
  
  unmapped.chr.transcripts <- genedf[genedf$contig %>% is.na,]
  unmapped.chr.transcripts$contig.name <- genedf[genedf$contig %>% is.na,] %>% extract2('transcript') %>%
    lapply(function(x){
      genome.cdna %>% extract2(x) %>% getAnnot %>% str_split(':') %>% vpluck(4)
    }) %>% unlist
  # unmapped.chr.transcripts$contig.name %>% lapply(function(x){(grc.issues$genbank == x) %>% any}) %>% unlist %>% any
  # ==> FALSE
  unmapped.chr.table <- unmapped.chr.transcripts$contig.name %>%
    lapply(function(x){
      grc.locations[which(grc.locations$MappedSeqInfo.GenBankID == x),]
    }) %>%
    do.call(what="rbind", args=.)
  # unmapped.chr.table %>% extract2('id') %>% sort %>% unique %>% length
  # ==> 20 unique unmappable IDs
  # unmapped.chr.table %>% extract2('id') %>% sort %>% unique
  # ==> "HG-1013" "HG-1100" "HG-1103" "HG-1104" "HG-1284" "HG-1285" "HG-1659" "HG-1662" "HG-1737" "HG-1803" "HG-1804" "HG-1847" "HG-1959" "HG-2011" "HG-255"  "HG-630"  "HG-660" 
  #     "HG-661"  "HG-789"  "HG-791" 
  
  # It's possible there were flaws in the way I handled chromosome assignment for the ones mapped through 'backup' function PatchIdToChr
  # which figuring out the representation of the correct contig placement per assembly (i.e. of issue encoding) may clarify, i.e. these
  # edge cases may improve the whole set, so worth doing rather than just throwing them away...
  
  # i.e. the method above may simply discard the Un option for pairs (e.g. "chr4 or chrUn"), when it may be the case GRCh37.p13 had removed the chromosome 4 mapping.
  # Clarifying the specification format per assembly version would verify it (some "chr4" cases in the above example may become "chrUn"), and avoid this problem.
  
  # Alternatively, it's possibly not decipherable directly from the given tables, i.e. if the issues do not record this systematically,
  # in which case confirming the 9 GenBank IDs by hand would be fine... but then not possible to do for other releases...
  
  # unmapped.chr.table$MappedSeqInfo.GenBankID %>% sort %>% unique
  # ==> "GL000193.1" "GL000194.1" "GL000213.1" "GL000215.1" "GL000218.1" "GL000219.1" "GL000221.1" "GL000222.1" "GL000223.1"
  
  # To compare the location data associated with a single GenBank accession [per transcript at issue]
  # echo "$(head -1 human_locations.tsv; grep 'GL000194.1' human_locations.tsv)" | cutf 1,2,3,4,6,7,8,13,16 | tsvview
  
  # For all unmapped GenBank IDs (9):
  # ids=("GL000193.1" "GL000194.1" "GL000213.1" "GL000215.1" "GL000218.1" "GL000219.1" "GL000221.1" "GL000222.1" "GL000223.1")
  # echo "$(head -1 human_locations.tsv; for gbid in ${ids[*]}; do grep $gbid human_locations.tsv; done | sort -k4,4 -k1,1)" | cutf 1,2,3,4,6,7,8,13,16 | tsvview
  
  # For all unmapped patches (20):
  # patchids=("HG-1013" "HG-1100" "HG-1103" "HG-1104" "HG-1284" "HG-1285" "HG-1659" "HG-1662" "HG-1737" "HG-1803" \
  #           "HG-1804" "HG-1847" "HG-1959" "HG-2011" "HG-255"  "HG-630"  "HG-660"  "HG-661"  "HG-789"  "HG-791")
  # ...View issue's longform locations:
  # echo "$(head -1 human_locations.tsv; for patchid in ${ids[*]}; do grep $patchid human_locations.tsv; done | sort -k4,4 -k1,1)" | cutf 1,2,3,4,6,7,8,13,16 | tsvview
  # ...View issue with shortform locations:
  # echo "$(head -1 human.tsv; for patchid in ${ids[*]}; do grep $patchid human.tsv; done)" | tsvview
  
  # For all locations:
  # echo "$(head -1 human_locations.tsv; tail --lines=+2 human_locations.tsv | sort -k4,4 -k1,1)" | cutf 1,2,3,4,6,7,8,13,16 | tsvview
  
  # Looking at multi-chromosome mapping IDs
  # manually get contig.map$genbank (no unplaced/unloc. contigs) onto clipboard, removing double quotes and iterating through bash string array:
  # gbids=`xclip -o`
  # for gbid in ${gbids[*]}; do echo "$(head -1 human_locations.tsv; grep "$gbid" human_locations.tsv)"  | tsvview; done
  
  # got a list from MANUAL INSPECTION(...) of list... for those with more than one chromosome
  
  # bad_patch_ids=(JH720453.7 JH636052.4 JH636053.3 KE332496.1 GL383561.2 KE332502.1)
  # for gbid in ${bad_patch_ids[*]}; do echo "$(head -1 human_locations.tsv; grep "$gbid" human_locations.tsv)"  | tsvview; done
  
  # in R:
  # grc.issues[which(grc.issues$genbank %in% c('JH720453.7', 'JH636052.4', 'JH636053.3', 'KE332496.1', 'GL383561.2', 'KE332502.1')),]
  # grc.locations[which(grc.locations$MappedSeqInfo.GenBankID %in% c('JH720453.7', 'JH636052.4', 'JH636053.3', 'KE332496.1', 'GL383561.2', 'KE332502.1')),]
  
  # properly in R, for each group of genedf transcript rows made by matching against GenBank ID:
  # list.of.bad.patch.transcript.locations <- lapply(contig.map$genbank, function(x) {
  #   patch.locations <- grc.locations[which(grc.locations$MappedSeqInfo.GenBankID == x),]
  #   if ((grc.locations$chr %>% sort %>% unique %>% length) > 1) {
  #     return(patch.locations)
  #   }
  # })
  # list.of.bad.patch.transcript.locations <- list.of.bad.patch.transcript.locations[lapply(list.of.bad.patch.transcript.locations, function(x){!is.null(x)}) %>% unlist]
  # list.of.bad.patch.transcript.locations[1]
  # ==> [[1]]
  #           id chr SuccessfullyMapped MappedSeqInfo.GenBankID MappedSeqInfo.RefSeqID MappedSeqInfo.SequenceType ChrStart  ChrEnd MappedVersionNumbers MappedVersionAccessions
  # 3645  HG-581   1             MAPPED              JH720453.1         NW_003871100.1               ALT_SCAFFOLD  1350798 1378514                  2,2       BX682528,BX682528
  # 6802 HG-1425   X             MAPPED              JH720453.1         NW_003871100.1               ALT_SCAFFOLD        1 1452651                  2,3       AC233982,CR407552
  #      Accession1Method Accession2Method AssemblyName GenBankAssemblyAcc RefSeqAssemblyAcc AssemblyStatus
  # 3645        alignment        alignment   GRCh37.p13   GCA_000001405.14  GCF_000001405.25            old
  # 6802        component        alignment   GRCh37.p13   GCA_000001405.14  GCF_000001405.25            old
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