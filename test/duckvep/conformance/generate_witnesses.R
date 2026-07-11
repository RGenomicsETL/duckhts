#!/usr/bin/env Rscript
# Generate variants at transcript boundaries and coding positions where rare
# consequence states occur.
#
# Ported from /root/duckvep-c@9f922c8. It reads
# a transcript's exon/CDS structure from a GFF and the reference bases from a FASTA range,
# then covers the positions where consequence bugs hide — each exon/intron
# boundary (donor / acceptor / 5th base / region), the start and stop codons, exon edges
# and interiors — crossed with allele shapes (three SNVs, one-base deletion/insertion,
# two-base MNV, two/three-base deletion, and short length-changing replacements). It
# only generates the VCF; the differential runner compares every result with VEP and
# preserves unsupported or unresolved cases.
#
# CLASS is a coarse grouping by genomic geometry (donor = intron 5' side, acceptor =
# intron 3' side), not the strand-resolved SO term. The actual
# SO term is decided by VEP and DuckVEP during comparison; CLASS only groups witnesses.
#
# air-formatted; strings via glue (project standard).
suppressMessages({
  library(optparse)
  library(duckdb)
  library(DBI)
  library(glue)
})
options(rlang_backtrace_on_error = "none")

root <- tryCatch(
  system2(
    "git",
    c("rev-parse", "--show-toplevel"),
    stdout = TRUE,
    stderr = FALSE
  ),
  error = function(e) "."
)
op <- OptionParser(
  usage = "%prog [options]  (defaults target the in-tree minimal fixture)"
)
op <- add_option(
  op,
  "--gff",
  default = file.path(root, "test/data/duckvep/minimal.gff3")
)
op <- add_option(
  op,
  "--fasta",
  default = file.path(root, "test/data/duckvep/minimal.fa")
)
op <- add_option(
  op,
  "--tx",
  default = "DUCK1-201",
  help = "target transcript id [%default]"
)
op <- add_option(
  op,
  "--out",
  default = file.path(
    root,
    "test/duckvep/conformance/results/witnesses.vcf"
  )
)
op <- add_option(
  op,
  "--ext",
  default = Sys.getenv(
    "DUCKHTS_EXT",
    file.path(root, "build/release/duckhts.duckdb_extension")
  )
)
op <- add_option(
  op,
  "--check",
  action = "store_true",
  default = FALSE,
  help = "self-test on the minimal fixture: assert expected witnesses, then exit"
)
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)

con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true")))
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
invisible(dbExecute(con, glue("LOAD '{opt$ext}'")))

# 1. transcript features (exon/CDS) for the target tx. Transcript line carries ID=transcript:X;
#    exon/CDS lines carry Parent=transcript:X.
feats <- dbGetQuery(
  con,
  glue(
    "SELECT seqname AS chrom, feature, start, \"end\" AS endp, strand,
          regexp_extract(attributes,'Parent=transcript:([^;]*)',1) AS parent_tx,
          regexp_extract(attributes,'ID=transcript:([^;]*)',1) AS id_tx
   FROM read_gff('{opt$gff}', scan_mode := 'sequential')"
  )
)
tx_line <- feats[feats$id_tx == opt$tx, ]
geom <- feats[feats$parent_tx == opt$tx & feats$feature %in% c("exon", "CDS"), ]
if (nrow(geom) == 0L) {
  die("no exon/CDS features for transcript '{opt$tx}' in {opt$gff}")
}
chrom <- geom$chrom[1]
strand <- if (nrow(tx_line)) tx_line$strand[1] else geom$strand[1]

exons <- geom[geom$feature == "exon", c("start", "endp")]
exons <- exons[order(exons$start), , drop = FALSE]
cds <- geom[geom$feature == "CDS", c("start", "endp")]
cds_lo <- if (nrow(cds)) min(cds$start) else NA_integer_
cds_hi <- if (nrow(cds)) max(cds$endp) else NA_integer_

# 2. class-relevant genomic POSITIONS (the covering-array offsets).
pos <- integer(0)
lab <- character(0)
add <- function(p, l) {
  pos <<- c(pos, as.integer(p))
  lab <<- c(lab, rep(l, length(p)))
}
for (i in seq_len(nrow(exons))) {
  s <- exons$start[i]
  e <- exons$endp[i]
  add(c(s, s + 1L, e - 1L, e), "exon_edge")
  if (i < nrow(exons)) {
    add(c(e + 1L, e + 2L, e + 5L), "donor")
  } # intron 5' (genomic)
  if (i > 1L) {
    add(c(s - 1L, s - 2L, s - 8L), "acceptor")
  } # intron 3' (genomic)
  add(as.integer(round((s + e) / 2)), "exon_mid")
}
if (!is.na(cds_lo)) {
  start_codon <- if (strand == "+") {
    cds_lo:(cds_lo + 2L)
  } else {
    (cds_hi - 2L):cds_hi
  }
  stop_codon <- if (strand == "+") {
    (cds_hi - 2L):cds_hi
  } else {
    cds_lo:(cds_lo + 2L)
  }
  add(start_codon, "start_codon")
  add(stop_codon, "stop_codon")
}

# dedup positions keeping the MOST SPECIFIC label (lower rank wins).
pos <- pmax(pos, 1L)
rnk <- c(
  start_codon = 1,
  stop_codon = 2,
  donor = 3,
  acceptor = 4,
  exon_edge = 5,
  exon_mid = 6
)
ord <- order(pos, rnk[lab])
pos <- pos[ord]
lab <- lab[ord]
keep <- !duplicated(pos)
pos <- pos[keep]
lab <- lab[keep]

# 3. reference bases over the position span (FASTA range; build the .fai if missing).
if (!file.exists(glue("{opt$fasta}.fai"))) {
  invisible(dbGetQuery(
    con,
    glue("SELECT * FROM fasta_index('{opt$fasta}')")
  ))
}
p0 <- min(pos)
span <- glue("{chrom}:{p0}-{max(pos) + 4L}")
seqstr <- toupper(dbGetQuery(
  con,
  glue(
    "SELECT SEQUENCE AS s FROM read_fasta('{opt$fasta}', region := '{span}') LIMIT 1"
  )
)$s[1])
nuc <- function(p, n) substr(seqstr, p - p0 + 1L, p - p0 + n)
base_at <- function(p) {
  b <- nuc(p, 1L)
  if (b %in% c("A", "C", "G", "T")) b else NA_character_
}
other_base <- function(...) {
  setdiff(c("A", "C", "G", "T"), c(...))[1]
}

# Tile allele shapes at every selected position. The differential runner decides scope.
W <- list()
emit <- function(p, ref, alt, class, shape) {
  W[[length(W) + 1L]] <<- data.frame(
    pos = p,
    ref = ref,
    alt = alt,
    class = class,
    shape = shape,
    stringsAsFactors = FALSE
  )
}
for (k in seq_along(pos)) {
  p <- pos[k]
  rb <- base_at(p)
  if (is.na(rb)) {
    next
  }
  for (a in setdiff(c("A", "C", "G", "T"), rb)) {
    emit(p, rb, a, lab[k], "snv")
  }
  tb <- nuc(p, 2L)
  if (grepl("^[ACGT]{2}$", tb)) {
    emit(p, tb, substr(tb, 1L, 1L), lab[k], "del1") # 1bp del (frameshift)
    emit(p, rb, glue("{rb}T"), lab[k], "ins1") # 1bp ins (frameshift)
    emit(p, rb, glue("{rb}ATG"), lab[k], "ins3") # 3bp in-frame insertion
    emit(p, tb, chartr("ACGT", "TGCA", tb), lab[k], "mnv2") # 2bp MNV (both flip)
    emit(
      p,
      tb,
      other_base(substr(tb, 1L, 1L), substr(tb, 2L, 2L)),
      lab[k],
      "delins2to1"
    )
    emit(
      p,
      rb,
      glue("{other_base(rb)}{other_base(rb, other_base(rb))}"),
      lab[k],
      "delins1to2"
    )
  }
  fb <- nuc(p, 3L)
  if (grepl("^[ACGT]{3}$", fb)) {
    emit(p, fb, substr(fb, 1L, 1L), lab[k], "del2")
  } # 2bp del
  qb <- nuc(p, 4L)
  if (grepl("^[ACGT]{4}$", qb)) emit(p, qb, substr(qb, 1L, 1L), lab[k], "del3") # 3bp del (in-frame)
}
w <- do.call(rbind, W)
w <- w[order(w$pos, w$ref, w$alt), , drop = FALSE]
w$variant_id <- glue_data(w, "{chrom}:{pos}:{ref}:{alt}")

if (opt$check) {
  stopifnot(
    "no witnesses generated" = nrow(w) > 0L,
    "missing donor class" = any(w$class == "donor"),
    "missing exon_edge class" = any(w$class == "exon_edge"),
    "missing start_codon class (coding tx)" = is.na(cds_lo) ||
      any(w$class == "start_codon"),
    "all SNV refs single ACGT base" = all(
      w$ref[w$shape == "snv"] %in% c("A", "C", "G", "T")
    ),
    "SNV alt differs from ref" = all(with(w[w$shape == "snv", ], ref != alt))
  )
  cat(glue(
    "OK: {nrow(w)} witnesses for {opt$tx} across {length(unique(w$class))} classes\n"
  ))
  quit(status = 0)
}

# Write the witness VCF. CLASS/TX/SHAPE are INFO fields; IDs provide stable joins.
dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
vc <- file(opt$out, "w")
writeLines(
  c(
    "##fileformat=VCFv4.2",
    glue("##source=duckvep/generate_witnesses.R;tx={opt$tx}"),
    "##INFO=<ID=CLASS,Number=1,Type=String,Description=\"coarse genomic-geometry class\">",
    "##INFO=<ID=TX,Number=1,Type=String,Description=\"target transcript id\">",
    "##INFO=<ID=SHAPE,Number=1,Type=String,Description=\"allele shape (snv/del1/ins1/ins3/mnv2/del2/del3/delins2to1/delins1to2)\">",
    glue("##contig=<ID={chrom}>"),
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  ),
  vc
)
writeLines(
  glue_data(
    w,
    "{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t.\t.\tCLASS={class};TX={opt$tx};SHAPE={shape}"
  ),
  vc
)
close(vc)
cat(
  glue(
    "wrote {nrow(w)} witnesses ({length(unique(w$class))} classes, {length(unique(w$shape))} shapes) -> {opt$out}"
  ),
  "\n",
  sep = ""
)
