suppressMessages(suppressWarnings({
    library(dada2)
}))

# NULL-coalescing helper — defined early so it's available in all functions below
`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)

# Helper to extract a single named argument
get_arg <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0) args[idx + 1] else default
}

# Helper to extract a multi-value named argument (collects until next -- flag)
rget_args_multi <- function(args, flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(NULL)
  start <- idx + 1
  end <- start
  while (end <= length(args) && !startsWith(args[end], "--")) {
    end <- end + 1
  }
  args[start:(end - 1)]
}

input_fwd_R1 <- rget_args_multi(args, "--input_fwd_R1")
input_fwd_R2 <- rget_args_multi(args, "--input_fwd_R2")
input_rev_R1 <- rget_args_multi(args, "--input_rev_R1")
input_rev_R2 <- rget_args_multi(args, "--input_rev_R2")
nproc        <- as.integer(get_arg(args, "--nproc",        "1"))
truncQ       <- as.integer(get_arg(args, "--truncQ",       "2"))
truncLen_R1  <- as.integer(get_arg(args, "--truncLen_R1",  "250"))
truncLen_R2  <- as.integer(get_arg(args, "--truncLen_R2",  "250"))
maxEE_R1     <- as.double(get_arg(args,  "--maxEE_R1",     "2"))
maxEE_R2     <- as.double(get_arg(args,  "--maxEE_R2",     "5"))
maxMismatch  <- as.integer(get_arg(args, "--maxMismatch",  "0"))
minOverlap   <- as.integer(get_arg(args, "--minOverlap",   "12"))

# Minimum gzip file size below which a file is treated as an empty stub.
# Valid gzip streams are >= 20 bytes; the Nextflow process writes 8-byte stubs.
MIN_GZ_BYTES <- 20L

is_nonempty_gz <- function(paths) {
  file.exists(paths) & (file.info(paths)$size >= MIN_GZ_BYTES)
}

getN <- function(x) sum(getUniques(x))

# Strip orientation + R-file suffix from a basename to recover the sample ID.
# Handles: sample_fwd_R1.fastq.gz  and  sample_rev_R1.fastq.gz
basename_to_sampleid <- function(paths, tag) {
  sub(paste0("_", tag, "_R[12]\\.fastq\\.gz$"), "", basename(paths))
}

# ── DADA2 pipeline for one orientation ────────────────────────────────────────
# Returns a list with $seqtab, $out, $dadaFs, $dadaRs, $mergers, $sample_names
# — or NULL if no reads exist for this orientation at all.
run_one_orientation <- function(fnFs, fnRs, tag) {

  if (is.null(fnFs) || length(fnFs) == 0) {
    message("[", tag, "] No input files provided — skipping orientation.")
    return(NULL)
  }

  sample_names <- basename_to_sampleid(fnFs, tag)

  # Drop samples whose input files are empty stubs
  keep_input <- is_nonempty_gz(fnFs) & is_nonempty_gz(fnRs)
  if (!any(keep_input)) {
    message("[", tag, "] No non-empty input files — skipping orientation.")
    return(NULL)
  }
  fnFs         <- fnFs[keep_input]
  fnRs         <- fnRs[keep_input]
  sample_names <- sample_names[keep_input]

  filtFs <- file.path("filtered", paste0(sample_names, "_", tag, "_F_filt.fastq.gz"))
  filtRs <- file.path("filtered", paste0(sample_names, "_", tag, "_R_filt.fastq.gz"))
  names(filtFs) <- sample_names
  names(filtRs) <- sample_names

  out <- filterAndTrim(
    fnFs, filtFs, fnRs, filtRs,
    truncLen = c(truncLen_R1, truncLen_R2),
    maxN = 0, maxEE = c(maxEE_R1, maxEE_R2),
    truncQ = truncQ, rm.phix = TRUE,
    compress = TRUE, multithread = nproc
  )

  # Drop samples that produced zero reads after filtering
  keep_filt <- is_nonempty_gz(filtFs) & is_nonempty_gz(filtRs)
  filtFs_ok <- filtFs[keep_filt]
  filtRs_ok <- filtRs[keep_filt]

  if (length(filtFs_ok) == 0) {
    message("[", tag, "] All samples filtered to zero reads — skipping denoising.")
    return(list(seqtab = NULL, out = out,
                dadaFs = NULL, dadaRs = NULL, mergers = NULL,
                sample_names = sample_names))
  }

  errF <- learnErrors(filtFs_ok, multithread = nproc, nbases = 1e9)
  errR <- learnErrors(filtRs_ok, multithread = nproc, nbases = 1e9)

  dadaFs  <- dada(filtFs_ok, err = errF, multithread = nproc)
  dadaRs  <- dada(filtRs_ok, err = errR, multithread = nproc)
  mergers <- mergePairs(
    dadaFs, filtFs_ok, dadaRs, filtRs_ok,
    verbose = TRUE, maxMismatch = maxMismatch, minOverlap = minOverlap
  )

  seqtab <- makeSequenceTable(mergers)

  list(seqtab = seqtab, out = out,
       dadaFs = dadaFs, dadaRs = dadaRs, mergers = mergers,
       sample_names = sample_names)
}

# ── Run both orientations ──────────────────────────────────────────────────────

fwd <- run_one_orientation(input_fwd_R1, input_fwd_R2, "fwd")
rev <- run_one_orientation(input_rev_R1, input_rev_R2, "rev")

if (is.null(fwd) && is.null(rev)) {
  stop("No reads in either orientation. Cannot produce ASV table.")
}
if (is.null(fwd$seqtab) && is.null(rev$seqtab)) {
  stop("All reads filtered out in both orientations. Cannot produce ASV table.")
}

# ── Merge seqtabs and remove chimeras ─────────────────────────────────────────
# CUTADAPT_DADA2_ORIENT already swapped R1/R2 for reverse-oriented pairs, so
# both seqtabs contain ASVs in the same (forward) direction.
# mergeSequenceTables(..., repeats="sum") collapses identical ASVs that appear
# in both orientation groups for the same sample.

if (is.null(fwd$seqtab)) {
  seqtab_all <- rev$seqtab
} else if (is.null(rev$seqtab)) {
  seqtab_all <- fwd$seqtab
} else {
  seqtab_all <- mergeSequenceTables(fwd$seqtab, rev$seqtab, repeats = "sum")
}

seqtab.nochim <- removeBimeraDenovo(
  seqtab_all, method = "consensus", multithread = nproc, verbose = FALSE
)

# ── ASV table ─────────────────────────────────────────────────────────────────

final_seqtab <- as.data.frame(t(seqtab.nochim))
final_seqtab <- cbind(SeqID = rownames(final_seqtab), final_seqtab)
write.table(final_seqtab, file = "asv_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ── Read tracking ─────────────────────────────────────────────────────────────

# Build a per-sample tracking data.frame for one orientation result.
# Returns NULL if the orientation was entirely skipped (no input files).
build_track <- function(orient) {
  if (is.null(orient)) return(NULL)

  sample_names <- orient$sample_names
  out_mat      <- orient$out

  # filterAndTrim rownames are input file paths — detect tag from first row
  tag <- if (grepl("_fwd_", rownames(out_mat)[1])) "fwd" else "rev"
  rownames(out_mat) <- basename_to_sampleid(rownames(out_mat), tag)

  # Samples that made it all the way through mergePairs
  surviving <- rownames(orient$seqtab %||% matrix(nrow = 0, ncol = 0))

  # Per-dada-object read counts indexed by sample name
  dada_counts <- function(dada_list, sids) {
    if (is.null(dada_list)) return(setNames(rep(0L, length(sids)), sids))
    setNames(sapply(dada_list, getN), names(dada_list))[sids]
  }

  df <- data.frame(
    row.names = sample_names,
    input     = out_mat[sample_names, "reads.in"],
    filtered  = out_mat[sample_names, "reads.out"],
    denoisedF = ifelse(sample_names %in% surviving,
                       dada_counts(orient$dadaFs, names(orient$dadaFs))[sample_names], 0L),
    denoisedR = ifelse(sample_names %in% surviving,
                       dada_counts(orient$dadaRs, names(orient$dadaRs))[sample_names], 0L),
    merged    = ifelse(sample_names %in% surviving,
                       dada_counts(orient$mergers, names(orient$mergers))[sample_names], 0L),
    stringsAsFactors = FALSE
  )
  df[is.na(df)] <- 0L
  df
}

track_fwd <- build_track(fwd)
track_rev <- build_track(rev)

# Union of all sample names seen in either orientation
all_samples <- union(rownames(track_fwd), rownames(track_rev))

tracking_cols <- c("input", "filtered", "denoisedF", "denoisedR", "merged")

# Returns one row of tracking data for sample sid from df, or a zero row if absent
safe_row <- function(df, sid) {
  if (!is.null(df) && sid %in% rownames(df)) {
    df[sid, tracking_cols]
  } else {
    as.data.frame(setNames(as.list(rep(0L, length(tracking_cols))), tracking_cols))
  }
}

# Combine per-sample tracking: sum fwd + rev for each metric.
# fwd_input / rev_input columns let users see the orientation split per sample.
track_list <- lapply(all_samples, function(sid) {
  f <- safe_row(track_fwd, sid)
  r <- safe_row(track_rev, sid)
  nochim_count <- if (sid %in% rownames(seqtab.nochim)) rowSums(seqtab.nochim)[sid] else 0L
  data.frame(
    SampleID  = sid,
    fwd_input = f$input,
    rev_input = r$input,
    input     = f$input     + r$input,
    filtered  = f$filtered  + r$filtered,
    denoisedF = f$denoisedF + r$denoisedF,
    denoisedR = f$denoisedR + r$denoisedR,
    merged    = f$merged    + r$merged,
    nonchim   = as.integer(nochim_count),
    stringsAsFactors = FALSE
  )
})

final_track <- do.call(rbind, track_list)
write.table(final_track, file = "track_control.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ── FASTA for taxonomic assignment ────────────────────────────────────────────

asv_seqs    <- colnames(seqtab.nochim)
asv_headers <- paste0(">", asv_seqs)
writeLines(c(rbind(asv_headers, asv_seqs)), "ASV_sequences.fasta")
