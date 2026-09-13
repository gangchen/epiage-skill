#!/usr/bin/env Rscript
# Local human methylation IDAT -> masked, prefix-collapsed beta matrix.
# Sourceable for tests; dependencies are loaded only when main() executes.

usage <- function() cat(
"Usage:
  Rscript preprocess_idat.R --input <directory|prefix|channel.idat[.gz]> --output-dir <new-dir>
  Rscript preprocess_idat.R --sample-sheet <samples.csv> --output-dir <new-dir>
  Rscript preprocess_idat.R --input <IDAT directory|prefix> --inspect
  Rscript preprocess_idat.R --check-deps [--cache <directory>]
  Rscript preprocess_idat.R --init-cache [--cache <directory>]

Options:
  --sample-sheet CSV with sample_id,idat_prefix; relative paths resolve beside CSV.
  --cache       ExperimentHub cache directory (use the same path at setup and run).
  --detection-p Detection p-value cutoff (default 0.05; mask p > cutoff).
  --min-beads   Minimum bead count (default 1, matching openSesame).
  --save-probe-qc  Also write per-sample probe beta/mask/p-value CSV.gz files.

Runs QCDPB with platform autodetection (HM450, EPIC, EPICv2, MSA when supported
by installed SeSAMe/data). One platform per run. No package install or cache
download is performed during processing. --init-cache explicitly downloads
public SeSAMe reference data; it does not read or upload sample data.
Outputs: betas.csv, qc.csv, samples.csv, run_info.txt, session_info.txt.
", sep = "")

parse_args <- function(args) {
    out <- list(input = NULL, sample_sheet = NULL, output_dir = NULL,
                cache = NULL, detection_p = 0.05, min_beads = 1L,
                save_probe_qc = FALSE, check_deps = FALSE, init_cache = FALSE,
                inspect = FALSE, help = FALSE)
    flags <- c("--save-probe-qc", "--check-deps", "--init-cache", "--inspect", "--help")
    values <- c("--input", "--sample-sheet", "--output-dir", "--cache",
                "--detection-p", "--min-beads")
    i <- 1L
    while (i <= length(args)) {
        arg <- args[i]
        key <- gsub("-", "_", sub("^--", "", arg))
        if (arg %in% flags) {
            out[[key]] <- TRUE
        } else if (arg %in% values && i < length(args)) {
            i <- i + 1L
            out[[key]] <- args[i]
        } else stop("Unknown option or missing value: ", arg, call. = FALSE)
        i <- i + 1L
    }
    out$detection_p <- suppressWarnings(as.numeric(out$detection_p))
    out$min_beads <- suppressWarnings(as.numeric(out$min_beads))
    if (!is.finite(out$detection_p) || out$detection_p <= 0 || out$detection_p >= 1)
        stop("--detection-p must be between 0 and 1.", call. = FALSE)
    if (!is.finite(out$min_beads) || out$min_beads < 1 || out$min_beads %% 1 != 0)
        stop("--min-beads must be a positive integer.", call. = FALSE)
    out
}

idat_prefix <- function(path) {
    path <- path.expand(path)
    path <- sub("_(Grn|Red)\\.idat(\\.gz)?$", "", path)
    file.path(normalizePath(dirname(path), mustWork = TRUE), basename(path))
}

resolve_pair <- function(prefix) {
    vapply(c("Grn", "Red"), function(channel) {
        candidates <- paste0(prefix, "_", channel, c(".idat", ".idat.gz"))
        found <- candidates[file.exists(candidates) & !dir.exists(candidates)]
        if (length(found) != 1L)
            stop("Expected exactly one ", channel, " IDAT for ", prefix,
                 "; found ", length(found), ". Check missing/duplicate compressed files.", call. = FALSE)
        if (is.na(file.info(found)$size) || file.info(found)$size == 0)
            stop("Empty IDAT: ", found, call. = FALSE)
        normalizePath(found, mustWork = TRUE)
    }, character(1))
}

discover_samples <- function(input = NULL, sample_sheet = NULL) {
    if (is.null(input) == is.null(sample_sheet))
        stop("Provide exactly one of --input or --sample-sheet.", call. = FALSE)
    if (!is.null(sample_sheet)) {
        sample_sheet <- normalizePath(path.expand(sample_sheet), mustWork = TRUE)
        samples <- read.csv(sample_sheet, colClasses = "character", check.names = FALSE)
        if (!all(c("sample_id", "idat_prefix") %in% names(samples)) || anyDuplicated(names(samples)))
            stop("Sample sheet requires unique columns sample_id,idat_prefix.", call. = FALSE)
        samples <- samples[, c("sample_id", "idat_prefix"), drop = FALSE]
        if (anyNA(samples) || any(trimws(samples$idat_prefix) == ""))
            stop("Sample sheet contains missing sample IDs or prefixes.", call. = FALSE)
        samples$idat_prefix <- vapply(samples$idat_prefix, function(p) {
            p <- path.expand(p)
            if (!grepl("^(/|[A-Za-z]:[/\\\\])", p)) p <- file.path(dirname(sample_sheet), p)
            idat_prefix(p)
        }, character(1))
    } else {
        input <- path.expand(input)
        if (dir.exists(input)) {
            files <- list.files(input, pattern = "\\.idat(\\.gz)?$",
                                full.names = TRUE, recursive = TRUE)
            if (any(!grepl("_(Grn|Red)\\.idat(\\.gz)?$", files)))
                stop("IDAT filenames must end in _Grn.idat or _Red.idat (optionally .gz).", call. = FALSE)
            prefixes <- sort(unique(vapply(files, idat_prefix, character(1))))
        } else prefixes <- idat_prefix(input)
        samples <- data.frame(sample_id = basename(prefixes), idat_prefix = prefixes,
                              stringsAsFactors = FALSE)
    }
    if (nrow(samples) == 0) stop("No IDAT files found.", call. = FALSE)
    if (any(trimws(samples$sample_id) == "") || anyDuplicated(samples$sample_id))
        stop("Sample IDs must be nonempty and unique; use --sample-sheet to name duplicate basenames.", call. = FALSE)
    if (any(samples$sample_id %in% c("CpG", "Beta_value")))
        stop("Sample IDs CpG and Beta_value are reserved; choose another sample_id.", call. = FALSE)
    if (anyDuplicated(samples$idat_prefix))
        stop("Same IDAT prefix listed more than once.", call. = FALSE)
    pairs <- lapply(samples$idat_prefix, resolve_pair)
    samples$green_idat <- vapply(pairs, `[[`, character(1), "Grn")
    samples$red_idat <- vapply(pairs, `[[`, character(1), "Red")
    rownames(samples) <- NULL
    samples
}

load_sesame <- function(cache = NULL, offline = TRUE) {
    bundled <- FALSE
    if (is.null(cache) && offline) {
        script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
        archive <- file.path(dirname(script), "..", "data", "sesame-reference-cache.tar.gz")
        if (is.na(script) || !file.exists(archive))
            stop("Bundled SeSAMe references are missing. Supply --cache with an offline cache, ",
                 "or ask the user whether to download the missing public references; see references/idat-sesame.md.")
        cache <- tempfile("epiage-reference-cache-")
        dir.create(cache)
        utils::untar(archive, exdir = cache)
        bundled <- TRUE
    }
    if (!is.null(cache)) {
        cache <- path.expand(cache)
        dir.create(cache, recursive = TRUE, showWarnings = FALSE)
        cache <- normalizePath(cache, mustWork = TRUE)
        Sys.setenv(EXPERIMENT_HUB_CACHE = cache)
    }
    failures <- vapply(c("sesame", "sesameData"), function(p) {
        tryCatch({ loadNamespace(p); "" }, error = function(e) conditionMessage(e))
    }, character(1))
    if (any(nzchar(failures))) stop("Missing or incompatible R dependencies:\n",
        paste(paste0("  - ", names(failures)[nzchar(failures)], ": ", failures[nzchar(failures)]), collapse = "\n"),
        ". Ask the user whether to download and install the missing dependencies; ",
        "see references/idat-sesame.md.", call. = FALSE)
    if (bundled) {
        versions <- read.csv(file.path(cache, "bundle_versions.csv"), colClasses = "character")
        for (i in seq_len(nrow(versions))) {
            if (as.character(packageVersion(versions$package[i])) != versions$version[i])
                stop("Bundled annotations were tested with ", versions$package[i], " ", versions$version[i],
                     ". Install the matching package versions, or prepare a matching --cache for your ",
                     "installed version after asking the user whether downloads are allowed.")
        }
    }
    needed <- c("readIDATpair", "prepSesame", "pOOBAH", "noob", "getBetas", "sesameQC_calcStats")
    if (!all(needed %in% getNamespaceExports("sesame")) ||
        !all(c("collapseToPfx", "collapseMethod") %in% names(formals(sesame::getBetas))) ||
        !"min_beads" %in% names(formals(sesame::readIDATpair)))
        stop("Installed SeSAMe lacks required APIs; upgrade SeSAMe and sesameData together.", call. = FALSE)
    # No remote fallback, metadata refresh, or implicit resource download at runtime.
    options(SESAMEDATA_USE_ALT = FALSE)
    ExperimentHub::setExperimentHubOption("LOCAL", offline)
    ExperimentHub::setExperimentHubOption("ASK", FALSE)
    AnnotationHub::setAnnotationHubOption("LOCAL", offline)
    AnnotationHub::setAnnotationHubOption("ASK", FALSE)
    invisible(NULL)
}

required_data_titles <- function() {
    titles <- sesameData::sesameDataList()$Title
    basic <- c("idatSignature", paste0(c("HM450", "EPIC", "EPICv2", "MSA"), ".address"))
    masks <- titles[grepl("^KYCG\\.(HM450|EPIC|EPICv2|MSA)\\.Mask\\.", titles)]
    if (!all(basic %in% titles)) stop("SeSAMe data package lacks supported human array resources; upgrade it.")
    unique(c(basic, masks))
}

check_resources <- function(titles) {
    catalog <- sesameData::sesameDataList()
    failed <- character()
    for (title in titles) {
        error <- tryCatch({ invisible(sesameData::sesameDataGet(title)); NULL }, error = identity)
        if (!is.null(error)) {
            id <- catalog$EHID[match(title, catalog$Title)]
            failed <- c(failed, paste0("  - ", title, " (", id, ")"))
        }
    }
    if (length(failed)) stop("Missing or unreadable local SeSAMe reference resources:\n",
        paste(failed, collapse = "\n"), "\nNo download was attempted. Ask the user whether downloads ",
        "are allowed; then run --init-cache --cache DIR, or import a prepared offline cache. ",
        "The bundled resource inventory and sizes are in data/sesame-resources.json.", call. = FALSE)
    invisible(TRUE)
}

read_idat_identity <- function(path, fields) {
    # SeSAMe 1.24 reads signal fields only, leaving Barcode/ChipType NULL.
    # Read length-prefixed identity strings at the offsets it already validated.
    con <- if (grepl("\\.gz$", path)) gzfile(path, "rb") else file(path, "rb")
    on.exit(close(con))
    read_string <- function(field) {
        if (!field %in% rownames(fields)) return(NA_character_)
        seek(con, fields[field, "byteOffset"], origin = "start")
        n <- 0
        for (i in 0:4) {
            byte <- readBin(con, "integer", size = 1, signed = FALSE, n = 1)
            if (length(byte) != 1L) stop("Truncated IDAT metadata: ", path)
            n <- n + (byte %% 128) * 2^(7 * i)
            if (n > 65536) stop("Invalid IDAT identity string length: ", path)
            if (byte < 128) break
            if (i == 4) stop("Invalid IDAT metadata length encoding: ", path)
        }
        if (!n) return(NA_character_)
        value <- readBin(con, "raw", n = n)
        if (length(value) != n) stop("Truncated IDAT identity string: ", path)
        rawToChar(value)
    }
    c(barcode = read_string("Barcode"), chip_type = read_string("ChipType"), stripe = read_string("MostlyA"))
}

inspect_pair <- function(green, red, sample_id) {
    # SeSAMe's IDAT reader validates the binary magic/version. The signature
    # detector matches actual Illumina addresses, not filenames or probe counts.
    # These two internal helpers are API-checked and covered by the smoke test.
    ns <- asNamespace("sesame")
    if (!all(vapply(c("readIDAT", "inferPlatformFromTango"), exists,
                    logical(1), envir = ns, inherits = FALSE)))
        stop("Installed SeSAMe lacks the tested IDAT inspection helpers.")
    read <- get("readIDAT", ns)
    infer <- get("inferPlatformFromTango", ns)
    g <- read(green)
    r <- read(red)
    if (!identical(rownames(g$Quants), rownames(r$Quants)))
        stop("Red/green probe addresses differ or are reordered: ", sample_id)
    g_identity <- read_idat_identity(green, g$fields)
    r_identity <- read_idat_identity(red, r$fields)
    if (!identical(g_identity, r_identity))
        stop("Red/green barcode, chip type or array position mismatch: ", sample_id)
    platform <- infer(r)
    if (length(platform) != 1L || !platform %in% c("HM450", "EPIC", "EPICv2", "MSA"))
        stop("IDAT is not a recognized supported human methylation array: ", sample_id)
    data.frame(sample_id = sample_id, platform = platform,
        chip_type = unname(r_identity["chip_type"]), barcode = unname(r_identity["barcode"]),
        array_position = unname(r_identity["stripe"]),
        pair_identity_available = !anyNA(r_identity[c("barcode", "stripe")]),
        n_addresses = nrow(r$Quants), stringsAsFactors = FALSE)
}

inspect_samples <- function(samples) {
    # The only annotation required for identification is the small signature set.
    check_resources("idatSignature")
    info <- do.call(rbind, lapply(seq_len(nrow(samples)), function(i)
        inspect_pair(samples$green_idat[i], samples$red_idat[i], samples$sample_id[i])))
    if (length(unique(info$platform)) != 1L)
        stop("Mixed platforms detected; process each platform in a separate run.")
    info
}

check_annotations <- function(platform) {
    titles <- required_data_titles()
    check_resources(c(paste0(platform, ".address"),
                      titles[startsWith(titles, paste0("KYCG.", platform, ".Mask."))]))
    masks <- sesame::getMask(platform)
    if (!length(masks)) warning("No recommended design mask in this SeSAMe/data version for ",
                               platform, "; review this limitation before interpreting clocks.")
    length(masks) > 0L
}

write_csv <- function(x, path) {
    con <- if (grepl("\\.gz$", path)) gzfile(path, "wt") else file(path, "wt")
    on.exit(close(con))
    write.csv(x, con, row.names = FALSE, na = "NA")
}

beta_matrix <- function(betas, sample_ids) {
    # Explicit name-based union prevents positional misalignment across samples.
    ids <- sort(unique(unlist(lapply(betas, names), use.names = FALSE)))
    result <- matrix(NA_real_, nrow = length(ids), ncol = length(betas),
                     dimnames = list(ids, sample_ids))
    for (i in seq_along(betas)) {
        if (anyDuplicated(names(betas[[i]]))) stop("Duplicate beta IDs after SeSAMe collapse.")
        result[match(names(betas[[i]]), ids), i] <- betas[[i]]
    }
    result
}

process_sample <- function(prefix, sample_id, opts, probe_qc_file = NULL) {
    raw <- sesame::readIDATpair(prefix, min_beads = opts$min_beads)
    platform <- attr(raw, "platform")
    if (length(platform) != 1L || !platform %in% c("HM450", "EPIC", "EPICv2", "MSA"))
        stop("Unsupported or unidentified human platform: ", paste(platform, collapse = ", "),
             ". Upgrade SeSAMe/data if this is a newer human array; do not force another manifest.")
    # Avoid the default 'betas' QC routine: it performs another preprocessing pass.
    raw_qc <- as.data.frame(sesame::sesameQC_calcStats(raw,
        c("intensity", "numProbes", "channel", "dyeBias")))
    names(raw_qc) <- paste0("raw_", names(raw_qc))
    sdf <- sesame::prepSesame(raw, "QCD")
    pval <- sesame::pOOBAH(sdf, return.pval = TRUE)
    pval <- pval[match(sdf$Probe_ID, names(pval))]
    # P before B: noob changes the out-of-band background used for detection.
    sdf <- sesame::pOOBAH(sdf, pval.threshold = opts$detection_p)
    # Older releases may return NA detection p-values. Never treat these as observed.
    sdf$mask <- is.na(sdf$mask) | sdf$mask | !is.finite(pval)
    sdf <- sesame::noob(sdf)
    full <- sesame::getBetas(sdf, mask = TRUE, collapseToPfx = FALSE)
    betas <- sesame::getBetas(sdf, mask = TRUE, collapseToPfx = TRUE, collapseMethod = "mean")
    betas[!is.finite(betas)] <- NA_real_
    # Keep cg and non-CpG ch features used by existing models; drop SNP/control rows.
    betas <- betas[grepl("^(cg[0-9]+$|ch[.])", names(betas))]
    if (!length(betas) || !any(is.finite(betas))) stop("No usable methylation betas after QC: ", sample_id)
    if (any(betas < 0 | betas > 1, na.rm = TRUE)) stop("Out-of-range beta values: ", sample_id)
    fraction <- mean(is.na(betas))
    qc <- data.frame(sample_id = sample_id, platform = platform, prep = "QCDPB",
        design_mask_available = opts$design_mask_available,
        detection_p_cutoff = opts$detection_p, min_beads = opts$min_beads,
        n_probes = nrow(sdf), n_masked_probes = sum(sdf$mask),
        fraction_masked_probes = mean(sdf$mask),
        n_detection_failed = sum(!is.finite(pval) | pval > opts$detection_p),
        n_beta_loci = length(betas), n_observed_loci = sum(!is.na(betas)),
        n_missing_loci = sum(is.na(betas)), fraction_missing_loci = fraction,
        beta_median = median(betas, na.rm = TRUE),
        qc_status = if (!opts$design_mask_available) "review_design_mask_unavailable" else
            if (fraction >= 0.20) "review_high_missingness" else "review_metrics",
        stringsAsFactors = FALSE)
    qc <- cbind(qc, raw_qc)
    if (!is.null(probe_qc_file)) write_csv(data.frame(
        Probe_ID = sdf$Probe_ID, detection_p = unname(pval),
        masked = sdf$mask, beta = unname(full[sdf$Probe_ID])), probe_qc_file)
    list(betas = betas, qc = qc, platform = platform)
}

run_preprocessing <- function(samples, opts, inspection) {
    target <- path.expand(opts$output_dir)
    if (file.exists(target) || dir.exists(target))
        stop("Output path already exists; choose a new --output-dir: ", target, call. = FALSE)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    parent <- normalizePath(dirname(target), mustWork = TRUE)
    target <- file.path(parent, basename(target))
    staging <- tempfile(".sesame-", tmpdir = parent)
    dir.create(staging)
    on.exit(unlink(staging, recursive = TRUE), add = TRUE)
    all_betas <- vector("list", nrow(samples))
    all_qc <- vector("list", nrow(samples))
    platform <- NULL
    for (i in seq_len(nrow(samples))) {
        message(sprintf("[%d/%d] %s", i, nrow(samples), samples$sample_id[i]))
        probe_file <- if (opts$save_probe_qc) file.path(staging, sprintf("probe_qc_%03d.csv.gz", i)) else NULL
        result <- process_sample(samples$idat_prefix[i], samples$sample_id[i], opts, probe_file)
        if (!is.null(platform) && platform != result$platform)
            stop("Mixed platforms detected; process each platform in a separate run.")
        platform <- result$platform
        all_betas[[i]] <- result$betas
        all_qc[[i]] <- result$qc
    }
    matrix <- beta_matrix(all_betas, samples$sample_id)
    write_csv(inspection, file.path(staging, "input_inspection.csv"))
    write_csv(data.frame(CpG = rownames(matrix), matrix, check.names = FALSE), file.path(staging, "betas.csv"))
    write_csv(do.call(rbind, all_qc), file.path(staging, "qc.csv"))
    samples$green_md5 <- unname(tools::md5sum(samples$green_idat))
    samples$red_md5 <- unname(tools::md5sum(samples$red_idat))
    if (opts$save_probe_qc) samples$probe_qc_file <- sprintf("probe_qc_%03d.csv.gz", seq_len(nrow(samples)))
    write_csv(samples, file.path(staging, "samples.csv"))
    writeLines(c(paste("completed_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
        paste("sesame:", packageVersion("sesame")), paste("sesameData:", packageVersion("sesameData")),
        paste("ExperimentHub:", packageVersion("ExperimentHub")),
        paste("cache:", ExperimentHub::getExperimentHubOption("CACHE")),
        paste("platform:", platform), "prep: QCDPB", paste("detection_p_cutoff:", opts$detection_p),
        paste("min_beads:", opts$min_beads), "collapse: SeSAMe getBetas mean after masking",
        "QC missingness >=20%: review heuristic, not a validated sample pass/fail criterion",
        "Missing values remain NA; no imputation in IDAT preprocessing."), file.path(staging, "run_info.txt"))
    writeLines(capture.output(sessionInfo()), file.path(staging, "session_info.txt"))
    if (file.exists(target) || dir.exists(target)) stop("Output appeared during processing; refusing to overwrite: ", target)
    if (!file.rename(staging, target)) stop("Could not finalize output directory: ", target)
    message("Wrote ", target, ". Review qc.csv before calculating clocks.")
    invisible(target)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
    opts <- parse_args(args)
    if (opts$help) { usage(); return(invisible(NULL)) }
    if (opts$init_cache || opts$check_deps) {
        if (!is.null(opts$input) || !is.null(opts$sample_sheet) || !is.null(opts$output_dir))
            stop("Run --init-cache/--check-deps separately from sample processing.")
        if (opts$init_cache && is.null(opts$cache))
            stop("--init-cache requires an explicit --cache directory. Use only after download consent.")
        load_sesame(opts$cache, offline = !opts$init_cache)
        if (opts$init_cache) {
            sesameData::sesameDataCache(required_data_titles())
            # Verify actual local retrieval; some sesameData releases swallow cache errors.
            ExperimentHub::setExperimentHubOption("LOCAL", TRUE)
            check_resources(required_data_titles())
        }
        if (opts$check_deps) check_resources(required_data_titles())
        cat("sesame", as.character(packageVersion("sesame")), "| sesameData",
            as.character(packageVersion("sesameData")), "\n")
        return(invisible(NULL))
    }
    if (is.null(opts$output_dir) && !opts$inspect) stop("--output-dir is required. Use --help for examples.")
    samples <- discover_samples(opts$input, opts$sample_sheet)
    load_sesame(opts$cache)
    inspection <- inspect_samples(samples)
    print(inspection[, c("sample_id", "platform", "chip_type", "n_addresses")], row.names = FALSE)
    if (opts$inspect) return(invisible(inspection))
    opts$design_mask_available <- check_annotations(inspection$platform[1])
    run_preprocessing(samples, opts, inspection)
}

if (sys.nframe() == 0L) tryCatch(main(), error = function(e) {
    message("ERROR: ", conditionMessage(e)); quit(status = 1L)
})
