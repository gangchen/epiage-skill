# Pure IO/assembly tests require only base R; no patient files or annotation download.
source("epigenetic-clocks/scripts/preprocess_idat.R")

expect_error <- function(expr, pattern) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
temporary <- tempfile("idat-tests-")
dir.create(temporary)
run_tests <- function() {
    on.exit(unlink(temporary, recursive = TRUE))
    pair <- function(prefix) {
        dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
        for (channel in c("Grn", "Red")) writeBin(charToRaw("IDAT"), paste0(prefix, "_", channel, ".idat"))
    }
    prefix <- file.path(temporary, "space and 中文", "sample_A")
    pair(prefix)
    samples <- discover_samples(temporary)
    stopifnot(nrow(samples) == 1L, samples$sample_id == "sample_A")
    stopifnot(identical(discover_samples(paste0(prefix, "_Red.idat")), samples))
    # Mixed compression is supported; duplicate representations are rejected.
    file.copy(paste0(prefix, "_Red.idat"), paste0(prefix, "_Red.idat.gz"))
    expect_error(discover_samples(prefix), "exactly one Red")
    unlink(paste0(prefix, "_Red.idat"))
    stopifnot(grepl("gz$", discover_samples(prefix)$red_idat))
    unlink(paste0(prefix, "_Red.idat.gz"))
    expect_error(discover_samples(temporary), "exactly one Red")
    pair(prefix)
    duplicate <- file.path(temporary, "second", "sample_A")
    pair(duplicate)
    expect_error(discover_samples(temporary), "unique")
    sheet <- file.path(temporary, "samples.csv")
    write.csv(data.frame(sample_id = c("Person 1", "Person 2"),
        idat_prefix = c("space and 中文/sample_A", "second/sample_A")), sheet, row.names = FALSE)
    selected <- discover_samples(sample_sheet = sheet)
    stopifnot(identical(selected$sample_id, c("Person 1", "Person 2")), nrow(selected) == 2L)
    expect_error(discover_samples(temporary, sheet), "exactly one")
    expect_error(parse_args(c("--min-beads", "0")), "positive integer")
    expect_error(parse_args(c("--detection-p", "1")), "between 0 and 1")
    expect_error(parse_args(c("--input")), "missing value")
    # IDAT identity strings use variable-length lengths; keep metadata checks
    # independent of SeSAMe's signal-only reader and validate gzip too.
    identity_file <- file.path(temporary, "identity.bin")
    fields <- cbind(byteOffset = c(0, 4, 8))
    rownames(fields) <- c("Barcode", "ChipType", "MostlyA")
    payload <- c(as.raw(3), charToRaw("123"), as.raw(3), charToRaw("abc"), as.raw(3), charToRaw("R01"))
    writeBin(payload, identity_file)
    expected <- c(barcode = "123", chip_type = "abc", stripe = "R01")
    stopifnot(identical(read_idat_identity(identity_file, fields), expected))
    compressed <- gzfile(paste0(identity_file, ".gz"), "wb")
    writeBin(payload, compressed)
    close(compressed)
    stopifnot(identical(read_idat_identity(paste0(identity_file, ".gz"), fields), expected))
    m <- beta_matrix(list(c(cg2 = NA_real_, cg1 = 0.2), c(cg3 = 0.9, cg1 = 0.8)), c("A", "B"))
    stopifnot(identical(rownames(m), c("cg1", "cg2", "cg3")), m["cg1", "A"] == 0.2,
              m["cg1", "B"] == 0.8, is.na(m["cg2", "B"]), is.na(m["cg3", "A"]))
    cat("IDAT discovery, pairing, sample-sheet, argument and matrix tests passed.\n")
}
run_tests()
