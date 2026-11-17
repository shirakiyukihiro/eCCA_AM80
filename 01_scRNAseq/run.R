# =============================================================================
# Step runner for 01_mouse_biliary_fibroblast_scRNAseq.R
#
# The analysis lives in one script so that it can be published as one file, but
# during development you almost never want to run all of it. This runner lets
# you name the steps you want:
#
#   source("run.R")
#   steps()              # list the steps and which of them have been run
#   run_steps(9)         # config + helpers, then step 9
#   run_steps(9:16)      # ... then steps 9 to 16
#   run_steps(9, reload = FALSE)   # keep the in-memory objects (see below)
#
# Config and helper functions (steps 0 and 1) are ALWAYS run first, so a change
# to CONFIG or to a helper never has to be hunted for and executed by hand.
#
# Objects a step needs are RE-READ FROM THE CHECKPOINT every time, because some
# steps modify their own input: step 9 drops contaminated cells from sc_fib, so
# running it twice on the in-memory object re-clusters an object that has
# already been filtered, and the second result is not reproducible from the
# script. Re-reading costs a few seconds and removes that whole class of error.
# reload = FALSE keeps whatever is in the session, which is only safe for steps
# that do not modify their inputs.
#
# Nothing here is part of the analysis. The published script does not need it.
# =============================================================================

SCRIPT <- "01_mouse_biliary_fibroblast_scRNAseq.R"

# Objects each step needs, and the checkpoint each is restored from. Only used
# when the object is missing from the session (or reload = TRUE).
# A step must never read the checkpoint it writes. Step 9 removes contaminated
# cells, so if it read and wrote the same file, every re-run would start from its
# own previous output and produce a different answer. It reads
# 03_fibroblasts.rds (written by step 8) and writes 03c_fibroblast_subsets.rds.
#
# ONE EXCEPTION, added 2026-09-29. Step 11 now reads AND writes
# 03c_fibroblast_subsets.rds, because the programme scores have to reach
# 02_figure_panels_Fig3_S1.R and 03c is the only file it loads. This is safe because the
# step is idempotent: it changes no counts, no embedding and no cell, it only
# writes the programme_* metadata columns, and it purges any programme_* columns
# it finds before scoring, so re-running it on its own output reproduces exactly
# the same columns. Do not extend this exception to any step that changes the
# cells, the matrix or the clustering.
#
# There is no automatic enforcement of the rule above; a sentence here claiming
# otherwise was wrong and has been removed (no such function ever existed).
STEP_INPUTS <- list(
  "6"  = c(sc     = "01_allcells_clustered.rds"),
  "7"  = c(sc     = "02_allcells_annotated.rds"),
  "8"  = c(sc     = "02_allcells_annotated.rds"),
  "9"  = c(sc_fib = "03_fibroblasts.rds"),
  "10" = c(sc_fib = "03c_fibroblast_subsets.rds"),
  "11" = c(sc_fib = "03c_fibroblast_subsets.rds"),
  "12" = c(sc_fib = "03c_fibroblast_subsets.rds",
           sc_ds  = "03b_fibroblasts_downsampled.rds"),
  "13" = c(sc_fib = "03c_fibroblast_subsets.rds"),
  # sc_ds added 2026-09-29: step 14 now builds its ranked lists one versus
  # rest from the DOWNSAMPLED object, so the depth asymmetry is removed at
  # source rather than by excluding gene families afterwards.
  "14" = c(sc_fib = "03c_fibroblast_subsets.rds",
           sc_ds  = "03b_fibroblasts_downsampled.rds"),
  "15" = c(sc_fib = "03c_fibroblast_subsets.rds")
)

# --- parse the script into steps ----------------------------------------------
.parse_steps <- function(path = SCRIPT) {
  src <- readLines(path, warn = FALSE)
  at  <- grep("^#@ STEP [0-9]+$", src)
  if (!length(at)) stop("no '#@ STEP n' markers found in ", path)
  num <- as.integer(sub("^#@ STEP ", "", src[at]))
  end <- c(at[-1] - 1L, length(src))
  title <- vapply(at, function(i) sub("^# [0-9]+\\. ", "", src[i + 2L]),
                  character(1))
  list(src = src,
       tbl = data.frame(step = num, title = title, from = at, to = end,
                        n_lines = end - at + 1L, stringsAsFactors = FALSE))
}

# --- list the steps ------------------------------------------------------------
steps <- function(path = SCRIPT) {
  p <- .parse_steps(path)
  ck <- if (exists("obj", mode = "function")) {
    vapply(p$tbl$step, function(s) {
      f <- unname(STEP_INPUTS[[as.character(s + 1L)]])
      if (is.null(f)) "" else paste(ifelse(file.exists(obj(f)), "ok", "MISSING"),
                                    collapse = ",")
    }, character(1))
  } else rep("", nrow(p$tbl))
  out <- p$tbl[, c("step", "title", "from", "to", "n_lines")]
  out$writes_checkpoint <- ck
  print(out, row.names = FALSE)
  invisible(out)
}

# --- run ----------------------------------------------------------------------
run_steps <- function(which_steps, path = SCRIPT, reload = TRUE,
                      with_setup = TRUE, echo = FALSE) {
  p   <- .parse_steps(path)
  tbl <- p$tbl
  ask <- as.integer(which_steps)
  unknown <- setdiff(ask, tbl$step)
  if (length(unknown))
    stop("no such step(s): ", paste(unknown, collapse = ", "),
         ". Run steps() to list them.")

  setup <- if (with_setup) intersect(c(0L, 1L), tbl$step) else integer(0)
  todo  <- sort(unique(c(setup, ask)))

  message("Running step(s): ", paste(todo, collapse = ", "))

  for (s in todo) {
    r <- tbl[tbl$step == s, ]

    # restore the objects this step needs, if they are not already here
    need <- STEP_INPUTS[[as.character(s)]]
    if (!is.null(need)) {
      for (nm in names(need)) {
        if (reload || !exists(nm, envir = globalenv())) {
          f <- obj(need[[nm]])
          if (!file.exists(f))
            stop("step ", s, " needs '", nm, "' and the checkpoint is missing:\n  ",
                 f, "\nRun the earlier step that writes it (see steps()).")
          message("  loading ", nm, " <- ", basename(f))
          assign(nm, readRDS(f), envir = globalenv())
        }
      }
    }

    # mtimes of this step's inputs, to detect a step overwriting what it read
    in_files <- if (is.null(need)) character(0) else unname(vapply(need, obj, ""))
    before   <- file.mtime(in_files)

    t0 <- Sys.time()
    message(sprintf("--- step %d: %s (lines %d-%d)", s, r$title, r$from, r$to))
    code <- p$src[r$from:r$to]
    # "object 'deg_robust' not found" names the symbol and not the cause. The
    # cause is always the same: it is made in a step that was not run, and in a
    # fresh session that is every step except the one you asked for. Finding
    # where it is assigned takes one grep, so do the grep.
    eval_step <- function() eval(parse(text = paste(code, collapse = "\n")),
                                 envir = globalenv())
    tryCatch(eval_step(), error = function(e) {
        m <- regmatches(conditionMessage(e),
                        regexpr("object '[^']+' not found", conditionMessage(e)))
        if (length(m)) {
          sym <- sub("^object '([^']+)' not found$", "\\1", m)
          hit <- grep(paste0("^\\s*", sym, "\\s*(<-|=[^=])"), p$src)
          where <- vapply(hit, function(i) {
            k <- tbl$step[tbl$from <= i & tbl$to >= i]
            if (length(k)) k[1] else NA_integer_
          }, integer(1))
          where <- sort(unique(where[!is.na(where) & where != s]))
          if (length(where))
            message("\n  '", sym, "' is created in step(s) ",
                    paste(where, collapse = ", "),
                    ". Run run_steps(c(", paste(c(where, s), collapse = ", "),
                    ")) in this session.")
        }
      stop(e)
    })
    message(sprintf("--- step %d done in %s", s,
                    format(round(Sys.time() - t0, 1))))

    if (length(in_files)) {
      touched <- in_files[!is.na(before) & file.mtime(in_files) != before]
      if (length(touched))
        warning("step ", s, " overwrote a checkpoint it reads:\n  ",
                paste(basename(touched), collapse = "\n  "),
                "\nRe-running it will start from its own previous output ",
                "instead of from the same input, and the result will drift. ",
                "Give the step its own output file.", call. = FALSE)
    }
  }
  invisible(todo)
}

# --- constants that cross a step boundary ---------------------------------------
# run_steps() guarantees sections 0 and 1 and nothing else, so a constant defined
# in section 6 and read in section 9 makes run_steps(9) fail in a fresh session
# with "object not found" -- an error that names the symbol but not the cause,
# and that only appears after a restart, which is when it is least expected.
#
# The check is deliberately crude: top-level assignments to ALL-CAPS names, and
# whether that name appears in a later section. It cannot see anything else, so a
# clean report is not a guarantee. What it does catch is the exact mistake that
# has now happened once, at the moment run.R is loaded rather than ten minutes
# into a re-run.
check_step_constants <- function(path = SCRIPT) {
  p <- .parse_steps(path)
  tbl <- p$tbl[order(p$tbl$step), ]
  defs <- list()
  for (i in seq_len(nrow(tbl))) {
    body <- p$src[tbl$from[i]:tbl$to[i]]
    nm <- sub("^([A-Z][A-Z0-9_]{2,})[[:space:]]*<-.*$", "\\1",
              grep("^[A-Z][A-Z0-9_]{2,}[[:space:]]*<-", body, value = TRUE))
    for (n in unique(nm))
      if (is.null(defs[[n]])) defs[[n]] <- tbl$step[i]
  }
  bad <- character(0)
  for (n in names(defs)) {
    s0 <- defs[[n]]
    if (s0 <= 1) next                       # always executed, so always present
    later <- tbl$step[tbl$step > s0]
    hit <- vapply(later, function(s) {
      b <- p$src[tbl$from[tbl$step == s]:tbl$to[tbl$step == s]]
      any(grepl(paste0("\\b", n, "\\b"), b))
    }, logical(1))
    if (any(hit))
      bad <- c(bad, sprintf("  %s  defined in step %d, used in step(s) %s",
                            n, s0, paste(later[hit], collapse = ", ")))
  }
  if (length(bad))
    warning("these constants are defined in a step that run_steps() will not ",
            "necessarily run first:\n", paste(bad, collapse = "\n"),
            "\nMove them to section 1, or those steps will fail in a fresh ",
            "session.", call. = FALSE, immediate. = TRUE)
  invisible(bad)
}
check_step_constants()

message("run.R loaded. steps() to list, run_steps(n) to run.")
