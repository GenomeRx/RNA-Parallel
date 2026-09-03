# The rails that keep the other suites honest, and the one that keeps the machine alive.
# Each of these guards a failure mode that reports green while doing nothing.

test_that("the gate list the suite zeroes matches the gates the code reads", {
  # setup-parallel.R zeroes gates by hand. A gate added to R/ and forgotten here would send
  # its companion down the serial fallback in every test, silently, still passing.
  # Source files only. An installed package ships R/ as a compiled .rdb, and reading that as
  # text yields non-empty binary junk that matches nothing, which would fail rather than skip.
  files <- list.files(testthat::test_path("..", ".."), pattern = "\\.R$",
                      recursive = TRUE, full.names = TRUE)
  files <- files[dirname(files) == file.path(testthat::test_path("..", ".."), "R")]
  skip_if(!length(files), "package sources not available, this is an installed check")
  src <- unlist(lapply(files, readLines), use.names = FALSE)

  in_code <- sort(unique(regmatches(
    paste(src, collapse = "\n"),
    gregexpr('combat\\.min\\.[a-z.]+', paste(src, collapse = "\n")))[[1]]))
  zeroed <- sort(names(Filter(function(v) identical(v, 0),
                              options()[grep("^combat\\.min\\.", names(options()))])))
  expect_setequal(in_code, zeroed)
})

test_that("combat_reap kills what a call created and spares what it did not", {
  skip_on_os("windows")
  reap <- rnaparallel:::combat_reap
  kids <- rnaparallel:::combat_children

  # a worker the "user" already had running: it must survive
  outsider <- parallel::mcparallel(Sys.sleep(30))
  spare <- kids()
  expect_true(outsider$pid %in% spare)

  # workers the call creates: these are what reap must take
  made <- replicate(3, parallel::mcparallel(Sys.sleep(30)), simplify = FALSE)
  expect_length(setdiff(kids(), spare), 3L)

  started <- proc.time()[["elapsed"]]
  expect_identical(reap(spare), 3L)
  expect_lt(proc.time()[["elapsed"]] - started, 3)
  expect_setequal(kids(), spare)
  expect_true(outsider$pid %in% kids())        # the caller's own worker is untouched

  expect_gte(reap(integer()), 1L)              # with nothing spared, the outsider is in scope
  expect_length(kids(), 0L)
})

test_that("the shipped report and NEWS belong to the version in DESCRIPTION", {
  root <- testthat::test_path("..", "..")
  desc <- file.path(root, "DESCRIPTION")
  skip_if(!file.exists(desc), "repository-only release files are unavailable")
  version <- read.dcf(desc, fields = "Version")[[1L]]

  expect_identical(readLines(file.path(root, "NEWS.md"), n = 1L, warn = FALSE),
                   paste("# rnaparallel", version))

  # Assert the ARTIFACT, not the script that reads it. The rail this replaces checked that a
  # path string appeared inside tools/doccheck.R, which passes just as happily against a
  # report rendered two versions ago on another machine -- and did, for the whole of 0.4.5.
  # The report prints its own sessionInfo, so the built package version is in the HTML.
  html <- file.path(root, "docs", "index.html")
  skip_if(!file.exists(html), "docs/index.html is not in this tree")
  expect_match(paste(readLines(html, warn = FALSE), collapse = "\n"),
               paste0("rnaparallel_", version), fixed = TRUE,
               info = "docs/index.html was rendered against a different package version")
})

test_that("the roxygen comments in R/ and the committed man pages agree", {
  root <- testthat::test_path("..", "..")
  skip_on_cran()
  skip_if_not_installed("roxygen2")
  skip_if(!dir.exists(file.path(root, "man")), "repository-only release files are unavailable")
  declared <- read.dcf(file.path(root, "DESCRIPTION"),
                       fields = "Config/roxygen2/version")[[1L]]
  skip_if(!is.na(declared) &&
            utils::compareVersion(as.character(utils::packageVersion("roxygen2")),
                                  declared) < 0,
          "installed roxygen2 is older than the version this package declares")

  # In a SUBPROCESS, deliberately. roxygenise() loads the package it documents, so running it
  # here swapped this session's namespace for one rooted in a tempdir; once that tempdir was
  # removed, packageVersion() could no longer find a DESCRIPTION and returned NA, which broke a
  # later test in a different file. A test that rewrites the namespace under the suite is worse
  # than the drift it checks for.
  tmp <- withr::local_tempdir()
  for (d in c("DESCRIPTION", "NAMESPACE", "R", "man")) {
    file.copy(file.path(root, d), tmp, recursive = TRUE)
  }
  script <- tempfile(fileext = ".R")
  writeLines(sprintf('suppressMessages(roxygen2::roxygenise(%s))', deparse(tmp)), script)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
    stdout = TRUE, stderr = TRUE,
    env = c(paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep)))),
    timeout = 120))
  status <- attr(out, "status"); if (is.null(status)) status <- 0L
  expect_identical(status, 0L, info = paste(out, collapse = "\n"))

  for (rd in list.files(file.path(root, "man"), pattern = "\\.Rd$")) {
    expect_identical(readLines(file.path(tmp, "man", rd), warn = FALSE),
                     readLines(file.path(root, "man", rd), warn = FALSE),
                     info = paste(rd, "is out of date; run roxygen2::roxygenise()"))
  }
})

test_that("rp_arrayweights_uniform reads the block-leading rows and nothing else", {
  f <- rnaparallel:::rp_arrayweights_uniform
  w <- matrix(1, 10, 4)
  attr(w, "arrayweights") <- TRUE

  expect_true(f(w, 1:4))                       # uniform
  expect_true(f(NULL, 1:4))                    # nothing to check
  bare <- matrix(1, 10, 4); bare[3, ] <- 2
  expect_true(f(bare, 1:4))                    # no arrayweights attribute, not our problem

  w2 <- w; w2[3L, ] <- 2
  expect_false(f(w2, 1:4))                     # row 3 leads a block, so it must be seen
  expect_true(f(w2, c(1L, 2L, 4L)))            # row 3 leads nothing, so it cannot change a qr
})

test_that("the dupcor tail gate refuses a limma whose post-loop run stops decomposing", {
  skip_if_not_installed("limma")
  # Everything limma runs between its per-gene loop and the pooled tail executes inside every
  # block on that block's own rho. Against 3.62.2 those statements decompose; the gate exists
  # so a future limma that made them depend on which genes share a block is refused instead of
  # returning a quietly different consensus. Without a mutation this test asserts nothing.
  dec <- rnaparallel:::rp_tail_decomposes
  expect_true(dec(limma::duplicateCorrelation, 2L, 4L))

  inject <- function(mut) {
    f <- limma::duplicateCorrelation
    b <- as.list(body(f))
    fori <- max(which(vapply(b, function(s) is.call(s) &&
                               identical(as.character(s[[1]]), "for"), logical(1))))
    body(f) <- as.call(append(b, list(mut), after = fori))
    f
  }
  # each of these makes the result depend on the other genes in the block
  expect_false(dec(inject(quote(rho <- rho - min(rho, na.rm = TRUE))), 2L, 4L))
  expect_false(dec(inject(quote(rho <- rho/max(abs(rho), na.rm = TRUE))), 2L, 4L))
  expect_false(dec(inject(quote(rho <- rho - mean(rho, na.rm = TRUE))), 2L, 4L))
  expect_false(dec(inject(quote(rho <- sort(rho))), 2L, 4L))
  # and one that does not
  expect_true(dec(inject(quote(rho <- rho * 1)), 2L, 4L))
})


test_that("each least-squares branch reads its own size gate", {
  # The Windows merge routed both lmFit branches through one function whose first act was to
  # return combat.min.ls.cells when it was set, so raising the weightless gate silently raised
  # the voom/weighted one from 2e4 to the same value and switched off a split measured at
  # 2.52x-3.39x. The whole suite is blind to it by construction: setup-parallel.R sets every
  # combat.min.* to 0, and 0 is returned from the first line for both branches.
  #
  # The `fork_default` path this test exercises when the option is unset goes through
  # rp_copy_free(), which is intentionally OS-gated: mclapply cannot fork on Windows, so a
  # Windows caller correctly reads Inf (the split never runs) rather than the fork_default
  # value a Unix caller would get. Same reason the neighbouring
  # "the break-even gates follow the backend's copy behaviour, not the OS" test skips on
  # Windows for this exact code path; this test asserts the OS-correct value on each platform
  # instead of skipping outright, since the option-set branch (the actual regression this
  # test exists to catch) is identical on every OS.
  ls_gate <- rnaparallel:::rp_ls_min_cells
  fork_default_ls <- if (rnaparallel:::rp_copy_free("mclapply")) 6e6 else Inf
  fork_default_cells <- if (rnaparallel:::rp_copy_free("mclapply")) 2e4 else Inf
  withr::with_options(list(combat.min.ls.cells = 6e7, combat.min.cells = NULL), {
    expect_identical(ls_gate("combat.min.ls.cells", 6e6, "mclapply"), 6e7)
    expect_identical(ls_gate("combat.min.cells", 2e4, "mclapply"), fork_default_cells)  # NOT 6e7
  })
  withr::with_options(list(combat.min.ls.cells = NULL, combat.min.cells = 5e5), {
    expect_identical(ls_gate("combat.min.cells", 2e4, "mclapply"), 5e5)
    expect_identical(ls_gate("combat.min.ls.cells", 6e6, "mclapply"), fork_default_ls)
  })
})

test_that("the break-even gates follow the backend's copy behaviour, not the OS", {
  # These thresholds were tuned where a worker INHERITS the matrix. Keyed to the OS, a macOS
  # caller on foreach got the inherited gate and lmFit measured 0.24x against the original.
  # The question is the COPY, not fork(): foreach builds a FORK cluster on Unix and is still
  # slow, because doParallel's cluster form serialises every task. So foreach is asked rather
  # than assumed -- doParallelMC drives mclapply and copies nothing, doParallelSNOW does not.
  fk <- rnaparallel:::rp_copy_free
  skip_on_os("windows")
  withr::with_options(list(combat.fork = TRUE), {
    expect_true(fk("mclapply"))
    expect_true(fk("BiocParallel"))
    expect_false(fk("serial"))
    expect_true(fk(function(idx, f, workers) lapply(idx, f)))  # custom keeps the inherited answer
  })
  # foreach depends on what is registered, which is the whole point
  skip_if_not_installed("doParallel")
  withr::with_options(list(combat.fork = TRUE), {
    foreach::registerDoSEQ()
    expect_false(fk("foreach"))              # unregistered: it gets this package's own cluster
    doParallel::registerDoParallel(cores = 2)
    expect_true(fk("foreach"))               # doParallelMC drives mclapply, nothing is copied
    cl <- parallel::makeCluster(2, type = "PSOCK")
    doParallel::registerDoParallel(cl)
    expect_false(fk("foreach"))              # doParallelSNOW serialises every task
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  })
  # the escape hatch turns every dispatch serial, so no gate should read as forking
  withr::with_options(list(combat.fork = FALSE), {
    expect_false(fk("mclapply"))
    expect_false(fk(function(idx, f, workers) lapply(idx, f)))
  })

  withr::with_options(list(combat.min.ls.cells = NULL, combat.min.norm.cells = NULL,
                           combat.min.order.cells = NULL), {
    expect_identical(rnaparallel:::rp_ls_min_cells("combat.min.ls.cells", 6e6, "mclapply"), 6e6)
    foreach::registerDoSEQ()
    expect_identical(rnaparallel:::rp_ls_min_cells("combat.min.ls.cells", 6e6, "foreach"), Inf)
    expect_identical(rnaparallel:::rp_norm_min_cells("mclapply"), 2e5)
    expect_identical(rnaparallel:::rp_norm_min_cells("foreach"), 2e6)
    expect_identical(rnaparallel:::rp_order_min_cells("mclapply"), 4e6)
    expect_identical(rnaparallel:::rp_order_min_cells("foreach"), 4e7)
  })
})
