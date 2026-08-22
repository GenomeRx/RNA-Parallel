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

  # the contract is which pids are in scope, not how fast waitpid clears the table: a killed
  # child stays listed by children() until it is collected, so asserting on that list races.
  expect_identical(reap(spare), 3L)
  expect_true(outsider$pid %in% kids())        # the caller's own worker is untouched

  expect_gte(reap(integer()), 1L)              # with nothing spared, the outsider is in scope
  invisible(parallel::mccollect(wait = FALSE))
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
