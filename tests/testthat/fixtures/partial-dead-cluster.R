# A child process cannot see the parent's devtools::load_all(), so a bare library() call here
# would test whatever rnaparallel is INSTALLED while the suite believes it is testing the
# working tree. Breaking combat_cluster_teardown in the source and re-running left this green.
local({
  src <- Sys.getenv("RNAPARALLEL_SRC", "")
  if (nzchar(src) && requireNamespace("pkgload", quietly = TRUE)) {
    suppressMessages(pkgload::load_all(src, quiet = TRUE))
  } else {
    library(rnaparallel)
  }
})

alive <- function(pid) isTRUE(tools::pskill(pid, 0L))
wait_gone <- function(pids, timeout = 4) {
  deadline <- proc.time()[["elapsed"]] + timeout
  while (any(vapply(pids, alive, logical(1))) &&
         proc.time()[["elapsed"]] < deadline) Sys.sleep(0.01)
  !any(vapply(pids, alive, logical(1)))
}

known <- integer()
on.exit({
  try(rnaparallel::combat_cluster_stop(), silent = TRUE)
  for (pid in known) if (alive(pid)) tools::pskill(pid, tools::SIGKILL)
}, add = TRUE)

rnaparallel::combat_cluster_stop()
cl <- rnaparallel:::combat_cluster(2L, "FORK")
old_pids <- unlist(parallel::clusterCall(cl, Sys.getpid), use.names = FALSE)
known <- c(known, old_pids)
cat("OLD_PIDS=", paste(old_pids, collapse = ","), "\n", sep = "")
flush.console()

trace_file <- tempfile("rnaparallel-post-node-")
assign(".rnaparallel_trace_file", trace_file, envir = .GlobalEnv)
assign(".rnaparallel_old_nodes", cl, envir = .GlobalEnv)
suppressMessages(trace(
  "postNode", where = asNamespace("parallel"), print = FALSE,
  tracer = quote({
    old_node <- any(vapply(.rnaparallel_old_nodes, identical, logical(1), y = con))
    if (old_node) cat("POST\n", file = .rnaparallel_trace_file, append = TRUE)
  })
))
on.exit(suppressMessages(untrace("postNode", where = asNamespace("parallel"))), add = TRUE)

tools::pskill(old_pids[1L], tools::SIGKILL)
stopifnot(wait_gone(old_pids[1L]))

# Force the first kill wait to time out after the node handles have been closed. Rebuild must
# continue, and later cleanup may use only the preserved PID identities.
ns <- asNamespace("rnaparallel")
old_wait_pids <- get("combat_wait_pids", envir = ns)
wait_calls <- 0L
assignInNamespace("combat_wait_pids", function(...) {
  wait_calls <<- wait_calls + 1L
  if (wait_calls == 1L) FALSE else old_wait_pids(...)
}, ns = "rnaparallel")
on.exit(assignInNamespace("combat_wait_pids", old_wait_pids, ns = "rnaparallel"), add = TRUE)

rebuilt <- rnaparallel:::combat_cluster(2L, "FORK")
post_writes <- if (file.exists(trace_file)) length(readLines(trace_file)) else 0L
cat("POST_WRITES_AFTER_SUSPICION=", post_writes, "\n", sep = "")
flush.console()
stopifnot(identical(post_writes, 0L))
# The surviving worker is either already gone -- closing its node socket is enough to end a
# FORK worker -- or held by PID for signal-only cleanup. Asserting it is ALWAYS held encoded
# one implementation's timing; what has to be true is that it is not still running and not
# still reachable through a connection.
retired <- rnaparallel:::.combat_clusters$retired
held <- if (length(retired)) unlist(lapply(retired, `[[`, "pids")) else integer()
stopifnot(!alive(old_pids[2L]) || old_pids[2L] %in% held)
cat("RETIRED_PIDS_AFTER_TIMEOUT=", paste(held, collapse = ","),
    " SURVIVOR_ALIVE=", alive(old_pids[2L]), "\n", sep = "")
flush.console()
stopifnot(identical(unlist(parallel::clusterCall(rebuilt, function() 2L + 2L)), c(4L, 4L)))
new_pids <- unlist(parallel::clusterCall(rebuilt, Sys.getpid), use.names = FALSE)
known <- c(known, new_pids)
cat("NEW_PIDS=", paste(new_pids, collapse = ","), "\n", sep = "")
flush.console()

rnaparallel::combat_cluster_stop()
stopifnot(wait_gone(c(old_pids, new_pids)))
post_writes <- if (file.exists(trace_file)) length(readLines(trace_file)) else 0L
stopifnot(identical(post_writes, 0L))
cat("PARTIAL_DEAD_CLUSTER_OK\n")
