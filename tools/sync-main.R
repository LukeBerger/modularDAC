#!/usr/bin/env Rscript
#
# Regenerate the `main` branch from `dev`.
#
# main is the user-facing branch: it carries only what someone needs to run
# modularDAC on their own experimental data. The benchmarking code -- anything
# that requires a known ground-truth graph, or that exists to score a result
# against one -- stays on dev.
#
# The split has two parts:
#
#   1. DROPPED entirely. Functions whose only purpose is benchmarking or
#      reporting on a benchmark, together with their tests and man pages.
#
#   2. KEPT BUT UNEXPORTED. The graph/data simulators, the ground-truth module
#      constructors and calc_F1. Every test fixture in the suite builds its data
#      with these, so removing them would take the whole test suite with them.
#      On main they are marked @keywords internal instead of @export, so they
#      still back the tests but never appear in the user-facing API.
#
# This script is the single source of truth for that transformation: it is
# re-applied in full on every sync, so the two branches cannot drift and merge
# conflicts on the transformed files never need hand-resolution.
#
# All main-side work happens in a temporary git worktree, so the checkout this
# runs from stays on dev and is never modified. (Switching branches in place
# would delete this very file mid-run, which aborts the merge.)
#
# Usage, from a clean checkout of the package root:
#     Rscript tools/sync-main.R           # build main locally
#     Rscript tools/sync-main.R --push    # build main and push to origin
#
# Pushing is opt-in. Without --push the script stops after committing so the
# diff can be reviewed.

args <- commandArgs(trailingOnly = TRUE)
do.push <- "--push" %in% args

# ---------------------------------------------------------------- the split --

# dropped from main (source, tests and generated docs)
DROP_R <- c(
  "R/modular_plot.R",        # visNetwork plot of a graph coloured by module
  "R/module_contiguity.R",   # how self-contained modules are within a graph
  "R/module_correlation.R",  # within-module correlation distributions
  "R/module_match.R"         # agreement between two module sets
)
DROP_TESTS <- c(
  "tests/testthat/test-module_contiguity.R",
  "tests/testthat/test-module_correlation.R",
  "tests/testthat/test-module_match.R"
)
DROP_MAN <- c(
  "man/modular_plot.Rd", "man/module_contiguity.Rd",
  "man/module_correlation.Rd", "man/module_match.Rd", "man/dot-match_modules.Rd"
)
# dev-only tooling. The example dataset it produces (data/lfr_example.rda, with
# its R/lfr_example.R docs) ships on both branches; the script that regenerates
# it does not, since it calls the simulators.
DROP_OTHER <- c(
  "data-raw/make_example_data.R"
)

# kept on main but unexported: test infrastructure, not user-facing API.
# calc_F1 is here rather than in DROP_R because test-learn_SILGGM_graph.R uses
# it as a sanity check on the learned graph.
UNEXPORT <- c(
  "R/make_modular_graph.R",
  "R/make_lfr.R",
  "R/sim_graph_data.R",
  "R/true_modules.R",
  "R/true_fuzzy.R",
  "R/calc_F1.R"
)

# Suggests entries no longer reachable once the dropped set is gone
DROP_SUGGESTS <- c("ggplot2", "visNetwork")

# files that must exist on main after a correct merge; guards against the
# transformation being applied to stale content if the merge silently fails
SENTINELS <- c("R/divide_and_conquer.R", "R/learn_ARACNE_graph.R",
               "R/module-class.R", "tests/testthat/helper-fixtures.R")

# ---------------------------------------------------------------- utilities --

# stderr is left to flow to the console rather than captured, so git's own
# messages stay visible; the exit status is what decides success.
run <- function(..., allow.fail = FALSE) {
  cmd <- paste(...)
  out <- suppressWarnings(system(cmd, intern = TRUE))
  st <- attr(out, "status")
  if (is.null(st)) st <- 0L
  if (!allow.fail && st != 0) {
    stop("command failed (status ", st, "): ", cmd, call. = FALSE)
  }
  attr(out, "status") <- st
  out
}

step <- function(...) cat("==> ", ..., "\n", sep = "")

if (!file.exists("DESCRIPTION")) {
  stop("run this from the package root (the directory holding DESCRIPTION).", call. = FALSE)
}
if (length(run("git status --porcelain")) > 0) {
  stop("working tree is not clean; commit or stash first.", call. = FALSE)
}

root <- normalizePath(".", winslash = "/")
wt <- file.path(normalizePath(tempdir(), winslash = "/"), "modularDAC-main")

# ------------------------------------------------------ main worktree setup --

step("fetching")
run("git fetch origin")

step("preparing a main worktree at ", wt)
run("git worktree remove --force", shQuote(wt), allow.fail = TRUE)  # leftover from a failed run
run("git worktree prune")
run("git worktree add", shQuote(wt), "main")
on.exit({
  setwd(root)
  run("git worktree remove --force", shQuote(wt), allow.fail = TRUE)
  run("git worktree prune", allow.fail = TRUE)
}, add = TRUE)

setwd(wt)

run("git merge --ff-only origin/main", allow.fail = TRUE)  # no-op when in sync

step("merging dev")
# -X theirs settles content conflicts in dev's favour; the transformation below
# is re-applied in full afterwards, so main's own edits never need preserving.
merge.out <- run("git merge -X theirs --no-edit dev", allow.fail = TRUE)
if (attr(merge.out, "status") != 0) {
  conflicts <- run("git diff --name-only --diff-filter=U", allow.fail = TRUE)
  if (!length(conflicts)) {
    stop("merge of dev failed and left no conflicts to resolve; see git's ",
         "output above. Nothing has been committed.", call. = FALSE)
  }
  # the only expected conflicts are modify/delete on files main drops
  for (f in conflicts) {
    if (f %in% c(DROP_R, DROP_TESTS, DROP_MAN, DROP_OTHER)) {
      run("git rm -q -f", shQuote(f))
    } else {
      run("git checkout --theirs", shQuote(f)); run("git add", shQuote(f))
    }
  }
  run("git commit --no-edit")
}

missing <- SENTINELS[!file.exists(SENTINELS)]
if (length(missing)) {
  stop("main does not look like dev after the merge; missing:\n  ",
       paste(missing, collapse = "\n  "),
       "\nAborting before any files are changed.", call. = FALSE)
}
behind <- as.integer(run("git rev-list --count HEAD..dev")[1])
if (behind > 0) {
  stop("main is still ", behind, " commit(s) behind dev after the merge; aborting.",
       call. = FALSE)
}

# ---------------------------------------------------------- apply the split --

step("removing benchmarking sources, tests, man pages and dev-only tooling")
for (f in c(DROP_R, DROP_TESTS, DROP_MAN, DROP_OTHER)) {
  if (file.exists(f)) run("git rm -q -f", shQuote(f))
}

step("unexporting the simulators, ground-truth constructors and calc_F1")
for (f in UNEXPORT) {
  if (!file.exists(f)) next
  txt <- readLines(f, warn = FALSE)
  hit <- grepl("^#' @export\\s*$", txt)
  if (!any(hit)) next            # already transformed; the script is idempotent
  txt[hit] <- paste0(
    "#' @keywords internal\n",
    "#' @note Benchmarking helper, kept on this branch only to back the test\n",
    "#'   suite. It is not part of the user-facing API; the exported version\n",
    "#'   lives on the dev branch.")
  writeLines(txt, f)
}

step("trimming Suggests: ", paste(DROP_SUGGESTS, collapse = ", "))
d <- readLines("DESCRIPTION", warn = FALSE)
drop.line <- vapply(d, function(l) any(vapply(DROP_SUGGESTS, function(p)
  grepl(paste0("^\\s+", p, ",?\\s*$"), l), logical(1))), logical(1))
writeLines(d[!drop.line], "DESCRIPTION")

step("regenerating NAMESPACE and man/")
roxygen2::roxygenise(".")

step("running the test suite on main")
res <- as.data.frame(testthat::test_local(".", reporter = "silent"))
cat("    tests: ", nrow(res), "  passed: ", sum(res$passed),
    "  failed: ", sum(res$failed), "  errors: ", sum(res$error), "\n", sep = "")
if (sum(res$failed) + sum(res$error) > 0) {
  stop("main's test suite does not pass; nothing has been committed.", call. = FALSE)
}

step("committing")
run("git add -A")
# Re-running with nothing new to sync leaves nothing staged, and `git commit`
# would fail. Skip the commit in that case so the script stays idempotent and
# `--push` still works on a second run.
if (!length(run("git diff --cached --name-only"))) {
  step("nothing to commit; main already matches dev")
} else {
  msg <- paste(
    "Rebuild main from dev without benchmarking code",
    "",
    "main carries only what a user needs to run modularDAC on their own data.",
    "Dropped: modular_plot, module_contiguity, module_correlation, module_match",
    "(and their tests and docs). The simulators, the ground-truth module",
    "constructors and calc_F1 are kept but unexported, since every test fixture",
    "depends on them. Generated by tools/sync-main.R.",
    sep = "\n")
  run("git commit -q -m", shQuote(msg))
}

n.scripts <- length(run("git ls-tree --name-only HEAD R/"))
n.exports <- length(grep("^export", readLines("NAMESPACE", warn = FALSE)))
step("done. main has ", n.scripts, " R/ scripts and ", n.exports, " exports.")

if (do.push) {
  step("pushing to origin/main")
  run("git push origin main")
} else {
  step("not pushed. Review with:  git log --oneline origin/main..main")
  step("                          git diff origin/main..main --stat")
  step("then push with:           git push origin main")
}
