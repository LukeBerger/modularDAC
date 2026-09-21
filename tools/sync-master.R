#!/usr/bin/env Rscript
#
# Regenerate the `master` branch from `dev`.
#
# master is the user-facing branch: it carries only what someone needs to run
# modularDAC on their own experimental data. The benchmarking code -- anything
# that requires a known ground-truth graph, or that exists to score a result
# against one -- stays on dev.
#
# The split has two parts:
#
#   1. DROPPED entirely. Functions whose only purpose is benchmarking or
#      reporting on a benchmark, together with their tests and man pages.
#
#   2. KEPT BUT UNEXPORTED. The graph/data simulators and the ground-truth
#      module constructors. Every test fixture in the suite builds its data with
#      these, so removing them would take the whole test suite with them. On
#      master they are marked @keywords internal instead of @export, so they
#      still back the tests but never appear in the user-facing API.
#
# This script is the single source of truth for that transformation: it is
# re-applied in full on every sync, so the two branches cannot drift and merge
# conflicts on the transformed files do not need hand-resolution.
#
# Usage, from a clean checkout of the package root:
#     Rscript tools/sync-master.R           # build master locally
#     Rscript tools/sync-master.R --push    # build master and push to origin
#
# Pushing is opt-in. Without --push the script stops after committing so the
# diff can be reviewed.

args <- commandArgs(trailingOnly = TRUE)
do.push <- "--push" %in% args

# ---------------------------------------------------------------- the split --

# dropped from master (source, tests and generated docs)
DROP_R <- c(
  "R/modular_plot.R",        # visNetwork plot of a graph coloured by module
  "R/calc_F1.R",             # scores a predicted graph against a true one
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
  "man/modular_plot.Rd", "man/calc_F1.Rd", "man/module_contiguity.Rd",
  "man/module_correlation.Rd", "man/module_match.Rd", "man/dot-match_modules.Rd"
)

# kept on master but unexported: test infrastructure, not user-facing API
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

# ---------------------------------------------------------------- utilities --

run <- function(...) {
  cmd <- paste(...)
  out <- suppressWarnings(system(cmd, intern = TRUE))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
    stop("command failed: ", cmd, "\n", paste(out, collapse = "\n"), call. = FALSE)
  }
  invisible(out)
}

step <- function(...) cat("==> ", ..., "\n", sep = "")

if (!file.exists("DESCRIPTION")) {
  stop("run this from the package root (the directory holding DESCRIPTION).", call. = FALSE)
}
if (length(suppressWarnings(system("git status --porcelain", intern = TRUE))) > 0) {
  stop("working tree is not clean; commit or stash first.", call. = FALSE)
}

# ------------------------------------------------------- merge dev to master --

step("checking out master and merging dev")
run("git checkout master")
run("git fetch origin")
# fast-forward master to whatever origin has, then take dev's content wholesale.
# -X theirs settles content conflicts in dev's favour; modify/delete conflicts on
# files this script removes anyway are resolved by the removal step below.
suppressWarnings(system("git merge --ff-only origin/master", intern = TRUE))
merged <- suppressWarnings(system("git merge -X theirs --no-edit dev", intern = TRUE))
cat(paste(merged, collapse = "\n"), "\n")
# any remaining conflict is a modify/delete on a file master drops: take the delete
conflicts <- suppressWarnings(system("git diff --name-only --diff-filter=U", intern = TRUE))
for (f in conflicts) {
  if (f %in% c(DROP_R, DROP_TESTS, DROP_MAN)) {
    run("git rm -q -f", shQuote(f))
  } else {
    run("git checkout --theirs", shQuote(f))
    run("git add", shQuote(f))
  }
}
if (length(conflicts)) run("git commit --no-edit")

# --------------------------------------------------------- apply the split --

step("removing benchmarking sources, tests and man pages")
for (f in c(DROP_R, DROP_TESTS, DROP_MAN)) {
  if (file.exists(f)) run("git rm -q -f", shQuote(f))
}

step("unexporting the simulators and ground-truth constructors")
for (f in UNEXPORT) {
  if (!file.exists(f)) next
  txt <- readLines(f, warn = FALSE)
  hit <- grepl("^#' @export\\s*$", txt)
  if (!any(hit)) next
  txt[hit] <- paste0(
    "#' @keywords internal\n",
    "#' @note Benchmarking helper. Present on this branch only to back the test\n",
    "#'   suite; it is not part of the user-facing API. The exported version\n",
    "#'   lives on the dev branch."
  )
  writeLines(txt, f)
}

step("trimming Suggests: ", paste(DROP_SUGGESTS, collapse = ", "))
d <- readLines("DESCRIPTION", warn = FALSE)
drop.line <- vapply(d, function(l) {
  any(vapply(DROP_SUGGESTS, function(p)
    grepl(paste0("^\\s+", p, ",?\\s*$"), l), logical(1)))
}, logical(1))
writeLines(d[!drop.line], "DESCRIPTION")

step("regenerating NAMESPACE and man/")
roxygen2::roxygenise(".")

step("committing")
run("git add -A")
msg <- paste(
  "Rebuild master from dev without benchmarking code",
  "",
  "master carries only what a user needs to run modularDAC on their own data.",
  "Dropped: modular_plot, calc_F1, module_contiguity, module_correlation,",
  "module_match (and their tests/docs). The simulators and ground-truth module",
  "constructors are kept but unexported, since every test fixture depends on",
  "them. Generated by tools/sync-master.R.",
  sep = "\n")
run("git commit -q -m", shQuote(msg))

step("done. master now has ",
     length(system("git ls-tree --name-only HEAD R/", intern = TRUE)), " R/ scripts and ",
     length(grep("^export", readLines("NAMESPACE", warn = FALSE))), " exports.")

if (do.push) {
  step("pushing to origin/master")
  run("git push origin master")
} else {
  step("not pushed. Review with:  git diff origin/master..master --stat")
  step("then push with:           git push origin master")
}
