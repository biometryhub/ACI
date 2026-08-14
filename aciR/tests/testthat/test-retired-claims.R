# -- claims that were retired, and must not come back --------------------------
#
# Twice in review this package was found asserting something that had ceased to
# be true: documentation saying the online smoother and causal influence range
# were scalar-only, months after they shipped, and a `@param` promising `NA`
# where the code had begun returning a censored bound. Both times the claim had
# been corrected on the surfaces the author happened to touch and left standing
# on the others.
#
# The failure is not the individual sentences. It is marking an item done from
# the surfaces one edited rather than from a search over the surfaces that could
# carry it -- a check that costs one grep and was not run. A habit of running it
# is exactly the kind of rule that fails under load, so it is a test instead.
#
# Each entry is a RETIRED SENTENCE, not a keyword. That is deliberate: the
# package is allowed, and often obliged, to discuss a rule it no longer uses.
# What it may not do is assert it in the present tense.

.retired_claims <- list(
  list(
    phrase = "remain scalar-only",
    why = paste(
      "The online smoother and causal influence range dispatch on the",
      "components exactly as the core does, and have since the vector work",
      "landed."
    )
  ),
  list(
    phrase = "still does not implement",
    why = paste(
      "Used to introduce the vector and conditional cases, both of which",
      "exist."
    )
  ),
  list(
    phrase = "Not implemented, so not graded",
    why = paste(
      "The vector online smoother and CIR are implemented; their grounding is",
      "collapse onto the scalar path plus an independent transcription, which",
      "is a scope statement rather than an absence."
    )
  ),
  list(
    phrase = "ranges returned as `NA`",
    why = paste(
      "A censored time returns a lower bound. Only insufficient status",
      "yields NA."
    )
  ),
  list(
    phrase = "the grid of the reference implementation",
    why = paste(
      "The default epsilon grid has 129 points; the reference uses 513. The",
      "default is this package's cheaper grid and is documented as such."
    )
  ),
  list(
    phrase = "closed with the\n#' Simpson 3/8 rule",
    why = paste(
      "The odd-interval closure follows the reference: the quadratic through",
      "the last three samples over the final interval."
    )
  )
)

# Gather every text surface that could carry a claim, from wherever it is
# reachable. Under `R CMD check` the source tree is gone and only the installed
# package exists, so the two roots are unioned rather than chosen between --
# and the manual is read through `Rd_db()`, since installed help is a database
# rather than files. A guard that quietly scans nothing is worse than no guard,
# so the count is asserted before the search.
.claim_surfaces <- function() {
  roots <- unique(c(
    normalizePath(testthat::test_path("..", ".."), mustWork = FALSE),
    system.file(package = "aciR")
  ))
  roots <- roots[nzchar(roots) & dir.exists(roots)]

  paths <- character(0)
  for (root in roots) {
    paths <- c(
      paths,
      list.files(file.path(root, "R"), "\\.R$", full.names = TRUE),
      list.files(file.path(root, "man"), "\\.Rd$", full.names = TRUE),
      list.files(file.path(root, "vignettes"), "\\.Rmd$", full.names = TRUE),
      list.files(file.path(root, "doc"), "\\.(Rmd|html)$", full.names = TRUE),
      list.files(file.path(root, "inst", "extdata"), "\\.yml$",
                 full.names = TRUE),
      list.files(file.path(root, "extdata"), "\\.yml$", full.names = TRUE),
      file.path(root, c("DESCRIPTION", "README.Rmd", "README.md", "NEWS.md"))
    )
  }
  paths <- unique(paths[file.exists(paths)])

  texts <- vapply(paths, function(f) {
    paste(readLines(f, warn = FALSE), collapse = "\n")
  }, character(1L))
  names(texts) <- basename(paths)

  # The installed manual is an Rd database, not a directory of files.
  manual <- tryCatch(tools::Rd_db("aciR"), error = function(e) list())
  if (length(manual) > 0L) {
    rd <- vapply(manual, function(entry) {
      paste(as.character(entry), collapse = " ")
    }, character(1L))
    names(rd) <- paste0("Rd:", names(manual))
    texts <- c(texts, rd)
  }
  texts
}

test_that("retired claims have not reappeared anywhere in the package", {
  texts <- .claim_surfaces()

  # The floor holds in both environments: a source tree offers dozens, an
  # installed package offers DESCRIPTION, NEWS, the manual database and the
  # oracle manifest. Zero means the search found nothing to search.
  expect_gt(length(texts), 5L)

  offences <- character(0)
  for (claim in .retired_claims) {
    hit <- names(texts)[vapply(texts, grepl, logical(1L),
                               pattern = claim$phrase, fixed = TRUE)]
    for (where in hit) {
      offences <- c(offences, sprintf(
        "%s\n    contains retired claim: \"%s\"\n    %s",
        where, claim$phrase, claim$why
      ))
    }
  }
  expect_identical(
    offences, character(0),
    info = paste0("\n  ", paste(offences, collapse = "\n  "))
  )
})
