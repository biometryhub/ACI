# render the side-by-side parity document --------------------------------------
#
# Generated from the manifests and the measured reports, never hand-written, so
# the document cannot drift from what was actually run. The MATLAB excerpts are
# read from the reference by the line ranges the manifest declares; the R
# excerpts are read from the package by function name.

.esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

# Pull an R function body by name: from its definition to the closing brace in
# column one. The codebase is uniform enough for this to be exact.
.r_function <- function(file, name) {
  src <- readLines(file, warn = FALSE)
  start <- grep(sprintf("^%s <- function", name), src)[1L]
  if (is.na(start)) return(sprintf("## %s not found in %s", name, file))
  end <- start + which(src[seq.int(start, length(src))] == "}")[1L] - 1L
  paste(src[seq.int(start, end)], collapse = "\n")
}

.matlab_lines <- function(script, spec) {
  src <- sub("\r$", "", readLines(file.path("matlab_reference", script),
                                  warn = FALSE))
  parts <- trimws(strsplit(spec, ",", fixed = TRUE)[[1L]])
  idx <- unlist(lapply(parts, function(p) {
    e <- as.integer(trimws(strsplit(p, "-", fixed = TRUE)[[1L]]))
    seq.int(e[1L], e[2L])
  }))
  paste(src[idx], collapse = "\n")
}

.pane <- function(title, sub, code) sprintf(
  '<div class="pane"><div class="pane-h"><span class="pane-t">%s</span><span class="pane-s">%s</span></div><pre><code>%s</code></pre></div>',
  .esc(title), .esc(sub), .esc(code))

# `companion` is the href to the development ledger, which differs by output
# location: the in-tree record sits four directories below the repository root,
# the delivery copy sits beside its companion. Rendering rather than copying is
# what keeps the two from drifting apart.
#
# The delivery copy goes to `aciR/pkgdown/assets/ledgers/`, which pkgdown copies
# verbatim into the built site. That directory is the single home for these
# documents: they are HTML, so GitHub serves them as source rather than as
# pages, and the rendered site is the only place they are actually readable.
# Use `render_parity_delivery()` rather than retyping the pair by hand, so that
# a correct output path cannot depend on remembering this comment.
render_parity <- function(out = "tools/oracle/parity/reports/parity.html",
                          companion = paste0("../../../../aciR/pkgdown/",
                                             "assets/ledgers/",
                                             "development_ledger.html")) {
  ex <- read_manifest("tools/oracle/parity/manifest/extracts.dcf")
  # Every results file in reports/ is tabulated, not just the scalar one. The
  # predator-prey comparison was invisible on this page for exactly as long as
  # the list was hard-coded to one file.
  res_files <- c("parity_scalar.csv", "parity_predprey.csv")
  res_files <- file.path("tools", "oracle", "parity", "reports", res_files)
  res <- do.call(rbind, lapply(res_files[file.exists(res_files)],
                               utils::read.csv, stringsAsFactors = FALSE))

  pairs <- list(
    list(m = "ref_filter_scalar", r = "aci_filter", f = "aciR/R/aci-core.R",
         cap = "The forward CGNS filter. Note the reference's two initial conditions are literals inside the range: the mean starts at the true latent y(1), the covariance at a hardcoded 0.1. aciR takes both as arguments."),
    list(m = "ref_aci_metric_scalar", r = "aci_metric", f = "aciR/R/aci-core.R",
         cap = "The ACI metric as a relative entropy of the smoother from the filter. Four lines in MATLAB; aciR routes through a shared implementation so that aci() and aci_metric() cannot drift apart."),
    list(m = "ref_pp_filter_dir1", r = "aci_filter", f = "aciR/R/aci-core.R",
         cap = "Predator-prey, direction x(t) to y. The same algebra as the dyad with the roles of x and y exchanged, and the same variable names as its own direction two, which is why that script cannot be run from top to bottom.")
  )

  blocks <- vapply(pairs, function(p) {
    i <- which(ex[, "Function"] == p$m)[1L]
    script <- ex[i, "Script"]; spec <- ex[i, "Lines"]
    sprintf('<section class="pair"><h3>%s</h3><p class="cap">%s</p><div class="grid">%s%s</div></section>',
      # The arrow is written as an entity, not as a literal, so the rendered
      # page stays pure ASCII and cannot be mis-decoded by a browser that
      # guesses the encoding. It sits outside .esc(), which would escape the &.
      paste0(.esc(p$m), " &harr; ", .esc(paste0(p$r, "()"))), .esc(p$cap),
      .pane("MATLAB (reference)", sprintf("%s:%s", script, spec),
            .matlab_lines(script, spec)),
      .pane("R (aciR)", sprintf("%s :: %s()", basename(p$f), p$r),
            .r_function(p$f, p$r)))
  }, character(1L))

  rows <- paste(apply(res, 1L, function(r) sprintf(
    '<tr><td>%s</td><td class="q">%s</td><td class="n">%s</td><td class="n">%s</td><td class="n">%s</td><td><span class="v %s">%s</span></td></tr>',
    .esc(r[["dataset"]]), .esc(r[["quantity"]]), .esc(r[["n"]]),
    .esc(format(as.numeric(r[["max_abs"]]), digits = 3, scientific = TRUE)),
    # "n/a" rather than a dash: a bare "--" reaches the page literally, since
    # nothing here processes typographic markup the way pandoc would.
    .esc(if (is.na(as.numeric(r[["headroom"]]))) "n/a"
         else format(as.numeric(r[["headroom"]]), digits = 3)),
    if (r[["verdict"]] == "ok") "ok" else "bad", .esc(r[["verdict"]]))),
    collapse = "\n")

  # Token substitution rather than sprintf: the template is far past
  # sprintf's 8192-character format limit, and its CSS contains percent signs.
  template <- paste(readLines("tools/oracle/parity/tools/parity_template.html",
                              warn = FALSE), collapse = "\n")
  html <- sub("{{PAIRS}}", paste(blocks, collapse = "\n"), template,
              fixed = TRUE)
  html <- sub("{{ROWS}}", rows, html, fixed = TRUE)
  html <- sub("{{COMPANION}}", companion, html, fixed = TRUE)

  # The headline statistics are computed here, not typed into the template.
  # The code panes were already pulled at render time while these cards were
  # hardcoded, so the page could claim it cannot drift while carrying a figure
  # that could. Deriving them from the manifest and the gate reports closes
  # that gap: a stale card is now impossible rather than merely unlikely.
  chk <- check_extraction("tools/oracle/parity/manifest/extracts.dcf")
  g1_files <- list.files(file.path("tools", "oracle", "parity", "reports"),
                         pattern = "^g1_", full.names = TRUE)
  g1 <- do.call(rbind, lapply(g1_files, utils::read.csv,
                              stringsAsFactors = FALSE))
  html <- sub("{{G1}}", sprintf("%d / %d", sum(g1$same == 1L), nrow(g1)),
              html, fixed = TRUE)
  html <- sub("{{VERBATIM_LINES}}",
              format(sum(chk$lines_verbatim), big.mark = ","), html,
              fixed = TRUE)
  html <- sub("{{N_EXTRACTS}}", as.character(nrow(chk)), html, fixed = TRUE)

  writeLines(html, out)
  out
}

# ---- delivery copy ---------------------------------------------------------

# The published parity ledger. It lands beside the development ledger in the
# pkgdown asset directory, so both are copied into the built site and their
# reciprocal links are plain siblings. Rendering it through this function is
# what stops the published copy and the in-tree record from diverging, and what
# stops a hand-typed relative path from being wrong.
render_parity_delivery <- function(
    dir = file.path("aciR", "pkgdown", "assets", "ledgers")) {
  if (!dir.exists(dir)) {
    stop("ledger directory not found: ", dir, call. = FALSE)
  }
  companion <- file.path(dir, "development_ledger.html")
  if (!file.exists(companion)) {
    stop("companion ledger not found: ", companion, call. = FALSE)
  }
  render_parity(out = file.path(dir, "parity_ledger.html"),
                companion = "development_ledger.html")
}
