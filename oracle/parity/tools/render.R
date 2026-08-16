# -- render the side-by-side parity document -----------------------------------
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
# location: the in-tree record sits three levels below `public/`, the delivery
# copy sits beside it. Rendering rather than copying is what keeps the two from
# drifting apart.
render_parity <- function(out = "oracle/parity/reports/parity.html",
                          companion = paste0("../../../public/",
                                             "aciR_development_ledger_2026-08-15.html")) {
  ex <- read_manifest("oracle/parity/manifest/extracts.dcf")
  # Every results file in reports/ is tabulated, not just the scalar one. The
  # predator-prey comparison was invisible on this page for exactly as long as
  # the list was hard-coded to one file.
  res_files <- c("parity_scalar.csv", "parity_predprey.csv")
  res_files <- file.path("oracle", "parity", "reports", res_files)
  res <- do.call(rbind, lapply(res_files[file.exists(res_files)],
                               utils::read.csv, stringsAsFactors = FALSE))

  pairs <- list(
    list(m = "ref_filter_scalar", r = "aci_filter", f = "aciR/R/aci-core.R",
         cap = "The forward CGNS filter. Note the reference's two initial conditions are literals inside the range: the mean starts at the true latent y(1), the covariance at a hardcoded 0.1. aciR takes both as arguments."),
    list(m = "ref_aci_metric_scalar", r = "aci_metric", f = "aciR/R/aci-core.R",
         cap = "The ACI metric as a relative entropy of the smoother from the filter. Four lines in MATLAB; aciR routes through a shared implementation so that aci() and aci_metric() cannot drift apart."),
    list(m = "ref_pp_filter_dir1", r = "aci_filter", f = "aciR/R/aci-core.R",
         cap = "Predator-prey, direction x(t) to y. The same algebra as the dyad with the roles of x and y exchanged -- and the same variable names as its own direction two, which is why that script cannot be run top to bottom.")
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
    .esc(if (is.na(as.numeric(r[["headroom"]]))) "--"
         else format(as.numeric(r[["headroom"]]), digits = 3)),
    if (r[["verdict"]] == "ok") "ok" else "bad", .esc(r[["verdict"]]))),
    collapse = "\n")

  # Token substitution rather than sprintf: the template is far past
  # sprintf's 8192-character format limit, and its CSS contains percent signs.
  template <- paste(readLines("oracle/parity/tools/parity_template.html",
                              warn = FALSE), collapse = "\n")
  html <- sub("{{PAIRS}}", paste(blocks, collapse = "\n"), template,
              fixed = TRUE)
  html <- sub("{{ROWS}}", rows, html, fixed = TRUE)
  html <- sub("{{COMPANION}}", companion, html, fixed = TRUE)
  writeLines(html, out)
  out
}
