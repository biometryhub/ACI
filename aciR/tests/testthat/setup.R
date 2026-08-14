# Partial matching of `$` silently rescued a test that reached for `aux$E` on
# an object whose field is `E_j`. It passed for the wrong reason, and would
# have turned into what looked like a numerical regression the moment a field
# named `E` appeared. Warnings are errors under testthat, so this makes the
# class un-reintroducible rather than merely fixed once.
options(warnPartialMatchDollar = TRUE)
