# Repository Guidelines

## Project Structure & Module Organization
BAGEL is an R package. Core functions live in `R/` (e.g., `calculate_copynumber.R`, `run_biscut.R`) and expose exported APIs via roxygen tags. Datasets and breakpoints are stored as `.rda` objects in `data/`, with reproducible generation scripts in `data-raw/`. Tests reside in `tests/testthat/` alongside `testthat.R`. User-facing docs and vignettes are in `docs/` and `Tutorial.Rmd`; sample run scripts live in `runMe.R` and assets under `inst/`.

## Build, Test, and Development Commands
- `Rscript -e "devtools::document()"`: regenerate roxygen2 docs before committing API changes.
- `Rscript -e "devtools::load_all()"`: load the package locally for quick iteration.
- `Rscript -e "devtools::test()"`: run the testthat suite in `tests/testthat`.
- `Rscript -e "devtools::check()"` or `R CMD check .`: full package check (build, lint, examples).
- `R CMD build .`: produce a release tarball (mirrors `BAGEL_2.0.0.tar.gz`).

## Coding Style & Naming Conventions
Follow tidyverse style: 4-space indentation, spaces around `=` in named arguments, and pipe-friendly, readable verbs. Exported functions use `snake_case` (`calculateCopyNumber_fixed` is legacy; prefer `snake_case` for new code). Keep files focused by topic. Annotate every exported function with roxygen2 comments, including `@param`, `@return`, and `@examples`. Place helper utilities in `R/utils.R` and keep internal objects prefixed with a leading dot if hidden.

## Testing Guidelines
Write unit tests with `testthat::test_that()` under `tests/testthat/`; mirror the file being tested (e.g., `test-calculate_copynumber.R`). Include representative TCGA cancer types using fixtures from `inst/extdata`. Ensure new features update or extend snapshot expectations. Run `Rscript -e "devtools::test()"` locally and attach coverage notes (via `covr::report()` if feasible) to PRs touching analytical logic.

## Commit & Pull Request Guidelines
Recent history is terse (`MAJOR UPDATES!!`); switch to imperative, descriptive commits such as `Add arm summarization helper`. Group related changes and avoid bundling regenerated docs unless code changed. PRs should list motivation, key changes, test evidence (`devtools::check` output), and any data dependencies or new files. Tag linked issues and add screenshots when modifying reports or plots.

## Data & Configuration Tips
Keep raw data generation scripts in `data-raw/` reproducible; rerun them before updating `.rda` data objects. When introducing new configuration defaults, expose them via exported functions, document them in `man/`, and surface user-facing knobs in `run_biscut.R` examples.
