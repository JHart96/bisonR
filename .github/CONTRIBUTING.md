# Contributing to bisonR

Thanks for helping improve bisonR! This guide explains how changes get from an idea into a release.

## Branches

* **`main`** is the latest release. It's protected: nothing is pushed to it directly, and every change arrives through a pull request (PR) that passes all checks. Anyone installing with `remotes::install_github("JHart96/bisonR")` gets `main`, so it must always work.
* **`release/<version>`** (e.g. `release/0.5.0`) collects the work for the next release. When the release is ready, a PR from this branch into `main` is reviewed and merged, and a GitHub release is tagged.
* **Feature and fix branches** (`feature/<short-name>`, `fix/<short-name>`) hold individual changes. Branch from the current release branch and open your PR against it. Small urgent fixes to a release can branch from `main` and target `main` directly, as long as they bump the patch version.

## Checks on every pull request

GitHub Actions runs these on every PR. A PR into `main` can only be merged when they all pass:

| Check | What it does |
|---|---|
| `R-CMD-check` | `R CMD check` on macOS, Windows and Linux, including the test suite and vignettes. |
| `test-coverage` | Measures test coverage and reports it to Codecov. |
| `pkgdown` | Builds the documentation site, so broken examples or articles are caught before merging. The site at [bisonr.org](https://bisonr.org) is only deployed from `main`. |
| `version-check` | PRs into `main` must increase the `Version` in `DESCRIPTION` and add a matching `# bisonR <version>` heading at the top of `NEWS.md`. |

## Versions

bisonR uses [semantic versioning](https://semver.org/) (`major.minor.patch`):

* **patch** (0.5.0 → 0.5.1): bug fixes that don't change how existing code is used.
* **minor** (0.5.0 → 0.6.0): new features, or changes that may alter results or require small changes to user code. These must be listed in `NEWS.md`.
* **major**: reserved for 1.0.0 and later breaking changes.

Bump the version with `usethis::use_version("patch")` (or `"minor"`), which also adds a NEWS heading for you.

## Making a change

1. Open an issue (or comment on an existing one) to discuss larger changes first.
2. Create a branch from the current release branch: `git switch -c feature/my-change origin/release/0.5.0`.
3. Make the change, with tests in `tests/testthat/` for any new behaviour or bug fix.
4. Run `devtools::document()` if you changed roxygen comments, then `devtools::check()`.
5. Add a bullet to `NEWS.md` under the version being worked on.
6. Push the branch and open a PR. Fill in the PR template checklist.

## Code style

* Follow the [tidyverse style guide](https://style.tidyverse.org/) for new code (`<-` for assignment, `snake_case` names).
* Call functions from other packages with `pkg::fn()`, and never use `library()` or `require()` inside package code.
* Give errors that tell the user what to do, e.g. ``stop("Column `x` not found in `data`.", call. = FALSE)``.

## Running the tests

Tests that fit Stan models need CmdStan installed and are skipped otherwise. Run everything with:

```r
devtools::test()
```

Compiled Stan models are cached in `tools::R_user_dir("bisonR", "cache")`. Delete that folder to force a recompile.
