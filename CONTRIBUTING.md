# Contributing to bisonR

This guide describes how changes to bisonR are made and released.

## Branches

- **`main`** is the latest release, and is what
  `remotes::install_github("JHart96/bisonR")` installs. It is protected:
  changes can only be made through pull requests (PRs) that pass the
  required checks.
- **`release/<version>`** (e.g. `release/1.1.0`) collects the work for
  the next release. When the release is ready, a PR from this branch
  into `main` is reviewed and merged, and a GitHub release is tagged.
- **Feature and fix branches** (`feature/<short-name>`,
  `fix/<short-name>`) hold individual changes. Branch from the current
  release branch and open your PR against it. Small urgent fixes to a
  release can branch from `main` and target `main` directly, as long as
  they bump the patch version.

## Checks on every pull request

GitHub Actions runs these checks on every PR. All except `test-coverage`
must pass before a PR can be merged into `main`:

| Check | What it does |
|----|----|
| `R-CMD-check` | `R CMD check` on macOS, Windows and Linux, including the tests and examples. Articles are not rebuilt here (the `pkgdown` check renders them), and the brms tests only run on Linux. |
| `test-coverage` | Measures test coverage and reports it to Codecov. |
| `pkgdown` | Builds the documentation site, so broken examples or articles are caught before merging. The site at [bisonr.org](https://bisonr.org) is only deployed from pushes to `main` and from published releases. |
| `version-check` | PRs into `main` must increase the `Version` in `DESCRIPTION` and add a matching `# bisonR <version>` heading at the top of `NEWS.md`. |

## Versions

bisonR uses [semantic versioning](https://semver.org/)
(`major.minor.patch`):

- **patch** (1.0.0 → 1.0.1): bug fixes that don’t change how existing
  code is used.
- **minor** (1.0.0 → 1.1.0): new features that don’t break existing
  code.
- **major** (1.0.0 → 2.0.0): breaking changes, such as removing or
  renaming arguments or functions, or changing results for the same
  input.

All changes should be listed in `NEWS.md`.

The version can be bumped with `usethis::use_version("patch")` (or
`"minor"` or `"major"`), which also adds a heading to `NEWS.md`.

## Making a change

1.  For larger changes, open an issue (or comment on an existing one)
    first.
2.  Create a branch from the current release branch:
    `git switch -c feature/my-change origin/release/<version>`.
3.  Make the change, with tests in `tests/testthat/` for any new
    behaviour or bug fix.
4.  Run `devtools::document()` if you changed roxygen comments, then
    `devtools::check()`.
5.  Add a bullet to `NEWS.md` under the version being worked on.
6.  Push the branch and open a PR. Fill in the PR template checklist.

## Code style

- Follow the [tidyverse style guide](https://style.tidyverse.org/) for
  new code (`<-` for assignment, `snake_case` names).
- Call functions from other packages with `pkg::fn()`, and never use
  [`library()`](https://rdrr.io/r/base/library.html) or
  [`require()`](https://rdrr.io/r/base/library.html) inside package
  code.
- Give errors that tell the user what to do,
  e.g. `` stop("Column `x` not found in `data`.", call. = FALSE) ``.

## Running the tests

Tests that fit Stan models need CmdStan installed and are skipped
otherwise. Run everything with:

``` r

devtools::test()
```

Compiled Stan models are cached in
`tools::R_user_dir("bisonR", "cache")`. Delete that folder to force a
recompile.
