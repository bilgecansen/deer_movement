# Repo conventions

## Script style

Two hard rules for every script here.

**1. Never use `;` as a statement separator.** One statement per line.

```r
# no
library(mgcv); library(amt)
id <- args[1]; season <- args[2]

# yes
library(mgcv)
library(amt)
id <- args[1]
season <- args[2]
```

A semicolon inside a string, or as English punctuation in a comment, is fine —
leave those alone. To find the real ones, use R's parser rather than grep:

```r
getParseData(parse(file, keep.source = TRUE))  # statement separators are ';' tokens
```

That never matches a semicolon inside a comment or a string literal.

In bash, `; then` and `; do` are required idiom and stay. Compressed one-liners
like `if ...; then x; else y; fi` should be expanded.

**2. Keep lines to 80 characters** where possible.

- Long messages split across several string arguments — `stop()`, `paste()` and
  `sprintf()` concatenate, so behaviour is unchanged.
- Long calls break at commas.
- Trailing comments move *above* the code rather than being padded out to the
  right margin.
- When rewrapping a long comment, reflow the whole paragraph. Splitting a single
  line in isolation leaves an orphaned word alone on the next line and reads
  worse than what you started with.

Irreducible exceptions: a line that is one long URL or file path.

`air.toml` pins `line-width = 80` for the formatter. Air has no semicolon rule,
so rule 1 is not enforced automatically.

## Layout

- `scripts/helpers/` — the helper library, split by theme. `helper_functions.R`
  is a thin aggregator that sources all of them; every script sources only that.
- `scripts/gam/`, `scripts/amt/` — one folder per modelling path, named for the
  fitting engine (mgcv GAM vs amt conditional logistic). Both paths fit an iSSF;
  the engine is what distinguishes them. Keep the `_gam` / `_amt` suffix on
  filenames even inside those folders; the redundancy is a deliberate safeguard.
  Note `amt::fit_issf` and `amt::make_issf_model` are the package's own API and
  keep their names.
- `scripts/retired/` — scripts kept for reference but not runnable as-is.
- `scripts/checks/` — correctness checks. `run_checks_all.sh` runs them all and
  exits non-zero on failure. `walkthrough_gam.R` is a straight-line, loop-free
  inlining of the GAM path for manual line-by-line review. `check_wrangle.R`
  re-derives every column of `data/tracks/` from the source rasters and the
  stored coordinates; run it as `check_wrangle.R all both` for the full cohort.
- `docs/gam_decision_inventory.md` — every non-obvious choice in the GAM path,
  attributed, with the alternative that was rejected.

Scripts run from the repository root; `source()` paths are relative to it.
