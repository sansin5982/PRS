# Verification record

## Completed in this build

- Twenty chapter sources are present, each with YAML and references.
- Pandoc parsed every chapter's Markdown structure successfully.
- 138 R chunk labels were checked for uniqueness across the course.
- A lexical delimiter check passed for R chunks and shared scripts. This is not an R parser or execution test.
- README/index chapter targets and shared helper paths were checked.
- Selected worked arithmetic, participant joins, missingness denominators and synthetic role counts were checked independently.
- Scientific review distinguishes simulation, published findings, software templates, assumptions and unassessed questions.
- There are 33 figure-captioned R chunks. Their plotting code is supplied; images are not rendered in this build.
- Primary method papers, PGS Catalog records and official PLINK, LDpred2, PRS-CS, PRS-CSx and R Markdown documentation were consulted. Some publisher pages blocked direct retrieval; a working identifier is not a guarantee of unrestricted full-text access.

## Not executed here

No R executable is available. R chunks, package-specific API calls, knitting, figure rendering, PLINK commands, PRS-CS/PRS-CSx runs and the end-to-end runner have NOT been executed in this environment. Numerical output from those analyses is not fabricated. Optional external inputs are not bundled or claimed downloaded.

## Local completion gates

1. Run `Rscript scripts/run_analysis.R` from a fresh session at the course root.
2. Inspect its output tables, warnings and audit. A successful process does not prove clinical validity.
3. Install rmarkdown/knitr and run `Rscript scripts/render_course.R`.
4. Inspect every generated figure and page; adjust label spacing if your device differs.
5. For optional genomic software exercises, record exact versions, stage documented inputs and compare the PLINK exercise with expected_scores.tsv.

## Scope and changes

This edition incorporates the previously reviewed foundations, corrects Chapter 7's export instructions and readiness labels, strengthens scoring-table checks, adds an unmatched-ID example and adds Chapters 8–20 with shared synthetic data and a separate analysis runner. Earlier standalone chapters remain separate artifacts.

The first seven chapters are broad lecture notes with additional mini-exercises. Later chapters are focused method or evaluation lessons. The course does not claim to cover every PRS algorithm or provide clinical validation.
