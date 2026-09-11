Chapter 20: Reporting and Reproducibility: The Complete Teaching
Analysis
================

# Could another researcher reproduce your conclusion?

A figure alone does not reveal which people, weights, alleles or model
choices produced it. A reproducible report connects the question to the
inputs, processing, evaluation and limitations. Reproducibility does not
prove validity, but it allows errors and assumptions to be examined.

**Objectives:** run the complete fixed-score teaching analysis from a
fresh session, identify its outputs, write a restrained interpretation
and distinguish completed checks from unresolved scientific questions.
Prerequisites: Chapters 6–10 and 15–19.

# 1. The final project question

In the fictional, fully observed five-year CAD cohort, does adding the
fixed synthetic score to age and SBP improve prediction in independent
test participants?

The model is developed in 2,000 participants and evaluated in 2,000
separate participants. The other simulation groups demonstrate separate
discovery, tuning and LD roles in method chapters. This final script
does not select the best method across those chapters. It evaluates one
prespecified fixed score.

All genotypes, weights and outcomes are synthetic. There are no real
patients, no clinical recommendations and no evidence of ancestry
portability. The outcome-generating model includes a genetic
contribution, so this exercise demonstrates correct workflow rather than
testing whether genetics matters in real CAD.

# 2. Run from a fresh session

The supplied `run_analysis.R` sources the helper, reconstructs the fixed
synthetic data, fits the two models, calculates metrics and uncertainty,
and writes a new output directory. It does not depend on previously
knitted chapters.

Run from the course root in the Terminal:

``` bash
Rscript scripts/run_analysis.R
```

R is required. Base R suffices for the analysis script; knitting
chapters additionally requires rmarkdown, knitr and Pandoc. The script
refuses to overwrite its default completed run. If intentionally
rerunning, supply a new output directory:

``` bash
Rscript scripts/run_analysis.R results/fixed_score_run_02
```

The chapter itself only performs an in-memory check below; it does not
run the export script during knitting.

``` r
bundle <- fit_course_models()
knitr::kable(data.frame(Parameter=names(coef(bundle$extended)),
                        Estimate=unname(coef(bundle$extended))),digits=4)
```

| Parameter   | Estimate |
|:------------|---------:|
| (Intercept) |  -2.4889 |
| age10       |   0.5506 |
| sbp10       |   0.1424 |
| PRS_z       |   0.4514 |

``` r
stopifnot(length(intersect(bundle$train$IID,bundle$test$IID))==0)
```

# 3. Know what each output means

| Output in the chosen run directory | Content | What it establishes |
|----|----|----|
| `model_bundle.rds` | Fitted models and reference transformation | Model implementation can be restored |
| `test_predictions.tsv` | Synthetic IDs, outcomes and fixed predictions | Evaluation inputs are explicit |
| `metrics.tsv` | AUC and Brier for each model | Point estimates on the same test sample |
| `delta_auc.tsv` | Paired bootstrap interval and failed resample count | Conditional evaluation uncertainty |
| `association.tsv` | Test association OR per development SD | Adjusted association, not causality |
| `calibration.tsv` | Defined calibration estimates | Probability agreement diagnostics |
| `decision_curve.tsv` | Net benefit over illustrative thresholds | Decision calculation under stated assumptions |
| `audit.tsv` | Automated checks, assumptions and unassessed items | Scope of verification |
| `session_info.txt` | R/platform details | Computational context |
| `input_checksums.tsv` | Helper and runner fingerprints | Exact script provenance |
| `report.txt` | Computed teaching interpretation | Numerical summary without clinical claims |

An `.rds` file preserves R objects. It is convenient for this course,
but a deployed model also needs explicit variable definitions,
preprocessing and software dependencies. A CSV or TSV alone does not
preserve those modelling rules.

# 4. Report design before the headline result

State the phenotype, horizon, eligible population, score identity and
data roles. Explain matching, QC, missing-data handling, relatedness and
discovery overlap. Give the model formula and transformation reference.
Report event counts and uncertainty, not just a P value or best AUC. The
PRS reporting standards developed by Wand and colleagues provide a
structured basis for these disclosures. \[1\]

A results paragraph should say that the analysis used synthetic data,
which models were compared, what the paired performance difference was,
and what remains unassessed. The supplied script constructs this
paragraph from computed results. It does not prefill favorable numerical
conclusions.

Do not interpret a nonsignificant difference as proof of equivalence. Do
not interpret a positive AUC difference as proof that adding PRS is
worth genotyping, changes treatment appropriately or prevents disease.

# 5. Three levels of evidence

| Level | Course example |
|----|----|
| Automatically checked when executed | Disjoint participant IDs; finite predictions; defined binary outcome; successful model fitting |
| Supplied by the simulation design | No censoring, competing death, ancestry structure or relatedness |
| Not assessed | Clinical effectiveness, real-world portability, real genotype QC and genotype consent/access |

A program cannot prove the source’s clinical units simply because the
values are numeric. A clean report must preserve that distinction.

# 6. Publication and GitHub

Keep source Rmd, shared scripts, small shareable teaching data,
references and verification records together. Knit using the supplied
render script. Commit the generated Markdown and figure paths needed by
your chosen GitHub rendering route. Rmd files do not execute
automatically merely because they are pushed to GitHub.

Avoid uploading participant-level clinical or genetic data unless their
sharing terms explicitly permit it. Our synthetic outputs can support
teaching, but real workflows should keep restricted locations and
credentials out of public configuration. A `.gitignore` is useful but is
not a substitute for checking what was already committed.

Published scoring resources such as the PGS Catalog make model
identification and reuse easier. A reusable scoring file still does not
provide all information required to reproduce a new clinical validation.
Retain both model provenance and study design. \[2\]

# 7. What counts as course completion?

You should now be able to distinguish scoring from prediction modelling,
identify data roles, audit participant and variant inputs, explain major
method families, prevent outcome-driven leakage, evaluate several
performance domains and report uncertainty honestly.

Further study can cover individual-level penalized models, more advanced
survival prediction, local ancestry, multi-trait methods, formal
sample-size simulation and prospective implementation. These are
extensions, not reasons to overclaim what this beginner project
established.

# Practice

**Does reproducing the same number prove the model is scientifically
valid?** No; reproducible mistakes remain mistakes.

**Which output establishes a clinical benefit?** None of these synthetic
outputs.

**What should accompany a reported PRS OR?** Its reference scale,
outcome definition, covariates, sample and event counts, uncertainty and
study limitations.

**Why save both scripts and session information?** Code specifies the
procedure; software versions help explain implementation differences.

# References

1.  Wand H, Lambert SA, Tamburro C, et al. Improving reporting standards
    for polygenic scores in risk prediction studies. *Nature*.
    2021;591:211–219. <https://doi.org/10.1038/s41586-021-03243-6>
2.  Lambert SA, Gil L, Jupp S, et al. The Polygenic Score Catalog as an
    open database for reproducibility and systematic evaluation. *Nature
    Genetics*. 2021;53:420–425.
    <https://doi.org/10.1038/s41588-021-00783-5>
