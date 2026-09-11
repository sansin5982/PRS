# Course plan

Twenty separate chapters; core lessons precede optional advanced methods and interpretation. Chapters 1–7 retain reviewed foundational explanations with corrected course conventions. Chapters 8–20 use the shared project or explicitly labelled methodological mini-exercises.

## 1. What Is a Polygenic Risk Score

Purpose: Understand what a score summarizes. Prerequisites: None.

Concepts: Susceptibility; weighted sum; limits. Practical: Hand calculation and reference plot.

Output: Arithmetic and interpretation. Handover: Genetic foundations. Level: Core.

## 2. Genetic Foundations for PRS

Purpose: Read genetic representations. Prerequisites: 1.

Concepts: Alleles; inheritance; LD; dosage. Practical: Genotype and dosage exercises.

Output: Checked definitions. Handover: GWAS and harmonization. Level: Core.

## 3. From GWAS to PRS

Purpose: Understand where weights originate. Prerequisites: 1–2.

Concepts: GWAS; effect scale; uncertainty. Practical: Inspect example summary statistics.

Output: GWAS schema. Handover: Calculation and methods. Level: Core.

## 4. How Is a PRS Calculated

Purpose: Calculate a score correctly. Prerequisites: 1–3.

Concepts: Weighted sums; matching; scaling. Practical: R matrix scoring.

Output: Toy individual scores. Handover: Interpretation and fixed scoring. Level: Core.

## 5. Interpreting a PRS

Purpose: Interpret a score without overclaiming. Prerequisites: 1–4.

Concepts: Percentiles; odds; probability; performance. Practical: Reference and evaluation mini-exercises.

Output: Interpretation examples. Handover: Planning and evaluation. Level: Core.

## 6. Planning a PRS Analysis

Purpose: Specify a defensible study. Prerequisites: 1–5.

Concepts: Data roles; independence; precision. Practical: Planning and leakage demonstrations.

Output: Teaching protocol. Handover: Project preparation. Level: Core.

## 7. Preparing the PRS Project and Input Files

Purpose: Audit inputs and identities. Prerequisites: 2,4,6.

Concepts: File formats; keys; units; provenance. Practical: Input checks and joins.

Output: Readiness report. Handover: QC and harmonization. Level: Core.

## 8. Genotype Quality Control

Purpose: Investigate measurement quality. Prerequisites: 2,6–7.

Concepts: Missingness; HWE; relatedness; batch. Practical: Synthetic missingness audit.

Output: QC count table and flags. Handover: Harmonization. Level: Core.

## 9. Genome Builds and Allele Harmonization

Purpose: Align the same variants and alleles. Prerequisites: 2,4,7–8.

Concepts: Build; strand; swap; coverage. Practical: Exact/swap/unresolved matching.

Output: Alignment report. Handover: Scoring. Level: Core.

## 10. Applying a Fixed Score

Purpose: Apply fixed weights reproducibly. Prerequisites: 4,7–9.

Concepts: Sum vs average; reference scale. Practical: R scoring; optional PLINK exercise.

Output: Scores and reference parameters. Handover: Methods and evaluation. Level: Core.

## 11. Clumping and Thresholding

Purpose: Learn C+T selection. Prerequisites: 3,6,8–10.

Concepts: Clumping; thresholds; tuning. Practical: Discovery GWAS and greedy clumping.

Output: Candidate and selected weights. Handover: LD-aware methods. Level: Core.

## 12. LDpred2

Purpose: Learn LDpred2 inputs and assumptions. Prerequisites: 3,9,11.

Concepts: Priors; LD; inf/grid/auto. Practical: Input preparation; optional software fit.

Output: Input audit and candidate weights. Handover: Alternative methods. Level: Advanced method.

## 13. PRS CS

Purpose: Learn continuous shrinkage. Prerequisites: 3,9,12.

Concepts: Global/local shrinkage; phi; auto. Practical: Schema checks; external test protocol.

Output: PRS-CS input schema. Handover: Multi-ancestry methods. Level: Advanced method.

## 14. Multi Ancestry Methods

Purpose: Choose methods across populations. Prerequisites: 11–13.

Concepts: Shared evidence; score combinations. Practical: Synthetic combination exercise.

Output: Method-choice comparison. Handover: Portability. Level: Advanced method.

## 15. Association Analysis

Purpose: Estimate an adjusted association. Prerequisites: 5,6,10.

Concepts: OR; covariates; reference SD. Practical: Logistic association.

Output: Association table. Handover: Prediction comparison. Level: Core.

## 16. Independent Prediction Evaluation

Purpose: Evaluate incremental prediction. Prerequisites: 5,6,10,15.

Concepts: AUC; Brier; paired uncertainty. Practical: Independent prediction and bootstrap.

Output: Metrics and interval. Handover: Calibration. Level: Core.

## 17. Calibration and Absolute Risk

Purpose: Assess probability agreement. Prerequisites: 5,16.

Concepts: Calibration; updating; absolute risk. Practical: Calibration shift and diagnostic plots.

Output: Calibration table. Handover: Decision analysis. Level: Core.

## 18. Clinical Utility

Purpose: Connect predictions to decisions. Prerequisites: 16–17.

Concepts: Thresholds; PPV; net benefit. Practical: Decision curves and uncertainty.

Output: Threshold and net-benefit tables. Handover: Sensitivity analysis. Level: Advanced interpretation.

## 19. Portability and Sensitivity Analysis

Purpose: Investigate robustness and transfer. Prerequisites: 6,14,16–18.

Concepts: Portability; measurement; subgroup precision. Practical: Measurement-noise stress test.

Output: Sensitivity comparison. Handover: Reporting. Level: Advanced interpretation.

## 20. Reporting and Reproducibility

Purpose: Complete a reproducible analysis. Prerequisites: 6–10,15–19.

Concepts: Provenance; outputs; claims. Practical: Fresh-session runner.

Output: Versioned run outputs and report. Handover: Independent application. Level: Core.

