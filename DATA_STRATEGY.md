# Dataset strategy

## Main recurring project: explicitly synthetic five-year CAD

Question: does a fixed synthetic PRS improve prediction beyond age and systolic blood pressure in independent test people?

The generator creates 120 biallelic autosomal variants and 11,000 unrelated people: 5,000 discovery; 2,000 development; 1,000 tuning; 2,000 test; 1,000 LD reference. LD-reference outcomes are omitted. Alleles are A/G with genotypes counting A. Positions and IDs are fictional, labelled SIMULATION_ONLY. There are 20 local blocks; this is not human genome-wide LD.

Age is in years and SBP in mmHg. CAD5 is 0/1 with complete five-year observation, no competing death, no censoring, no population structure and no technical batch effects. These are assumptions supplied by the generator, not empirical QC findings. The fixed score is noisy relative to the disclosed generating effects, but it is not a fitted published score.

The integer-like genotype matrix has 1.32 million entries; even as R doubles its raw payload is approximately 10.6 MB, before object copies and other memory use. This is a teaching-scale workload, not a genome-wide runtime estimate. Recreating the helper needs no network access or permission-controlled data. Exact random draws depend on R's RNG implementation; record session information.

## Public real-data resources

- [1000 Genomes public data](https://www.internationalgenome.org/data/) provide genotype/reference resources and population metadata. Public downloads do not require the clinical-cohort application process, but consult the resource's use conditions. Whole-genome files can be large; inspect release sizes and use a documented subset. These files alone do not supply the prospective CAD outcome needed by our project.
- [PGS Catalog PGS000013](https://www.pgscatalog.org/score/PGS000013/) supplies a real CAD score record and linked scoring downloads. Its stated terms distinguish academic research and commercial use. The original record lists hg19; choose any harmonized build deliberately and record the exact file. Genotype scoring is possible with compatible inputs, but performance validation needs appropriate outcomes.
- [Official LDpred2 educational data](https://privefl.github.io/bigsnpr/articles/LDpred2.html) are presented as fake educational data. Do not relabel their phenotypes clinical observations.
- [PRS-CS](https://github.com/getian107/PRScs) and [PRS-CSx](https://github.com/getian107/PRScsx) document public test inputs and LD resources. Their test examples verify a software workflow, not our synthetic or any new clinical cohort. Reference downloads may dominate storage requirements.

The course uses these real resources for inspection and optional documented software exercises. It does not bundle them, claim successful downloads, or invent clinical phenotypes for their real participants. Restricted cohorts require their own authorization and do not block the self-contained teaching route.

## Interpretation boundary

Separate arithmetic demonstration, software smoke test, association, prediction validation and clinical implementation. No synthetic result in this package establishes a clinical benefit or ancestry portability.
