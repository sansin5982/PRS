Chapter 15: Association Analysis: Does PRS Track the Outcome?
================

# An association is a useful result, but not the final answer

We have calculated a fixed score. The next question is whether higher
scores are associated with more CAD events in the evaluated sample. This
is distinct from whether probabilities are accurate or whether using the
score improves healthcare.

**Objectives:** fit an adjusted association model, interpret an odds
ratio per specified SD, report uncertainty, and distinguish association
from prediction and causation. Prerequisites: Chapters 5, 6 and 10. Base
R and the course helper are sufficient.

# 1. State the comparison before fitting

We analyze the independent test participants with the fixed synthetic
score, adjusting for age and SBP. The score is standardized using
development parameters. We estimate an association in these test data;
this is not training a prediction model and then claiming its in-sample
performance is independent.

The logistic model relates predictors to **log odds**. Exponentiating a
PRS coefficient gives the odds ratio for a one-unit increase in that
predictor, holding the other included predictors fixed. Here the unit is
one development-sample SD. The coefficient depends on model
specification; it is not a biological constant.

``` text
Odds = probability / (1 - probability)
OR per reference SD = exp(PRS coefficient)
```

An OR of 1.5 does not mean a 50 percentage-point increase in risk. The
same odds ratio corresponds to different probability differences at
different baseline risks. A significance test asks about evidence under
a statistical model, not whether the effect is clinically important.

# 2. Fit the association model

``` r
bundle <- fit_course_models()
d <- bundle$test
association <- glm(CAD5~PRS_z+age10+sbp10,data=d,family=binomial())
stopifnot(association$converged,all(is.finite(coef(association))))
tab <- coef(summary(association))
b <- tab["PRS_z","Estimate"]; se <- tab["PRS_z","Std. Error"]
association_result <- data.frame(
  Comparison="One development SD increase in fixed synthetic PRS",
  OR=exp(b),Lower=exp(b-1.96*se),Upper=exp(b+1.96*se),
  P=tab["PRS_z","Pr(>|z|)"],N=nrow(d),Cases=sum(d$CAD5))
knitr::kable(association_result,digits=4)
```

| Comparison | OR | Lower | Upper | P | N | Cases |
|:---|---:|---:|---:|---:|---:|---:|
| One development SD increase in fixed synthetic PRS | 1.5309 | 1.3114 | 1.7871 | 0 | 2000 | 179 |

The interval is a large-sample Wald interval. Sparse events, separation
and misspecification can make this approximation unsuitable. The
convergence check is necessary but not a full diagnostic. For real data,
inspect fitting warnings, functional form, influential observations and
sampling structure.

The simulation was constructed with a genetic contribution, so an
association is expected on average. The precise output is generated when
you run the chapter; no unexecuted result is presented as an observed
finding.

# 3. Adjustment has a purpose

Age and SBP are included because they are part of the prespecified
teaching question. The baseline simulation has no population structure,
relatedness or batch effect. Their absence is a declared assumption, not
something this regression demonstrates.

For real genetic association, ancestry PCs and technical factors may be
relevant. Avoid copying a fixed number of PCs without investigating the
dataset. Relatedness may require mixed models or appropriate
dependence-aware inference. A covariate that occurs after disease onset
may be unsuitable for a prospective prediction question even if it
correlates with disease.

Adding covariates can change an odds ratio even without conventional
confounding because odds ratios are non-collapsible. Consequently, a
difference between adjusted and unadjusted ORs should not automatically
be interpreted as “the percentage explained by confounding.” \[1\]

# 4. Why not compare top versus bottom groups only?

Extreme-group comparisons can be intuitive, but lose information and
depend on chosen cutoffs. If you report them, define the reference
group, percentile source, threshold selection and uncertainty. The top
5% versus everyone else is different from top 5% versus bottom 5%.

``` r
# Descriptive categories using development cutoffs, not test-selected cutoffs.
cuts <- quantile(bundle$train$PRS_z,c(.2,.4,.6,.8),names=FALSE)
d$group <- cut(d$PRS_z,breaks=c(-Inf,cuts,Inf),labels=paste0("Q",1:5))
group_summary <- do.call(rbind,lapply(split(d,d$group),function(x)
  data.frame(N=nrow(x),Cases=sum(x$CAD5),Observed_fraction=mean(x$CAD5))))
knitr::kable(group_summary,digits=3)
```

|     |   N | Cases | Observed_fraction |
|:----|----:|------:|------------------:|
| Q1  | 433 |    19 |             0.044 |
| Q2  | 410 |    26 |             0.063 |
| Q3  | 425 |    41 |             0.096 |
| Q4  | 350 |    39 |             0.111 |
| Q5  | 382 |    54 |             0.141 |

These observed group fractions are descriptive, with finite-sample
uncertainty. They are not personalized risk estimates. Test group sizes
need not be exactly equal because cutoffs came from development data.

# 5. Published example and causal restraint

Khera et al. reported CAD association and risk stratification for a
genome-wide score in their study population. Such results motivate
evaluation elsewhere but do not supply the coefficient or baseline
probability for every new cohort. \[2\]

A PRS association can reflect tagged biological effects, population
structure, indirect genetic influences or selection processes. It does
not satisfy Mendelian-randomization assumptions by itself. Do not claim
that changing a person’s PRS would cause the estimated change in disease
risk.

# Practice and handover

**What is the unit of our OR?** One development-sample SD of the fixed
synthetic PRS.

**Is a small P value enough to establish useful prediction?** No. Sample
size can make a modest association precise without making it useful for
decisions.

**Are group fractions clinical probabilities?** They describe the
evaluated group and sampling design; individual prediction needs a
suitable model.

Save an association table with N, event count, reference SD, covariates,
OR, interval and model limitations. Chapter 16 evaluates predictions
from models fitted entirely outside the test sample.

# References

1.  Greenland S, Robins JM, Pearl J. Confounding and collapsibility in
    causal inference. *Statistical Science*. 1999;14:29–46.
    <https://doi.org/10.1214/ss/1009211805>
2.  Khera AV, Chaffin M, Aragam KG, et al. Genome-wide polygenic scores
    for common diseases identify individuals with risk equivalent to
    monogenic mutations. *Nature Genetics*. 2018;50:1219–1224.
    <https://doi.org/10.1038/s41588-018-0183-z>
