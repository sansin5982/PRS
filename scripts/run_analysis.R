# Run at the course root: Rscript scripts/run_analysis.R [new_output_directory]
stopifnot(file.exists("scripts/course_helpers.R"))
source("scripts/course_helpers.R")
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else "results/fixed_score_run_01"
if(dir.exists(out)) stop("Output directory already exists; choose a new run directory")
dir.create(out,recursive=TRUE)
write_tsv <- function(x,name) write.table(x,file.path(out,name),sep="\t",
                                        row.names=FALSE,quote=FALSE,na="NA")
demo <- make_course_data(); fit <- fit_course_models(demo); d <- fit$test
stopifnot(!anyDuplicated(d$IID),length(intersect(fit$train$IID,d$IID))==0,
          all(is.finite(d$p_base)),all(is.finite(d$p_prs)),
          all(d$p_base>0 & d$p_base<1),all(d$p_prs>0 & d$p_prs<1))
metrics <- data.frame(Model=c("Age_SBP","Age_SBP_PRS"),
 AUC=c(auc(d$CAD5,d$p_base),auc(d$CAD5,d$p_prs)),
 Brier=c(mean((d$CAD5-d$p_base)^2),mean((d$CAD5-d$p_prs)^2)))
b <- paired_auc_bootstrap(d$CAD5,d$p_base,d$p_prs,B=2000)
delta <- data.frame(Delta_AUC=b$delta,Lower=b$CI[1],Upper=b$CI[2],
                    Failed=b$failed,Resamples=2000)
assoc <- glm(CAD5~PRS_z+age10+sbp10,data=d,family=binomial())
a <- coef(summary(assoc))["PRS_z",]
association <- data.frame(OR=exp(a[1]),Lower=exp(a[1]-1.96*a[2]),
 Upper=exp(a[1]+1.96*a[2]),P=a[4],N=nrow(d),Cases=sum(d$CAD5))
lp <- qlogis(d$p_prs)
c0 <- glm(d$CAD5~1+offset(lp),family=binomial())
c1 <- glm(d$CAD5~lp,family=binomial())
stopifnot(assoc$converged,c0$converged,c1$converged,
 all(is.finite(coef(assoc))),all(is.finite(coef(c0))),all(is.finite(coef(c1))))
cal <- data.frame(Measure=c("Intercept_slope_fixed_1","Joint_intercept","Joint_slope"),
                  Estimate=c(coef(c0)[1],coef(c1)[1],coef(c1)[2]))
t <- seq(.05,.25,.01)
decision <- data.frame(Threshold=t,
 Baseline=vapply(t,function(x) net_benefit(d$CAD5,d$p_base,x),numeric(1)),
 Extended=vapply(t,function(x) net_benefit(d$CAD5,d$p_prs,x),numeric(1)),
 Act_all=mean(d$CAD5)-(1-mean(d$CAD5))*t/(1-t),Act_none=0)
audit <- data.frame(Item=c("Disjoint development/test","Finite predictions",
 "Binary outcome","Unrelated / complete follow-up","Real-world clinical utility"),
 Status=c("CHECKED","CHECKED","CHECKED","SIMULATION ASSUMPTION","NOT ASSESSED"))
saveRDS(list(base=fit$base,extended=fit$extended,reference=fit$reference,
             simulation=demo$metadata),file.path(out,"model_bundle.rds"))
write_tsv(d,"test_predictions.tsv");write_tsv(metrics,"metrics.tsv")
write_tsv(delta,"delta_auc.tsv");write_tsv(association,"association.tsv")
write_tsv(cal,"calibration.tsv");write_tsv(decision,"decision_curve.tsv")
write_tsv(audit,"audit.tsv")
writeLines(capture.output(sessionInfo()),file.path(out,"session_info.txt"))
paths <- c("scripts/course_helpers.R","scripts/run_analysis.R")
write_tsv(data.frame(File=paths,MD5=unname(tools::md5sum(paths))),"input_checksums.tsv")
report <- c("SYNTHETIC TEACHING ANALYSIS: no real patients or clinical evidence.",
 sprintf("Test sample: %d people, %d events. Development sample: %d people.",
         nrow(d),sum(d$CAD5),nrow(fit$train)),
 sprintf("AUC baseline %.4f; extended %.4f. Paired difference %.4f (95%% percentile interval %.4f to %.4f).",
 metrics$AUC[1],metrics$AUC[2],b$delta,b$CI[1],b$CI[2]),
 "Intervals are conditional on fitted models; calibration estimates are reported separately.",
 "No claim of clinical utility, causality, real genotype validity or ancestry portability is supported.")
writeLines(report,file.path(out,"report.txt"))
message("Completed synthetic analysis in ",out)
