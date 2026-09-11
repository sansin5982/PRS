# Run at course root. Produces SYNTHETIC hard-call input, not real genomic data.
stopifnot(file.exists("scripts/course_helpers.R"))
source("scripts/course_helpers.R")
out <- "results/plink_demo"
if(dir.exists(out)) stop("Refusing to overwrite an existing PLINK demo directory")
dir.create(out,recursive=TRUE)
demo <- make_course_data(); ix <- which(demo$people$role=="Test")[1:50]
G <- demo$G[ix,,drop=FALSE]; w <- demo$weights
ped <- t(vapply(seq_along(ix),function(k) {
 g <- G[k,]
 alleles <- as.vector(rbind(ifelse(g>=1,"A","G"),ifelse(g==2,"A","G")))
 c("0",demo$people$IID[ix[k]],"0","0","1","-9",alleles)
},character(6+2*ncol(G))))
write.table(ped,file.path(out,"target.ped"),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(data.frame(CHR=w$CHR,ID=w$ID,CM=0,POS=w$POS),
 file.path(out,"target.map"),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(w[c("ID","A1","weight")],file.path(out,"weights.tsv"),
 quote=FALSE,row.names=FALSE,sep="\t")
write.table(data.frame(IID=demo$people$IID[ix],R_sum=as.numeric(G%*%w$weight)),
 file.path(out,"expected_scores.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
message("Wrote 50 fictional participants and 120 fictional loci in ",out)
