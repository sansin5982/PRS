# Run from course root: Rscript scripts/render_course.R
stopifnot(file.exists("README.md"),dir.exists("chapters"))
if(!requireNamespace("rmarkdown",quietly=TRUE) || !requireNamespace("knitr",quietly=TRUE))
 stop("Install rmarkdown and knitr first")
if(!rmarkdown::pandoc_available()) stop("Pandoc is required; RStudio normally supplies it")
files <- c("index.Rmd",sort(list.files("chapters",pattern="\\.Rmd$",full.names=TRUE)))
for(f in files) {
 message("Rendering ",f)
 rmarkdown::render(f,envir=new.env(parent=globalenv()),quiet=FALSE)
}
