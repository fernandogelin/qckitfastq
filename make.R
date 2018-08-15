library("Rd2md")
library("knitr")

path = "man/"
files = dir(path)

for(file in files){
  infile = file.path(path, file)
  outfile = file.path("docs", gsub(".Rd", ".md", file))
  Rd2markdown(infile, outfile, append = FALSE)
}

vignette_file = dir("vignettes/")
in_ = file.path("vignettes", vignette_file)
out_ = file.path("docs", gsub(".Rmd", ".md", vignette_file))
knit(in_, out_)


file_rename <- function(from, to) {
    todir <- dirname(to)
    if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
    file.rename(from = from,  to = to)
}

figures = dir("figure/")
for (f in figures){
  file_rename(from = file.path("figure", f),
              to = file.path("docs", "figure", f))
}

file_rename(from = file.path("README.md"),
            to = file.path("docs", "README.md"))
