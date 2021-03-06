% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{qual_score_per_read}
\alias{qual_score_per_read}
\title{calculate mean quality per read}
\usage{
qual_score_per_read(infile)
}
\arguments{
\item{infile}{A string giving the path for the fastqfile}
}
\value{
mean quality per read
}
\description{
Calculate the mean quality score per read of the FASTQ gzipped file
}
\examples{
infile <- system.file("extdata", "10^5_reads_test.fq.gz", package = "qckitfastq")
qual_score_per_read(infile)
}
