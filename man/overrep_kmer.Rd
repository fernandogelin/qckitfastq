% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/over_rep_kmer.R
\name{overrep_kmer}
\alias{overrep_kmer}
\title{Generate overrepresented kmers from all kmer counts results}
\usage{
overrep_kmer(path, k, nc, nr, writefile = FALSE, prefix)
}
\arguments{
\item{path}{the path to the gz file}

\item{k}{the length of the sequence looking for}

\item{nc}{number of positions}

\item{nr}{number of reads}

\item{writefile}{TRUE for writing the output to csv file}

\item{prefix}{prefix of output file if writefile is TRUE}
}
\value{
the index of reads that has overrepresented kmers
}
\description{
Generate overrepresented kmers from all kmer counts results
}
\examples{

path <-system.file("extdata", "10^5_reads_test.fq.gz", package = "qckitfastq")
overrep_kmer(path,k=4,nc=100,nr=25000)

}
