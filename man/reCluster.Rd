% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AStrap.R
\name{reCluster}
\alias{reCluster}
\title{Adjusting the cluster.}
\usage{
reCluster(Align, cdhit.result, re.id = 0.7, re.co = 0.7)
}
\arguments{
\item{Align}{A data frame contation the information of sequence alignments.}

\item{cdhit.result}{A data frame holds isoform clusters.}

\item{re.id}{An numeric value indicating the
threshold of identity of alignment for the second-round
mapping (default: 0.7). Used only if recluster = TRUE.}

\item{re.co}{An numeric value indicating the
threshold of coverage of alignment for the second-round
mapping (default: 0.7). Used only if recluster = TRUE.}
}
\value{
A list contation the pairwise sequence alignments
 in a cluster, including two elements: alignment and cluster.
 alignment is a data frame storage the pairwise alignment
 of isoform in a cluster; cluster is a data frame storage
 the new clustering result after second-round mapping.
}
\description{
To obtain all potential clusters with more
  than one isoform for further analysis, a second-round mapping
  can be performed on clusters with single isoform. Isoforms passing
  the second-round mapping were reassigned to the .
}
\examples{
##Loading the example data
load(system.file("data","example_recluster.Rdata",package = "AStrap"))
result <- reCluster(gmap.align,cdhit.result,re.id=0.7,re.co=0.7)
}
