% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AStrap.R
\name{AStrap}
\alias{AStrap}
\title{Identification AS events and classification AS types.}
\usage{
AStrap(cluster.alignment, transcriptSeq, trainModel, identity = 0.7,
  coverage = 0.7, bias = 0, ASlength = 0)
}
\arguments{
\item{cluster.alignment}{A data frame holds pairwise alignment of isoforms in a cluster.}

\item{transcriptSeq}{An XStringSet object holds the transcript sequence.}

\item{trainModel}{A model object for which prediction is desired.}

\item{identity}{AS detection is performed if identity above a given threshold (default: 0.7).}

\item{coverage}{AS detection is performed if coverage above a given threshold (default: 0.7).}

\item{bias}{Maximum number of mismatches in pairwise sequence alignment (default: 0).}

\item{ASlength}{AS detection is performed if AS length above a given threshold (default: 0).}
}
\value{
A list of AS events and AS types. This object contains three elements:
 ASevent, feature, predict. ASevent is a data frame stroage the identified AS events;
 feature is a feature data frame for AS events based on transcript, including 511 features;
 predict is the calssification of AS events. and the form of the value returned by predict
 depends on the class of its argument.
}
\description{
This function can identify AS events and classify AS events from
transcripts without a reference.Two classification models trained on collected
AS data from rice (MSU7) and human (GRCH38) were integrated in AStrap, which
could be directly applied for distinguishing AS types for other species.
Certainly, users can also train a specific model on their own data sets
using several functions provided in AStrap, see the function of buildTrainModel.
}
\examples{
##Loading transcript sequence
trSequence.path <- system.file("extdata","example_TRsequence.fasta",
                              package = "AStrap")
trSequence <-  readDNAStringSet(trSequence.path,format = "fasta")

##Loadinding example data output by CD-HIT-EST
cdhit.path <- system.file("extdata","example_cdhitest.clstr",
                         package = "AStrap")
raw.cluster <- readCDHIT(cdhit.path)

##Loading example data output by GMAP
gmap.path <- system.file("extdata","example_gmap.gff3",package = "AStrap")

##getting pairwise alignment of isoforms of the same cluster
cluster.align <- readGMAP(gmap.path,raw.cluster, recluster = TRUE,
                          recluster.identity = 0.7,
                          recluster.coverage = 0.7)
alignment <- cluster.align$alignment
rew.cluster <- cluster.align$cluster

##Loading rice model
load(system.file("data","rice_model.Rdata",package = "AStrap"))

##Identification and prediction based on rice RF-based model
result <- AStrap(alignment,trSequence,rice_RFmodel)
}
