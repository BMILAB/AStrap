#' @title Read the clusters from a file.
#'
#' @description  Read a text file of list of clusters output
#'   by CD-HIT-EST and creates a data frame form it.
#'
#' @param cdhit.path A string for the name of a text file of
#'   list of clusters output by CD-HIT-EST.
#' @return A data frame contation the information of isoform
#'   clusters with four columns: ClusterID, seqID, seqLen, seqNum.
#'   ClusterID is the name of the cluster; seqID is the name of the
#'   transcript; seqLen is the length of the transcript; seqNum is
#'   the number of transcripts contatined in this class.
#' @examples
#' ##load the exapmle data
#' cdhit.path <- system.file("extdata","example_cdhitest.clstr",
#'                            package = "AStrap")
#' raw.cluster <- readCDHIT(cdhit.path)
readCDHIT <- function(cdhit.path) {
  options(stringsAsFactors = FALSE)
  cdhit <- read.table(cdhit.path, header = FALSE, sep = "\n", quote = NULL,
                      stringsAsFactors = FALSE)
  cluster.idx <- which(substr(cdhit$V1, 1, 1) == ">")
  num.reads <- nrow(cdhit) - length(cluster.idx)
  clusters <- data.frame(cluster = rep("", num.reads), id = rep("", num.reads),
                         len = rep(0, num.reads), stringsAsFactors = FALSE)
  count <- 0
  for (i in 1:length(cluster.idx)) {
    if (i == length(cluster.idx)) {
      end.idx <- nrow(cdhit)
    } else {
      end.idx <- cluster.idx[i + 1] - 1
    }
    for (j in (cluster.idx[i] + 1):end.idx) {
      count <- count + 1
      len <- str_extract(cdhit$V1[j], "(?<=\\t)\\d+(?=nt||aa)")
      id <- str_extract(cdhit$V1[j], "(?<=\\>)\\S+?(?=\\.\\.\\.)")
      clusters$cluster[count] <- i
      clusters$id[count] <- id
      clusters$len[count] <- len
    }
  }
  cluster.sumary <- as.data.frame(table(clusters$cluster))
  colnames(cluster.sumary) <- c("cluster", "number")
  cluster.sumary$cluster <- as.character(cluster.sumary$cluster)
  Cluster <- merge(clusters, cluster.sumary, by = "cluster")

  Cluster.order <- Cluster[order(as.numeric(Cluster$cluster)), ]
  colnames(Cluster.order) <- c("ClusterID", "seqID", "seqLen", "seqNum")
  Cluster.order$ClusterID <- as.numeric(Cluster.order$ClusterID)
  Cluster.order$seqLen <- as.numeric(Cluster.order$seqLen)
  row.names(Cluster.order)<-c(1:nrow(Cluster.order))
  return(Cluster.order)
}



#' @title Alignment trimming
#'
#' @description This function is to treat redundant sequence alignment
#'   and remove the invalid sequence alignments.
#'
#' @param gmap.file A string for the name of a GFF3 file of pairwise
#'   sequence alignments output by GMAP.
#' @return A data frame holds trimmed sequence alignment.
#' @examples
#' ##Load example data
#' gmap.path <- system.file("extdata","example_gmap.gff3",package = "AStrap")
#' trimmed.ali <- trimmGMAP(gmap.path)
trimmGMAP <- function(gmap.file) {
  options(stringsAsFactors = F)
  gmap <- import.gff(gmap.file, format = "gff3",
                     feature.type = c("mRNA", "exon"))
  gmap <- as.data.frame(gmap)
  gmap$seqnames <- as.character(gmap$seqnames)
  gmap$type <- as.character(gmap$type)
  gmap$coverage <- as.numeric(gmap$coverage)
  gmap$identity <- as.numeric(gmap$identity)
  gmap$strand <- as.character(gmap$strand)
  gmap <- gmap[which(gmap$strand == "+"), ]
  gmap <- gmap[which(gmap$seqnames != gmap$Name), ]
  mRNA.id <- which(gmap$type == "mRNA")
  align <- gmap[mRNA.id, ]
  align$num <- 0
  align$Subject <- ""
  align$Query <- ""
  for (i in 1:length(mRNA.id)) {
    if (i == length(mRNA.id)) {
      end.idx <- nrow(gmap)
    } else {
      end.idx <- mRNA.id[i + 1] - 1
    }
    exon.num <- end.idx - mRNA.id[i]
    align$num[i] <- exon.num
    for (j in (mRNA.id[i] + 1):end.idx) {
      start <- str_extract(gmap$Target[j], "(?<=\\s)\\d+(?=\\s\\d+)")
      end <- str_extract(gmap$Target[j], "(?<=\\s)\\d+(?=\\s\\+)")
      Query.align <- paste(start, end, sep = "-")
      align$Query[i] <- paste(align$Query[i], Query.align, sep = ":")

      Subject.align <- paste(gmap$start[j], gmap$end[j], sep = "-")
      align$Subject[i] <- paste(align$Subject[i], Subject.align, sep = ":")
    }
  }
  gmap.align <- data.frame(Sid = as.character(align$seqnames), Qid = align$Name,
                           Coverage = align$coverage, identity = align$identity,
                           Alinum = align$num, Salign = align$Subject,
                           Qalign = align$Query)
  gmap.align$NAME <- ""
  index1 <- gmap.align$Sid > gmap.align$Qid
  gmap.align$NAME[index1] <- paste(gmap.align$Sid[index1],
                                   gmap.align$Qid[index1], sep = "")
  index2 <- gmap.align$Sid < gmap.align$Qid
  gmap.align$NAME[index2] <- paste(gmap.align$Qid[index2],
                                   gmap.align$Sid[index2], sep = "")

  gmap.rep <- as.data.frame(table(gmap.align$NAME))
  colnames(gmap.rep) <- c("id", "fre")
  gmap.rep$id <- as.character(gmap.rep$id)
  gmap.rep$fre <- as.integer(gmap.rep$fre)
  one.name <- as.character(gmap.rep$id[which(gmap.rep$fre == 1)])
  index.one <- which((gmap.align$NAME) %in% one.name)
  dup.name <- as.character(gmap.rep$id[which(gmap.rep$fre > 1)])
  index.need <- c(index.one)
  for (i in 1:length(dup.name)) {
    index.tem <- which(gmap.align$NAME == dup.name[i])
    k <- index.tem[1]
    for (j in 2:length(index.tem)) {
      m <- index.tem[j]
      if (gmap.align$Alinum[m] > gmap.align$Alinum[k]) {
        k <- m
      } else if (gmap.align$Alinum[m] == gmap.align$Alinum[k]) {
        if (gmap.align$Coverage[m] > gmap.align$Coverage[k]) {
          k <- m
        }
      }
    }
    index.need <- c(index.need, k)
  }
  gmap.result <- gmap.align[index.need, ]
  return(gmap.result)
}



#' @title Adjusting the cluster.
#'
#' @description To obtain all potential clusters with more
#'   than one isoform for further analysis, a second-round mapping
#'   can be performed on clusters with single isoform. Isoforms passing
#'   the second-round mapping were reassigned to the .
#'
#' @param Align A data frame contation the information of sequence alignments.
#' @param cdhit.result A data frame holds isoform clusters.
#' @param re.id An numeric value indicating the
#'   threshold of identity of alignment for the second-round
#'   mapping (default: 0.7). Used only if recluster = TRUE.
#' @param re.co An numeric value indicating the
#'   threshold of coverage of alignment for the second-round
#'   mapping (default: 0.7). Used only if recluster = TRUE.
#' @return A list contation the pairwise sequence alignments
#'  in a cluster, including two elements: alignment and cluster.
#'  alignment is a data frame storage the pairwise alignment
#'  of isoform in a cluster; cluster is a data frame storage
#'  the new clustering result after second-round mapping.
#' @examples
#' ##Loading the example data
#' load(system.file("data","example_recluster.Rdata",package = "AStrap"))
#' result <- reCluster(gmap.align,cdhit.result,re.id=0.7,re.co=0.7)
reCluster <- function(Align, cdhit.result, re.id = 0.7, re.co = 0.7) {
  indexSone <- which(Align$SseqNum == 1)
  indexQone <- which(Align$QseqNum == 1)
  OneName <- unique(c(Align$Sid[indexSone], Align$Qid[indexQone]))
  Cluster1 <- cdhit.result
  ClusterOne <- data.frame(seqId = OneName, coverage = rep(0, length(OneName)),
                           identity = rep(0, length(OneName)),
                           Target = rep(0, length(OneName)))
  for (i in 1:length(OneName)) {
    index.s <- which(Align$Sid == OneName[i])
    index.q <- which(Align$Qid == OneName[i])
    if (length(index.q) > 0) {
      for (j in 1:length(index.q)) {
        k <- index.q[j]
        if (Align$identity[k] > ClusterOne$identity[i]) {
          ClusterOne$identity[i] <- Align$identity[k]
          ClusterOne$coverage[i] <- Align$Coverage[k]
          ClusterOne$Target[i] <- as.character(Align$Sid[k])
        } else if (Align$identity[k] == ClusterOne$identity[i]) {
          if (Align$Coverage[k] > ClusterOne$coverage[i]) {
            ClusterOne$identity[i] <- Align$identity[k]
            ClusterOne$coverage[i] <- Align$Coverage[k]
            ClusterOne$Target[i] <- Align$Sid[k]
          }

        }
      }
    }

    if (length(index.s) > 0) {
      for (j in 1:length(index.s)) {
        k <- index.s[j]
        if (Align$identity[k] > ClusterOne$identity[i]) {
          ClusterOne$identity[i] <- Align$identity[k]
          ClusterOne$coverage[i] <- Align$Coverage[k]
          ClusterOne$Target[i] <- Align$Qid[k]
        } else if (Align$identity[k] == ClusterOne$identity[i]) {
          if (Align$Coverage[k] > ClusterOne$coverage[i]) {
            ClusterOne$identity[i] <- Align$identity[k]
            ClusterOne$coverage[i] <- Align$Coverage[k]
            ClusterOne$Target[i] <- Align$Qid[k]
          }

        }
      }
    }
  }

  index.c <- which(ClusterOne$coverage >= re.co * 100)
  index.i <- which(ClusterOne$identity >= re.id * 100)
  index.ci <- intersect(index.c, index.i)
  ClusterOne <- ClusterOne[index.ci, ]
  index.dell <- which(Cluster1$seqID %in% ClusterOne$seqId)
  Cluster1 <- Cluster1[-index.dell, ]
  max.clusterid <- max(as.numeric(cdhit.result$ClusterID))
  data <- data.frame(id1 = ClusterOne$seqId, id2 = ClusterOne$Target)
  g <- graph_from_data_frame(data, directed = FALSE)
  connect.igraph <- clusters(g)$membership
  connect.id <- unique(connect.igraph)
  for (i in 1:length(connect.id)) {
    re.name <- names(connect.igraph)[which(connect.igraph %in% connect.id[i])]
    index <- which(Cluster1$seqID %in% re.name)
    overlap.num <- length(index)
    re.num <- length(re.name)
    if (overlap.num > 0) {
      cluster.id <- as.numeric(unique(Cluster1$ClusterID[index]))
      index.id <- which(Cluster1$ClusterID %in% cluster.id)
      Cluster1$seqNum[index.id] <- Cluster1$seqNum[index.id] + re.num - overlap.num
      Cluster1 <- Cluster1[-index, ]
      new.cluster <- data.frame(ClusterID = as.numeric(rep(cluster.id, re.num)),
                                seqID = as.character(re.name),
                                seqLen = as.numeric(rep(0,re.num)),
                                seqNum = as.integer(rep(Cluster1$seqNum[index.id[1]], re.num)))
      Cluster1 <- rbind(Cluster1, new.cluster)
    } else {
      max.clusterid <- max.clusterid + 1
      new.cluster <- data.frame(ClusterID = as.numeric(rep(max.clusterid, re.num)),
                                seqID = as.character(re.name),
                                seqLen = as.numeric(rep(0,re.num)),
                                seqNum = as.integer(rep(re.num, re.num)))
      Cluster1 <- rbind(Cluster1, new.cluster)

    }
  }
  Cluster1 <- Cluster1[order(Cluster1$ClusterID), ]
  Cluster1$seqLen <- cdhit.result$seqLen[match(Cluster1$seqID, cdhit.result$seqID)]
  index.q <- match(Align$Qid, Cluster1$seqID)
  Align$QClusterID <- Cluster1$ClusterID[index.q]
  Align$QseqLen <- Cluster1$seqLen[index.q]
  Align$QseqNum <- Cluster1$seqNum[index.q]
  index.s <- match(Align$Sid, Cluster1$seqID)
  Align$SClusterID <- Cluster1$ClusterID[index.s]
  Align$SseqLen <- Cluster1$seqLen[index.s]
  Align$SseqNum <- Cluster1$seqNum[index.s]
  result <- list(alignment = Align, cluster = Cluster1)
  return(result)
}


#' @title Network graph showing the distribution of isoforms in a cluster.
#'
#' @description This function plots a clustering network form a data frames.
#'  It shows the distribution of isoforms in a cluster.
#'
#' @param cluster A data frame contation the information of isoform clusters.
#' @param cluster.id An interger vector indicating the ID of cluster.
#' @return This function returns a network graph object invisibly.
#' @examples
#' ##Loading example data
#' load(system.file("data","example_cluster.Rdata",package = "AStrap"))
#' g<- plotCluster(cluster.align$cluster,cluster.id=c("7","5"))
#' plot(g)
plotCluster <- function(cluster, cluster.id = c("2")) {
  Cluster <- as.data.frame(cluster)
  index <- which(Cluster$ClusterID %in% cluster.id)
  if (length(index) == 0) {
    stop(c("no exist this Cluster.id: ", cluster.id))
  }
  data <- data.frame(from = Cluster$ClusterID[index], to = Cluster$seqID[index])
  index.no <- setdiff(cluster.id, unique(data$from))
  if (length(index.no) > 0) {
    warning(c("no exist the following Cluster.id: ", index.no))
  }
  g <- graph_from_data_frame(data, directed = FALSE)
  V(g)$label <- V(g)$name
  V(g)$color <- rainbow(6, alpha = 0.6)[4]
  index.no <- which(as.character(names(V(g))) %in% as.character(unique(data$from)))
  V(g)[index.no]$label <- paste0("Cluster", V(g)[index.no]$label)
  V(g)[index.no]$color <- rainbow(1, alpha = 0.6)[1]
  V(g)[index.no]$size <- 13
  V(g)[index.no]$label.font <- 2
  V(g)[-index.no]$size <- 10
  V(g)$label.dist <- 0.4
  V(g)$label.color <- "#333300"
  V(g)$label.cex <- 0.6
  V(g)[index.no]$label.cex <- 0.8
  V(g)[index.no]$label.color <- "#333333"
  E(g)$width <- 2
  E(g)$color <- "grey"
  g$layout <- layout.fruchterman.reingold
  V(g)$frame.color <- NA
  return(g)
}


#' @title Network graph showing identities of pairwise alignments in a cluster.
#'
#' @description This function plot pairwise alignments of isoforms of the
#'   same cluster network graph. It shows the relationship between isoforms
#'   of the same cluster, where the identity of a pairwise alignment is
#'   denoted as a line with different thickness and color. The higher
#'   the identity, the thicker the line. Red line denotes the identity of 1.
#'
#' @param alignment A data frame contation the information of sequence alignments.
#' @param cluster.id An interger vector indicating the ID of cluster.
#' @return This function returns a network graph object invisibly.
#' @examples
#' ##Loading example data
#' load(system.file("data","example_cluster.Rdata",package = "AStrap"))
#' g <- plotAlign(cluster.align$alignment,cluster.id=c("7","5"))
#' plot(g)
plotAlign <- function(alignment, cluster.id = c("2")) {
  alignment <- as.data.frame(alignment)
  index <- which(alignment$SClusterID %in% cluster.id)
  if (length(index) == 0) {
    stop(c("no exist this Cluster.id: ", cluster.id))
  }
  data <- data.frame(from = alignment$Qid[index], to = alignment$Sid[index])
  index.no <- setdiff(cluster.id, unique(alignment$SClusterID[index]))
  if (length(index.no) > 0) {
    warning(c("no exist the following Cluster.id: ", index.no))
  }
  g <- graph_from_data_frame(data, directed = FALSE)
  V(g)$label <- V(g)$name
  V(g)$frame.color <- NA
  V(g)$size <- 10
  V(g)$label.dist <- 0.4
  V(g)$label.cex <- 0.6
  V(g)$label.color <- "#333333"
  E(g)$color <- "grey"
  g$layout <- layout.fruchterman.reingold
  mel.col <- rainbow(length(unique(V(g)$name)), alpha = 0.6)
  V(g)$color <- mel.col
  E(g)$width <- as.numeric(alignment$identity[index]/100)
  E(g)$color[which(E(g)$width == 1)] <- "lightcoral"
  return(g)
}



#' @title Read the pairwise sequence alignments from a file.
#'
#' @description This functions is a container for storing pairwise
#'   sequence alignments in a cluster.
#'
#' @param gmap.path A string for the name of a GFF3 file of pairwise
#'   sequence alignments output by GMAP.
#' @param cdhit.result A data frame holds isoform clusters.
#' @param recluster To obtain all potential clusters with more
#'   than one isoform for further analysis, a second-round mapping
#'   can be performed on clusters with single isoform. If TRUE (default),
#'   isoforms passing the second-round mapping were reassigned
#'   to the corresponding new clusters.
#' @param recluster.identity An numeric value indicating the
#'   threshold of identity of alignment for the second-round
#'   mapping (default: 0.7). Used only if recluster = TRUE.
#' @param recluster.coverage An numeric value indicating the
#'   threshold of coverage of alignment for the second-round
#'   mapping (default: 0.7). Used only if recluster = TRUE.
#' @return A list contation the pairwise sequence alignments
#'  in a cluster, including two elements: alignment and cluster.
#'  alignment is a data frame storage the pairwise alignment
#'  of isoforms in a cluster; cluster is a data frame storage
#'  the new clustering result after second-round mapping.
#' @examples
#' ##Load the example data
#' cdhit.path <- system.file("extdata","example_cdhitest.clstr",package = "AStrap")
#' gmap.path <- system.file("extdata","example_gmap.gff3",package = "AStrap")
#'
#' ##Read raw isoform clusters
#' raw.cluster <- readCDHIT(cdhit.path)
#'
#' ##Pairwise sequence alignments of isoforms of the same cluster.
#' cluster.align <- readGMAP(gmap.path,raw.cluster, recluster = TRUE,
#'                           recluster.identity = 0.7,
#'                           recluster.coverage = 0.7)
#' ##Pairwise sequence alignments
#' alignment <- cluster.align$alignment
#' head(alignment[,1:4])
#' ##Adjusted cluster
#' rew.cluster <- cluster.align$cluster
readGMAP <- function(gmap.path, cdhit.result, recluster = TRUE, recluster.identity = 0.7, recluster.coverage = 0.7) {
  gmap.align <- trimmGMAP(gmap.path)
  S.infor <- cdhit.result
  colnames(S.infor) <- paste("S", colnames(S.infor), sep = "")
  gmap.align <- merge(gmap.align, S.infor, by.x = "Sid", by.y = "SseqID")
  Q.infor <- cdhit.result
  colnames(Q.infor) <- paste("Q", colnames(Q.infor), sep = "")
  gmap.align <- merge(gmap.align, Q.infor, by.x = "Qid", by.y = "QseqID")
  if (recluster) {
    result <- reCluster(gmap.align, cdhit.result, re.id = recluster.identity, re.co = recluster.coverage)
    gmap.align <- result$alignment
    cdhit.result <- result$cluster
  }
  gmap.align$sameCluster <- 0
  gmap.align$sameCluster[gmap.align$SClusterID == gmap.align$QClusterID] <- 1
  gmap.align <- gmap.align[which(gmap.align$same == 1), ]
  result <- list(alignment = gmap.align, cluster = cdhit.result)
  return(result)
}



#' @title Identification AS events and classification AS types.
#'
#' @description This function can identify AS events and classify AS events from
#' transcripts without a reference.Two classification models trained on collected
#' AS data from rice (MSU7) and human (GRCH38) were integrated in AStrap, which
#' could be directly applied for distinguishing AS types for other species.
#' Certainly, users can also train a specific model on their own data sets
#' using several functions provided in AStrap, see the function of buildTrainModel.
#'
#' @param cluster.alignment A data frame holds pairwise alignment of isoforms in a cluster.
#' @param transcriptSeq An XStringSet object holds the transcript sequence.
#' @param trainModel A model object for which prediction is desired.
#' @param identity AS detection is performed if identity above a given threshold (default: 0.7).
#' @param coverage AS detection is performed if coverage above a given threshold (default: 0.7).
#' @param bias Maximum number of mismatches in pairwise sequence alignment (default: 0).
#' @param ASlength AS detection is performed if AS length above a given threshold (default: 0).
#' @return A list of AS events and AS types. This object contains three elements:
#'  ASevent, feature, predict. ASevent is a data frame stroage the identified AS events;
#'  feature is a feature data frame for AS events based on transcript, including 511 features;
#'  predict is the calssification of AS events. and the form of the value returned by predict
#'  depends on the class of its argument.
#' @examples
#' ##Loading transcript sequence
#' trSequence.path <- system.file("extdata","example_TRsequence.fasta",
#'                               package = "AStrap")
#' trSequence <-  readDNAStringSet(trSequence.path,format = "fasta")
#'
#' ##Loadinding example data output by CD-HIT-EST
#' cdhit.path <- system.file("extdata","example_cdhitest.clstr",
#'                          package = "AStrap")
#' raw.cluster <- readCDHIT(cdhit.path)
#'
#' ##Loading example data output by GMAP
#' gmap.path <- system.file("extdata","example_gmap.gff3",package = "AStrap")
#'
#' ##getting pairwise alignment of isoforms of the same cluster
#' cluster.align <- readGMAP(gmap.path,raw.cluster, recluster = TRUE,
#'                           recluster.identity = 0.7,
#'                           recluster.coverage = 0.7)
#' alignment <- cluster.align$alignment
#' rew.cluster <- cluster.align$cluster
#'
#' ##Loading rice model
#' load(system.file("data","rice_model.Rdata",package = "AStrap"))
#'
#' ##Identification and prediction based on rice RF-based model
#' result <- AStrap(alignment,trSequence,rice_RFmodel)
AStrap <- function(cluster.alignment, transcriptSeq, trainModel, identity = 0.7, coverage = 0.7, bias = 0, ASlength = 0) {
  bias <- bias + 1
  ASevent <- data.frame(Name = "", Target = 1, Qname = "", Sname = "", Qstart = 1, Qend = 1, Sstart = 1, Send = 1)
  ASevent$Name <- as.character(ASevent$Name)
  ASevent$Target <- as.numeric(ASevent$Target)
  ASevent$Qname <- as.character(ASevent$Qname)
  ASevent$Sname <- as.character(ASevent$Sname)
  ASevent$Qstart <- as.numeric(ASevent$Qstart)
  ASevent$Qend <- as.numeric(ASevent$Qend)
  ASevent$Sstart <- as.numeric(ASevent$Sstart)
  ASevent$Send <- as.numeric(ASevent$Send)
  as.align <- cluster.alignment[which(cluster.alignment$Alinum > 1), ]
  index.need <- intersect(which(as.align$identity >= identity * 100), which(as.align$Coverage >= coverage * 100))
  as.align <- as.align[index.need, ]
  as.align$Sali <- str_extract_all(as.align$Salign, "\\d+")
  as.align$Qali <- str_extract_all(as.align$Qalign, "\\d+")
  for (i in 1:nrow(as.align)) {
    S.align <- as.numeric(as.align$Sali[[i]])
    Q.align <- as.numeric(as.align$Qali[[i]])
    align.num <- as.numeric(as.align$Alinum[i])
    for (j in 1:(align.num - 1)) {
      S.dis <- S.align[2 * j + 1] - S.align[2 * j]
      Q.dis <- Q.align[2 * j + 1] - Q.align[2 * j]
      if ((Q.dis <= bias) && (S.dis > 1)) {
        TEMP <- data.frame(Name = as.align$NAME[i], Target = 1,
                           Qname = as.align$Qid[i], Sname = as.align$Sid[i],
                           Qstart = Q.align[2 * j], Qend = Q.align[2 * j + 1],
                           Sstart = S.align[2 * j], Send = S.align[2 * j + 1])
        ASevent <- rbind(ASevent, TEMP)

      } else if ((S.dis <= bias) && (Q.dis > 1)) {
        TEMP <- data.frame(Name = as.align$NAME[i], Target = 0,
                           Qname = as.align$Sid[i], Sname = as.align$Qid[i],
                           Qstart = S.align[2 * j],Qend = S.align[2 * j + 1],
                           Sstart = Q.align[2 * j], Send = Q.align[2 * j + 1])
        ASevent <- rbind(ASevent, TEMP)

      }
    }
  }
  ASevent <- ASevent[-1, ]
  rew.ASevent <- cbind(ASevent, as.align[match(as.character(ASevent$Name), as.character(as.align$NAME)), ])
  rew.ASevent <- extract_IsoSeq_tr(rew.ASevent, transcriptSeq)
  row.names(rew.ASevent) <- c(1:nrow(rew.ASevent))
  testfeature <- getFeature(rew.ASevent)
  pred <- predict(trainModel, testfeature)
  rew.ASevent$spliceSeq <- paste0(subseq(rew.ASevent$seq, start = 1, width = 2), "-", subseq(rew.ASevent$seq, end = rew.ASevent$length, width = 2))
  rew.ASevent$Predict <- as.character(pred)
  ASevent <- data.frame(Qid = rew.ASevent$Qname, Qlength = rew.ASevent$QseqLen,
                        Sid = rew.ASevent$Sname, Slength = rew.ASevent$SseqLen,
                        Clusterid = rew.ASevent$QClusterID,CluterSeqNum = rew.ASevent$QseqNum,
                        Qstart = rew.ASevent$Qstart, Qend = rew.ASevent$Qend,
                        Sstart = rew.ASevent$Sstart, Send = rew.ASevent$Send,
                        identity = rew.ASevent$identity, coverage = rew.ASevent$Coverage,
                        ASlength = rew.ASevent$length,Qalign = rew.ASevent$Qalign,
                        Salign =rew.ASevent$Salign,
                        prediction = rew.ASevent$Predict, spliceSeq = rew.ASevent$spliceSeq)
  index <- which(rew.ASevent$Target == 0)
  ASevent$Qlength[index] <- rew.ASevent$SseqLen[index]
  ASevent$Slength[index] <- rew.ASevent$QseqLen[index]
  ASevent$Qalign[index] = rew.ASevent$Qalign[index]
  ASevent$Salign[index] = rew.ASevent$Salign[index]
  ASevent <- ASevent[ASevent$ASlength > ASlength, ]
  return(list(ASevent = ASevent, feature = testfeature, predict = pred))
}



#' @title Extract the sequence around splice sites based on transcript.
#'
#' @description A function for extracting a set of subsequences surrounding splice sites
#' from a sequence container like XStraingSet.
#'
#' @param ASdata A data frame holds the coordinate of splice sites with as least
#'   three columns: Sname, Sstart, Send. Sname is the name of subject in the pariwise
#'   sequence alignment; Sstart is the start of subject in the pairwise sequence
#'   alignment; Send is the end of subject in the pairwise sequence;
#'
#' @param transcriptSeq An XStringSet object holds the transcript sequence.
#' @return A alternative splicing data frame with sequence information.
#' @examples
#' ##Loading pairwise alignment data
#' load(system.file("data","sample_Aligndata.Rdata",package = "AStrap"))
#' head(Aligndata)
#'
#' ##Loading transcript sequence
#' trSequence.path <- system.file("extdata","example_TRsequence.fasta",
#'                               package = "AStrap")
#' trSequence <-  readDNAStringSet(trSequence.path,format = "fasta")
#'
#' ##Extract sequence around splice site based on the transcript sequence
#' Aligndata <- extract_IsoSeq_tr(Aligndata,trSequence )
#' colnames(Aligndata)
#' head(Aligndata$Ddown10)
#'
#' head(Aligndata$Aup10)
extract_IsoSeq_tr <- function(ASdata, transcriptSeq) {
  ASdata <- as.data.frame(ASdata)
  ASdata$num <- as.numeric(c(1:nrow(ASdata)))
  ASdata$length <- ASdata$Send - ASdata$Sstart - 1
  start <- as.integer(ASdata$Sstart + 1)
  end <- as.integer(ASdata$Send - 1)
  ASdata$seq <- as.character(DNAStringSet(subseq(getSeq(transcriptSeq, ASdata$Sname), start = start, end = end)))
  TEMP <- ASdata$SseqLen - start + 1
  Ddown10 <- rep(10, nrow(ASdata))
  Ddown10[which(TEMP < 10)] <- TEMP[which(TEMP < 10)]
  Ddown20 <- rep(20, nrow(ASdata))
  Ddown20[which(TEMP < 20)] <- TEMP[which(TEMP < 20)]
  Aup10 <- rep(10, nrow(ASdata))
  Aup10[which(end < 10)] <- end[which(end < 10)]
  Aup20 <- rep(20, nrow(ASdata))
  Aup20[which(end < 20)] <- end[which(end < 20)]
  Aup30 <- rep(30, nrow(ASdata))
  Aup30[which(end < 30)] <- end[which(end < 30)]
  ASdata$Ddown10 <- as.character(DNAStringSet(subseq(getSeq(transcriptSeq, ASdata$Sname), start = start, width = Ddown10)))
  ASdata$Ddown20 <- as.character(DNAStringSet(subseq(getSeq(transcriptSeq, ASdata$Sname), start = start, width = Ddown20)))
  ASdata$Aup10 <- as.character(DNAStringSet(subseq(getSeq(transcriptSeq, ASdata$Sname), end = end, width = Aup10)))
  ASdata$Aup20 <- as.character(DNAStringSet(subseq(getSeq(transcriptSeq, ASdata$Sname), end = end, width = Aup20)))
  ASdata$Aup30 <- as.character(DNAStringSet(subseq(getSeq(transcriptSeq, ASdata$Sname), end = end, width = Aup30)))
  ASdata$donorSeq <- as.character(DNAStringSet(subseq(getSeq(transcriptSeq, ASdata$Sname), start = as.integer(ASdata$Sstart), width = 6)))
  ASdata$acceptorSeq <- as.character(DNAStringSet(subseq(getSeq(transcriptSeq, ASdata$Sname), end = as.integer(ASdata$Send), width = 6)))
  return(ASdata)
}


featureGeLength <- function(ASdata) {
  ASdata <- as.data.frame(ASdata)
  ASdata$num <- as.numeric(c(1:nrow(ASdata)))
  if (is.null(ASdata$length)) {
    AltD <- ASdata[grep("AltD", ASdata$type, ignore.case = TRUE), ]
    AltD$length <- abs(AltD$coord3 - AltD$coord1)
    AltA <- ASdata[grep("AltA", ASdata$type, ignore.case = TRUE), ]
    AltA$length <- abs(AltA$coord2 - AltA$coord4)
    ES <- ASdata[grep("ES", ASdata$type, ignore.case = TRUE), ]
    ES$length <- abs(ES$coord3 - ES$coord2) + 1
    IR <- ASdata[grep("IR", ASdata$type, ignore.case = TRUE), ]
    IR$length <- abs(IR$coord1 - IR$coord2) - 1
    AS <- rbind(AltA, AltD, IR, ES)
    AS <- AS[order(ASdata$num), ]
    AS$num <- NULL
    return(AS)
  } else {
    warning("there is already a column called length,
            Please make sure that the column value represents the length of the AS ")
  }
  }


featureTHREE <- function(ASlength) {
  ASlength <- as.numeric(ASlength)
  ASthree <- ASlength%%3
  ASthree[which(ASthree != 0)] <- 1
  return(ASthree)
}


#' @title Extract the sequence around splice sites based on genome.
#'
#' @description A function for extracting a set of subsequences surrounding splice sites
#' from a sequence container like a BSgenome or XStraingSet.
#'
#' @param ASdata A data frame holds an alternative splicing database,
#'   including intron retention (IR), exon skipping (ES), alternative
#'   donor sites (AltD), and alternative acceptor sites (AltA).
#' @param genome A BSgenome or XStringSet object.
#' @return A alternative splicing data frame with sequence information.
#' @examples
#' ##Load the rice genome (MSU)
#' library("BSgenome.Osativa.MSU.MSU7")
#'
#' ##Load the alternative splicing database
#' path <- system.file("extdata","sample_riceAS.txt",
#'                      package = "AStrap")
#' rice_ASdata <-read.table(path,sep="\t",head = TRUE,
#'                          stringsAsFactors = FALSE)
#'
#' ##Extract the sequence around splice sites based on genome
#' AS.data <- extract_IsoSeq_ge(rice_ASdata,Osativa)
extract_IsoSeq_ge <- function(ASdata, genome) {
  ASdata <- as.data.frame(ASdata)
  if (is.null(ASdata$length)) {
    ASdata <- featureGeLength(ASdata)
  } else {
    warning("there is already a column called length,
            Please make sure that the column value represents the length of the AS ")
  }
  ASdata$num <- as.numeric(c(1:nrow(ASdata)))
  ASdata$seq <- "N"
  ASdata$up <- "N"
  ASdata$down <- "N"
  AltD <- ASdata[grep("AltD", ASdata$type, ignore.case = TRUE), ]
  AltA <- ASdata[grep("AltA", ASdata$type, ignore.case = TRUE), ]
  ES <- ASdata[grep("ES", ASdata$type, ignore.case = TRUE), ]
  IR <- ASdata[grep("IR", ASdata$type, ignore.case = TRUE), ]
  ASdata = ASdata[1, ]
  if (nrow(AltD) > 0) {
    altd <- getGeAltDsequence(AltD, genome)
    ASdata <- rbind(ASdata, altd)
  }
  if (nrow(AltA) > 0) {
    alta <- getGeAltAsequence(AltA, genome)
    ASdata <- rbind(ASdata, alta)
  }
  if (nrow(ES) > 0) {
    ex <- getGeESsequence(ES, genome)
    ASdata <- rbind(ASdata, ex)
  }
  if (nrow(IR) > 0) {
    ir <- getGeIRsequence(IR, genome)
    ASdata <- rbind(ASdata, ir)
    PWMir <- ir[!duplicated(paste0(ir$up, ir$seq, ir$down)), ]
    PWMir$gt <- subseq(PWMir$up, start = 50, width = 6)  #XGTXXX
    PWMir$ag <- subseq(PWMir$down, start = 47, width = 6)  #XXXAGX
    pfm_ir_gt <- consensusMatrix(PWMir$gt)
    PWM_Donor <<- PWM(pfm_ir_gt)
    pfm_ir_ag <- consensusMatrix(PWMir$ag)
    PWM_Acceptor <<- PWM(pfm_ir_ag)
  }
  ASdata <- ASdata[-1, ]
  ASdata$Ddown10 <- as.character(subseq(ASdata$up, start = 51, width = 10))
  ASdata$Ddown20 <- as.character(subseq(ASdata$up, start = 51, width = 20))
  ASdata$Aup10 <- as.character(subseq(ASdata$down, start = 42, width = 10))
  ASdata$Aup20 <- as.character(subseq(ASdata$down, start = 32, width = 20))
  ASdata$Aup30 <- as.character(subseq(ASdata$down, start = 22, width = 30))
  ASdata$donorSeq <- as.character(subseq(ASdata$up, start = 50, width = 6))
  ASdata$acceptorSeq <- as.character(subseq(ASdata$down, start = 47, width = 6))
  ASdata <- ASdata[order(ASdata$num), ]
  ASdata$num <- NULL
  return(ASdata)
  }


getGeIRsequence <- function(IR, genome) {
  ir <- IR[1, ]
  pos <- which(IR$strand == "+")
  if (length(pos) > 0) {
    irp <- IR[pos, ]
    start <- as.integer(irp$coord1 + 1)
    end <- as.integer(irp$coord2 - 1)
    index <- which(irp$coord1 > irp$coord2)
    if (length(index) > 0) {
      c1 <- irp$coord1[index]
      c2 <- irp$coord2[index]
      irp$coord1[index] <- c2
      irp$coord2[index] <- c1
    }
    irp$seq <- as.character(getSeq(genome, irp$chr, start, end))
    irp$up <- as.character(getSeq(genome, irp$chr, start - 50, start + 50))
    irp$down <- as.character(getSeq(genome, irp$chr, end - 50, end + 50))
    ir <- rbind(ir, irp)
  }
  neg <- which(IR$strand == "-")
  if (length(neg) > 0) {
    irn <- IR[neg, ]
    start <- as.integer(irn$coord2 + 1)
    end <- as.integer(irn$coord1 - 1)
    index <- which(irn$coord1 < irn$coord2)
    if (length(index) > 0) {
      c1 <- irn$coord1[index]
      c2 <- irn$coord2[index]
      irn$coord1[index] <- c2
      irn$coord2[index]<- c1
    }
    irn$seq <- as.character(reverseComplement(getSeq(genome, irn$chr, start, end)))
    irn$up <- as.character(reverseComplement(getSeq(genome, irn$chr, end - 50, end + 50)))
    irn$down <- as.character(reverseComplement(getSeq(genome, irn$chr, start - 50, start + 50)))
    ir <- rbind(ir, irn)
  }
  ir <- ir[-1, ]
  return(ir)
}


getGeAltAsequence <- function(AltA, genome) {
  alta <- AltA[1, ]
  if (length(which(AltA$coord1 != AltA$coord3)) > 0) {
    stop("error of input data:
         for AltA, the Coordinate1 and Coordinate3 are  different")
  }
  pos <- which(AltA$strand == "+")
  altap <- AltA[pos, ]
  if (length(pos) > 0) {
    index <- which(altap$coord2 > altap$coord4)  #just no
    if (length(index) > 0) {
      c2 <- altap$coord2[index]
      c4 <- altap$coord4[index]
      altap$coord2[index] <- c4
      altap$coord4[index] <- c2
    }

    if (length(which(altap$coord1 > altap$coord2)) > 0) {
      stop("error of input data:
           for AltA, the Coordinate1 is behind Coordinate2")
    }
    start <- as.integer(altap$coord2)
    end <- as.integer(altap$coord4 - 1)
    start1 <- as.integer(altap$coord1)
    altap$seq <- as.character(getSeq(genome, altap$chr, start, end))
    altap$up <- as.character(paste0(as.character(getSeq(genome, altap$chr, start1 - 49, start1)),
                                    as.character(getSeq(genome, altap$chr, start,start + 50))))
    altap$down <- as.character(getSeq(genome, altap$chr, end - 50, end + 50))
    alta <- rbind(alta, altap)
    }
  neg <- which(AltA$strand == "-")
  if (length(neg) > 0) {
    altan <- AltA[neg, ]
    index <- which(altan$coord2 < altan$coord4)  #just no
    if (length(index) > 0) {
      c2 <- altan$coord2[index]
      c4 <- altan$coord4[index]
      altan$coord2[index] <- c4
      altan$coord4[index] <- c2
    }

    if (length(which(altan$coord1 < altan$coord2)) > 0) {
      stop("error of input data:
           for AltA, the Coordinate1 is behind Coordinate2")
    }
    start <- as.integer(altan$coord4 + 1)
    end <- as.integer(altan$coord2)
    start1 <- as.integer(altan$coord1)
    altan$seq <- as.character(reverseComplement(getSeq(genome, altan$chr, start, end)))
    altan$up <- as.character(reverseComplement(DNAStringSet(paste0(as.character(getSeq(genome, altan$chr, end - 50, end)),
                                                                   as.character(getSeq(genome, altan$chr, start1, start1 + 49))))))
    altan$down <- as.character(reverseComplement(DNAStringSet(getSeq(genome, altan$chr, start - 50, start + 50))))
    alta <- rbind(alta, altan)
    }
  alta <- alta[-1, ]
  less <- which(alta$length <= 50)
  if (length(less) > 0) {
    alta$down[less] <- as.character(paste0(subseq(alta$up[less], alta$length[less], 50),
                                           alta$seq[less], subseq(alta$down[less], 52, 101)))
  }
  return(alta)
  }


getGeAltDsequence <- function(AltD, genome) {
  altd <- AltD[1, ]
  if (length(which(AltD$coord2 != AltD$coord4))) {
    stop("error of input data:
         for AltD, the Coordinate2 and Coordinate4 are  different")
  }
  pos <- which(AltD$strand == "+")
  if (length(pos) > 0) {
    altdp <- AltD[pos, ]
    index <- which(altdp$coord1 > altdp$coord3)
    if (length(index) > 0) {
      c1 <- altdp$coord1[index]
      c3 <- altdp$coord3[index]
      altdp$coord1[index] <- c3
      altdp$coord3[index] <- c1
    }
    if (length(which(altdp$coord3 > altdp$coord2))) {
      stop("error of input data:
           for AltD, the Coordinate3 is behind Coordinate2")
    }
    start <- as.integer(altdp$coord1 + 1)
    end <- as.integer(altdp$coord3)
    start1 <- as.integer(altdp$coord2)
    altdp$seq <- as.character(getSeq(genome, altdp$chr, start, end))
    altdp$up <- as.character(getSeq(genome, altdp$chr, start - 50, start + 50))
    altdp$down <- as.character(DNAStringSet(paste0(as.character(getSeq(genome, altdp$chr, end - 50, end)),
                                                   as.character(getSeq(genome, altdp$chr, start1, start1 + 49)))))
    altd <- rbind(altd, altdp)
    }

  neg <- which(AltD$strand == "-")
  if (length(neg) > 0) {
    altdn <- AltD[neg, ]
    index <- which(altdn$coord3 > altdn$coord1)
    if (length(index) > 0) {
      c3 <- altdn$coord3[index]
      c1 <- altdn$coord1[index]
      altdn$coord3[index] <- c1
      altdn$coord1[index] <- c3
    }
    if (length(which(altdn$coord3 < altdn$coord2))) {
      stop("error of input data:
           for AltD, the Coordinate3 is behind Coordinate2")
    }
    start <- as.integer(altdn$coord3)
    end <- as.integer(altdn$coord1 - 1)
    start1 <- as.integer(altdn$coord2)
    altdn$seq <-as.character(reverseComplement(DNAStringSet(getSeq(genome, altdn$chr, start, end))))
    altdn$up <- as.character(reverseComplement(DNAStringSet(getSeq(genome, altdn$chr, end - 50, end + 50))))
    altdn$down <- as.character(reverseComplement(DNAStringSet(paste0(as.character(getSeq(genome, altdn$chr, start1 - 49, start1)),
                                                                     as.character(getSeq(genome, altdn$chr, start, start + 50))))))
    altd <- rbind(altd, altdn)
    }
  altd <- altd[-1, ]
  less <- which(altd$length <= 50)
  if (length(less) > 0) {
    altd$up[less] <- as.character(paste0(subseq(altd$up[less], 1, 50), altd$seq[less],
                                         subseq(altd$down[less], 52, 102 - altd$length[less])))
  }
  return(altd)
  }


getGeESsequence <- function(ES, genome) {
  ex <- ES[1, ]
  pos <- which(ES$strand == "+")
  if (length(pos) > 0) {
    exp <- ES[pos, ]
    if (length(which(exp$coord2 < exp$coord1))) {
      stop("error of input data:
           for ES, the Coordinate1 is behind Coordinate2")
    }
    if (length(which(exp$coord3 < exp$coord2))) {
      stop("error of input data:
           for ES, the Coordinate2 is behind Coordinate3")
    }
    if (length(which(exp$coord4 < exp$coord3))) {
      stop("error of input data:
           for ES, the Coordinate3 is behind Coordinate4")
    }
    start <- as.integer(exp$coord2)
    end <- as.integer(exp$coord3)
    exp$seq <- as.character(getSeq(genome, exp$chr, start, end))
    start1 <- as.integer(exp$coord1)
    end4 <- as.integer(exp$coord4)
    exp$up <- as.character(paste0(as.character(getSeq(genome, exp$chr, start1 - 49, start1)),
                                  as.character(getSeq(genome, exp$chr, start, start + 50))))
    exp$down <- as.character(paste0(as.character(getSeq(genome, exp$chr, end - 50, end)),
                                    as.character(getSeq(genome, exp$chr, end4, end4 + 49))))

    ex <- rbind(ex, exp)
    }

  neg <- which(ES$strand == "-")
  if (length(neg > 0)) {
    exn <- ES[neg, ]
    if (length(which(exn$coord2 > exn$coord1))) {
      stop("error of input data:
           for AS type ES, the Coordinate1 is behind Coordinate2")
    }
    if (length(which(exn$coord3 > exn$coord2))) {
      stop("error of input data:
           for AS type ES, the Coordinate2 is behind Coordinate3")
    }
    if (length(which(exn$coord4 > exn$coord3))) {
      stop("error of input data:
           for AS type ES, the Coordinate3 is behind Coordinate4")
    }
    start <- as.integer(exn$coord3)
    end <- as.integer(exn$coord2)
    start1 <- as.integer(exn$coord4)
    end4 <- as.integer(exn$coord1)
    exn$seq <- as.character(reverseComplement(DNAStringSet(getSeq(genome, exn$chr, start, end))))
    exn$up <- as.character(reverseComplement(DNAStringSet(paste0(as.character(getSeq(genome, exn$chr, end - 50, end)),
                                                                 as.character(getSeq(genome, exn$chr, end4, end4 + 49))))))
    exn$down <- as.character(reverseComplement(DNAStringSet(paste0(as.character(getSeq(genome, exn$chr, start1 - 49, start1)),
                                                                   as.character(getSeq(genome, exn$chr, start, start + 50))))))

    ex <- rbind(ex, exn)
    }
  ex <- ex[-1, ]
  less <- which(ex$length <= 50)
  if (length(less) > 0) {
    ex$up[less] <- as.character(paste0(subseq(ex$up[less], 1, 50), ex$seq[less],
                                       subseq(ex$down[less], 52, 102 - ex$length[less])))
    ex$down[less] <- as.character(paste0(subseq(ex$up[less], ex$length[less], 50),
                                         ex$seq[less], subseq(ex$down[less], 52, 101)))
  }
  return(ex)
    }


#' @title Classification models.
#'
#' @description This function builds various classification models,
#'   including support vector machine (SVM), random forests (RF),
#'   and adaptive boosting (AdaBoost).
#'
#' @param ASdata A data frame including the coordinates of splice sites,
#'   class label and the sequence around splice sites. The "type" column is
#'   a vector of class label comprising of "AltA","AltD","ES" and "IR".
#' @param chooseNum A interger for the number of AS events from each AS type
#'   for building classification model.
#' @param proTrain The proportion of training dataset using random sampling.
#' @param proTest The proportion of testing dataset using random sampling.
#' @param ASlength AS data is trimmed if AS length below a given threshold.
#' @param classifier A string for the classification method. This must be one of
#'   the string "svm", "rf", "adaboost", not case sensitive.
#' @param use.all Whether to use all alternative splicing dataset for building
#'   classificaiton model (default: FALSE).
#' @return This function returns a fitted model with eight elements, including
#'   trainset, testset, model, predict, accuracy, confusion, evaluate, ROC.
#'   trainset is the training data set;
#'   testset is the testing data set;
#'   model is the fitted model;
#'   predict is the predicted classification results;
#'   accuracy is the  prediction accuracy;
#'   confusion is the confusion matrix of the prediction;
#'   evaluate is the evaluation matrix of the classification,
#'   including precition, sp, recall, f1;
#'   ROC: A ROC curve.
#' @examples
#' ##Loading example alternative splicing data
#' path <- system.file("extdata","sample_riceAS.txt",package = "AStrap")
#' rice_ASdata <-read.table(path,sep="\t",head = TRUE,stringsAsFactors = FALSE)
#' head(rice_ASdata)
#'
#' ##Loading geneome using the package of BSgenome
#' library("BSgenome.Osativa.MSU.MSU7")
#' rice_ASdata<- extract_IsoSeq_ge(rice_ASdata,Osativa)
#' names(rice_ASdata)
#'
#' ##Classification model building based on random forest method
#' library(randomForest)
#' library(ROCR)
#' library(ggplot2)
#' model <- buildTrainModel(rice_ASdata, chooseNum = 100,
#'                        proTrain = 2/3, proTest = 1/3,ASlength =0,
#'                        classifier = "rf", use.all = FALSE)
#' ##Performance evaluation
#' names(model)
#' model$evaluate
#' model$confusion
#' model$accuracy
#'
#' ##Or classification model building based on  SVM method
#' library(e1071)
#' library(ROCR)
#' library(ggplot2)
#' model <- buildTrainModel(rice_ASdata, chooseNum = 100,
#'                        proTrain = 2/3, proTest = 1/3,ASlength =0,
#'                        classifier = "svm", use.all = FALSE)
#'
#' ##Or classification model building based on  AdaBoost method
#' library(adabag)
#' library(ROCR)
#' library(ggplot2)
#' model <- buildTrainModel(rice_ASdata, chooseNum = 100,
#'                        proTrain = 2/3, proTest = 1/3,ASlength =0,
#'                        classifier = "adaboost", use.all = FALSE)
buildTrainModel <- function(ASdata, chooseNum = 1000, proTrain = 2/3, proTest = 1/3, ASlength =0, classifier = "rf", use.all = FALSE) {
  ASdata<-ASdata[which(ASdata$length>=ASlength),]
  classifier <- tolower(classifier)
  ASdata <- as.data.frame(ASdata)
  ASdata$type <- tolower(ASdata$type)
  ASdata$type[which(ASdata$type=="alta")] <- "AltA"
  ASdata$type[which(ASdata$type=="altd")] <- "AltD"
  ASdata$type[which(ASdata$type=="es")] <- "ES"
  ASdata$type[which(ASdata$type=="ir")] <- "IR"
  AltD <- ASdata[grep("AltD", ASdata$type, ignore.case = TRUE), ]
  AltA <- ASdata[grep("AltA", ASdata$type, ignore.case = TRUE), ]
  ES <- ASdata[grep("ES", ASdata$type, ignore.case = TRUE), ]
  IR <- ASdata[grep("IR", ASdata$type, ignore.case = TRUE), ]
  if (!use.all) {
    if (nrow(AltD) > chooseNum) {
      NAltD <- AltD[sample(1:nrow(AltD), chooseNum, replace = F), ]
    } else {
      NAltD <- AltD
    }
    if (nrow(AltA) > chooseNum) {
      NAltA <- AltA[sample(1:nrow(AltA), chooseNum, replace = F), ]
    } else {
      NAltA <- AltA
    }
    if (nrow(ES) > chooseNum) {
      NES <- ES[sample(1:nrow(ES), chooseNum, replace = F), ]
    } else {
      NES <- ES
    }
    if (nrow(IR) > chooseNum) {
      NIR <- IR[sample(1:nrow(IR), chooseNum, replace = F), ]
    } else {
      NIR <- IR
    }
  } else {
    NAltD <- AltD
    NAltA <- AltA
    NES <- ES
    NIR <- IR
  }

  ind <- sample(2, nrow(NAltD), replace = TRUE, prob = c(proTrain, proTest))  #2/3%为训练集 1/3为测试集
  trainAD <- NAltD[ind == 1, ]
  testAD <- NAltD[ind == 2, ]
  ind <- sample(2, nrow(NAltA), replace = TRUE, prob = c(proTrain, proTest))  #2/3%为训练集 1/3为测试集
  trainAA <- NAltA[ind == 1, ]
  testAA <- NAltA[ind == 2, ]
  ind <- sample(2, nrow(NES), replace = TRUE, prob = c(proTrain, proTest))  #2/3%为训练集 1/3为测试集
  trainES <- NES[ind == 1, ]
  testES <- NES[ind == 2, ]
  ind <- sample(2, nrow(NIR), replace = TRUE, prob = c(proTrain, proTest))  #2/3%为训练集 1/3为测试集
  trainIR <- NIR[ind == 1, ]
  testIR <- NIR[ind == 2, ]
  trainset <- rbind(trainAD, trainAA, trainES, trainIR)
  row.names(trainset) <- c(1:nrow(trainset))
  trainfeature <- getFeature(trainset)
  train <- data.frame(trainfeature, class = as.factor(trainset$type))
  testset <- rbind(testAD, testAA, testES, testIR)
  row.names(testset) <- c(1:nrow(testset))
  testfeature <- getFeature(testset)
  test <- data.frame(testfeature, class = as.factor(testset$type))
  return(chooseMethod(train = train, test = test, classifier = classifier))
}



chooseMethod <- function(train, test, classifier = "rf") {
  if (classifier == "svm") {
    model <- svm(class ~ ., data = train, kernel = "radial", cross = 10, probability = TRUE)
    if (nrow(test) > 0) {
      pred <- predict(model, test[, -ncol(test)], probability = TRUE, decision.values = TRUE)
      proba <- as.data.frame(attr(pred, "probabilities"))
      ROC <- plotROC(proba, test$class)
      result <- as.matrix(table(pred = pred, true = test$class))
      accuracy <- sum(diag(result))/sum(result)
      tp <- diag(result)
      fp <- rowSums(result) - diag(result)
      fn <- colSums(result) - diag(result)
      tn <- sum(result) - tp - fn - fp
      sp <- tn/(tn + fp)
      precision <- diag(result)/rowSums(result)
      recall <- (diag(result)/colSums(result))
      f1 <- 2 * precision * recall/(precision + recall)
      macroPrecision <- mean(precision)
      macroRecall <- mean(recall)
      macroF1 <- mean(f1)
      macroSP <- mean(sp)
      temp <- rbind(data.frame(precision, sp, recall, f1), mean = c(macroPrecision, macroSP, macroRecall, macroF1))
      return(list(trainSet = train, testSet = test, svmModel = model,
                  svmPre = pred, accuracy = accuracy,
                  confusion = result, evaluate = temp, ROC = ROC))
    } else {
      return(list(trainSet = train, testSet = test, svmModel = model))
    }

  } else if (classifier == "rf") {
    model <- randomForest(class ~ ., data = train, proximity = TRUE, type = "classification", importance = TRUE)
    if (nrow(test) > 0) {
      pred <- predict(model, test[, -ncol(test)])
      proba <- as.data.frame(predict(model, test[, -ncol(test)], type = "prob"))
      ROC <- plotROC(proba, test$class)
      result <- as.matrix(table(pred = pred, true = test$class))
      accuracy <- sum(diag(result))/sum(result)
      tp <- diag(result)
      fp <- rowSums(result) - diag(result)
      fn <- colSums(result) - diag(result)
      tn <- sum(result) - tp - fn - fp
      sp <- tn/(tn + fp)
      precision <- diag(result)/rowSums(result)
      recall <- (diag(result)/colSums(result))
      f1 <- 2 * precision * recall/(precision + recall)
      macroPrecision <- mean(precision)
      macroRecall <- mean(recall)
      macroF1 <- mean(f1)
      macroSP <- mean(sp)
      temp <- rbind(data.frame(precision, sp, recall, f1), mean = c(macroPrecision, macroSP, macroRecall, macroF1))
      return(list(trainSet = train, testSet = test, randomforestModel = model,
                  randomforestPre = pred, accuracy = accuracy, confusion = result,
                  evaluate = temp, ROC = ROC))
    } else {
      return(list(trainSet = train, testSet = test, randomforestModel = model))
    }

  } else if (classifier == "adaboost") {
    model <- boosting(class ~ ., data = train, pro = TRUE)
    #model <- boosting.cv(class ~ ., data = train,v=10)
    if (nrow(test) > 0) {
      pred <- predict(model, test[, -ncol(test)])
      proba <- as.data.frame(pred$prob)
      names(proba) <- c("AltA", "AltD", "ES", "IR")
      ROC <- plotROC(proba, test$class)
      result <- as.matrix(table(pred = pred$class, true = test$class))
      accuracy <- sum(diag(result))/sum(result)
      tp <- diag(result)
      fp <- rowSums(result) - diag(result)
      fn <- colSums(result) - diag(result)
      tn <- sum(result) - tp - fn - fp
      sp <- tn/(tn + fp)
      precision <- diag(result)/rowSums(result)
      recall <- (diag(result)/colSums(result))
      f1 <- 2 * precision * recall/(precision + recall)
      macroPrecision <- mean(precision)
      macroRecall <- mean(recall)
      macroF1 <- mean(f1)
      macroSP <- mean(sp)
      temp <- rbind(data.frame(precision, sp, recall, f1), mean = c(macroPrecision, macroSP, macroRecall, macroF1))
      return(list(trainSet = train, testSet = test, adaboosModel = model,
                  adaboostPre = pred, accuracy = accuracy, confusion = result,
                  evaluate = temp, ROC = ROC))
    } else {
      return(list(trainSet = train, testSet = test, adaboosModel = model))
    }

  } else {
    stop("the input method is incorrect,
         please input one of the following classification methods:svm,randomforest,adaboost")

  }
}



#' @title Feature construction.
#'
#' @description Extract features from the sequences around splice sites,
#'   A total of 511 features.
#'
#' @param ASdata A data frame with nine columns: donorSeq, acceptorSeq,
#'   seq, length, Ddown10, Ddown20, Aup10, Aup20, Aup30.
#'   donorSeq is the sequence of of the [-2, +3] region of donor sites;
#'   acceptorSeq is the sequence of the [-2, +3] region of acceptor site;
#'   seq is the sequence of the splice region;
#'   length is the length of the splice region;
#'   Ddown10 is the sequence of the downstream 10bp region of donor site;
#'   Ddown20 is the sequence of the downstream 20bp region of donor site;
#'   Aup10 is the sequence of the upstream 10bp region of acceptor site;
#'   Aup20 is the sequence of the upstream 20bp region of acceptor site;
#'   Aup30 is the sequence of the upstream 30bp region of acceptor site;
#' @return A list of features.
#' @examples
#' ##Example 1: based on pairwise sequence alignment and transcript
#' ##sequences.
#' ##Loading pairwise alignment data
#' load(system.file("data","sample_Aligndata.Rdata",package = "AStrap"))
#' head(Aligndata)
#'
#' ##Loading transcript sequence
#' trSequence.path <- system.file("extdata","example_TRsequence.fasta",
#'                               package = "AStrap")
#' trSequence <-  readDNAStringSet(trSequence.path,format = "fasta")
#'
#' ##Extract sequence around splice site based on the transcript sequence
#' Aligndata <- extract_IsoSeq_tr(Aligndata,trSequence )
#' colnames(Aligndata)
#' head(Aligndata$Ddown10)
#' head(Aligndata$Aup10)
#' ##Loading the consensus matrix of the sequences of the [-2,+3]
#' ##region of acceptor sites
#' load(system.file("data","example_PWM_acceptor.Rdata",
#'                  package = "AStrap"))
#' ##Loading the consensus matrix of the sequences of the [-2,+3]
#' ##region of donor sites
#' load(system.file("data","example_PWM_donor.Rdata",
#'                  package = "AStrap"))
#' ##Construction the feature space
#' feature <- getFeature(Aligndata)
#' ncol(feature)
#'
#' ##Example2: based on the alternative splicing database and genome.
#' path <- system.file("extdata","sample_riceAS.txt",package = "AStrap")
#' rice_ASdata <-read.table(path,sep="\t",head = TRUE,
#'                          stringsAsFactors = FALSE)
#' head(rice_ASdata)
#' ##Loading geneome using the package of BSgenome
#' library("BSgenome.Osativa.MSU.MSU7")
#' rice_ASdata<- extract_IsoSeq_ge(rice_ASdata,Osativa)
getFeature <- function(ASdata) {
  ASdata <- as.data.frame(ASdata)
  DonorMotif <- featureDonorMotif(ASdata$donorSeq, PWM_Donor)
  AcceptorMotif <- featureAcceptorMotif(ASdata$acceptorSeq, PWM_Acceptor)
  GT <- rep(0, nrow(ASdata))
  GT[grep("^GT", ASdata$seq)] <- 1
  GC <- rep(0, nrow(ASdata))
  GC[grep("^GC", ASdata$seq)]<- 1
  AG <- rep(0, nrow(ASdata))
  AG[grep("AG$", ASdata$seq)] <- 1
  AT <- rep(0, nrow(ASdata))
  AT[grep("^AT", ASdata$seq)] <- 1
  AC <- rep(0, nrow(ASdata))
  AC[grep("AC$", ASdata$seq)] <- 1
  GCpro <- letterFrequency(DNAStringSet(ASdata$seq), "GC", as.prob = TRUE)
  Stopco <- vcountPattern("TAG", DNAStringSet(ASdata$seq)) +
    vcountPattern("TAA", DNAStringSet(ASdata$seq)) +
    vcountPattern("TGA", DNAStringSet(ASdata$seq))
  three <- featureTHREE(ASdata$length)
  ASnuFre <- featureASFrequency(ASdata)
  ASCTD <- featureASCTD(ASdata)

  summary.feature <- data.frame(length = ASdata$length, Muthree = three,
                                donorGT = GT, donorGC = GC, donorAT = AT,
                                acceptorAG = AG, acceptorAC = AC,
                                donorMotif = DonorMotif, acceptorMotif = AcceptorMotif,
                                ASnuFre , ASCTD, GC = GCpro, Stopdocon = Stopco)
  return(summary.feature)

}


featureASFrequency <- function(ASdata) {
  up1 <- oligonucleotideFrequency(DNAStringSet(ASdata$Ddown10), width = 3, step = 1,
                                  as.prob = TRUE, as.array = FALSE, fast.moving.side = "right",
                                  with.labels = TRUE, simplify.as = "matrix")
  colnames(up1) <- paste0("Ddown10", colnames(up1))

  up2 <- oligonucleotideFrequency(DNAStringSet(ASdata$Ddown20), width = 3, step = 1,
                                  as.prob = TRUE, as.array = FALSE, fast.moving.side = "right",
                                  with.labels = TRUE, simplify.as = "matrix")
  colnames(up2) <- paste0("Ddown20", colnames(up2))

  down1 <- oligonucleotideFrequency(DNAStringSet(ASdata$Aup10), width = 3, step = 1,
                                    as.prob = TRUE, as.array = FALSE, fast.moving.side = "right",
                                    with.labels = TRUE, simplify.as = "matrix")
  colnames(down1) <- paste0("Aup10", colnames(down1))

  down2 <- oligonucleotideFrequency(DNAStringSet(ASdata$Aup20), width = 3, step = 1,
                                    as.prob = TRUE, as.array = FALSE, fast.moving.side = "right",
                                    with.labels = TRUE, simplify.as = "matrix")
  colnames(down2) <- paste0("Aup20", colnames(down2))

  down3 <- oligonucleotideFrequency(DNAStringSet(ASdata$Aup30), width = 3, step = 1,
                                    as.prob = TRUE, as.array = FALSE, fast.moving.side = "right",
                                    with.labels = TRUE, simplify.as = "matrix")
  colnames(down3) <- paste0("Aup30", colnames(down3))
  row.names(up1) <- c(1:nrow(ASdata))
  row.names(up2) <- c(1:nrow(ASdata))
  row.names(down1) <- c(1:nrow(ASdata))
  row.names(down2) <- c(1:nrow(ASdata))
  row.names(down3) <- c(1:nrow(ASdata))
  ASnuFre <- data.frame(up1, up2, down1, down2, down3)
  row.names(ASnuFre) <- c(1:nrow(ASdata))
  return(ASnuFre)
}


featureASCTD <- function(ASdata) {
  CTD_seq <- featureCTD(ASdata$seq, class = elements("rnaBase"))
  CTD_up1 <- featureCTD(ASdata$Ddown10, class = elements("rnaBase"))
  CTD_up2 <- featureCTD(ASdata$Ddown20, class = elements("rnaBase"))
  CTD_down1 <- featureCTD(ASdata$Aup10, class = elements("rnaBase"))
  CTD_down2 <- featureCTD(ASdata$Aup20, class = elements("rnaBase"))
  CTD_down3 <- featureCTD(ASdata$Aup30, class = elements("rnaBase"))
  colnames(CTD_seq) <- paste0("CTD_seq", colnames(CTD_seq))
  colnames(CTD_up1) <- paste0("CTD_Ddown10", colnames(CTD_up1))
  colnames(CTD_up2) <- paste0("CTD_Ddown20", colnames(CTD_up2))
  colnames(CTD_down1) <- paste0("CTD_Aup10", colnames(CTD_down1))
  colnames(CTD_down2) <- paste0("CTD_Aup20", colnames(CTD_down2))
  colnames(CTD_down3) <- paste0("CTD_Aup30", colnames(CTD_down3))
  row.names(CTD_seq) <- c(1:nrow(ASdata))
  row.names(CTD_up1) <- c(1:nrow(ASdata))
  row.names(CTD_up2) <- c(1:nrow(ASdata))
  row.names(CTD_down1) <- c(1:nrow(ASdata))
  row.names(CTD_down2) <- c(1:nrow(ASdata))
  row.names(CTD_down3) <- c(1:nrow(ASdata))
  ASCTD <- data.frame(CTD_seq, CTD_up1, CTD_up2, CTD_down1, CTD_down2, CTD_down3)
  row.names(ASCTD) <- c(1:nrow(ASdata))
  return(ASCTD)
}



featureDonorMotif <- function(donorSeq, PWMdonor) {
  donorSeq <- as.character(donorSeq)
  donorScore <- rep(0, length(donorSeq))
  for (j in 1:length(donorSeq)) {
    seq <- donorSeq[j]
    x <- strsplit(x = seq, split = "")
    seq_score <- vector()
    for (i in 1:length(x[[1]])) {
      seq_score[i] <- PWMdonor[x[[1]][i], i]
    }
    donorScore[j] <- sum(seq_score)
  }
  return(donorScore)
}


featureAcceptorMotif <- function(acceptorSeq, PWMacceptor) {
  acceptorSeq <- as.character(acceptorSeq)
  acceptorScore <- rep(0, length(acceptorSeq))
  for (j in 1:length(acceptorSeq)) {
    seq <- acceptorSeq[j]
    x <- strsplit(x = seq, split = "")
    seq_score <- vector()
    for (i in 1:length(x[[1]])) {
      seq_score[i] <- PWMacceptor[x[[1]][i], i]
    }
    acceptorScore[j] <- sum(seq_score)
  }
  return(acceptorScore)
}


#' @title Attribute Evaluators.
#'
#' @description  Using function InfoGainAttributeEval from Rweka  to
#'  evaluate the worth of an attribute by measuring the information
#'  gain with respect to the class.
#'
#' @param trainset The training data set.
#' @param num The number of top features based on the entropy value can
#'   be display in figure.
#' @returns A figure about attribute evaluatiors
#' @examples
#' ##Loading example data
#' rice_model<- load(system.file("data","rice_model.Rdata",package = "AStrap"))
#' library(RWeka)
#' library(ggplot2)
#' plotGain(trainset,20)
plotGain <- function(trainset, num) {
  infor <- InfoGainAttributeEval(class ~ ., data = trainset)
  df <- data.frame(features = names(infor), information = infor)
  df <- df[order(df$information, decreasing = TRUE), ]
  # df$features <- factor(df$features,levels=df$features) datasub <- df[(nrow(df)-19):nrow(df),]
  p <- ggplot(data = df[1:num, ], aes(x = features, y = information, fill = features)) +
    geom_bar(stat = "identity") + labs(y = "Information gain",x = NULL) +
    coord_flip() + theme_bw() + guides(fill = FALSE)
  return(p)
}



#' @title Plotting a ROC curve
#'
#' @description Plotting a ROC curve based on ggplots.
#'
#' @param probality A vector, matrix, list, or data frame containing
#'   the predictions
#' @param labels A vector, matrix, list, or data frame containing
#'   the true class labels. Must have the same dimensions as 'predictions
#' @return A ROC curve.
#' @examples
#' ##Loading example data
#' load(system.file("data","example_svm.Rdata",package = "AStrap"))
#' proba <- as.data.frame(attr(example_svm$svmPre, "probabilities"))
#' testset <- example_svm$testSet
#' label <- testset$class
#' library(ggplot2)
#' plotROC(proba,label)
plotROC <- function(probalility, labels) {
  name <- colnames(probalility)
  Plotdata <- data.frame(x = 1, y = 2, type = "aa")
  Auc <- c()
  for (i in 1:length(name)) {
    label <- as.character(labels)
    label[which(labels == name[i])] <- 1
    label[-which(labels == name[i])] <- 0
    pred <- prediction(probalility[, i], label)
    perf <- ROCR::performance(pred, "tpr", "fpr")
    auc <- ROCR::performance(pred, "auc")
    auc <- unlist(slot(auc, "y.values"))
    Auc <- c(Auc, auc)
    x <- unlist(perf@x.values)
    y <- unlist(perf@y.values)
    plotdata <- data.frame(x = x, y = y)
    plotdata$type <- name[i]
    Plotdata <- rbind(Plotdata, plotdata)
  }
  names(Auc) <- name
  Plotdata <- Plotdata[-1, ]
  Auc <- round(Auc, 3)
  names(Plotdata) <- c("fpr", "tpr", "type")
  g <- ggplot(Plotdata, aes(x = fpr, y = tpr, colour = type)) + geom_path(size = 1) +
    scale_color_manual(breaks = c("AltA", "AltD", "ES", "IR"),
                       values = c("lightcoral", "steelblue", "hotpink", "mediumturquoise")) +
    labs(x = "False positive rate FPR=FP/N", y = "Ture positive rate TPR=TP/P") +
    theme(plot.title = element_text(face = "bold", size = 15)) + theme_bw()


  g <- g + annotate("text", x = 0.35, y = 0.85, label = paste("AUC", Auc[which(names(Auc) == "AltA")]), color = "lightcoral", size = 6) +
    annotate("text", x = 0.35, y = 0.75, label = paste("AUC", Auc[which(names(Auc) == "AltD")]), color = "steelblue", size = 6) +
    annotate("text", x = 0.35, y = 0.65, label = paste("AUC", Auc[which(names(Auc) == "ES")]), color = "hotpink", size = 6) +
    annotate("text", x = 0.35, y = 0.55, label = paste("AUC", Auc[which(names(Auc) == "IR")]), color = "mediumturquoise", size = 6) +
    theme(legend.position = "right") + guides(color = guide_legend(title = NULL))
  print(g)
  return(list(AUC = Auc, ROC = Plotdata))

}


#' @title Plotting a AS event.
#'
#' @description Visualization of AS events of different AS types, including AltA, AltD,
#'   ES, IR.
#'
#' @param ASevent A list of AS events.
#' @param id The ID of AS events.
#' @return This function retures a picture of AS event.
#' @examples
#' ##Load example data
#' load(system.file("data","sample_ASresult.Rdata",package = "AStrap"))
#' ##Visualization of AS events of different AS types
#' library(Gviz)
#' plotAS(result$ASevent, id = 1)
#' plotAS(result$ASevent, id = 7)
#' plotAS(result$ASevent, id = 13)
#' plotAS(result$ASevent, id = 21)
plotAS <- function(ASevent, id = 1) {
  Sali <- unlist(str_extract_all(ASevent$Salign[id], "\\d+"))
  alignNum <- length(Sali)
  Slength <- ASevent$Slength[id]
  Astart <- ASevent$Sstart[id]
  Aend <- ASevent$Send[id]
  Qname <- ASevent$Qid[id]
  Sname <- ASevent$Sid[id]
  AStype <- ASevent$prediction[id]
  Seq <- unlist(str_extract_all(ASevent$spliceSeq[id], "\\w+"))
  aTrack1 <- AnnotationTrack(start = c(Sali[seq(1, alignNum, 2)]),
                             end = c(Sali[seq(2, alignNum, 2)]),
                             chromosome = "chr1", strand = "+", col = "gray",
                             shape = c("smallArrow", "arrow"), fill = "deepskyblue1",
                             group = Qname, id = "Query", genome = "hg19", name = "Query")
  aTrack2 <- AnnotationTrack(start = c(1),
                             end = c(Slength),
                             chromosome = "chr1", strand = "+", col = "gray",
                             shape = c("smallArrow","arrow"), fill = "deepskyblue1",
                             group = Sname, id = "Subject", genome = "hg19", name = "Subject")
  aTrack3 <- AnnotationTrack(start = c(Astart, Aend),
                             end = c(Astart, Aend),
                             chromosome = "chr1", strand = "*", col = "gray",
                             shape = c("smallArrow","arrow"), fill = "gray",
                             group = AStype, id = c(Seq[1], Seq[2]), genome = "hg19", name = AStype)
  plotTracks(list(aTrack1, aTrack3, aTrack2), groupAnnotation = "group",
             shape = "box", just.group = "left", fontcolor.group = "black",
             fontsize.group = 15, background.title = "lightcoral", featureAnnotation = "id",
             fontcolor.feature = "deepskyblue1", fontsize.feature = 8)

}

