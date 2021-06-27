\name{MiATDS}
\alias{MiATDS}
\title{MiATDS}
\usage{
MiATDS(y, otu.tab, tree=NULL, cov=NULL, model = c("gaussian", "binomial"), pow=c(1:5), comp=FALSE, CLR=FALSE, opt.ncl=30, n.perm=5000)
}
\description{
 This function "MiATDS" is an adaptive microbiome-based association test integrating the microbiome higher criticism analysis (MiHC) and adaptive weighted sum of powered score tests (aWSPU). MiATDS has a good performance for detecting microbial association signals with diverse sparsity levels.
}

\arguments{
  \item{y}{
 Host phenotype of interest (Continuous or binary).
}

  \item{otu.tab}{
 A matrix of OTU count table.
}
  \item{tree}{
 A rooted phylogenetic tree.
}

\item{cov}{
 Covariate (e.g., age, gender).
}

  \item{model}{
 "gaussian" for Gaussian outcomes, "binomial" for Binomial outcomes.
}

  \item{comp}{
 An indicator if the OTU table contains absolute abundances or relative abundances. Default is comp=FALSE for absolute abundances.
}
  \item{CLR}{
 The centered log-ratio (CLR) transformation. Default is CLR=FALSE for no CLR. transformation.
}
  \item{opt.ncl}{
 A upper limit to find the optimal number of clusters. Default is opt.ncl=30.
}
  \item{n.perm}{
 A number of permutations. Default is n.perm=5000.
}
}

\value{
  pd.rank: The ranking for probability degree.



  awSPU.pvs: The p-value for the wSPU test and awSPU test.

  aWSPU.pvs: The p-value for the WSPU test and aWSPU test.

  omnibus.pvs: The p-value for the MiATDS test.

}



\examples{
# Import requisite R packages
library(cluster)
library(compositions)
library(phyloseq)
library(permute)
library(MiHC)


# Import example microbiome data

load("D:/R_pro/MiATDS/data/obesity_data.rda")
otu.tab <- obesity_data@otu_table
tree <- obesity_data@phy_tree
y<- obesity_data@sam_data$label
cov <- as.matrix(obesity_data@sam_data$x1)


# Fit MiATDS
set.seed(123)
out <- MiATDS(y, cov=cov, otu.tab=otu.tab, tree=tree, model="binomial", n.perm=5000)
out

}
