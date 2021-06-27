\name{OTUs_simulated}
\alias{OTUs_simulated}
\title{OTUs_simulated}
\usage{
OTUs_simulated(data, nSam, nOTU, n_repeat, mu, size)
}
\description{
We used a phylogenetic tree of OTUs from a real throat microbiome data set, which consists of 856 OTUs after discarding singleton OTUs. Then based on a complicated statistical meodel, we generated the OTU counts for each individual to simulate the features of real microbiome data.
}

\arguments{
  \item{data}{
Real data.
}
  \item{nSam}{Sample size.

}
  \item{nOTU}{
The number of OTUs.
}
  \item{n_repeat}{
The number of repeat.

}
  \item{mu}{
The mean of the negative binomial distribution.
}
  \item{size}{
The size of the negative binomial distribution.
}

}

\value{
OTU_simulated: OTU counts table simulated based on real data.

}



\examples{
# Import requisite R packages
library(MiSPU)

#Simulated OTUs
data("throat.otu.tab")
otu.tab <- round(OTUs_simulated(data=throat.otu.tab, nSam=100, nOTU=100, n_repeat=10, mu=1000, size=25)$OTU_simulated)

}
