# MiATDS

Type: Package

Title: A powerful adaptive microbiome-based association test for microbial association signals with diverse sparsity levels
Version: 1.0

Author: Xiaoyun Huang

Maintainer: Xiaoyun Huang <huangxiaoyun@mails.ccnu.edu.cn>

Imports: phyloseq, cluster, compositions, permute, vegan, ape, aSPU, MiSPU, devtools, MiHC 

Description: An adaptive microbiome-based association test for detecting microbial association signals with diverse sparsity levels

License: GPL-2

Encoding: UTF-8

LazyData: true

URL: https://github.com/XiaoyunHuang33/MiATDS

## Introduction

`MiATDS`, an adaptive microbiome-based association test combining the microbiome higher criticism analysis (MiHC) and adaptive weighted sum of powered score tests (aWSPU), has a good performance for detecting microbial association signals with diverse sparsity levels.

## Installation 

phyloseq:
```
BiocManager::install("phyloseq")
```

cluster:
```
install.packages("cluster")
```

compositions:
```
install.packages("compositions")
```

permute:
```
install.packages("permute")
```

vegan:
```
install.packages("vegan")
```

ape:
```
install.packages("ape")
```

aSPU:
```
install.packages("aSPU")
```

MiSPU:
```
install.packages("MiSPU")
```

devtools:
```
install.packages("devtools")
```

MiHC:
```
devtools::install_github("hk1785/MiHC", force=T)
```


You may install `MiATDS` from GitHub using the following code: 

```
devtools::install_github("XiaoyunHuang33/MiATDS", force=T)
```

---------------------------------------------------------------------------------------------------------------------------------------

## Description
The function `MiATDS` tests association between microbiome and a host phenotype. The association signals of association test can be at diverse sparsity levels. And host phenotype must be continuous or binary.


## Usage
```
MiATDS(y, otu.tab, cov=NULL, tree=NULL, model = c("gaussian", "binomial"), pow=c(1:5), comp=FALSE, CLR=FALSE, opt.ncl=30, n.perm=5000)
```

## Arguments
* _y_ - host phenotype of interest (Continuous or binary).
* _otu.tab_ - OTU count table.
* _covs_ - covariate (e.g., age, gender).
* _tree_ - A rooted phylogenetic tree.
* _model_ - "gaussian" for Gaussian outcomes, "binomial" for Binomial outcomes.
* pow - A set of gamma values. Default is c(1:5).
* _comp_ - An indicator if the OTU table contains absolute abundances or relative abundances. Default is comp=FALSE for absolute abundances.
* _CLR_ - The centered log-ratio (CLR) transformation. Default is CLR=FALSE for no CLR transformation.
* _opt.ncl_ - A upper limit to find the optimal number of clusters. Default is opt.ncl=30.
* _n.perm_ - A number of permutations. Default is n.perm=5000. 

## Values
_$pd.rank_ - The ranking for probability degree.

_$awSPU.pvs_  - The p-value for the wSPU test and awSPU test.

_$aWSPU.pvs_  - The p-value for the WSPU test and aWSPU test.

_$omnibus.pvs_  - The p-value for the MiATDS test.


## Example

Import requisite R packages: 

```
library(cluster)
library(permute)
library(phyloseq)
library(MiHC)
library(MiATDS)  
```


Import example microbiome data:

```
data(obesity_data)
otu.tab <- obesity_data@otu_table
tree <- obesity_data@phy_tree
y<- obesity_data@sam_data$label
cov <- as.matrix(obesity_data@sam_data$x1)
```

Fit GEEMiHC:

```
set.seed(123)
out <- MiATDS(y, otu.tab, cov=cov, tree, model = "binomial")
out
```

## References
* Pan, W., Kim, J., Zhang, Y., Shen, X., Wei, P., 2014. A powerful and adaptive association test for rare variants. Genetics 197, 1081-1095.

* Koh, H., Blaser, M.J., Li, H., 2017. A powerful microbiome-based association test and a microbial taxa discovery framework for comprehensive association mapping. Microbiome 5, 45.

* Koh, H., Zhao, N., 2020. A powerful microbial group association test based on the higher criticism analysis for sparse microbial association signals. Microbiome 8, 63.


---------------------------------------------------------------------------------------------------------------------------------------

# Other Resources

## OTUs_simulated function

### Description
We generate the OTUs count table simulated based on the Dirichlet-multinomial model according to real data.

### Usage
```
OTUs_simulated(data, nSam, nOTU, n_repeat, mu, size)
```

### Arguments
* _data_ - real data.

* _nSam_ - Sample size.

* _nOTU_ - The number of OTUs.

* _n\_repeat_ - The number of repeat.

* _mu_ - The mean of the negative binomial distribution.

* _size_ - The size of the negative binomial distribution.


### Values

_$OTU\_simulated_ - OTU counts table simulated based on real data.

### Example

```
data("throat.otu.tab", package ="MiSPU")
otu.tab <- round(OTUs_simulated(data=throat.otu.tab, nSam=100, nOTU=100, n_repeat=10, mu=1000, size=25)$OTU_simulated)
```

### References

* Wu, C., Chen, J., Kim, J., Pan,W., 2016. An adaptive association test for microbiome data. Genome Med. 8, 56.

* La Rosa, P.S., Brooks, J.P., Deych, E., Boone, E.L., Edwards, D.J., Wang, Q., Sodergren, E., Weinstock, G., Shannon, W.D., 2012. Hypothesis Testing and Power Calculations for Taxonomic-Based Human Microbiome Data. PLoS ONE 7, e52078.

## The datasets of case applications

(1) Association analysis between gut microbiome and obesity：

the information of processed data: ```data(obesity_data)```.

(2) Association analysis between gut microbiome and colorectal cancer:

the information of processed data: ```data(crc_data)```.

(3) Association analysis between gut microbiome and autism:
 
the information of processed data: ```data(autism_data)```.


### References

* Escobar, J.S., Klotz, B., Valdes, B.E., Agudelo, G.M., 2014. The gut microbiota of Colombians differs from that of Americans, Europeans and Asians. BMC Microbiol. 14, 311.

* Wang, T., Cai, G., Qiu, Y., Fei, N., Zhang, M., Pang, X., Jia, W., Cai, S., Zhao, L., 2012. Structural segregation of gut microbiota between colorectal cancer patients and healthy volunteers. ISME J. 6, 320-329.

* Kang, D.-W., Park, J.G., Ilhan, Z.E., Wallstrom, G., LaBaer, J., Adams, J.B., Krajmalnik-Brown, R., 2013. Reduced incidence of Prevotella and Other fermenters in intestinal microflora of autistic children. PLoS ONE 8, e68322.

## Statement

Our code mainly refers to R packages, _MiHC_, _OMiAT_, and _MiSPU_, where _MiATDS_ function refers to _MiHC_ and _OMiAT_, the generation of OTU table (i.e., _OTUs\_simulated_ function) refers to _MiSPU_.


## Citation

If you use this code for you research, please cite our paper.
```
@article{
title={A powerful adaptive microbiome-based association test for microbial association signals with diverse sparsity levels},
author={Han Sun, Xiaoyun Huang, Lingling Fu, Ban Huo, Tingting He, Xingpeng Jiang},
journal={Journal of Genetics and Genomics},
volume={48(9)},
pages={851–859},
year={2021},
url={https://linkinghub.elsevier.com/retrieve/pii/S1673852721002599},
doi = {10.1016/j.jgg.2021.08.002}
}
```
