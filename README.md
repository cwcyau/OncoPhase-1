# OncoPhase
OncoPhase: Quantification of somatic mutation cellular prevalence using phase information.

Developed in the Ovarian Cancer Cell Laboratory, University of Oxford, United Kingdom.

This package offers a direct method to quantify the cellular
prevalence of single nucleotide variants (SNVs) using phase information. The
method utilizes three sources of information: the phasing information, the copy
number variation, and the allele counts. The method is demonstrated to bring
more capabilities in Cancer Genomic and allows computing the cell prevalence of
a mutation in various cancer contexts.


To Download and install the package.

The library devtools offer a command to
quickly install any R package deposited in Github.
To install OncoPhase directly from this repository using the library devtools, process as follow:


1- If you dont already have devtools installed, proceed to its installation : 

install.packages("devtools")

library(devtools)


2- Download and install OncoPhase :

install_github("chedonat/OncoPhase")

library(OncoPhase)

3- Have a look at the package help

help(OncoPhase)




