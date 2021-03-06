% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OncoPhase_methods.R
\name{getPrevalenceMultiSamples}
\alias{getPrevalenceMultiSamples}
\title{Somatic mutations cellular prevalence using haplotype phasing on a multi sample study.}
\usage{
getPrevalenceMultiSamples(snp_allelecount_df, ref_allelecount_df,
  major_copynumber_df, minor_copynumber_df, mode = "PhasedSNP",
  cnv_fraction = NULL, phasing_association_df = NULL,
  NormalcellContamination_df = NULL, tumoursamples = NULL,
  nbFirstColumns = 3, region = NULL, detail = TRUE, LocusRadius = 10000,
  SameTumour = TRUE, LocusCoverage = 1, ProgressOutputs = T)
}
\arguments{
\item{snp_allelecount_df}{A data frame containing for each mutation the  allelic 
counts of the variant at each tumor samples. The data frame should contains at least the following three columns among its firsts columns: Chrom (The mutation
 chromosome) , End (The mutation position) and IsGermline (is the mutation a germline
  or somatic mutation).}

\item{ref_allelecount_df}{A data frame containing for each mutation the allelic count
of the reference at each tumor sample. The data frame should contains at least the following three columns among its firsts columns:  Chrom (The mutation
  or Somatic mutation)}

\item{major_copynumber_df}{A data frame containing for each mutation, its  major}

\item{minor_copynumber_df}{A data frame containing for each mutation the minor
chromosomal copy number at each tumor samples.}

\item{mode}{The mode under which the prevalence is computed  (default : PhasedSNP , alternatives methods  are FlankingSNP, OptimalSNP,and SNVOnly).  Can also be provided as a numeric 0=SNVOnly, 1= PhasedSNP, 2=FlankingSNP and 3 = OptimalSNP
#@param formula The formula used to compute the prevalence. can be either "matrix" for the linear equations or "General" for the exact allele count cases. Default : Matrix}

\item{phasing_association_df}{A data frame containing for each somatic mutation, 
a colon separated list of germline SNP phased to it.}

\item{tumoursamples}{: The list of tumor samples to consider for the prevalence
  snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df
  and  CNVfraction_df. If not provided, the headers from nbFirstColumns + 1 to 
  the last column of snp_allelecount_df is retrieved and its intersection with the
other inputted data frames headers is considered.}

\item{nbFirstColumns}{Number of first columns in snp_allelecount_df to reproduce in
the output dataframe e.g: Chrom, Pos, Vartype. Columns from  nbFirstColumns +1 to the last column should contains the information needed for the prevalence computation at each tumour sample}

\item{region}{The region of the genome to consider for the prevalence computation  in the format chrom:start-end 
e.g "chr22:179800-98767}

\item{detail}{when set to TRUE, a detailed output is generated containing, the context and the detailed prevalence for each group of cells (germline cells, cells affected by one of the two genomic alterations (SNV or CNV) but not both, cells affected by  by both copynumber alteration and SNV ). Default : TRUE.}

\item{LocusCoverage}{when set to true, lambda_S and mu_S might be adjusted if necessary so that they meet the rules lambda_S <= lambda_G. mu_S >= mu_G and lambda_S + mu_S = lambda_G + mu_G. Not used if mode=SNVOnly,  Default = FALSE}

\item{CNVfraction_df, }{If provided, represents a data frame containing for each mutation,  the fraction of
cells affected by the copy number alteration. Used only if the mode is "PhasedSNP" and formula is "General".}

\item{min_cells}{Minimum number of cells (default 2). In case the estimated number of cells sequenced at the locus of the mutation is less than min_cells, NA is returned.}

\item{min_alleles}{Minimum number of alleles. (default 4). In case the estimated number of alleles sequenced at the locus of the mutation is less than min_alleles, NA is returned.}
}
\value{
A data frame containing :
 \describe{
       \item{}{Column 1 to NbFirstcolumn of the input data frame snp_allelecount_df. 
       This will generally include the chromosome and the position of the mutation plus
       any other columns to report in the prevalence dataframe (e.g REF, ALL, ...) }
        \item{}{One column per tumour sample reporting the prevalence of the mutation 
        at each samples}
     }
}
\description{
This is a generic function to compute the cellular prevalence of somatic mutations in
 cancer using haplotype phasing.  The function applies the model to a range of mutations located at a given genomic region or at the whole genome scale. The model computes the prevalence of a somatic
  mutation relatively to close and eventually phased germline mutations. It uses three sources
   of information as input : The allelic counts, the phasing information and the 
   copy number alteration.  Multiple tumor samples can be provided for the prevalence computation.
}
\examples{

#Example 1: Loading a simple example data set with two somatic mutations, 5 germlines SNP, and 3 tumor samples
data(simpleExample2)
se=simpleExample2
prevalence_df=getPrevalenceMultiSamples(se$snp_allelecount_df, se$ref_allelecount_df, se$major_copynumber_df,se$minor_copynumber_df,phasing_association_df=se$phasing_association_df, )
print(prevalence_df)

#Chrom     End IsGermline  Tumour1        Tumour2        Tumour3
#mutation2  chr2 3003000          0 C2:0|0|1 C2:0.15|0|0.85 C2:0.12|0|0.88
#mutation6  chr2 4008000          0 C1:1|0|0       C1:1|0|0 C2:0|0.24|0.76

#Example 2 : Computing somatic mutation cellular prevalence on chromosome 15 of  patient 11152 (data retrieved from a parallel study)

data("chr22_1152")
ds=chr22_11152
masterprevalence_df=getPrevalenceMultiSamples(ds$snp_allelecount_df, ds$ref_allelecount_df,  ds$major_copynumber_df,ds$minor_copynumber_df,phasing_association_df = ds$phasing_association_df, cnv_fraction=ds$CNVFraction_df,nbFirstColumns=6,detail=FALSE)
print(head(masterprevalence_df))

data("chr10_11152")
df=chr10_11152
masterprevalence_df=getPrevalenceMultiSamples(df$snp_allelecount_df, df$ref_allelecount_df, df$major_copynumber_df,df$minor_copynumber_df,phasing_association_df=df$phasing_association_df, cnv_fraction=df$CNVFraction_df,nbFirstColumns=6, region="chr10:50000000-180000000")
print(head(masterprevalence_df))



}
\seealso{
\code{\link{getPrevalence}}
}

