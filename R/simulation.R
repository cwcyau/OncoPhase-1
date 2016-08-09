


#' The Phylogeny class
#'
#' An S4 base class representing a phylogeny
#' @slot Name Object of class \code{\link{character}} representing the name of the phylogeny
#' @slot NbSNV Object of class integer, attribute of the class Phylogeny representing the number of single nucleotide variations (SNVs).
#' @slot NbSCNAs Object of class integer, attribute of the class Phylogeny representing the number of somatic copy number alterations (SCNAs)
# #' @slot NbSNVClusters Object of class integer, attribute of the class Phylogeny representing the number of SNVs clusters.
#' @slot snv_ids Object of class list, attribute of the class Phylogeny representing the identifiers or names of the NbSNVs SNVs clusters.
# #' @slot snvclusters_ids Object of class list, attribute of the class Phylogeny representing the identifiers or names of the NbSNVClusters SNV clusters.
# #' @slot snv_clutsers Object of class list, attribute of the class Phylogeny giving the clusters assigned to each SNV
#' @slot scna_list Object of class list, attribute of the class Phylogeny representing the list of SCNA given in the form list(scna_1_name=scna_1_attibutes, scna_2_name=scna_2_attibutes,...,scna_NbSCNAs_name=scna_NbSCNAs_attibutes ). 
#' Each SCNA attribute is a list containing 2 fields :
#' \describe{
#' \item{CN}{A pair of integer representing the major and minor copies numbers of the SCNA given in the form c(major, minor)}
#' \item{LOC}{ The location or genomic region affected by the SCNA represented in term of the list of SNVs spanned by the SCNA. LOC is given inform of a vector (a_1, a_2, ..., a_NbSNV\]). Each  a_i takes values in {0,1,2,3,4} as follow:
#' \describe{
#' \item{a_i=0}{the scna do not span the locus of the SNV i}
#' \item{a_i=1}{the SCNA span the SNV i locus, and the SNV is harbored by all the copies of the major copy number chromosome}
#' \item{a_i=2}{the SCNA span the SNV i locus, and the SNV is harbored by all the copies of the minor copy number chromosome}
#' \item{a_i=3}{ the SCNA span the SNV i locus, and the SNV is harbored by one copy of the major copy number chromosome}
#' \item{a_i=4}{the SCNA span the SNV i locus, and the SNV is harbored by one copy of the minor copy number chromosome}
#' }
#' }
#' }
#' @slot Clones Object of class list, attribute of the Class Phylogeny containing the list of the clones given in the form list(clone_1_name=clone_1_attibutes, clone_2_name=clone_2_attibutes,...,clone_NbClones_name=clone_NbClones_attibutes). The Germline clones do not need to be list, it will be deduced. Each clone's attribute is a list containing the following fields :
#' \describe{
#' \item{snv}{ID list of the SNVs harbored by cells of the clone}
#' \item{scna}{ID list of the SCNAs affecting the  cells of the clone}
#' \item{prev}{Cellular prevalence of the clone}
#' }
#' 
#' 
#' @seealso \code{\link{simulation.Phylogeny}}, 
#' @export Phylogeny
#' @exportClass Phylogeny
#' 
Phylogeny <- setClass("Phylogeny", 
                         slots = c(
                           Name = "character",
                           NbSNVs = "numeric",
                           NbSCNAs = "numeric",
                           NbSNVClusters = "numeric",
                           snv_ids = "character",
                          # snvclusters_ids = "list",
                          # snv_clusters = "list",
                           scna_list = "list",
                           Clones = "list"),
                         prototype = list(Name = "NormalSample",
                                          NbSNV = 0,
                                          NbSCNAs = 0,
                                          NbSNVClusters = 0,
                                          snv_ids = NULL,
                                         # snvclusters_ids = list(),
                                         # snv_clusters = list(),
                                          scna_list = NULL,
                                          Clones = NULL
                                          ),
                         validity=function(object)
                         {
                           is.wholenumber <-
                             function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
                           

                           if(!is.wholenumber(NbSNVs) || !is.wholenumber(NbSCNAs) #|| !is.integer(NbSNVClusters) 
                              ||NbSNVs<0 || NbSCNAs <0 #|| NbSNVClusters <0
                              )
                             return("The attributes  NbSNVs and NbSCNA should all be positive integer. Check their values ")
                           
                           if(length(snv_ids) != NbSNVs)
                             return(NbSNVs, "SNV IDs expected but ",  lenght(snv_ids) , "provided. Check your parameter snv_ids")
                           #if(length(snvclusters_ids) != NbSNVClusters)
                          #   return(NbSNVClusters, "SNV IDs expected but ",  lenght(snvclusters_ids) , "provided. Check your parameter snvclusters_ids")
                         #  if(length(snv_clusters) != NbSNVs)
                          #   return(NbSNVs, "clusters IDs expected but ",  lenght(snv_clutsers) , "provided. Check your parameter snv_clutsers")
                           if(length(scna_list) != NbSCNAs)
                             return(NbSCNAs, "SNV IDs expected but ",  lenght(scna_list) , " provided. Check your parameter snv_ids")
                           
                         #  if(!(unique(snv_clusters) %in% snvclusters_ids))
                         #     return(setdiff(unique(snv_clusters), snvclusters_ids), " are unknown cluster IDs")
                              
                           for(scna in scna_list){
                             if(!("CN" %in% names(scna)))
                               return("Parameter CN (copy numbers) is missing for this SCNA.")
                             if(!("LOC" %in% names(scna)))
                               return("Parameter LOC (Genomic Location) is missing for this SCNA.")
                           }
                           
                           for(clone in Clones){
                             if(sum(!(names(clone) %in% c("snv","scna","prev"   )))>1)
                               return(paste("Unknown parameters to clones : ", setdiff(names(clone), c("snv","scna","prev"   ))))
                             
                            # if(!("scna" %in% names(clone)))
                           #    return("Parameter scna is missing for this clone.")
                             if(!is.null(clone$snv))
                             if(sum(!(clone$snv %in% snv_ids ))>0)
                               return( "Unknown snv provided :",setdiff(cluster,snvclusters_ids  ))
                             if(!is.null(clone$scna))
                             if(sum(!(clone$scna %in% names(scna_list) ))>0)
                               return( "Unknown scna provided :",setdiff(cna,names(scna_list)   ))
                           }
                           
                           

                         }
)








#Examples
#Examples

#load a configuration
phylogeny="phylogeny11"
#Number of SNVs and Number of SCNAs
NbSNVs= 5
NbSCNAs= 2
#NbSNVClusters = 5
snv_ids = c("M1","M2","M3","M4","M5")
#snvclusters_ids = c("M1","M2","M3","M4","M5")
#snv_clusters = c("M1","M2","M3","M4","M5")
scna_list = list(
  "M6"=list(CN=c(2,1),LOC=c(0,0,1,0,0) ),
  "M7"=list(CN=c(2,0),LOC=c(0,0,0,0,1))
)
Clones = list(
  #  "Germline" = list(snv=c(0,0,0,0,0),0p=0.0 ),
  "CloneA" = list(snv="M1",
                  prev=0.1),
  "CloneB" = list(snv=c("M1","M2"),
                  prev=0.3),
  "CloneC" = list(snv=c("M1","M2","M3"),
                  prev=0.1),
  "CloneD" = list(snv=c("M1","M2","M3"),
                  scna="M6",
                  prev=0.15),
  "CloneE" = list(snv=c("M1","M2","M3","M4"),
                  scna="M6",
                  prev=0.15),
  "CloneF" = list(snv=c("M1","M2","M5"),
                  scna="M7",
                  prev=0.20)
)




phylogeny11=Phylogeny(Name="phylogeny11",NbSNVs=NbSNVs,NbSCNAs=NbSCNAs, snv_ids=snv_ids, scna_list=scna_list, Clones=Clones)



