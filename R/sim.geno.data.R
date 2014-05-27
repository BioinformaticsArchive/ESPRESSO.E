#' 
#' @title Generates genotypes for a genetic variant
#' @description Generates two alleles and combines them to form the genotype of 
#' a SNP under a binary or additive genetic model.
#' @param num.obs number of observations to generate.
#' @param geno.model genetic model; binary=0 and additive=1.
#' @param MAF minor allele frequency of the SNP (in ESPRESSO this is the frequency of the 'at risk' allele)
#' @return a dataframe that contains the following variables:
#' \code{allele.A} major allele
#' \code{allele.B} minor allele
#' \code{genotype} genotype
#' @keywords internal
#' @author Gaye A.
#' 
sim.geno.data <- function(num.obs=10000, geno.model=0, MAF=0.1){
  
  numobs <- num.obs
  geno.mod <- geno.model
  geno.maf <- MAF

  # CORRECTION TERM FOR MEAN CENTERING FOR ADDITIVE 
  mean.geno.additive <- (2*geno.maf*(1-geno.maf)+2*(geno.maf^2))

  # CORRECTION TERM FOR MEAN CENTERING FOR BINARY GENE
  mean.geno.binary <- (2*geno.maf*(1-geno.maf)+(geno.maf^2))

  # CREATE, CENTRE AND ROUND AN ADDITIVE GENETIC GENOTYPE COVARIATE WITH
  # APPROPRIATE MAF 
  allele.A <- rbinom(numobs,1,geno.maf)
  allele.B <- rbinom(numobs,1,geno.maf)

  # ACTUAL GENOTYPE IS SUM OF ALLELES IN ADDITIVE GENETIC MODEL
  # AND IS 1 IF SUM OF ALLELES IS 1 OR GREATER IN THE BINARY MODEL 
  genotype <- allele.A+allele.B

		if(geno.model==0){
				geno.U <- genotype > 0
				genotype <- geno.U - mean.geno.binary
		}
		if(geno.model==1){
				geno.U <- genotype
				genotype <- geno.U - mean.geno.additive
		}

  # return the data as a dataframe
  dtframe <- data.frame(allele.A, allele.B, genotype)
}

