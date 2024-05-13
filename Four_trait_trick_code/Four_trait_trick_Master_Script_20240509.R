#Set the working directory

#The first time you run this, please be sure to install
# AlphaSimR and simplePHENOTYPES



#############################################################################################
##############################################################
# The steps immeidately below are temporary. I envision that this would be ran
# as a function, and applied to every single trait that is simulated in
# every single popluation

#Set your working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project/")
home.dir <- getwd()


#Open all libraries
library(package = "AlphaSimR")
library(package = "simplePHENOTYPES")

#Read in a simulated rep so that you have something to work with
load("2.FactorA..0.05.FactorB..0.05.FactorC..0.05.FactorD..3.Rep.Rdata")

#Obtain the parameters used by simplePHENOTYPES
source("Simulating_Omnigenic_Genetic_Architecture/setInputParametersforSimplePhentoypes.R")

#Obtain the genotypic data
this.input.SNPs <- hapmap.file.of.founder.QTLs




#################################################
#Partition out the which loci and effects for each of the four traits
#Trait 1: the whole thing
#Trait 2: only core effects
#Trait 3: only peripheral effects
#Trait 4: only core-by-peripheral epistatic effects



#Simulate the four traits in simplePHENOTYPES


#Simulate the genetic values of the four traits in simplePHENOTYPES


#Use popVar to calculate the variance-covariance matrix of the traits

#Use popVar to calculate the variance-covariance matrix of the genetic values


