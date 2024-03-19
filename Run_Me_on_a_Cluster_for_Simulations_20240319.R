#Set the working directory

#The first time you run this, please be sure to install
# AlphaSimR and simplePHENOTYPES

####Begin temporary code
this.factor.A <- 2
this.factor.B <- 0.05
this.factor.C <- 0.05
this.factor.D <- 0.05
this.rep <- 2
this.home.dir<- getwd()



factor.A <- this.factor.A 
factor.B <- this.factor.B 
factor.C <- this.factor.C
factor.D <- this.factor.D 
rep <- this.rep
home.dir <- this.home.dir
this.output.directory.name <-  paste("Factor.A.", factor.A,
                                      "Factor.B.", factor.B,
                                      "Factor.C.", factor.C,
                                      "Factor.D.", factor.D,
                                      "Rep.", this.rep,
                                      sep = "")


#Set your working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/Sabbatical_Project")
this.home.dir.name <- getwd()


#Open all libraries
library(package = "AlphaSimR")
library(package = "simplePHENOTYPES")

#Source in all files
source("Simulating_Omnigenic_Genetic_Architecture/Simulate_Omnigenic_Genetic_Architecture_as_a_Function_20240219.R")
source("Simulating_Omnigenic_Genetic_Architecture/Simulate_Multiple_Generations_of_Selection.R")
source("Functions_to_Make_Life_Easier/get.me.my.SNPs.in.hapmap.format.R")

#Simulate a founder population in AlphaSimR
source("Simulate_a_founder_Population_20240319.R")

###For loop through all settings of the factorial experiment
### The end product should be a whole bunch of R workspaces
### Not to mention output directories of simulated traits
for(i in c(1,2,4)){
  for(j in c(0.05,0.5,1,2)){
    for(k in c(0.05,0.5,1,2)){
       for(el in c(0.05,0.5,1,2)){
         for(rep in 1:3){
           #Run the function
           this.factor.A <- i
           this.factor.B <- j
           this.factor.C <- k
           this.factor.D <- el
           this.rep <- rep
           one.replicate.of.simulation.run(factor.A = this.factor.A, 
                                           factor.B = this.factor.B,
                                           factor.C = NULL, 
                                           factor.D = NULL,
                                           rep = NULL,
                                           home.dir.name = NULL)
         }#End for(rep in 1:3)
      }#End for(el in 1:4)
    }#End for(k in 1:4)
  }#End for(j in 1:4)
}#End for(i in 1:3)

these.factorial.levels <- data.frame(factor.A, factor.B, factor.C, factor.D)

setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/Housekeeping_notes_for_factorial_Experiment")

write.csv(these.factorial.levels, "Tally_of_Experimental_Levels_Ran.csv", row.names = FALSE)
