
#Set your working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim")
home.dir <- getwd()


#Initiate the objects of: vc.traits, vc.genetic.values, and vc.breeding.values
vc.traits <- list()
vc.genetic.values <- list()
vc.breeding.values <- list()

count <- 1
for(rep in 1:3){
  for(i in c(1,2,4)){
    for(j in c(0.05,0.5,1,2)){
       for(k in c(0.05,0.5,1,2)){
         for(el in c(0.05,0.5,1,2)){
           #Run the function
           factor.A <- i
           factor.B <- j
           factor.C <- k
           factor.D <- el 
           this.rep <- rep
           #Read in the appropriate .Rdata file
           load(paste(i,".FactorA..",j, ".FactorB..",
                      k,".FactorC..",el,".FactorD..",
                      rep, ".Rep.Rdata", sep = ""))
           #Append the appropriate objects to the list
           vc.traits[[count]] <- these.trait.var.covar
           vc.genetic.values[[count]] <- these.genetic.value.var.covar
           vc.breeding.values[[count]] <- these.breeding.value.var.covar
           count <- count+1
         }#End for(rep in 1:3)
      }#End for(el in 1:4)
    }#End for(k in 1:4)
  }#End for(j in 1:4)
}#End for(i in 1:3)

#End product: 
# vc.traits: object containing every single variance-covariance matrix of every single trait
# vc.genetic.values: object containing every single variance-covariance matrix of every single genetic value
# vc.breeding.values: object containing every single variance-covariance matrix of every single breeding value
