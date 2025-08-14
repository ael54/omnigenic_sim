



factor.A <- 1
factor.B <- 0.05
factor.C <- 0.05
factor.D <- 0.05 
this.rep <- 1

i <- 1
j <- 0.05
k <- 0.05
el <- 0.05 
rep <- 1

this.setting <-  paste(".Factor.A.", i,
                       ".Factor.B.", j,
                       ".Factor.C.", k,
                       ".Factor.D.", el,
                       ".Rep.", this.rep,
                       sep = "")
vc.traits.this.setting <- these.trait.var.covar
vc.genetic.values.this.setting <- these.genetic.value.var.covar
vc.breeding.values.this.setting <- these.breeding.value.var.covar


trait.1.variance <- NULL
trait.2.variance <- NULL
trait.3.variance <- NULL
trait.4.variance <- NULL
genetic.value.1.variance <- NULL
genetic.value.2.variance <- NULL
genetic.value.3.variance <- NULL
genetic.value.4.variance <- NULL
breeding.value.1.variance <- NULL
breeding.value.2.variance <- NULL
breeding.value.3.variance <- NULL
breeding.value.4.variance <- NULL
factor.A.vector <- NULL
factor.B.vector <- NULL
factor.C.vector <- NULL
factor.D.vector <- NULL
rep.vector <- NULL
subpopulation.vector <- NULL
selection.number.vector <- NULL
generation.vector <- NULL

#Get info for the founder population
#Initiate everything for the founding population
factor.A.vector <- c(factor.A.vector, factor.A)
factor.B.vector <- c(factor.B.vector, factor.B)
factor.C.vector <- c(factor.C.vector, factor.C)
factor.D.vector <- c(factor.D.vector, factor.D)
rep.vector <- c(rep.vector, this.rep)  
subpopulation.vector <- c(subpopulation.vector, "Founder")
selection.number.vector <- c(selection.number.vector, "Founder")
generation.vector <- c(generation.vector, 0)


this.trait.var.covar <- matrix(unlist(vc.traits.this.setting[which((grepl("Founder.pop", names(vc.traits.this.setting)))&
                                                                    (grepl(paste("Rep.", this.rep,sep = ""), names(vc.traits.this.setting))))]),nrow = 4)
this.genetic.value.var.covar <- matrix(unlist(vc.genetic.values.this.setting[which((grepl("Founder.pop", names(vc.traits.this.setting)))&
                                                                                     (grepl(paste("Rep.", this.rep,sep = ""), names(vc.traits.this.setting))))]),nrow = 4)
this.breeding.value.var.covar <- matrix(unlist(vc.breeding.values.this.setting[which((grepl("Founder.pop", names(vc.traits.this.setting)))&
                                                                                       (grepl(paste("Rep.", this.rep,sep = ""), names(vc.traits.this.setting))))]),nrow = 4)

trait.1.variance <- c(trait.1.variance, this.trait.var.covar[1,1])
trait.2.variance <- c(trait.2.variance, this.trait.var.covar[2,2])
trait.3.variance <- c(trait.3.variance, this.trait.var.covar[3,3])
trait.4.variance <- c(trait.4.variance, this.trait.var.covar[4,4]) 
genetic.value.1.variance <- c(genetic.value.1.variance, this.genetic.value.var.covar[1,1])
genetic.value.2.variance <- c(genetic.value.2.variance, this.genetic.value.var.covar[2,2])
genetic.value.3.variance <- c(genetic.value.3.variance, this.genetic.value.var.covar[3,3])
genetic.value.4.variance <- c(genetic.value.4.variance, this.genetic.value.var.covar[4,4])
breeding.value.1.variance <- c(breeding.value.1.variance, this.breeding.value.var.covar[1,1])
breeding.value.2.variance <- c(breeding.value.2.variance, this.breeding.value.var.covar[2,2])
breeding.value.3.variance <- c(breeding.value.3.variance, this.breeding.value.var.covar[3,3])
breeding.value.4.variance <- c(breeding.value.4.variance, this.breeding.value.var.covar[4,4])        



for(this.selection.type in c("Direct","Disruptive","Stabilizing")){
 
  
  for(this.selection.number in c("112", "224")){
    #vc.traits.this.subpop <- vc.traits.this.setting[which(grepl(this.selection.type, names(vc.traits.this.setting)))]
    vc.traits.this.subpop <- vc.traits.this.setting[which((grepl(this.selection.type, names(vc.traits.this.setting))) & 
                                                          (grepl(this.selection.number, names(vc.traits.this.setting))) )]
    vc.genetic.values.this.subpop <- vc.genetic.values.this.setting[which((grepl(this.selection.type, names(vc.traits.this.setting))) & 
                                                                            (grepl(this.selection.number, names(vc.traits.this.setting))) )]  
    vc.breeding.values.this.subpop <- vc.breeding.values.this.setting[which((grepl(this.selection.type, names(vc.traits.this.setting))) & 
                                                                             (grepl(this.selection.number, names(vc.traits.this.setting))) )]  
    
    
    for(this.generation in 1:10){
      
      factor.A.vector <- c(factor.A.vector, factor.A)
      factor.B.vector <- c(factor.B.vector, factor.B)
      factor.C.vector <- c(factor.C.vector, factor.C)
      factor.D.vector <- c(factor.D.vector, factor.D)
      rep.vector <- c(rep.vector, this.rep)  
      subpopulation.vector <- c(subpopulation.vector, this.selection.type)
      selection.number.vector <- c(selection.number.vector, this.selection.number)
      generation.vector <- c(generation.vector, this.generation)
      
      this.trait.var.covar <- matrix(unlist(vc.traits.this.subpop[which(grepl(paste("Gen.",this.generation, ".Factor", sep = ""), names(vc.traits.this.subpop)))]),nrow = 4)
      this.genetic.value.var.covar <- matrix(unlist(vc.genetic.values.this.subpop[which(grepl(paste(this.generation), names(vc.genetic.values.this.subpop)))]),nrow = 4)
      this.breeding.value.var.covar <- matrix(unlist(vc.breeding.values.this.subpop[which(grepl(paste(this.generation),names(vc.breeding.values.this.subpop)))]),nrow = 4)
      
      trait.1.variance <- c(trait.1.variance, this.trait.var.covar[1,1])
      trait.2.variance <- c(trait.2.variance, this.trait.var.covar[2,2])
      trait.3.variance <- c(trait.3.variance, this.trait.var.covar[3,3])
      trait.4.variance <- c(trait.4.variance, this.trait.var.covar[4,4]) 
      genetic.value.1.variance <- c(genetic.value.1.variance, this.genetic.value.var.covar[1,1])
      genetic.value.2.variance <- c(genetic.value.2.variance, this.genetic.value.var.covar[2,2])
      genetic.value.3.variance <- c(genetic.value.3.variance, this.genetic.value.var.covar[3,3])
      genetic.value.4.variance <- c(genetic.value.4.variance, this.genetic.value.var.covar[4,4])
      breeding.value.1.variance <- c(breeding.value.1.variance, this.breeding.value.var.covar[1,1])
      breeding.value.2.variance <- c(breeding.value.2.variance, this.breeding.value.var.covar[2,2])
      breeding.value.3.variance <- c(breeding.value.3.variance, this.breeding.value.var.covar[3,3])
      breeding.value.4.variance <- c(breeding.value.4.variance, this.breeding.value.var.covar[4,4])
      
      
    }#End for loop through the generations:for(this.generation in 1:10)
  }#End for loop through the number of individuals to select
} #End for loop through the different types of selection: for(this.selection type in c("Direct","Disruptive","Stabilizing"))

data.for.vc.scatterplot <- data.frame(trait.1.variance, trait.2.variance,
                                      trait.3.variance, trait.4.variance,
                                      genetic.value.1.variance, genetic.value.2.variance,
                                      genetic.value.3.variance, genetic.value.4.variance,
                                      breeding.value.1.variance, breeding.value.2.variance,
                                      breeding.value.3.variance, breeding.value.4.variance,
                                      factor.A.vector, factor.B.vector, factor.C.vector,
                                      factor.D.vector, rep.vector,subpopulation.vector, selection.number.vector, 
                                      generation.vector)


#Temporary code to make one plot for one level of one factor

this.factor <- which(grepl("factor",colnames(data.for.vc.scatterplot)))[1]
this.level <- unique(data.for.vc.scatterplot[,this.factor])[1]

#This will be inside of the loops that goes through each factor and each factor levels
data.for.vc.scatterplot.this.factor.level <- data.for.vc.scatterplot[which(data.for.vc.scatterplot[,this.factor] == this.level),]

pdf(paste("VC_Plots_20240820_Temp/VC.Factor",substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),"Level.",this.level,".plot.trait.genetic.breeding.value.three.selection.types.pdf", sep = ""), width = 50)
#For loop through the different kind of selection/subpopulation levels; index on this.selection.type 
for(this.selection.type in unique(data.for.vc.scatterplot.this.factor.level$subpopulation.vector)[-1]){
  #- there will be a different page for each selection/subpopulation level ####
  for(this.selection.number in unique(data.for.vc.scatterplot.this.factor.level$selection.number.vector)[-1]){
    #Source in the code below that will help make these plots
    input.scatter.plot.data <- data.for.vc.scatterplot.this.factor.level[which((data.for.vc.scatterplot.this.factor.level$subpopulation.vector == "Founder")|
                                                                                 ((data.for.vc.scatterplot.this.factor.level$subpopulation.vector == this.selection.type)&
                                                                                    data.for.vc.scatterplot.this.factor.level$selection.number.vector == this.selection.number) ),] 
    y.axis.label <- paste("log(Variance)",sep = "")
    this.min.y.axis <- log(min(data.for.vc.scatterplot.this.factor.level[,which(grepl("variance",
                                                                                      colnames(input.scatter.plot.data)))]))
    this.max.y.axis <- log(max(data.for.vc.scatterplot.this.factor.level[,which(grepl("variance",
                                                                                      colnames(input.scatter.plot.data)))]))
    #Source in some gplot code that will the first row plots; this time for the four trait values
    source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim/Functions_to_Make_Life_Easier/Make_VC_Scatter_Plots_Traits_20240715.R")
    fa.row.1 <- fa
    fb.row.1 <- fb
    fc.row.1 <- fc
    fd.row.1 <- fd
    
    #Source in some gplot code that will the first row plots; this time for the four genetic values
    source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim/Functions_to_Make_Life_Easier/Make_VC_Scatter_Plots_Genetic_Values_20240716.R")
    fa.row.2 <- fa
    fb.row.2 <- fb
    fc.row.2 <- fc
    fd.row.2 <- fd
    
    #Source in some gplot code that will the first row plots; this time for the four breeding values
    source("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim/Functions_to_Make_Life_Easier/Make_VC_Scatter_Plots_Breeding_Values_20240716.R")
    fa.row.3 <- fa
    fb.row.3 <- fb
    fc.row.3 <- fc
    fd.row.3 <- fd
    
    #Use the code below as a starting point for making the figures
    
    print(plot_grid(fa.row.1, fb.row.1, fc.row.1, fd.row.1,
                    fa.row.2, fb.row.2, fc.row.2, fd.row.2,
                    fa.row.3, fb.row.3, fc.row.3, fd.row.3,          
                    nrow = 3, ncol = 4))
    
    
  }#End  for(this.selection.number in unique(data.for.vc.scatterplot.this.factor.level$selection.number.vector)[-1])
}#End for(this.selection.type in unique(data.for.vc.scatterplot.this.factor.level$subpopulation.vector)[-1]) 
dev.off()


