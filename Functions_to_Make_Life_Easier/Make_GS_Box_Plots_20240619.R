

#####This code is in the process of being written

val.pop <- 1 #The index of the for loop through validation populations
#Next idea: have separate graphs per each training set
data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop <- melt(data = data.for.boxplot.this.validation.pop,
                                                                           measure.vars = c("predictive.ability.GBLUP",
                                                                                            "predictive.ability.Multi.Kern.add",
                                                                                            "predictive.ability.Multi.Kern.epi"))
variable.truncated <- substr(data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop$variable,start = 20, stop = 3000)
data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop <- data.frame(data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop, variable.truncated)


  #Make a box plot for Factor A  
  
  #Make a box plot for Factor B  
  
  #Make a box plot for Factor C  
  
  #Make a box plot for Factor D  