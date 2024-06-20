#Delete the code below
val.pop <- 1 #The index of the for loop through validation populations


###Format the data for boxplot so that a plot can be made out of them
data.for.boxplot.this.validation.pop <- data.for.boxplot[which(data.for.boxplot$validation.set.vector == val.pop),]
data.for.boxplot.this.validation.pop <- na.omit(data.for.boxplot.this.validation.pop)


data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop <- melt(data = data.for.boxplot.this.validation.pop,
                                                                           measure.vars = c("predictive.ability.GBLUP",
                                                                                            "predictive.ability.Multi.Kern.add",
                                                                                            "predictive.ability.Multi.Kern.epi"))
variable.truncated <- substr(data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop$variable,start = 20, stop = 3000)
data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop <- data.frame(data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop, variable.truncated)


#Create a separate data set for each training population combination
data.for.boxplot.this.validation.pop.melted.first.pop.train <- data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop[which(
  (data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop$training.set.1.boolean == 1)&
    (data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop$training.set.2.boolean == 0)),]

data.for.boxplot.this.validation.pop.melted.second.pop.train<- data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop[which(
  (data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop$training.set.1.boolean == 0)&
    (data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop$training.set.2.boolean == 1)),]

data.for.boxplot.this.validation.pop.melted.both.pop.train<- data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop[which(
  (data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop$training.set.1.boolean == 1)&
    (data.for.boxplot.this.validation.pop.melted.sep.for.each.train.pop$training.set.2.boolean == 1)),]
