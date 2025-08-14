#Make a box plot for Factor A  
title.label = "Trait 1 -whole shebang"
fa <- ggplot(input.scatter.plot.data, aes(x=generation.vector)) + 
             geom_point(col ="red", aes(y=log(trait.1.variance))) + ylim(this.min.y.axis,this.max.y.axis) +
  if(!is.na(title.label)) {
    labs(title = title.label, x=paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, sep = ""), 
         y = y.axis.label, size = 25)
  }else{
    labs(x= paste("Factor = ", this.factor,", Level = ",  this.level, ", Selection = ", this.selection.type, sep = ""), y = y.axis.label, size = 25)
  }


title.label = "Trait 2 - core component"
fb <- ggplot(input.scatter.plot.data, aes(x=generation.vector)) + 
  geom_point(col ="red", aes(y=log(trait.2.variance))) +  ylim(this.min.y.axis,this.max.y.axis) +
  if(!is.na(title.label)) {
    labs(title = title.label, x=paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, sep = ""), 
         y = y.axis.label, size = 25)
  }else{
    labs(x= paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, sep = ""), y = y.axis.label, size = 25)
  }



title.label = "Trait 3 - peripheral component"
fc <- ggplot(input.scatter.plot.data, aes(x=generation.vector)) + 
  geom_point(col ="red", aes(y=log(trait.3.variance))) +  ylim(this.min.y.axis,this.max.y.axis) +
  if(!is.na(title.label)) {
    labs(title = title.label, x=paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, sep = ""), 
         y = y.axis.label, size = 25)
  }else{
    labs(x= paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, sep = ""), y = y.axis.label, size = 25)
  }


title.label = "Trait 4 - between core and perph"
fd <- ggplot(input.scatter.plot.data, aes(x=generation.vector)) + 
  geom_point(col ="red", aes(y=log(trait.4.variance))) +  ylim(this.min.y.axis,this.max.y.axis) +
  if(!is.na(title.label)) {
    labs(title = title.label, x=paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, sep = ""), 
         y = y.axis.label, size = 25)
  }else{
    labs(x= paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, sep = ""), y = y.axis.label, size = 25)
  }


  