#Make a box plot for Factor A  
title.label = NA
fa <- ggplot(input.scatter.plot.data, aes(x=generation.vector)) + 
             geom_point(col ="blue", aes(y=log(genetic.value.1.variance))) + ylim(this.min.y.axis,this.max.y.axis) +
  if(!is.na(title.label)) {
    labs(title = title.label, x=paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, ", No. Inds. Selected = ", this.selection.number, sep = ""), 
         y = y.axis.label, size = 25)
  }else{
    labs(x= paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, ", No. Inds. Selected = ", this.selection.number, sep = ""), y = y.axis.label, size = 25)
  }

title.label = NA
fb <- ggplot(input.scatter.plot.data, aes(x=generation.vector)) + 
  geom_point(col ="blue", aes(y=log(genetic.value.2.variance))) +  ylim(this.min.y.axis,this.max.y.axis) +
  if(!is.na(title.label)) {
    labs(title = title.label, x=paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, ", No. Inds. Selected = ", this.selection.number, sep = ""), 
         y = y.axis.label, size = 25)
  }else{
    labs(x= paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, ", No. Inds. Selected = ", this.selection.number, sep = ""), y = y.axis.label, size = 25)
  }



title.label = NA
fc <- ggplot(input.scatter.plot.data, aes(x=generation.vector)) + 
  geom_point(col ="blue", aes(y=log(genetic.value.3.variance))) +  ylim(this.min.y.axis,this.max.y.axis) +
  if(!is.na(title.label)) {
    labs(title = title.label, x=paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, ", No. Inds. Selected = ", this.selection.number, sep = ""), 
         y = y.axis.label, size = 25)
  }else{
    labs(x= paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, ", No. Inds. Selected = ", this.selection.number, sep = ""), y = y.axis.label, size = 25)
  }


title.label = NA
fd <- ggplot(input.scatter.plot.data, aes(x=generation.vector)) + 
  geom_point(col ="blue", aes(y=log(genetic.value.4.variance))) +  ylim(this.min.y.axis,this.max.y.axis) +
  if(!is.na(title.label)) {
    labs(title = title.label, x=paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, ", No. Inds. Selected = ", this.selection.number, sep = ""), 
         y = y.axis.label, size = 25)
  }else{
    labs(x= paste("Factor = ", substr(colnames(data.for.vc.scatterplot)[this.factor],start = 7,stop = 9),", Level = ",  this.level, ", Selection = ", this.selection.type, ", No. Inds. Selected = ", this.selection.number, sep = ""), y = y.axis.label, size = 25)
  }



 
  