#Make a box plot for Factor A  
title.label = "Trait 1 -whole shebang"
fa <- ggplot(input.scatter.plot.data, aes(x=generation.vector)) + 
             geom_point(col ="red", aes(y=trait.1.variance)) + ylim(this.min.y.axis,this.max.y.axis) +
  if(!is.na(title.label)) {
    labs(title = title.label, x=paste("Factor 1 = ", i, " Selection = ", this.selection.type, sep = ""), 
         y = y.axis.label, size = 25)
  }else{
    labs(x= paste("Factor ", i, sep = ""), y = y.axis.label, size = 25)
  }

title.label = "Trait 2 - core component"
fb <- ggplot(input.scatter.plot.data, aes(x=generation.vector)) + 
  geom_point(col ="red", aes(y=trait.2.variance)) +  ylim(this.min.y.axis,this.max.y.axis) +
  if(!is.na(title.label)) {
    labs(title = title.label, x=paste("Factor 1 = ", i, " Selection = ", this.selection.type, sep = ""), 
         y = y.axis.label, size = 25)
  }else{
    labs(x= paste("Factor ", i, sep = ""), y = y.axis.label, size = 25)
  }



title.label = "Trait 3 - peripheral component"
fc <- ggplot(input.scatter.plot.data, aes(x=generation.vector)) + 
  geom_point(col ="red", aes(y=trait.3.variance)) +  ylim(this.min.y.axis,this.max.y.axis) +
  if(!is.na(title.label)) {
    labs(title = title.label, x=paste("Factor 1 = ", i, " Selection = ", this.selection.type, sep = ""), 
         y = y.axis.label, size = 25)
  }else{
    labs(x= paste("Factor ", i, sep = ""), y = y.axis.label, size = 25)
  }


title.label = "Trait 4 - between core and perph"
fd <- ggplot(input.scatter.plot.data, aes(x=generation.vector)) + 
  geom_point(col ="red", aes(y=trait.4.variance)) +  ylim(this.min.y.axis,this.max.y.axis) +
  if(!is.na(title.label)) {
    labs(title = title.label, x=paste("Factor 1 = ", i, " Selection = ", this.selection.type, sep = ""), 
         y = y.axis.label, size = 25)
  }else{
    labs(x= paste("Factor ", i, sep = ""), y = y.axis.label, size = 25)
  }

#Temporary code
pdf("VC.experimental.plot.pdf", width = 50)
print(plot_grid(fa,fb,fc,fd,
                nrow = 1, ncol = 4))
dev.off()
#End temporary code

####Old code, which can be updated

#Make a box plot for Factor B  
fb <- ggplot(input.box.plot.data, aes(x=as.factor(factor.B.vector), 
                                y=value, fill=variable.truncated)) + 
  geom_boxplot() + ylim(this.min.y.axis,this.max.y.axis) + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25),
                         plot.title = element_text(size = 30), legend.position = "none")+
if(!is.na(title.label)) {
  labs(title =  title.label, x="Cor Epi vs Cor Add", y = y.axis.label, size = 25)
}else{
  labs(x="Cor Epi vs Cor Add", y = y.axis.label, size = 25)
}

#Make a box plot for Factor C  
fc <- ggplot(input.box.plot.data, aes(x=as.factor(factor.C.vector), 
                                y=value, fill=variable.truncated)) + 
  geom_boxplot() + ylim(this.min.y.axis,this.max.y.axis) +  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25),
                         plot.title = element_text(size = 30), legend.position = "none")+
  
if(!is.na(title.label)) {
  labs(title =  title.label, x="Per Epi vs Per Add", y = y.axis.label, size = 25)
}else{
  labs(x="Per Epi vs Per Add", y = y.axis.label, size = 25)
}
#Make a box plot for Factor D 
fd <- ggplot(input.box.plot.data, aes(x=as.factor(factor.D.vector), 
                                y=value, fill=variable.truncated)) + 
  geom_boxplot() + ylim(this.min.y.axis,this.max.y.axis) + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25),
                         plot.title = element_text(size = 30), legend.text = element_text(size = 20))+
if(!is.na(title.label)) {
  labs(title =  title.label, x="Btw Epi vs Per Add", y = y.axis.label, size = 25)
}else{
  labs(x="Btw Epi vs Per Add", y = y.axis.label, size = 25)
}


 
  