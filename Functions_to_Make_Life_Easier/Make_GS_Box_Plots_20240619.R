#Make a box plot for Factor A  
fa <- ggplot(input.box.plot.data, aes(x=as.factor(factor.A.vector), 
                                y=value, fill=variable.truncated)) + 
  geom_boxplot() + ylim(this.min.y.axis,this.max.y.axis) + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25),
                         plot.title = element_text(size = 30), legend.position = "none")+
  if(!is.na(title.label)) {
    labs(title = title.label, x="Cor Add vs Per Add", y = y.axis.label, size = 25)
  }else{
    labs(x="Cor Add vs Per Add", y = y.axis.label, size = 25)
  }

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


 
  