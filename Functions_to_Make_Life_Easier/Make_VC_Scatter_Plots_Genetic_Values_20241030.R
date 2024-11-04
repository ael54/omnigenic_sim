#Move these booleans into the main code once you are done pressure testing this



the.plot.genetic.value.vc <- the.plot.genetic.value.vc + ylab("log(Genetic Value Variance)")
#Initiate the ggplot objects
if(plot.trait.1){
  the.plot.genetic.value.vc <- the.plot.genetic.value.vc  + geom_point(data= input.scatter.plot.data, aes(x = generation.vector, y= log(genetic.value.1.variance)), alpha = 0.01, col = "black")+
    geom_smooth(data = input.scatter.plot.data, mapping = aes(x = generation.vector, y= log(genetic.value.1.variance)), size = 0.5, linetype = this.line.type, method = "loess", col = "black")
}#end if(plot.trait.1)
if(plot.trait.2){
  the.plot.genetic.value.vc <- the.plot.genetic.value.vc+geom_point(data= input.scatter.plot.data, aes(x = generation.vector, y= log(genetic.value.2.variance)), alpha = 0.01, col = "red")+
    geom_smooth(data = input.scatter.plot.data, mapping = aes(x = generation.vector, y= log(genetic.value.2.variance)), size = 0.5, linetype = this.line.type, method = "loess", col = "red")
}#end if(plot.trait.2)
if(plot.trait.3){
  the.plot.genetic.value.vc <- the.plot.genetic.value.vc + geom_point(data= input.scatter.plot.data, aes(x = generation.vector,y= log(genetic.value.3.variance)), alpha = 0.01,  col = "blue" )+
    geom_smooth(data = input.scatter.plot.data, mapping = aes(x = generation.vector, y= log(genetic.value.3.variance)), size = 0.5, linetype = this.line.type,method = "loess", col = "blue") 
}#end if(plot.trait.3)
if(plot.trait.4){
  the.plot.genetic.value.vc <- the.plot.genetic.value.vc + geom_point(data= input.scatter.plot.data, aes(x = generation.vector, y= log(genetic.value.4.variance)), alpha = 0.01, col = "purple")+
    geom_smooth(data = input.scatter.plot.data, mapping = aes(x = generation.vector, y= log(genetic.value.4.variance)), size = 0.5, linetype = this.line.type, method = "loess", col = "purple")
}#end if(plot.trait.4)




  