#Move these booleans into the main code once you are done pressure testing this



the.plot.trait.vc <- the.plot.trait.vc + ylab("log(Trait Variance)")
#Initiate the ggplot objects
if(plot.trait.1){
  the.plot.trait.vc <- the.plot.trait.vc  + geom_point(data= input.scatter.plot.data, aes(x = generation.vector, y= log(trait.1.variance)), alpha = 0.10, col = "black")+
    geom_smooth(data = input.scatter.plot.data, mapping = aes(x = generation.vector, y= log(trait.1.variance)), size = 0.5, linetype = this.line.type, method = "loess", col = "black")
}#end if(plot.trait.1)
if(plot.trait.2){
  the.plot.trait.vc <- the.plot.trait.vc+geom_point(data= input.scatter.plot.data, aes(x = generation.vector, y= log(trait.2.variance)), alpha = 0.10, col = "red")+
    geom_smooth(data = input.scatter.plot.data, mapping = aes(x = generation.vector, y= log(trait.2.variance)), size = 0.5, linetype = this.line.type, method = "loess", col = "red")
}#end if(plot.trait.2)
if(plot.trait.3){
  the.plot.trait.vc <- the.plot.trait.vc + geom_point(data= input.scatter.plot.data, aes(x = generation.vector,y= log(trait.3.variance)), alpha = 0.10,  col = "blue" )+
    geom_smooth(data = input.scatter.plot.data, mapping = aes(x = generation.vector, y= log(trait.3.variance)), size = 0.5, linetype = this.line.type,method = "loess", col = "blue") 
}#end if(plot.trait.3)
if(plot.trait.4){
  the.plot.trait.vc <- the.plot.trait.vc + geom_point(data= input.scatter.plot.data, aes(x = generation.vector, y= log(trait.4.variance)), alpha = 0.10, col = "purple")+
    geom_smooth(data = input.scatter.plot.data, mapping = aes(x = generation.vector, y= log(trait.4.variance)), size = 0.5, linetype = this.line.type, method = "loess", col = "purple")
}#end if(plot.trait.4)




  