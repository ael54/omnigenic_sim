#Move these booleans into the main code once you are done pressure testing this



the.plot.trait.breeding.value <- the.plot.trait.breeding.value + ylab("Breeding Value")
#Initiate the ggplot objects
if(plot.trait.1){
  tthe.plot.trait.breeding.value<- the.plot.trait.breeding.value  + 
    geom_smooth(data = input.scatter.plot.data, mapping = aes(x = generation.vector, y= breeding.value.1), size = 0.5, linetype = this.line.type, method = "loess", col = "black")
}#end if(plot.trait.1)
if(plot.trait.2){
  the.plot.trait.breeding.value <- the.plot.trait.breeding.value+
    geom_smooth(data = input.scatter.plot.data, mapping = aes(x = generation.vector, y= breeding.value.2), size = 0.5, linetype = this.line.type, method = "loess", col = "red")
}#end if(plot.trait.2)
if(plot.trait.3){
  the.plot.trait.breeding.value <- the.plot.trait.breeding.value +
    geom_smooth(data = input.scatter.plot.data, mapping = aes(x = generation.vector, y= breeding.value.3), size = 0.5, linetype = this.line.type,method = "loess", col = "blue") 
}#end if(plot.trait.3)
if(plot.trait.4){
  the.plot.trait.breeding.value <- the.plot.trait.breeding.value+ 
    geom_smooth(data = input.scatter.plot.data, mapping = aes(x = generation.vector, y= breeding.value.4), size = 0.5, linetype = this.line.type, method = "loess", col = "purple")
}#end if(plot.trait.4)

#The code below will put the points back in. Add this to each line within the if(plot.trait.i) statements
#geom_point(data= input.scatter.plot.data, aes(x = generation.vector, y= breeding.value.1), alpha = 0.10, col = "black")+
#geom_point(data= input.scatter.plot.data, aes(x = generation.vector, y= breeding.value.2), alpha = 0.10, col = "red")+
#geom_point(data= input.scatter.plot.data, aes(x = generation.vector,y= breeding.value.3), alpha = 0.10,  col = "blue" )
#geom_point(data= input.scatter.plot.data, aes(x = generation.vector, y= breeding.value.4), alpha = 0.10, col = "purple")+
  