

factor.A <- NULL
factor.B <- NULL
factor.C <- NULL
factor.D <- NULL
for(i in 1:3){
  for(j in 1:4){
    for(k in 1:4){
       for(el in 1:4){
         factor.A <- c(factor.A, rep(i,1))
         factor.B <- c(factor.B, rep(j,1))
         factor.C <- c(factor.C, rep(k,1))
         factor.D <- c(factor.D, rep(el,1))
      }#End for(el in 1:4)
    }#End for(k in 1:4)
  }#End for(j in 1:4)
}#End for(i in 1:3)

these.factorial.levels <- data.frame(factor.A, factor.B, factor.C, factor.D)

setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/Housekeeping_notes_for_factorial_Experiment")

write.csv(these.factorial.levels, "Tally_of_Experimental_Levels_Ran.csv", row.names = FALSE)
