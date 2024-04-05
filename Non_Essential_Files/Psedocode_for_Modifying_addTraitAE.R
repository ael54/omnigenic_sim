#THIS code is modified from https://github.com/gaynorr/AlphaSimR/blob/a77f3a274f785c53ce06716d460a50dc530ae0b1/R/Class-SimParam.R
#THIS is Pseudocode, and will not actually work.
#The purpose of this is to show potentially
# where the things we need for manually entering 
# Locations and positions of additive/epistatic QTLs

######AEL: Include the following as input parameters:
 # List positions of additive QTL
 # List of additive effects corresponding to each QTL
 # List/matrix of positions of epistatic QTL
 # List/matrix of positions of epistatic effects corresponding to each QTL
addTraitAEOmnigenic = function(nQtlPerChr,mean=0,var=1,relAA=0,corA=NULL,
                      corAA=NULL,useVarA=TRUE,gamma=FALSE,shape=1,force=FALSE,
                      name=NULL){
  if(!force){
    private$.isRunning()
  }
  if(length(nQtlPerChr)==1){
    nQtlPerChr = rep(nQtlPerChr,self$nChr)
  }
  nTraits = length(mean)
  relAA = relAA*4
  if(length(gamma)==1) gamma = rep(gamma,nTraits)
  if(length(shape)==1) shape = rep(shape,nTraits)
  if(length(relAA)==1) relAA = rep(relAA,nTraits)
  if(is.null(corA)) corA=diag(nTraits)
  if(is.null(corAA)) corAA=diag(nTraits)
  if(is.null(name)){
    name = paste0("Trait",1:nTraits+self$nTraits)
  }
  stopifnot(length(mean)==length(var),
            isSymmetric(corA),
            isSymmetric(corAA),
            length(relAA)==length(mean),
            length(mean)==nrow(corA),
            (sum(nQtlPerChr)%%2L)==0L,
            length(mean)==length(name))
 ###############AEL: Here is where specific loci can be hard coded in
   qtlLoci = private$.pickLoci(nQtlPerChr)#Replace with the list of additive QTL
  addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                      corr=corA,gamma=gamma,shape=shape)#Replace with list of additive QTL effects
  #New object: epiQtlLoci - this will be list/matrix of all a*a epistatic loci postions
  epiEff = sampEpiEff(qtlLoci=qtlLoci,nTraits=nTraits,
                      corr=corA,gamma=gamma,shape=shape,
                      relVar=relAA) #Replace with list of epistatic QTL effects
 
  ###############AEL: No more modifications were made to this function
  
   E = matrix(sample.int(sum(nQtlPerChr),sum(nQtlPerChr)),ncol=2)
  for(i in 1:nTraits){
    trait = new("TraitAE",
                qtlLoci,
                addEff=addEff[,i],
                epiEff=cbind(E,epiEff[,i]),
                intercept=0,
                name=name[i])
    tmp = calcGenParam(trait, self$founderPop,
                       self$nThreads)
    if(useVarA){
      scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
    }else{
      scale = sqrt(var[i])/sqrt(popVar(tmp$gv)[1])
    }
    trait@addEff = trait@addEff*scale
    trait@epiEff[,3] = trait@epiEff[,3]*scale
    trait@intercept = mean[i]-mean(tmp$gv*scale)
    if(useVarA){
      private$.addTrait(trait,var[i],popVar(tmp$gv*scale)[1])
    }else{
      private$.addTrait(trait,popVar(tmp$bv*scale)[1],var[i])
    }
  }
  invisible(self)
}