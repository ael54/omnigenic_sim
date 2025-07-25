###############################################################################   
############################################################################### 
###############################################################################   
# The vast majority of this script was written by Chris Gaynor. Alex Lipka
# made some minor, trivial tweaks so that this can be ran with the rest
# of the simulation pipeline.
##### Written approximately May 11-15, 2024

get.me.my.breeding.values <- function(SP.within.function = NULL,
                                      current.generation = NULL,
                                      traits = NULL,
                                      core.genes = NULL,
                                      peripheral.genes = NULL,
                                      core.core.epistasis = NULL,
                                      peri.peri.epistasis = NULL,
                                      core.peri.epistasis = NULL,
                                      calculate.epi.dev = FALSE){
  # Extract QTL genotypes
  M = pullQtlGeno(pop=current.generation, trait = 1, asRaw = FALSE, 
                  simParam = SP.within.function)
  
  # Center genotypes for calculation of breeding value, this sucks off the mean
  M_centered = scale(M, scale=FALSE)
  
  #Loop through each of the traits.  
  breeding.values <- matrix(0,nrow = nrow(traits), ncol = 4)
  rownames(breeding.values) = traits$`<Trait>`
  colnames(breeding.values) = colnames(traits)[2:5]
  for(trait.i in 1:4){
    print(paste("----------Calculating breeding value for trait ", trait.i,"---------", sep = ""))
    #Pick out the trait of interest. 
    this.trait <- traits[, (trait.i+1)]
    if(trait.i == 1){
      #Pull out the all additive and epistatic effects
      colnames(core.genes) <- c("QTL", "Add.effect")
      colnames(peripheral.genes) <- c("QTL", "Add.effect")
      all.additive.effects <- rbind(core.genes, peripheral.genes) 
      
      #Put the core epistatic effects into three objects
      loc.1 <- NULL
      loc.2 <- NULL
      epi.effect <- NULL
      for(trait.j in 1:nrow(core.core.epistasis)){
        if(!is.na(core.core.epistasis$epi.effects[trait.j])){
          loc.1 <- c(loc.1, core.core.epistasis$core.core.epi[trait.j])
          epi.effect <- c(epi.effect, core.core.epistasis$epi.effects[trait.j])
        }else{
          loc.2 <- c(loc.2, core.core.epistasis$core.core.epi[trait.j])
        }
      }#end for(trait.j in 1:nrow(core.core.epistasis))
      
      #Add the peripheral epistatic effects into these three objects
      
      for(trait.j in 1:nrow(peri.peri.epistasis)){
        if(!is.na(peri.peri.epistasis$epi.effects[trait.j])){
          loc.1 <- c(loc.1, peri.peri.epistasis$peri.peri.epi[trait.j])
          epi.effect <- c(epi.effect, peri.peri.epistasis$epi.effects[trait.j])
        }else{
          loc.2 <- c(loc.2, peri.peri.epistasis$peri.peri.epi[trait.j])
        }
      }#end for(trait.j in 1:nrow(peri.peri.epistasis))
      
      #Add the between core-peripheral epistatic effects into these three objects
      for(trait.j in 1:nrow(core.peri.epistasis)){
        if(!is.na(core.peri.epistasis$epi.effects[trait.j])){
          loc.1 <- c(loc.1, core.peri.epistasis$core.peri.epi[trait.j])
          epi.effect <- c(epi.effect, core.peri.epistasis$epi.effects[trait.j])
        }else{
          loc.2 <- c(loc.2, core.peri.epistasis$core.peri.epi[trait.j])
        }
      }#end for(trait.j in 1:nrow(core.peri.epistasis))
      
      
      
      all.epi.effects <- data.frame(loc.1,loc.2, epi.effect)
      
    }# end if(trait.i == 1)
    if(trait.i == 2){
      #Pull out  all core additive and epistatic effects
      
      colnames(core.genes) <- c("QTL", "Add.effect")
      all.additive.effects <- core.genes
      
      loc.1 <- NULL
      loc.2 <- NULL
      epi.effect <- NULL
      for(trait.j in 1:nrow(core.core.epistasis)){
        if(!is.na(core.core.epistasis$epi.effects[trait.j])){
          loc.1 <- c(loc.1, core.core.epistasis$core.core.epi[trait.j])
          epi.effect <- c(epi.effect, core.core.epistasis$epi.effects[trait.j])
        }else{
          loc.2 <- c(loc.2, core.core.epistasis$core.core.epi[trait.j])
        }
      }#end for(trait.j in 1:nrow(core.core.epistasis))
    
      
      all.epi.effects <- data.frame(loc.1,loc.2, epi.effect)
      
      
    }#end  if(trait.i == 2)
    if(trait.i == 3){
      #Pull out all peripheral additive and epistatic effects
 
      colnames(peripheral.genes) <- c("QTL", "Add.effect")
      all.additive.effects <- peripheral.genes 
      
      loc.1 <- NULL
      loc.2 <- NULL
      epi.effect <- NULL
    
      #Add the peripheral epistatic effects into these three objects
      
      for(trait.j in 1:nrow(peri.peri.epistasis)){
        if(!is.na(peri.peri.epistasis$epi.effects[trait.j])){
          loc.1 <- c(loc.1, peri.peri.epistasis$peri.peri.epi[trait.j])
          epi.effect <- c(epi.effect, peri.peri.epistasis$epi.effects[trait.j])
        }else{
          loc.2 <- c(loc.2, peri.peri.epistasis$peri.peri.epi[trait.j])
        }
      }#end for(trait.j in 1:nrow(peri.peri.epistasis))
      
      all.epi.effects <- data.frame(loc.1,loc.2, epi.effect)
      
  
      
    } #end if(trait.i == 3)
    if(trait.i == 4){
      #Pull out only the epistatic effects between the core and peripheral
      all.additive.effects <- NULL
      
      loc.1 <- NULL
      loc.2 <- NULL
      epi.effect <- NULL
  
      #Add the between core-peripheral epistatic effects into these three objects
      for(trait.j in 1:nrow(core.peri.epistasis)){
        if(!is.na(core.peri.epistasis$epi.effects[trait.j])){
          loc.1 <- c(loc.1, core.peri.epistasis$core.peri.epi[trait.j])
          epi.effect <- c(epi.effect, core.peri.epistasis$epi.effects[trait.j])
        }else{
          loc.2 <- c(loc.2, core.peri.epistasis$core.peri.epi[trait.j])
        }
      }#end for(trait.j in 1:nrow(core.peri.epistasis))
      
      
      
      all.epi.effects <- data.frame(loc.1,loc.2, epi.effect)
      
    } #end if(trait.i == 4)
  
  
  ############ Calculate breeding values
  # Get the column numbers of in M_Centered that correspond to each QTL
  if(!is.null(all.additive.effects)){
    index.of.additive.effects <- NULL
    for(index.i in 1:nrow(all.additive.effects)){
      index.of.additive.effects <- c(index.of.additive.effects,
            which(all.additive.effects$QTL[index.i]  == colnames(M_centered)))
    }#end for(index.i in 1:nrow(all.additive.effects))
    
    BV = M_centered[,index.of.additive.effects] %*% all.additive.effects$Add.effect
    mu = mean((M[,index.of.additive.effects]-1)%*%all.additive.effects$Add.effect) + mean(this.trait)

  }else{
    BV = 0
    mu = mean(this.trait)
  } #end else 
  
  # Calculate frequencies for all genotype classes
  # This is needed for calculating the average effect of an allele substitution
  freq = matrix(0, 
                nrow=3, # Three genotype classes for diploids (0,1,2)
                ncol=ncol(M))
  
  for(i.epi in seq_len(3)){
    freq[i.epi,] = colSums(M == (i.epi-1)) / nrow(M)
  }
  
  colnames(freq) = colnames(M)
  
  # Create a variable to hold the additive-by-additive epistatic deviations
  # Making it a matrix to match BV above...ALEX SAYS::::I am changing BV 
  #to trait so that this will still work for trait 4, which has no additive
  # effects

  AA = matrix(0, nrow=nrow(traits))
  
  # Loop over all pairs of interacting loci
  #
  # I'll calculate the average effect for each locus to find breeding values
  # I'll then calculate the additive-by-additive epistatic deviations using the 
  # lack-of-fit after accounting for the breeding values 
  # This approach works because there is no other form of epistasis in the simulation.
  for(i.epi in seq_len(nrow(all.epi.effects))){
      if(i.epi %% 1000 == 0){ 
       print(paste("Calculating epistatic contribution to breeding value for ", i.epi,"th locus pair", sep = ""))
      }#end if(i.epi %% 1000 == 0)
      # Extract positions for loci and the epistatic effect
      loc1 = which(colnames(freq) == all.epi.effects[i.epi,1]) # Locus 1
      loc2 = which(colnames(freq) == all.epi.effects[i.epi,2]) # Locus 2
      eff = all.epi.effects[i.epi,3] # Effect
      
      
      
      # Construct centered genotypes for both loci
      x = c(-1,0,1)
      x1 = x - sum(freq[,loc1]*x)
      x2 = x - sum(freq[,loc2]*x)
      
      
  
      # Calculate marginal means for locus 1 and center
      # Values correspond to genotypes 0, 1, and 2
      y1 = c((freq[1,loc2]-freq[3,loc2])*eff,
             0,
             (freq[3,loc2]-freq[1,loc2])*eff)
      y1 = y1 - sum(y1*freq[,loc1])
      
      # Calculate marginal means for locus 2 and center
      # Values correspond to genotypes 0, 1, and 2
      y2 = c((freq[1,loc1]-freq[3,loc1])*eff,
             0,
             (freq[3,loc1]-freq[1,loc1])*eff)
      y2 = y2 - sum(y2*freq[,loc2])
      
      
      # Calculate average effects for both loci
      # Average effects are just regression coefficients
      alpha = c(
        sum(freq[,loc1]*y1*x1) / sum(freq[,loc1]*x1^2),
        sum(freq[,loc2]*y2*x2) / sum(freq[,loc2]*x2^2)
      )
      
      # Account for a potential divide by zero case
      if(any(is.na(alpha))) alpha[is.na(alpha)] = 0
      
      # Calculate breeding values
      tmpBv = M_centered[, c(loc1, loc2)] %*% alpha
      
      # Add breeding values to the global breeding values
      # At this point the values could be kept separate to partition variance due 
      # to different types of QTL
      BV = BV + tmpBv
      
      if(calculate.epi.dev){
        # Calculate epistatic deviations using lack-of-fit
        # This calculation could be skipped if you don't care about epistatic variance
        # It takes considerably more time than the previous calculations
        
        # Construct a 3x3 table for the epistatic deviations
        # The rows will represent the 0, 1, and 2 genotypes for locus 1
        # The columns will represent the 0, 1, and 2 genotypes for locus 2
        # This will be done in three parts
        # 1. Construct initial matrix using genetic values
        # 2. Subtract average genetic value using observed genotype frequencies (includes LD)
        # 3. Subtract breeding values
        devTab = tcrossprod(x) * eff
        freqTab = matrix(0, nrow=3, ncol=3)
        for(j.epi in seq_len(nrow(M))){
          freqTab[M[j.epi ,loc1]+1, M[j.epi ,loc2]+1] = freqTab[M[j.epi ,loc1]+1, M[j.epi ,loc2]+1] + 1
        }
        freqTab = freqTab / nrow(M)
        tmpMu = sum(freqTab*devTab)
        devTab = devTab - tmpMu
        devTab = sweep(devTab, 1, x1*alpha[1], "-")
        devTab = sweep(devTab, 2, x2*alpha[2], "-")
        
        # Determine which element of the table to access for each individual
        # This is being done to avoid using a double loop
        # To make sense of this trick, note that R saves data in a column major order
        # Thus, a value of 3 would be column 1, row 3 and a value of 4 would be column 2, row 1
        ind = M[,loc1] + 3*M[,loc2] + 1
        
        # Get each individual's epistatic deviation
        tmpAa = devTab[ind]
        
        # Add deviations to global value
        AA = AA + tmpAa
        
        # Add contribution to trait mean
        mu = mu + tmpMu
      } # End if(calculate.epi.dev)
    }#End for(i.epi in seq_len(nrow(all.epi.effects)))
    #Add the breeding value of the ith trait to the ith column of BV
    breeding.values[,trait.i] = BV
  }#End for(trait.i in 1:4)

  
  
  ##########This return code is here as a placeholder. It will get updated
  ###once this part of the code is written
  return(breeding.values = breeding.values)

} #end get.me.my.breeding.values()
