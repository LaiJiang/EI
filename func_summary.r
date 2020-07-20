

perpersonsummary <- function(fn.ind.data, 
                             functioncalls ) {
  # var.data global parameter
  Genes <- table(fn.ind.data$gene.name)
  #only keep these genes appears at least once
  Genes <- Genes[Genes>0]
  for (ii in (1:length(Genes) )) {

        one.data <- fn.ind.data[fn.ind.data$gene.name==names(Genes)[ii],]
    
    sone.data <- NULL
    for(func_ii in 1:length(functioncalls)){ 
      sone.data <- c(sone.data, functioncalls[[func_ii]]( fn.ind.data,one.data))  
      #print( c(func_ii,functioncalls[[func_ii]](one.data)) )
      #print(func_ii)
    }            
    
    if(length(sone.data)!=length(functioncalls))print( c( ii, length(sone.data)) )
    
    if (ii==1) {summar.data <- sone.data}
    if (ii>1) {summar.data <- rbind(summar.data, sone.data)}
    
  }
  
  #colnames(summar.data)[2] <- vlabel
  return(summar.data)
}


aggregate.snps  <- function(x1,x2, inter.way, aggre.way){
  
  if(inter.way=="M")vector.snp <- mapply(function(x,y)x*y, x1, x2)
  
  if(inter.way=="D")vector.snp <- mapply(function(x,y)x/y, x1, x2)
  
  aggre.way(vector.snp)
  
}

#####################
#for the DHS regsion thats nearest to this gene, is this gene also the nearest gene (among all neighbors) to this DHS
fn.gene.closest.DHS.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( unique(one.data$gene.name)==snp.id.overlap.dat$gene.name[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  else return(0)
}

fn.gene.closest.DHS <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( unique(one.data$gene.name)==snp.id.overlap.dat$gene.name[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  else return(0)
}


#for the DHS regsion thats nearest to this gene, the gwas beta for matched SNPs in the gene nearest to this DHS region
fn.gene.closest.DHS.beta.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$beta[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  else return(0)
}

fn.gene.closest.DHS.beta <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$beta[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  else return(0)
}


fn.gene.closest.DHS.lead.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$is.leadsnp[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  else return(0)
}


fn.gene.closest.DHS.lead <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$is.leadsnp[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  else return(0)
}




fn.gene.closest.DHS.tss.dist.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$snp.tss.dist[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  else return( median(fn.ind.data$snp.tss.dist) )
  
}


fn.gene.closest.DHS.tss.dist <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$snp.tss.dist[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  else return( median(fn.ind.data$snp.tss.dist) )
}


fn.gene.closest.DHS.tss.dist.inv.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( 1/ (0.01+ snp.id.overlap.dat$snp.tss.dist[which.min(  snp.id.overlap.dat$snp.gene.dist)]))
  else return(0)
}

fn.gene.closest.DHS.tss.dist.inv <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( 1/ (0.01+ snp.id.overlap.dat$snp.tss.dist[which.min(  snp.id.overlap.dat$snp.gene.dist)]))
  else return(0)
}

fn.gene.closest.DHS.tss.dist.inv2.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( 1/ (0.01+ snp.id.overlap.dat$snp.tss.dist[which.min(  snp.id.overlap.dat$snp.gene.dist)])^2)
  else return(0)
}

fn.gene.closest.DHS.tss.dist.inv2 <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( 1/ (0.01+ snp.id.overlap.dat$snp.tss.dist[which.min(  snp.id.overlap.dat$snp.gene.dist)])^2)
  else return(0)
}



fn.gene.closest.DHS.gene.dist.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return(  ( snp.id.overlap.dat$snp.gene.dist[which.min(  snp.id.overlap.dat$snp.gene.dist)]))
  else return(median(fn.ind.data$snp.gene.dist))
  
}


fn.gene.closest.DHS.gene.dist.inv.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( 1/ (0.01+ snp.id.overlap.dat$snp.gene.dist[which.min(  snp.id.overlap.dat$snp.gene.dist)]))
  else return(0)
  
}

fn.gene.closest.DHS.gene.dist.inv2.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( 1/ (0.01+ snp.id.overlap.dat$snp.gene.dist[which.min(  snp.id.overlap.dat$snp.gene.dist)])^2)
  else return(0)
  
}

fn.gene.closest.DHS.overlap.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$is.gene.overlap[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
  
  
}


fn.gene.closest.DHS.nearest.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$is.nearest.to.tss[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
  
}


fn.gene.closest.DHS.maf.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$maf[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
  
}


fn.gene.closest.DHS.se.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$se[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  else return(median(fn.ind.data$se))
  
}

fn.gene.closest.DHS.z.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$z[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  else return(0)
  
}

fn.gene.closest.DHS.prob.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$prob[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
}

fn.gene.closest.DHS.bf.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$log10bf[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
  
}

fn.gene.closest.DHS.dbsnp.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$is.dbsnp.delit[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
  
}



fn.gene.closest.DHS.snpeff.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$is.snpeff.delit[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return( 0 )
  
}


fn.gene.closest.DHS.counts.trait <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.trait.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$snp.in.trait.DHS[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
  
}



fn.gene.closest.DHS.gene.dist <- function(fn.ind.data, one.data)  {                                             
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return(  ( snp.id.overlap.dat$snp.gene.dist[which.min(  snp.id.overlap.dat$snp.gene.dist)]))
  else return(median(fn.ind.data$snp.gene.dist))
  
}


fn.gene.closest.DHS.gene.dist.inv <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( 1/ (0.01+ snp.id.overlap.dat$snp.gene.dist[which.min(  snp.id.overlap.dat$snp.gene.dist)]))
  else return(0)
  
}

fn.gene.closest.DHS.gene.dist.inv2 <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( 1/ (0.01+ snp.id.overlap.dat$snp.gene.dist[which.min(  snp.id.overlap.dat$snp.gene.dist)])^2)
  else return(0)
  
}

fn.gene.closest.DHS.overlap <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$is.gene.overlap[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
  
  
}


fn.gene.closest.DHS.nearest <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$is.nearest.to.tss[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
  
}


fn.gene.closest.DHS.maf <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$maf[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
  
}


fn.gene.closest.DHS.se <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$se[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  else return(median(fn.ind.data$se))
  
}

fn.gene.closest.DHS.z <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$z[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  else return(0)
  
}

fn.gene.closest.DHS.prob <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$prob[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
}

fn.gene.closest.DHS.bf <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$log10bf[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
  
}

fn.gene.closest.DHS.dbsnp <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$is.dbsnp.delit[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
  
}




fn.gene.closest.DHS.snpeff <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$is.snpeff.delit[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return( 0 )
  
}


fn.gene.closest.DHS.counts <- function(fn.ind.data, one.data)  {
  snp.id.overlap.dat  <- fn.ind.data[which( fn.ind.data$SNP %in% one.data$snp.name[which(one.data$nearest.DHS.from.gene==1)] ),]
  if(nrow(snp.id.overlap.dat)>0)return( snp.id.overlap.dat$snp.in.DHS[which.min(  snp.id.overlap.dat$snp.gene.dist)])
  
  else return(0)
  
}



#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################



fn.no.SNPs <- function(fn.ind.data, one.data)  {
  return(nrow(one.data))  }


fn.max.length <- function(fn.ind.data, one.data)  {
  return(max((one.data$gene.length)))  }

fn.max.length.inv <- function(fn.ind.data, one.data)  {
  return(max(1/( max(one.data$gene.length,0.01) )))  }


fn.max.length.inv2 <- function(fn.ind.data, one.data)  {
  return(max(1/( max(one.data$gene.length,0.01) )^2))  }


fn.strand.pos <- function(fn.ind.data, one.data)  {
  return(sum((one.data$gene.strand)=="+"))  }


fn.strand.neg <- function(fn.ind.data, one.data)  {
  return(sum((one.data$gene.strand)=="-"))  }


fn.max.BF <- function(fn.ind.data, one.data)  {
  return(max((one.data$log10bf)))  }


fn.mean.BF <- function(fn.ind.data, one.data)  {
  return(mean((one.data$log10bf)))  }



fn.max.maf <- function(fn.ind.data, one.data)  {
  return(max((one.data$maf)))  }


fn.mean.maf <- function(fn.ind.data, one.data)  {
  return(mean((one.data$maf)))  }



fn.max.beta <- function(fn.ind.data, one.data)  {
  return(max((one.data$beta)))  }


fn.mean.beta <- function(fn.ind.data, one.data)  {
  return(mean((one.data$beta)))  }


fn.max.se <- function(fn.ind.data, one.data)  {
  return(max((one.data$se)))  }


fn.mean.se <- function(fn.ind.data, one.data)  {
  return(mean((one.data$se)))  }



fn.max.prob <- function(fn.ind.data, one.data)  {
  return(max((one.data$prob)))  }


fn.mean.prob <- function(fn.ind.data, one.data)  {
  return(mean((one.data$prob)))  }



fn.max.entropy <- function(fn.ind.data, one.data)  {
  return(max(( -one.data$prob *(log(  max(one.data$prob,0.0001)  )) )))  }


fn.mean.entropy <- function(fn.ind.data, one.data)  {
  return(mean(( -one.data$prob *(log( max(one.data$prob,0.0001)  )) )))  }


#the location features are removed!!!
if(FALSE){
  fn.max.gene.tss <- function(fn.ind.data, one.data)  {
    return(max((one.data$gene.tss)))  }
  
  fn.mean.gene.tss <- function(fn.ind.data, one.data)  {
    return(mean((one.data$gene.tss)))  }
}




fn.max.snp.tss.dist <- function(fn.ind.data, one.data)  {
  return(max((one.data$snp.tss.dist)))  }

fn.min.snp.tss.dist <- function(fn.ind.data, one.data)  {
  return(min((one.data$snp.tss.dist)))  }


fn.max.snp.tss.dist.inv <- function(fn.ind.data, one.data)  {
  return(max(1/(  max( 0.01, one.data$snp.tss.dist)   )))  }

fn.max.snp.tss.dist.inv2 <- function(fn.ind.data, one.data)  {
  return(max(1/( max( 0.01, one.data$snp.tss.dist)  )^2))  }


fn.mean.snp.tss.dist <- function(fn.ind.data, one.data)  {
  return(mean((one.data$snp.tss.dist)))  }

fn.mean.snp.tss.dist.inv <- function(fn.ind.data, one.data)  {
  return(mean(1/( max(0.01,one.data$snp.tss.dist) )))  }

fn.mean.snp.tss.dist.inv2 <- function(fn.ind.data, one.data)  {
  return(mean(1/(max(0.01,one.data$snp.tss.dist))^2))  }





fn.max.snp.gene.dist <- function(fn.ind.data, one.data)  {
  return(max((one.data$snp.gene.dist)))  }

fn.max.snp.gene.dist.inv <-  function(fn.ind.data, one.data) {aggregate.snps(1, max(0.01,one.data$snp.gene.dist), "D", max)}
fn.max.snp.gene.dist.inv2 <-  function(fn.ind.data, one.data) {aggregate.snps(1, max(0.01,one.data$snp.gene.dist)^2, "D", max)}



fn.min.snp.gene.dist <- function(fn.ind.data, one.data)  {
  return(min((one.data$snp.gene.dist)))  }

fn.min.snp.gene.dist.inv <-  function(fn.ind.data, one.data) {aggregate.snps(1, max(0.01,one.data$snp.gene.dist), "D", min)}
fn.min.snp.gene.dist.inv2 <-  function(fn.ind.data, one.data) {aggregate.snps(1, max(0.01,one.data$snp.gene.dist)^2, "D", min)}



fn.mean.snp.gene.dist <- function(fn.ind.data, one.data)  {
  return(mean((one.data$snp.gene.dist)))  }

fn.mean.snp.gene.dist.inv <-  function(fn.ind.data, one.data) {aggregate.snps(1, max(0.01,one.data$snp.gene.dist), "D", mean)}
fn.mean.snp.gene.dist.inv2 <-  function(fn.ind.data, one.data) {aggregate.snps(1, max(0.01,one.data$snp.gene.dist)^2, "D", mean)}



fn.gene.overlap.snp <- function(fn.ind.data, one.data)  {
  return(sum((one.data$is.gene.overlap)))  }


fn.sum.is.nearest.to.tss <- function(fn.ind.data, one.data)  {
  return(sum((one.data$is.nearest.to.tss)))  }




fn.max.z <- function(fn.ind.data, one.data)  {
  return(max((one.data$z)))  }

fn.mean.z <- function(fn.ind.data, one.data)  {
  return(mean((one.data$z)))  }




fn.max.log10bf_group <- function(fn.ind.data, one.data)  {
  return(max((one.data$log10bf_group)))  }

fn.mean.log10bf_group <- function(fn.ind.data, one.data)  {
  return(mean((one.data$log10bf_group)))  }




fn.max.prob_group <- function(fn.ind.data, one.data)  {
  return(max((one.data$prob_group)))  }

fn.mean.prob_group <- function(fn.ind.data, one.data)  {
  return(mean((one.data$prob_group)))  }





fn.sum.is.dbsnp.delit <- function(fn.ind.data, one.data)  {
  return(sum((one.data$is.dbsnp.delit)))  }




fn.sum.is.snpeff.delit <- function(fn.ind.data, one.data)  {
  return(sum((one.data$is.snpeff.delit)))  }




fn.max.snp.in.trait.DHS <- function(fn.ind.data, one.data)  {
  return(max((one.data$snp.in.trait.DHS)))  }


fn.mean.snp.in.trait.DHS <- function(fn.ind.data, one.data)  {
  return(mean((one.data$snp.in.trait.DHS)))  }



fn.sum.nearest.trait.DHS.from.gene <- function(fn.ind.data, one.data)  {
  return(sum((one.data$nearest.trait.DHS.from.gene)))  }


########################################
#two way interactions between potential importan features


fn.max.bf.dist <-  function(fn.ind.data, one.data)  {
  return( max((one.data$log10bf /  max(one.data$snp.gene.dist , 0.01) )))  }


fn.max.bf.dist.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$log10bf, max(0.01,one.data$snp.gene.dist)^2, "D", max)}



fn.mean.bf.dist <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$log10bf /  max(one.data$snp.gene.dist,0.01) )))  }

fn.mean.bf.dist.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$log10bf, max(0.01,one.data$snp.gene.dist)^2, "D", mean)}


#beta X distance
fn.max.beta.dist <-  function(fn.ind.data, one.data)  {
  return( max((one.data$beta /  max(one.data$snp.gene.dist,0.01) )))  }

fn.max.beta.dist.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$beta, max(0.01,one.data$snp.gene.dist)^2, "D", max)}


fn.mean.beta.dist <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$beta /  max(one.data$snp.gene.dist,0.01) )))  }

fn.mean.beta.dist.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$beta, max(0.01,one.data$snp.gene.dist)^2, "D", mean)}

#BF X overlap 

fn.max.bf.overlap <-  function(fn.ind.data, one.data)  {
  return( max((one.data$log10bf *  (one.data$is.gene.overlap) )))  }


fn.mean.bf.overlap <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$log10bf *  (one.data$is.gene.overlap) )))  }

#beta X overlap

fn.max.beta.overlap <-  function(fn.ind.data, one.data)  {
  return( max((one.data$beta *  (one.data$is.gene.overlap) )))  }


fn.mean.beta.overlap <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$beta *  (one.data$is.gene.overlap) )))  }

##BF X chrm
fn.max.bf.snp.in.trait.DHS <-  function(fn.ind.data, one.data)  {
  return( max((one.data$log10bf *  (one.data$snp.in.trait.DHS) )))  }

fn.mean.bf.snp.in.trait.DHS <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$log10bf *  (one.data$snp.in.trait.DHS) )))  }

##Beta X chrm
fn.max.beta.snp.in.trait.DHS <-  function(fn.ind.data, one.data)  {
  return( max((one.data$beta *  (one.data$snp.in.trait.DHS) )))  }

fn.mean.beta.snp.in.trait.DHS <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$beta *  (one.data$snp.in.trait.DHS) )))  }

fn.max.bf.snp.in.DHS <-  function(fn.ind.data, one.data)  {
  return( max((one.data$log10bf *  (one.data$snp.in.DHS) )))  }

fn.mean.bf.snp.in.DHS <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$log10bf *  (one.data$snp.in.DHS) )))  }

##Beta X chrm
fn.max.beta.snp.in.DHS <-  function(fn.ind.data, one.data)  {
  return( max((one.data$beta *  (one.data$snp.in.DHS) )))  }

fn.mean.beta.snp.in.DHS <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$beta *  (one.data$snp.in.DHS) )))  }


#dist X overlap

fn.mean.overlap.dist <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$is.gene.overlap /  max(one.data$snp.gene.dist , 0.01) )))  }

fn.mean.overlap.dist.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, max(0.01 , one.data$snp.gene.dist)^2, "D", mean)}


fn.max.overlap.dist <-  function(fn.ind.data, one.data)  {
  return( max((one.data$is.gene.overlap /  max(one.data$snp.gene.dist , 0.01) )))  }

fn.max.overlap.dist.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, max(0.01 , one.data$snp.gene.dist)^2, "D", max)}

#dist X chrm

fn.mean.snp.in.trait.DHS.dist <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$snp.in.trait.DHS /  max(one.data$snp.gene.dist,0.01) )))  }

fn.mean.snp.in.trait.DHS.dist.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, max(0.01,one.data$snp.gene.dist)^2, "D", mean)}

fn.max.snp.in.trait.DHS.dist <-  function(fn.ind.data, one.data)  {
  return( max((one.data$snp.in.trait.DHS /  max(one.data$snp.gene.dist,0.01) )))  }

fn.max.snp.in.trait.DHS.dist.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, max(0.01,one.data$snp.gene.dist)^2, "D", max)}

fn.mean.snp.in.DHS.dist <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$snp.in.DHS /  max(one.data$snp.gene.dist,0.01) )))  }

fn.mean.snp.in.DHS.dist.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, max(0.01,one.data$snp.gene.dist)^2, "D", mean)}

fn.max.snp.in.DHS.dist <-  function(fn.ind.data, one.data)  {
  return( max((one.data$snp.in.DHS /  max(one.data$snp.gene.dist,0.01) )))  }

fn.max.snp.in.DHS.dist.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, max(0.01,one.data$snp.gene.dist)^2, "D", max)}



#overlap X chrm

fn.mean.overlap.snp.in.trait.DHS <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$snp.in.trait.DHS *  (one.data$is.gene.overlap) )))  }


fn.max.overlap.snp.in.trait.DHS <-  function(fn.ind.data, one.data)  {
  return( max((one.data$snp.in.trait.DHS *  (one.data$is.gene.overlap) )))  }

fn.mean.overlap.snp.in.DHS <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$snp.in.DHS *  (one.data$is.gene.overlap) )))  }


fn.max.overlap.snp.in.DHS <-  function(fn.ind.data, one.data)  {
  return( max((one.data$snp.in.DHS *  (one.data$is.gene.overlap) )))  }


########################################
#lead X distance

fn.max.lead.dist <-  function(fn.ind.data, one.data)  {
  return( max((one.data$is.leadsnp /  max(one.data$snp.gene.dist,0.01) )))  }

fn.max.lead.dist.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.leadsnp, max(0.01,one.data$snp.gene.dist)^2, "D", max)}


fn.mean.lead.dist <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$is.leadsnp /  max(one.data$snp.gene.dist,0.01) )))  }

fn.mean.lead.dist.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.leadsnp, max(0.01,one.data$snp.gene.dist)^2, "D", mean)}


#lead X overlap
fn.max.lead.overlap <-  function(fn.ind.data, one.data)  {
  return( max((one.data$is.leadsnp *  (one.data$is.gene.overlap) )))  }


#lead X overlap
fn.sum.lead.overlap <-  function(fn.ind.data, one.data)  {
  return( sum((one.data$is.leadsnp *  (one.data$is.gene.overlap) )))  }


#lead X nearest 
fn.max.lead.nearest <-  function(fn.ind.data, one.data)  {
  return( max((one.data$is.leadsnp *  (one.data$is.nearest.to.tss) )))  }

fn.sum.lead.nearest <-  function(fn.ind.data, one.data)  {
  return( sum((one.data$is.leadsnp *  (one.data$is.nearest.to.tss) )))  }

#lead X maf

fn.max.lead.maf <-  function(fn.ind.data, one.data)  {
  return( max((one.data$is.leadsnp *  (one.data$maf) )))  }

fn.sum.lead.maf <-  function(fn.ind.data, one.data)  {
  return( sum((one.data$is.leadsnp *  (one.data$maf) )))  }

#lead X beta

fn.max.lead.beta <-  function(fn.ind.data, one.data)  {
  return( max((one.data$is.leadsnp *  (one.data$beta) )))  }

fn.mean.lead.beta <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$is.leadsnp *  (one.data$beta) )))  }


#lead X Z

fn.max.lead.z <-  function(fn.ind.data, one.data)  {
  return( max((one.data$is.leadsnp *  (one.data$z) )))  }

fn.mean.lead.z <-  function(fn.ind.data, one.data)  {
  return( mean((one.data$is.leadsnp *  (one.data$z) )))  }


#lead x prob



fn.max.lead.prob <-  function(fn.ind.data, one.data)  {
  return( max((one.data$is.leadsnp *  (one.data$prob) )))  }


fn.sum.lead.prob <-  function(fn.ind.data, one.data)  {
  return( sum((one.data$is.leadsnp *  (one.data$prob) )))  }

#####################other interactions#####################
###############################################################
###############################################################
###############################################################
####################################################################################


######################################
#lead snps
fn.sum.lead.bf <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.leadsnp, one.data$log10bf, "M", sum)}
fn.max.lead.bf <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.leadsnp, one.data$log10bf, "M", max)}
fn.sum.lead.dbsnp <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.leadsnp, one.data$is.dbsnp.delit, "M", sum)}
fn.sum.lead.snpeff <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.leadsnp, one.data$is.snpeff.delit, "M", sum)}

fn.sum.lead.snp.in.trait.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.leadsnp, one.data$snp.in.trait.DHS, "M", sum)}
fn.sum.lead.near.snp.in.trait.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.leadsnp, one.data$nearest.trait.DHS.from.gene, "M", sum)}

fn.sum.lead.snp.in.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.leadsnp, one.data$snp.in.DHS, "M", sum)}
fn.sum.lead.near.snp.in.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.leadsnp, one.data$nearest.DHS.from.gene, "M", sum)}

######################################
#distance
fn.sum.dist.nearest.trait.DHS.from.gene <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$snp.tss.dist, "M", sum)}
fn.max.dist.nearest.trait.DHS.from.gene.inv <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, max(one.data$snp.tss.dist,0.01), "D", max)}
fn.max.dist.nearest.trait.DHS.from.gene.inv2 <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, max(one.data$snp.tss.dist,0.01)^2, "D", max)}
fn.sum.dist.nearest.DHS.from.gene <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$snp.tss.dist, "M", sum)}
fn.max.dist.nearest.DHS.from.gene.inv <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, max(one.data$snp.tss.dist,0.01), "D", max)}
fn.max.dist.nearest.DHS.from.gene.inv2 <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, max(one.data$snp.tss.dist,0.01)^2, "D", max)}

fn.mean.dist.prob <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$prob, max(one.data$snp.tss.dist,0.01), "D", mean)}
fn.mean.dist.prob.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$prob, max(one.data$snp.tss.dist,0.01)^2, "D", mean)}
fn.max.dist.prob <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$prob, max(one.data$snp.tss.dist,0.01), "D", max)}
fn.max.dist.prob.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$prob, max(one.data$snp.tss.dist,0.01)^2, "D", max)}
######################################
#snp.gene.distance: dist2
fn.sum.dist2.nearest.trait.DHS.from.gene <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$snp.gene.dist, "M", sum)}
fn.max.dist2.nearest.trait.DHS.from.gene.inv <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, max(one.data$snp.gene.dist,0.01), "D", max)}
fn.max.dist2.nearest.trait.DHS.from.gene.inv2 <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, max(one.data$snp.gene.dist,0.01)^2, "D", max)}
fn.sum.dist2.nearest.DHS.from.gene <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$snp.gene.dist, "M", sum)}
fn.max.dist2.nearest.DHS.from.gene.inv <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, max(one.data$snp.gene.dist,0.01), "D", max)}
fn.max.dist2.nearest.DHS.from.gene.inv2 <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, max(one.data$snp.gene.dist,0.01)^2, "D", max)}

fn.mean.dist2.prob <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$prob, max(one.data$snp.gene.dist,0.01), "D", mean)}
fn.mean.dist2.prob.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$prob, max(one.data$snp.gene.dist,0.01)^2, "D", mean)}

fn.max.dist2.prob <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$prob, max(one.data$snp.gene.dist,0.01), "D", max)}
fn.max.dist2.prob.square <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$prob, max(one.data$snp.gene.dist,0.01)^2, "D", max)}

######################################
#overlap 

fn.sum.overlap.bf <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, one.data$log10bf, "M", sum)}
fn.max.overlap.bf <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, one.data$log10bf, "M", max)}
fn.sum.overlap.dbsnp <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, one.data$is.dbsnp.delit, "M", sum)}
fn.sum.overlap.snpeff <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, one.data$is.snpeff.delit, "M", sum)}

fn.sum.overlap.snp.in.trait.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, one.data$snp.in.trait.DHS, "M", sum)}
fn.sum.overlap.near.nearest.trait.DHS.from.gene <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, one.data$nearest.trait.DHS.from.gene, "M", sum)}

fn.sum.overlap.snp.in.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, one.data$snp.in.DHS, "M", sum)}
fn.sum.overlap.near.nearest.DHS.from.gene <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, one.data$nearest.DHS.from.gene, "M", sum)}


######################################
fn.sum.nearest.nearest.trait.DHS.from.gene.bf <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$log10bf, "M", sum)}
fn.max.nearest.nearest.trait.DHS.from.gene.bf <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$log10bf, "M", max)}
fn.sum.nearest.nearest.trait.DHS.from.gene.dbsnp <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$is.dbsnp.delit, "M", sum)}
fn.sum.nearest.nearest.trait.DHS.from.gene.snpeff <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$is.snpeff.delit, "M", sum)}
fn.sum.nearest.nearest.trait.DHS.from.gene.snp.in.trait.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$snp.in.trait.DHS, "M", sum)}
fn.sum.nearest.nearest.trait.DHS.from.gene.p <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$prob, "M", sum)}
fn.max.nearest.nearest.trait.DHS.from.gene.p <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$prob, "M", max)}
fn.sum.nearest.nearest.trait.DHS.from.gene.beta <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$beta, "M", sum)}
fn.max.nearest.nearest.trait.DHS.from.gene.beta <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$beta, "M", max)}

fn.sum.nearest.nearest.DHS.from.gene.bf <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$log10bf, "M", sum)}
fn.max.nearest.nearest.DHS.from.gene.bf <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$log10bf, "M", max)}
fn.sum.nearest.nearest.DHS.from.gene.dbsnp <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$is.dbsnp.delit, "M", sum)}
fn.sum.nearest.nearest.DHS.from.gene.snpeff <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$is.snpeff.delit, "M", sum)}
fn.sum.nearest.nearest.DHS.from.gene.snp.in.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$snp.in.DHS, "M", sum)}
fn.sum.nearest.nearest.DHS.from.gene.p <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$prob, "M", sum)}
fn.max.nearest.nearest.DHS.from.gene.p <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$prob, "M", max)}
fn.sum.nearest.nearest.DHS.from.gene.beta <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$beta, "M", sum)}
fn.max.nearest.nearest.DHS.from.gene.beta <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$beta, "M", max)}

######################################
#the chromatin counts interactions
fn.sum.snp.in.trait.DHS.count.bf <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, one.data$log10bf, "M", sum)}
fn.max.snp.in.trait.DHS.count.bf <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, one.data$log10bf, "M", max)}
fn.sum.snp.in.trait.DHS.count.dbsnp <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, one.data$is.dbsnp.delit, "M", sum)}
fn.sum.snp.in.trait.DHS.count.snpeff <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, one.data$is.snpeff.delit, "M", sum)}
fn.sum.snp.in.trait.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, 1, "M", sum)}
fn.sum.snp.in.trait.DHS.count.p <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, one.data$prob, "M", sum)}
fn.max.snp.in.trait.DHS.count.p <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, one.data$prob, "M", max)}
fn.sum.snp.in.trait.DHS.count.beta <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, one.data$beta, "M", sum)}
fn.max.snp.in.trait.DHS.count.beta <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, one.data$beta, "M", max)}

fn.nearest.trait.DHS.from.gene.dist <- function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.tss.dist, one.data$nearest.trait.DHS.from.gene, "M", max)}
fn.nearest.trait.DHS.from.gene.dist.inv <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, max(one.data$snp.tss.dist,0.01), "D", max)}
fn.nearest.trait.DHS.from.gene.dist.inv2 <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, max(one.data$snp.tss.dist,0.01)^2, "D", max)}


fn.sum.snp.in.DHS.count.bf <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, one.data$log10bf, "M", sum)}
fn.max.snp.in.DHS.count.bf <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, one.data$log10bf, "M", max)}
fn.sum.snp.in.DHS.count.dbsnp <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, one.data$is.dbsnp.delit, "M", sum)}
fn.sum.snp.in.DHS.count.snpeff <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, one.data$is.snpeff.delit, "M", sum)}
fn.sum.snp.in.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, 1, "M", sum)}
fn.sum.snp.in.DHS.count.p <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, one.data$prob, "M", sum)}
fn.max.snp.in.DHS.count.p <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, one.data$prob, "M", max)}
fn.sum.snp.in.DHS.count.beta <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, one.data$beta, "M", sum)}
fn.max.snp.in.DHS.count.beta <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, one.data$beta, "M", max)}

fn.nearest.DHS.from.gene.dist <- function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.tss.dist, one.data$nearest.DHS.from.gene, "M", max)}
fn.nearest.DHS.from.gene.dist.inv <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, max(one.data$snp.tss.dist,0.01), "D", max)}
fn.nearest.DHS.from.gene.dist.inv2 <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, max(one.data$snp.tss.dist,0.01)^2, "D", max)}


########################################################################
########################################################################
######################################################################## Vince new features

#if this gene is the nreaest gene to the trait =matche dopen chromatin SNP.
fn.max.snp.in.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, 1, "M", max)}
fn.sum.snp.in.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, 1, "M", sum)}
fn.max.snp.in.trait.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, 1, "M", max)}
fn.sum.snp.in.trait.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, 1, "M", sum)}

fn.max.nearest.DHS.from.gene <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, 1, "M", max)}
fn.sum.nearest.DHS.from.gene <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, 1, "M", sum)}
fn.max.nearest.gene.from.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.gene.from.DHS, 1, "M", max)}
fn.sum.nearest.gene.from.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.gene.from.DHS, 1, "M", sum)}

fn.max.nearest.trait.DHS.from.gene <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, 1, "M", max)}
fn.sum.nearest.trait.DHS.from.gene <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, 1, "M", sum)}
fn.max.nearest.gene.from.trait.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.gene.from.trait.DHS, 1, "M", max)}
fn.sum.nearest.gene.from.trait.DHS <-  function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.gene.from.trait.DHS, 1, "M", sum)}
################################################################################################

#number of genes in the same locus
fn.locus.no.genes <-  function(fn.ind.data, one.data) {
  length(unique(fn.ind.data$gene.name[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)]))
}

#number of SNPs in the same locus
fn.locus.no.SNPs <-  function(fn.ind.data, one.data) {
  length(unique(fn.ind.data$snp.name[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)]))
}

#number of lead SNPs =1 in any locus


#max of GWWAS beta in the same locus
fn.max.locus.beta <-  function(fn.ind.data, one.data) {
  max(fn.ind.data$beta[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)])
}

##########locus specific maf
fn.max.locus.maf <-  function(fn.ind.data, one.data) {
  max(fn.ind.data$maf[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)])
}

##########locus specific se
fn.max.locus.se <-  function(fn.ind.data, one.data) {
  max(fn.ind.data$se[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)])
}


##########locus specific z
fn.max.locus.z <-  function(fn.ind.data, one.data) {
  max(fn.ind.data$z[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)])
}

##########locus specific prob
fn.max.locus.prob <-  function(fn.ind.data, one.data) {
  max(fn.ind.data$prob[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)])
}


##########locus specific log10bf
fn.max.locus.bf <-  function(fn.ind.data, one.data) {
  max(fn.ind.data$log10bf[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)])
}


#mean of GWWAS beta in the same locus
fn.mean.locus.beta <-  function(fn.ind.data, one.data) {
  mean(unique(fn.ind.data$beta[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)]))
}

##########locus specific maf
fn.mean.locus.maf <-  function(fn.ind.data, one.data) {
  mean(unique(fn.ind.data$maf[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)]))
}

##########locus specific se
fn.mean.locus.se <-  function(fn.ind.data, one.data) {
  mean(unique(fn.ind.data$se[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)]))
}


##########locus specific z
fn.mean.locus.z <-  function(fn.ind.data, one.data) {
  mean(unique(fn.ind.data$z[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)]))
}

##########locus specific prob
fn.mean.locus.prob <-  function(fn.ind.data, one.data) {
  mean(unique(fn.ind.data$prob[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)]))
}


##########locus specific log10bf
fn.mean.locus.bf <-  function(fn.ind.data, one.data) {
  mean(unique(fn.ind.data$log10bf[which(fn.ind.data$locus.name %in% one.data$locus.name ==1)]))
}



###########################eqtl score

fn.sum.gtex <- function(fn.ind.data, one.data)  {
  return(sum((one.data$in.gtex)))  }


fn.mean.gtex <- function(fn.ind.data, one.data)  {
  return(mean((one.data$in.gtex)))  }

fn.max.gtex <- function(fn.ind.data, one.data)  {
  return(max((one.data$in.gtex)))  }




###########################

############################
#if the gene overlap a DHS region
fn.gene.ovelap.DHS <- function(fn.ind.data, one.data)  {
  
  #if gene overlap SNP
  if(sum(one.data$is.gene.overlap)>0){
    
    snp.overlap <- one.data$snp.name[which(one.data$is.gene.overlap==1)]
    return ( max(fn.ind.data$nearest.DHS.from.gene[match(snp.overlap,fn.ind.data$snp.name)]) )
  }
  
  else { return(0) }
  
}


############################
############################
############################
############################
#snp_effect variables
fn.max.1m.snpeff.impact.none <- function(fn.ind.data, one.data) {max(1-one.data$snpeff.impact.none)}
fn.sum.1m.snpeff.impact.none <- function(fn.ind.data, one.data) {sum(1-one.data$snpeff.impact.none)}

fn.sum.snpeff.impact.modifier <- function(fn.ind.data, one.data) {sum(one.data$snpeff.impact.modifier)}
fn.max.snpeff.impact.modifier <- function(fn.ind.data, one.data) {max(one.data$snpeff.impact.modifier)}

fn.sum.snpeff.impact.low <- function(fn.ind.data, one.data) {sum(one.data$snpeff.impact.low)}
fn.max.snpeff.impact.low <- function(fn.ind.data, one.data) {max(one.data$snpeff.impact.low)}

fn.sum.snpeff.impact.moderate <- function(fn.ind.data, one.data) {sum(one.data$snpeff.impact.moderate)}
fn.max.snpeff.impact.moderate <- function(fn.ind.data, one.data) {max(one.data$snpeff.impact.moderate)}

fn.sum.snpeff.impact.high <- function(fn.ind.data, one.data) {sum(one.data$snpeff.impact.high)}
fn.max.snpeff.impact.high <- function(fn.ind.data, one.data) {max(one.data$snpeff.impact.high)}


##interaction between overlapping SNPs and impacts
fn.sum.overlap.1m.snpeff.one <- function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, 1-one.data$snpeff.impact.none, "M", sum)}
fn.max.overlap.1m.snpeff.one <- function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, 1-one.data$snpeff.impact.none, "M", max)}

##interaction between overlapping SNPs and modifiers
fn.sum.overlap.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, one.data$snpeff.impact.modifier, "M", sum)}
fn.max.overlap.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$is.gene.overlap, one.data$snpeff.impact.modifier, "M", max)}



###################interatction between modifers and DHS


fn.sum.snp.in.trait.DHS.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, one.data$snpeff.impact.modifier, "M", sum)}
fn.max.snp.in.trait.DHS.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.trait.DHS, one.data$snpeff.impact.modifier, "M", max)}

fn.sum.nearest.trait.DHS.from.gene.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$snpeff.impact.modifier, "M", sum)}
fn.max.nearest.trait.DHS.from.gene.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene, one.data$snpeff.impact.modifier, "M", max)}


fn.sum.nearest.gene.from.trait.DHS.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.gene.from.trait.DHS, one.data$snpeff.impact.modifier, "M", sum)}
fn.max.nearest.gene.from.trait.DHS.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.gene.from.trait.DHS, one.data$snpeff.impact.modifier, "M", max)}


fn.sum.snp.in.DHS.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, one.data$snpeff.impact.modifier, "M", sum)}
fn.max.snp.in.DHS.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$snp.in.DHS, one.data$snpeff.impact.modifier, "M", max)}


fn.sum.nearest.DHS.from.gene.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$snpeff.impact.modifier, "M", sum)}
fn.max.nearest.DHS.from.gene.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene, one.data$snpeff.impact.modifier, "M", max)}

fn.sum.nearest.gene.from.DHS.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.gene.from.DHS, one.data$snpeff.impact.modifier, "M", sum)}
fn.max.nearest.gene.from.DHS.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.gene.from.DHS, one.data$snpeff.impact.modifier, "M", max)}

fn.sum.nearest.DHS.from.gene.dist.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene/(one.data$snp.gene.dist+0.01), one.data$snpeff.impact.modifier, "M", sum)}
fn.max.nearest.DHS.from.gene.dist.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.DHS.from.gene/(one.data$snp.gene.dist+0.01), one.data$snpeff.impact.modifier, "M", max)}

fn.sum.nearest.trait.DHS.from.gene.dist.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene/(one.data$snp.gene.dist+0.01), one.data$snpeff.impact.modifier, "M", sum)}
fn.max.nearest.trait.DHS.from.gene.dist.modifier <- function(fn.ind.data, one.data) {aggregate.snps(one.data$nearest.trait.DHS.from.gene/(one.data$snp.gene.dist+0.01), one.data$snpeff.impact.modifier, "M", max)}





