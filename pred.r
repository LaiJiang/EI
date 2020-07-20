#for each (one of 12) training trait, predict test trait


#############
#load test trait data to be predicted
path.main <- "/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/lai.jiang/Brent/Github/" 

path.save <- paste0(path.main,"save/")


path.dat <- paste0(path.main,"dat/")


#load test traits

load(paste0(path.dat, "covid19.nejm.2.RData"),verbose=TRUE)

#############
path.trained <- "/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/lai.jiang/Brent/snpeff/save_rm/locus/xg/"

trained.models <- list.files(path.trained)
  
#the predicted probabilities by trained traits
all.probs <- NULL
for(ii.model in 1:length(trained.model)){


load(paste0(path.trained, trained.models[ii.model]))

#############################################################
#first re-order test.x to be the same format of train.x (by selected features)

test.x <- snp.collapse[,match(colnames(train.x),colnames(snp.collapse))]

#the default 
test.y <- rep(1:nrow(test.x)) 

#standardization
test.x <- sweep(test.x, 2, mn.x, FUN = "-")
test.x <- sweep(test.x, 2, sd.x, FUN = "/")

library(xgboost)
library(mlr)
library(pROC)


testset<-data.frame(cbind(test.x,test.y))
testset$test.y<-as.factor(testset$test.y)
colnames(testset)<-c(colnames(test.x),"Response")
testtask <- makeClassifTask (data = testset,target = "Response")
xgpred <- predict(xgmodel,testtask)
pred.prob <- xgpred[["data"]][["prob.1"]]

all.probs <- cbind(all.probs,pred.prob)
}


colnames(all.probs) <-c( "calcium" ,   "dbilirubin", "dbp",       
                          "ebmd",       "glucose" ,   "height"   , 
                          "ldl" ,       "lowtsh"   ,  "rbc"    ,   
                          "sbp" ,       "t2d"   ,     "tg"   )
  
rownames(all.probs) <- rownames(test.x)

save(all.probs, file=paste0(path.save,"Covid19_pred.RData"))

