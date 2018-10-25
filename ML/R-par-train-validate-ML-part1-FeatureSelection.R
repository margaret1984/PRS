#################This script is to write ALL steps of ML in one script, staring from data split, filtering and finally model fit. File for all with SNPs filtered by phred and/or plink files for GWAS are neeeded 
#########http://appliedpredictivemodeling.com/blog/2017/9/2/njdc83d01pzysvvlgik02t5qnaljnd ###nested CV in a lool 
library(doParallel)
library(parallel)
library(caret)
library(sda)
library(glmnet)
library(data.table)
###############
####Read data and set parametres####
###############
dat <- fread("snp-data.txt",	header= TRUE ,	sep= "\t")
###############
####get training/testing files: single run or loop ####
###############
set.seed(100) ##make sure results are replicable
trainIndex <- createDataPartition(dat$Response_class, p = .70, list = FALSE, times=1) ### do it once, to get training and holdout set, based on predefined class, e.g. responde to treatment yes/no
dat.train <- as.data.frame(dat[ trainIndex,])
dat.holdout <- as.data.frame(dat[ -trainIndex,])
set.seed(100)
rep.folds<-createMultiFolds(dat.train$Response_class, k = 5, times = 10) ###for repeated cross-validation
#folds <- createFolds(dat.train$Response_class, k = 5,returnTrain = TRUE) ####for cross validation
cl <- makeCluster(5) ###let's go parallel
registerDoParallel(cl)
outputfile.int.cv ="cross-validation-performance.txt" ##see how methods are working
write.table(t("Selection_method	Fold	Accuracy	Sensitivity	Specificity	Alg	AccuracyPValue	AccuracyLower	AccuracyUpper	McnemarPValue"), fileoutputfile.int.cv,quote=F, row.names=F,col.names=F,sep="\t",append=T)
###############
###FEATURE SELECTION###
###############
result <- foreach(i = 1:50,.packages  =c("glmnet","sda","caret")) %dopar% {
    comp <- tryCatch({
       rep.cv.training <- as.data.frame(dat.train[rep.folds[[i]],])
       rep.cv.testing <- as.data.frame(dat.train[-rep.folds[[i]],])
       snps.train<- as.matrix (rep.cv.training  [56:length(rep.cv.training)]) 
        cv.fitControl <- trainControl(method = "cv",  number = 5, classProbs = TRUE,   summaryFunction = twoClassSummary, savePredictions = TRUE)
        cv.ml.methods=c("rpart","C5.0","gbm","naive_bayes","rf","glm","svmRadial","svmLinear")
####################
###CAT###
#################
        set.seed(100)
        ranking.DDA <- sda.ranking(Xtrain = snps.train, L = rep.cv.training$Response_class, fdr = TRUE)
        best.idx <- ranking.DDA[which(ranking.DDA[,"lfdr"] < 0.8)]
        snps.sda <- snps.train[,best.idx]
        sda.train.set =  rep.cv.training[c("Response_class",colnames(snps.sda))]
        sda.test.set =  rep.cv.testing[c("Response_class",colnames(snps.sda))]
        output.file.sda <- paste("post-CAT0.8-n=",i,"-Response_class.txt",sep="") ###save SNPs found in each fold
        write.table(names(best.sda), output.file.sda, row.names=F,  quote=F, sep="\t")
        for (alg in 1:length(cv.ml.methods)){
            set.seed(100)
            model.train <- train(Response_class ~., data = sda.train.set, method = cv.ml.methods[alg],tuneLength = 10,  trControl = cv.fitControl,metric = "ROC")
            model.pred <- predict(model.train, newdata = sda.test.set)
            mc <- confusionMatrix(data =  model.pred, as.factor(sda.test.set[,1]))
            measures.cv.test <- cbind ("CAT",i,mc$overall["Accuracy"],mc$byClass["Sensitivity"], mc$byClass["Specificity"] , cv.ml.methods[alg],mc$overall["AccuracyPValue"],mc$overall["AccuracyLower"],mc$overall["AccuracyUpper"],mc$overall["McnemarPValue"] )
            write.table(measures.cv.test , file=outputfile.int.cv,quote=F, row.names=F,col.names=F,sep="\t",append=T) 
        }
        


#################
####LASSO###
################
        set.seed(100)
        fit<-glmnet(snps.train, rep.cv.training$Response_class, family="binomial", alpha=1,  standardize=TRUE) ####alpha=1 = LASSO, =0 Ridge, =.5 Elastic net
        coefs=coef(fit,s=min(fit$lambda))
        lasso.coeff.best=cbind(coefs@Dimnames[[1]][which(coefs != 0 ) ], coefs[ which(coefs != 0 ) ]) 
        output.file.lasso=paste(ml.path,"post-glmnet-LASSO-n=",i,"-Response_class-",parameter,".txt",sep="")
        write.table(lasso.coeff.best, output.file.lasso, row.names=F,  quote=F, sep="\t")
        lasso.train.set =  rep.cv.training[c("Response_class",lasso.coeff.best[2:length(lasso.coeff.best[,1])])]
        lasso.test.set =  rep.cv.testing[c("Response_class",lasso.coeff.best[2:length(lasso.coeff.best[,1])])]

        for (alg in 1:length(cv.ml.methods)){
            set.seed(100)
            model.train <- train(Response_class ~., data = lasso.train.set, method = cv.ml.methods[alg],tuneLength = 10,  trControl = cv.fitControl,metric = "ROC")
            model.pred <- predict(model.train, newdata = lasso.test.set)
            mc <- confusionMatrix(data =  model.pred, as.factor(lasso.test.set[,1]))
            measures.cv.test <- cbind ("LASSO",i,mc$overall["Accuracy"],mc$byClass["Sensitivity"], mc$byClass["Specificity"] , cv.ml.methods[alg],mc$overall["AccuracyPValue"],mc$overall["AccuracyLower"],mc$overall["AccuracyUpper"],mc$overall["McnemarPValue"] )
            write.table(measures.cv.test , file=outputfile.int.cv,quote=F, row.names=F,col.names=F,sep="\t",append=T) 
        }
              

        
    }, error=function(e) NULL)
    
}
stopCluster (cl)
####################
####get SNPs into one file and select those present in certain no of folds###
################
###use cat to set all SNPs into one
cat1 <- paste("cat post-CAT0.8*.*phred=2.txt >> rep10cv5-CAT0.8-snps-ld0.8-phred=2.txt") ###use cat to set all SNPs into one
system(cat1)
cat2 <- paste("cat  repcv5-post-glmnet-LASSO-*.*phred=2.txt >> rep10cv5-glmnet-LASSO-snps-ld0.8-phred=2.txt") 
system(cat2)
snps.lasso <- read.table("rep10cv5-glmnet-LASSO-snps-ld0.8-phred=2.txt",header=F)
sel.snps.lasso <- which(table(snps.lasso[1])>24) ####LASSO, select SNPS, present in certain number of initial folds
snps.cat <- read.table("rep10cv5-CAT0.8-snps-ld0.8-phred=2.txt",header=F)
sel.snps.cat <- table(snps.cat)[table(snps.cat)>24] ####CAT
# ###break line ###
fil.par <- "rep10cv5-CAT0.8-snps-ld0.8-phred=2-present=25"
write.table(sel.snps.cat,file=paste(ml.path, fil.par,"-pre.txt",sep=""),sep="\n",quote = F, append=FALSE)
grepp <- paste("grep rs ", ml.path,fil.par, "-pre.txt >>",ml.path, fil.par,".txt",sep="")
system(grepp)
