#################This script is to write ALL steps of ML in one script, staring from data split, filtering and finally model fit. File for all with SNPs filtered by phred and/or plink files for GWAS are neeeded 
#########http://appliedpredictivemodeling.com/blog/2017/9/2/njdc83d01pzysvvlgik02t5qnaljnd ###nested CV in a lool 
library(doParallel)
library(parallel)
library(caret)
library(sda)
library(glmnet)
library(data.table)
library(pROC)
library(ggplot2)
library(plotROC)
source("crossValidationSet.R") ######Anne-Christin script 
###############
####Read data and set parametres####
###############
dat <- fread("./IRL_training/irl-random-lables-ld0.8-maf0.05-phred=2.txt",	header= TRUE ,	sep= "\t") ####use for CADD-filtered SNPS only
irl.path <- "/Users/malgorzata_maciukiewicz/Documents/camh/IRL-Grey/analyses/GWAS/"
stard.path <- "/Users/malgorzata_maciukiewicz/Documents/camh/STAR*D/GWAS/"
ml.path <- "./IRL_training/response/"
parameter <- "ld0.8-maf0.05-phred=2"
###############
####get training/testing files: single run or loop ####
###############
set.seed(100)
trainIndex <- createDataPartition(dat$Response_50HAMD, p = .70, list = FALSE, times=1) ### do it once, to get training and holdout set; repeat PCA and all 100 times
dat.train <- as.data.frame(dat[ trainIndex,])
dat.holdout <- as.data.frame(dat[ -trainIndex,])
#write.table(file="irl-trainig-response-clin.txt",dat.train[1:55],sep="\t",quote=F,row.names=F)
set.seed(100)
rep.folds<-createMultiFolds(dat.train$Response_50HAMD, k = 5, times = 10)
#folds <- createFolds(dat.train$Response_50HAMD, k = 5,returnTrain = TRUE)
cl <- makeCluster(5) 
registerDoParallel(cl)
outputfile.int.cv = paste(ml.path,"random-labels-internal-rep10cv5-CAT0.8-trained-irl",parameter,"-only-genet.txt",sep="")
write.table(t("Fold	Accuracy	Sensitivity	Specificity	Alg	AccuracyPValue	AccuracyLower	AccuracyUpper	McnemarPValue"), file=outputfile.int.cv,quote=F, row.names=F,col.names=F,sep="\t",append=T)
###############
###FEATURE SELECTION ONLY###
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
        ranking.DDA <- sda.ranking(Xtrain = snps.train, L = rep.cv.training$Response_50HAMD, fdr = TRUE)
        best.idx <- ranking.DDA[which(ranking.DDA[,"lfdr"] < 0.8)]
        best.sda <- which(ranking.DDA[,"lfdr"] < 0.8)
        snps.sda <- snps.train[,best.idx]
        sda.train.set =  rep.cv.training[c("Response_50HAMD",colnames(snps.sda))]
        sda.test.set =  rep.cv.testing[c("Response_50HAMD",colnames(snps.sda))]
#         output.file.sda <- paste(ml.path,"post-CAT0.8-n=",i,"-Response_50HAMD-",parameter,".txt",sep="")
#         write.table(names(best.sda), output.file.sda, row.names=F,  quote=F, sep="\t")
        for (alg in 1:length(cv.ml.methods)){
            set.seed(100)
            model.train <- train(Response_50HAMD ~., data = sda.train.set, method = cv.ml.methods[alg],tuneLength = 10,  trControl = cv.fitControl,metric = "ROC")
            model.pred <- predict(model.train, newdata = sda.test.set)
            mc <- confusionMatrix(data =  model.pred, as.factor(sda.test.set[,1]))
            measures.cv.test <- cbind (i,mc$overall["Accuracy"],mc$byClass["Sensitivity"], mc$byClass["Specificity"] , cv.ml.methods[alg],mc$overall["AccuracyPValue"],mc$overall["AccuracyLower"],mc$overall["AccuracyUpper"],mc$overall["McnemarPValue"] )
            write.table(measures.cv.test , file=outputfile.int.cv,quote=F, row.names=F,col.names=F,sep="\t",append=T) 
        }
        
        
        

# ################
####LASSO###
################
#         set.seed(100)
#         fit<-glmnet(snps.train, rep.cv.training$Response_50HAMD, family="binomial", alpha=1,  standardize=TRUE) ####alpha=1 = LASSO, =0 Ridge, =.5 Elastic net
#         coefs=coef(fit,s=min(fit$lambda))
#         lasso.coeff.best=cbind(coefs@Dimnames[[1]][which(coefs != 0 ) ], coefs[ which(coefs != 0 ) ]) 
#         output.file.lasso=paste(ml.path,"post-glmnet-LASSO-n=",i,"-Response_50HAMD-",parameter,".txt",sep="")
#         write.table(lasso.coeff.best, output.file.lasso, row.names=F,  quote=F, sep="\t")
#         lasso.train.set =  rep.cv.training[c("Response_50HAMD",lasso.coeff.best[2:length(lasso.coeff.best[,1])])]
#         lasso.test.set =  rep.cv.testing[c("Response_50HAMD",lasso.coeff.best[2:length(lasso.coeff.best[,1])])]
# 
#         for (alg in 1:length(cv.ml.methods)){
#             set.seed(100)
#             model.train <- train(Response_50HAMD ~., data = lasso.train.set, method = cv.ml.methods[alg],tuneLength = 10,  trControl = cv.fitControl,metric = "ROC")
#             model.pred <- predict(model.train, newdata = lasso.test.set)
#             mc <- confusionMatrix(data =  model.pred, as.factor(lasso.test.set[,1]))
#             measures.cv.test <- cbind (i,mc$overall["Accuracy"],mc$byClass["Sensitivity"], mc$byClass["Specificity"] , cv.ml.methods[alg],mc$overall["AccuracyPValue"],mc$overall["AccuracyLower"],mc$overall["AccuracyUpper"],mc$overall["McnemarPValue"] )
#             write.table(measures.cv.test , file=outputfile.int.cv,quote=F, row.names=F,col.names=F,sep="\t",append=T) 
#         }
#        
# 
#         
    }, error=function(e) NULL)
    
}
stopCluster (cl)
####################
####get SNPs into one file and select those present in certain no of folds###
################
# cat1 <- paste("cat ", ml.path,"post-CAT0.8*.*phred=2.txt >>",ml.path, "rep10cv5-CAT0.8-snps-ld0.8-phred=2.txt",sep="") ###use cat to set all SNPs into one
# system(cat1)
# cat2 <- paste("cat ", ml.path,"repcv5-post-glmnet-LASSO-*.*phred=2.txt >>",ml.path, "rep10cv5-glmnet-LASSO-snps-ld0.8-phred=2.txt",sep="") ###use cat to set all SNPs into one
# system(cat2)
# snps.lasso <- read.table(paste(ml.path, "rep10cv5-glmnet-LASSO-snps-ld0.8-phred=2.txt",sep=""),header=F)
# sel.snps.lasso <- which(table(snps.lasso[1])>24) ####LASSO
# snps.cat <- read.table(paste(ml.path, "rep10cv5-CAT0.8-snps-ld0.8-phred=2.txt",sep=""),header=F)
# sel.snps.cat <- table(snps.cat)[table(snps.cat)>24] ####CAT
# ###break line ###
# fil.par <- "rep10cv5-CAT0.8-snps-ld0.8-phred=2-present=25"
# write.table(sel.snps.cat,file=paste(ml.path, fil.par,"-pre.txt",sep=""),sep="\n",quote = F, append=FALSE)
# grepp <- paste("grep rs ", ml.path,fil.par, "-pre.txt >>",ml.path, fil.par,".txt",sep="")
# system(grepp)
####LASSO###
# train.newdat=dat.train[c("Response_50HAMD","HAMD1_BL","HAMD2_BL","HAMD13_BL", "SEX",names(sel.snps.lasso[3:length(sel.snps.lasso)-1]))]
# houldout.newdat=dat.holdout[c("Response_50HAMD","HAMD1_BL","HAMD2_BL","HAMD13_BL", "SEX",names(sel.snps.lasso[3:length(sel.snps.lasso)-1]))]
####CAT###
# train.newdat=dat.train[c("Response_50HAMD","HAMD1_BL","HAMD2_BL","HAMD13_BL", "SEX",names(sel.snps.cat[1:length(sel.snps.cat)-1]))]
# houldout.newdat=dat.holdout[c("Response_50HAMD","HAMD1_BL","HAMD2_BL","HAMD13_BL", "SEX",names(sel.snps.cat[1:length(sel.snps.cat)-1]))]
###############
###get STAR*D SNPs####
###############
# extract.sel <- paste("plink1.9 --bfile ", stard.path,"stard-eur-irl-ch --extract " ,ml.path, fil.par,".txt --recode A  --tab --out ",  stard.path, fil.par, sep="")
# system(extract.sel)
# system("perl -pi -w -e 's/_[ACGT12]\t/\t/g;' /Users/malgorzata_maciukiewicz/Documents/camh/STAR*D/GWAS/*.*raw")
# system("perl -pi -w -e 's/_[ACGT12]\n/\n/g;' /Users/malgorzata_maciukiewicz/Documents/camh/STAR*D/GWAS/*.*raw") 
# system("perl -pi -w -e 's/NA/1/g;' /Users/malgorzata_maciukiewicz/Documents/camh/STAR*D/GWAS/*.*raw") 
# stard.file <- paste(stard.path,fil.par,".raw", sep="")
# stard.snps <- read.table(stard.file,header=T)
# stard.pheno.file <- paste(stard.path,"pheno-hamd-ml-mds.txt", sep="")
# stard.pheno <- read.table(stard.pheno.file,header=T,sep="\t")
# dat.val <- merge(stard.pheno, stard.snps, by=c("FID","IID"))
# old.idx <- which(dat.val$Age_at_entry_to_study > 60)
# dat.val.old <- dat.val[old.idx,]
# colnames(dat.val.old)[which(names(dat.val.old) == "SEX.x")] <- "SEX"
################
###prepare files for ML###
###############
#######clinical variables for irl-responsee "HAMD1_BL","HAMD2_BL","HAMD13_BL", "SEX.x",
#######clinical variables for irl-remissionn"HAMD1_BL","HAMD2_BL","HAMD3_BL","HAMD7_BL", "SEX.x",
# val.newdat <- dat.val.old[c( "Response_50HAMD","HAMD1_BL","HAMD2_BL","HAMD13_BL", "SEX",names(dat.val.old[64:length(dat.val.old)]))]
# outputfile = paste(ml.path,"model-rep10cv5trained-holdout-irl",parameter,fil.par,"-clin-genet.txt",sep="")
# outputfile.val =paste(ml.path,"model-rep10cv5trained-stard-irl",parameter,fil.par,"-clin-genet.txt",sep="")
# write.table(t("Fold	Accuracy	Sensitivity	Specificity	Alg	AccuracyPValue	AccuracyLower	AccuracyUpper	McnemarPValue"), file=outputfile,quote=F, row.names=F,col.names=F,sep="\t",append=T)
# write.table(t("Fol	dAccuracy	Sensitivity	Specificity	Alg	AccuracyPValue	AccuracyLower	AccuracyUpper	McnemarPValue"), file=outputfile.val,quote=F, row.names=F,col.names=F,sep="\t",append=T)
###############
###checks###
###############
# train.snps=train.newdat[6:length(train.newdat)]
# corr=cor(train.snps)
# highlyCorDescr=findCorrelation(corr, cutoff = 0.3)
# train.newdat<- train.newdat[,-highlyCorDescr] 
# houldout.newdat<- houldout.newdat[,-highlyCorDescr]
# outputfile.snps = paste(ml.path,"snps-no-corr-irl",parameter,fil.par,"-.txt",sep="")
# write.table(names(train.newdat), file=outputfile.snps,quote=F, row.names=F,col.names=F,sep="\t",append=F)
# val.newdat<- val.newdat[,-highlyCorDescr]
#comboInfo <- findLinearCombos(filteredCorr )
###############
###fit model###
###############
# fitControl <- trainControl(method = "repeatedcv",  number = 5, classProbs = TRUE,   summaryFunction = twoClassSummary, savePredictions = TRUE, repeats=20)###classic
# ml.methods=c("rpart","C5.0","gbm","naive_bayes","rf","glm","svmRadial","svmLinear")
#     for (alg in 1 : length(ml.methods)){ 
#         model.train <- train(Response_50HAMD ~., data = train.newdat , method = ml.methods[alg],tuneLength = 10,  trControl = fitControl,metric = "ROC")
#         model.pred.hold <- predict(model.train, newdata = houldout.newdat)
#         mc.hold <- confusionMatrix(data =  model.pred.hold, as.factor(houldout.newdat[,1]))
#         model.pred.val <- predict(model.train, newdata = val.newdat)
#         mc.val <- confusionMatrix(data = model.pred.val, val.newdat[,1])
#         measures.hold <- cbind (mc.hold$overall["Accuracy"],mc.hold$byClass["Sensitivity"], mc.hold$byClass["Specificity"] , ml.methods[alg],mc.hold$overall["AccuracyPValue"],mc.hold$overall["AccuracyLower"],mc.hold$overall["AccuracyUpper"],mc.hold$overall["McnemarPValue"] )
#         write.table(measures.hold , file=outputfile,quote=F, row.names=F,col.names=F,sep="\t",append=T)
#         measures.val <- cbind (mc.val$overall["Accuracy"],mc.val$byClass["Sensitivity"], mc.val$byClass["Specificity"] , ml.methods[alg],mc.val   $overall["AccuracyPValue"],mc.val$overall["AccuracyLower"],mc.val$overall["AccuracyUpper"],mc.val$overall["McnemarPValue"] )
#         write.table(measures.val, file=outputfile.val,quote=F, row.names=F,col.names=F,sep="\t",append=T)          
#         jpeg(paste(ml.path,"caret-model-tuning",parameter,fil.par,ml.methods[alg],"-ROC.jpg", sep=""),width=1700, height=1700, units="px",res=300)
#         plot.roc(as.numeric(model.train$pred$obs),as.numeric(model.train$pred$pred))
#         dev.off()
        #write.table( model.train$results,file=paste(ml.path,"caret-model-tuning",parameter,fil.par,ml.methods[alg],".txt", sep=""), quote=F, row.names=F,sep="\t",append=T)
        #g <- ggplot(model.train$pred, aes(m=R, d=factor(obs, levels = c("R", "RN")))) + geom_roc(n.cuts=0) + coord_equal() +style_roc()

# }

#g <- ggplot(model.train$pred, aes(m=R, d=factor(obs, levels = c("RN", "R")))) + geom_roc(n.cuts=0) + coord_equal() +style_roc()
###############
###compare resulls from CV on a real world, get BEST model from previous step###
###############
##model.train <- train(Response_50HAMD ~., data = train.newdat , method = "svmLinear",tuneLength = 10,  trControl = fitControl,metric = "ROC")
###############
###clean a space###
###############
#rm(list=ls())
###double loop
# # # # # result <- foreach(i = 1:5,.packages  =c("glmnet","sda","caret")) %dopar% {
# # # # #     cv.training <- as.data.frame(train.newdat[folds[[i]],])
# # # # #     cv.testing <- as.data.frame(train.newdat[-folds[[i]],])
# write.table(t("Fold	Accuracy	Sensitivity	Specificity	Alg	AccuracyPValue	AccuracyLower	AccuracyUpper	McnemarPValue"), file=outputfile.cv,quote=F, row.names=F,col.names=F,sep="\t",append=T)
# outputfile.cv = paste(ml.path,"model-cv5trained-innercv-irl",parameter,fil.par,"-clin-genet.txt",sep="")