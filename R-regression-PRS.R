#!	/bin/sh
Rargs="--no-save	-q"
library(car)
library(QuantPsyc)
dat.clin=read.table("clinical-data.txt", 	header=TRUE,	sep="\t")
thresholds=c("0.1","0.05","0.01","0.00005","0.00000005") ### thresholds like in profies
outfile="linear-regression-results.txt"
outfile.log="logistic-regression-results.txt"
write.table(t(c("Thres.score" ,"Estimate", "Std.Error",   "t-value",   "P","adj.r.squared","BETA" )), file=outfile, row.names=F,  col.names=F,quote=F, sep="\t", append=T)####linear
write.table(t(c("Thres.score" ,"Estimate", "SE", "z", "P","OR","2.5%", "97.5%")), file=outfile.log, row.names=F,  col.names=F,quote=F, sep="\t", append=T) ###logistic
for (i in 1:length(thresholds)){
    dat.genet = read.table (paste("profile",thresholds[i],".profile",sep=""),header=T)
    dat=merge(dat.clin, dat.genet, by=c("FID","IID"))
    dat[dat=="NA"]=NA
    set.seed(100)
    a.lm<-lm(continues_phenotype~scale(SCORE)+covariate1+covariate2,data=dat, na.action=na.omit)
    reg.sum=summary(a.lm)
    coef1 <- lm.beta(a.lm)
    write.table(t(c(thresholds[i], reg.sum$coefficients["scale(SCORE)",c(1,2,3,4)], reg.sum$adj.r.squared,coef1["scale(SCORE)"])), file=outfile, row.names=F,   col.names=F, quote=F, sep="\t", append=T)######linear
    a.glm<-glm(binary_phenotype ~ cale(SCORE)+covariate1+covariate2,data=dat, family="binomial",na.action=na.omit)
    glm.sum = summary(a.glm)
    or=exp(cbind(OR = coef( a.glm), confint( a.glm)))
    write.table(t(c(thresholds[i],glm.sum$coefficients["scale(SCORE)",c(1,2,3,4)],or["scale(SCORE)",c(1,2,3)])), file=outfile.log, row.names=F, col.names=F, quote=F, sep="\t", append=T)
}