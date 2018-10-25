#!	/bin/sh
library(ggplot2)
library(psych)
###############
####Read data and set parametres####
###############
dat<-read.table	("./IRL_training/response/random-labels-internal-rep10cv5-CAT0.8-trained-irlld0.8-maf0.05-phred=2-only-genet.txt",	header=T,	sep="\t") 
dat[dat=="NA"]=NA
dat[dat==" "]=NA
###############
###external function for multiple plots of ggplot2###
################
###############
###this is a code downloaded from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/###
###############
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
#   # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {  
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
 }
###############
###plot summary of repeated CV###
################
p.acc <- ggplot(dat, aes(x=Alg, y=Accuracy, fill=Alg))+ 
geom_boxplot(width=0.5)+
theme_bw() + 
theme(legend.position="none",text = element_text(size=14))+ 
xlab("")+ 
scale_fill_brewer(palette="Blues")
p.sens<- ggplot(dat, aes(x=Alg, y=Sensitivity, fill=Alg))+ 
geom_boxplot(width=0.5)+
theme_bw() + 
theme(legend.position="none",text = element_text(size=14))+
xlab("")+ 
scale_fill_brewer(palette="Blues")+
ylim(c(0, 1))
p.spec <-ggplot(dat, aes(x=Alg, y=Specificity, fill=Alg))+ 
geom_boxplot(width=0.5)+
theme_bw() + 
theme(legend.position="none",text = element_text(size=14))+ 
xlab("")+ylim(c(0, 1))+ 
scale_fill_brewer(palette="Blues")
jpeg("./IRL_training/response/random-labels-internal-rep10cv5-CAT0.8-trained-irlld0.8-maf0.05-phred=2-only-genet.jpg", width=4400, height=2200,units="px",res=300)
multiplot(p.acc,p.sens,p.spec,cols=2)
dev.off()
###############
###get summary stats###
################
a = describeBy(dat[,2:4], dat$Alg,digits=3)
a2 <- do.call("rbind",a)
write.table(a2,file="./IRL_training/response/random-labels-internal-rep10cv5-CAT0.8-trained-irlld0.8-maf0.05-phred=2-only-genet-summary-stats.txt", quote=F, sep="\t")
