#################################################
#                                               #
#   Process Bayenv2 environmental association   #
#              analysis results                 #
#                                               #
#################################################

##--- LIBRARIES ---##

library(VennDiagram)

##--- FUNCTION ---##

avg.bayenv<-function(dir.prefix,fileIN){

  ## Calculate average of Bayenv2 enviromental association runs ##
  # dir.prefix = directory prefix pattern for series of Bayenv2 runs, including wildcard where appropriate
  # fileIN = name of Bayenv2 output file (assumes same name for all runs)
  
  # create list of files to be averaged, based on directories and input file name given
  files<-Sys.glob(file.path(dir.prefix,fileIN))
  
  # import each result file and add together
  for(i in 1:length(files)){
    assign(paste("tmp",i,sep=""),read.table(files[i],row.names=1,header=F,stringsAsFactors = F))
    if(i==1){
      newdf<-tmp1
    } else
      if(i>1){
        newdf<-newdf+get((paste("tmp",i,sep="")))
      }
  }
  
  # get average by dividing by the number of results files
  avgdf<-newdf/(length(files))
  
  return(avgdf)
}


#----------------------------------------------#
#   Check consistency of covariance matrices   #
#----------------------------------------------#

file.list<-dir(pattern="BALX.2ndfilt.neutral.*.txt")
for(i in 1:length(file.list)){
  tmp<-as.matrix(read.table(paste(file.list[i]),stringsAsFactors = F))
  assign(sub(".txt","",paste(file.list[i])),tmp)
}
 
name="BALX.2ndfilt.neutral."
endings<-c('.covmat.99k','.covmat.99-5k','.covmat.100k')
png("CovMat_Oct2016.png",width=6,height=6,units="in",res=150)
par(mfcol=c(3,3))
for(rep in 1:3){
  for(i in 1:length(endings)){
    tmp=paste(name,rep,endings[i],sep="")
    image(cov2cor(get(tmp)),main=paste(tmp))
  }
}
dev.off()

#- check mean covariance matrix created from average of last 5 draws of run #1

BALX.neutral.mean<-as.matrix(read.table('BALX.2ndfilt.neutral.1.covmat.mean',stringsAsFactors = F))
image(cov2cor(BALX.neutral.mean))

#- check mean covariance matrix created from average of last 5 draws of run #3

BALX.neutral.3.mean<-as.matrix(read.table('BALX.2ndfilt.neutral.3.covmat.mean',stringsAsFactors = F))
image(cov2cor(BALX.neutral.3.mean))

#- create mean covariance matrix from final matrix of 3 runs

BALX.neutral.3.mean <- (BALX.2ndfilt.neutral.1.covmat.100k + BALX.2ndfilt.neutral.2.covmat.100k + BALX.2ndfilt.neutral.3.covmat.100k)/3
image(cov2cor(BALX.neutral.3.mean))

#--------------------------------------------------------------#
#   Get average environmental association results across run   #
#--------------------------------------------------------------#

run1.avg<-avg.bayenv("run1-*","bf_environ.enviro1-subset2.txt")
write.table(run1.avg,file="run1.bf_environ.enviro1-subset2.avg.txt",
            quote=F,row.names=T,col.names=F,sep="\t")


#------------------------------------------------------#
#    Venn Diagram of loci by environmental variable    #
#------------------------------------------------------#

# Fst outlier (2+/4 tests) + strong (BF > 20) environmental association

files<-dir(pattern='run1.*.loci')
run1.snps<-list()
for(i in 1:length(files))
{
  name.tmp<-sub('.loci','',paste(files[i]))
  name<-sub('run1.','',name.tmp)
  tmp<-read.table(files[i],header=F,stringsAsFactors = F)
  run1.snps[[paste(name)]]<-c(tmp[,1])
}

run1.strong<-list(Temperature=run1.snps$strongtemp,Precipitation=run1.snps$strongprec,Aridity=run1.snps$strongarid)
venn.diagram(run1.strong, 'run1-strong.tiff', cex=2.5, cat.cex=2.5, margin=0.2, lwd=4, cat.dist=c(0.1,0.1,0.07), cat.pos = c(335,25,180))

