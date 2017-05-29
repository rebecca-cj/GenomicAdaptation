##############################################
#                                            #
#   Perform correction for multiple testing  #
#    on Arlequin outlier (FDIST2) results    #
#                                            #
#    Rebecca Jordan July 2016                #
#                                            #
##############################################

library(qvalue)

arl.data<-read.table('mcarpa2.2ndfilt.fdist2_ObsOut.txt',header=T,stringsAsFactors = F)

arl.BH.adj<-qvalue(p=arl.data$FSTP.value,lambda=0,fdr.level=0.05)
length(subset(arl.BH.adj$qvalues,arl.BH.adj$qvalues<0.05))
length(subset(arl.BH.adj$qvalues,arl.BH.adj$qvalues<0.1))

arl.data.qval<-cbind(arl.data,BHqval=arl.BH.adj$qvalue)
write.table(arl.data.qval,file="mcarpa2.2ndfilt.arl.qval.txt",quote=F,row.names=F,sep="\t")

