##########################################################
#                                                        #
#           FDR correction of Lositan output             #
#        using "correct.pval.dataframe" function         #
#                 from getPvalues.R                      #
#   by Lotterhos & Whitlock 2014 (Mol Ecol 23:2178-92)   #
#                                                        #
##########################################################

##--- LIBRARIES ---##
library(qvalue)

##--- FUNCTIONS ---##

#-- From Lotterhos & Whitlock 2014 (Mol Ecol 23:2178-92) getPvalue.R
correct.pval.dataframe <- function(dataframe, p.colName="p.val.cum", p.colNum){
  ##########################################
  ### Correct a list of cumulative p-values to indicate tails,
  ###### and indicate significance at the Bonferroni, FDR=0.05, and FDR=0.01 levels
  ###### Depends on Storey's qvalue.R function
  ### infilepath is path to the file, if it is not yet loaded as a dataframe
  ### assumes the column in the dataframe is named p.val.cum
  ### if the column has another name, please specify name and number of the column 
  ### the p-value column should be based on cumulative probabilities 
  ###  (i.e. starting at 0 in the left tail and ending at 1 in the right tail)
  ### In the output, "L" indicates left-tail, and "R" indicates right-tail	
  ### if write.outfile==TRUE
  ### will write to a file in outfilepath, but with the extension ".Cpval"
  ### if outfilepath is NOT specified, will write to infilepath with the extension ".Cpval"
  ### if neither outfile path or infilepath is specified, will give an error
  
  if (p.colName!="p.val.cum"){
    p.val <- assign(p.colName, dataframe[,p.colNum])
  }else{
    p.val <- dataframe$p.val.cum
  }
  
  L.p <- p.val
  R.p <- 1-L.p
  num.obs <- length(L.p)
  Bonf.p <- 0.05/num.obs
  
  ### Get qvalues for left hand side of distribution
  q2 <- qvalue(L.p)
  L.q <- q2$qvalues
  
  ### Get qvalues for right hand side of distribution	
  q3 <- qvalue(R.p)
  R.q <- q3$qvalues
  
  ### Convert p-values and q-values to the right and left sides of the distribution
  p.val.tail <- L.p
  p.val.tail[L.p>0.5] <- R.p[L.p>0.5] 
  
  qval <- L.q
  maxq <- min(q2$pi0, q3$pi0) 
  #if pi0 is different for the left and right sides, take the minimum
  qval[L.q<maxq] <- L.q[L.q<maxq]
  qval[R.q<maxq] <- R.q[R.q<maxq]
  
  Tail <- L.p
  Tail[L.p>0.5] <- "R"
  Tail[L.p<=0.5] <- "L"
  tail <- as.factor(Tail)
  
  Bonf <- rep(FALSE, num.obs)
  Bonf[p.val.tail<Bonf.p] <- TRUE
  
  FDR.01 <- rep(FALSE, num.obs)
  FDR.01[qval<0.01] <-  TRUE
  
  FDR.05 <- rep(FALSE, num.obs) 
  FDR.05[qval<0.05] <-  TRUE
  
  out <- data.frame(dataframe, tail=tail, p.val.tail=p.val.tail, qval.tail=qval, 
                    Bonf=Bonf, FDR.01=FDR.01, FDR.05=FDR.05)
  
  return(out)
}

##--- IMPORT DATA ---###

test1<-read.table("test1_mcarpa2_2ndfilt_loci.txt",header=T,stringsAsFactors = F)
test2<-read.table("test2_mcarpa2_2ndfilt_loci.txt",header=T,stringsAsFactors = F)
test3<-read.table("test3_mcarpa2_2ndfilt_loci.txt",header=T,stringsAsFactors = F)

##--- FDR CORRECTION ---##

corrected.test1<-correct.pval.dataframe(test1,p.colName="p.val.cum")
write.table(corrected.test1,file="mcarpa2.2ndfilt_Test1_loci_plus_q.txt",quote=F,row.names=F)

corrected.test2<-correct.pval.dataframe(test2,p.colName="p.val.cum")
write.table(corrected.test2,file="mcarpa2.2ndfilt_Test2_loci_plus_q.txt",quote=F,row.names=F)

corrected.test3<-correct.pval.dataframe(test3,p.colName="p.val.cum")
write.table(corrected.test3,file="mcarpa2.2ndfilt_Test3_loci_plus_q.txt",quote=F,row.names=F)


