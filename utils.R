fdrQvalue <- function(results, FDRthreshold){
  results = results[!is.na(results$P),]
  totalN = length(results$P)
  data = data.frame(Rank=as.numeric(1:totalN), P=sort(results$P))

  # calculate q-values from p-values: qvalue = pi0*m*p/v
  data$Q = (totalN*data$P)/data$Rank

  # FDR Method:
  # find  Q value under FDR threshold, then find corresponding maximum P value
  # find number of snps with p value under maxium
  # calculate percent of significant SNPs
  maxP = max(data[data$Q < FDRthreshold,]$P)
  NSigSnp=length(data[data$P <= maxP, ]$P)
  PerSigSnp=round((NSigSnp/totalN)*100,4)
  output=data.frame(NSigSnps=NSigSnp,
                    NTotalSnps=totalN,
                    PerSigSnps=PerSigSnp)
  output
}

qqunif <- function(p,BH=T,CI=T,FDRthres=0.05, title=""){
   nn = length(p)
   xx =  -log10((1:nn)/(nn+1))
   dat<-cbind(sort(p),1:nn)
   q<-(nn*dat[,1])/dat[,2] # calculate q-values from p-values
   dat<-cbind(dat,q)

   plot(xx,  -sort(log10(p)), cex.lab=1, mgp=c(2,1,0), pch=20, xlab="Expected -log(10)", ylab="Observed -log(10)", main=title)

   if(CI) {
     c95 <- rep(0,nn)
     c05 <- rep(0,nn)

     for(i in 1:nn)
     {
       c95[i] <- qbeta(0.95,i,nn-i+1)
       c05[i] <- qbeta(0.05,i,nn-i+1)
     }
     polygon(c(xx, rev(xx)), c(-log10(c95), rev(-log10(c05))),col = "grey", border = NA)

   }

   abline(0,1,col='red')
   if(BH) {
     abline(-log10(0.05),1, col='black',lty=2)
     abline(-log10(0.10),1, col='black',lty=3)
     abline(-log10(0.25),1, col='black',lty=4)
     abline(-log10(FDRthres),1, col='blue',lty=5)
     abline(h=-log10(0.05/nn),col="blue") ## bonferroni
     legend('bottomright', c("FDR = 0.05","FDR = 0.10","FDR = 0.25", paste("FDR =",FDRthres)),col=c('black','black','black', 'blue'),lty=2:5, cex=0.7)
   }

   points(xx,  -sort(log10(p)), cex.lab=1, mgp=c(2,1,0), pch=20)

   #FDR
   nsnps<-round((sum(p<=max(dat[dat[,3]<FDRthres,1]))/nn)*100,4)
   y<-max(-log10(p))
   text(0,y,paste(nsnps,"% of SNPs have a q-value <= ",FDRthres,sep=""),pos=4)
   print(paste(nsnps,"% of SNPs have a q-value <= ",FDRthres,sep=""))
 }

MAFfilter <- function(alldata,  MAFthreshold, MAFfile){
  #alldata = dataframe object with required column "SNP"
  #MAFfile = .frq file for dataset
  #NOTE: returns all SNPs >= MAFthreshold
mafDf=read.table(file=MAFfile, header=T, colClasses=c("NULL", "character", "NULL", "NULL", "numeric","NULL"))
combinedDf=merge(alldata,mafDf, by="SNP")
output = combinedDf[combinedDf$MAF >= MAFthreshold,]
output
}

qqplot2datasets <- function(df1, df2, FDR, title){
  #data 1
  results1 = data.frame(obs=-log10(sort(df1$P, decreasing=F)),
                        exp=-log10( 1:length(df1$P)/length(df1$P)))
  r1FDR = fdrQvalue(df1, FDR)
  #data 2
  results2 = data.frame(obs=-log10(sort(df2$P, decreasing=F)),
                        exp=-log10( 1:length(df2$P)/length(df2$P)))
  r2FDR = fdrQvalue(df2, FDR)
  p <- (ggplot() +
          geom_point(results1, mapping=(aes(x=exp, y=obs, color="Data 1"))) +
          geom_point(results2, mapping=(aes(x=exp, y=obs, color="Data 2"))) +
          geom_abline(slope=1, intercept=0, color="black") +
          labs(title=title) +
          xlab("Observed -log(p)") +
          ylab("Expected -log(p)") +
          guides(colour = guide_legend(title = "Data Sets")) +
          geom_text(aes(label = paste(r1FDR$PerSigSnps,"% of SNPs <= ",FDR," FDR",sep=""),
                   x = 2, y = ceiling(max(results1$obs)),color="Data 1")) +
          geom_text(aes(label = paste(r2FDR$PerSigSnps,"% of SNPs <= ",FDR," FDR",sep=""),
                    x = 2, y = ceiling(max(results1$obs)-1),color="Data 2"))
)
  p
}

##added 12/22/2015
histSigSNPs <- function(df, title, originalper, FDR){

      n = length(df[df$PerSigSnps >= originalper, 1])
      hist(df$PerSigSnps, 
           col="darkgray",
           main=paste(title,"\n",n,"perm sets >= True proportion", originalper),
           xlab= paste("% of Significant SNPs (FDR", FDR,")"),
           breaks=20
      )
      abline(v=originalper, col="red", lty=2)
}
