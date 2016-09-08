###function for pulling variant windows
plucking <- function(input_file, topHits_file, output_file, dataset) {

 
 
  #########################################################################
  ### RANK SNPS BY THEIR FIXED EFFECT P-VALUE (OR OTHER TEST STATISTIC) ###
  #########################################################################
 
  #####garden <- read.table(input_file, header=T)
  ###REPLACED above line with input format for metasoft results
  garden <-  read.table(input_file, header=F, skip=1, stringsAsFactors=FALSE)
  colnames(garden)[1:16] =  c("RSID", "NUM_STUDY", "PVALUE_FE", "BETA_FE","STD_FE", "PVALUE_RE", "BETA_RE", "STD_RE", "PVALUE_RE2", "STAT1_RE2", "STAT2_RE2", "PVALUE_BE", "I_SQUARE", "Q", "PVALUE_Q", "TAU_SQUARE")
  
  garden <- as.data.frame(garden)
  garden <- garden[,c("RSID", "PVALUE_FE")] ##added
  garden$RSID <- as.character(garden$RSID)
  garden$PVALUE_FE <- as.numeric(garden$PVALUE_FE)
 
  noduplicates.garden <- garden[!duplicated(garden[,1:1]),]
  garden <- noduplicates.garden
 
  growth <- order(garden$PVALUE_FE)
  growth
  garden <- garden[growth,]
  garden
 
  rownames(garden) = c(1:nrow(garden))
 
  ###replace rownames(garden), which has the positions of the former rows they were
  ###sorted in, to instead have a linear order of 1 through nrows of garden.
 
  garden

 
  flowers <- read.csv(topHits_file, header=T) ##changed from read.table
  flowers <- flowers[,c("RSID", "PVALUE_FE")] ###ADDED
  flowers <- as.data.frame(flowers)
  flowers$RSID <- as.character(flowers$RSID)
  flowers$PVALUE_FE <- as.numeric(flowers$PVALUE_FE)
 
  noduplicates.flowers <- flowers[!duplicated(flowers[,1:1]),]
  flowers <- noduplicates.flowers
 
  growth <- order(flowers$PVALUE_FE)
  growth
  flowers <- flowers[growth,]
  flowers
 
  rownames(flowers) = c(1:nrow(flowers))
 
  flowers
 
  for (i in 1:nrow(flowers)) {
   
    blossom <- flowers[i,]
    flowerInGarden <- garden[garden$RSID == blossom$RSID,]
   
    if( i == 1 ) {
     
      location <- rownames(flowerInGarden)
      location <- as.numeric(location)
      location
     
      neighboringFlowers <- garden[(max((location-75),1):min((location+75),nrow(garden))),]  
      plucked <- neighboringFlowers
      plucked
      plucked[,"DATASET"] <- c(dataset)     

     
      write.table(plucked, paste(output_file,".txt",sep=""), append = TRUE, row=F, quote=F, col.names=T)
     
    }
   
    if (i > 1)
    {
      location <- rownames(flowerInGarden)
      location <- as.numeric(location)
      location
     
      neighboringFlowers <- garden[(max((location-75),1):min((location+75),nrow(garden))),]  
      plucked <- neighboringFlowers
      plucked[,"DATASET"] <- c(dataset) 
     
     
      write.table(plucked, paste(output_file, ".txt",sep=""), append = TRUE, row=F, quote=F, col.names=F)
     
    }
   
   
  }
}
