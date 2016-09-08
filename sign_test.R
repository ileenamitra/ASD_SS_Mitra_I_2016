getTopAutosomesSNPs <- function(bimDF, DIR, data, nSNPs){
		#read in results
		topSnps = read.table(paste(DIR, data,"_ASD_results_11202015.output", sep=""), header=F, skip=1, stringsAsFactors=FALSE)
		colnames(topSnps)[1:16] =  c("SNP", "NUM_STUDY", "PVALUE_FE", "BETA_FE","STD_FE", "PVALUE_RE", "BETA_RE", "STD_RE", "PVALUE_RE2", "STAT1_RE2", "STAT2_RE2", "PVALUE_BE", "I_SQUARE", "Q", "PVALUE_Q", "TAU_SQUARE")
		topSnps = topSnps[,c("SNP", "PVALUE_FE", "BETA_FE")]
		# remove NAs
		topSnps = topSnps[!is.na(topSnps$PVALUE_FE),]
		topSnps = topSnps[topSnps$PVALUE_FE != "NAN",]
		topSnps$SNP = as.character(topSnps$SNP)
		topSnps$PVALUE_FE = as.numeric(topSnps$PVALUE_FE)
		topSnps$BETA_FE = as.numeric(topSnps$BETA_FE)
    	
 		#filter for chromosomes
		topSnps = merge(topSnps, bimDF, by.x="SNP", all.x=T)
		topSnps = topSnps[topSnps$CHR < 23,]

		###filter for top SNPs
		topSnps = topSnps[order(topSnps$PVALUE_FE),]
		topSnps = topSnps[1:nSNPs,]
		topSnps
}

getTopChrXSNPs <- function(bimDF, DIR, data, nSNPs){
		#read in results
		topSnps = read.table(paste(DIR, data,"_ASD_results_11202015.output", sep=""), header=F, skip=1, stringsAsFactors=FALSE)
		colnames(topSnps)[1:16] =  c("SNP", "NUM_STUDY", "PVALUE_FE", "BETA_FE","STD_FE", "PVALUE_RE", "BETA_RE", "STD_RE", "PVALUE_RE2", "STAT1_RE2", "STAT2_RE2", "PVALUE_BE", "I_SQUARE", "Q", "PVALUE_Q", "TAU_SQUARE")
		topSnps = topSnps[,c("SNP", "PVALUE_FE", "BETA_FE")]
		# remove NAs
		topSnps = topSnps[!is.na(topSnps$PVALUE_FE),]
		topSnps = topSnps[topSnps$PVALUE_FE != "NAN",]
		topSnps$SNP = as.character(topSnps$SNP)
		topSnps$PVALUE_FE = as.numeric(topSnps$PVALUE_FE)
		topSnps$BETA_FE = as.numeric(topSnps$BETA_FE)
    	
 		#filter for chromosomes
		topSnps = merge(topSnps, bimDF, by.x="SNP", all.x=T)
		topSnps = topSnps[topSnps$CHR == 23,]

		###filter for top 600
		topSnps = topSnps[order(topSnps$PVALUE_FE),]
		topSnps = topSnps[1:nSNPs,]
		topSnps
}

getCompareSNPs <-function (DIR, data){
		compareSNPs = read.table(paste(DIR, data,"_ASD_results_11202015.output", sep=""), header=F, skip=1, stringsAsFactors=FALSE)
		colnames(compareSNPs)[1:16] =  c("SNP", "NUM_STUDY", "PVALUE_FE", "BETA_FE","STD_FE", "PVALUE_RE", "BETA_RE", "STD_RE", "PVALUE_RE2", "STAT1_RE2", "STAT2_RE2", "PVALUE_BE", "I_SQUARE", "Q", "PVALUE_Q", "TAU_SQUARE")
		compareSNPs = compareSNPs[,c("SNP", "PVALUE_FE", "BETA_FE")]
		# remove NAs
		compareSNPs = compareSNPs[!is.na(compareSNPs$PVALUE_FE),]
		compareSNPs = compareSNPs[compareSNPs$PVALUE_FE != "NAN",]
		compareSNPs$SNP = as.character(compareSNPs$SNP)
		compareSNPs$PVALUE_FE = as.numeric(compareSNPs$PVALUE_FE)
		compareSNPs$BETA_FE = as.numeric(compareSNPs$BETA_FE)
    	compareSNPs
}

compareSignBetweenSetsTest <- function(topSNPs, compareSNPs){
    	#change col names to not be confused
    	colnames(topSNPs)[1:3] = c("SNP", "PVALUE_FE_1", "BETA_FE_1")
    	colnames(compareSNPs) = c("SNP", "PVALUE_FE_2", "BETA_FE_2")
    	#compare to TopSNPs
    	betaResults = merge(topSNPs[,c("SNP", "BETA_FE_1")], compareSNPs[,c("SNP", "BETA_FE_2")], by.x="SNP", all.X=T)
    	betaResults$SIGN1 = ifelse(betaResults$BETA_FE_1 < 0, 0, 1)
    	betaResults$SIGN2 = ifelse(betaResults$BETA_FE_2 < 0, 0, 1)
    	betaResults$COMPARE = betaResults$SIGN1==betaResults$SIGN2
    	##output
		out = data.frame(numSameBeta = sum(betaResults$COMPARE), perSameBeta = (sum(betaResults$COMPARE)/nrow(betaResults)*100))
		out

}

bim = read.table("/path/plink.bim", header=F, colClasses = c("numeric", "character", "NULL", "numeric", "NULL", "NULL"))
colnames(bim) = c("CHR", "SNP", "BP")


DIR="/dir_path/"

###AUTOSOMES
###TOP MALE Results Compared to females results
mTop = getTopAutosomesSNPs(bim, DIR, "MALE", 500)
fAll = getCompareSNPs(DIR, "FEMALE")
mTopSignTest = compareSignBetweenSetsTest(mTop, fAll)
print("Autosomes TOP MALE Results Compared to females results")
print(mTopSignTest)
####TOP FEMALE Results Compared to male results	
fTop = getTopAutosomesSNPs(bim, DIR, "FEMALE", 500)
mAll = getCompareSNPs(DIR, 'MALE')
fTopSignTest = compareSignBetweenSetsTest(fTop, mAll)
print("Autosomes TOP FEMALE Results Compared to male results")
print(fTopSignTest)

###CHR X
###TOP MALE Results Compared to females results
mTop = getTopChrXSNPs(bim, DIR, "MALE", 50)
fAll = getCompareSNPs(DIR, "FEMALE")
mTopSignTest = compareSignBetweenSetsTest(mTop, fAll)
print("CHR X TOP MALE Results Compared to females results")
print(mTopSignTest)
####TOP FEMALE Results Compared to male results	
fTop = getTopChrXSNPs(bim, DIR, "FEMALE", 50)
mAll = getCompareSNPs(DIR, "MALE")
fTopSignTest = compareSignBetweenSetsTest(fTop, mAll)
print("CHR X TOP FEMALE Results Compared to male results")
print(fTopSignTest)
