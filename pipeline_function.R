source("/path/utils.R")##if needed change this to the correct path.

createGenomeFileSnpCount <- function(refGenesFile, kbpRange, mafFile, MAFthresh, bimFile, outputGenomeFile){
    ###INPUT:
    #refGenomeFile: tab delimited file that MUST have columns "GENE" "CHR" "START" "END", where the CHR, START, END columns are only numeric and are the true values.
    #kbpRange is the number of bases to + and - from start and end of genes.
    #mafFile is the .frq file produced by plink --freq. Leave blank ("") and set MAFthresh to 0, if do not want to use filter.
    #bimFile is the .bim file for the dataset in use
    #outputGenomeFile will be written as tab delimited file with only genes that have atleast 1 SNP in dataset, and new START, END, LENGTH,and SNPCOUNT.
    refGenes = read.table(refGenesFile, header=T, fill=F)
    refGenes$START <- refGenes$START - kbpRange
    refGenes$END <- refGenes$END + kbpRange
    refGenes$LENGTH <-  refGenes$END - refGenes$START

    resultSnps = read.table(file=bimFile, header=F, fill=F, row.names=NULL, col.names=c("CHR", "SNP", "CM", "BP", "A1", "A2"))
    if(MAFthresh > 0) {
        resultSnps = MAFfilter(resultSnps,  MAFthresh, mafFile)
    }

    refGenes$SNPCOUNT = numeric(length=length(refGenes$GENE))
    for (i in seq(from=1, to=length(refGenes$GENE), by=1)) {
        refGenes$SNPCOUNT[i] = length(resultSnps[(resultSnps$CHR == refGenes$CHR[i] & refGenes$START[i] <= resultSnps$BP & resultSnps$BP <= refGenes$END[i]), "SNP"])
    }
    refGenes = refGenes[refGenes$SNPCOUNT > 0, ]
    write.table(refGenes, file=outputGenomeFile, quote=F,  row.names=F, sep="\t")
    refGenes
}

getOverlapSnps <- function(resultsSnps, gseaSnps, kbp_range){
    gseaSnps$L1 <- gseaSnps$BP-kbp_range
    gseaSnps$L2 <- gseaSnps$BP+kbp_range

    overlap <- data.frame()
    for (i in seq(from=1, to=length(gseaSnps$SNP), by=1)) {
        #filter for snps in same chromosome
        df = resultsSnps[resultsSnps$CHR == gseaSnps$CHR[i],]
        #filter for snps within position range
        df = df[(gseaSnps$L1[i] <= df$BP & df$BP <= gseaSnps$L2[i]), ]
        overlap <- rbind(overlap, df)
    }
    overlap <- with(overlap, overlap[order(overlap$SNP),])
    overlap <- overlap[!duplicated(overlap$SNP), ]
    overlap
}

getOverlapGenes <- function(resultsSnps, gseaGenes, kbp_range){
    gseaGenes$L1 <- gseaGenes$START - kbp_range
    gseaGenes$L2 <- gseaGenes$END + kbp_range

    overlap <- data.frame()
    for (i in seq(from=1, to=length(gseaGenes$GENE), by=1)) {
        df = resultsSnps[resultsSnps$CHR == gseaGenes$CHR[i],]
        df = df[(gseaGenes$L1[i] <= df$BP & df$BP <= gseaGenes$L2[i]), ]
        df$GENE = character(length=length(df$SNP))
        df$GENE = rep(gseaGenes$GENE[i], length(df$SNP))
        overlap <- rbind(overlap, df)
    }
    overlap = with(overlap, overlap[order(SNP),])
    overlap = overlap[!duplicated(overlap$SNP), ]
    rownames(overlap) = NULL
    overlap
}

genePermutationsByLength <- function(permN, refGenesFile, setGenesFile, outputPermFile){
    ###INPUT:
    #permN = number of permutations
    #refGenesFiles = txt file produced by createGenomeFileSnpCount
    #setGenesFile = txt listing genes of interest
    #outputPermFile = file of output text permutation gene files ending in "_" (such as "/home/gsea/myGenePerms_")
    refGenes <- read.table(refGenesFile, header=T, row.names=NULL)
    refGenes <- refGenes[order(refGenes$LENGTH),]
    refGenes <- refGenes[!(duplicated(refGenes$GENE)),]
    refGenesN <- nrow(refGenes)
    rownames(refGenes) <- 1:refGenesN #row.name = rank

    setGenes <- read.table(setGenesFile, header=T, row.names=NULL) #must have header with "GENE"
    setGenes <- unique(setGenes$GENE)
    setGenesN <- length(setGenes)

    for (n in 1:permN){
        permGenes <- data.frame()
        for (gene in setGenes){
            pos <- as.numeric(rownames(refGenes[refGenes$GENE == paste(gene),]))
            lowerPos <- max(c((pos-50), 1))
            upperPos <- min(c((pos+50), refGenesN))
            neighbors <- refGenes[lowerPos:upperPos, ]
            neighbors <- neighbors[neighbors$GENE != paste(gene),]
            permGene <- sample(neighbors$GENE, 1, replace=FALSE)
            permGenes <- rbind(permGenes, refGenes[refGenes$GENE == paste(permGene), ])
        }
        write.table(permGenes, file=paste(outputPermFile, n, ".txt",sep=""), quote=F,  row.names=F, sep="\t")
    }
}

genePermutationsBySNPCount <- function(permN, refGenesFile, setGenesFile, outputPermFile){
    ###INPUT:
    #permN = number of permutations
    #refGenesFiles = txt file produced by createGenomeFileSnpCount
    #setGenesFile = txt listing genes of interest
    #outputPermFile = file of output text permutation gene files ending in "_" (such as "/home/gsea/myGenePerms_")
    refGenes <- read.table(refGenesFile, header=T, row.names=NULL)
    refGenes <- refGenes[order(refGenes$SNPCOUNT),]
    refGenes <- refGenes[!(duplicated(refGenes$GENE)),]
    refGenesN <- nrow(refGenes)
    rownames(refGenes) <- 1:refGenesN #row.name = rank

    setGenes <- read.table(setGenesFile, header=F, row.names=NULL, col.names=c("GENE"))
    setGenes <- unique(setGenes$GENE)
    setGenesN <- length(setGenes)

    for (n in 1:permN){
        permGenes <- data.frame()
        for (gene in setGenes){
            pos <- as.numeric(rownames(refGenes[refGenes$GENE == paste(gene),]))
            lowerPos <- max(c((pos-50), 1))
            upperPos <- min(c((pos+50), refGenesN))
            neighbors <- refGenes[lowerPos:upperPos, ]
            neighbors <- neighbors[neighbors$GENE != paste(gene),]
            permGene <- sample(neighbors$GENE, 1, replace=FALSE)
            permGenes <- rbind(permGenes, refGenes[refGenes$GENE == paste(permGene), ])
        }
        write.table(permGenes, file=paste(outputPermFile, n, ".txt",sep=""), quote=F,  row.names=F, sep="\t")
    }
}


permGeneEnrichAnalysis <- function(startPermN, endPermN, genePermFile, FDRthresh, resultsDataFrame, outCSVFile){
    ###INPUT:
    #startPermN, endPermN = the numbers that genePermFile ends in (such as 1, 100)
    #genePermFile = same as outputPermFile in genePermutationsByLength function
    #FDRthresh = number threshold for calculating percent significant SNPs.
    #resultsDataFrame = data frame object with columns SNP (text) CHR (numeric) BP (numeric) P (numeric)
    #outCSVFile = .csv file with significance results
    #create output df
    outputdf = data.frame(NSigSNPs=numeric(), NTotalSNPs=numeric(), PerSigSNPs=numeric())
    for (n in startPermN:endPermN) {
        #read in genes list
        gseaGenes = read.table(file=paste(genePermFile,n,".txt",sep=""), header=T, fill=F,row.names=NULL)
        #find SNPs from results that overlap with perm Genes
        overlap = getOverlapGenes(resultsDataFrame, gseaGenes, 0)
        #calculate % Significant SNPs at FDR threshold
        fdrOut = fdrQvalue(overlap, FDRthresh)
        rownames(fdrOut) = n
        outputdf = rbind(outputdf, fdrOut)
    }
    write.table(outputdf, file = outCSVFile, sep=",", quote=F, row.names=T)
    outputdf
}

processMetaRE2ResultsDF <- function(metaResultsFile, snpMapFile, MAFthresh, MAFFile){
    ###Output dataframe will have columns SNP, CHR, BP, P, MAF
    #read in Metasoft .mmap file
    snpMap=read.table(file=snpMapFile, header=F, fill=F, row.names=NULL, colClasses=c("character", rep("numeric", 2), rep("NULL", 2)))
    colnames(snpMap)=c("SNP", "CHR", "BP")
    #read in Metasoft RE2 results
    metaResults = read.table(file=metaResultsFile, header=F, skip=1, fill=T)
    metaResults = metaResults[,c(1,9)] #SNP and RE2 Pvalue cols
    colnames(metaResults) = c("SNP", "P")
    rownames(metaResults) = NULL
    #MAF filter
    metaResults=MAFfilter(metaResults,  MAFthresh, MAFFile)
    metaResults=metaResults[!is.na(metaResults$P),]
    metaResults=merge(metaResults, snpMap, by="SNP")
    metaResults
}
