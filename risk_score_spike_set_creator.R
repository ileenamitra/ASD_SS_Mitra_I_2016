OUTDIR="/outdir_path/"
###TAB FILE WITH FID, IID, SEX, PHENO, BATCH header
peopleFile = "/path/ASD_SS_HeritabilityPerm_input_people_12082015.info"

people = read.table(file = peopleFile, header=T, colClasses=c(rep("character",2), rep("numeric", 2)))

femaleCases = people[(people$SEX == 2 & people$PHENO == 2),]
NfemaleCases = nrow(femaleCases) 

femaleControls = people[(people$SEX == 2 & people$PHENO == 1),]
NfemaleControls = nrow(femaleControls) 

maleCases = people[(people$SEX == 1 & people$PHENO == 2),]
NmaleCases = nrow(maleCases)

maleControls = people[(people$SEX == 1 & people$PHENO == 1),]
NmaleControls = nrow(maleControls) 

print(paste("# Total People", nrow(people), "; #Female Cases", NfemaleCases,"; Female Controls", (nrow(people) - NfemaleCases - NmaleCases -  NmaleControls), "; #Male Cases", NmaleCases,"; #Male Controls", NmaleControls))

####CREATE RISK SCORE SET ####
###for each female case - find male case froms same batch.
###to avoid technical batch issues
###[1] "AGP_imputed_QC+", "SS_imputed_QC+", "CHARGE_imputed_QC+_update_pheno", "SSC1M_imputed_QC+""SSCDuo_imputed_QC+", "SSCOmni_imputed_QC+""AGRE_imputed_QC+", "CHOP_imputed_QC+"
testMaleCases= data.frame()
for (batch in unique(people$BATCH)){
	#in batchA , N number of female cases
	permCasesBatchN = nrow(femaleCases[femaleCases$BATCH==batch,])
	print(paste(batch, permCasesBatchN, sep=": "))
	###add N male cases from batchA as females
	permMaleCasesIID = sample(maleCases[maleCases$BATCH==batch ,"IID"], permCasesBatchN, replace = FALSE)
	temp =  maleCases[((maleCases$IID %in% permMaleCasesIID) & !(maleCases$IID %in% testMaleCases$IID)), ]
	testMaleCases <- rbind(testMaleCases, temp)
}

testMaleControlsIID = gsub("_T$", "_U", testMaleCases$IID)
testMaleControls = maleControls[maleControls$IID %in% testMaleControlsIID, ]

if ((nrow(testMaleCases) == NfemaleCases) & (nrow(testMaleControls)==NfemaleControls)) {
		
		testMalesSet = rbind(testMaleCases, testMaleControls)
		write.table(testMalesSet[,c("FID", "IID", "PHENO")],
					file=paste(OUTDIR, "RiskScore/","testMalesSet_12092015_TotalN",nrow(testMalesSet),".pheno", sep=""), 
					quote=F, row.names=F, col.names = F)

		discMalesSet = people[(!(people$IID %in% testMalesSet$IID) & people$SEX == 1), ]
		write.table(discMalesSet[,c("FID", "IID", "PHENO")],
					file=paste(OUTDIR, "RiskScore/","discoveryMalesSet_12092015_TotalN",nrow(discMalesSet),".pheno", sep=""), 
					quote=F, row.names=F, col.names = F)

		testFemalesSet = people[people$SEX == 2, ]
		write.table(testFemalesSet[,c("FID", "IID", "PHENO")],
					file=paste(OUTDIR, "RiskScore/","testFemalesSet_12092015_TotalN",nrow(testFemalesSet),".pheno", sep=""), 
					quote=F, row.names=F, col.names = F)		
	for (NpermFemaleCases in  c(seq(0,NfemaleCases, by=10),NfemaleCases)){
			#randomly select N female cases
			permFemaleCasesIID = sample(femaleCases$IID, NpermFemaleCases, replace = FALSE)
			permFemaleCases = femaleCases[femaleCases$IID %in% permFemaleCasesIID,]
			##add N male cases from test set - with matching batch as female case
			permMaleCases= data.frame()
			###to avoid technical batch issues
			for (batch in unique(people$BATCH)){
				#in batchA , N number of female cases
				permCasesBatchN = nrow(permFemaleCases[permFemaleCases$BATCH==batch,])
				###add N male cases from batchA as females
				permMaleCasesIID = sample(testMaleCases[testMaleCases$BATCH==batch ,"IID"], permCasesBatchN, replace = FALSE)
				temp =  testMaleCases[((testMaleCases$IID %in% permMaleCasesIID) & !(testMaleCases$IID %in% permMaleCases$IID)), ]
				permMaleCases <- rbind(permMaleCases, temp)
			}
			##add pseudo controls matching cases
			##create male and female sets of spiked perm set 
			permMaleControlsIID = gsub("_T$", "_U", permMaleCases$IID)
			permMaleControls = maleControls[maleControls$IID %in% permMaleControlsIID, ]
			finalPermMale = rbind(permMaleCases, permMaleControls)
			finalPermMale = rbind(finalPermMale, discMalesSet)

			permFemaleControlsIID = gsub("_T$", "_U", permFemaleCases$IID)
			permFemaleControls = femaleControls[femaleControls$IID %in% permFemaleControlsIID, ]
			finalPermFemale = rbind(permFemaleCases, permFemaleControls)
			finalPermFemale = rbind(finalPermFemale, discMalesSet)
			
			###check final spiked perm set
			if (nrow(finalPermMale) == ((2*NpermFemaleCases)+(nrow(discMalesSet))) & 
				nrow(finalPermFemale) == ((2*NpermFemaleCases)+(nrow(discMalesSet)))){
				if ((sum(duplicated(finalPermMale$IID))==0) & (sum(duplicated(finalPermFemale$IID))==0)){
					#write perm
					write.table(finalPermMale[,c("FID", "IID", "PHENO")],
					file=paste(OUTDIR, "SpikedMaleH2Perm/", "GCTA_spiked_perm_Nmales_", nrow(permMaleCases), ".pheno", sep=""), 
					quote=F, row.names=F, col.names = F)

					write.table(finalPermFemale[,c("FID", "IID", "PHENO")],
					file=paste(OUTDIR,"SpikedFemaleH2Perm/", "GCTA_spiked_perm_Nfemales_", nrow(permFemaleCases), ".pheno", sep=""), 
					quote=F, row.names=F, col.names = F)
				}
				else{ 
					print("Duplicate IIDs found") 
					break
				}
		
			}
			else{
				print("ERROR total number of permuted people do not match.")
				print(paste("Total PERM female", nrow(finalPermFemale)))
				print(paste("Total PERM male", nrow(finalPermMale)))
				print(paste("Expected number:", (2*NpermFemaleCases)+(nrow(discMalesSet))))
				break 	
			}

	}
}else{print(paste("ERROR:", "Number of test male cases", nrow(testMaleCases), "does not equal", NfemaleCases,
	"and Number of test male controls",nrow(testMaleControls), "does not equal",NfemaleControls))}
