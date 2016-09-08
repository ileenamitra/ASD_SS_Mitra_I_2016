###Protocol:
###Create sex-perms Trios:
###Input R args from shell script:
### Trio .fam file
TrioInput[1]=AGP_Trio_pheno_input_sexPerm.fam
TrioInput[2]=AGRE_Trio_pheno_input_sexPerm.fam
TrioInput[3]=CHARGE_Trio_pheno_input_sexPerm.fam
TrioInput[4]=CHOP_Trio_pheno_input_sexPerm.fam
TrioInput[5]=SS_groups_Trio_pheno_input_sexPerm.fam
TrioInput[6]=SSC1M_Trio_pheno_input_sexPerm.fam
TrioInput[7]=SSCDuo_Trio_pheno_input_sexPerm.fam
TrioInput[8]=SSCOmni_Trio_pheno_input_sexPerm.fam
####R script:
args <- commandArgs(trailingOnly = TRUE)
TrioInput <- args[1]   ### a list of probands and sexes: FamID, IndID, Sex, "2"
Trio_inds <- read.table(file=paste("/path/",TrioInput,sep=""),header=F)
colnames(Trio_inds) = c("FID", "IID", "PID", "MID", "SEX", "PHENO")
Trio_Case <- Trio_inds[Trio_inds$PHENO == 2,]
Trio_Case_n <- nrow(Trio_Case)
Trio_fem_Case <- Trio_Case[Trio_Case$SEX == 2,]
Trio_fem_Case_n <- nrow(Trio_fem_Case)
####ensure correct pheno being written out 
Perm = 100
for ( i in 1:Perm )	{
	Trio_fem_casePerm_n <- sample(Trio_Case_n,Trio_fem_Case_n)
	Trio_fem_Perm <-Trio_Case[Trio_fem_casePerm_n,c("FID","IID","PHENO")]
	Trio_male_Perm <-Trio_Case[-Trio_fem_casePerm_n,c("FID","IID","PHENO")]
	write.table(Trio_fem_Perm, file=paste("/path/",TrioInput,"_fem_perm_",i,sep=""),col.names=F,row.names=F,quote=F)
	write.table(Trio_male_Perm, file=paste("/path/",TrioInput,"_male_perm_",i,sep=""),col.names=F,row.names=F,quote=F)
}


###Create sex-perms Case Control:
###Input R args from shell script:
CCinput[1]=AGP_CC_pheno_input_sexPerm.fam
CCinput[2]=AGRE_CC_pheno_input_sexPerm.fam
CCinput[3]=CHARGE_CC_pheno_input_sexPerm.fam
CCinput[4]=CHOP_CC_pheno_input_sexPerm.fam
CCinput[5]=EMA_CC_pheno_input_sexPerm.fam
CCinput[6]=SEED_CC_pheno_input_sexPerm.fam
CCinput[7]=SS_groups_CC_pheno_input_sexPerm.fam
CCinput[8]=SSC1M_CC_pheno_input_sexPerm.fam
CCinput[9]=SSCDuo_CC_pheno_input_sexPerm.fam
CCinput[10]=SSCOmni_CC_pheno_input_sexPerm.fam
####R script:
args <- commandArgs(trailingOnly = TRUE)
CCInput <- args[1]      ### .fam file with FamID, IndID, FatID, MotID, Sex & affection
CC_inds <- read.table(file=paste("/path/",CCInput, sep=""),header=F)
colnames(CC_inds) = c("FID", "IID", "PID", "MID", "SEX", "PHENO")
CC_Control <- CC_inds[CC_inds$PHENO==1,]
CC_fem_Control_n <- nrow(CC_Control[CC_Control$SEX==2,])
CC_Case <- CC_inds[CC_inds$PHENO==2,]
CC_fem_Case_n <- nrow(CC_Case[CC_Case$SEX==2,])
CC_Case_n <- nrow(CC_Case)
CC_Control_n <- nrow(CC_Control)
####ensure correct pheno being written out for CC
Perm = 100
for ( i in 1:Perm )	{
	CC_fem_casePerm_n <- sample(CC_Case_n,CC_fem_Case_n)
	CC_fem_ctrlPerm_n <- sample(CC_Control_n,CC_fem_Control_n)
	CC_fem_casePerm <- CC_Case[CC_fem_casePerm_n,c("FID","IID","PHENO")]
	CC_fem_ctrlPerm <- CC_Control[CC_fem_ctrlPerm_n,c("FID","IID","PHENO")]
	CC_male_casePerm <- CC_Case[-CC_fem_casePerm_n,c("FID","IID","PHENO")]
	CC_male_ctrlPerm <- CC_Control[-CC_fem_ctrlPerm_n,c("FID","IID","PHENO")]
	CC_fem_perm <- merge(CC_fem_casePerm,CC_fem_ctrlPerm,all=TRUE)
	CC_male_perm <- merge(CC_male_casePerm,CC_male_ctrlPerm,all=TRUE)
	write.table(CC_fem_perm,file=paste("/path/",CCInput, "_fem_perm_",i,sep=""),col.names=F,row.names=F,quote=F)
	write.table(CC_male_perm,file=paste("/path/",CCInput,"_male_perm_",i,sep=""),col.names=F,row.names=F,quote=F)
}
