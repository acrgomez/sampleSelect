setwd("~/Dropbox/Sample Selection/April19/")
require(dplyr)
require(lubridate)
options(stringsAsFactors = FALSE)


emily <- read.csv("Emily_Whitmer_CSL Urine samples for PCR 2.11.19_KP.updated.csv")

correctedcols <- read.csv("feedback/Copy of Rehab CSL serum MAT Results_CorrectedDate.csv")
correctedcols$DateSampledCorrected <- as.Date(correctedcols$DateSampledCorrected, "%m/%d/%y") 


merg <- read.csv("~/Google Drive/CSL_Leptospirosis/Current_Data/Merged_Strand_WCSL_Chem_MAT_PCR_Culture_Master.csv")
merg$StrandDate <- as.Date(merg$StrandDate)
#merg <- filter(merg,year(StrandDate)>2013)
bloodData <- read.csv("feedback/2017-2019 Bloodwork Data 5-7-19.csv")
bloodData$Specimens..Sample_Date <-as.Date(bloodData$Specimens..Sample_Date,"%m/%d/%Y")
bloodData$pasty <- paste(bloodData$Animal..TMMC_ID,bloodData$Specimens..Sample_Date)

strandData <- read.csv("2016-2019 CSL Strand Data 4-1-19.csv")
tmeData <- read.csv("2016-2019 CSL TME Data 4-1-19.csv")
specimen <- read.csv("2016-2019 Lepto Specimen Data 4-1-19.csv") %>% filter(aux_Species..Common_Name=="California Sea Lion")
specimen$Index <- seq_along(specimen$Animal..TMMC_ID)
specimen$AdmitDate <- strandData$Admit_Date[match(specimen$Animal..TMMC_ID,strandData$Animal..TMMC_ID)] %>% as.Date("%m/%d/%Y")
specimen$Sample_Date <- as.Date(specimen$Sample_Date,"%m/%d/%Y")
specimen$pasty <- paste(specimen$Animal..TMMC_ID,specimen$Sample_Date)
specimen$InEmilyList <- specimen$Animal..TMMC_ID %in% emily$ID

specimenav <- filter(specimen, specimen$Specimens_Box..Name!="",
                     specimen$Collection_Method !="unable to collect",
                     specimen$Storage_Method !="+30 C (incubator)",
                     specimen$Storage_Method != "culture media")
specimenav$SinceAdmit <- specimenav$Sample_Date - specimenav$AdmitDate



merg$SampleDateChem <- as.Date(merg$SampleDateChem)
merg$chemPasty <- paste(merg$ID,merg$SampleDateChem)
mergchem <- select(merg,BUN,CA,CL,CREAT,PHOS,K,NA.) %>%
            rowSums()
pastymergchem <- merg$chemPasty[!is.na(mergchem)]

chemcols <- select(bloodData,BUN,Calcium,Chloride,Creatinine,Phosphorus,Potassium,Sodium) %>%
            rowSums()
hasChem <- filter(bloodData,!is.na(chemcols))
allchempasty <- c(hasChem$pasty,pastymergchem)


specimen$hasChem <- specimen$pasty %in% allchempasty

unavailable <- filter(specimen, specimen$Specimens_Box..Name=="",
                      specimen$Collection_Method =="unable to collect",
                      specimen$Storage_Method =="+30 C (incubator)")



serum <- filter(specimenav,Tissue == "serum")
#Do any serum have MAT already?
merg$pasty <- paste(merg$ID,merg$SampleDateMAT)
serum <- filter(serum,!(pasty%in%merg$pasty)) ###removes 13 samples
serum <- serum %>% filter(!duplicated(pasty)) %>% arrange(Animal..TMMC_ID,Sample_Date)
wanted <- !duplicated(serum$Animal..TMMC_ID)
serum$wanted <- "Maybe"
serum$wanted[wanted] = "Yes"


#Does any kid, pel, or uri have same day PCR already
pastypcr <- c(paste(merg$ID,merg$SampleDatePCR.kidney),
              paste(merg$ID,merg$SampleDatePCR.urine),
              paste(merg$ID,merg$SampleDatePCR.pellet),
              paste(merg$ID,merg$SampleDatePCR.all),
              paste(merg$ID,merg$SampleDatePCR.first),
              paste(merg$ID,merg$SampleDatePCR.last))

specimenav$SameDayPCRresult <- specimenav$pasty %in% pastypcr


kidneys <- filter(specimenav,Tissue == "kidney",Storage_Method %in% c("whirlpak","cryovial"))
kidneys$priority <-3
pellets <- filter(specimenav,Tissue == "urine pellet")
pellets$priority <-1
urine <- filter(specimenav,Tissue %in% c("urine"))
urine$priority <- 2

allPCR <- rbind(kidneys,pellets,urine)

allPCR <- arrange(allPCR,allPCR$Animal..TMMC_ID,allPCR$Sample_Date,priority)

pcrsamples <- filter(allPCR,!duplicated(pasty),!SameDayPCRresult)

pcrsamples$SameDaySerumMAT <- pcrsamples$pasty %in% c(serum$pasty,merg$pasty)


#chem within 2 days
pcrsamples$ChemWithin2days <- pcrsamples$pasty %in% allchempasty |
                        paste(pcrsamples$Animal..TMMC_ID,pcrsamples$Sample_Date-2)   %in% allchempasty |
                        paste(pcrsamples$Animal..TMMC_ID,pcrsamples$Sample_Date-1)   %in% allchempasty |
                        paste(pcrsamples$Animal..TMMC_ID,pcrsamples$Sample_Date+1)   %in% allchempasty |
                        paste(pcrsamples$Animal..TMMC_ID,pcrsamples$Sample_Date+2)   %in% allchempasty

hasMAT <- filter(merg,!is.na(LogPom))
hasMAT$SampleDateMAT <- as.Date(hasMAT$SampleDateMAT)
hasMATmini <- select(hasMAT,ID,SampleDateMAT)

corrmini <- data.frame(ID=correctedcols$TMMC.ID,SampleDateMAT=correctedcols$DateSampledCorrected)
hasmatplus <- rbind(hasMATmini,corrmini)
hasmatplus$pasty <- paste(hasmatplus$ID,hasmatplus$SampleDateMAT)

pcrsamples$MATwithin21 <- NA
pcrsamples$Samplewithin21 <- NA

for(ii in seq_along(pcrsamples$Animal..TMMC_ID)){
    sample <- pcrsamples[ii,]
    matdone <- filter(hasmatplus,ID==sample$Animal..TMMC_ID)
    
    if(dim(matdone)[1]>0){
        difs <- abs(matdone$SampleDateMAT - sample$Sample_Date)
        if(any(difs<=21)){
            pcrsamples$MATwithin21[ii] <- TRUE
        } else {
            pcrsamples$MATwithin21[ii] <- FALSE
        }
    } else {
        pcrsamples$MATwithin21[ii] <- FALSE
    }
    
    seravai <- filter(serum,Animal..TMMC_ID==sample$Animal..TMMC_ID)
    if(dim(seravai)[1]>0){
        difs <- abs(seravai$Sample_Date - sample$Sample_Date)
        if(any(difs<=21)){
            pcrsamples$Samplewithin21[ii] <- TRUE
        } else {
            pcrsamples$Samplewithin21[ii] <- FALSE
        }
    } else {
        pcrsamples$Samplewithin21[ii] <- FALSE
    }
}



wantedpcrsamples <- filter(pcrsamples,Samplewithin21 | MATwithin21 | SameDaySerumMAT | InEmilyList | ChemWithin2days)
dim(pcrsamples)
dim(wantedpcrsamples)
#####



serum <- filter(serum,!(wanted=="Maybe" & !(Animal..TMMC_ID %in% wantedpcrsamples$Animal..TMMC_ID)))
serum$SameDayPCRsampleAvailable <- serum$pasty %in% wantedpcrsamples$pasty
serum$diffDayPCRsampleAvailable <- (serum$Animal..TMMC_ID %in% wantedpcrsamples$Animal..TMMC_ID) & !serum$SameDayPCRsampleAvailable


negPCRlist <- filter(pcrsamples,!(Index %in% wantedpcrsamples$Index))


write.csv(serum,"~/Desktop/serumSamples.csv",row.names = F)
write.csv(wantedpcrsamples,"~/Desktop/PCRsamples.csv",row.names = F)
write.csv(negPCRlist,"~/Desktop/UnwantedPCRsamples.csv",row.names=F)


##############

serum <- read.csv("~/Desktop/serumSamples.csv")
serum$Sample_Date <- as.Date(serum$Sample_Date)

serum$DaysToClosest <- NA
for(ii in seq_along(serum$Animal..TMMC_ID)){
    sample <- serum[ii,]
    all <- filter(serum[-ii,], Animal..TMMC_ID == sample$Animal..TMMC_ID)
    if (dim(all)[1] >= 1){
        diff <- all$Sample_Date - sample$Sample_Date
        serum$DaysToClosest[ii] <-diff[which.min(abs(diff))]
    }
}

write.csv(serum,"~/Desktop/serumSamples_newcol.csv",row.names = F)



minipcrres <- select(merg,ID,SampleDatePCR.kidney,SampleDatePCR.urine,
                     SampleDatePCR.pellet,SampleDatePCR.all,SampleDatePCR.first,
                     SampleDatePCR.last)

pcralldates <- melt(minipcrres,id.vars = 1)
colnames(pcralldates)=c("ID","var","Date")
pcralldates$Date <- as.Date(pcralldates$Date)

haspcr <- filter(merg, !is.na(PCR.ever))
chemis <- as.data.frame(matrix(unlist(strsplit(allchempasty," ")),ncol=2,byrow=T))
colnames(chemis)=c("ID","Date")
chemis$Date <- as.Date(chemis$Date)

pcr <- read.csv("~/Desktop/PCRsamples_all.csv")
pcr$Sample_Date <- as.Date(pcr$Sample_Date,"%m/%d/%Y")


pcr$DaysToClosestSerum <- NA
pcr$DaysToClosestMAT <- NA
pcr$DaysToClosestPCR <- NA
pcr$DaysToClosestChem <- NA
pcr$DaysToClosestPCRsample <- NA


for (ii in seq_along(pcr$Animal..TMMC_ID)){
    samplepcr <- pcr[ii,]
    serums <- filter(serum, Animal..TMMC_ID == samplepcr$Animal..TMMC_ID)
    chems <- filter(chemis, ID == samplepcr$Animal..TMMC_ID)
    mats <- filter(hasmatplus,ID == samplepcr$Animal..TMMC_ID)
    pcrs <- filter(pcralldates,ID == samplepcr$Animal..TMMC_ID)
    

    pcr$DaysToClosestSerum[ii] <- min(abs(serums$Sample_Date - samplepcr$Sample_Date),na.rm=T)
    pcr$DaysToClosestMAT[ii] <- min(abs(mats$SampleDateMAT - samplepcr$Sample_Date),na.rm=T)
    pcr$DaysToClosestPCR[ii] <- min(abs(pcrs$Date - samplepcr$Sample_Date),na.rm=T)
    pcr$DaysToClosestChem[ii] <- min(abs(chems$Date - samplepcr$Sample_Date),na.rm=T)
    
    otherpcrs <- filter(pcr[-ii,],Animal..TMMC_ID==samplepcr$Animal..TMMC_ID)
    if(dim(otherpcrs)[1]>=1){
        diff <- samplepcr$Sample_Date - otherpcrs$Sample_Date
        pcr$DaysToClosestPCRsample[ii] <- diff[which.min(abs(diff))]
    }

}




write.csv(pcr,"~/Desktop/finalPCRlist.csv",row.names=F)

























filter(serum,Animal..TMMC_ID == "CSL-13888")

filter(bloodData,Animal..TMMC_ID == "CSL-13836")

filter(merg,ID=="CSL-8079")


filter(bloodData,Animal..TMMC_ID == "CSL-10073")


require(ggplot2)

mini$SinceAdmit <- as.numeric(as.Date(mini$SampleDateChem)-as.Date(mini$AdmitDate))
mini2 <- filter(mini,SinceAdmit <30)
ggplot(mini2)+
    geom_point(aes(x=BUN,y=CREAT,color=LeptoCode),pch=3)+
    ylim(c(0,30))+
    xlim(c(0,600))+
    #scale_color_viridis()+
    theme_classic()
