setwd("~/Dropbox/Sample Selection/")

###This code is for a specific sample pull request by Katie Prager
require(dplyr)
require(lubridate)
options(stringsAsFactors = FALSE)

urineSent <- read.csv("samples_to_WU_April_2015_to_CF_4.30.15_sent 06-10-2015_kp_update_7.10.15.csv")
serumSent <- read.csv("~/Dropbox/Sample Selection/CDC_SHIPMENT_SERUM2015.csv")
savedSerum <- read.csv("Admit_Serum_Set_Aside_KP_to_Jan2016.csv")

allMAT <- read.csv("ALL_CSL_MAT_13_May_2016.csv")
allSamples <- read.csv("CSL_samples_to_05.20.2016.csv")
bloodwork <-read.csv("CSL_FM_bloodwork_new_fromKatie.csv",na.strings=c("","NA"," "))

pks <- read.csv("all2016PKSAnimals.csv")


###function that adds repeated entries to columns that have 
###blank cells on filemaker, and whose value is only implied
populate <- function(list, indexlist){
    len = length(list)
    ii=1
    while(ii < len){
        
        if(!indexlist[ii+1]) list[ii+1] <- list[ii]
        ii=ii+1
    }
    return(list)
}

###Just to make sure that we don't carry over a value to an empty cell
###that is actually supposed to be empty
master <- allSamples



filledIndexes <- allSamples$MMC_ID != ''
for (col in seq_along(colnames(allSamples))){
    master[,col][filledIndexes&as.character(master[,col])==''] = NA             
}



#####populate and convert dates#####

master$Admit_Date <- as.Date(populate(master$Admit_Date,filledIndexes), "%m/%d/%Y")
master$Age_Class <- populate(master$Age_Class,filledIndexes)
master$Disposition <- populate(master$Disposition,filledIndexes)
master$Studies..Study <- populate(master$Studies..Study, filledIndexes)
master$Studies..Study[is.na(master$Studies..Study)] = ''
master$Disposition_Date <- populate(master$Disposition_Date,filledIndexes)
master$MMC_ID <- populate(master$MMC_ID,filledIndexes)
master$Name <- populate(master$Name,filledIndexes)
master$Sex <- populate(master$Sex,filledIndexes)
master$Strand_Date <- as.Date(populate(master$Strand_Date,filledIndexes),
                              "%m/%d/%Y")
master$Post_Mortem_Date = as.Date(populate(master$Post_Mortem_Date,
                                           filledIndexes), "%m/%d/%Y")
master$Sample_Date <- as.Date(master$Specimens..date_sampled, "%m/%d/%Y")

master$PKS <- (master$MMC_ID %in% pks$MMC_ID) | (master$Studies..Study=="Prospective Kidney Study")
master$DaysSinceAdmit <- master$Sample_Date - master$Admit_Date
master$DaysToPostMortem <- master$Post_Mortem_Date - master$Sample_Date
master$Tissue <- master$Specimens..specimen_tissue
master$Cancer <- master$Studies..Study == "Cancer case control"

cancer <- unique(master$MMC_ID[master$Cancer])


master <- arrange(master, MMC_ID, Admit_Date)
master$LastAtTMMC <- NA  ###placeholder

master= filter(master,Name!="Godsealla")
master=filter(master,Name!="Ellye")

serumSent$DateSampled= as.Date(serumSent$DateSampled, "%m/%d/%Y")
urineSent$DateSampled= as.Date(urineSent$DateSampled, "%m/%d/%Y")

#####################################################
#end of setup
#####################################################


#####################################################



###############################################################################
#urine work##
###############

lapply(strsplit(as.character(test),split=" "),FUN=unlist)

master$Tissue <- tolower(master$Specimens..specimen_tissue)


master$urine <- grepl(x=master$Tissue,pattern = "urine")

urinemaster0 <- filter(master, urine)



strsplit(master$Specimens..specimen_tissue,split = c(" "))

urinemaster2 = subset(master, master$Specimens..specimen_tissue 
                      %in% c('urine','urine pellet') & 
                          master$Specimens..collection_method !='pen floor')

urinemaster1 = filter(urinemaster2,Specimens..freezer!='')

urinemaster  =  arrange(urinemaster1,MMC_ID,DaysSinceAdmit,
                        desc(Specimens..specimen_tissue),
                        grepl(pattern = "Prager",Specimens..box,ignore.case = T))



besturine = filter(urinemaster,!duplicated(MMC_ID))
senturine2 <- filter(urinemaster2,Specimens..freezer ==''&
                         Specimens..project =="Leptospira PCR")

notsentbest=filter(besturine,!(MMC_ID%in%senturine2$MMC_ID))



pksurine = filter(urinemaster1,MMC_ID %in% pks$MMC_ID)
pksurine = arrange(pksurine,MMC_ID,abs(DaysToPostMortem),
                   desc(Specimens..specimen_tissue))
postmpks <- filter(pksurine,!duplicated(MMC_ID))
postmpks <- filter(postmpks,!(interaction(MMC_ID,Sample_Date) %in%
                                  interaction(besturine$MMC_ID,besturine$Sample_Date)))

finalurine <- rbind(notsentbest,postmpks)
oldurine = filter(finalurine,DaysSinceAdmit >14)




ffinalurine <- read.csv("~/Dropbox/final/urineforpdf.csv")
ffinalurine$Date <- as.Date(ffinalurine$Date, "%m/%d/%Y")
ffinalurine$Date[117:119] = ffinalurine$Date[117:119]+years(2000)
index <- match(interaction(ffinalurine$MMC_ID,ffinalurine$Date),
               interaction(finalurine$MMC_ID,finalurine$Sample_Date))


ffinalurine$Collection_method <- finalurine$Specimens..collection_method[index]

#write.csv(ffinalurine,"~/Dropbox/final/urineforpdf.csv",row.names = F)


###############################################################################
#serum work##
###############


serummaster2 = subset(master, master$Specimens..specimen_tissue  == 'serum')

serummaster1 = subset(serummaster2, grepl("CSL", serummaster2$Specimens..box))    #available serum is in a CSL box
serummaster1$LastAtTMMC <- !(duplicated(serummaster1$MMC_ID) |
                                 duplicated(serummaster1$MMC_ID,fromLast = T))    #last at tmmc means not duplicated, but duplicated leaves first occurance


serummaster <- arrange(serummaster1,MMC_ID,DaysSinceAdmit,Specimens..freezer,
                       Specimens..box)

wantedsera1 = filter(serummaster, !(MMC_ID %in% c(savedSerum$ID,
                                                  allMAT$MMCID,
                                                  serumSent$MMC_ID)))

wantedsera2 = filter(wantedsera1, !LastAtTMMC | PKS)
wantedsera2 = arrange(wantedsera2,MMC_ID,DaysSinceAdmit,Specimens..freezer,
                      Specimens..box)
wantedsera = filter(wantedsera2, !duplicated(MMC_ID))

wantedboyce = filter(serummaster2, grepl("Boyce",Specimens..recipient),
                     !(MMC_ID %in% c(savedSerum$ID,
                                     allMAT$MMCID,
                                     serumSent$MMC_ID,
                                     wantedsera$MMC_ID)))
wantedboyce = arrange(wantedboyce,DaysSinceAdmit)
wantedboyce = filter(wantedboyce, !duplicated(MMC_ID))
wantedboyce$Boyce="Yes"



oldurineSera = filter(wantedsera2, MMC_ID %in% oldurine$MMC_ID)
index = match(oldurineSera$MMC_ID,oldurine$MMC_ID)
oldurineSera$diff=abs(oldurineSera$Sample_Date-oldurine$Sample_Date[index])
oldurineSera <- filter(oldurineSera,diff<=5)
oldurineSera <- arrange(oldurineSera,MMC_ID,diff,Specimens..freezer)
oldurineSera = filter(oldurineSera, !duplicated(MMC_ID))
oldurineSera$diff=NULL

pkspostmsera <- filter(serummaster1,PKS,abs(DaysSinceAdmit) > abs(DaysToPostMortem))

pkspostmsera <- arrange(pkspostmsera,MMC_ID,DaysToPostMortem)

pkspostmsera <- filter(pkspostmsera,!duplicated(MMC_ID))

pkspostmsera <- filter(pkspostmsera,!(interaction(MMC_ID,Specimens..bag) %in%
                                          interaction(MMC_ID,Specimens..bag)))

finalSera <- rbind(wantedsera,pkspostmsera,oldurineSera)

finalSera$Boyce=''

finalSera <- rbind(finalSera,wantedboyce)


###############################################
###pups


pups <- filter(finalSera,Age_Class=="Pup", !(MMC_ID %in% c(urineSent$MMC_ID,finalurine$MMC_ID,senturine2$MMC_ID)))


pups$cummulativeMon <- month(pups$Admit_Date) + (year(pups$Admit_Date) - 
                        range(year(pups$Admit_Date))[1])*12
boys <- filter(pups, Sex=="1. Male")
girls <- filter(pups, Sex=="2. Female")
selected <-data.frame()
for (mon in unique(boys$cummulativeMon)){
    temp <- filter(boys, cummulativeMon == mon)
    if (length(temp[,1]) > 5)
        temp <- sample_n(temp,5)
    selected <- rbind(selected,temp)
}

for (mon in unique(girls$cummulativeMon)){
    temp <- filter(girls, cummulativeMon == mon)
    if (length(temp[,1]) > 5)
        temp <- sample_n(temp,5)
    selected <- rbind(selected,temp)
}

finalseranopup <- filter(finalSera,!(interaction(MMC_ID,Admit_Date) %in% 
                             interaction(pups$MMC_ID, pups$Admit_Date)))
selected$cummulativeMon<- NULL

finalsera <- rbind(finalseranopup,selected)

#####################################################
#####no urine without serum

finalurine <- filter(finalurine,(MMC_ID %in% c(finalSera$MMC_ID,serumSent$MMC_ID,allMAT$MMCID,savedSerum$ID)))

#temp=filter(finalurine,!(MMC_ID %in% c(finalSera$MMC_ID,serumSent$MMC_ID,allMAT$MMCID,savedSerum$ID)))
#filter(serummaster,MMC_ID %in% temp$MMC_ID)
#unique(temp$MMC_ID)
#####################################################



#tmmcsera <- filter(master, Tissue == "serum",!is.na(Specimens..bag))
# 
# savedSerum$ID[savedSerum$ID=="WCSL"]="WCSL-117-15"
# savedSerum <- filter(savedSerum,To == '')
# savedSerum$Cancer <- savedSerum$ID %in% cancer
# savedSerum$NoOtherAtTMMC <- !(savedSerum$ID %in% tmmcsera$MMC_ID)
# savedSerum$Sample_Date <- as.Date(savedSerum$Date.Drawn,"%m/%d/%Y")
#Saved serum already sent?
# 
# filter(savedSerum,ID %in% serumSent$MMC_ID)   #lakzmi
# savedSerum <- filter(savedSerum,!(ID %in% allMAT$MMCID))
# write.csv(savedSerum,"savedSerumFinal.csv",row.names=F)
######################Chem#####################################################
bloodwork$MMC_ID <- bloodwork$Stranding..MMC_ID

chem <- filter(bloodwork,Sample_Type=="Chemistry")

chem$Sample_Date <- as.Date(chem$Sampled_TS,"%m/%d/%Y")

chemWant <- select(chem,MMC_ID,Sample_Date,BUN,Calcium,Chloride,Creatinine,     #all we need done
                   Phosphorus,Potassium,Sodium)

for (ii in seq_along(chemWant$MMC_ID)){                                         #this is way more complicated than it needs to be, but a pretty flag for chem
    if (sum(is.na(chemWant[ii,]))){
        chemWant$Flag[ii] <- paste(c("Needs:", 
                                     colnames(chemWant)[is.na(chemWant[ii,])]), 
                                   collapse=' ')
    } else {chemWant$Flag[ii] = "Done"}
}


chemWant$Flag <-'' 

chemWant <- arrange(chemWant,MMC_ID,Flag)

# savedSerum$chemAnyDate <- chemWant$Flag[match(savedSerum$ID,chemWant$MMC_ID)]
# savedSerum$chemSameDate <- chemWant$Flag[match(
#                             interaction(savedSerum$ID,savedSerum$Sample_Date),
#                             interaction(chemWant$MMC_ID,chemWant$Sample_Date))]
# savedSerum$chemSameDate[is.na(savedSerum$chemSameDate)] ="NO"

##############################################################################################
##############################################################################################





finalsera$chemAnyDate <- chemWant$Flag[match(finalsera$MMC_ID,chemWant$MMC_ID)]
finalsera$chemSameDate <- chemWant$Flag[match(
    interaction(finalsera$MMC_ID,finalsera$Sample_Date),
    interaction(chemWant$MMC_ID,chemWant$Sample_Date))]
finalsera$chemSameDate[is.na(finalsera$chemSameDate)] ="NO"
finalsera$chemAnyDate[is.na(finalsera$chemAnyDate)] ="NO"


##############################################################################################
##############################################################################################


finalseraclean <- select(finalsera,MMC_ID,Name,Sample_Date,Specimens..collection_method,Specimens..bag,
                         Specimens..box,Specimens..freezer,PKS,Cancer,LastAtTMMC,Boyce,chemSameDate,chemAnyDate)
finalseraclean=arrange(finalseraclean,Cancer,Boyce,PKS,Specimens..freezer,Specimens..box,Specimens..bag)

finalurineclean <- select(finalurine,MMC_ID,Name,Sample_Date,Specimens..bag,
                          Specimens..box,Specimens..freezer,Tissue,Specimens..tissue_quantity,PKS,Cancer)
finalurineclean=arrange(finalurineclean,Cancer,PKS,Specimens..freezer,Specimens..box,Specimens..bag)

write.csv(finalseraclean,"Selected_Serum_Clean.csv",row.names=F)
write.csv(finalurineclean,"Selected_Urine_Clean.csv",row.names=F)


##############################################################################################
##############################################################################################
##############################################################################################
#Extra verification

temp=subset(finalurine, !(finalurine$MMC_ID %in% finalsera$MMC_ID))
filter(temp,MMC_ID %in% allMAT$MMCID)
filter(finalurine,MMC_ID %in% savedSerum$ID)


subset(finalsera,LastAtTMMC)



##############################################################################################
##############################################################################################
##############################################################################################

finalserum <- read.csv("~/Dropbox/final/serum.csv")

filter(finalserum,MMC_ID %in% allMAT$MMCID)
filter(finalserum,MMC_ID %in% serumSent$MMC_ID)


finalpks = filter(finalserum,MMC_ID %in% pks$MMC_ID & MMC_ID %in% serumSent$MMC_ID)


finalnosent <- filter(finalserum, !(MMC_ID %in% serumSent$MMC_ID))
                      
                      
                      
filter(finalnosent,MMC_ID %in% allMAT$MMCID)
filter(finalnosent,MMC_ID %in% serumSent$MMC_ID)                  
        
write.csv(finalnosent,"~/Dropbox/final/Serum_revised_2016.csv",row.names = F)
                      
                      
                      
                      
                      
                      
                      

