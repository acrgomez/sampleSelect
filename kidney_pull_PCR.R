require(dplyr)
require(lubridate)
options(stringsAsFactors = F)


savedSerum <- read.csv("~/Dropbox/Sample Selection/new/setAside.csv")
allSamples <- read.csv("~/Dropbox/Sample Selection/new/allSamples.csv")
bloodwork <- read.csv("~/Dropbox/Sample Selection/new/bloodwork.csv",na.strings=c("","NA"," "))

# recentSamples <- read.csv("~/Google Drive/allSamples_dec2014todec2016.csv")
# 
# 
# #savedSerum$ID <- paste("CSL-",savedSerum$ID,sep='')
# 
# 
# 
# 
# ###Just to make sure that we don't carry over a value to an empty cell
# ###that is actually supposed to be empty
# # colnames(master) <- c("Admit_date","Disposition","MMC_ID","Name","bag","box",
# #                       "Collection_method","Sample_Date",
# #                       "Freezer","Project",
# #                       "Recipient","Result","Sent","Tissue","Media","Volume")
# 
# #write.csv(master,"~/Google Drive/allSamples_dec2014todec2016.csv",row.names=F)
# 
# #####populate and convert dates#####
# master <- allSamples
# 
# 
# 
# filledIndexes <- allSamples$MMC_ID != ''
# for (col in seq_along(colnames(allSamples))){
#   master[,col][filledIndexes&as.character(master[,col])==''] = NA             
# }
# 
# ###blank cells on filemaker, and whose value is only implied
# populate <- function(list, indexlist){
#   len = length(list)
#   ii=1
#   while(ii < len){
#     
#     if(!indexlist[ii+1]) list[ii+1] <- list[ii]
#     ii=ii+1
#   }
#   return(list)
# }
# 
# 
# 
# 
# master$Admit_date <- as.Date(populate(master$Admit_date,filledIndexes), "%m/%d/%Y")
# master$Age_Class <- populate(master$Age_Class,filledIndexes)
# master$Disposition <- populate(master$Disposition,filledIndexes)
# master$Project <- populate(master$Project, filledIndexes)
# master$MMC_ID <- populate(master$MMC_ID,filledIndexes)
# master$Name <- populate(master$Name,filledIndexes)
# master$Sample_Date <- as.Date(master$Sample_Date, "%m/%d/%Y")
# master$DaysSinceAdmit <- master$Sample_Date - master$Admit_date
# master <- arrange(master, MMC_ID, Admit_date)
# master$LastAtTMMC <- NA  ###placeholder
# 
# 
# write.csv(master,"~/Google Drive/prepared_allSamples_dec2014todec2016.csv",row.names=F)


recentSamples <- read.csv("~/Google Drive/prepared_allSamples_dec2014todec2016.csv")
recentSamples[is.na(recentSamples)] = ''



ids <- unique(c(master$MMC_ID,savedSerum$ID))

haskidney <- filter(master,Tissue=="kidney")
haskidneyids <- unique(haskidney$MMC_ID)

hasserum <- unique(c(savedSerum$ID, filter(master,Tissue=="serum")$MMC_ID))
hasurine <- unique(c(filter(master,Tissue %in% c("urine pellet","urine"))$MMC_ID))

hasall <- filter(haskidney, MMC_ID %in% hasserum & MMC_ID %in% hasurine)

unsent_hasall = filter(hasall, Sent=='')

sentkidneyPCR <- filter(haskidney,Sent !='',Recipient=="Hollings Marine Lab")

tosend_hasall <- filter(unsent_hasall,!(unsent_hasall$MMC_ID %in% sentkidneyPCR$MMC_ID),
                        Media =="whirlpak",Recipient != "Discarded")
tosend_hasall$LastAtTMMC<- NULL


write.csv(tosend_hasall,"~/Google Drive/kidneys_tosend.csv",row.names = F)



#####recent urine


urine <- filter(recentSamples, Tissue %in% c("urine","urine pellet"),!(Collection_method %in% 
                  c("pen floor","free catch","culture media")))

sent_pcr_urine <- filter(urine, Recipient =="Hollings Marine Lab", Sent !="")

notsenturine <- filter(urine, !(MMC_ID %in% sent_pcr_urine$MMC_ID))

notsenturine <- arrange(notsenturine, MMC_ID, DaysSinceAdmit)
uninotsenturine <- filter(notsenturine,!duplicated(MMC_ID))

serum <- filter(recentSamples, Tissue =="serum", 
                Recipient %in% c("CDC",""," ","CAHFS"))
temp <- filter(serum,duplicated(MMC_ID))
serum$last <- !(serum$MMC_ID %in% temp$MMC_ID) & serum$Recipient == ""

temp2 <- filter(serum, last)
urinewithserum <- filter(uninotsenturine, MMC_ID %in% serum$MMC_ID)
recent_urine_tosent <- filter(urinewithserum, !(MMC_ID %in% temp2$MMC_ID))


write.csv(recent_urine_tosent,"~/Google Drive/recenturine_tosent.csv",row.names = F)





