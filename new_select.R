require(dplyr)
require(lubridate)
options(stringsAsFactors = F)

master1 <- read.csv("~/Google Drive/1970to2017_niceandclean.csv")

savedSerum <- read.csv("~/Dropbox/Sample Selection/new/setAside.csv")
bloodwork <- read.csv("~/Dropbox/Sample Selection/new/bloodwork.csv",na.strings=c("","NA"," "))
matresults <- read.csv("~/Google Drive/CSL_Leptospirosis/Current_Data/MAT_Master.csv")

#savedSerum$ID <- paste("CSL-",savedSerum$ID,sep='')
#write.csv(savedSerum,"~/Dropbox/Sample Selection/new/setAside.csv",row.names=F)

master <- select(master1,Admit_Date,
                 MMC_ID,Name,Sex,Age_Class,Disposition,Specimens..amount,
                 Specimens..collection_method,Specimens..date_sampled,
                 Specimens..specimen_tissue,Specimens..storage,Specimens..tissue_quantity,Studies..Date_Sent,
                 Specimens..freezer,Specimens..box,Specimens..bag,Studies..Study,
                 Specimens..sent,Specimens..project,Specimens..recipient,
                 Specimens..results)

master <- filter(master,!is.na(MMC_ID)&!is.na(Specimens..date_sampled)&!is.na(Admit_Date))

master$ID <- sub(x=master$MMC_ID,pattern = "-R1",replacement = "")
master$ID <- sub(x=master$ID,pattern = "-F",replacement = "")

#####populate and convert dates#####

master$Admit_Date <- as.Date(master$Admit_Date)
master$Studies..Study[is.na(master$Studies..Study)] = ''
master$Specimens..date_sampled <- as.Date(master$Specimens..date_sampled,"%m/%d/%Y")
master$DaysSinceAdmit <- master$Specimens..date_sampled - master$Admit_Date
master <- arrange(master, MMC_ID, Admit_Date)
master$PKSflagged <- agrepl(x = master$Studies..Study,pattern = "Prospective Kidney Study",ignore.case = T)

####

sera <- filter(master,agrepl(x=master$Specimens..specimen_tissue,pattern ="serum"))
kidneys <- filter(master,agrepl(x=master$Specimens..specimen_tissue,pattern ="kidney"))
urines <- filter(master,agrepl(x=master$Specimens..specimen_tissue,pattern ="urine"))
supernatants <- filter(master,agrepl(x=master$Specimens..specimen_tissue,pattern ="supernatant"))


###
sera$hasresults <- sera$ID %in% matresults$ID
sera$setaside <- sera$ID %in% savedSerum$ID
sera$atTMMC <- grepl("-80",sera$Specimens..freezer)
sera$restrand <- agrepl("restrand",sera$Name)
tmmcsera <- filter(sera,atTMMC)
tmmcsera$last <- !duplicated(tmmcsera$ID) & !duplicated(tmmcsera$ID,fromLast = T)
nottmmcsera <- filter(sera,!atTMMC)
nottmmcsera$last <- NA
sera <- rbind(tmmcsera,nottmmcsera)


wantsera <- filter(sera,atTMMC & !hasresults & !setaside & !last & Admit_Date > "2016-01-01")




  
haskidneyids <- unique(haskidney$MMC_ID)


ids <- unique(c(master$MMC_ID,savedSerum$ID))

haskidney <- filter(master,agrepl(x=master$Tissue,pattern ="kidney"))
haskidneyids <- unique(haskidney$MMC_ID)

hasserum <- unique(c(savedSerum$ID, filter(master,Tissue=="serum")$MMC_ID))
hasurine <- unique(c(filter(master,Tissue %in% c("urine pellet","urine"))$MMC_ID))

hasall <- filter(haskidney, MMC_ID %in% hasserum & MMC_ID %in% hasurine)


notpks_hasall <- filter(hasall, !(MMC_ID %in% pks$MMC_ID))


