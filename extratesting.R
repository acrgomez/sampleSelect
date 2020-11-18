#2017 final list consistency test

urilist <- read.csv("~/Google Drive/Sample_select/april2017/urine.csv")


templist <- paste(urilist$MMC_ID,urilist$Specimens..bag)
tempall <- paste(urines$MMC_ID,urines$Specimens..bag)



urilist$collection <- urines$Specimens..collection_method[match(templist,tempall)]
write.csv(urilist,"~/Google Drive/Sample_select/april2017/urine.csv",row.names=F)



extrasera <- read.csv("~/Google Drive/Sample_select/april2017/extrasera_extractedDNA_KP.csv")


serlist <- read.csv("~/Google Drive/Sample_select/april2017/2017.4.15.SeraPullAll.csv")




filter(extrasera,!(MMC_ID %in% serlist$MMC_ID))



extdna <-read.csv("~/Google Drive/Sample_select/april2017/selected_extracted_DNA_nadc.csv")



cdc <- filter(matresults,Lab != "CAHFS")
cahfs <- filter(matresults,Lab == "CAHFS")

extdna$cdcMAT <- extdna$MMC_ID2 %in% cdc$ID
extdna$cahfsonly <- extdna$MMC_ID2 %in% cahfs$ID
extdna$newlist <- extdna$MMC_ID2 %in% serlist$MMC_ID
