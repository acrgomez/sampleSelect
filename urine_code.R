selsera <- read.csv("~/Google Drive/Sample_select/april2017/2017.4.15.SeraPullAll.csv")
matresults <- read.csv("~/Google Drive/CSL_Leptospirosis/Current_Data/MAT_Master.csv")
pcrresults <- read.csv("~/Google Drive/CSL_Leptospirosis/Current_Data/PCR_Master.csv")
fromprev <- read.csv("~/Google Drive/Sample_select/april2017/a_final/fromprevious.csv")
uri <- filter(urines,Specimens..collection_method %in% c("catheter","cystocentesis","necropsy","post mortem"))

tmmcurine <- filter(uri, !is.na(Specimens..bag))
senturine <- filter(uri, is.na(Specimens..bag), grepl("Leptospira PCR",Specimens..project))


###not sent
candidateurine <- filter(tmmcurine, !(MMC_ID %in% senturine$MMC_ID))

temp <- candidateurine$Specimens..specimen_tissue
temp[grepl("pellet",temp)] = "a"

candidateurine <- arrange(candidateurine,MMC_ID,temp,DaysSinceAdmit,desc(Specimens..box))

urine <- filter(candidateurine,!duplicated(MMC_ID))



urine$hasMAT <- matresults$Lab[match(urine$MMC_ID,matresults$ID)]
urine$newlist <- urine$MMC_ID %in% selsera$MMC_ID
urine$oldlist <- urine$MMC_ID %in% fromprev$MMC_ID

temp3 <-filter(pcrresults,ID %in% urine$MMC_ID)
urine$haspcr <- urine$MMC_ID %in% pcrresults$ID

write.csv(urine,"~/Google Drive/Sample_select/temp.csv",row.names = F)

write.csv(temp3, "~/Google Drive/Sample_select/temp3.csv",row.names = F)

seluri <- filter(urine,MMC_ID %in%
                     c(selsera$MMC_ID,matresults$MMC_ID,fromprev$MMC_ID))

seluri <- filter(seluri,Name!="Elias")






cafsres <- filter(sera,agrepl("CAHFS",sera$Specimens..recipient))



list1 <- read.csv("~/Google Drive/Sample_select/april2017/temp1.csv")


cahfsinlist <- filter(cafsres,MMC_ID %in% list1$MMC_ID)
cahfsres1 <- filter(cahfsinlist,duplicated(MMC_ID))
cahfsres2 <- filter(cahfsres1,duplicated(MMC_ID))
cahfsres3 <- filter(cahfsres2,duplicated(MMC_ID))
cahfsres4 <- filter(cahfsres3,duplicated(MMC_ID))
cahfsres5 <- filter(cahfsres4,duplicated(MMC_ID))


idmatch <- match(list1$MMC_ID,cahfsinlist$MMC_ID)
list1$cahfsdate1 <- cahfsinlist$Specimens..date_sampled[idmatch]
list1$cahfsres1 <- cahfsinlist$Specimens..results[idmatch]


idmatch <- match(list1$MMC_ID,cahfsres1$MMC_ID)
list1$cahfsdate2 <- cahfsres1$Specimens..date_sampled[idmatch]
list1$cahfsres2 <- cahfsres1$Specimens..results[idmatch]


idmatch <- match(list1$MMC_ID,cahfsres2$MMC_ID)
list1$cahfsdate3 <- cahfsres2$Specimens..date_sampled[idmatch]
list1$cahfsres3 <- cahfsres2$Specimens..results[idmatch]


idmatch <- match(list1$MMC_ID,cahfsres3$MMC_ID)
list1$cahfsdate4 <- cahfsres3$Specimens..date_sampled[idmatch]
list1$cahfsres4 <- cahfsres3$Specimens..results[idmatch]


idmatch <- match(list1$MMC_ID,cahfsres4$MMC_ID)
list1$cahfsdate5 <- cahfsres4$Specimens..date_sampled[idmatch]
list1$cahfsres5 <- cahfsres4$Specimens..results[idmatch]

idmatch <- match(list1$MMC_ID,cahfsres5$MMC_ID)
list1$cahfsdate6 <- cahfsres5$Specimens..date_sampled[idmatch]
list1$cahfsres6 <- cahfsres5$Specimens..results[idmatch]


write.csv(list1,"~/Google Drive/Sample_select/april2017/temp1.csv",row.names=F)



