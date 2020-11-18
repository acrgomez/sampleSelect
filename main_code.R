keeps$logIct <- log2(keeps$Icterohaemorrhagiae.RGA/100) +1
keeps$logIct[keeps$logIct <0] = 0

keeps$logPom <- log2(keeps$Pomona.P/100) +1
keeps$logPom[keeps$logPom <0] = 0

keeps3 <- filter(keeps,TMMC.ID %in% ids)


keeps2$dayssince[keeps2$TMMC.ID=="Anakin"] = keeps2$Sample.Date[keeps2$TMMC.ID=="Anakin"] - 
                                                keeps2$Sample.Date[keeps2$TMMC.ID=="Anakin"][1]  

keeps2$dayssince[keeps2$TMMC.ID=="Yoda"] = keeps2$Sample.Date[keeps2$TMMC.ID=="Yoda"] - 
    keeps2$Sample.Date[keeps2$TMMC.ID=="Yoda"][1] 


keeps2 <- filter(keeps2,TMMC.ID %in% c("Yoda","Anakin"))

keeps$Admit.Date <- as.Date(keeps$Admit.Date,"%m/%d/%Y")
keeps$dayssince <- keeps$Sample.Date-keeps$Admit.Date
pdf("~/Google Drive/ictpomLK.pdf")
ggplot(keeps2,aes(x=dayssince))+
    geom_line(aes(y=logIct,color="Ictero"))+
    geom_line(aes(y=logPom,color="Pomona"))+
    facet_wrap(~TMMC.ID,nrow = 2,ncol=1)

ggplot(keeps3,aes(x=dayssince))+
    geom_line(aes(y=logIct,color="Ictero"))+
    geom_line(aes(y=logPom,color="Pomona"))+
    facet_wrap(~TMMC.ID)


dev.off()







###boyce      
                   
                   
boycesera <- filter(sera, grepl(pattern = "Boyce", x = Specimens..recipient))

pulllist <- read.csv("~/Google Drive/Sample_select/a_final/list.csv")
fromprev <- read.csv("~/Google Drive/Sample_select/a_final/fromprevious.csv")

boycesera$currentpull <-    boycesera$MMC_ID %in% pulllist$MMC_ID
boycesera$notsentpreviouspull <- boycesera$MMC_ID %in% fromprev$MMC_ID
boycesera$leptosentsera <- boycesera$MMC_ID %in% leptosentsera$MMC_ID
boycesera$hasresults <- boycesera$MMC_ID %in% matdata$ID
                   
                   write.csv(boycesera,"~/Google Drive/Sample_select/boyce_new_all.csv")
                   
                   

require(ggplot2)
require(dplyr)
require(lubridate)

options(stringsAsFactors = F)

master1 <- read.csv("~/Google Drive/data/1970to2017_niceandclean.csv")

savedSerum <- read.csv("~/Dropbox/Sample Selection/new/setAside.csv")
bloodwork <- read.csv("~/Dropbox/Sample Selection/new/bloodwork.csv",na.strings=c("","NA"," "))
matresults <- read.csv("~/Google Drive/CSL_Leptospirosis/Current_Data/MAT_Master.csv")

pcrresults <- read.csv("~/Google Drive/CSL_Leptospirosis/Current_Data/PCR_Master.csv")


master <- select(master1,Admit_Date,
                 MMC_ID,Name,Sex,Age_Class,Disposition,Specimens..amount,
                 Specimens..collection_method,Specimens..date_sampled,
                 Specimens..specimen_tissue,Specimens..storage,Specimens..tissue_quantity,Studies..Date_Sent,
                 Specimens..freezer,Specimens..box,Specimens..bag,Studies..Study,
                 Specimens..sent,Specimens..project,Specimens..recipient,
                 Specimens..results)

master <- filter(master,!is.na(MMC_ID)&!is.na(Specimens..date_sampled)&!is.na(Admit_Date),
                 Specimens..box !="Missing",grepl(MMC_ID,pattern = "CSL-[0-9]"))

master$ID <- sub(x=master$MMC_ID,pattern = "-[A-Z][0-9]",replacement = "-F")
master$ID <- sub(x=master$ID,pattern = "-[A-Z]",replacement = "")

rm(master1)


master$Admit_Date <- as.Date(master$Admit_Date)
master$Studies..Study[is.na(master$Studies..Study)] = ''
master$Specimens..date_sampled <- as.Date(master$Specimens..date_sampled,"%m/%d/%Y")
master$DaysSinceAdmit <- master$Specimens..date_sampled - master$Admit_Date
master <- arrange(master, MMC_ID, Admit_Date)
master$PKSflagged <- agrepl(x = master$Studies..Study,pattern = "Prospective Kidney Study",ignore.case = T)


###agrepl uses fuzzy matching for the tissues##

sera <- filter(master,agrepl(x=master$Specimens..specimen_tissue,pattern ="serum"))
kidneys <- filter(master,agrepl(x=master$Specimens..specimen_tissue,pattern ="kidney"))
urines <- filter(master,grepl(x=master$Specimens..specimen_tissue,pattern ="urine"))


###Sera selection goes first, as others rely on it###
sera$hasresults <- sera$ID %in% matresults$ID
sera$setaside <- sera$ID %in% savedSerum$ID
sera$atTMMC <- !is.na(sera$Specimens..bag)
#sera$restrand <- agrepl("restrand",sera$Name)
tmmcsera <- filter(sera,atTMMC)
tmmcsera$last <- !duplicated(tmmcsera$ID) & !duplicated(tmmcsera$ID,fromLast = T)
nottmmcsera <- filter(sera,!atTMMC)
nottmmcsera$last <- NA
sera <- rbind(tmmcsera,nottmmcsera)
leptosentsera <- filter(nottmmcsera,grepl(Specimens..recipient,pattern = "CDC"))
basicsera <- filter(tmmcsera,!last,!hasresults, !(MMC_ID %in% leptosentsera$MMC_ID),
                    !setaside,Specimens..date_sampled>"2016-01-01",
                    DaysSinceAdmit <=14)

basicsera <- arrange(basicsera,MMC_ID,DaysSinceAdmit)
selectsera <- filter(basicsera,!duplicated(MMC_ID))

temp <- filter(whirl,grepl(pattern = "last",x = whirl$admitsera))
temp2<- filter(whirl,postmortsera == "yes")
extra <- filter(tmmcsera, MMC_ID %in% temp$MMC_ID)
extra2 <- arrange(filter(tmmcsera, MMC_ID %in% temp2$MMC_ID),DaysSinceAdmit)
extra2 <- filter(extra2, !duplicated(MMC_ID,fromLast = T))

selsera <- rbind(selectsera,extra,extra2)



###urine
uri <- filter(urines,Specimens..collection_method %in% c("catheter","cystocentesis","necropsy","post mortem"))
tmmcurine <- filter(uri, !is.na(Specimens..bag))
senturine <- filter(uri, is.na(Specimens..bag), grepl("Leptospira PCR",Specimens..project))

candidateurine <- filter(tmmcurine, !(MMC_ID %in% senturine$MMC_ID))

temp <- candidateurine$Specimens..specimen_tissue
temp[grepl("pellet",temp)] = "a"

candidateurine <- arrange(candidateurine,MMC_ID,temp,DaysSinceAdmit,desc(Specimens..box))

urine <- filter(candidateurine,!duplicated(MMC_ID))

seluri <- filter(urine,MMC_ID %in% c(selsera$MMC_ID,matresults$MMC_ID,leptosentsera$MMC_ID),
                 year(Admit_Date)>2015)

seluri <- filter(seluri,Name!="Elias")

write.csv(seluri,"~/Google Drive/Sample_select/Urine_Pull.csv",row.names=F)






###bloodwork

bloodwork1 <- read.csv("~/Google Drive/Untitled.csv")
bloodwork <- filter(bloodwork1,BloodValues..Sample_Type=="Chemistry")


filledIndexes <- bloodwork$MMC_ID != ''
for (col in seq_along(colnames(bloodwork))){
    bloodwork[,col][filledIndexes&as.character(bloodwork[,col])==''] = NA             
}


bloodwork$Admit_Date <- as.Date(populate(bloodwork$Admit_Date,filledIndexes), "%m/%d/%Y")
bloodwork$MMC_ID <- populate(bloodwork$MMC_ID,filledIndexes)
bloodwork$Name <- populate(bloodwork$Name,filledIndexes)
bloodwork$SampleDate <- as.Date(bloodwork$BloodValues..Sampled_Date,"%m/%d/%Y")



bloodwork$DaysSinceAdmit <- bloodwork$SampleDate- bloodwork$Admit_Date

bloodwork <- filter(bloodwork, !is.na(BloodValues..Albumin))

selsera$bloodwork = NA

for (sample in seq_along(selsera$MMC_ID)){
    if (!(selsera$MMC_ID[sample] %in% bloodwork$MMC_ID)) {
        selsera$bloodwork[sample] <- "NEEDED"
    } else {
        temp1 <- filter(bloodwork,MMC_ID == selsera$MMC_ID[sample])
        if (sum(abs(temp1$DaysSinceAdmit - 
                    selsera$DaysSinceAdmit[sample]) < 4) > 0) {
            selsera$bloodwork[sample] <- "OK"
        } else {
            selsera$bloodwork[sample] <- "NEEDED (DATE)"
        }
    }
}
                    
    



write.csv(selsera,file="~/Google Drive/Sample_select/Sera_Pull.csv",row.names=F)


####
####
####
###


nadc <- filter(master,grepl(x = Specimens..recipient,pattern = "NADC"))
nadc$clean <- nadc$Specimens..collection_method %in% c("catheter",
                                                       "cystocentesis",
                                                       "necropsy","post mortem")
nadc$anymatresult <- nadc$MMC_ID %in% matresults$ID
nadc$dupli <- duplicated(nadc$MMC_ID) | duplicated(nadc$MMC_ID,fromLast = T)
matresults$dupli <- duplicated(matresults$ID) | duplicated(matresults$ID, 
                                                           fromLast = T)


index <- match(nadc$MMC_ID,matresults$ID)
nadc$matresdate <- as.Date(matresults$SampleDateMAT[index])

matres2 <- filter(matresults,duplicated(ID))
index2 <- match(nadc$MMC_ID,matres2$ID)
nadc$matresdate2 <- as.Date(matres2$SampleDateMAT[index2])

matres3 <- filter(matres2,duplicated(ID))
index3 <- match(nadc$MMC_ID,matres3$ID)
nadc$matresdate3 <- as.Date(matres3$SampleDateMAT[index3])

matres4 <- filter(matres3,duplicated(ID))
index4 <- match(nadc$MMC_ID,matres4$ID)
nadc$matresdate4 <- as.Date(matres4$SampleDateMAT[index4])

###next has nothing, so I can stop here

temp1 <- nadc[,c(9,9,9,9)] - nadc[,28:31]


temp1$min <- apply(temp1,MARGIN = 1,FUN=min,na.rm=T)


nadc <- bind(nadc, temp1)


nadc$diffdays <- as.numeric(gsub(pattern = '[^0-9]','',nadc$min))


nadc$mat_within_4_days <- nadc$diffdays < 5

nadc2 <- filter(nadc, grepl(pattern = "PCR",x=nadc$Specimens..project),clean)

write.csv(nadc2,"~/Google Drive/Sample_select/nadc_samples.csv",row.names=F)



hollings <- filter(master,grepl(x = Specimens..recipient,pattern = "olling"))

bla <- filter(nadc2, MMC_ID %in% hollings$MMC_ID)

bla2 <- filter(hollings, MMC_ID %in% nadc2$MMC_ID)


