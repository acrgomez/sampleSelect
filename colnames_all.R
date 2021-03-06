##colnames pulled directly from the old system
##many of them were maintained but this is not up to date and should not be used w/o verification

colss<-c("Admit_Date","AdmitLatLong_ActEst","Age_Class","ATOS","Body_of_Water","BornInRehab","Common_Name","Disposition",
"FirstObs","Genus","Lat_Long_Valid","LatLong_Determinded","Lock","LockDisplay","Mass_Stranding","MassStrandRef",
"MMC_ID","Mother_ID","Name","Name_Definition_Explanation","Name_Fundraising_Only","NMFS_Disp_Form","NMFS_NDB",
"NMFS_Strand_Form","No_Report","Number_of_animals","Other_Markings","PhotoDispo","PhotosVideos",
"PhotosVideos_PreAdmit","Recovery_Area","Remarks","Restrand","SamplesCollected","Sex","Species",
"SpecimenDisposition","SpecimensCollected","Strand_Date","Strand_Latitude_new","Strand_Longitude_new",
"Stranding_City","Stranding_County","Stranding_Details","Stranding_Group","Stranding_Locality","TC",
"TransferFrom","zz_created_By","zz_created_TS","zz_modified_By","zz_modified_TS","Other_ID..ID","Specimens..amount",
"Specimens..bag","Specimens..box","Specimens..cleanup","Specimens..collection_method","Specimens..date_sampled",
"Specimens..details","Specimens..flag_onReport","Specimens..freezer","Specimens..initials","Specimens..project",
"Specimens..recipient","Specimens..report_num","Specimens..results","Specimens..results_details","Specimens..sent",
"Specimens..specimen_tissue","Specimens..storage","Specimens..sum_RecCount","Specimens..tissue_quantity",
"Specimens..zz_rw_tmpFlag_missingPostMortDate","Studies..Date_Recieved","Studies..Date_Sent",
"Studies..Facility_Address","Studies..flg_TempTransfer","Studies..Primary_Investigator",
"Studies..Primary_Investigator_Affiliation","Studies..Primary_Investigator_Permit_Number","Studies..Study",
"Tags..Color","Tags..Date","Tags..ID","Tags..Placement","Tags..TagApplied","Tags..Type")


pullcols<-c("Admit_Date",
            "MMC_ID","Name","Sex","Age_Class","Disposition","Specimens..amount",
            "Specimens..collection_method","Specimens..date_sampled",
            "Specimens..specimen_tissue","Specimens..storage","Specimens..tissue_quantity","Studies..Date_Sent",
            "Specimens..freezer","Specimens..box","Specimens..bag","Studies..Study",
            "Specimens..sent","Specimens..project","Specimens..recipient",
            "Specimens..results",
            "Specimens..results_details",
            "Specimens..details")