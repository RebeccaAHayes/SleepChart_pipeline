args=commandArgs(trailingOnly = TRUE)

if (length(args)!=2){
  stop("incorrect n arguments supplied")	#don't run if you don't have all arguments
} else {
  analysis_name=args[1]
  lstfile=args[2]
}

library(tidyverse)
library(data.table)
library(lubridate)

pdir = "/lab-share/Psych-Jalbrzikowski-PNR-e2/Public/Sandbox/LunaC_pipeline/"
studydir=paste0(pdir,"/derivatives/LunaR_analyses_",analysis_name)
datadir<-paste0(pdir,"/derivatives/formatted_data")

#########################
#SAMPLE LIST
#########################
print("updating lists")

setwd(datadir)

system(paste0("luna --build ",datadir,"/final -ext=.annot > muse.lst")) #Generate final subject list of rereferenced edfs and eannot files
system(paste0("luna --build ",datadir,"/final -ext=.eannot > muse_eannot.lst")) #Generate final subject list of rereferenced edfs and eannot files
read_tsv("muse_eannot.lst",col_names = c("sub","edf","eannot")) %>% filter(!str_detect(sub,"sub-T0009")) %>% write_tsv("muse_eannot.lst",col_names=FALSE)

read_tsv("muse.lst",col_names = c("sub","edf","eannot")) %>% mutate(edf=str_replace(edf,"final","denoised_final")) %>% write_tsv("muse_denoised.lst",col_names=FALSE)
read_tsv("muse_eannot.lst",col_names = c("sub","edf","eannot")) %>% mutate(edf=str_replace(edf,"final","denoised_final")) %>% write_tsv("muse_eannot_denoised.lst",col_names=FALSE)

print("lists updated")

#########################
#STAGING AND ANALYSES
#########################
setwd(studydir)

#output directory for compiled data
todaysdate= (Sys.Date())
outputdir = paste0(studydir, "/output/")
if (!dir.exists(outputdir)){dir.create(outputdir,recursive=T,mode="0775")} else {print("outputdir exists")}

#-------------------------------------------------------------------------------
# BEGIN ANALYSIS
#-------------------------------------------------------------------------------
#Read in sample list
setwd(studydir) 
getwd()

#sl <- lsl(paste0("/project/derivatives/formatted_data/",lstfile))
if (!dir.exists("spindles")){dir.create("spindles",mode="0775")} else {print("spindles dir exists")} #Make directory for spindle annots
vars<-read.table("cmd/vars.txt") %>% 	#Get list of channels we'll be looping through from vars.txt aliases
      filter(str_detect(V1,"alias")) %>% 
      separate(V2,into=c("varname","signame"),sep="[|]")

sublist<-read_tsv(paste0(datadir,"/",lstfile,".lst"),col_names=c("ID","edf","annot"))

#-------------------------------------------------------------------------------
#sleep stages in all subjects
#-------------------------------------------------------------------------------

#build luna command
lunacmd<-paste0(
  "luna ",
  datadir,"/",lstfile,".lst", #lst file
  " @cmd/vars.txt ",
  "sig=", vars$varname[1],
  " -t ",outputdir,"/staging -s 'HYPNO epoch'"
)

system(lunacmd)

#GENERATE SLEEP ARCHITECTURE AND STAGING INFORMATION

sleep_arch<-rbindlist(lapply(list.files(path=paste0(outputdir,"/staging"),pattern="HYPNO.txt",recursive=T,full.names=T),fread))
sleep_stages<-rbindlist(lapply(list.files(path=paste0(outputdir,"/staging"),pattern="HYPNO_SS.txt",recursive=T,full.names=T),fread))
sleep_cycles<-rbindlist(lapply(list.files(path=paste0(outputdir,"/staging"),pattern="HYPNO_C.txt",recursive=T,full.names=T),fread))
ss_epoch<-rbindlist(lapply(list.files(path=paste0(outputdir,"/staging"),pattern="HYPNO_E.txt",recursive=T,full.names=T),fread))

ss_minutes<-sleep_stages[ sleep_stages$SS %in% c("N1","N2","N3","R") ,c("ID","SS","MINS")] #Extract out the minute estimates
ss_minutes_w<-pivot_wider(ss_minutes, id_cols = "ID", names_from="SS", names_prefix="MINS.",values_from="MINS") %>% #Pivot wider 
		mutate(TST=MINS.N1+MINS.N2+MINS.N3+MINS.R, #Calculate TST based on only N1/N2/N3/R and get percentages
		       PERC.N1=100*MINS.N1/TST,
		       PERC.N2=100*MINS.N2/TST,
		       PERC.N3=100*MINS.N3/TST,
		       PERC.R=100*MINS.R/TST)

#store data
write.csv(ss_minutes_w,paste0(outputdir,"/sleep-stages_", todaysdate, ".csv"),row.names=FALSE )
write.csv(sleep_arch, paste0(outputdir,"/sleep-arch_", todaysdate, ".csv"),row.names=FALSE)
write.csv(sleep_cycles, paste0(outputdir,"/stages-by-cycle_", todaysdate, ".csv"),row.names=FALSE)
write.csv(ss_epoch, paste0(outputdir,"/stages-by-epoch_", todaysdate, ".csv"),row.names=FALSE)

remove(list=c("ss_minutes","ss_minutes_w","sleep_arch","sleep_cycles"))

#-------------------------------------------------------------------------------
#spindle and power info in all subjects
#-------------------------------------------------------------------------------


for (v in vars$varname){ #Now, loop through each of our channels
  
  #GENERATE SPINDLE INFORMATION FROM EACH OF OUR COILS
  lunacmd<-paste0(
    "luna ",
    datadir,"/",lstfile,".lst", #lst file
    " @cmd/vars.txt ",
    "sig=",v," chan=",v,
    " -t ",outputdir,"/nrem_analyses_",v,
    " < cmd/nrem_power_spectra_so.txt"
  )
  
  system(lunacmd)
  #----------------------------------------------------------------------------
  #NREM analyses
  #----------------------------------------------------------------------------
  
  nrem_power_df<-rbindlist(lapply(list.files(path=paste0(outputdir,"/nrem_analyses_",v),pattern="PSD_F_CH.txt",recursive=T,full.names=T),fread)) %>%
    filter(F>=0.5)
  write.csv(nrem_power_df,paste0(outputdir,"/",v,"_NREM_PSD_", todaysdate, ".csv"),row.names=FALSE) #Save that data
  
  ep2n<-rbindlist(lapply(list.files(path=paste0(outputdir,"/nrem_analyses_",v),pattern="DUMP-MASK_E.txt",recursive=T,full.names=T),fread))
  write.csv(ep2n,paste0(outputdir,"/",v,"_NREM_masked_epochs_",todaysdate,".csv"),row.names=FALSE) #Save that data
  names(ep2n)<-c("ID","E",paste0(v,"_EMASK_N")) #Rename it so we can save all the channels' info in an Rdata file at the end of the script
  
  nrem_bands<-rbindlist(lapply(list.files(path=paste0(outputdir,"/nrem_analyses_",v),pattern="PSD_B_CH.txt",recursive=T,full.names=T),fread))
  nrem_bands_wide<-pivot_wider(nrem_bands, id_cols = c("ID", "CH"), names_from = "B", values_from=c("PSD","RELPSD"),names_sep = ".")  #Reshape and save
  write.csv(nrem_bands_wide,paste0(outputdir,"/",v,"_NREM_BANDS_", todaysdate, ".csv"),row.names=FALSE)

  
  nrem_bands_by_epoch<-rbindlist(lapply(list.files(path=paste0(outputdir,"/nrem_analyses_",v),pattern="PSD_E_B_CH.txt",recursive=T,full.names=T),fread))
  write.csv(nrem_bands_by_epoch, paste0(outputdir, "/",v,"NREM_BANDS_BY_EPOCH_", todaysdate, ".csv"), row.names=FALSE)
  
  nrem_spectra_df<-rbindlist(lapply(list.files(path=paste0(outputdir,"/nrem_analyses_",v),pattern="PSD_E_F_CH.txt",recursive=T,full.names=T),fread)) %>% 
    filter(F>=0.5) #Extract out the by-epoch spectral data for making spectrograms later
  
  nrem_epoch<-rbindlist(lapply(list.files(path=paste0(outputdir,"/nrem_analyses_",v),pattern="MASK_EMASK.txt",recursive=T,full.names=T),fread))
  write.csv(nrem_epoch,paste0(outputdir,"/",v,"_NREM_EPOCHS_", todaysdate, ".csv"),row.names=FALSE) #Extract out total counts of masked and unmasked epochs 
  
  #Extract and save spindle data
  nrem_spindle <- rbindlist(lapply(list.files(path=paste0(outputdir,"/nrem_analyses_",v),pattern="SPINDLES_F_CH.txt",recursive=T,full.names=T),fread))
  write.csv(nrem_spindle, paste0(outputdir,"/",v,"_spindles_", todaysdate, ".csv"),row.names=FALSE)

  #Extract and save SO data
  nrem_so <- rbindlist(lapply(list.files(path=paste0(outputdir,"/nrem_analyses_",v),pattern="SO_CH.txt",recursive=T,full.names=T),fread))
  write.csv(nrem_so, paste0(outputdir,"/",v,"_slowwave_", todaysdate, ".csv"),row.names=FALSE)
  
  #----------------------------------------------------------------------------
  #REM analyses
  #----------------------------------------------------------------------------

  #Everything repeats for REM
  
  #GENERATE SPINDLE INFORMATION FROM EACH OF OUR COILS
  lunacmd<-paste0(
    "luna ",
    datadir,"/",lstfile,".lst", #lst file
    " @cmd/vars.txt ",
    "sig=",v," chan=",v,
    " -t ",outputdir,"/rem_analyses_",v,
    " < cmd/nrem_power_spectra_so.txt"
  )
  
  system(lunacmd)
 
  rem_power_df<-rbindlist(lapply(list.files(path=paste0(outputdir,"/rem_analyses_",v),pattern="PSD_F_CH.txt",recursive=T,full.names=T),fread)) %>%
    filter(F>=0.5)
  write.csv(rem_power_df,paste0(outputdir,"/",v,"_REM_PSD_", todaysdate, ".csv"),row.names=FALSE) #Save that data
  
  ep2r<-rbindlist(lapply(list.files(path=paste0(outputdir,"/rem_analyses_",v),pattern="DUMP-MASK_E.txt",recursive=T,full.names=T),fread))
  write.csv(ep2r,paste0(outputdir,"/",v,"_REM_masked_epochs_",todaysdate,".csv"),row.names=FALSE) #Save that data
  names(ep2r)<-c("ID","E",paste0(v,"_EMASK_R")) #Rename it so we can save all the channels' info in an Rdata file at the end of the script
  
  rem_bands<-rbindlist(lapply(list.files(path=paste0(outputdir,"/rem_analyses_",v),pattern="PSD_B_CH.txt",recursive=T,full.names=T),fread))
  rem_bands_wide<-pivot_wider(rem_bands, id_cols = c("ID", "CH"), names_from = "B", values_from=c("PSD","RELPSD"),names_sep = ".")  #Reshape and save
  write.csv(rem_bands_wide,paste0(outputdir,"/",v,"_REM_BANDS_", todaysdate, ".csv"),row.names=FALSE)
  
  rem_bands_by_epoch<-rbindlist(lapply(list.files(path=paste0(outputdir,"/rem_analyses_",v),pattern="PSD_E_B_CH.txt",recursive=T,full.names=T),fread))
  write.csv(rem_bands_by_epoch, paste0(outputdir, "/",v,"REM_BANDS_BY_EPOCH_", todaysdate, ".csv"), row.names=FALSE)
  
  rem_spectra_df<-rbindlist(lapply(list.files(path=paste0(outputdir,"/rem_analyses_",v),pattern="PSD_E_F_CH.txt",recursive=T,full.names=T),fread)) %>%
    filter(F>=0.5) #Extract out the by-epoch spectral data for making spectrograms later
  
  rem_epoch<-rbindlist(lapply(list.files(path=paste0(outputdir,"/rem_analyses_",v),pattern="MASK_EMASK.txt",recursive=T,full.names=T),fread))
  write.csv(rem_epoch,paste0(outputdir,"/",v,"_REM_EPOCHS_", todaysdate, ".csv"),row.names=FALSE) #Extract out total counts of masked and unmasked epochs 
  
  #----------------------------------------------------------------------------
  #Cleanup
  #----------------------------------------------------------------------------

  #Rename things in the R environment so that we can save them in allplotdata.Rdata at the end
  assign(paste0(v,"_nrem_power_df"),nrem_power_df)
  assign(paste0(v,"_nrem_spectra_df"),nrem_spectra_df)
  assign(paste0(v,"_rem_power_df"),rem_power_df)
  assign(paste0(v,"_rem_spectra_df"),rem_spectra_df)
  assign(paste0(v,"_ep2n"),ep2n)
  assign(paste0(v,"_ep2r"),ep2r)
  
  #Clear out the old dataframe names so we can use them again in the next iteration without worrying about spillover
  remove(nrem_power_df,nrem_spectra_df,rem_power_df,rem_spectra_df,ep2n,ep2r,nrem_spindle,nrem_so,
         nrem_bands,nrem_bands_wide,nrem_epoch,nrem_epoch_wide,nrem_power,
         rem_bands,rem_bands_wide,rem_epoch,rem_epoch_wide,rem_power)
}

#########################
#SAVING RDATA
#########################

#Make a dataframe that combines by-epoch mask info for all channels
allep<-select(ss_epoch, ID, E, STAGE)
for (df in ls(pattern="*_ep2[nr]")){
  print(df)
  allep<-full_join(allep,get(df))
}

#Save data to allplotdata.Rdata
save(list=c("ss_epoch","allep",
            ls(pattern=".*_nrem_power_df"),
            ls(pattern=".*_nrem_spectra_df"),
            ls(pattern=".*_rem_power_df"),
            ls(pattern=".*_rem_spectra_df")),
     file=paste0(studydir, "/output/allplotdata.Rdata"))

print("end run_analyses.R")
