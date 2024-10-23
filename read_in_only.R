library(tidyverse)
library(data.table)
library(lubridate)
library(RSQLite)

pdir = "/lab-share/Psych-Jalbrzikowski-PNR-e2/Public/Sandbox/LunaC_pipeline"

#########################
#Data Processing
#########################

datadir<-paste0(pdir,"/derivatives/formatted_data")

bids <- paste0(pdir,"/bids")
erdir <- paste0(datadir,"/edf_raw")
if (!dir.exists(erdir)){dir.create(erdir,recursive=T,mode="0775")} else {print("dir exists")}

# find the files that you want
edfs <- list.files(path=bids, pattern="*.edf",recursive=T,full.names=T)
# copy the files to the new folder
file.copy(edfs, paste0(erdir,"/",str_remove(basename(edfs),"_eeg")),recursive=FALSE)
setwd(datadir)

#Make a list of annotation files in bids and move them to the derivatives folder
annots<-list.files(path=bids, pattern="events.csv$",recursive=T,full.names=T)
ardir <- paste0(datadir,"/annot_raw")
if (!dir.exists(ardir)){dir.create(ardir,recursive=T,mode="0775")} else {print("dir exists")}						 
file.copy(annots, paste0(ardir,"/",str_remove(basename(annots),"_events")),recursive=FALSE)

#direct to folder with raw annotation files exported from dreem
setwd(ardir)

#get list of file names
annotfile_list<-dir(pattern="*.csv")

#Rewrite annotations so that Luna can read them
for (af in annotfile_list){
  df<-read.csv(af) %>% select(sleep_stage,times) %>% #Select sleep stage and start time
    mutate(stop="...") %>% #Set end time to "beginning of next annotation"
    select("class"=sleep_stage,"start"=times,stop)
  write.table(df,file=str_replace(af,".csv",".annot"),quote=F,sep="\t",row.names=F,col.names=T) #Save as tab-separated annot file

}

annotfile_list<-list.files(ardir,pattern="*.annot")
#get subject info from each file
subinfo<-as_tibble(annotfile_list) %>% 
  mutate("studyid"=str_extract(value,"sub-.*?(?=_)"),
         "ses"=str_extract(value,"ses-[0-9]+"),
         "run"=str_extract(value,"run-[0-9]+"))


fdir<-paste0(datadir,'/final/')
if (!dir.exists(fdir)){dir.create(fdir,recursive=T,mode="0775")} else {print("dir exists")}
dfdir<-paste0(datadir,'/denoised_final/')
if (!dir.exists(dfdir)){dir.create(dfdir,recursive=T,mode="0775")} else {print("dir exists")}
#nfdir<-paste0(datadir,'/notched_final/')
#if (!dir.exists(nfdir)){dir.create(nfdir,recursive=T,mode="0775")} else {print("dir exists")}
#loop through each subject and recode sleep stages
for (i in 1:length(annotfile_list)){
  #new annot name
  newfile<-paste0(subinfo$studyid[i], "_",
                  subinfo$ses[i],"_",
                  subinfo$run[i],".annot")
  newdf<-read_delim(paste0(ardir,"/",subinfo$value[i])) %>%
    mutate("stop"="...")
  
  write.table(newdf, paste0(fdir,newfile),sep="\t",row.names=F,quote=F,col.names=T)
  write.table(newdf$class, paste0(fdir,str_replace(newfile,".annot",".eannot")),row.names=F,quote=F,col.names=F)
}

print(paste0("N edfs: ",length(dir(erdir,pattern="*.edf$"))))
#EDIT POINT FOR NEW PROJECTS
setwd(datadir)

file.remove(paste0("final/",dir("final",pattern="*.edf$")))
file.remove(paste0("denoised_final/",dir("denoised_final",pattern="*.edf$")))

system(paste0("luna --build ",datadir,"/edf_raw > rawedfs.lst"))
system("(luna rawedfs.lst < cmd/rereference.txt) 2>&1 | tee reref-out.txt")
system("(luna rawedfs.lst < cmd/rereference_denoise.txt) 2>&1 | tee reref-denoise-out.txt")
#system("(luna rawedfs.lst < cmd/rereference_notch.txt) 2>&1 | tee reref-notch-out.txt")
system(paste0("luna --build ",datadir,"/final > ref.lst"))
system(paste0("luna --build ",datadir,"/denoised_final > ref-den.lst"))
#system("luna --build /project/derivatives/formatted_data/notched_final > ref-notch.lst")

print(paste0("N re-referenced edfs: ",length(str_subset(readLines("ref.lst"),".*ref.*"))))

setwd(fdir)
#get list of rereferenced EDF files
edffile_list<-dir(pattern="*.edf")

#get subject info from edf files and generate new, shorter edf names
e_subinfo<-as_tibble(edffile_list) %>% 
  mutate("studyid"=str_extract(value,"sub-.*?(?=_)"),
         "ses"=str_extract(value,"ses-[0-9]+"),
         "run"=str_extract(value,"run-[0-9]+"))%>%
  mutate(newname=paste0(studyid,"_",ses,"_",run,".edf"))

#Rename files
file.rename(from = paste0(fdir,edffile_list),
            to = paste0(fdir, e_subinfo$newname))

setwd(dfdir)
#get list of rereferenced EDF files
edffile_list<-dir(pattern="*.edf")

#get subject info from edf files and generate new, shorter edf names
e_subinfo<-as_tibble(edffile_list) %>% 
  mutate("studyid"=str_extract(value,"sub-.*?(?=_)"),
         "ses"=str_extract(value,"ses-[0-9]+"),
         "run"=str_extract(value,"run-[0-9]+"))%>%
  mutate(newname=paste0(studyid,"_",ses,"_",run,".edf"))

#Rename files
file.rename(from = paste0(dfdir,edffile_list),
            to = paste0(dfdir, e_subinfo$newname))


setwd(datadir)

#Make folders for sleep stage predictions and confidence ratings
apdir <- paste0(datadir,"/annot_predicted")
if (!dir.exists(apdir)){dir.create(apdir,recursive=T,mode="0775")} else {print("dir exists")}						 
acdir <- paste0(datadir,"/annot_confidence")
if (!dir.exists(acdir)){dir.create(acdir,recursive=T,mode="0775")} else {print("dir exists")}		

#dapdir <- paste0(datadir,"/denoised_annot_predicted")
#if (!dir.exists(dapdir)){dir.create(dapdir,recursive=T,mode="0775")} else {print("dir exists")}						 
#dacdir <- paste0(datadir,"/denoisedannot_confidence")
#if (!dir.exists(dacdir)){dir.create(dacdir,recursive=T,mode="0775")} else {print("dir exists")}		

#napdir <- paste0(datadir,"/notched_annot_predicted")
#if (!dir.exists(napdir)){dir.create(napdir,recursive=T,mode="0775")} else {print("dir exists")}						 
#nacdir <- paste0(datadir,"/notchedannot_confidence")
#if (!dir.exists(nacdir)){dir.create(nacdir,recursive=T,mode="0775")} else {print("dir exists")}		

#Generate Luna-POPS predicted sleep stages
system(paste0("luna --build ",datadir,"/final -ext=.null>pops.lst"))
system("luna pops.lst chan=fp2d_fp1d -t annot_confidence/pops_fp2d_fp1d < cmd/pops-cmd.txt")
system("luna pops.lst chan=fp2d_tp7d -t annot_confidence/pops_fp2d_tp7d < cmd/pops-cmd.txt")
system("luna pops.lst chan=fp1d_tp8d -t annot_confidence/pops_fp1d_tp8d < cmd/pops-cmd.txt")

#system("luna --build /project/derivatives/formatted_data/denoised_final -ext=.null>pops_denoised.lst")
#system("luna pops_denoised.lst chan=fp2d_fp1d -t denoised_annot_confidence/pops_fp2d_fp1d < cmd/denoised_pops_cmd.txt")
#system("luna pops_denoised.lst chan=fp2d_tp7d -t denoised_annot_confidence/pops_fp2d_tp7d < cmd/denoised_pops_cmd.txt")
#system("luna pops_denoised.lst chan=fp1d_tp8d -t denoised_annot_confidence/pops_fp1d_tp8d < cmd/denoised_pops_cmd.txt")

#system("luna --build /project/derivatives/formatted_data/notched_final -ext=.null>pops_notched.lst")
#system("luna pops_notched.lst chan=fp2d_fp1d -t notched_annot_confidence/pops_fp2d_fp1d < cmd/notched_pops_cmd.txt")
#system("luna pops_notched.lst chan=fp2d_tp7d -t notched_annot_confidence/pops_fp2d_tp7d < cmd/notched_pops_cmd.txt")
#system("luna pops_notched.lst chan=fp1d_tp8d -t notched_annot_confidence/pops_fp1d_tp8d < cmd/notched_pops_cmd.txt")


