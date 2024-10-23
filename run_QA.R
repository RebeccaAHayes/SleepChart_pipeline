args=commandArgs(trailingOnly = TRUE)

if (length(args)!=1){
  stop("incorrect n arguments supplied")	#don't run if you don't have all arguments
} else {
  analysis_name=args[1]
}

library(kableExtra)
library(tidyverse)
library(data.table)
library(flextable)
library(magrittr)

pdir = "/lab-share/Psych-Jalbrzikowski-PNR-e2/Public/Sandbox/LunaC_pipeline/"
studydir=paste0(pdir,"/derivatives/LunaR_analyses_",analysis_name)
datadir<-paste0(pdir,"/derivatives/formatted_data")

#output directory for compiled data

load(paste0(studydir, "/output/allplotdata.Rdata"))
vars<-read.table(paste0(studydir,"/cmd/vars.txt")) %>% filter(str_detect(V1,"alias")) %>% separate(V2,into=c("varname","signame"),sep="[|]")

# loop through the file list to read in data and clean it up
ss_epoch_subset<-filter(allep, !is.na(STAGE))
ss_epoch_subset$STAGE<- factor(ss_epoch_subset$STAGE, levels=c("W", "N1", "N2", "N3", "R")) 
ss_epoch_subset$stage_num<-recode(ss_epoch_subset$STAGE, W=0, N1=1, N2=2, N3=3, R=5)
tmp<-get(paste0(vars$varname[1],"_nrem_power_df"))
ids <- unique(tmp$ID)

for(n in 1:length(ids)){  
#  n=1
  i=ids[n]  
  rmarkdown::render(input = "QA_report.rmd", 
                    output_format = "html_document",
                    output_file = paste0(i, ".html"),
                    output_dir = paste0(studydir, "/output/"))
  
}
