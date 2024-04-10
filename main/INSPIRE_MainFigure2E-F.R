library(tidyverse)
library(survminer)
library(survival)

### Input data preparation ----------------------------------------------------------------------------------------
data_path <- "https://raw.githubusercontent.com/pughlab/INSPIRE_TCR/main/Data"

### Clinical data:
clinical_data_fname <- "INSPIRE_ClinicalData.csv"

clinical_data <- readr::read_csv(
        file.path (data_path , clinical_data_fname))%>%
        dplyr::select(Patient_id , 
                      COHORT , 
                      `Best response longevity` ,
                      `OS (months)` ,
                      `PFS (months)`,
                      `Lost to follow up` ,
                      `PFS Status` , 
                      `OS Status`)

### Diversity FMM clusters:
FMM_data_fname <- "INSPIRE_Tumour_DiversityFMMClusters.csv"

FMM_data <- readr::read_csv(
        file.path (data_path , FMM_data_fname))%>%
        dplyr::select(Patient_id ,  `FMM Cluster`) %>%
        rename(FMM_Cluster = `FMM Cluster`)

### Dataframe for KM assessments:
km_assessment <- left_join(FMM_data ,
                           clinical_data,
                           by = "Patient_id") %>%
        mutate(FMM_Cluster = fct_relevel(
                FMM_Cluster , 
                "Diversified" ,
                "Stagnant diversity" ,
                "Persistently low diversity" ))
### OS ASSESSMENT -------------------------------------------------------------------------------------------------
os_km <- survminer::ggsurvplot (
        survfit (Surv(`OS (months)`, 
                      `OS Status` ,
                      type = "right") ~ FMM_Cluster, 
                 data = km_assessment  ) ,
        ###linewidth:
        size = 0.5 ,
        #linecolor:
        palette = c(  
                "#4F9900" ,  #"Diversified"
                "#BDB4B3" ,  #"Stagnant diversity" , 
                "#3B3029"    #"Persistently low diversity",
        ) ,
        #linetype:
        linetype = c("solid" , "solid" , "solid") ,
        
        #Adding p-value:
        pval = FALSE ,
        pval.size = 5 ,
        pval.coord = c(0 , 0.05) ,
        
        #Fine-tuning legend features:
        legend = "right" ,
        legend.title = "Diversity shift" ,
        legend.labs = c("Diversified" ,
                        "Stagnant diversity" ,
                        "Persistently low diversity" ),
        
        #Risk table:
        risk.table = TRUE,
        risk.table.col = "strata" ,
        risk.table.height = 0.25 ,
        risk.table.y.title = FALSE ) +
        ylab ("Overall survival probability") +
        xlab ("Months")



os_km$plot <- os_km$plot + 
        theme_minimal() +
        theme(
                panel.grid = element_blank() ,
                
                axis.line = element_line(color = "#000000" , linewidth = 0.1),
                axis.ticks = element_line(color = "#000000" ) ,
                
                axis.title = element_text ( size = 13, family = "Helvetica", face="plain" , colour = "#000000"),
                
                axis.text = element_text ( size = 13, family = "Helvetica", face="plain" ,
                                           hjust = 0.5, colour = "#000000"),
                
                legend.title = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 10  ),
                legend.text  = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 10) ,
                legend.position= "none" ,
                aspect.ratio = 0.65)



os_km$table <- os_km$table + 
        theme_minimal() +
        theme(panel.grid = element_blank() ,
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              
              axis.text.y = element_text ( size = 13, family = "Helvetica", face="plain" ,
                                           colour = "#000000"),
              axis.text.x = element_blank(),
              legend.position= "none" ,
              aspect.ratio = 0.176)



print(os_km)
### PFS ASSESSMENT ------------------------------------------------------------------------------------------------
pfs_km <- survminer::ggsurvplot (
        survfit (Surv(`PFS (months)`, 
                      `PFS Status` ,
                      type = "right") ~ FMM_Cluster, 
                 data = km_assessment  ) ,
        ###linewidth:
        size = 0.5 ,
        #linecolor:
        palette = c(  
                "#4F9900" ,  #"Diversified"
                "#BDB4B3" ,  #"Stagnant diversity" , 
                "#3B3029"    #"Persistently low diversity",
        ) ,
        #linetype:
        linetype = c("solid" , "solid" , "solid") ,
        
        #Adding p-value:
        pval = FALSE ,
        pval.size = 5 ,
        pval.coord = c(0 , 0.05) ,
        
        #Fine-tuning legend features:
        legend = "right" ,
        legend.title = "Diversity shift" ,
        legend.labs = c("Diversified" ,
                        "Stagnant diversity" ,
                        "Persistently low diversity" ),
        #Risk table:
        risk.table = TRUE,
        risk.table.col = "strata" ,
        risk.table.height = 0.25 ,
        risk.table.y.title = FALSE ) +
        ylab ("Progression-free survival\nprobability") +
        xlab ("Months")


pfs_km$plot <- pfs_km$plot + 
        theme_minimal() +
        theme(
                panel.grid = element_blank() ,
                
                axis.line = element_line(color = "#000000" , linewidth = 0.1),
                axis.ticks = element_line(color = "#000000" ) ,
                
                axis.title = element_text ( size = 13, family = "Helvetica", face="plain" , colour = "#000000"),
                
                axis.text = element_text ( size = 13, family = "Helvetica", face="plain" ,
                                           hjust = 0.5, colour = "#000000"),
                
                legend.title = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 10  ),
                legend.text  = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 10) ,
                legend.position= "none" ,
                aspect.ratio = 0.65)


pfs_km$table <- pfs_km$table + 
        theme_minimal() +
        theme(panel.grid = element_blank() ,
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              
              axis.text.y = element_text ( size = 13, family = "Helvetica", face="plain" ,
                                           colour = "#000000"),
              axis.text.x = element_blank(),
              legend.position= "none" ,
              aspect.ratio = 0.176)

print(pfs_km)
