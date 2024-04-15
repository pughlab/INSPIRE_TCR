library(tidyverse)
library(ggdist)

### Input data preparation ----------------------------------------------------------------------------------------
data_path <- "https://raw.githubusercontent.com/pughlab/INSPIRE_TCR/main/Data"



### cfDNA clonotypes:
cfdna_fname <- "INSPIRE_Plasma_TRBClonotypes.csv"

cfdna <- readr::read_csv( file.path(data_path , cfdna_fname)) %>%
        dplyr::select(Patient_id , Cycle , aaSeqCDR3) %>%
        unique()

cfdna_counts <- cfdna %>%
        group_by (Patient_id , Cycle) %>%
        summarise(Count = n()) %>%
        ungroup()


### PBMC clonotypes:
pbmc_fname <- "INSPIRE_PBMC_TRBClonotypes.csv"

pbmc <- readr::read_csv(file.path(data_path , pbmc_fname)) %>%
        dplyr::select(Patient_id , Cycle , aaSeqCDR3) %>%
        unique()

pbmc_counts <- pbmc %>%
        group_by (Patient_id , Cycle) %>%
        summarise(Count = n()) %>%
        ungroup()


### Tumour clonotypes:
tumour_fname <- "INSPIRE_Tumour_TRBClonotypes.csv"

tumour <- readr::read_csv(file.path(tumour_path , tumour_fname)) %>%
        filter(!(Patient_id == "INS-D-006" & Cycle == "EOTT")) %>%
        filter(! (Patient_id %in% c("INS-B-025" , "INS-E-028_L" , "INS-E-028_R") ) ) %>%
        dplyr::select(Patient_id , Cycle , aaSeqCDR3) %>%
        unique()

BaselineTumour_PatientList <- levels(as.factor( (tumour %>% filter (Cycle == "ST"))$Patient_id ))
OnICBTumour_PatientList <- levels(as.factor( (tumour %>% filter (Cycle != "ST"))$Patient_id ))

### ---------------------------------------------------------------------------------------------------------
### Calculating the commonality of clonotypes between tumour and PBMC samples: ------------------------------


BaselineTumour_PBMCCommonality <- tibble( Patient_id = BaselineTumour_PatientList )

BaselineTumour_PBMCCommonality [ , c("SB" , 
                                     "C3B" , "C6B" , "C9B" , 
                                     "C12B" , "C15B" , "C18B" ,
                                     "C21B" , "C24B" , "C27B" ,
                                     "C30B" , "C33B", "EOTB")] <- NA


for (i in BaselineTumour_PatientList) {
        
        tumourDF <- tumour %>%
                filter(Cycle == "ST" & Patient_id == i)
        
        for (j in colnames(BaselineTumour_PBMCCommonality) [-1]) {
                
                pbmcDF <- pbmc %>%
                        filter(Cycle == j & Patient_id == i)
                
                if (nrow (pbmcDF) > 0 ) {
                        
                        BaselineTumour_PBMCCommonality [[j]][BaselineTumour_PBMCCommonality$Patient_id == i] <- length(intersect(tumourDF$aaSeqCDR3 ,
                                                                                                                                 pbmcDF$aaSeqCDR3))
                        
                }
                
                rm (pbmcDF)
                
        }
        rm (tumourDF)
        
}



#------------------------------------------------------------

OnICBTumour_PBMCCommonality <- tibble( Patient_id = OnICBTumour_PatientList )

OnICBTumour_PBMCCommonality [ , c("SB" , 
                                  "C3B" , "C6B" , "C9B" , 
                                  "C12B" , "C15B" , "C18B" ,
                                  "C21B" , "C24B" , "C27B" ,
                                  "C30B" , "C33B", "EOTB")] <- NA


for (i in OnICBTumour_PatientList) {
        
        tumourDF <- tumour %>%
                filter(Cycle != "ST" & Patient_id == i)
        
        for (j in colnames(OnICBTumour_PBMCCommonality) [-1]) {
                
                pbmcDF <- pbmc %>%
                        filter(Cycle == j & Patient_id == i)
                
                if (nrow (pbmcDF) > 0 ) {
                        
                        OnICBTumour_PBMCCommonality [[j]][OnICBTumour_PBMCCommonality$Patient_id == i] <- length(intersect(tumourDF$aaSeqCDR3 ,
                                                                                                                           pbmcDF$aaSeqCDR3))
                        
                }
                
                rm (pbmcDF)
                
        }
        rm (tumourDF)
        
}


### ---------------------------------------------------------------------------------------------------------
### Calculating the commonality of clonotypes between tumour and Plasma samples: ----------------------------

BaselineTumour_PlasmaCommonality <- tibble( Patient_id = BaselineTumour_PatientList )

BaselineTumour_PlasmaCommonality [ , c("SB" , 
                                       "C3B" , "C6B" , "C9B" , 
                                       "C12B" , "C15B" , "EOTB")] <- NA


for (i in BaselineTumour_PatientList) {
        
        tumourDF <- tumour %>%
                filter(Cycle == "ST" & Patient_id == i)
        
        for (j in colnames(BaselineTumour_PlasmaCommonality) [-1]) {
                
                PlasmaDF <- cfdna %>%
                        filter(Cycle == j & Patient_id == i)
                
                if (nrow (PlasmaDF) > 0 ) {
                        
                        BaselineTumour_PlasmaCommonality [[j]][BaselineTumour_PlasmaCommonality$Patient_id == i] <- length(intersect(tumourDF$aaSeqCDR3 ,
                                                                                                                                     PlasmaDF$aaSeqCDR3))
                        
                }
                
                rm (PlasmaDF)
                
        }
        rm (tumourDF)
        
}


#------------------------------------------------------------


OnICBTumour_PlasmaCommonality <- tibble( Patient_id = OnICBTumour_PatientList )

OnICBTumour_PlasmaCommonality [ , c("SB" , 
                                    "C3B" , "C6B" , "C9B" , 
                                    "C12B" , "C15B" , "EOTB")] <- NA


for (i in OnICBTumour_PatientList) {
        
        tumourDF <- tumour %>%
                filter(Cycle != "ST" & Patient_id == i)
        
        for (j in colnames(OnICBTumour_PlasmaCommonality) [-1]) {
                
                PlasmaDF <- cfdna %>%
                        filter(Cycle == j & Patient_id == i)
                
                if (nrow (PlasmaDF) > 0 ) {
                        
                        OnICBTumour_PlasmaCommonality [[j]][OnICBTumour_PlasmaCommonality$Patient_id == i] <- length(intersect(tumourDF$aaSeqCDR3 ,
                                                                                                                               PlasmaDF$aaSeqCDR3))
                        
                }
                
                rm (PlasmaDF)
                
        }
        rm (tumourDF)
        
}


### ---------------------------------------------------------------------------------------------------------
### Constructing a dataframe of all commonolaties among tumour, pbmc, and plasma for downstream visualization: 

LocalSystemic <- 
        rbind(
                BaselineTumour_PBMCCommonality %>%
                        pivot_longer(cols = c(2:ncol(BaselineTumour_PBMCCommonality)) ,
                                     names_to = "Cycle" ,
                                     values_to = "Commonality") %>%
                        mutate(`Systemic Repertoire` = "PBMC" ,
                               `Tumour Timepoint` = "Baseline")%>%
                        filter(!is.na(Commonality)) %>%
                        left_join(pbmc_counts , by = c("Patient_id" , "Cycle")),
                
                
                OnICBTumour_PBMCCommonality %>%
                        pivot_longer(cols = c(2:ncol(OnICBTumour_PBMCCommonality)) ,
                                     names_to = "Cycle" ,
                                     values_to = "Commonality") %>%
                        mutate(`Systemic Repertoire` = "PBMC" ,
                               `Tumour Timepoint` = "On-ICB")%>%
                        filter(!is.na(Commonality)) %>%
                        left_join(pbmc_counts , by = c("Patient_id" , "Cycle")),
                
                
                BaselineTumour_PlasmaCommonality %>%
                        pivot_longer(cols = c(2:ncol(BaselineTumour_PlasmaCommonality)) ,
                                     names_to = "Cycle" ,
                                     values_to = "Commonality") %>%
                        mutate(`Systemic Repertoire` = "Plasma" ,
                               `Tumour Timepoint` = "Baseline")%>%
                        filter(!is.na(Commonality))%>%
                        left_join(cfdna_counts , by = c("Patient_id" , "Cycle")) ,
                
                
                OnICBTumour_PlasmaCommonality %>%
                        pivot_longer(cols = c(2:ncol(OnICBTumour_PlasmaCommonality)) ,
                                     names_to = "Cycle" ,
                                     values_to = "Commonality") %>%
                        mutate(`Systemic Repertoire` = "Plasma" ,
                               `Tumour Timepoint` = "On-ICB")%>%
                        filter(!is.na(Commonality))%>%
                        left_join(cfdna_counts , by = c("Patient_id" , "Cycle")) ) %>%
        filter(Cycle %in% c("SB" , "C3B" , "C6B" )) %>%
        mutate(Cycle = fct_relevel(Cycle ,
                                   c("SB" , "C3B" , "C6B" )))

### ---------------------------------------------------------------------------------------------------------
### Quick statistics of the commonalities: ------------------------------------------------------------------

LocalSystemic %>%
        filter(`Tumour Timepoint` == "Baseline") %>%
        filter(Cycle == "SB") %>%
        mutate(Fraction = Commonality / Count) %>%
        group_by(`Systemic Repertoire`) %>%
        summarise(Median = median(Fraction * 100))

# `Systemic Repertoire`         Median
# 1 PBMC                        1.72
# 2 Plasma                      7.11


LocalSystemic %>%
        filter(`Tumour Timepoint` == "On-ICB" ) %>%
        filter(Cycle == "C3B") %>%
        mutate(Fraction = Commonality / Count) %>%
        group_by(`Systemic Repertoire`) %>%
        summarise(Median = median(Fraction * 100))

# `Systemic Repertoire`         Median
# 1 PBMC                        1.85
# 2 Plasma                      16.2 

### ---------------------------------------------------------------------------------------------------------
### Visualization: ------------------------------------------------------------------------------------------


### Figure 4A

ggplot(data = LocalSystemic ,
       aes(
               x = Cycle ,
               y = 100 * Commonality / Count )) +
        ggdist::stat_pointinterval(
                aes( color = `Systemic Repertoire` ,
                     side = `Systemic Repertoire`),
                position = position_dodge(width = - 0.2) ,
                .width = c(0.125, 0.4),
                fatten_point = 3 ,
                
                point_color = "#000000" ,
                point_fill = "#000000" ,
                point_alpha = 0.5 ,
                shape = 16
        ) +
        
        ggdist::geom_weave(
                aes( fill = `Systemic Repertoire` ,
                     side = `Systemic Repertoire`)  ,
                slab_linewidth = 0.1 ,
                slab_color = "#000000" ,
                #stackratio = 0.5 ,
                position = position_dodge(width = - 0.5) ,
                dotsize = 2.25 ,
                #alpha = 0.75
        ) +
        # stat_spike(aes( color = `Systemic Repertoire` ) , 
        #            at = median) +
        
        scale_side_mirrored(guide = "none") +
        
        scale_color_manual( values = c("#FFD403" , "#A12A2A") ,
                            breaks = c("Plasma" , "PBMC")) +
        scale_fill_manual( values = c("#FFD403" , "#A12A2A") ,
                           breaks = c("Plasma" , "PBMC")) +
        scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), 
                           breaks = c(1, 10 , 100) )+
        ylab ("Percentage of systemic TRB CDR3s\noverlapping with tumour")+
        theme_minimal() +
        theme(
                panel.grid = element_blank() ,
                panel.spacing = unit(2 , units = "line") ,
                
                axis.line.y = element_line(color = "#000000" , linewidth = 0.1),
                axis.ticks = element_line(color = "#000000" , linewidth = 0.25) ,
                axis.ticks.length = unit(0.25 , "line") ,
                
                axis.title.x = element_blank(),
                axis.title.y = element_text ( size = 15, family = "Helvetica", face="plain" , colour = "#000000"),
                
                axis.text.y = element_text ( size = 15, family = "Helvetica", face="plain" ,
                                             hjust = 0.5, colour = "#000000"),
                axis.text.x = element_text ( size = 13, family = "Helvetica", face="plain" ,
                                             hjust = 0.5, colour = "#000000" ) ,
                legend.title = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 12  ),
                legend.text  = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 12) ,
                
                strip.text.y  = element_text(family = "Helvetica", face="plain" ,
                                             color = "#000000", size = 12 , angle = 0) ,
                aspect.ratio = 0.3
        ) +
        facet_grid(`Tumour Timepoint` ~ .)






### Figure 4B

ggplot(data = LocalSystemic ,
       aes(
               x = Cycle ,
               y = Commonality  )) +
        ggdist::stat_pointinterval(
                aes( color = `Systemic Repertoire` ,
                     side = `Systemic Repertoire`),
                position = position_dodge(width = - 0.2) ,
                .width = c(0.125, 0.4),
                fatten_point = 3 ,
                
                point_color = "#000000" ,
                point_fill = "#000000" ,
                point_alpha = 0.5 ,
                shape = 16
        ) +
        
        ggdist::geom_weave(
                aes( fill = `Systemic Repertoire` ,
                     side = `Systemic Repertoire`)  ,
                slab_linewidth = 0.1 ,
                slab_color = "#000000" ,
                #stackratio = 0.5 ,
                position = position_dodge(width = - 0.5) ,
                dotsize = 2.25 ,
                #alpha = 0.75
        ) +
        # stat_spike(aes( color = `Systemic Repertoire` ) , 
        #            at = median) +
        
        scale_side_mirrored(guide = "none") +
        
        scale_color_manual( values = c("#FFD403" , "#A12A2A") ,
                            breaks = c("Plasma" , "PBMC")) +
        scale_fill_manual( values = c("#FFD403" , "#A12A2A") ,
                           breaks = c("Plasma" , "PBMC")) +
        scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), 
                           breaks = c(1, 10 , 100) )+
        ylab ("Absolute count of systemic TRB CDR3s\noverlapping with tumour")+
        theme_minimal() +
        theme(
                panel.grid = element_blank() ,
                panel.spacing = unit(2 , units = "line") ,
                
                axis.line.y = element_line(color = "#000000" , linewidth = 0.1),
                axis.ticks = element_line(color = "#000000" , linewidth = 0.25) ,
                axis.ticks.length = unit(0.25 , "line") ,
                
                axis.title.x = element_blank(),
                axis.title.y = element_text ( size = 15, family = "Helvetica", face="plain" , colour = "#000000"),
                
                axis.text.y = element_text ( size = 15, family = "Helvetica", face="plain" ,
                                             hjust = 0.5, colour = "#000000"),
                axis.text.x = element_text ( size = 13, family = "Helvetica", face="plain" ,
                                             hjust = 0.5, colour = "#000000" ) ,
                legend.title = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 12  ),
                legend.text  = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 12) ,
                
                strip.text.y  = element_text(family = "Helvetica", face="plain" ,
                                             color = "#000000", size = 12 , angle = 0) ,
                aspect.ratio = 0.3
        ) +
        facet_grid(`Tumour Timepoint` ~ .)
