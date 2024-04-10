library(tidyverse)

### Input data preparation ----------------------------------------------------------------------------------------
data_path <- "https://raw.githubusercontent.com/pughlab/INSPIRE_TCR/main/Data"

### Clinical data:
clinical_data_fname <- "INSPIRE_ClinicalData.csv"

clinical_data <- readr::read_csv(
        file.path (data_path , clinical_data_fname))%>%
        dplyr::select(Patient_id , COHORT , `Best response longevity`)



### Longitudinal diversity FMM clusters:
FMM_data_fname <- "INSPIRE_PBMC_LongitudinalDiversityFMMClusters.csv"

FMM_data <- readr::read_csv(
        file.path (data_path , FMM_data_fname))%>%
        dplyr::select(Patient_id ,  `FMM Cluster`)

### Diversity indices:
div_fname <- "INSPIRE_PBMC_DiversityIndices.csv"

diversity_indices <- readr::read_csv(file.path (data_path , div_fname) ) %>% 
        filter (Locus == "CLONES_TRB") %>%
        filter(Order_q == 1 ) %>%
        dplyr::select(Patient_id , Cycle , Cycle_Code ,  Diversity) %>%
        left_join(clinical_data ,
                  by = "Patient_id")%>%
        mutate( COHORT = fct_relevel( COHORT , 
                                      "SCCHN", 
                                      "TNBC",
                                      "HGSOC",
                                      "MM", 
                                      "MST" )) %>%
        mutate( Cycle = fct_relevel(Cycle ,
                                    "SB" , "C3B" , "C6B" , "C9B" ,
                                    "C12B" , "C15B" , "C18B" , "C21B" ,
                                    "C24B" , "C27B" , "C30B" , "C33B" ,  
                                    "EOTB"))

### Patients with PBMC samples available at 5 or more timepoints:
Patient_with_continued_benefit <- (diversity_indices %>%
                                           group_by(Patient_id) %>%
                                           summarise(n = n()) %>%
                                           ungroup () %>%
                                           arrange(desc (n)) %>%
                                           filter(n > 4) )$Patient_id


diversity_indices <- diversity_indices %>%
        filter(Patient_id %in% Patient_with_continued_benefit) %>%
        left_join(FMM_data ,
                  by = "Patient_id")


### Visualization -------------------------------------------------------------------------------------------------
### Color palettes:
RECIST_color_pal <- tibble(response = 
                                   c ("CR", "PR",
                                      "SD, 6 < n ICB cycles",
                                      "SD, n < 6 ICB cycles",
                                      "PD" , "NE"),
                           color = 
                                   c ("#033483" , "#A7D2E9" ,
                                      "#7DB290" ,
                                      "#FEA500" ,
                                      "#C5231B" , "#AFAFAF" ))



### Figure 1E

ggplot(data = diversity_indices ) +
        geom_line(
                aes(
                        x = Cycle_Code * 3 ,
                        y = Diversity ,
                        group = Patient_id ,
                        color =  `Best response longevity`
                ),
                linewidth = 0.2
        ) +
        scale_color_manual(
                values = RECIST_color_pal$color ,
                breaks = RECIST_color_pal$response
        ) +
        ylab ("PBMC Shannon diversity") +
        xlab ("Time on-ICB (weeks)") +
        theme_minimal() +
        
        theme(
                panel.grid = element_blank(),
                panel.spacing = unit(1 , units = "line"),
                
                axis.line = element_line(colour = "#000000" , linewidth = 0.1),
                axis.ticks = element_line(color = "#000000" , linewidth = 0.25) ,
                axis.ticks.length = unit(0.25 , "line") ,
                
                axis.title = element_text ( vjust = 2, family = "Helvetica", face="plain" ,
                                            size = 15 , hjust = 0.5 ),
                
                
                axis.text = element_text ( size = 15, family = "Helvetica", face="plain" ,
                                           hjust = 0.5, colour = "#000000"),
                
                
                strip.text = element_text ( size = 15, family = "Helvetica", face="plain" ,
                                            hjust = 0, colour = "#000000"),
                
                legend.position = "none" ,
                
                
                aspect.ratio = 1/1.6 
        )+
        facet_grid(. ~ `FMM Cluster` ,
                   labeller = as_labeller(c(`Cluster I` = "Stable" ,
                                            `Cluster II` = "Sporadic fluctuations to\nstabilization transition" ,
                                            `Cluster III` = "Cyclic fluctuations")))
