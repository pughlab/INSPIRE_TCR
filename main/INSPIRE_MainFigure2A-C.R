library(tidyverse)
library (ggbeeswarm)

### Input data preparation ----------------------------------------------------------------------------------------
data_path <- "https://github.com/pughlab/INSPIRE_TCR/blob/main/Data"

### Clinical data:
clinical_data_fname <- "INSPIRE_ClinicalData.xlsx"

clinical_data <- readxl::read_xlsx(
        file.path (data_path , clinical_data_fname))%>%
        dplyr::select(Patient_id , 
                      COHORT , 
                      `Best response longevity` )

### Metastasis site data:
MetSite_data_fname <- "Metastasis_site.xlsx"
MetSite_data <- readxl::read_xlsx(file.path (data_path , MetSite_data_fname))

### Diversity indices:
div_fname <- "inspire_tumor_diversity_inices.csv"

diversity_indices <- readr::read_csv(
        file.path (data_path , div_fname)) %>% 
        filter (Locus == "CLONES_TRB") %>%
        filter(!(Patient_id == "INS-D-006" & Cycle == "EOTT")) %>%
        filter(Order_q == 1) %>%
        dplyr::select(Patient_id , Cycle , Diversity) %>%
        left_join(clinical_data ,
                  by = "Patient_id") %>%
        left_join(
                MetSite_data ,
                by = c("Patient_id" , "Cycle"))

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
                              
cohort_color_pal = tibble(
        COHORT = c ("SCCHN", "TNBC", "HGSOC", "MM", "MST"),
        color  = c ( "#441e11", "#FDDA23",  "#FAC0C0" , "#A7A7A7", "#B0BF1A") )
                              


### Figure 2A

ggplot( data = diversity_indices %>%
                filter(Cycle == "ST") %>%
                mutate(COHORT = fct_relevel(
                        COHORT ,
                        "SCCHN", 
                        "TNBC",
                        "HGSOC",
                        "MM",
                        "MST")),
        aes(
                x = COHORT ,
                y = Diversity )) +
        
        geom_boxplot(
                outlier.shape = NA ,
                width = 0.5 ,
                fill = "transparent" ,
                linewidth = 0.25) +
        ggbeeswarm::geom_beeswarm(
                aes(fill = `Best response longevity`),
                method = "center" ,
                cex = 4 ,
                shape = 21 ,
                size = 4.5 ,
                stroke = 0.05 ) +

        scale_fill_manual(values = RECIST_color_pal$color,
                          name = "Best response",
                          breaks = RECIST_color_pal$response ,
                          guide = "none") +
        
        scale_x_discrete(
                breaks = c(
                        "SCCHN", 
                        "TNBC",
                        "HGSOC",
                        "MM",
                        "MST" ) ,
                label = c(
                        "SCCHN\nN=5" ,
                        "TNBC\nN=5" ,
                        "HGSOC\nN=5" ,
                        "MM\nN=4" ,
                        "MST\nN=14" ),
                position = "bottom" )+
        scale_y_log10(limits = c(5 , 1000)) +
        ylab ("Baseline tumour Shannon diversity")+
        
        theme_minimal() +
        theme(
                panel.grid = element_blank() ,
                
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
                legend.position='bottom' ,
                aspect.ratio = 1.25 )



### Figure 2B

ggplot( data = diversity_indices %>%
                filter(Cycle == "ST")%>%
                mutate(Tissue = fct_relevel(
                        Tissue ,
                        c("Lymph node" ,
                          "Lung" ,
                          "Liver" ,
                          "Subcutaneous tissue" ,
                          "Tongue and sublingual" ,
                          "Peritoneal mass" ,
                          "Abdominal wall" ,
                          "Gluteal mass" ,
                          "Skin" ,
                          "Adnexal mass")
                )) ,
        aes(
                x = Tissue ,
                y = Diversity )) +
        
        geom_boxplot(
                outlier.shape = NA ,
                width = 0.5 ,
                fill = "transparent" ,
                linewidth = 0.25) +
        ggbeeswarm::geom_beeswarm(
                aes(
                        fill = COHORT
                ),
                method = "center" ,
                cex = 3 ,
                shape = 21 ,
                size = 4.5 ,
                stroke = 0.05 ) +
        
        scale_fill_manual(values = cohort_color_pal$color,
                          name = "Cohort",
                          breaks = cohort_color_pal$COHORT ,
                          guide = "none") +
        
        scale_y_log10(limits = c(5 , 1000)) +
        ylab ("Baseline tumour Shannon diversity")+
        
        theme_minimal() +
        theme(
                panel.grid = element_blank() ,
                
                axis.line.y = element_line(color = "#000000" , linewidth = 0.1),
                axis.ticks = element_line(color = "#000000" , linewidth = 0.25) ,
                axis.ticks.length = unit(0.25 , "line") ,
                
                axis.title.x = element_blank(),
                axis.title.y = element_text ( size = 15, family = "Helvetica", face="plain" , colour = "#000000"),
                
                axis.text.y = element_text ( size = 15, family = "Helvetica", face="plain" ,
                                             hjust = 0.5, colour = "#000000"),
                axis.text.x = element_text ( size = 13, family = "Helvetica", face="plain" ,
                                             hjust = 1, angle = 30 , colour = "#000000" ) ,
                legend.title = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 12  ),
                legend.text  = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 12) ,
                aspect.ratio = 1/1.25 )



### Figure 2C

ggplot( data = diversity_indices %>%
                filter(Cycle != "ST") %>%
                mutate(COHORT = fct_relevel(
                        COHORT ,
                        "SCCHN", 
                        "TNBC",
                        "HGSOC",
                        "MM",
                        "MST")),
        aes(
                x = COHORT ,
                y = Diversity )) +
        
        geom_boxplot(
                outlier.shape = NA ,
                width = 0.5 ,
                fill = "transparent" ,
                linewidth = 0.25) +
        ggbeeswarm::geom_beeswarm(
                aes( fill = `Best response longevity`),
                method = "center" ,
                cex = 4 ,
                shape = 21 ,
                size = 4.5 ,
                stroke = 0.05 ) +
        
        scale_fill_manual(values = RECIST_color_pal$color,
                          name = "Best response",
                          breaks = RECIST_color_pal$response ,
                          guide = "none") +

        
        scale_x_discrete(
                breaks = c(
                        "SCCHN", 
                        "TNBC",
                        "HGSOC",
                        "MM",
                        "MST" ) ,
                label = c(
                        "SCCHN\nN=2" ,
                        "TNBC\nN=2" ,
                        "HGSOC\nN=5" ,
                        "MM\nN=4" ,
                        "MST\nN=11" ),
                position = "bottom" )+
        scale_y_log10(limits = c(5 , 1000)) +
        ylab ( expression ("Cycle"[2/3]* " on-ICB tumour Shannon diversity") ) +
        
        theme_minimal() +
        theme(
                panel.grid = element_blank() ,
                
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
                legend.position='bottom' ,
                aspect.ratio = 1.25 )




