library(tidyverse)
library (ggbeeswarm)

### Input data preparation ----------------------------------------------------------------------------------------
data_path <- "https://raw.githubusercontent.com/pughlab/INSPIRE_TCR/main/Data"

### Clinical data:
clinical_data_fname <- "INSPIRE_ClinicalData.csv"

clinical_data <- readr::read_csv(
        file.path (data_path , clinical_data_fname))%>%
        dplyr::select(Patient_id , 
                      COHORT , 
                      `Best response longevity` )



### Diversity indices:
div_fname <- "INSPIRE_PBMC_DiversityIndices.csv"

diversity_indices <- readr::read_csv(
        file.path (data_path , div_fname)) %>% 
        filter (Locus == "CLONES_TRB") %>%
        filter(Order_q == 1) %>%
        dplyr::select(Patient_id , Cycle , Diversity) %>%
        left_join(clinical_data ,
                  by = "Patient_id") 

### Statistical tests ---------------------------------------------------------------------------------------------

kruskal.test( Diversity ~ COHORT, 
              data = diversity_indices %>%
                      filter(Cycle == "SB"))

pairwise.wilcox.test (
        (diversity_indices %>%
                 filter(Cycle == "SB"))$Diversity ,
        (diversity_indices %>%
                 filter(Cycle == "SB"))$COHORT ,
        p.adjust.method = "BH")

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


### Figure 1B

#Quick look at the sample counts:
diversity_indices %>%
        filter(Cycle == "SB") %>%
        group_by(COHORT) %>%
        summarise(n = n())



ggplot( data = diversity_indices %>%
                filter(Cycle == "SB") %>%
                mutate(COHORT = fct_relevel(
                        COHORT ,
                        "SCCHN", 
                        "TNBC",
                        "HGSOC",
                        "MM",
                        "MST")),
        aes(
                x = COHORT ,
                y = Diversity
        )) +
        
        geom_boxplot(
                outlier.shape = NA ,
                width = 0.5 ,
                fill = "transparent" ,
                linewidth = 0.25) +
        ggbeeswarm::geom_beeswarm(
                aes(
                        fill = `Best response longevity`
                ),
                method = "center" ,
                cex = 2.5 ,
                shape = 21 ,
                size = 4.5 ,
                stroke = 0.05 ) +
        ggpubr::stat_compare_means( 
                ref.group = ".all." ,
                hide.ns = TRUE ,
                label = "p.signif" ,
                method="wilcox.test" ,
                size = 8)+
        
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
                        "SCCHN\nN=14" ,
                        "TNBC\nN=12" ,
                        "HGSOC\nN=10" ,
                        "MM\nN=10" ,
                        "MST\nN=22" ),
                position = "bottom" )+
        ylim(c(40 , 3950)) +
        ylab ("Baseline PBMC Shannon diversity")+
        
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
