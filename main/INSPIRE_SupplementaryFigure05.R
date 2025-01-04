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
div_fname <- "INSPIRE_Tumour_DiversityIndices.csv"

diversity_indices <- readr::read_csv(
        file.path (data_path , div_fname)) %>% 
        filter (Locus == "CLONES_TRB") %>%
        filter(!(Patient_id == "INS-D-006" & Cycle == "EOTT")) %>%
        dplyr::select(Patient_id , Cycle , Order_q , Diversity) %>%
        filter(Order_q <= 127) %>%
        left_join(clinical_data ,
                  by = "Patient_id") %>%
        mutate(Cycle = ifelse(Cycle == "ST" , "ST" , "on-ICB")) %>%
        mutate(Cycle = fct_relevel(Cycle , "ST" , "on-ICB") ,
               COHORT = fct_relevel(COHORT , "SCCHN", "TNBC", "HGSOC", "MM", "MST"))

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



### Figure S5


ggplot(data = diversity_indices) +
        geom_point(
                aes(
                        x = Order_q ,
                        y = Diversity ,
                        fill = `Best response longevity`
                ) ,
                shape = 21 ,
                stroke = 0.01,
                size = 1.5) +

        geom_line(
                aes(x = Order_q ,
                    y = Diversity ,
                    group = Patient_id,
                    color = `Best response longevity`) ,
                linewidth = 0.1)+
        ggrepel::geom_text_repel(data = diversity_indices %>%
                                         filter(Order_q == 0) ,
                                 aes(x = Order_q ,
                                     y = Diversity ,
                                     label = Patient_id,
                                     color = `Best response longevity`),
                                 hjust = 3,
                                 direction = "y",
                                 size = 2.5,
                                 segment.size = .5,
                                 segment.linetype = "dotted",
                                 #segment.square = TRUE,
                                 #Always show all labels, even when they have too many overlaps:
                                 max.overlaps = Inf,
                                 #Do not repel labels from data points:
                                 point.size = NA,
                                 #Use min.segment.length = 0 to draw all line segments, no matter how short they are and min.segment.length = Inf to never show the lines:
                                 min.segment.length = Inf,
                                 box.padding = 0.2,
                                 nudge_x = 0.2,
                                 show.legend = FALSE)+
        ylab ("Hill Number") +
        scale_fill_manual(values = RECIST_color_pal$color,
                          name = "Best response",
                          breaks = RECIST_color_pal$response ) +
        scale_color_manual(values = RECIST_color_pal$color,
                           name = "Best response",
                           breaks = RECIST_color_pal$response ) +
        scale_x_continuous(trans = scales::pseudo_log_trans(base = 2) ,
                           breaks = c(0, 1, 2 , 127) ,
                           labels = c("Richness" ,
                                      "Shannon" ,
                                      "Simpson" ,
                                      "Effective number of\nexpanded clonotypes"),
                           guide = guide_axis(n.dodge = 1),
                           expand = c(0.3 , 0.2)) +
        theme_minimal() +
        theme(
                panel.grid = element_blank() ,
                
                axis.line.y = element_line(color = "#000000" , linewidth = 0.1),
                axis.ticks = element_line(color = "#000000" , linewidth = 0.25) ,
                axis.ticks.length = unit(0.25 , "line") ,
                
                axis.title.x = element_blank(),
                axis.title.y = element_text ( size = 7.5, family = "Helvetica", face="plain" , colour = "#000000"),
                
                axis.text.y = element_text ( size = 7.5, family = "Helvetica", face="plain" ,
                                             hjust = 0.5, colour = "#000000"),
                axis.text.x = element_text ( size = 7.5, family = "Helvetica", face="plain" ,
                                             hjust = 1, angle = 30 , colour = "#000000" ) ,
                legend.title = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 7.5  ),
                legend.text  = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 7.5)  ,
                legend.position = "bottom" ,
                
                strip.text.x = element_text ( size = 7.5, family = "Helvetica", face="plain" ,
                                              hjust = 0.5, colour = "#000000" ) ,
                strip.text.y = element_text ( size = 7.5, family = "Helvetica", face="plain" ,
                                              hjust = 0.5, angle = 0 , colour = "#000000" ))+
        facet_grid(Cycle ~ COHORT)
