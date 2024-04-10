library(tidyverse)

### Input data preparation ----------------------------------------------------------------------------------------
data_path <- "https://raw.githubusercontent.com/pughlab/INSPIRE_TCR/main/Data"

### Clinical data:
clinical_data_fname <- "INSPIRE_ClinicalData.csv"

clinical_data <- readr::read_csv(
        file.path (data_path , clinical_data_fname))%>%
        dplyr::select(Patient_id , 
                      COHORT , 
                      `Best response longevity` )



### Diversity FMM clusters:
FMM_data_fname <- "INSPIRE_PBMC_DiversityFMMClusters.csv"

FMM_data <- readr::read_csv(
        file.path (data_path , FMM_data_fname))%>%
        dplyr::select(Patient_id ,  `FMM Cluster`)

### Diversity indices:
div_fname <- "INSPIRE_PBMC_DiversityIndices.csv"

diversity_indices <- readr::read_csv(
        file.path (data_path , div_fname)) %>% 
        filter (Locus == "CLONES_TRB") %>%
        filter(Cycle %in% c ("SB" , "C3B")) %>%
        filter(Order_q == 1 ) %>%
        dplyr::select(Patient_id , Cycle , Diversity) %>%
        pivot_wider(
                names_from = Cycle ,
                values_from = Diversity
        ) %>%
        na.exclude() %>%
        left_join(clinical_data ,
                  by = "Patient_id")%>%
        left_join(FMM_data ,
                  by = "Patient_id")%>% 
        mutate( COHORT = fct_relevel( COHORT , 
                                      "SCCHN", 
                                      "TNBC",
                                      "HGSOC",
                                      "MM", 
                                      "MST" ))
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



### Figure 1C

ggplot(data = diversity_indices) +

### Baseline-Low, Diversified contours --------------------------------------------

        ggdensity::geom_hdr_lines( data = diversity_indices %>%
                                           filter(`FMM Cluster` == "Baseline-Low, Diversified" ) %>%
                                           dplyr::select(Patient_id ,
                                                         SB , C3B ,
                                                         `FMM Cluster`),
                                   aes(
                                           x = SB ,
                                           y = log2 (C3B / SB) ,
                                           color = `FMM Cluster` ),
                                   linewidth = 0.1 ,
                                   alpha = 1,
                                   #h = c(360 , 1.3) ,
                                   probs = seq(from = 0.01 ,
                                               to = 0.20 ,
                                               by = 0.025) ,
                                   linejoin = "round",
                                   lineend = "round") +
        
        ggdensity::geom_hdr_lines( data = diversity_indices %>%
                                           filter(`FMM Cluster` == "Baseline-Low, Diversified" ) %>%
                                           dplyr::select(Patient_id ,
                                                         SB , C3B ,
                                                         `FMM Cluster`),
                                   aes(
                                           x = SB ,
                                           y = log2 (C3B / SB) ,
                                           color = `FMM Cluster` ),
                                   linewidth = 0.1 ,
                                   alpha = 1,
                                   ylim = c(0.1 , 2) ,
                                   xlim = c(0 , 870) ,
                                   #h = c(360 , 1.3) ,
                                   probs = seq(from = 0.30 ,
                                               to = 0.99 ,
                                               by = 0.1) ,
                                   linejoin = "round",
                                   lineend = "round") +

### Baseline-Low, Stasis/Decline contours -----------------------------------------

        ggdensity::geom_hdr_lines( data = diversity_indices %>%
                                           filter(`FMM Cluster` == "Baseline-Low, Stasis/Decline" ) %>%
                                           dplyr::select(Patient_id ,
                                                         SB , C3B ,
                                                         `FMM Cluster`),
                                   aes(
                                           x = SB ,
                                           y = log2 (C3B / SB) ,
                                           color = `FMM Cluster` ),
                                   linewidth = 0.1 ,
                                   alpha = 1,
                                   ylim = c(-3 , 0.05) ,
                                   xlim = c(0 , 850) ,
                                   probs = seq(from = 0.01 , to = 0.99 , by = 0.05) ,
                                   linejoin = "round",
                                   lineend = "round") +
### Baseline-High, Stasis/Diversified contours ------------------------------------


        ggdensity::geom_hdr_lines( data = diversity_indices %>%
                                           filter(`FMM Cluster` == "Baseline-High, Stasis/Diversified" )%>%
                                           dplyr::select(Patient_id ,
                                                         SB , C3B ,
                                                         `FMM Cluster`),
                                   aes(
                                           x = SB ,
                                           y = log2 (C3B / SB) ,
                                           color = `FMM Cluster` ),
                                   linewidth = 0.1 ,
                                   alpha = 1,
                                   #h = c(360 , 1.3) ,
                                   probs = seq(from = 0.01 , to = 0.3 , by = 0.025) ,
                                   linejoin = "round",
                                   lineend = "round") +
        ggdensity::geom_hdr_lines( data = diversity_indices %>%
                                           filter(`FMM Cluster` == "Baseline-High, Stasis/Diversified" )%>%
                                           dplyr::select(Patient_id ,
                                                         SB , C3B ,
                                                         `FMM Cluster`),
                                   aes(
                                           x = SB ,
                                           y = log2 (C3B / SB) ,
                                           color = `FMM Cluster` ),
                                   linewidth = 0.1 ,
                                   alpha = 1,
                                   ylim = c(-0.25 , 3) ,
                                   xlim = c(860 , 3300) ,
                                   probs = seq(from = 0.3 , to = 0.99 , by = 0.1) ,
                                   linejoin = "round",
                                   lineend = "round") +
        
### Baseline-High, Diminution contours --------------------------------------------

ggdensity::geom_hdr_lines( data = diversity_indices %>%
                                   filter(`FMM Cluster` == "Baseline-High, Diminution" )%>%
                                   dplyr::select(Patient_id ,
                                                 SB , C3B ,
                                                 `FMM Cluster`),
                           aes(
                                   x = SB ,
                                   y = log2 (C3B / SB) ,
                                   color = `FMM Cluster` ),
                           linewidth = 0.1 ,
                           alpha = 1,
                           ylim = c(-3 , -0.1) ,
                           xlim = c(900 , 3300) ,
                           probs = seq(from = 0.01 , to = 0.99 , by = 0.1) ,
                           linejoin = "round",
                           lineend = "round") +
        
### Adding the scatter plot -------------------------------------------------------

        geom_hline(
                yintercept = 0,
                color = "#000000" ,
                linewidth = 0.1 ,
                linetype = "dashed")+

        
        
        geom_point(
                aes(
                        x = SB ,
                        y = log2 (C3B / SB) ,
                        fill = `Best response longevity`  ),
                size = 4.5 ,
                shape = 21 ,
                stroke = 0.1 ) +
        scale_fill_manual(values = RECIST_color_pal$color,
                          name = "Best response",
                          breaks = RECIST_color_pal$response ,
                          guid = "none")+
        scale_color_manual(
                values = c("#016401"  , #"Baseline-High, Stasis/Diversified" ,
                           "#C5C8C7" , #"Baseline-High, Diminution" ,         
                           "#CD4E00", #"Baseline-Low, Diversified" ,        
                           "#3B3029"   #"Baseline-Low, Stasis/Decline"
                ),
                name = "Diversity shift",
                breaks = c( "Baseline-High, Stasis/Diversified" ,
                            "Baseline-High, Diminution" ,         
                            "Baseline-Low, Diversified" ,        
                            "Baseline-Low, Stasis/Decline" ) ) +
        guides(color=guide_legend(ncol=1 ,
                                  override.aes = list(linewidth = 0.5 )))+
        ylab (expression(paste("log ( " ,
                               frac("Cycle"[3]*" on-ICB PBMC Shannon diversity", 
                                    "Baseline PBMC Shannon diversity") ,
                               " )" , sep = "")
        ))+
        xlab ("Baseline PBMC Shannon diversity")+
        facet_wrap( COHORT ~ . )  +
        theme_minimal()+
        theme (
                panel.grid = element_blank (),
                panel.spacing = unit(1 , units = "line"),
                
                axis.line.y = element_line(color = "#000000" , linewidth = 0.1),
                axis.ticks = element_line(color = "#000000" , linewidth = 0.25) ,
                axis.ticks.length = unit(0.25 , "line") ,
                
                axis.title = element_text ( vjust = 2, family = "Helvetica", face="plain" ,
                                            size = 15 , hjust = 0.5 ),
                
                
                axis.text = element_text ( size = 15, family = "Helvetica", face="plain" ,
                                           hjust = 0.5, colour = "#000000"),
                
                legend.title = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 12  ),
                legend.text  = element_text(family = "Helvetica", face="plain" ,
                                            color = "#000000", size = 12) ,
                
                strip.text = element_text (family = "Helvetica", face="plain" ,
                                           size = 15 , color = "#000000"),
                
                legend.position = "none" ,
                
                
                aspect.ratio = 1/2.1) 
