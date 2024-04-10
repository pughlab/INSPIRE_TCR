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

### Metastasis site data:
MetSite_data_fname <- "Metastasis_site.csv"
MetSite_data <- readr::read_csv(file.path (data_path , MetSite_data_fname))%>%
        dplyr::select(Patient_id , Tissue) %>%
        unique()

### Diversity FMM clusters:
FMM_data_fname <- "INSPIRE_Tumour_DiversityFMMClusters.csv"

FMM_data <- readr::read_csv(
        file.path (data_path , FMM_data_fname))%>%
        dplyr::select(Patient_id ,  `FMM Cluster`)

### Diversity indices:
div_fname <- "INSPIRE_Tumour_DiversityInices.csv"

diversity_indices <- readr::read_csv(
        file.path (data_path , div_fname)) %>% 
        filter (Locus == "CLONES_TRB") %>%
        filter(!(Patient_id == "INS-D-006" & Cycle == "EOTT")) %>%
        filter(Order_q == 1 ) %>%
        dplyr::select(Patient_id , Cycle , Diversity) %>%
        mutate(Cycle = ifelse(Cycle == "ST" , "PreICB" , "OnICB"))%>%
        pivot_wider(
                names_from = Cycle ,
                values_from = Diversity
        ) %>%
        na.exclude() %>%
        left_join(clinical_data ,
                  by = "Patient_id")%>%
        left_join(FMM_data ,
                  by = "Patient_id") %>%
        left_join(MetSite_data , by = "Patient_id" )

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
                              


### Figure 2D



ggplot(
        data = diversity_indices ) +
        geom_hline(
                yintercept = 0,
                color = "#000000" ,
                linewidth = 0.1 ,
                linetype = "dashed")+
        stat_density_2d( data = diversity_indices %>%
                                 filter(`FMM Cluster` == "Stagnant diversity"),
                         aes(
                                 x = PreICB ,
                                 y = log2 (OnICB / PreICB) ,
                                 color = `FMM Cluster` ),
                         linewidth = 0.1 ,
                         h = c(0.5 , 1) ,
                         bins = 15 ,
                         linejoin = "round",
                         lineend = "round") +
        
        stat_density_2d( data = diversity_indices %>%
                                 filter(`FMM Cluster` %in% c( "Persistently low diversity")),
                         aes(
                                 x = PreICB ,
                                 y = log2 (OnICB / PreICB) ,
                                 color = `FMM Cluster` ),
                         linewidth = 0.1 ,
                         bins = 10 ,
                         linejoin = "round" ,
                         lineend = "round") +
        stat_density_2d( data = diversity_indices %>%
                                 filter(`FMM Cluster` %in% c( "Diversified")),
                         aes(
                                 x = PreICB ,
                                 y = log2 (OnICB / PreICB) ,
                                 color = `FMM Cluster` ),
                         linewidth = 0.1 ,
                         h = c(0.2 , 0.5) ,
                         bins = 15 ,
                         linejoin = "round" ,
                         lineend = "round") +
        scale_color_manual(values = c(  "#3B3029" , #"Persistently low diversity" , 
                                        "#C5C8C7" , #"Stagnant diversity" ,
                                        "#4F9900"   #"Diversified"
        ),
        name = "Diversity shift",
        breaks = c( "Persistently low diversity" , 
                    "Stagnant diversity" ,
                    "Diversified") ) +
        guides(color=guide_legend(ncol=1 ,
                                  override.aes = list(linewidth = 1 )))+
        
        
        geom_point(
                aes(
                        x = PreICB ,
                        y = log2 (OnICB / PreICB) ,
                        
                        fill = `Best response longevity`
                ),
                shape = 21 ,
                size = 4.5 ,
                stroke = 0.1 ) +
        scale_fill_manual(values = RECIST_color_pal$color,
                          name = "Best response",
                          breaks = RECIST_color_pal$response ,
                          guid = "none") +
        
        ggdist::stat_pointinterval(data = diversity_indices  ,
                                   aes(
                                           x = PreICB ,
                                           y = -3
                                   ),
                                   position = position_dodge(width = -0.5, 
                                                             preserve = "single") ,
                                   point_interval = "median_qi" ,
                                   interval_size_range = c(0.1 ,0.75) ,
                                   #.width = c(.5, 0.95)
        ) +
        
        ggrepel::geom_text_repel(
                aes(
                        x = PreICB ,
                        y = log2 (OnICB / PreICB) ,
                        label = paste(COHORT , "; " , Tissue , sep = "")
                ) ,
                direction = "y" ,
                hjust = 0 ,
                #nudge_x = 4,
                nudge_y = 0 ,
                size = 3.5 ,
                min.segment.length = 0 ,
                segment.size = 0.25 ) +
        
        
        ylab (expression(paste("log ( " ,
                               frac("Cycle"[2/3]*" on-ICB tumour Shannon diversity", 
                                    "Baseline tumour Shannon diversity") ,
                               " )" , sep = "")
        ))+
        scale_x_log10(limits = c(5 , 1000))+
        xlab ("Baseline tumour Shannon diversity")+
        theme_minimal()+
        theme (
                panel.grid = element_blank() ,
                
                
                
                axis.line.y = element_line(colour = "#000000" , linewidth = 0.1),
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
                #legend.position='bottom' ,
                aspect.ratio = 1/1.4)

