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
        dplyr::select(Patient_id , Cycle , Order_q , Diversity) %>%
        left_join(clinical_data ,
                  by = "Patient_id") %>%
        mutate(COHORT = fct_relevel(COHORT ,
                                    "SCCHN" , "TNBC" , "HGSOC" , "MM" , "MST"),
               `Best response longevity` = fct_relevel(`Best response longevity` ,
                                             "CR" , "PR" ,
                                             "SD, 6 < n ICB cycles" ,
                                             "SD, n < 6 ICB cycles" ,
                                             "PD"))

### Aesthetics ----------------------------------------------------------------------------------------------------
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

theme <- theme_minimal() +
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


### Figure S4A ---------------------------------------------------------------------------------------------

kruskal.test( Diversity ~ COHORT, 
              data = diversity_indices %>%
                      filter(Order_q == 0) %>%
                      filter(Cycle == "SB"))

# Kruskal-Wallis rank sum test
# Kruskal-Wallis chi-squared = 30.322, df = 4, p-value = 4.208e-06

dunn.test::dunn.test (x = (diversity_indices %>%
                                   filter(Order_q == 0) %>%
                                   filter(Cycle == "SB"))$Diversity,
                      g = (diversity_indices %>%
                                   filter(Order_q == 0) %>%
                                   filter(Cycle == "SB"))$COHORT, 
                      kw = FALSE ,
                      method = "BH")

# (Benjamini-Hochberg)                              
# Col Mean-|
# Row Mean |      HGSOC         MM        MST      SCCHN
# ---------+--------------------------------------------
#       MM |   1.803705
# |            0.0509
# |
#      MST |   3.536907   1.421875
# |            0.0007*    0.0861
# |
#    SCCHN |   5.065558   3.117333   2.189110
# |            0.0000*    0.0023*    0.0286
# |
#     TNBC |   1.597483  -0.286425  -1.852819  -3.592644
# |            0.0688     0.3873     0.0533    0.0008*
#         
# alpha = 0.05
# Reject Ho if p <= alpha/2




ggplot( data = diversity_indices %>%
                filter(Cycle == "SB") %>%
                filter(Order_q == 0),
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
        ylab ("Baseline PBMC Richness")+
        theme




### Figure S4B ---------------------------------------------------------------------------------------------

kruskal.test( Diversity ~ COHORT, 
              data = diversity_indices %>%
                      filter(Order_q == 2) %>%
                      filter(Cycle == "SB"))

# Kruskal-Wallis rank sum test
# Kruskal-Wallis chi-squared = 11.372, df = 4, p-value = 0.02269

dunn.test::dunn.test (x = (diversity_indices %>%
                                   filter(Order_q == 2) %>%
                                   filter(Cycle == "SB"))$Diversity,
                      g = (diversity_indices %>%
                                   filter(Order_q == 2) %>%
                                   filter(Cycle == "SB"))$COHORT, 
                      kw = FALSE ,
                      method = "BH")

# (Benjamini-Hochberg)                              
# Col Mean-|
# Row Mean |      HGSOC         MM        MST      SCCHN
# ---------+--------------------------------------------
#       MM |   1.096903
# |            0.1704
# |
#      MST |   1.200645  -0.085588
# |            0.1916     0.5177
# |
#    SCCHN |   2.795339   1.610548   2.045948
# |            0.0130*    0.1341     0.0679
# |
#     TNBC |   0.021653  -1.124024  -1.250135  -2.918441
# |            0.4914     0.1864     0.2113    0.0176*
#         
# alpha = 0.05
# Reject Ho if p <= alpha/2





ggplot( data = diversity_indices %>%
                filter(Cycle == "SB") %>%
                filter(Order_q == 2),
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
        ylab ("Baseline PBMC Simpson diversity")+
        theme



### Figure S4C ---------------------------------------------------------------------------------------------

kruskal.test( Diversity ~ COHORT, 
              data = diversity_indices %>%
                      filter(Order_q == 0) %>%
                      filter(Cycle == "C3B"))

# Kruskal-Wallis rank sum test
# Kruskal-Wallis chi-squared = 26.578, df = 4, p-value = 2.419e-05

dunn.test::dunn.test (x = (diversity_indices %>%
                                   filter(Order_q == 0) %>%
                                   filter(Cycle == "C3B"))$Diversity,
                      g = (diversity_indices %>%
                                   filter(Order_q == 0) %>%
                                   filter(Cycle == "C3B"))$COHORT, 
                      kw = FALSE ,
                      method = "BH")

# (Benjamini-Hochberg)                              
# 
# Col Mean-|
# Row Mean |      HGSOC         MM        MST      SCCHN
# ---------+--------------------------------------------
#       MM |   1.235967
# |            0.1353
# |
#      MST |   2.796687   1.331450
# |            0.0052*    0.1307
# |
#    SCCHN |   4.864456   3.556430   2.877943
# |            0.0000*    0.0009*    0.0050*
# |
#     TNBC |   1.808733   0.572765  -0.652438  -2.950272
# |            0.0587     0.2834     0.2856    0.0053*
#         
# alpha = 0.05
# Reject Ho if p <= alpha/2



        ggplot( data = diversity_indices %>%
                        filter(Cycle == "C3B") %>%
                        filter(Order_q == 0),
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
                method = "wilcox.test" ,
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
                        "TNBC\nN=11" ,
                        "HGSOC\nN=11" ,
                        "MM\nN=11" ,
                        "MST\nN=26" ),
                position = "bottom" )+
        ylab ("Cycle 3 on-ICB PBMC Richness")+
        theme



### Figure S4D ---------------------------------------------------------------------------------------------

kruskal.test( Diversity ~ COHORT, 
              data = diversity_indices %>%
                      filter(Order_q == 1) %>%
                      filter(Cycle == "C3B"))

# Kruskal-Wallis rank sum test
# Kruskal-Wallis chi-squared = 20.467, df = 4, p-value = 0.0004037

dunn.test::dunn.test (x = (diversity_indices %>%
                                   filter(Order_q == 1) %>%
                                   filter(Cycle == "C3B"))$Diversity,
                      g = (diversity_indices %>%
                                   filter(Order_q == 1) %>%
                                   filter(Cycle == "C3B"))$COHORT, 
                      kw = FALSE ,
                      method = "BH")

# (Benjamini-Hochberg)                              
# Col Mean-|
# Row Mean |      HGSOC         MM        MST      SCCHN
# ---------+--------------------------------------------
#       MM |   0.974706
# |     0.2355
# |
#      MST |   1.940820   0.785308
# |            0.0523     0.2702
# |
#    SCCHN |   4.084350   3.052818   2.858413
# |            0.0002*    0.0038*    0.0053*
# |
#     TNBC |   0.713444  -0.261261  -1.095033  -3.329311
# |            0.2642     0.3969     0.2279    0.0022*
#         
# alpha = 0.05
# Reject Ho if p <= alpha/2




ggplot( data = diversity_indices %>%
                filter(Cycle == "C3B") %>%
                filter(Order_q == 1),
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
                method = "wilcox.test" ,
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
                        "TNBC\nN=11" ,
                        "HGSOC\nN=11" ,
                        "MM\nN=11" ,
                        "MST\nN=26" ),
                position = "bottom" )+
        ylab ("Cycle 3 on-ICB PBMC Shannon diversity")+
        theme



### Figure S4E ---------------------------------------------------------------------------------------------

kruskal.test( Diversity ~ COHORT, 
              data = diversity_indices %>%
                      filter(Order_q == 2) %>%
                      filter(Cycle == "C3B"))

# Kruskal-Wallis rank sum test
# Kruskal-Wallis chi-squared = 11.733, df = 4, p-value = 0.01945

dunn.test::dunn.test (x = (diversity_indices %>%
                                   filter(Order_q == 2) %>%
                                   filter(Cycle == "C3B"))$Diversity,
                      g = (diversity_indices %>%
                                   filter(Order_q == 2) %>%
                                   filter(Cycle == "C3B"))$COHORT, 
                      kw = FALSE ,
                      method = "BH")

# (Benjamini-Hochberg)                              
# Col Mean-|
# Row Mean |      HGSOC         MM        MST      SCCHN
# ---------+--------------------------------------------
#       MM |   0.512474
# |            0.3802
# |
#      MST |   0.869612   0.262074
# |            0.3204     0.4407
# |
#    SCCHN |   2.823420   2.281068   2.488124
# |            0.0119*    0.0282     0.0214*
# |
#     TNBC |  -0.040194  -0.552668  -0.917262  -2.865957
# |           0.4840     0.4146     0.3590     0.0208*
#         
# alpha = 0.05
# Reject Ho if p <= alpha/2



ggplot( data = diversity_indices %>%
                filter(Cycle == "C3B") %>%
                filter(Order_q == 2),
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
                method = "wilcox.test" ,
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
                        "TNBC\nN=11" ,
                        "HGSOC\nN=11" ,
                        "MM\nN=11" ,
                        "MST\nN=26" ),
                position = "bottom" )+
        ylab ("Cycle 3 on-ICB PBMC Simpson diversity")+
        theme
