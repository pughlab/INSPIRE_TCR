library(tidyverse)
library(survminer)
library(survival)
library(forestmodel)
library(patchwork)

### Input data preparation ----------------------------------------------------------------------------------------
data_path <- "https://raw.githubusercontent.com/pughlab/INSPIRE_TCR/main/Data"

### Clinical data:
clinical_data_fname <- "INSPIRE_ClinicalData.csv"

clinical_data <- readr::read_csv(
        file.path (data_path , clinical_data_fname))%>%
        dplyr::select(Patient_id , 
                      COHORT , 
                      `OS (months)` ,
                      `PFS (months)`,
                      `Lost to follow up` ,
                      `PFS Status` , 
                      `OS Status`) %>%
        mutate(CodedCOHORT = ifelse(COHORT == "SCCHN" , COHORT , "Other") )

### Diversity FMM clusters:
FMM_data_fname <- "INSPIRE_PBMC_DiversityFMMClusters.csv"

FMM_data <- as.data.frame (readr::read_csv(
        file.path (data_path , FMM_data_fname))%>%
        dplyr::select(Patient_id ,  `FMM Cluster`) %>%
        left_join(clinical_data,
                  by = "Patient_id") )


FMM_data$`FMM Cluster` <- as.factor(FMM_data$`FMM Cluster`)
FMM_data$`FMM Cluster` <- relevel(FMM_data$`FMM Cluster`, 
                                  ref = "Baseline-High, Diminution")



FMM_data$`CodedCOHORT` <- as.factor(FMM_data$`CodedCOHORT`)
FMM_data$`CodedCOHORT` <- relevel(FMM_data$`CodedCOHORT`, 
                                  ref = "SCCHN")
#-------------------------------------------------------------------------------
#Figure S5A

# Fit the univariate Cox model for FMM_Cluster
cox_univariate_fmm <- coxph(Surv(`PFS (months)`, `PFS Status`) ~ `FMM Cluster`, 
                            data = FMM_data)

summary(cox_univariate_fmm) # Check summary to confirm HRs and p-values


# This plot will show the HR for each level of FMM_Cluster compared to reference cluster
forest_univariate <- forest_model(cox_univariate_fmm ,
                                    format_options = forest_model_format_options(text_size = 3.25 ,
                                                                                 point_size = 2.5) ,

                                    panels =
                                            
                                            list (
                                                    
                                                    list(width = 0.75,
                                                         display = ~variable,
                                                         fontface = "bold",
                                                         heading = "Variable" ,
                                                         hjust = 0),

                                                    list(width = 1.1,
                                                         display = ~level),

                                                    list(width = 0.75,
                                                         display = ~n,
                                                         hjust = 0.5,
                                                         heading = "N" ,
                                                         fontface = "plain"),

                                                    list(width = 0.05,
                                                         item = "vline",
                                                         hjust = 0.5 ),
                                                    #-----------------------------------
                                                    list(
                                                            item = "forest",
                                                            width = 0.4 ,
                                                            hjust = 0.5 ,
                                                            heading = "Hazard ratio" ,
                                                            #The linetype for Hazard ratio vertical line
                                                            linetype = "dashed" ,

                                                            #The position of Hazard ratio vertical line
                                                            line_x = 0 ,
                                                            lwd = 1) ,
                                                    
                                                    list(width = 0.05,
                                                         item = "vline",
                                                         hjust = 0.5),
                                                    #-----------------------------------
                                                    list(width = 0.25,
                                                         display = ~ ifelse(reference,
                                                                            "Reference",
                                                                            sprintf("%0.2f (%0.2f, %0.2f)",
                                                                                    trans(estimate),
                                                                                    trans(conf.low),
                                                                                    trans(conf.high)) ),
                                                         display_na = NA ,
                                                         hjust = 0),

                                                    list(
                                                            width = 0.05,
                                                            display = ~ ifelse(reference, "", 
                                                                               paste0(format.pval(p.value, digits = 1, eps = 0.001), 
                                                                                      case_when (
                                                                                              p.value < 0.001 ~  "***" ,
                                                                                              p.value < 0.01  ~  "**" ,
                                                                                              p.value < 0.05  ~  "*" ,
                                                                                              TRUE ~ ""))),
                                                            display_na = NA, 
                                                            hjust = 0.5,
                                                            heading = "P"
                                                    ) ))
#-------------------------------------------------------------------------------
#Figure S5B

# Fit the multivariate Cox model, adjusting for COHORT
cox_multivariate_fmm_cohort <- coxph(Surv(`PFS (months)`, `PFS Status`) ~ `FMM Cluster` + CodedCOHORT, 
                                     data = FMM_data)
summary(cox_multivariate_fmm_cohort) # Check summary

# Generate and print the forest plot
ggforest_multivariate <- forest_model(cox_multivariate_fmm_cohort ,
                                      format_options = forest_model_format_options(text_size = 3.25 ,
                                                                                   point_size = 2.5) ,
                                      
                                      panels =
                                              
                                              list (
                                                      
                                                      list(width = 0.75,
                                                           display = ~variable,
                                                           fontface = "bold",
                                                           heading = "Variable" ,
                                                           hjust = 0),
                                                      
                                                      list(width = 1.1,
                                                           display = ~level),
                                                      
                                                      list(width = 0.75,
                                                           display = ~n,
                                                           hjust = 0.5,
                                                           heading = "N" ,
                                                           fontface = "plain"),
                                                      
                                                      list(width = 0.05,
                                                           item = "vline",
                                                           hjust = 0.5 ),
                                                      #-----------------------------------
                                                      list(
                                                              item = "forest",
                                                              width = 0.4 ,
                                                              hjust = 0.5 ,
                                                              heading = "Hazard ratio" ,
                                                              #The linetype for Hazard ratio vertical line
                                                              linetype = "dashed" ,
                                                              
                                                              #The position of Hazard ratio vertical line
                                                              line_x = 0 ,
                                                              lwd = 1) ,
                                                      
                                                      list(width = 0.05,
                                                           item = "vline",
                                                           hjust = 0.5),
                                                      #-----------------------------------
                                                      list(width = 0.25,
                                                           display = ~ ifelse(reference,
                                                                              "Reference",
                                                                              sprintf("%0.2f (%0.2f, %0.2f)",
                                                                                      trans(estimate),
                                                                                      trans(conf.low),
                                                                                      trans(conf.high)) ),
                                                           display_na = NA ,
                                                           hjust = 0),
                                                      
                                                      list(
                                                              width = 0.05,
                                                              display = ~ ifelse(reference, "", 
                                                                                 paste0(format.pval(p.value, digits = 1, eps = 0.001), 
                                                                                        case_when (
                                                                                                p.value < 0.001 ~  "***" ,
                                                                                                p.value < 0.01  ~  "**" ,
                                                                                                p.value < 0.05  ~  "*" ,
                                                                                                TRUE ~ ""))),
                                                              display_na = NA, 
                                                              hjust = 0.5,
                                                              heading = "P"
                                                      ) ))
#-------------------------------------------------------------------------------
design <- "A
B"

ggsave (filename = "Uni-and-MultiVariate_Survival.svg" ,
        device = "svg" ,
        dpi = 1000 ,
        units = "px" ,
        width = 7500 ,
        height = 5000)

print (ggforest_univariate + 
               ggforest_multivariate + 
               patchwork::plot_layout(design = design ,
                                      heights = c(0.9 , 1.1)))
dev.off()
