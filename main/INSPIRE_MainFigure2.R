library(tidyverse)
library (ggbeeswarm)

### Input data preparation ----------------------------------------------------------------------------------------

### Clinical data:
clinical_data_path <- "/Users/shirin/Desktop/Immunarch/downsamplign/Projetcs/INSPIRE/clinical_data/coding_ref"
clinical_data_fname <- "INSPIRE_ClinicalData_CodingRef.xlsx"

clinical_data <- readxl::read_xlsx(
        file.path (clinical_data_path , clinical_data_fname))%>%
        dplyr::select(Patient_id , 
                      COHORT , 
                      `Best response longevity` )

### Metastasis site data:
MetSite_data_path <- "/Users/shirin/Desktop/Immunarch/downsamplign/Projetcs/INSPIRE/clinical_data/coding_ref"
MetSite_data_fname <- "Metastasis_site.xlsx"
MetSite_data <- readxl::read_xlsx(file.path (MetSite_data_path , MetSite_data_fname))

### Diversity indices:
input_path  <- "~/Desktop/Immunarch/downsamplign/Projetcs/INSPIRE/dna/tumor/diversity_profiles"
input_fname <- "inspire_tumor_diversity_inices.csv"

diversity_indices <- readr::read_csv(
        file.path (input_path , input_fname)) %>% 
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
                              
