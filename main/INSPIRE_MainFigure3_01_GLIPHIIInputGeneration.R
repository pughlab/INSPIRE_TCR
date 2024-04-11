### This script prepares the input required for the GLIPHII algorithm run.
### The script takes patient-derived TRB CDR3 sequences, combines it with a set of external TCR sequences 
### and formats the combined TCR sequences into the standard input format required by the GLIPHII algorithm.
### ----------------------------------------------------------------------------------------

library(tidyverse)

data_path <- "https://raw.githubusercontent.com/pughlab/INSPIRE_TCR/main/Data"

### Clinical data:
clinical_data_fname <- "INSPIRE_ClinicalData.csv"
clinical_data <- readr::read_csv(
        file.path (data_path , clinical_data_fname))%>%
        dplyr::select(Patient_id , 
                      COHORT , 
                      `Best response longevity` ) %>%
        mutate(`Best response longevity` = case_when(
                `Best response longevity`   == "CR"  ~ "CR",
                `Best response longevity`   == "PR"  ~ "PR",
                `Best response longevity`   == "SD, 6 < n ICB cycles"  ~ "LongTremSD",
                `Best response longevity`   == "SD, n < 6 ICB cycles"  ~ "ShortTermSD",
                `Best response longevity`   == "PD"  ~ "PD",
                `Best response longevity`   == "NE"  ~ "NE" ))


### Tumour clonotypes:
TumourClonotypes_fname <- "INSPIRE_Tumour_TRBClonotypes.csv"
TumourClonotypes <- readr::read_csv(file.path( data_path , TumourClonotypes_fname )) %>%
        filter(Patient_id != "INS-B-025") %>%
        mutate(Patient_id = if_else(Patient_id == "INS-E-028_L" | Patient_id == "INS-E-028_R" , 
                                    "INS-E-028" , 
                                    Patient_id)) %>%
        mutate(TRBV = unlist(
                sapply(
                        strsplit(allVHitsWithScore, split = "*", fixed = TRUE), 
                        function(x) x[[1]][1], 
                        simplify=FALSE)) ,
               
               TRBJ = unlist(
                       sapply(
                               strsplit(allJHitsWithScore, split = "*", fixed = TRUE), 
                               function(x) x[[1]][1], 
                               simplify=FALSE)) ) %>%
        
        dplyr::select(Patient_id , Cycle , 
                      aaSeqCDR3 , TRBV , TRBJ ) %>%
        
        unique() %>%
        
        arrange(Patient_id , match(Cycle, 
                                   c("ST" , "C2T" , "C3T" , "EOTT"))) %>%
        
        group_by(Patient_id , aaSeqCDR3 ) %>%
        
        mutate(TRBV = str_replace_all( paste(sort (unique (TRBV) ) , collapse = "-or"), "(?!^)TRBV", ""),
               TRBJ = str_replace_all( paste(sort (unique (TRBJ) ) , collapse = "-or"), "(?!^)TRBJ", "")) %>%
        
        ungroup() %>%
        
        unique() %>%
        
        mutate(cloneCount = 1) %>%
        group_by(Patient_id , aaSeqCDR3 , TRBV , TRBJ) %>%
        pivot_wider(values_from = cloneCount ,
                    names_from = Cycle) %>%
        ungroup() %>%
        unique() %>%
        rowwise() %>%
        mutate(Cycle = paste(
                c("ST" , "C2T" , "C3T" , "EOTT") [
                        which (
                                is.na (c_across(c("ST" , "C2T" , "C3T" , "EOTT"))) == 0)] , 
                collapse = "-")) %>%
        left_join(
                clinical_data ,
                by = "Patient_id") %>%
        mutate(
                `subject:condition` = paste( "INSPIRE" , (paste(Patient_id ,
                                                                Cycle ,
                                                                COHORT ,
                                                                `Best response longevity` , sep = "_")) , 
                                             sep = ":") ,
                CDR3a = NA ,
                
                count = 1 ) %>%
        dplyr::select(aaSeqCDR3 , TRBV , TRBJ , CDR3a ,
                      `subject:condition` , count) %>%
        rename(CDR3b = aaSeqCDR3)


### NSCLC dataset from  PMID 33691136 :
mdavis_fname <- "MarkDavis_TumorEnrichedTRB-CDR3.tsv"
mdavis_spcf <- readr::read_tsv(file.path( data_path , mdavis_fname ))


### deorphanized set of TCRs :
DeOrphanizedTCRs_fname <- "TCRDB_PublicRelease.tsv"
DeOrphanizedTCRs <- readr::read_tsv( file.path( data_path , DeOrphanizedTCRs_fname )) %>%
        filter(!(is.na(TRBV))) %>%
        mutate(TRBV = unlist(sapply(strsplit(TRBV, split = "*", fixed = TRUE), function(x) x[[1]][1], simplify=FALSE)) ,
               TRBJ = unlist(sapply(strsplit(TRBJ, split = "*", fixed = TRUE), function(x) x[[1]][1], simplify=FALSE)) ,
               CDR3a = NA ,
               count = 1)%>% 
        unite("pMHC",
              c(`MHC` , Epitope , `Epitope gene`),
              sep = "__")%>%
        unite("subject:condition",
              c(`Epitope species` , pMHC),
              sep = ":")%>% 
        select(CDR3b ,
               TRBV ,
               TRBJ ,
               CDR3a ,
               `subject:condition` ,
               count)

### Pooling Patient-derived and external TCRs and generating the final table:
write.table(
        bind_rows(TumourClonotypes , 
                  mdavis_spcf , 
                  DeOrphanizedTCRs), 
            file = "INSPIRE_CumulativeTumour.txt" ,
            col.names = FALSE , row.names = FALSE , 
            quote = FALSE ,
            sep = "\t" )

