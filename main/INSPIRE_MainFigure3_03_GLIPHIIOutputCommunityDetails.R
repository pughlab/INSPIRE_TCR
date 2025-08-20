library(tidyverse)
library(igraph)
library(graphlayouts)

data_path <- "https://raw.githubusercontent.com/pughlab/INSPIRE_TCR/main/Data"

### ---------------------------------------------------------------------------------------------------------
### Reading in the network: ---------------------------------------------------------------------------------

network_fname <- "INSPIRE_Tumour_GLIPHIINetwork.rds"
network <- readr::read_rds(file.path( data_path , network_fname ))

### ---------------------------------------------------------------------------------------------------------
### Storing community specificities and their details in a dataframe for downstream analysis: ---------------


### Defining the variable we'd like to extract:
GLIPHII_Community_stats <- as.data.frame(table(V(network)$LeidenCommunity)) %>%
        set_names(c("Community_id", "Community_size")) %>%
        mutate(Component_id = NA ,
               Component_size = NA ,
               PatientDerived_Nodes = NA ,
               NumberOf_PatientDerived_Nodes = NA ,
               Number_of_Patients = NA ,
               Cycle = NA ,
               RECIST_Range = NA ,
               Cohort_Range = NA ,
               Number_of_Cohorts = NA ,
               Community_specificity = NA ,
               KnownExternalTCRs = NA)


### Looping through the communities to extract the features:
for (n in c(1:nrow(GLIPHII_Community_stats))) {
        
        community_id <- GLIPHII_Community_stats$Community_id [n]
        
        Component_id <- unique(components(network)$membership [V(network)$name [V(network)$LeidenCommunity == community_id]])
        
        GLIPHII_Community_stats$Component_id [n] <- Component_id
        GLIPHII_Community_stats$Component_size [n] <- components(network)$csize [Component_id]
        
        GLIPHII_Community_stats$Community_specificity [n] <- ifelse(sum (!grepl("INSPIRE" , 
                                                                                unique(V(network)$Source [V(network)$LeidenCommunity == community_id]))) > 0 ,
                                                                    paste(
                                                                            sort (
                                                                                    grep("INSPIRE" , 
                                                                                         unique(V(network)$Source [V(network)$LeidenCommunity == community_id]) , 
                                                                                         value = TRUE , 
                                                                                         invert = TRUE)) ,
                                                                            
                                                                            collapse = ",") ,
                                                                    "Patient_intrinsic") 
        
        
        GLIPHII_Community_stats$PatientDerived_Nodes [n] <- paste(
                sort (
                        unique(V(network)$name [V(network)$LeidenCommunity == community_id & 
                                                        V(network)$Source == "INSPIRE"])) ,
                collapse = ",")
        
        
        
        GLIPHII_Community_stats$NumberOf_PatientDerived_Nodes [n] <- length(unique(V(network)$name [V(network)$LeidenCommunity == community_id & 
                                                                                                            V(network)$Source == "INSPIRE"]))
        
        
        
        GLIPHII_Community_stats$Number_of_Patients [n] <- length(unique(V(network)$Patient_id [V(network)$LeidenCommunity == community_id & 
                                                                                                       V(network)$Source == "INSPIRE"]))
        
        
        GLIPHII_Community_stats$KnownExternalTCRs [n] <- paste(
                sort (
                        unique(V(network)$name [V(network)$LeidenCommunity == community_id & 
                                                        V(network)$Source != "INSPIRE" & V(network)$Source != "MDavis"])) ,
                collapse = ",")
        
        
        GLIPHII_Community_stats$RECIST_Range [n] <- unique(
                V(network)$RECIST [V(network)$LeidenCommunity == community_id & 
                                           V(network)$Source == "INSPIRE"]) %>% 
                {.[order(match(.,  c("CR" , "PR" , "LongTremSD" , "ShortTermSD" , "PD" , "NE")))]} %>%
                paste(collapse = ",")
        
        
        
        
        
        GLIPHII_Community_stats$Cohort_Range [n] <- unique(
                V(network)$COHORT [V(network)$LeidenCommunity == community_id & 
                                           V(network)$Source == "INSPIRE"]) %>% 
                {.[order(match(.,  c("SCCHN" , "TNBC" , "HGSOC" , "MM" , "MST" )))]} %>%
                paste(collapse = ",")
        
        
        GLIPHII_Community_stats$Number_of_Cohorts [n] <- length(unique(V(network)$COHORT [V(network)$LeidenCommunity == community_id & 
                                                                                                  V(network)$Source == "INSPIRE"]))
        
        GLIPHII_Community_stats$Cycle [n] <- unique(V(network)$Cycle [V(network)$LeidenCommunity == community_id & 
                                                                              V(network)$Source == "INSPIRE"]) %>% 
                strsplit(split = "-", fixed = TRUE) %>%
                unlist () %>% 
                str_replace_all("C2T|C3T|EOTT", "OnICB") %>%
                unique()%>%
                {.[order(match(.,  c("ST" , "OnICB")))]} %>%
                paste(collapse = ",")
        
        rm (community_id , Component_id)
        
}


### Defining the functions for specificity community de-orphanization:
Annotation_Ref <- tibble(
        ExternalSpecificity = c("Patient_intrinsic" ,
                                "MDavis" ,
                                "HomoSapiens" ,
                                
                                "CEF" ,  "CMV" , 
                                "EBV" ,  "HCV" , 
                                "HPV" ,  "MCPyV" ,
                                "Influenza" , 
                                
                                "S-pneumoniae" , "M.tuberculosis" , 
                                "DENV" , "HTLV-1" , "YFV") ,
        
        Abstract_Annotation = c("Patient_intrinsic" ,
                                "MDavis" ,
                                "Human tumour and auto-antigens" ,
                                
                                rep("Viral", 7), 
                                rep("Other pathogens", 5)) )



# Function to extract HomoSapiens genes
getGenes <- function(string) {
        fullStrings <- grep("HomoSapiens" , 
                            unlist (strsplit(string , ",")),
                            value = TRUE) 
        genes <- paste (unique(sapply(
                strsplit (fullStrings , "__"), 
                function(x) x[[3]])),
                collapse = ",")
        
        return(genes)
}




# Function for Abstract Annotation:
AbstractAnnotate <- function(string) {
        if (length(unlist (str_split(string , ","))) == 1) {
                annotation <- Annotation_Ref$Abstract_Annotation [Annotation_Ref$ExternalSpecificity == string]
        }
        
        else if (length(unlist (str_split(string , ","))) == 2 &
                 grepl("MDavis" , string) ) {
                
                annotation <- Annotation_Ref$Abstract_Annotation [Annotation_Ref$ExternalSpecificity == grep ("MDavis" , 
                                                                                                              unlist(str_split(string , ",")) ,
                                                                                                              value = TRUE ,
                                                                                                              invert = TRUE)]
        }
        
        else {
                annotation <- "Cross-species"
        }
        
        return(annotation)
}



# Function for Detailed Annotation:
DetailAnnotate <- function(string) {
        if (length(unlist (str_split(string , ","))) == 1) {
                annotation <- string
        }
        else {
                annotation <- paste( grep( "MDavis",
                                           unlist (str_split(string , ",")) ,
                                           value = TRUE ,
                                           invert = TRUE) , collapse = ",")
                
        }
        
        return(annotation)
}





Color_Ref <- tibble(
        Abstract_Annotation = c("Patient_intrinsic" ,
                                "MDavis" ,
                                "Human tumour and auto-antigens" ,
                                "Viral" ,
                                "Other pathogens" ,
                                "Cross-species") ,
        Color = c(          "#777777" ,
                            "#CDEDFA" ,
                            "#FF9F01" ,
                            "#AC261B" ,
                            "#BAC70D" ,
                            "#FFDF01"))


### Assigning the specificities:
GLIPHII_Community_stats <- GLIPHII_Community_stats %>%
        rowwise() %>%
        mutate(Abstract_Annotation = AbstractAnnotate (Community_specificity) ,
               Detailed_Annotation = DetailAnnotate (Community_specificity)) %>%
        mutate(Detailed_Annotation = str_replace(Detailed_Annotation , 
                                                 "HomoSapiens" , 
                                                 getGenes (KnownExternalTCRs))) %>%
        left_join(Color_Ref ,
                  by = "Abstract_Annotation")


### Screening the TRBVs for each community:

# communities_of_interest <- (GLIPHII_Community_stats %>% filter(Abstract_Annotation == "MDavis"))$Community_id 
communities_of_interest <- c()

for (community_id in communities_of_interest) {
        
        # Print a header for the community to keep the output organized
        cat(paste("\n--- Vertices in LeidenCommunity:", community_id, "---\n"))
        
        vertex_data <- data.frame(
                Patient_id = V(network)$Patient_id[V(network)$LeidenCommunity == community_id],
                Source = V(network)$Source[V(network)$LeidenCommunity == community_id],
                TcRb   = V(network)$TcRb[V(network)$LeidenCommunity == community_id],
                V      = V(network)$V[V(network)$LeidenCommunity == community_id]) %>%
                arrange(V)

        print(vertex_data)
                
}
### Saving the final dataframe:
write_csv(GLIPHII_Community_stats ,
          "INSPIRE_Tumour_GLIPHIICommunitiesDetails.csv" ))
