library(tidyverse)
library(igraph)
library(graphlayouts)

data_path <- "https://raw.githubusercontent.com/pughlab/INSPIRE_TCR/main/Data"

### ---------------------------------------------------------------------------------------------------------
### Reading in the network and the community features dataframe: --------------------------------------------

network_fname <- "INSPIRE_Tumour_GLIPHIINetwork.rds"
network <- readr::read_rds(file.path( data_path , network_fname ))


communityStats_fname <- "INSPIRE_Tumour_GLIPHIICommunitiesDetails.csv"
GLIPHII_Community_stats <- read_csv(file.path( data_path , communityStats_fname ))

### ---------------------------------------------------------------------------------------------------------
### Condensing the nodes participating in each community into supoer-nodes: ---------------------------------


GLIPHII_Community_stats$Mock_Community_id <- 1:nrow(GLIPHII_Community_stats)

network <- set_vertex_attr(network,
                           name = "Mock_Community_id",
                           value = GLIPHII_Community_stats$Mock_Community_id [match(V(network)$LeidenCommunity , 
                                                                                    GLIPHII_Community_stats$Community_id) ])


Abstract_Network <- igraph::contract.vertices(
        graph = network ,
        mapping = V(network)$Mock_Community_id)
        ###The mapping argument passed to contract function is a numeric vector that specifies the mapping. 
        ### We pass the community ids as the mapping argument. 
        ### However, because we've trimmed some communities, some community ids are missing.
        ### For those withdrawn communities, empty super nodes are generated, which are not favorable. 
        ### Therefore, we need to provide mock community ids to prevent empty supernode generation. 
        ### We add a new column to GLIPHII_Community_stats as Mock_Community_id to do this.

### Adding attributes to each super-node from GLIPHII_Community_stats dataframe:

Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "Community_Size", 
                                    value = GLIPHII_Community_stats$Community_size)

Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "Component_id", 
                                    value = GLIPHII_Community_stats$Component_id)

Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "Community_id", 
                                    value = GLIPHII_Community_stats$Community_id)


Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "Community_specificity", 
                                    value = GLIPHII_Community_stats$Community_specificity)

Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "Detailed_Annotation", 
                                    value = GLIPHII_Community_stats$Detailed_Annotation)

Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "Abstract_Annotation", 
                                    value = GLIPHII_Community_stats$Abstract_Annotation)

Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "Color", 
                                    value = GLIPHII_Community_stats$Color)

Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "Number_of_Patients", 
                                    value = GLIPHII_Community_stats$Number_of_Patients)

Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "NumberOf_PatientDerived_Nodes", 
                                    value = GLIPHII_Community_stats$NumberOf_PatientDerived_Nodes)

Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "RECIST_Range", 
                                    value = GLIPHII_Community_stats$RECIST_Range)

Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "Cohort_Range", 
                                    value = GLIPHII_Community_stats$Cohort_Range)

Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "Number_of_Cohorts", 
                                    value = GLIPHII_Community_stats$Number_of_Cohorts)

Abstract_Network <- set_vertex_attr(Abstract_Network, 
                                    name = "Cycle", 
                                    value = GLIPHII_Community_stats$Cycle)

### Removing multi-edges:
Abstract_Network <- simplify(Abstract_Network)

### ---------------------------------------------------------------------------------------------------------
### Super-node network visualization: -----------------------------------------------------------------------

### Defining the coordination of each super-node using layouts readily available in both igraph and 
### graphlayouts packages. This step is up to the user and based on user's visualization preferences,
### Any other packages or algorithms can be utilized:


fr_layout <- layout_with_fr (Abstract_Network)

metro_layout <- graphlayouts::layout_as_metromap (Abstract_Network , 
                                            xy = fr_layout , 
                                            gr = 0.01)


graphopt_layout <- layout_with_graphopt(
        Abstract_Network ,
        start = metro_layout ,
        niter = 10000 ,
        charge = 0.0025 )

### For extra, minimal manual tweaking, we can use tkid:

tkid <- tkplot(
        Abstract_Network ,
        
        ## Adjust the size of the super-nodes relevantto their community size:
        vertex.size = V (Abstract_Network)$Community_Size**0.5 ,   
        
        # Color the super-nodes based on community GLIPHII-defined specificity:
        vertex.color = V (Abstract_Network)$Color,  
        
 
        vertex.label = NA ,
        vertex.shape = "circle" ,
        edge.width = 0.5 ,
        edge.color = "#000000" ,
        rescale = FALSE ,
        ###To get all the available layouts: grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
        layout = norm_coords(graphopt_layout ,
                             ymin = -1, 
                             ymax =  1, 
                             xmin = -1, 
                             xmax =  1))

### Get the modified layout:
layout <- tkplot.getcoords(tkid)

### We can save the layout for later use or add the X and Y coordinations 
### as node attributes tothe super-node network:


# write_csv(as.data.frame(layout2) ,
#           file.path(
#                   community_path ,
#                   "INSPIRE_CumulativeTumour_GLIPHII_SpecificityLandscape_SuperNodes_VisualizationCoordinations.csv"))


# Abstract_Network <- set_vertex_attr(Abstract_Network,
#                                     name = "LayoutX",
#                                     value = layout [,1])
# 
# Abstract_Network <- set_vertex_attr(Abstract_Network,
#                                     name = "LayoutY",
#                                     value = layout [,2])


### Saving the final super-node network:
saveRDS(Abstract_Network,
        file.path(community_path ,
                  "INSPIRE_Tumour_GLIPHIISuperNodeNetwork.rds") )

### ---------------------------------------------------------------------------------------------------------
### Visualization with annotation of the Community ID for zoom-ins within each circle: ----------------------

plot(
        simplify (Abstract_Network),
        vertex.size = V (Abstract_Network)$Community_Size**0.51 ,             
        vertex.color = V (Abstract_Network)$Color,  
        vertex.label = V(Abstract_Network)$Community_id ,
        #vertex.label = NA ,
        vertex.label.color = "#000000" ,
        vertex.label.cex = 0.4 , 
        vertex.label.degree = 75 , 
        vertex.label.family = "Helvetica" , 
        vertex.frame.width = 0.1 ,
        vertex.shape = "circle" ,
        edge.width = 0.5 ,
        edge.color = "#BDBDBB" ,
        rescale = FALSE ,
        layout = norm_coords(graphlayouts::layout_rotate(graphlayouts::layout_mirror(layout) , 
                                                         angle = -50)    ,
                             ymin=-0.98, 
                             ymax=0.98, 
                             xmin=-0.98, 
                             xmax=0.98))
