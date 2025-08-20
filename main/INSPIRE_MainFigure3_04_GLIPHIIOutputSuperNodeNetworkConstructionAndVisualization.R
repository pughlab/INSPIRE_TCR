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


### ---------------------------------------------------------------------------------------------------------
### Improved Layout Design (Layout Version Designed For Thesis Defense): ------------------------------------

# This layout segregates components by community composition. Components with a single community are positioned
# in five peripheral circles, while multi-community components are centralized.
# Developed with assistance from Google's Gemini.


# ==========================================================
# ============== NEW: Parameters to Tune ===================
# ==========================================================
# 1. To expand the central components
core_scale_factor <- 1.6 # (60% larger)

# 2. To compress the peripheral circles
radius_increment_factor <- 0.09 # (Smaller number = closer circles)

# 3. To create partial (e.g., 4/5) circles, leaving a gap at the bottom
# The arc will go from start_angle to end_angle (in radians)
# A full circle is 2*pi (~6.28). These values create a 4/5 arc with a gap at the bottom.
arc_proportion <- 4/5
gap_proportion <- 1 - arc_proportion
start_angle <- (gap_proportion / 2) * 2 * pi 
end_angle <- (1 - gap_proportion / 2) * 2 * pi

# Number of circles for the periphery
num_circles <- 5
# ==========================================================


# Step 1: Calculate the Initial Backbone Layout
# =============================================
Backbone3 <- graphlayouts::layout_as_backbone(Abstract_Network, keep = 0.75)


# Step 2: Identify and Separate Components
# ========================================
comps <- components(Abstract_Network)
membership <- comps$membership
csize <- comps$csize

isolate_ids <- which(csize[membership] == 1)
main_component_ids <- which(csize[membership] > 1)
coords_main <- Backbone3$xy[main_component_ids, ]


# Step 3: Arrange Components with New Scaling and Arcs
# ====================================================
# --- Split the isolates into groups for each circle ---
isolate_groups <- vector("list", num_circles)
for (i in 1:num_circles) {
        isolate_groups[[i]] <- isolate_ids[seq(i, length(isolate_ids), by = num_circles)]
}

# --- Recenter and SCALE the main components ---
center_main <- apply(coords_main, 2, mean)
coords_main_centered <- sweep(coords_main, 2, center_main, "-")
# --- Scale the core components to make them larger
coords_main_scaled <- coords_main_centered * core_scale_factor

# --- Create an empty list for isolate coordinates ---
isolate_coords_list <- vector("list", num_circles)

# --- Determine the radii for the compressed circles ---
# NOTE: Use the newly scaled coordinates to calculate the starting radius
max_dist <- max(abs(coords_main_scaled))
base_radius <- max_dist * 1.20 # Start first circle a bit closer
radius_increment <- base_radius * radius_increment_factor # Use the new smaller factor

# --- Loop to calculate coordinates for each ARC of isolates ---
for (i in 1:num_circles) {
        current_isolate_ids <- isolate_groups[[i]]
        num_current_isolates <- length(current_isolate_ids)
        
        if (num_current_isolates > 0) {
                current_radius <- base_radius + (i - 1) * radius_increment
                
                # Use the start_angle and end_angle to create an arc
                angles <- seq(start_angle, end_angle, length.out = num_current_isolates)
                
                # Add a small offset to stagger points (looks better on arcs too)
                if (num_current_isolates > 1) {
                        angle_step <- (angles[2]-angles[1])
                        angle_offset <- (angle_step / num_circles) * (i-1) - (angle_step/2)
                        angles <- angles + angle_offset
                }
                
                new_coords <- matrix(NA, nrow = num_current_isolates, ncol = 2)
                new_coords[, 1] <- current_radius * cos(angles)
                new_coords[, 2] <- current_radius * sin(angles)
                
                isolate_coords_list[[i]] <- new_coords
        }
}


# Step 4: Combine into a Final Layout Matrix
# ===========================================
final_layout <- matrix(NA, nrow = vcount(Abstract_Network), ncol = 2)
# Use the SCALED coordinates for the main components
final_layout[main_component_ids, ] <- coords_main_scaled
# Populate the layout for the isolates
for (i in 1:num_circles) {
        final_layout[isolate_groups[[i]], ] <- isolate_coords_list[[i]]
}

### ---------------------------------------------------------------------------------------------------------


tkid <- tkplot(
        Abstract_Network,
        # vertex.size = ifelse(V (Abstract_Network)$Source == "HEALTHYcfDNA,INSPIRETumour" ,
        #                      V (Abstract_Network)$Community_Size**0.85 ,
        #                      0) ,             # Adjust the size of the super-nodes
        vertex.size =  V (Abstract_Network)$Community_Size**0.51 ,
        vertex.color = V (Abstract_Network)$Color ,
        # vertex.label = ifelse(V(Abstract_Network)$Detailed_Annotation %in% c("IOKIN_intrinsic" , "HNSCC") ,
        #                       NA ,
        #                       V(Abstract_Network)$Detailed_Annotation ),
        vertex.label = NA ,
        vertex.label.color = "#000000" ,
        # vertex.label.dist = ifelse(V(Abstract_NetworkPrime)$EdgeBetwCommunity_Size > 25  ,
        #                            1.7 , 1
        # ),
        # vertex.frame.color = ifelse(V (Abstract_Network)$ResistanceStatus_Range == "Ar,Pr" ,
        #                             "red" ,
        #                             ifelse(V (Abstract_Network)$ResistanceStatus_Range == "Ar" ,
        #                                    "blue" ,
        #                                    "#000000")) ,
        vertex.label.cex = 0.5 ,
        vertex.label.degree = 75 ,
        vertex.label.family = "Helvetica" ,
        vertex.shape = "circle" ,
        edge.width = 0.1 ,
        edge.color = "#000000" ,
        rescale = FALSE ,
        ###To get all the available layouts: grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
        layout = norm_coords(#layout [, c(1,2)] ,
                final_layout ,
                ymin=-0.9,
                ymax=0.9,
                xmin=-0.9,
                xmax=0.9))



ManualLayout <- tk_coords(tkid)


Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "LayoutX",
                                    value = CrappyManualLayout [,1])

Abstract_Network <- set_vertex_attr(Abstract_Network,
                                    name = "LayoutY",
                                    value = CrappyManualLayout [,2])


### Saving the final super-node network:
saveRDS(Abstract_Network, "INSPIRE_Tumour_GLIPHIISuperNodeNetwork.rds")


### ---------------------------------------------------------------------------------------------------------
### Visualization with no annotations : ---------------------------------------------------------------------

svg(filename = "INSPIRE_Tumour_SuperNodeNetwork_unfiltered_UnAnnotated.svg"),
    family = "Helvetica" ,
    symbolfamily = "Helvetica" ,
    width = 6.0, height = 6.0,
    onefile = TRUE ,
    bg = "transparent")


plot(
        simplify (Abstract_Network),
        vertex.size = V (Abstract_Network)$Community_Size**0.51 ,             
        vertex.color = V (Abstract_Network)$Color,  
        vertex.label = NA ,
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
        layout = norm_coords(ManualLayout    ,
                             ymin=-0.98, 
                             ymax=0.98, 
                             xmin=-0.98, 
                             xmax=0.98))


dev.off()

### ---------------------------------------------------------------------------------------------------------
### Visualization of TRBV convergent communities with no annotations : --------------------------------------

svg(filename = "INSPIRE_Tumour_SuperNodeNetwork_TRBVConvergent_UnAnnotated.svg"),
    family = "Helvetica" ,
    symbolfamily = "Helvetica" ,
    width = 6.0, height = 6.0,
    onefile = TRUE ,
    bg = "transparent")


plot(
        simplify (Abstract_Network),
        vertex.size = V (Abstract_Network)$Community_Size**0.51 ,             
        vertex.color = ifelse(V (Abstract_Network)$Community_id %in% (GLIPHII_Community_stats %>% filter(TRBVConvergence != "Divergent"))$Community_id ,
                              V (Abstract_Network)$Color,  
                              "transparent" ) ,
        
        vertex.frame.color = ifelse(V (Abstract_Network)$Community_id %in% (GLIPHII_Community_stats %>% filter(TRBVConvergence != "Divergent"))$Community_id ,
                                    "#000000",  
                                    "transparent") ,
        
        vertex.label = NA ,
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
        layout = norm_coords(as.matrix(data.frame(V(Abstract_Network)$LayoutX , V(Abstract_Network)$LayoutY))    ,
                             ymin=-0.98, 
                             ymax=0.98, 
                             xmin=-0.98, 
                             xmax=0.98))


dev.off()
