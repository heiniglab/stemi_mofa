# Visualization Parameters


# General Parameters: Sizes etc.
plot_config = theme_bw() + # background
              theme(
                  panel.grid.major = element_blank(), # background
                  panel.grid.minor = element_blank(), # background
                  panel.background = element_blank(), # background
                  panel.border = element_blank(),
                  axis.line = element_line(colour = "black", size = 1), # background
                
                  
                  # general text
                  text = element_text(family = "sans" ,face = "bold", colour = "black" ,  size = 8) ,# general text size
                  
                  # title text
                  plot.title = element_text(family ="sans" ,face = "bold", colour = "black" ,  size = 8, hjust = 0.5 ),  #hjust: position of title
                  
                  # x axis title
                  axis.title.x = element_text(family = "sans" ,face = "bold", colour = "black" ,  size = 8),
                  
                  # y axis title
                  axis.title.y = element_text(family = "sans" ,face = "bold", colour = "black" ,  size = 8),
                  
                  
                  # axis text x
                  axis.text.x = element_text(color = "black", size = 8, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
                  # axis text y 
                  axis.text.y = element_text(color = "black", size = 8, angle = 0, hjust = 0, vjust = 0, face = "bold"),
                  
                  # legend title
                  legend.title = element_blank(),
                  strip.text.y = element_text(size = 8), strip.text.x = element_text(size = 8)
              
              
              )

# Line Plot Parameters
line_size = 0.5

# Point Plot specific parametres
point_size = 1


# Colors 

## Good/ Bad Category
good_bad_color =   scale_color_manual(values = c("bad" = "#EBAB65",
                                "good"="#428D45"))


## ACS subtypes coloring
acs_subtypes_color =   scale_color_manual(values = c("acs_subacute" = "#94541D",
                                "acs_w_infection"="#FED579", "acs_w_o_infection" = "#E73420"))


## CCS, ACS, no_CCS color
ccs_acs_color =   scale_fill_manual(values = c("CCS" = "#7A4897",
                                "No CCS"="#4359A5", "ACS" = "#E23A3A"))


## ccs, TP colors
ccs_tp_fill =   scale_fill_manual(values = c("CCS" = "#844D9A", "No CCS"="#4359A5", "TP1" = "#E6321E", "TP2" = "#EB5A53", "TP3" = "#F29392", "TP4" = "#F6AFB0"))
ccs_tp_color =   scale_color_manual(values = c("CCS" = "#844D9A","No CCS"="#4359A5", "TP0" = "#844D9A", "TP1" = "#E6321E", "TP2" = "#EB5A53", "TP3" = "#F29392", "TP4" = "#F6AFB0"))

ccs_tp_fill2 =   scale_fill_manual(values = c( "TP0" = "#844D9A", "TP1" = "#E6321E", "TP2" = "#EB5A53", "TP3" = "#F29392", "TP4" = "#F6AFB0"))
ccs_tp_color =   scale_color_manual(values = c("CCS" = "#844D9A","No CCS"="#4359A5", "TP0" = "#844D9A", "TP1" = "#E6321E", "TP2" = "#EB5A53", "TP3" = "#F29392", "TP4" = "#F6AFB0"))

ccs_tp_color2 =   scale_fill_manual(values = c("TP0_ccs" = "#844D9A", "TP1_acs_w_o_infection" = "#E6321E", "TP2_acs_w_o_infection" = "#EB5A53", "TP3_acs_w_o_infection" = "#F29392", "TP4_acs_w_o_infection" = "#F6AFB0", "ccs" = "#844D9A", "no_ccs" = "#4359A5", "TP1_bad" = "#EBAB65", "TP1_good" ="#428D45" ,"TP2_bad" = "#EBAB65", "TP2_good" ="#428D45","TP3_bad" = "#EBAB65", "TP3_good" ="#428D45", "TP4_bad" = "#EBAB65", "TP4_good" ="#428D45"))


## Cluster colors TBD
cell_type_colors =   scale_fill_manual(values = c(
    "CD4+ T cells (Cluster 0)" = "#3c75af", 
    "CD8+ T cells (Cluster 1)" = "#e1833d", 
    "CD4+ T cells (Cluster 2)" = "#4e9b6d", 
    "NK cells (Cluster 3)" = "#c43b36", 
    "CD14(high) Monocytes (Cluster 4)" = "#6d5296", 
    "CD4+ T cells (Cluster 5)" = "#845950",
    "CD14(high) Monocytes (Cluster 6)" = "#c67dac",
    "CD14(high) Monocytes (Cluster 7)" = "#b6bc6d",
    "B cells (Cluster 8)" = "#59b4c6",
    "FCGR3A(high) Monocytes (Cluster 9)" = "#b1c5e4",
    "B cells (Cluster 10)" = "#efbc82",
    "CD4+ T cells (Cluster 11)" = "#a5c88e",
    "FCGR3A(high) Monocytes (Cluster 12)" = "#e69998",
    "Dendritic cells (Cluster 13)" = "#c0b0d2",
    "Clinical" = "#330FED",
    "Cytokine" = "#07580A",
    "Neutrophil" = "#9F8918",
    "Proteomics" = "#128890"
))


# General Parameters: Sizes etc.
plot_config_heatmap = theme_bw() + # background
              theme(
                  panel.grid.major = element_blank(), # background
                  panel.grid.minor = element_blank(), # background
                  #panel.background = element_blank(), # background
                  panel.border = element_blank(),
                  #axis.line = element_line(colour = "black", size = 1), # background
                  axis.ticks = element_blank(),
                
                  
                  # general text
                  text = element_text(family = "sans" ,face = "bold", colour = "black" ,  size = 8) ,# general text size
                  
                  # title text
                  plot.title = element_text(family ="sans" ,face = "bold", colour = "black" ,  size = 8, hjust = 0.5 ),  #hjust: position of title
                  
                  # x axis title
                  axis.title.x = element_text(family = "sans" ,face = "bold", colour = "black" ,  size = 8),
                  
                  # y axis title
                  axis.title.y = element_text(family = "sans" ,face = "bold", colour = "black" ,  size = 8),
                  
                  
                  # axis text x
                  axis.text.x = element_text(color = "black", size = 8, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
                  # axis text y 
                  axis.text.y = element_text(color = "black", size = 8, angle = 0, hjust = 0, vjust = 0, face = "bold"),
                  
                  # legend title
                  legend.title = element_blank(),
                  strip.text.y = element_text(size = 8), strip.text.x = element_text(size =8)
              
              
              )



plot_config_heatmap_supp = theme_bw() + # background
              theme(
                  panel.grid.major = element_blank(), # background
                  panel.grid.minor = element_blank(), # background
                  #panel.background = element_blank(), # background
                  panel.border = element_blank(),
                  #axis.line = element_line(colour = "black", size = 1), # background
                  axis.ticks = element_blank(),
                
                  
                  # general text
                  text = element_text(family = "sans" ,face = "bold", colour = "black" ,  size = 6) ,# general text size
                  
                  # title text
                  plot.title = element_text(family ="sans" ,face = "bold", colour = "black" ,  size = 6, hjust = 0.5 ),  #hjust: position of title
                  
                  # x axis title
                  axis.title.x = element_text(family = "sans" ,face = "bold", colour = "black" ,  size = 6),
                  
                  # y axis title
                  axis.title.y = element_text(family = "sans" ,face = "bold", colour = "black" ,  size = 6),
                  
                  
                  # axis text x
                  axis.text.x = element_text(color = "black", size = 6, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
                  # axis text y 
                  axis.text.y = element_text(color = "black", size = 6, angle = 0, hjust = 0, vjust = 0, face = "bold"),
                  
                  # legend title
                  legend.title = element_blank(),
                  strip.text.y = element_text(size = 6), strip.text.x = element_text(size =6)
              
              
              )
