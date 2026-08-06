theme_pub <- theme_classic(base_size = 9) +
         theme(
                  panel.border = element_blank(),
                  axis.line = element_line(colour = "black", linewidth = 0.5),
                  axis.text = element_text(colour = "black"),
                  legend.title = element_text(face = "bold"),
                  legend.position = "right",
                  plot.title = element_text(face = "bold", hjust = 0.5),
                  strip.background = element_blank(),
                  strip.text = element_text(face = "bold")
         )
cluster_cols <- c(
         "#0072B2","#E69F00","#009E73","#CC79A7","#56B4E9",
         "#D55E00","#F0E442","#000000","#882255","#44AA99",
         "#117733","#999933","#88CCEE","#AA4499","#DDCC77",
         "#332288","#661100","#6699CC","#AA4466","#228833"
)

circlize::colorRamp2(
         c(-2, 0, 2),
         c("#2C7BB6",
           "white",
           "#D7191C")
)
