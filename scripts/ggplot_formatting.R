##
library(here)
here::i_am("scripts/ggplot_formatting.R")

# color palette for species
species.palette <- c("alone" = "#FFFFFF", "JB7" = "#85e6fb", "JB5" = "#1451a2", "BC10" = "#00b153",
                     "BC9" = "#2ee863", "JB370" = "#a15100", "JBC" = "#a2a2a2",
                     "135E" = "#d7cb6c", "community" = "#000000")

# associate genus names with strain designation
names <- c(`135E`="Diutina", 
           `BC10`="S. xylosus", 
           `BC9`="S. equorum", 
           `JB5`="Brevibacterium", 
           `JB7`="Brachybacterium", 
           `JBC`="Penicillium", 
           `JB370`="Scopulariopsis", 
           `community`="community" , 
           `alone`="monoculture")

names.italic <- c(`135E`=expression(italic("Diutina")), 
                  `BC10`=expression(italic("S. xylosus")), 
                  `BC9`=expression(italic("S. equorum")), 
                  `JB5`=expression(italic("Brevibacterium")), 
                  `JB7`=expression(italic("Brachybacterium")), 
                  `JBC`=expression(italic("Penicillium")), 
                  `JB370`=expression(italic("Scopulariopsis")), 
                  `community`="community" , `alone`="alone")

# for relative abundance plots, e.g. of succession
barp.theme <- list(
  scale_y_continuous(labels = scales::percent, expand = c(0.01,0.01)),
  scale_fill_manual(values = species.palette, name = "Species", 
                    labels = c("Brachybacterium", "Brevibacterium", 
                               "Staph. xylosus", "Staph. equorum", 
                               "Scopulariopsis", "Penicillium", "Diutina")),
  labs(x = "days", y = "relative composition"),
  theme(panel.spacing.x = unit(0.5, "lines"),
        plot.title = element_text(size = 15), 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 12),
        plot.subtitle=element_text(size=13))
)

save(species.palette, names, names.italic, barp.theme, file = here("figures/plotting_objects.RData"))
