setwd("/Users/nicolemejia/Documents/WTSH_peerj/boxplots_normalized_counts")
install.packages("reshape2")
library(reshape2)
library(ggplot2)
#read the data
ankrd <- read.csv("WTSH_boxplot_gene_ankrd - Sheet1.csv", header = TRUE)
hsp <- read.csv("WTSH_boxplot_HSPH1_counts - Sheet1.csv", header = TRUE)
#reshape the data so it's in the long form
data_long <- melt(ankrd, variable.name = "Body_Weight_Category", value.name = "Standardized_Count")
#make the plot using ggplot
ank <- ggplot(data_long, aes(x = Body_Weight_Category, y = Standardized_Count, na.rm = TRUE)) 
ank2 <- ank + geom_boxplot() + labs(title = "Gene Ankrd11_1", x = "Group",y = "Normalized Count") + theme_minimal()
print(ank3)
ank3 <- ank2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_text(face = "bold", size = 14), axis.title.y = element_text(face = "bold", size = 14), 
                     axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
print(ank3)
#add points
ank3 + geom_jitter(width = 0.2, size = 3, color = "black")                                                                                                     

# Do the same for the other gene
data_long2 <- melt(hsp, variable.name = "Body_Weight_Category", value.name = "Standardized_Count")
hsp <- ggplot(data_long2, aes(x = Body_Weight_Category, y = Standardized_Count, na.rm = TRUE))
hsp2 <- hsp + geom_boxplot() + labs(title = "Gene Hsph1", x = "Group",y = "Normalized Count") + theme_minimal()
hsp3 <- hsp2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_text(face = "bold", size = 14), axis.title.y = element_text(face = "bold", size = 14)
                     , axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13) )
hsp3 + geom_jitter(width = 0.2, size = 3, color = "black")
print(hsp3)
