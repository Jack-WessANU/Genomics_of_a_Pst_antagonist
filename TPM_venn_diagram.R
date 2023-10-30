library(dplyr)

install.packages(c("tidyverse","hrbrthemes","tm","proustr"))
install.packages("VennDiagram")

countmatrix_abundance<-data.frame(Inplanta1$abundance,Inplanta2$abundance, Inplanta3$abundance, Invitro1$abundance,Invitro2$abundance, Invitro3$abundance)
row.names(countmatrix_abundance)<-row.names(Inplanta1$abundance)
plotMDS(countmatrix)

View(countmatrix_abundance)
countmatrixinvitro <- rowSums(countmatrix_abundance[,c(4,5,6)])
countmatrixinvitro_mean <- countmatrixinvitro[]/3
countmatrixinvitro_mean <- as.data.frame(countmatrixinvitro_mean)
countmatrixinvitro_expressed_100 <- countmatrixinvitro_mean %>% filter(countmatrixinvitro_mean > 1)

View(countmatrixinvitro_mean)

countmatrixinplanta <- rowSums(countmatrix_abundance[,c(1,2,3)])
countmatrixinplanta_mean <- countmatrixinplanta[]/3
countmatrixinplanta_mean <- as.data.frame(countmatrixinplanta_mean)
countmatrixinplanta_expressed_100 <- countmatrixinplanta_mean %>% filter(countmatrixinplanta_mean > 1)

View(countmatrixinplanta_mean)
inplantaonly <- countmatrixinplanta_expressed_100[!(countmatrixinplanta_expressed_100 %in% countmatrixinvitro_expressed_100)]

View(inplantaonly)
#Graphing venn Diagrams

library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)

# Load library
library(VennDiagram)

# Generate 3 sets of 200 words
invitro_explist <- row.names(countmatrixinvitro_expressed_100)
inplanta_explist <- row.names(countmatrixinplanta_expressed_100)

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- c("#2171B5", "#D94701")

# Chart
Venny <- venn.diagram(
  x = list(invitro_explist, inplanta_explist), category.names = c('axenic
growth','hyperparasitic
    growth'), filename = 'penx_expressed_genes_venn_diagram.png',output=TRUE, ext.text = FALSE,
  
  #Titles and such
  main = "Expressed genes
  Total = 10581", main.fontface = 0.5 ,main.fontfamily = "sans" ,main.cex = 0.9,
  
  euler.d = TRUE,
  scaled = TRUE,
  #Output features
  imagetype="png" ,
  height = 1100 , 
  width = 1100 , 
  resolution = 300,
  compression = "lzw",
  
  #Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  #Numbers
  cex = .9,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.7,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-45, 45),
  cat.dist = c(0.035, 0.035),
  cat.fontfamily = "sans",
  # rotation = 1
  
)

Venny
