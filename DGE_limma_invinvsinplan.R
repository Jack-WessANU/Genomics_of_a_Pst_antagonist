library(edgeR)
library(limma)
library(tximport)

Inplanta1 <-tximport("./salmon_pstp_wheat_both_cds_removed/salmon_cds_wheat_pstp_removed_inpla1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
Inplanta2 <-tximport("./salmon_pstp_wheat_both_cds_removed/salmon_cds_wheat_pstp_removed_inpla2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
Inplanta3 <-tximport("./salmon_pstp_wheat_both_cds_removed/salmon_cds_wheat_pstp_removed_inpla3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
Invitro1 <-tximport("./salmon_pstp_wheat_both_cds_removed/salmon_cds_wheat_pstp_removed_invin1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
Invitro2 <-tximport("./salmon_pstp_wheat_both_cds_removed/salmon_cds_wheat_pstp_removed_invin2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
Invitro3 <-tximport("./salmon_pstp_wheat_both_cds_removed/salmon_cds_wheat_pstp_removed_invin3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)

#combining the Salmon count information into a single matrix
countmatrix<-data.frame(Inplanta1$counts,Inplanta2$counts, Inplanta3$counts, Invitro1$counts,Invitro2$counts, Invitro3$counts)
row.names(countmatrix)<-row.names(Inplanta1$counts)
plotMDS(countmatrix)



#import design matrix (generated seperately in excel)
design<-read.csv("./design.csv")


#filtering genes by expression
dge <- DGEList(counts=countmatrix)
group <- as.factor(c("inplanta","inplanta","inplanta","invitro","invitro","invitro"))
keep <- filterByExpr(dge, group=group)
dge <- dge[keep,,keep.lib.sizes=FALSE]

#checking the results of filtering
dim(dge)
dim(design)


#Use voom to determine weights for each gene to be passed to limma, voom normalises the data
v<-voom(dge, design, plot=TRUE, normalize="quantile")

#fit linear model for each gene
fit <- lmFit(v,design)

#define the comparisons you want to examine
contrast.matrix <- makeContrasts(invitro - inplanta, levels = design)

#estimate contrasts for each gene
fit2 <- contrasts.fit(fit, contrast.matrix)

#Smoothing of standard error
fit2 <- eBayes(fit2)

#Generate results tables
pst_wheat_removed_penx_invitrovsinplanta <- topTable(fit2,coef=1, number=Inf)

pst_wheat_removed_penx_invitrovsinplanta

#Export transcriptome profiles
write.csv(pst_wheat_removed_penx_invitrovsinplanta, file = "pst_wheat_removed_penx_invitrovsinplanta.csv")

pst_wheat_removed_penx_invitrovsinplanta

## inplantavsinvitro is a file that may then be used to generate graphs to display data, such as the Volcano plot


