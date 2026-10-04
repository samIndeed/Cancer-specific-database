## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
conflicted::conflict_prefer("filter", "dplyr")

## ----createConfig-------------------------------------------------------------
library(prolfqua)

# ## ----LoadDataAndConfigure-----------------------------------------------------
# xx <- prolfqua::sim_lfq_data_protein_config(Nprot = 100)
# xx
# 



# Load data
cm=read.csv("~/Projects/RadNet/GlasgowCancer/Proteomics/NewDataCountMatrix/NewDataCountMatrix.csv",row.names =1)
cm=read.csv("~/Projects/RadNet/GlasgowCancer/Proteomics/NewDataCountMatrix/NewDataUniqueCountMatrix.csv",row.names =1)

rownames(cm)=cm$Genes
cm=cm[,-1]



cm[,"P006_T02"]=(cm[,"P006_T00"]+cm[,"P006_T06"])/2# interpolate for P006-T02



# cm0=cm[,grep("T00",colnames(cm))]

# cm0=cm0[-1,]

# cm0 <- cm0 %>% replace(is.na(.), 0)
# cm <- cm %>% replace(is.na(.), 0)

# df=data.frame(cm0)
df=data.frame(cm)

# df$protein_Id=rownames(cm0)
df$protein_Id=rownames(cm)

library("stringr")

# pat=str_replace_all(colnames(cm0),"_.*","")
pat=str_replace_all(colnames(cm),"_.*","")
patoutcome=str_replace_all(pat,c("P006"="Disease Free","V001" = "Disease Free", "V003" = "Disease Free", "P002" = "Cancer death/Metastasis","V002" = "Disease Free", "V015" = "Disease Free", "V013" = "Disease Free","V009" = "Disease Free", "P003" = "Cancer death/Metastasis"))
patoutcome2=str_replace_all(patoutcome,c("Disease Free"="A","Cancer death/Metastasis"="B"))

patoutcome2=str_replace_all(colnames(cm),".*_","")

# annot <- data.frame(Sample = colnames(cm0), Group = patoutcome2)
annot <- data.frame(Sample = colnames(cm), Group = patoutcome2,Subject=pat)


# convert into long format
table_long <- tidyr::pivot_longer(df, contains("_T"),names_to = "Sample", values_to = "Intensity")

table_long <- dplyr::inner_join(annot, table_long)

head(table_long)
# create TableAnnotation and AnalysisConfiguration

atable <- prolfqua::AnalysisTableAnnotation$new()
atable$fileName = "Sample"
atable$workIntensity = "Intensity"
atable$hierarchy[["protein_Id"]]    <-  "protein_Id"

# atable$hierarchy[["Genes"]]    <-  "Genes"
atable$factors[["Group"]] <- "Group"
atable$factors[["Subject"]] <- "Subject"

config <- prolfqua::AnalysisConfiguration$new(atable)

# Build LFQData object
analysis_data <- prolfqua::setup_analysis(table_long, config)
lfqdata <- prolfqua::LFQData$new(analysis_data, config)
lfqdata$hierarchy_counts()

## ----removeSmallIntensities---------------------------------------------------
# lfqdata <- prolfqua::LFQData$new(xx$data, xx$config)
lfqdata$remove_small_intensities()
lfqdata$factors()

## ----showWide, eval = TRUE----------------------------------------------------
lfqdata$to_wide()$data[1:3,1:7]

## ----getPlotter---------------------------------------------------------------
lfqplotter <- lfqdata$get_Plotter()
density_nn <- lfqplotter$intensity_distribution_density()

## ----makeMissingHeatmap, fig.cap="Heatmap where missing proteins (zero in case of MaxQuant reported intensities), black - missing protein intensities, white - present"----
lfqplotter$NA_heatmap()

## ----missignessPerGroup, fig.cap="# of proteins with 0,1,...N missing values"----
lfqdata$get_Summariser()$plot_missingness_per_group()

## ----missignessHistogram, fig.cap="Intensity distribution of proteins depending on # of missing values"----
lfqplotter$missigness_histogram()

## ----PlotCVDistributions, fig.cap="Violin plots of CVs in the different groups and among all groups"----
stats <- lfqdata$get_Stats()
stats$violin()
prolfqua::table_facade( stats$stats_quantiles()$wide, paste0("quantile of ",stats$stat ))

## ----plotCVsplitbyMedianIntensity, fig.cap="Distribution of CV's for top 50% and bottom 50% proteins by intensity."----
stats$density_median()


## ----normalizedata------------------------------------------------------------
lt <- lfqdata$get_Transformer()
transformed <- lt$log2()$robscale()$lfq
transformed$config$table$is_response_transformed


## ----genplotNorm, fig.cap="Normalized intensities."---------------------------
pl <- transformed$get_Plotter()
density_norm <- pl$intensity_distribution_density()

## ----showIntensityDistributions, fig.cap="Distribution of intensities before and after normalization."----
gridExtra::grid.arrange(density_nn, density_norm)

## ----plotScatterMatrix, fig.cap = "Scatterplot matrix"------------------------
pl$pairs_smooth()

## ----createHeatmap------------------------------------------------------------
p <- pl$heatmap_cor()

## ----plotHeatmap, fig.cap="Heatmap, Rows - proteins, Columns - samples", fig.align=5, fig.height=5----
p

## ----lookatfactors------------------------------------------------------------
transformed$factors()


## ----defineModelAndContrasts--------------------------------------------------
# formula_Condition <-  strategy_lm("transformedIntensity ~ group_")
formula_Condition <-  strategy_lm("transformedIntensity ~ Subject+Group")

# specify model definition
modelName  <- "Model"
# Contrasts <- c("T06vsT00" = "GroupT06 - GroupT00",
#                "T02vsT00" = "GroupT02 - GroupT00")
# Contrasts <- c("AvsB" = "GroupA - GroupB")
Contrasts <- c("AvsB" = "GroupT06 - GroupT00")


## ----buildModel---------------------------------------------------------------
mod <- prolfqua::build_model(
  transformed$data,
  formula_Condition,
  subject_Id = transformed$config$table$hierarchy_keys() )


## ----showANOVA, fig.cap="Distribtuion of adjusted p-values (FDR)"-------------
mod$anova_histogram("FDR")

## ----filterDataForFDR---------------------------------------------------------
aovtable <- mod$get_anova()
head(aovtable)
dim(aovtable)
xx <- aovtable |> dplyr::filter(FDR < 0.1)
signif <- transformed$get_copy()
signif$data <- signif$data |> dplyr::filter(protein_Id %in% xx$protein_Id)
hmSig <- signif$get_Plotter()$heatmap()

aovtable[aovtable$p.value<0.00001,]

## ----showSigHeatmap, fig.cap="Heatmap for proteins with FDR < 0.2 in the analysis of variance"----
hmSig


## ----computeContrasts---------------------------------------------------------
contr <- prolfqua::Contrasts$new(mod, Contrasts)
v1 <- contr$get_Plotter()$volcano()

## ----computeContrastsModerated------------------------------------------------
contr <- prolfqua::ContrastsModerated$new(contr)
contrdf <- contr$get_contrasts()

## ----plotVolcanos, fig.cap="Volcano plot, Left panel - no moderation, Right panel - with moderation."----
plotter <- contr$get_Plotter()
v2 <- plotter$volcano()
gridExtra::grid.arrange(v1$FDR,v2$FDR, ncol = 1)


## ----showMAplot, fig.cap="MA plot showing the dependency of mean abuncance with respect to the difference"----
plotter$ma_plotly()

## ----checkProteinsInOutput----------------------------------------------------
#myProteinIDS <- c("sp|Q12246|LCB4_YEAST",  "sp|P38929|ATC2_YEAST",  "sp|Q99207|NOP14_YEAST")
myProteinIDS <- c("RIGI",  "HLA-F",  "A_HLA")
dplyr::filter(contrdf, protein_Id %in% myProteinIDS)


EnhancedVolcano(contrdf,
                lab = contrdf$protein_Id,
                # selectLab= grep("HLA",plotter$contrastDF$protein_Id,value = T),
                x = 'diff',
                y = 'p.value',legendPosition = 'none',
                pCutoff =0.001,
                subtitle = "Proteotypic and Razor peptides",
                title = "Disease-free vs Non Disease-free")

## ----computeMissing, eval=TRUE------------------------------------------------
mC <- ContrastsMissing$new(lfqdata = transformed, contrasts = Contrasts)
colnames(mC$get_contrasts())


## ----mergeResults, fig.cap="Volcano plots for the two contrasts with missing value imputation from the group_average model."----

merged <- prolfqua::merge_contrasts_results(prefer = contr,add = mC)$merged
plotter <- merged$get_Plotter()
tmp <- plotter$volcano()
tmp$FDR
tmp$p.value
plotter$contrastDF[plotter$contrastDF$p.value<0.001,]$protein_Id
plotter$contrastDF[plotter$contrastDF$protein_Id%in%"RIGI",]$p.value
plotter$contrastDF[plotter$contrastDF$protein_Id%in%"A_HLA",]$p.value
plotter$contrastDF[plotter$contrastDF$protein_Id%in%"HLA-F",]$p.value
plotter$contrastDF[plotter$contrastDF$protein_Id%in%grep("HLA",plotter$contrastDF$protein_Id,value = T),]$p.value


EnhancedVolcano(plotter$contrastDF,
                lab = plotter$contrastDF$protein_Id,
                # selectLab= grep("HLA",plotter$contrastDF$protein_Id,value = T),
                x = 'diff',
                y = 'p.value',legendPosition = 'none',
                pCutoff =0.001,
                subtitle = "Proteotypic and Razor peptides",
                title = "Disease-free vs Non Disease-free")


library(EnhancedVolcano)
pdf("~/Projects/RadNet/GlasgowCancer/Proteomics/CompareTimePoints/VolcProteo.pdf")
EnhancedVolcano(plotter$contrastDF,
                lab = plotter$contrastDF$protein_Id,
                # selectLab= grep("HLA",plotter$contrastDF$protein_Id,value = T),
                x = 'diff',
                y = 'p.value',legendPosition = 'none',
                pCutoff =0.001,
                subtitle = "Proteotypic peptides",
                title = "T06 vs T00")

EnhancedVolcano(plotter$contrastDF,
                lab = plotter$contrastDF$protein_Id,
                selectLab= grep("HLA",plotter$contrastDF$protein_Id,value = T),
                x = 'diff',
                y = 'p.value',legendPosition = 'none',
                pCutoff =0.001,
                subtitle = "Proteotypic peptides",
                title = "T06 vs T00")
dev.off()

pdf("~/Projects/RadNet/GlasgowCancer/Proteomics/CompareTimePoints/VolcProteoAndRaz.pdf")
EnhancedVolcano(plotter$contrastDF,
                lab = plotter$contrastDF$protein_Id,
                # selectLab= grep("HLA",plotter$contrastDF$protein_Id,value = T),
                x = 'diff',
                y = 'p.value',legendPosition = 'none',
                pCutoff =0.001,
                subtitle = "Proteotypic and Razor peptides",
                title = "T06 vs T00")

EnhancedVolcano(plotter$contrastDF,
                lab = plotter$contrastDF$protein_Id,
                selectLab= grep("HLA",plotter$contrastDF$protein_Id,value = T),
                x = 'diff',
                y = 'p.value',legendPosition = 'none',
                pCutoff =0.001,
                subtitle = "Proteotypic and Razor peptides",
                title = "T06 vs T00")
dev.off()

plotter$contrastDF
## ----mergedMore---------------------------------------------------------------
merged <- prolfqua::merge_contrasts_results(prefer = contr,add = mC)

moreProt <- transformed$get_copy()
moreProt$data <- moreProt$data |> dplyr::filter(protein_Id %in% merged$more$contrast_result$protein_Id)
moreProt$get_Plotter()$raster()

# here we do not get anything because there is nothing imputed!

## ----prepForGSEA--------------------------------------------------------------
#evalAll <- require("clusterProfiler") & require("org.Sc.sgd.db") & require("prora")
# evalAll <- require("clusterProfiler") & require("org.Sc.sgd.db2") & require("prora")


# ----clusterProfiler, eval=evalAll--------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
# library(prora)

bb <- prolfqua::get_UniprotID_from_fasta_header(merged$merged$get_contrasts(),
                                             idcolumn = "protein_Id")

write.csv(bb,"~/Projects/RadNet/GlasgowCancer/Proteomics/CompareTimePoints/T06vsT00Proteo.csv")
write.csv(bb,"~/Projects/RadNet/GlasgowCancer/Proteomics/CompareTimePoints/T06vsT00ProteoAndRaz.csv")

# bb <- prora::map_ids_uniprot(bb)
ranklist <- bb$statistic
# names(ranklist) <- bb$P_ENTREZGENEID
names(ranklist) <- bb$protein_Id
library(biomaRt)

resn <- getBM(attributes = c(
                             'external_gene_name',"entrezgene_id"),
             filters = 'external_gene_name', 
             values = bb$protein_Id,
             mart = mart)
# resn

bbb=merge(bb,resn,by.x="protein_Id",by.y="external_gene_name")

bbbb=bbb[!(is.na(bbb$entrezgene_id)),]

ranklist <- bbbb$statistic
names(ranklist)=bbbb$entrezgene_id
# names(ranklist) <- resn[!(is.na(resn$entrezgene_id)),]$entrezgene_id

bbbbb=bbbb[!duplicated(names(ranklist)),]


# resnn=resn[!(is.na(resn$entrezgene_id)),]
# resnnn=resnn[!duplicated(names(ranklist)),]

ranklist <- bbbbb$statistic
names(ranklist)=bbbbb$entrezgene_id

# names(ranklist) <- resnnn$entrezgene_id


res <- clusterProfiler::gseGO(
  sort(ranklist, decreasing = TRUE),
  OrgDb = org.Hs.eg.db,
  ont = "ALL",eps=0)

bbbbb[names(ranklist)%in%c("10159","4179","4345","55612"),]

bbbbb[bbbbb$protein_Id%in%c("RIGI"),]


# RIGI not signif because all NAs...

grep("10159",names(ranklist))
## ----ridgeplot, fig.cap="ridgeplot", eval = evalAll---------------------------
ridgeplot( res )

## ----dotplot, fig.cap = "Dotplot", eval = evalAll-----------------------------
pdf("~/Projects/RadNet/GlasgowCancer/Proteomics/CompareTimePoints/GSEADotProteo.pdf")
dotplot(res , showCategory = 16)
dev.off()
## ----upsetplot, fig.cap="Upset Plot", eval = evalAll--------------------------
enrichplot::upsetplot(res,n=15)

res@result[res@result$Description%in%"negative regulation of non-canonical NF-kappaB signal transduction",]
res@result[res@result$Description%in%"clathrin vesicle coat",]
str_split(res@result[res@result$Description%in%"negative regulation of non-canonical NF-kappaB signal transduction",]$core_enrichment,"/")[[1]]

bbbbb[names(ranklist)%in%str_split(res@result[res@result$Description%in%"single-stranded DNA binding",]$core_enrichment,"/")[[1]],]
bbbbb[names(ranklist)%in%str_split(res@result[res@result$Description%in%"complement activation",]$core_enrichment,"/")[[1]],]
bbbbb[names(ranklist)%in%str_split(res@result[res@result$Description%in%"ribosome",]$core_enrichment,"/")[[1]],]
bbbbb[names(ranklist)%in%str_split(res@result[res@result$Description%in%"extracellular matrix",]$core_enrichment,"/")[[1]],1]

pdf("~/Projects/RadNet/GlasgowCancer/Proteomics/CompareTimePoints/VolcECMProteo.pdf")
EnhancedVolcano(plotter$contrastDF,
                lab = plotter$contrastDF$protein_Id,
                selectLab= bbbbb[names(ranklist)%in%str_split(res@result[res@result$Description%in%"extracellular matrix",]$core_enrichment,"/")[[1]],1],
                x = 'diff',
                y = 'p.value',legendPosition = 'none',
                pCutoff =0.0001,
                subtitle = "ECM - Proteotypic peptides",
                title = "T06 vs T00")
dev.off()
pdf("~/Projects/RadNet/GlasgowCancer/Proteomics/CompareTimePoints/VolcRibosomeProteo.pdf")
EnhancedVolcano(plotter$contrastDF,
                lab = plotter$contrastDF$protein_Id,
                selectLab= bbbbb[names(ranklist)%in%str_split(res@result[res@result$Description%in%"ribosome",]$core_enrichment,"/")[[1]],1],
                x = 'diff',
                y = 'p.value',legendPosition = 'none',
                pCutoff =0.0001,
                subtitle = "Ribosome - Proteotypic peptides",
                title = "T06 vs T00")
dev.off()

write.csv(bbbbb[names(ranklist)%in%str_split(res@result[res@result$Description%in%"negative regulation of non-canonical NF-kappaB signal transduction",]$core_enrichment,"/")[[1]],],
          "~/Projects/RadNet/GlasgowCancer/Proteomics/CompareTimePoints/negative regulation of non-canonical NF-kappaB signal transduction_ProtAndRaz.csv")



## -----------------------------------------------------------------------------
sessionInfo()

