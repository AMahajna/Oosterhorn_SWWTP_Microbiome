#To do: 
# 1. bicluster according to 10.4.1 
# cluster is pathway class and heatmap is sample id. 
################################################################################
##Create folder for project organization

if(!dir.exists("input_data")){dir.create("input_data")}
if(!dir.exists("output_data")){dir.create("output_data")}
if(!dir.exists("figures")){dir.create("figures")}
if(!dir.exists("scripts")){dir.create("scripts")}

################################################################################
##load packages 

source(file = "scripts/install_load_packages.r")

################################################################################
##Load data

#Microbiome Data for assay and rowData in TSE 
table_raw <- read_csv(file = "input_data/4717_groupby_filtered_all_summary_pathways_taxonomy_abs.csv", show_col_types = FALSE)

#Relevant process data to be added to column data in TSE
process_data <- read_csv(file = "input_data/relevant_process_data.csv", show_col_types = FALSE)

#Check the class of each column to be numeric 
#sapply(process_data,class)

#Sample data overview added to column data in TSE
samdat = read_excel("input_data/WWTP_overview_samples_20180125.xlsx") %>%
  filter(NGS_performed == "YES") %>% 
  dplyr::select(sampleid = MonsterCodeBC, everything()) 

#removal_efficiency respective to sample to be added to column data in TSE
removal_efficiency  = read_excel("input_data/Removal_Efficiency.xlsx")

################################################################################
##Tidying and subsetting data

################################################################################
#Assay: these are the reads/counts for the microbiome
#Merge uniprot name and code to have unique names 
table_raw$uniprot_species_name = paste(table_raw$species_strain,table_raw$uniprot_name,table_raw$uniprot_acc, sep="_")

#change the sample name in column names in the assay data to match that 
#in the sample data overview
assay_data <-
  table_raw %>%  
  dplyr::select(uniprot_species_name, starts_with("47"))%>% 
  rename_with(~str_sub(., end = 8))

colnames(assay_data)[1] <- "uniprot_species_name"

#Cleaning "<DL" values= replace less than detection level with zero
assay_data =  data.frame(lapply(assay_data, function(x) {
  gsub("<DL", 0, x)
}))

#Checking the class of the values in the assay data
#sapply(assay_data,class)

#Changing the class of the assay data to numeric
assay_data[,2:33] <- sapply(assay_data[,2:33],as.numeric)

#Checking class and dimensionality
#sapply(assay_data,class)
#dim(assay_data)

#convert assay data to datafrme format 
assay_data <- 
  assay_data %>% 
  as_tibble() %>% 
  rename_with(~gsub("X","",.x)) %>% 
  as.data.frame()



#make row name the unique entity is for row data that contains gene data 
rownames(assay_data) = assay_data$uniprot_species_name
#Here we take only numeric data and convert to matrix for assay data
#The unique entity data is row names in the matrix 
assay_df = assay_data
assay_data = assay_data[,-1] %>% as.matrix()

################################################################################
#taxonomic data _ goes in row data in TSE 
row_data_tax <-
  table_raw %>% 
  dplyr::select(uniprot_species_name,taxid:domain)

#change the taxonomic ranking name to match mia package 
names(row_data_tax)[names(row_data_tax) == 'sub_domain'] <- "kingdom"
names(row_data_tax)[names(row_data_tax) == 'species'] <- "species_missing"
names(row_data_tax)[names(row_data_tax) == 'species_strain'] <- "species"
#check
#colnames(row_data_tax)
#reorder taxonomic ordering 
row_data_tax_reorder = select(row_data_tax,uniprot_species_name, domain, kingdom, phylum, class, order, family, genus,species, everything())
#capitalize taxonomy name 
colnames(row_data_tax_reorder) <- c('uniprot_species_name', 'Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species',"taxid","species_missing")

################################################################################
#Functional data of enzymes _ goes in row data in TSE 
row_data_genes <-
  table_raw %>% 
  dplyr::select(uniprot_species_name,kegg_reaction_number:go_number) %>% 
  mutate_all(funs(str_replace(., "_","NA"))) %>% 
  mutate_all(funs(str_replace(., "-$","NA"))) 

################################################################################
#sample data _ goes in column data in TSE including the  

process_data <- subset(process_data, select = -rarity)
#Standardize environmental data - reflect on standardization method choice 
process_data <- decostand(process_data, method ="log")

samdat <- samdat  %>%
  mutate(Season = case_when(
    month(Sample_Date) %in% c(3, 4, 5) ~ "Spring",
    month(Sample_Date) %in% c(6, 7, 8) ~ "Summer",
    month(Sample_Date) %in% c(9, 10, 11) ~ "Fall",
    month(Sample_Date) %in% c(12, 1, 2) ~ "Winter"
  ))

#for coloring in groups and clustering 
samdat$Year_Sample<-as.character(samdat$Year_Sample)
samdat$Sample_Date<-as.character(samdat$Sample_Date)

samdat = data.frame(cbind(samdat,process_data,(removal_efficiency[ ,3:6])/100))

#check dim of samdat it should be: 32*33
#dim(samdat)

################################################################################
##plot missing data for both gene and taxa

missing_tax_plot = row_data_tax %>% 
  dplyr::select('domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus','species_missing','species')%>% 
  rename("spcies_strain"="species","species"="species_missing")

png(filename="figures/missing_taxonomy.png" ,units = 'in',width=9, height=6, res=1000)
plot_missing(missing_tax_plot, title = "missing data profile for taxonomy")
dev.off()

missing_genes_plot <- row_data_genes %>%
  select(-uniprot_species_name)%>%
  mutate_all(~ifelse(. == "NA", NA, 1))

png(filename="figures/missing_functions.png" ,units = 'in',width=9, height=6, res=1000)
plot_missing(missing_genes_plot, title = "missing data profile for functional data")
dev.off()

################################################################################
#Multi-Assay Experiment 
# 1.  "kegg_type_pathway"  
# 2.  "kegg_pathway" 

#The choice of the data structure as multi-assay experiment sperms from the fact
#that missing_data in the functional genes dataset is too high to use reliably 
#thus, the pathway data gets grouped and these is dimensionality reduction 

#Note: KEGG metabolic pathways are relevant for the process 


pathway_class = row_data_genes %>% 
  dplyr::select(kegg_class_pathway,kegg_type_pathway, kegg_pathway)
#check chould be TRUE FALSE 
#dim(assay_df) == dim(pathway_class)
pathway_class = cbind(pathway_class,assay_df)
rownames(pathway_class) <- NULL

#Extract Metabolism related pathway data 
pathway_metabolism =  pathway_class[pathway_class$kegg_class_pathway == "A09100 Metabolism", ]

#Group by metabolism pathway and metabolism pathway type and sum reads respectively  
pathway_metabolism <- pathway_metabolism %>%
  group_by(kegg_pathway,kegg_type_pathway, kegg_class_pathway) %>%
  summarise(across(where(is.numeric), sum)) %>%
  arrange(kegg_type_pathway)

# We sum up all other metabolic pathways
pathway_other = colSums(pathway_class[pathway_class$kegg_class_pathway != "A09100 Metabolism", -c(1, 2,3)])
pathway_other$kegg_class_pathway ="Other pathway class"
pathway_other$kegg_type_pathway = "Other pathway type"
pathway_other$kegg_pathway ="Other pathway "
pathway_other %>% 
  dplyr::select(kegg_pathway,kegg_type_pathway, kegg_class_pathway,  everything())%>% 
  data.frame()

#We add it to total assay 
pathway_metabolism_other = rbind(pathway_metabolism,pathway_other)

#Check the sum of reads in the new dataframe 
sum(colSums(pathway_metabolism_other[ ,-c(1,2,3)]) == colSums(assay_df)) 
#should equal 32 n samples = if it's equal to 32 meaning grouping went right

pathway_assay = as.data.frame(pathway_metabolism_other[ ,4:ncol(pathway_metabolism_other)]) 
rownames(pathway_assay) = pathway_metabolism_other$kegg_pathway

#Check 
length(colnames(pathway_assay))
#should be 32 
nrow(pathway_assay) == nrow(pathway_metabolism_other$kegg_pathway)
#Should be true 

#set rownames to pathway name to convert into assay for TSE_pathway
#rownames(pathway_assay) = pathway_metabolism_other$kegg_pathway
row_data_pathway = pathway_metabolism_other[ ,1:2]
#row_data_pathway = as.data.frame(row_data_pathway)
#rownames(row_data_pathway) = row_data_pathway$kegg_pathway
#row_data_pathway = row_data_pathway[ ,2]
#pathway_type_df <- t(pathway_type[, -1])  
#colnames(pathway_type_df) <- pathway_type$kegg_type_pathway

################################################################################
#Assay for pathway_type

#Group by metabolism type
pathway_metabolism_type <- pathway_metabolism %>%
  group_by(kegg_type_pathway) %>%
  summarise(across(where(is.numeric), sum)) %>%
  arrange(kegg_type_pathway)

#We add it to total assay 
pathway_other = pathway_other %>%
  as.data.frame()%>%
  dplyr::select(-kegg_class_pathway,-kegg_pathway)%>%
  select(kegg_type_pathway, everything())


pathway_metabolism_type_other = rbind(pathway_metabolism_type,pathway_other)

#Check the sum of reads in the new dataframe 
sum(colSums(pathway_metabolism_other[ ,-c(1,2,3)]) == colSums(assay_df)) 
#should equal 32 n samples = if it's equal to 32 meaning grouping went right

pathway_type_assay = as.data.frame(pathway_metabolism_other[ ,4:ncol(pathway_metabolism_other)]) 
rownames(pathway_type_assay) = pathway_metabolism_type_other$kegg_type_pathway

#Check 
length(colnames(pathway_assay))
#should be 32 
nrow(pathway_assay) == nrow(pathway_metabolism_other$kegg_pathway)
#Should be true 

#set rownames to pathway name to convert into assay for TSE_pathway
#rownames(pathway_assay) = pathway_metabolism_other$kegg_pathway
row_data_pathway = pathway_metabolism_other[ ,1:2]
#row_data_pathway = as.data.frame(row_data_pathway)
#rownames(row_data_pathway) = row_data_pathway$kegg_pathway
#row_data_pathway = row_data_pathway[ ,2]
#pathway_type_df <- t(pathway_type[, -1])  
#colnames(pathway_type_df) <- pathway_type$kegg_type_pathway




################################################################################














#Check the grouping is correct 
#sum(colSums(pathway_type[ ,-1]) == colSums(assay_df[ , -1])) should equal 32 n samples

#pathway = row_data_genes %>% 
#  dplyr::select(kegg_pathway)
#pathway = cbind(pathway,assay_df[ ,-1])

#pathway <- pathway %>%
#  group_by(kegg_pathway) %>%
#  summarise(across(where(is.numeric), sum))

#pathway_assay = as.data.frame(pathway[ ,2:33])
#rownames(pathway_assay) = pathway$kegg_pathway

#pathway_df <- t(pathway[, -1])  
#colnames(pathway_df) <- pathway$kegg_pathway

#Check the grouping is correct 
#sum(colSums(pathway_type[ ,-1]) == colSums(assay_df[ , -1])) should equal 32 n samples
#taxonomyRanks(tse)

################################################################################
#Create TSE 

assays = SimpleList(counts = assay_data)
colData = data.frame(samdat)
rowData = data.frame(row_data_tax_reorder)

tse<- TreeSummarizedExperiment(assays = assays,
                               colData = colData,
                               rowData = rowData
)

#taxonomyRanks(tse)
################################################################################
#tax limited TSE not including the enzymes 
tse_species <- mergeFeaturesByRank(tse, "Species", na.rm=TRUE)

################################################################################
#Phylum level limited TSE 
tse_phylum <- mergeFeaturesByRank(tse, "Phylum", na.rm=TRUE)

################################################################################
#pathway TSE
#pathway_assay = (counts = pathway_assay)
#colData = data.frame(samdat)
#colData$Year_Sample<-as.character(colData$Year_Sample)

tse_pathway <- TreeSummarizedExperiment(assays = SimpleList(counts = as.matrix(pathway_assay)) ,
                                        colData = data.frame(samdat),
                                        rowData = row_data_pathway)

################################################################################
#pathway_type_assay = (counts = pathway_type_assay)
#colData = data.frame(samdat)
#colData$Year_Sample<-as.character(colData$Year_Sample)

tse_pathway_type <- TreeSummarizedExperiment(assays = SimpleList(counts = as.matrix(pathway_type_assay)) ,
                                        colData = data.frame(samdat),
                                        rowData = DataFrame(pathway_type[ ,1]))


################################################################################
#Relative Abundance Assays
tse <- transformAssay(tse, method = "relabundance")
tse_species <- transformAssay(tse_species, method = "relabundance")
tse_phylum <- transformAssay(tse_phylum, method = "relabundance")
tse_pathway <- transformAssay(tse_pathway, method = "relabundance")
tse_pathway_type <- transformAssay(tse_pathway_type, method = "relabundance")

################################################################################
#1. Can group by pathway and create a different TSE and analyze as MAE
# Sum all samples, and then group by Phylum and pathway_type 
#The do cross-correlation of other plots
# Create an ExperimentList that includes experiments
experiments <- ExperimentList(microbiome = tse,
                              tax = tse_species, 
                              phylum = tse_phylum, 
                              pathways = tse_pathway,
                              type = tse_pathway_type)

# Create a MAE
mae <- MultiAssayExperiment(experiments = experiments)
################################################################################
################################################################################
#Exploration 

################################################################################
#Abundance Denisty: jitter plot 
png(filename="figures/abundance_phylum.png" ,units = 'in',width=9, height=6, res=1000)
plotAbundanceDensity(tse_phylum, layout = "jitter", assay.type = "relabundance",
                     n = 40, colour_by="Sample_Date", point_size=2, point_shape=19) + 
  scale_x_log10(label=scales::percent)
dev.off()

png(filename="figures/abundance_species.png" ,units = 'in',width=9, height=6, res=1000)
plotAbundanceDensity(tse_species, layout = "jitter", assay.type = "relabundance",
                     n = 40, colour_by="Sample_Date", point_size=2, point_shape=19) + 
  scale_x_log10(label=scales::percent)
dev.off()

#We see that the core community is limited and most of the species/phylum are
#less prevalent in our community 

################################################################################
#Abundance Density: density plot 

#The relative abundance values for the top-5 taxonomic features can be visualized
#as a density plot over a log scaled axis,
png(filename="figures/density_species.png" ,units = 'in',width=9, height=6, res=1000)
plotAbundanceDensity(tse_species, layout = "density", assay.type = "relabundance",
                     n = 5, colour_by="Year_Sample", point_alpha=1/10) +
  scale_x_log10()
dev.off()

png(filename="figures/density_phylum.png" ,units = 'in',width=9, height=6, res=1000)
plotAbundanceDensity(tse_phylum, layout = "density", assay.type = "relabundance",
                     n = 5, colour_by="Year_Sample", point_alpha=1/10) +
  scale_x_log10()
dev.off()

#Abundance: of metabolic pathway and metabolic pathway type 
png(filename="figures/density_pathway.png" ,units = 'in',width=9, height=6, res=1000)
plotAbundanceDensity(tse_pathway, layout = "density", assay.type = "relabundance",
                     n = 10, colour_by="Year_Sample", point_alpha=1/10) + 
  scale_x_log10(label=scales::percent)
dev.off()

png(filename="figures/density_pathway_type.png" ,units = 'in',width=9, height=6, res=1000)
plotAbundanceDensity(tse_pathway_type, layout = "density", assay.type = "relabundance",
                     n = 10, colour_by="Year_Sample", point_alpha=1/10) + 
  scale_x_log10(label=scales::percent)
dev.off()
################################################################################
################################################################################
#Prevalence 

#Prevalence quantifies the frequency of samples where certain microbes were
#detected (above a given detection threshold).
################################################################################
#Core community on Phylum level
core_community_phylum = getPrevalence(tse, rank = "Phylum", detection = 0, prevalence = 99/100, sort = TRUE, as_relative = TRUE)
core_community_phylum = data.frame(core_community_phylum)
colnames(core_community_phylum) <- c("prevalence")
core_community_phylum$core_community_phylum = rownames(core_community_phylum)
core_community_phylum_names = core_community_phylum$core_community_phylum[core_community_phylum$prevalence ==  1]
print("Core community is present is ")
core_community_phylum_names


core_community_species = getPrevalence(tse, rank = "Species", detection = 0, prevalence = 99/100, sort = TRUE, as_relative = TRUE)
core_community_species = data.frame(core_community_species)
colnames(core_community_species) <- c("prevalence")
core_community_species$core_community_species = rownames(core_community_species)
core_community_species_names = core_community_species$core_community_species[core_community_species$prevalence ==  1]
print("on species level:")
core_community_species_names
#plot prevalence of different phylum across most samples-
#show prevalence of different phylum within the activated sludge microbiome of the activated sludge in SWWTP

################################################################################
#Plotting prevalence 
#The prevalence of different life domains across all samples
tse_domain <- mergeFeaturesByRank(tse, "Domain", na.rm=TRUE)
altExp(tse, "Domain") <- tse_domain

rowData(altExp(tse,"Domain"))$prevalence <- 
  getPrevalence(altExp(tse,"Domain"), detection = 0, prevalence =100/100, sort = FALSE,
                assay.type = "counts", as_relative = TRUE)

png(filename="figures/prevalence_domain.png" ,units = 'in',width=9, height=6, res=1000)
plotRowData(altExp(tse,"Domain"), "prevalence", colour_by = "Domain", point_size = 10)
dev.off()
#Here, I was trying to visualize the prevalence for all phylum present in 100% of 
#the samples. However, there are way too many phylum in the plot
#transposed_vector <- (t(core_community_phylum))
#colnames(transposed_vector) <- transposed_vector[2, ]
#transposed_vector <- transposed_vector[-2, ]
#rowData(altExp(tse, "Phylum"))$prevalence <- transposed_vector
#plotRowData(altExp(tse,"Phylum"), "prevalence", colour_by = "Phylum")
################################################################################
#Prevalence tree
altExps(tse) <- splitByRanks(tse)
altExps(tse) <-
  lapply(altExps(tse),
         function(y){
           rowData(y)$prevalence <- 
             getPrevalence(y, detection = 0, prevalence =100/100, sort = FALSE,
                           assay.type = "counts", as_relative = TRUE)
           y
         })

top_phyla_mean <- getTopFeatures(altExp(tse,"Phylum"),
                                 method="mean",
                                 top=5L,
                                 assay.type="counts")
x <- unsplitByRanks(tse, ranks = taxonomyRanks(tse)[1:6])
x <- addTaxonomyTree(x)

#The most abundance feature in the samples 
#Prevalence of top phyla as judged by mean abundance
png(filename="figures/abundant_phyla_prevalence_top5.png" ,units = 'in',width=9, height=6, res=1000)
plotRowTree(x[rowData(x)$Phylum %in% top_phyla_mean,],
                  edge_colour_by = "Phylum",
                  tip_colour_by = "prevalence",
                  node_colour_by = "prevalence")
dev.off() 

################################################################################
#Rare taxa 

# getRareFeatures returns the inverse
#rare <- getRareFeatures(tse,
#                        rank = "Phylum",
#                        detection = 0,
#                        prevalence = 10/100,
#                        as_relative = TRUE)
#head(rare)

# Gets a subset of object that includes rare taxa
#altExp(tse, "rare") <- subsetByRareFeatures(tse,
#                                            rank = "Class",
#                                            detection = 0.001,
#                                            prevalence = 0.001,
#                                            as_relative = TRUE)
#altExp(tse, "rare")     

################################################################################
################################################################################
#Quality Control

# Pick the top taxa
top_features <- getTopFeatures(tse_species, method="median", top=10)
# Check the information for these
rowData(tse_species)[top_features, taxonomyRanks(tse_species)]

#Library Size/Read Count 
#The total counts/sample can be calculated using perCellQCMetrics/addPerCellQC
#from the scater package. The former one just calculates the values,
#whereas the latter one directly adds them to colData.

perCellQCMetrics(tse)
tse <- addPerCellQC(tse)
perCellQCMetrics(tse_species)
tse_species <- addPerCellQC(tse_species)
#colData(tse)


p1 <- ggplot(as.data.frame(colData(tse_species))) +
  geom_histogram(aes(x = sum), color = "black", fill = "gray", bins = 30) +
  labs(x = "Library size", y = "Frequency (n)") + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), # Removes the grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) # Adds y-axis


df <- as.data.frame(colData(tse_species)) %>%
  arrange(sum) %>%
  mutate(index = 1:n())
p2 <- ggplot(df, aes(y = index, x = sum/1e6)) +
  geom_point() +  
  labs(x = "Library size (million reads)", y = "Sample index") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(), # Removes the grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) # Adds y-axis

#The distribution of calculated library sizes can be visualized as a histogram (left)
#or by sorting the samples by library size (right).

LibrarySize =p1 + p2

png(filename="figures/library_size.png" ,units = 'in',width=9, height=6, res=1000)
#pdf("figures/LibrarySize.pdf")
print(LibrarySize)
dev.off() 

################################################################################
################################################################################
#Community Similarity

#the environmental (process data has been normalized)

# Perform RDA
tse_species <- runRDA(tse_species,
              assay.type = "relabundance",
              formula = assay ~ INF_Cl_mg_per_l + INF_COD_mg_O2_per_l + INF_Nkj_mg_N_per_l + INF_PO4o_mg_P_per_l+  INF_SO4_µg_per_l + INF_TSS_mg_per_l + Glycerol_kg +Return_sludge_m3_per_h+ Inf_Flow_m3_per_h + Capacity_blowers_. +DW_AT_g_per_l+SVI_10 + T_avg_C,
              distance = "bray",
              na.action = na.exclude)

tse <- runRDA(tse,
              assay.type = "relabundance",
              formula = assay ~ INF_Cl_mg_per_l + INF_COD_mg_O2_per_l + INF_Nkj_mg_N_per_l + INF_PO4o_mg_P_per_l+  INF_SO4_µg_per_l + INF_TSS_mg_per_l + Glycerol_kg +Return_sludge_m3_per_h+ Inf_Flow_m3_per_h + Capacity_blowers_. +DW_AT_g_per_l+SVI_10 + T_avg_C,
              distance = "bray",
              na.action = na.exclude)
#full equatiom: overly redundant 
#formula = assay ~ INF_Cl_mg_per_l + INF_COD_mg_O2_per_l + INF_Nkj_mg_N_per_l + INF_PO4o_mg_P_per_l+  INF_SO4_µg_per_l + INF_TSS_mg_per_l + Glycerol_kg +Return_sludge_m3_per_h+ Inf_Flow_m3_per_h + Capacity_blowers_. + DW_AT_g_per_l + SVI_10 + T_avg_C,


# Store results of PERMANOVA test
rda_info <- attr(reducedDim(tse, "RDA"), "significance")

#Print out rda_info
rda_info$permanova |>
  knitr::kable()

#To ensure that the homogeneity assumption holds, we retrieve the corresponding 
#information from the results of RDA. In this case, none of the p-values is lower 
#than the significance threshold, and thus homogeneity is observed.
rda_info$homogeneity |>
  knitr::kable()

#visualize the weight and significance of each variable on the similarity
#between samples with an RDA plot
# Load packages for plotting function


# Generate RDA plot coloured by clinical status
plotRDA(tse, "RDA", colour_by = "Season")
#plotRDA(tse, "RDA")
#From the plot above, we can see that only age significantly describes differences 
#between the microbial profiles of different samples 
#Statistically significant (P < 0.05)


################################################################################
library(stringr)

## ???? Isolate metabolic pathways 

##
# Drop off those bacteria that do not include information in Phylum or lower levels
mae[[1]] <- mae[[1]][!is.na(rowData(mae[[1]])$Phylum), ]
# Clean taxonomy data, so that names do not include additional characters
rowData(mae[[1]]) <- DataFrame(apply(rowData(mae[[1]]), 2, 
                                     str_remove, pattern = "._[0-9]__"))


# Agglomerate microbiome data at family level
mae[[1]] <- mergeFeaturesByPrevalence(mae[[1]], rank = "Phylum")
# Does log10 transform for microbiome data
mae[[1]] <- transformAssay(mae[[1]], method = "log10", pseudocount = TRUE)

# Give unique names so that we do not have problems when we are creating a plot
rownames(mae[[1]]) <- getTaxonomyLabels(mae[[1]])

# Cross correlates data sets
correlations <- testExperimentCrossCorrelation(mae, 
                                               experiment1 = 5,
                                               experiment2 = 3,
                                               assay.type1 = "relabundance", 
                                               assay.type2 = "relabundance",
                                               method = "spearman", 
                                               p_adj_threshold = NULL,
                                               cor_threshold = NULL,
                                               # Remove when mia is fixed
                                               mode = "matrix",
                                               sort = TRUE,
                                               show_warnings = FALSE)
BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap) 

# Create a heatmap and store it
plot <- Heatmap(correlations$cor,
                # Print values to cells
                cell_fun = function(j, i, x, y, width, height, fill) {
                  # If the p-value is under threshold
                  if( !is.na(correlations$p_adj[i, j]) & correlations$p_adj[i, j] < 0.05 ){
                    # Print "X"
                    grid.text(sprintf("%s", "X"), x, y, gp = gpar(fontsize = 10, col = "#1dff00"))
                  }
                },
                heatmap_legend_param = list(title = "", legend_height = unit(5, "cm"))
)
plot

################################################################################


#add functional gene tse_species#add functional gene data to the TSE
# modify the Description entries
#colData(tse)$Description <- paste(colData(tse)$Description, "modified description")
# view modified variable
#head(tse$Description)
# simulate new data
#new_data <- runif(ncol(tse))
# store new data as new variable in colData
#colData(tse)$NewVariable <- new_data
# view new variable
#head(tse$NewVariable)





unique(tse$Year_Sample)
tse$Year_Sample%>% table()

unique(rowData(tse)$Phylum)
rowData(tse)$Phylum %>% table()


#For several custom analysis and visualization packages, 
#such as those from tidyverse, the SE data can be converted to a long data.frame
#format with meltAssay
molten_tse <- mia::meltAssay(tse,
                             add_row_data = TRUE,
                             add_col_data = TRUE,
                             assay.type = "relabundance")
#molten_tse
#dim(tse) 

# Not agglomerated data
tse_subset_by_feature <- tse[rowData(tse)$Phylum %in% c("Acidobacteria","Actinobacteria") & !is.na(rowData(tse)$Phylum), ]

# Show dimensions
#dim(tse_subset_by_feature)

# Agglomerate by phylum
tse_phylum <- tse %>% mergeFeaturesByRank(rank = "Phylum")

# Subset by feature and remove NAs
tse_phylum_subset_by_feature <- tse_phylum[rowData(tse_phylum)$Phylum %in% c("Acidobacteria","Actinobacteria") & !is.na(rowData(tse_phylum)$Phylum), ]

# Show dimensions
dim(tse_phylum_subset_by_feature)

# Subset by sample and feature and remove NAs
tse_subset_by_sample_feature <- tse[rowData(tse)$Phylum %in% c("Acidobacteria","Actinobacteria") & !is.na(rowData(tse)$Phylum), tse$Year_Sample %in% c("2015", "2016")]

# Show dimensions
dim(tse_subset_by_sample_feature)

#split by ranks
altExps(tse) <- splitByRanks(tse)
altExps(tse)

#split on other matters
tse_year = splitOn(tse, "Year_Sample")
tse_year

#Abundance: The relative abundance values for the top-5 taxonomic features can be 
#visualized as a density plot over a log scaled axis
plotAbundanceDensity(tse, layout = "density", assay.type = "relabundance",
                     n = 5, colour_by="Year_Sample", point_alpha=1/10) +
  scale_x_log10()

#Prevalence:
#Prevalence quantifies the frequency of samples where certain microbes were 
#detected (above a given detection threshold). The prevalence can be given as 
#sample size (N) or percentage (unit interval).
head(getPrevalence(tse, detection = 1/100, sort = TRUE, as_relative = TRUE))

head(getPrevalence(tse, detection = 1, sort = TRUE, assay.type = "counts",
                   as_relative = FALSE))

#Prevalence Analysis: 
head(getPrevalence(tse, rank = "Phylum", detection = 1/100, sort = TRUE,
                   assay.type = "counts", as_relative = TRUE))

#If you only need the names of the prevalent taxa, getPrevalentFeatures is 
#available. This returns the taxa that exceed the given prevalence 
#and detection thresholds.
getPrevalentFeatures(tse, detection = 0, prevalence = 80/100)
prev <- getPrevalentFeatures(tse, detection = 0, prevalence = 80/100,
                             rank = "Phylum", sort = TRUE)
prev

#Rare taxa
#more about that in the PDF
getRareFeatures(tse, detection = 0, prevalence = 80/100)
subsetByRareFeatures(tse, detection = 0, prevalence = 80/100)

#Plot prevalence 
tse_domain <- mergeFeaturesByRank(tse, "Domain", na.rm=TRUE)
altExp(tse, "Domain") <- tse_domain

rowData(altExp(tse,"Domain"))$prevalence <- 
  getPrevalence(altExp(tse,"Domain"), detection = 1/100, sort = FALSE,
                assay.type = "counts", as_relative = TRUE)




#Chapter: Taxonomic information 
#we will refer to co-abundant groups as CAGs, 
#which are clusters of taxa that co-vary across samples.
#checkTaxonomy(tse)
#Since the rowData can contain other data, 
#taxonomyRanks will return the columns mia assumes to contain the taxonomic information.
#taxonomyRanks(tse)
#rowData(tse)[, taxonomyRanks(tse)]
#all(!taxonomyRankEmpty(tse, rank = "Species"))
#table(taxonomyRankEmpty(tse, rank = "Species"))
#table(taxonomyRankEmpty(tse, rank = "Genus"))
#head(getTaxonomyLabels(tse,with_rank = TRUE))



#assay(altExp(tse, "Family"), "relabundance")[1:5, 1:7]

#Rare: 
#Rare taxa can also be aggregated into a single group “Other” instead of 
#filtering them out. A suitable function for this is mergeFeaturesByPrevalence.
altExp(tse, "Species_byPrevalence") <- mergeFeaturesByPrevalence(tse, 
                                                                 rank = "Species", 
                                                                 other_label = "Other", 
                                                                 prevalence = 5 / 100, 
                                                                 detection = 0, 
                                                                 as_relative = T)
altExp(tse, "Species_byPrevalence")

assay(altExp(tse, "Species_byPrevalence"), "relabundance")[88:92, 1:7]

#Taxa Clustering 
#library(bluster)
# The result of the CLR transform is stored in the assay clr
#tse <- transformAssay(tse, method = "clr", pseudocount = 1)

#tse <- transformAssay(tse, assay.type = "clr", method = "z", 
#                      MARGIN = "features")

# Cluster (with euclidean distance) on the features of the z assay
#tse <- cluster(tse,
#               assay.type = "z",
#               clust.col = "hclustEuclidean",
#               MARGIN = "features",
#               HclustParam(dist.fun = stats::dist, method = "ward.D2"))

# Declare the Kendall dissimilarity computation function
#kendall_dissimilarity <- function(x) {
#  as.dist(1 - cor(t(x), method = "kendall"))
#}

# Cluster (with Kendall dissimilarity) on the features of the z assay
#tse <- cluster(tse,
#               assay.type = "z",
#               clust.col = "hclustKendall",
#               MARGIN = "features",            
#               HclustParam(dist.fun = kendall_dissimilarity, method = "ward.D2"))

# Checking the clusters
clusters_euclidean <- rowData(tse)$hclustEuclidean
head(clusters_euclidean, 10)
clusters_kendall <- rowData(tse)$hclustKendall
head(clusters_kendall, 10)

#plotting clusters 
library(ggplot2)
library(patchwork) # TO arrange several plots as a grid
plot1 <- ggplot(as.data.frame(rowData(tse)), aes(x = clusters_euclidean)) +
  geom_bar() +
  labs(title = "CAG size distribution (Euclidean distance)",
       x = "Clusters", y = "Feature count (n)")
plot2 <- ggplot(as.data.frame(rowData(tse)), aes(x = clusters_kendall)) +
  geom_bar() +
  labs(title = "CAG size distribution (1 - tau)",
       x = "Clusters", y = "Feature count (n)")
TaxaClusters =  plot1 + plot2 + plot_layout(ncol = 2)
TaxaClusters
pdf("figures/TaxaClusters.pdf")
print(TaxaClusters)
dev.off() 

#Data Transformation 
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance", pseudocount = 1)
tse <- transformAssay(x = tse, assay.type = "relabundance", method = "clr", 
                      pseudocount = 1, name = "clr")
#head(assay(tse, "clr"))

#In ‘pa’ transformation, abundance table is converted to present/absent table.
tse <- transformAssay(tse, method = "pa")
#head(assay(tse, "pa"))

# list of abundance tables that assays slot contains
assays(tse)

tse <- mia::estimateRichness(tse, 
                             assay.type = "counts", 
                             index =   c("ace", "chao1", "hill", "observed"), 
                             name=  c("ace", "chao1", "hill", "observed"))

#head(tse$observed)

#Figure 7.1: Richness estimates plotted grouped by sampleid colored by year
library(scater)
plotColData(tse, 
            "observed", 
            "sampleid", 
            colour_by = "Year_Sample") +
  theme(axis.text.x = element_text(angle=45,hjust=1)) + 
  ylab(expression(Richness[Observed]))

tse <- mia::estimateDiversity(tse, 
                              assay.type = "counts",
                              index =  c("coverage", "fisher", "gini_simpson", "inverse_simpson",
                                         "log_modulo_skewness", "shannon"), 
                              name =  c("coverage", "fisher", "gini_simpson", "inverse_simpson",
                                        "log_modulo_skewness", "shannon"))
#head(tse$log_modulo_skewness)

library(ggsignif)
library(ggplot2)
library(patchwork)
library(ggsignif)

# Subsets the data. Takes only those samples that are from feces, skin, or tongue,
# and creates data frame from the collected data
df <- as.data.frame(colData(tse)[tse$Year_Sample)

# Changes old levels with new levels
df$Year_Sample <- factor(df$Year_Sample)

# For significance testing, all different combinations are determined
comb <- split(t(combn(levels(df$Year_Sample), 2)), 
              seq(nrow(t(combn(levels(df$Year_Sample), 2)))))

ggplot(df, aes(x = Year_Sample, y = shannon)) +
  # Outliers are removed, because otherwise each data point would be plotted twice; 
  # as an outlier of boxplot and as a point of dotplot.
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = comb, map_signif_level = FALSE) +
  theme(text = element_text(size = 10))
#Option to plot 

#Faith phylogenetic diversity
tse <- mia::estimateFaith(tse,
                          assay.type = "counts")
head(tse$faith)
#Evenness 
tse <- estimateEvenness(tse, 
                        assay.type = "counts", 
                        index=c("camargo", "pielou", "simpson_evenness", "evar", "bulla"),
                        name = c("camargo", "pielou", "simpson_evenness", "evar", "bulla"))
head(tse$camargo)

#Donminance
tse <- estimateDominance(tse, 
                         assay.type = "counts", 
                         index=c("absolute", "dbp", "core_abundance", "gini", "dmn", "relative",
                                 "simpson_lambda"), 
                         name = c("absolute", "dbp", "core_abundance", "gini", "dmn", "relative",
                                  "simpson_lambda"))

head(tse$relative)

#Rarity 
tse <- mia::estimateDiversity(tse, 
                              assay.type = "counts",
                              index = "log_modulo_skewness")

install.packages("remotes")
remotes::install_github("microbiome/microbiome")

# getRareFeatures returns the inverse
rare <- getRareFeatures(tse,
                        rank = "Phylum",
                        detection = 1/100,
                        prevalence = 50/100,
                        as_relative = TRUE)
head(rare)

#head(getRareTaxa(tse, prevalence = 50/100, include_lowest = FALSE))



# By default, reference is median of all samples. The name of column where results
# is "divergence" by default, but it can be specified. 
tse <- estimateDivergence(tse)

# The method that are used to calculate distance in divergence and 
# reference can be specified. Here, euclidean distance and dist function from 
# stats package are used. Reference is the first sample.
tse <- estimateDivergence(tse, name = "divergence_first_sample", 
                          reference = assays(tse)$counts[,1], 
                          FUN = stats::dist, method = "euclidean")

# Reference can also be median or mean of all samples. 
# By default, divergence is calculated by using median. Here, mean is used.
tse <- estimateDivergence(tse, name = "divergence_average", reference = "mean")

colData(tse)

#A plot comparing all the diversity measures calculated above and stored in 
#colData can then be constructed directly.
plots <- lapply(c("observed", "shannon", "simpson_evenness", "relative", "faith", "log_modulo_skewness"),
                plotColData,
                object = tse,
                x = "Year_Sample",
                colour_by = "Year_Sample")

plots <- lapply(plots, "+", 
                theme(axis.text.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.ticks.x = element_blank()))

((plots[[1]] | plots[[2]] | plots[[3]]) / 
    (plots[[4]] | plots[[5]] | plots[[6]])) +
  plot_layout(guides = "collect")

#Better to compare each group of measures on its own with correlation analysis 

#Beta Diversity 
#Comparing communities by beta diversity analysis
# Load package to plot reducedDim
library(scater)

# Run PCoA on relabundance assay with Bray-Curtis distances
tse <- runMDS(tse,
              FUN = vegan::vegdist,
              method = "bray",
              assay.type = "relabundance",
              name = "MDS_bray")

# Create ggplot object
p <- plotReducedDim(tse, "MDS_bray")

# Calculate explained variance
e <- attr(reducedDim(tse, "MDS_bray"), "eig")
rel_eig <- e / sum(e[e > 0])

# Add explained variance for each axis
p <- p + labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
              y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = ""))

#MDS plot based on the Bray-Curtis distances on the GlobalPattern dataset.
p

#Factors in ordination are: 	Assay type (relabundance/absolute)
# & the Beta diversity metric (Bray-Curtis, Aitchison,... Jaccard,...)
#the choice of these three factors can affect the resulting lower-dimensional data
# Run NMDS on relabundance assay with Bray-Curtis distances
tse <- runNMDS(tse,
               FUN = vegan::vegdist,
               method = "bray",
               assay.type = "relabundance",
               name = "NMDS_bray")

# Run MDS on clr assay with Aitchison distances
tse <- runMDS(tse,
              FUN = vegan::vegdist,
              method = "euclidean",
              assay.type = "clr",
              name = "MDS_aitchison")

# Run NMDS on clr assay with Euclidean distances
tse <- runNMDS(tse,
               FUN = vegan::vegdist,
               method = "euclidean",
               assay.type = "clr",
               name = "NMDS_aitchison")

# Load package for multi-panel plotting
library(patchwork)

# Generate plots for all 4 reducedDims
plots <- lapply(c("MDS_bray", "MDS_aitchison",
                  "NMDS_bray", "NMDS_aitchison"),
                plotReducedDim,
                object = tse)

# Generate multi-panel plot
#Comparison of MDS and NMDS plots based on the Bray-Curtis or 
#Aitchison distances on the GlobalPattern dataset.
wrap_plots(plots) +
  plot_layout(guides = "collect")

#PCA
tse <- runPCA(tse,
              name = "PCA",
              assay.type = "counts",
              ncomponents = 10)

plotReducedDim(tse, "PCA")

#UMAP for ordination
tse <- runUMAP(tse,
               name = "UMAP",
               assay.type = "counts",
               ncomponents = 3)

plotReducedDim(tse, "UMAP",
               ncomponents = c(1:3))

#Explained variance
# “stress” function measures the difference in pairwise similarities
# between the data points in the original
# Load vegan package
library(vegan)

# Quantify dissimilarities in the original feature space
x <- assay(tse, "relabundance") # Pick relabunance assay separately
d0 <- as.matrix(vegdist(t(x), "bray"))

# PCoA Ordination = MDS 
#Confusingly, PCoA is also abbreviated PCO, and is also known as 
#metric MultiDimensional Scaling (MDS)
pcoa <- as.data.frame(cmdscale(d0, k = 2))
names(pcoa) <- c("PCoA1", "PCoA2")

# Quantify dissimilarities in the ordination space
dp <- as.matrix(dist(pcoa))

# Calculate stress i.e. relative difference in the original and
# projected dissimilarities
stress <- sum((dp - d0)^2) / sum(d0^2)

#A Shepard plot visualizes the original versus the ordinated dissimilarity
# between the observations. 
ord <- order(as.vector(d0))
df <- data.frame(d0 = as.vector(d0)[ord],
                 dmds = as.vector(dp)[ord])

ggplot(df, aes(x = d0, y = dmds)) +
  geom_smooth() +
  geom_point() +    
  labs(title = "Shepard plot",
       x = "Original distance",
       y = "MDS distance",   
       subtitle = paste("Stress:", round(stress, 2))) +
  theme_bw()




#Here





#Example with metabolites 
library(mia)
data(HintikkaXOData, package="mia")
mae <- HintikkaXOData 

#Example with Clustering of data 
library(mia)
data("enterotype", package = "mia")
tse_cluster <- enterotype


#colData is for process data

#Multiexperiments is for metatranscriptomics 

#assay(tse, "counts")[1:5,1:7]
#tse <- relAbundanceCounts(tse)
#assays(tse)
#assay(tse, "relabundance")
#colData(tse)
#rowData(tse)

'''
transformAssay(
  method = c("alr", "chi.square", "clr", "frequency", "hellinger", "log", "log10",
             "log2", "max", "normalize", "pa", "range", "rank", "rclr", "relabundance", "rrank",
             "standardize", "total", "z"))

#assay(tse, "counts")[1:5,1:7]
#tse <- relAbundanceCounts(tse)
#assays(tse)
#assay(tse, "relabundance")
#colData(tse)
#rowData(tse)
'''
#top_features <- getTopFeatures(tse, method="median", top=10)
#rowData(tse)[top_features, taxonomyRanks(tse)]


'''



```{r Agglomerate data using the Alternative experiment transformation in the Summarized experiment data container}

# Agglomerate the data to Phylym level
tse_phylum <- agglomerateByRank(tse, "Phylum")
# both have the same number of columns (samples)
dim(tse)
#Decreasing the dimensionality of the assay matrix due to agglomerating by phylum rank
dim(tse_phylum)

#Question, why are not they equal in dimensionality between both? 107 vs. 103 (grouping is the reason)
length(unique(row_data_tax$phylum))

#assay(se_phylum)

# Add the newtable as an alternative experiment
altExp(tse, "phylum") <- tse_phylum
altExpNames(tse)

# Pick a sample subset: this acts on both altExp and assay data
tse[,1:10]
dim(altExp(tse[,1:10],"phylum"))
```




```{r transcriptomic data }

dim(row_data_genes)
colnames(row_data_genes)

#Multiple experiments relate to complementary measurement types, such as transcriptomic or metabolomic profiling of #the microbiome or the host. 
#if the samples can be matched directly 1:1, then transcriptomic data can be stores as Separate altExp
# Here, must learn how to use another formating becuase this is not the case use MultiAssay anot altExp 
#Go back and read more about MultiAssay 

#To construct a MultiAssayExperiment object, just combine multiple TreeSE data containers. Here we import metabolite #data from the same study.

```

```{r transcript}

assays = SimpleList(counts = assay_data)
colData = DataFrame(samdat)
rowData = DataFrame(row_data_genes)

tse_transcript<- TreeSummarizedExperiment(assays = assays,
                                          colData = colData,
                                          rowData = rowData
)

tse_transcript

```


```{r mae - MultiAssayExperiment}
# Create an ExperimentList that includes experiments
experiments <- ExperimentList(microbiome = tse, 
                              transcript = tse_transcript)

# Create a MAE
mae <- MultiAssayExperiment(experiments = experiments)

mae

```


```{r melt SE into long dataframe}

#the SE data can be converted to long data.frame format with meltAssay
tse <- transformSamples(tse, method="relabundance")

molten_tse <- meltAssay(tse,
                             add_row_data = TRUE,
                             add_col_data = TRUE,
                             abund_values  = "relabundance")
molten_tse



#help("meltAssay")


```

```{r subsetting se by column }

#Subsetting data helps to draw the focus of analysis on particular sets of samples and / or features. When dealing with large data sets, the subset of interest can be extracted and investigated separately. This might improve performance and reduce the computational load.

dim(tse)

# Note: when subsetting by sample, expect the number of columns to decrease; when subsetting by feature, expect the number of rows to decrease.

head(colData)

# inspect possible values for domain of life of the micro-organism 
unique(tse$Year_Sample)

# show recurrence for each value
tse$Year_Sample %>% table()

# subset by by year 2015 and 2016 (28 smaples)
tse_subset_by_year <- tse[ , tse$Year_Sample %in% c("2015", "2016")]

# show dimensions
dim(tse_subset_by_year)

```

``` {r  Subset by feature (row-wise)}

# inspect possible values for phylum
unique(rowData(tse)$phylum)

# show recurrence for each value
rowData(tse)$phylum %>% table()

#Non-agglomerated data

# subset by feature
tse_subset_by_feature <- tse[rowData(tse)$phylum %in% c("Actinobacteria","Proteobacteria") & !is.na(rowData(tse)$phylum), ]

# show dimensions
dim(tse_subset_by_feature)

#Agglomerated data 

# agglomerate by Phylum
tse_phylum <- tse %>% agglomerateByRank(rank = "Phylum")

# subset by feature and get rid of NAs
tse_phylum_subset_by_feature <- tse_phylum[rowData(tse_phylum)$phylum %in% c("Actinobacteria","Proteobacteria") & !is.na(rowData(tse_phylum)$phylum), ]

# show dimensions
dim(tse_phylum_subset_by_feature)


#Alternatively, the code below returns the not agglomerated version of the data.
# store features of interest into phyla
phyla <- c("Phylum:Actinobacteria", "Phylum:Proteobacteria")
phyla1 <- c("Actinobacteria", "Proteobacteria")
# subset by feature
tse_phylum_subset_by_feature <- tse_phylum[rowData(tse_phylum)$phylum %in% phyla1 & !is.na(rowData(tse_phylum)$phylum), ]
# show dimensions
dim(tse_subset_by_feature)
```

``` {r Subset by sample and feature}

# subset by year and feature and get rid of NAs
tse_subset_by_sample_feature <- tse[rowData(tse)$phylum %in% c("Actinobacteria","Proteobacteria") & !is.na(rowData(tse)$phylum), tse$Year_Sample %in% c("2015", "2016")]

# show dimensions
dim(tse_subset_by_sample_feature)


#Resulting dimensionality is correct 
```
``` { r splitting the data } 

altExps(tse) <- splitByRanks(tse)
altExps(tse)



```



```{r miaTime install} 



#install.packages("devtools", dependencies = TRUE)
#library(devtools)
#devtools::install_github("microbiome/miaTime")
#library(miaTime)

```

```{r A Jitter plot based on relative abundance data} 



# Add relative abundances
tse <- transformSamples(tse, method = "relabundance")

plotAbundance(tse[rowData(tse)$Phylum %in% Phylum],
              rank = "Phylum",
              abund_values = "relabundance")

```


```{r prevalence }

head(getPrevalence(tse, detection = 1/100, sort = TRUE, as_relative = TRUE))

head(getPrevalence(tse, detection = 1, sort = FALSE, assay_name = "counts",
                   as_relative = FALSE))

# Agglomerate taxa abundances to Phylum level, and add the new table to the altExp slot
altExp(tse,"Phylum") <- agglomerateByRank(tse, "Phylum")
# Check prevalence for the Phylum abundance table from the altExp slot
getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = TRUE,
              assay_name = "counts", as_relative = TRUE)
#This is better
getPrevalence(tse, rank = "Phylum", detection = 1/100, sort = TRUE,
              assay_name = "counts", as_relative = TRUE)

prev <- getPrevalentTaxa(tse, detection = 0, prevalence = 50/100,
                         rank = "Phylum", sort = TRUE)
prev


getRareTaxa(tse)
```

```{r installing  scuttle}
#BiocManager::install("scuttle")
```

```{r plotting prevalence }

rowData(altExp(tse,"Phylum"))$prevalence <- 
  getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = FALSE,
                assay_name = "counts", as_relative = TRUE)

library(scater)
plotRowData(altExp(tse,"Phylum"), "prevalence", colour_by = "Phylum")
```

```{r install ggsignif}
devtools::install_github("const-ae/ggsignif")
```

```{r diversity/richness/evenness/dominance}
#Richness
se <- mia::estimateRichness(se, 
                            abund_values = "counts", 
                            index = "observed", 
                            name="observed")

colData(se)$observed


#Diversity 
se <- mia::estimateDiversity(se, 
                             abund_values = "counts",
                             index = "shannon", 
                             name = "shannon")
head(colData(se)$shannon)



#if( !require(ggsignif) ){
#  install.packages(ggsignif)
#}
#library(ggplot2)
#library(ggsignif)

# Subsets the data. Takes only those samples that are from feces, skin, or tongue,
# and creates data frame from the collected data
df <- as.data.frame(colData(se)[colData(se)$Year_Sample %in% 
                                  c("2014", "2015", "2016"), ])

#Take all samples 
df1 <- as.data.frame(colData(se))
# Changes old levels with new levels
df$Year_Sample <- factor(df$Year_Sample)


# For significance testing, all different combinations are determined
comb <- split(t(combn(levels(df$Year_Sample), 2)), 
              seq(nrow(t(combn(levels(df$Year_Sample), 2)))))


ggplot(df, aes(x = Year_Sample, y = shannon)) +
  # Outliers are removed, because otherwise each data point would be plotted twice; 
  # as an outlier of boxplot and as a point of dotplot.
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = comb, map_signif_level = FALSE) +
  theme(text = element_text(size = 10))
ggsave("figures/Diversity1.png")


```

```{r estimate diversity }
#Use this for pre-processing 

#The available indices include the ‘Coverage’, ‘Faith's phylogenetic diversity’, ‘Fisher alpha’, ‘Gini-Simpson’, ‘Inverse Simpson’, ‘log-modulo skewness’, and ‘Shannon’ diversity indices


#Richness
tse <- mia::estimateRichness(tse, 
                             assay_name = "counts", 
                             index = "observed", 
                             name="observed")

sample_description = as.data.frame(colData(tse)$observed)
setnames(sample_description, old = colnames(sample_description), new = "observed")

#Diversity 
tse <- mia::estimateDiversity(tse, 
                              assay_name = "counts",
                              index =  "shannon", 
                              name = "shannon")

sample_description$shannon <- paste(colData(tse)$shannon)


#Dominance 
tse <- estimateDominance(tse, 
                         assay_name = "counts", 
                         index="relative")

sample_description$relative <- paste(colData(tse)$relative)


#Rarity 
tse <- mia::estimateDiversity(tse, 
                              assay_name = "counts",
                              index ="log_modulo_skewness")

sample_description$log_modulo_skewness <- paste(colData(tse)$log_modulo_skewness)

#Evenness  
tse <- estimateEvenness(tse, 
                        abund_values = "counts", 
                        index="simpson")

sample_description$simpson <- paste(colData(tse)$simpson)

#Divergence
tse <- mia::estimateDivergence(tse,
                               assay_name = "counts",
                               reference = "median",
                               FUN = vegan::vegdist)

sample_description$divergence <-paste(colData(tse)$divergence)

#Fischer's alpha diversity 

tse <- mia::estimateDiversity(tse, 
                              assay_name = "counts",
                              index = "fisher")

sample_description$fisher <- paste(colData(tse)$fisher)


#Gini-Simpson diversity
tse <- mia::estimateDiversity(tse, 
                              assay_name = "counts",
                              index = "gini_simpson")

sample_description$gini_simpson <- paste(colData(tse)$gini_simpson)


#Inverse Simpson diversity
tse <- mia::estimateDiversity(tse, 
                              assay_name = "counts",
                              index = "inverse_simpson")

sample_description$inverse_simpson <- paste(colData(tse)$inverse_simpson)

install.packages("writexl")
library("writexl")
write_xlsx(sample_description,"C:/Users/MWP-WKS050217/Desktop/RUG_Asala/PhD/Code/microbiome_practice/output_data/Diversity_Measures.xlsx")

colnames(sample_description)
view(sample_description)
#write.csv(sample_description,"output_data\\AlphaDiversity.csv", row.names = FALSE)
#Sorting dataframe by diversity values 
sample_description[
  with(sample_description, order("shannon", "fisher","gini_simpson", "inverse_simpson")),
]

sample_description %>% 
  as_tibble %>% 
  mutate_if(is.character, as.numeric) %>% 
  
  ggpairs()

cormat = data.frame()
cormat <- round(cor(sample_description),2)
head(cormat)

```



``` {r plot relative abundance } 


tse <- transformCounts(tse, method = "relabundance")
plotAbundance(tse, rank = "phylum", abund_values = "relabundance")

prev_phylum <- getPrevalentTaxa(tse, rank = "phylum",
                                detection = 0.01)
getPrevalentTaxa
library(patchwork)
plots <- plotAbundance(tse[rowData(tse)$Phylum %in% prev_phylum],
                       features = "Year_Sample",
                       rank = "Phylum",
                       abund_values = "relabundance")
plots$abundance / plots$Year_Sample +
  plot_layout(heights = c(9, 1))

??plotAbundance
ggsave("figures/RelAbundance_PrevPhylum.png")

plots <- plotAbundance(tse[rowData(tse)$Phylum %in% prev_phylum],
                       features = "Year_Sample",
                       rank = "Phylum",
                       abund_values = "relabundance")
plots$abundance / plots$Year_Sample +
  plot_layout(heights = c(9, 1))

```
```{r plot prevalence taxa }
plotTaxaPrevalence(tse, rank = "Phylum",
                   detections = c=(0. 0.001, 0.1, 0.2)
                   
                   ggsave("figures/Prev_phylum.png")
                   ```
                   ```{r relative abundance vs prevalence}
                   
                   plotPrevalentAbundance(tse, rank = "Family",
                                          colour_by = "Phylum")+ 
                     scale_x_log10()
                   
                   ggsave("figures/RelAbundance_Prev.png")
                   ```
                   
                   ```{r }
                   library(scater)
                   
                   plots <- lapply(c("shannon", "faith"),
                                   plotColData,
                                   object = tse, colour_by = "Sample_Date")
                   plots[[1]] + plots[[2]] +
                     plot_layout(guides = "collect")
                   
                   ```
                   
                   
                   