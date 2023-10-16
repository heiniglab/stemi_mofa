# Load Seurat object and create pseudobulk for MOFA analysis

#############################################
# Prerequisites - Load Libraries


library(dplyr)

library(Seurat, quietly = TRUE, verbose = FALSE)
library(SeuratDisk, quietly = TRUE, verbose = FALSE)

library(muscat)

library(reshape2)

library(SummarizedExperiment)

library(stringr)

###############################################
# Preqrequisites Configurations & Parameters

result_path = '../results/current/Reproduction'

seurat_input_data_name = 'G4_Seurat_Input_Replication.h5seurat'

name_save = 'V_AZIMUTH_REPRODUCTION'


## Define columns in seurat object containing sample-id and cluster annotations

sample_column = 'sample_id' # to be sample-id

cluster_column =  'cluster_id' # to be cluster_id; cell-type annotation


##################################################

### Should quantile normalization be applied?

quantile_normalization_single_cell = TRUE

standardize = FALSE

set_zero_na = FALSE

quantile_norm_feat = TRUE

# Functions

### Function for quantile normalization

quantile_normalization = function(X){
  set.seed(42)
  ranks = apply(X, 2, rank, ties.method = 'min')  # determine ranks of each entry
  
  sorted = data.frame(apply(X, 2, sort)) # sort the entries
  means = apply(sorted, 1, mean) # calculate the means
  
  normalized_data = apply(ranks, 2 ,function(x){ means[x]}) # substitute the means into ranks matrix
}


### Gene wise quantile normalization

stdnorm <- function(x) {
  set.seed(42)
  r = rank(x[!is.na(x)], ties.method="average")
  x[!is.na(x)] = qnorm(r / (length(x[!is.na(x)]) + 1))
  return(x)
}



# Load Data 

## Load Seurat object

###### Load the generated seurat objects

source_text = paste( result_path, '/', seurat_input_data_name , sep = '')
print(source_text)
print(file.info(source_text)$mtime)
rna_seurat_data = LoadH5Seurat(source_text, assays = "RNA", quietly = TRUE )




## Should contain raw counts
## After QC and Pre-processing
## annotations:
#### 'sample_id' - identification of sample/ patient incl. timepoint 
### 'cluster_id' - cell-type annotation 

colnames(rna_seurat_data[[]])


rna_seurat_data

rna_seurat_data_subset = rna_seurat_data # rename 

# Data Processing (Pseudobulk)

obs = rna_seurat_data_subset@meta.data

obs$cell = rownames(obs)

### add a dummy group column
obs$group_id = '1' ## group-id neede in DE analysis, here not --> DUMMY Variable

### assign sample-id: TBD remove library part for script sharing
obs$sample_id   = obs[[sample_column]]

### add cell-type assignment
obs$cluster_id = as.character(obs[[cluster_column]])

rownames(obs) = obs$cell

nrow(obs)

head(obs,2)

### Adjust B-cell mapping/ aggregate Azimuth cell-types
obs$cluster_id = str_replace(obs$cluster_id, 'B_intermediate|B_memory|B_naive', 'B cell')

sort(unique(obs$cluster_id))



## Add to seurat dataset

## group-id  

rna_seurat_data_subset = AddMetaData(object = rna_seurat_data_subset, metadata = obs[,'group_id', drop = FALSE], col.name = 'group_id')

## cluster-id

rna_seurat_data_subset = AddMetaData(object = rna_seurat_data_subset, metadata = obs[,'cluster_id', drop = FALSE], col.name = 'cluster_id')

## sample-id

rna_seurat_data_subset = AddMetaData(object = rna_seurat_data_subset, metadata = obs[,'sample_id', drop = FALSE], col.name = 'sample_id')

head(rna_seurat_data_subset@meta.data,2)

## Convert to SCE 

rna_sce = as.SingleCellExperiment(rna_seurat_data_subset)

rna_sce  # rows = genes; columns = cells



### Check amount of cells per sample and cluster

colSums(table(rna_sce$cluster_id,rna_sce$cluster_id ))

cells_per_sample_cluster = t(table(rna_sce$cluster_id, rna_sce$sample_id))

cells_per_sample_cluster = data.frame(cells_per_sample_cluster)

colnames(cells_per_sample_cluster) = c('Sample', 'Cluster_Cell_Type', 'amount_cells')

head(cells_per_sample_cluster,2)

## Analyze and calculate gene expression percentages per cluster

gene_list = list()

gene_cell_expr = list()

clusters = unique(rna_sce$cluster_id)
#clusters = unique(rna_sce$cell_type_Scanorama)

clusters

 for(i in clusters){
#    print(i)
    
    # subset data on cluster
    rna_sce_subset = rna_sce[,rna_sce$cluster_id== i] # cluster
    
    amount_cells = dim(rna_sce_subset)[2]
    
    # Calculate percentage of cells expressing gene
    amount_cells_expressing_gene = rowSums(assay(rna_sce_subset) > 0 )
    perc_cells_expressing_gene = (amount_cells_expressing_gene/ amount_cells) * 100
    
    
    gene_cell_expr[[i]] = data.frame(perc_cells_expressing_gene = perc_cells_expressing_gene, total_amount_cells_expressing_gene = amount_cells_expressing_gene)
    
    }

 ### Resulting amount of genes per cluster

gene_cell_expr_data = data.frame()

 for(i in names(gene_cell_expr)){
    data = gene_cell_expr[[i]]
    data$gene = rownames(gene_cell_expr[[i]])
    data$cluster = i
    gene_cell_expr_data = rbind( gene_cell_expr_data, data)
    }

head(gene_cell_expr_data,2)

## Add cluster, group and sample columns for aggregation

#### Add cluster_id, sample_id and group_id columns
(rna_sce <- prepSCE(rna_sce, 
    kid = 'cluster_id', # subpopulation assignments
    gid = 'group_id',  # group IDs (ctrl/stim)   # sample_id; using dummy sample id which corresponds to cluster columns
    sid = 'sample_id',   # sample IDs (ctrl/stim.1234)
    drop = FALSE))  # drop all other colData columns

nk <- length(kids <- levels(rna_sce$cluster_id))
ns <- length(sids <- levels(rna_sce$sample_id))
names(kids) <- kids; names(sids) <- sids

nk # amount of cluster

ns # amount of samples

kids  # cluster ids

length(kids) # amount cluster-id

## Aggregate single cell to pseudo-bulk data

pb <- aggregateData(rna_sce,
    assay = "counts", fun = "mean",
    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation

pb

#sum(colSums(assay(pb)))

### Save aggregated data

#save(  pb , file = paste0(result_path, '/G0_aggregated_RNA_input_correlations_all.RDS'))

# Normalization 

## RNA-Single-Seq

### Cell / gene expression data filtering

cell_perc_cluster = gene_cell_expr_data

head(cell_perc_cluster,2)

nrow(cell_perc_cluster)

length(unique(cell_perc_cluster$cluster))

##### Condition for filtering genes
cell_perc_cluster =  cell_perc_cluster[((cell_perc_cluster$perc_cells > 50) & (cell_perc_cluster$total_amount_cells_expressing_gene > 1200)) | ((cell_perc_cluster$perc_cells > 40) & (cell_perc_cluster$total_amount_cells_expressing_gene > 3000)) ,] 


nrow(cell_perc_cluster)

### Pseudobulk 

pb

all_genes = rownames(pb)

head(all_genes)

length(all_genes)

### Pre-Process

#### Remove Clusters (TBD)

names(assays(pb))

assay(pb, 'platelet') = NULL

assay(pb, 'plasmablast') = NULL

assay(pb, 'pDC') = NULL

assay(pb, 'Nkbright') = NULL

assay(pb, 'NK_Proliferating') = NULL

assay(pb, 'ILC') = NULL

assay(pb, 'HSPC') = NULL

assay(pb, 'Eryth') = NULL

assay(pb, 'Doublet') = NULL

assay(pb, 'doublet') = NULL

assay(pb, 'dnT') = NULL

assay(pb, 'cDC1') = NULL

assay(pb, 'CD8_TCM') = NULL

assay(pb, 'CD8_Proliferating') = NULL

assay(pb, 'CD4_Proliferating') = NULL

assay(pb, 'ASDC') = NULL

length(names(assays(pb)))

names(assays(pb))

#### Prepare gene-cluster dataframe + normalize

final_data = data.frame(samples = colnames(pb))

rownames(final_data) = final_data$samples

final_data_vis = data.frame(samples = colnames(pb))

rownames(final_data_vis) = final_data_vis$samples

genes_subset = cell_perc_cluster

name_save


for(i in unique(genes_subset$cluster)){
    data = assay(pb, i)


    ##### Normalize counts per sample (library size) - currently only for non-scanorama functions

    if(is.na(str_extract(name_save, 'scano')) == TRUE){
        scaling_factor = colSums(data) /mean(colSums(data))

        for (j in 1:ncol(data)){
            if(scaling_factor[j] != 0){
                data[,j] = data[,j]/ scaling_factor[j]
                }
            }
        }

    ### Subset data on genes with minimum expression in cluster
    data = data[rownames(data) %in% genes_subset$gene[genes_subset$cluster == i],]

    ### Alternative - cluster independent subsetting
    #data = data[rownames(data) %in% genes_subset,]

    ##### TBD pre-processing stepd

    if(is.na(str_extract(name_save, 'scano')) == TRUE){
        data = log2(data+1) # logarithmize count values (optional!)
        }

    #### Quantile normalization (TBD maybe also on complete dataset?)

    if(quantile_normalization_single_cell == TRUE){
        data_rows = rownames(data)
        data  = quantile_normalization(data ) 
        rownames(data) = data_rows
        }

    rownames(data) = paste0(i, '__' ,rownames(data))

    data = data.frame(t(data))

    expr_mean = data.frame( mean_expr = rowMeans(data))
    colnames(expr_mean) = i
    rownames(expr_mean) = rownames(data)

    final_data = merge(final_data, data, by = 0)
    final_data_vis = merge(final_data_vis, expr_mean, by = 0)

    rownames(final_data) =  final_data$Row.names
    rownames(final_data_vis) = final_data_vis$Row.names
    final_data$Row.names = NULL
    final_data_vis$Row.names = NULL
    }

   

head(final_data,2)

ncol(final_data)

nrow(genes_subset)

dim(final_data)

final_data$samples = str_replace(final_data$samples, '-.*', '')

head(final_data,2)



#### Filter genes

### Remove mitochondrial & ribosomal genes

head(final_data,2)

ncol(final_data)

final_data = final_data[, !colnames(final_data) %in% (colnames(final_data)[!is.na(str_extract(colnames(final_data), '__MT.*|__RPL.*|__RPS.*'))])]

ncol(final_data)   # minus sample + sample_id column --> 11.831

head(final_data,2)

## Genes with high variance

head(final_data,2)

final_data$samples = NULL

final_data$sample_id = NULL

ncol(final_data)

final_data$samples = rownames(final_data)

head(final_data,2)



#### Prepare long format

final_data_long = melt(final_data)

### Decide what to do with duplicates

head(final_data_long,2)

final_data_long$type = 'single_cell'

final_data_long$samples = str_replace(final_data_long$samples, '-.*', '')

final_data_long = final_data_long %>% group_by(samples, type, variable) %>% summarise(value = mean(value))  # take average in case same samples measured multiple times

length(unique(final_data_long$variable))

final_data_long$sample_id = final_data_long$samples
final_data_long$samples = NULL

# Quantile Normalization

head(final_data_long,2)

data_long = final_data_long

head(data_long,2)

nrow(data_long)

length(unique(data_long$sample_id))

### Normalization & wide format

### Standardize values

standardize

if(standardize == TRUE){
    data_long = merge(data_long, data_long %>% group_by(variable) %>% summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE)))
    
    data_long[data_long == 0] = NA
    
    data_long = data_long[(data_long$sd != 0) & (!is.na(data_long$sd)),]
    
    data_nas = data.frame(is.na(data_long))
    data_long$value = (data_long$value - data_long$mean)/data_long$sd
    
    #data_long = data.frame(data_long)
    data_long$mean = NULL
    data_long$sd = NULL
    data_long$value[data_nas$value] = NA
    }

unique(data_long$type)

## Prepare wide format for correlations

data_long$ident = paste0(data_long$type, '_0_', data_long$variable)

nrow(unique(data_long[,c('sample_id', 'ident')]))

nrow(data_long)

### Transform to wide

final_data = dcast(data_long, sample_id ~ ident , value.var = "value") # ! with this merging there might be NA values for some samples on some data types

head(final_data,2)

rownames(final_data) = final_data$sample_id

final_data$sample_id = NULL

ncol(final_data)

nrow(final_data)

### Deal with NA - Set NA for 0 observation + remove samples with only NA

head(final_data,2)

set_zero_na

if(set_zero_na == TRUE){
    final_data[final_data == 0] = NA
    }



### Remember NA's

data_nas = is.na(final_data)

ncol(final_data)

keep_samples = names(rowSums(data_nas))[rowSums(data_nas) != ncol(final_data)]

final_data = final_data[keep_samples,]

data_nas = data_nas[keep_samples,]

### Apply feature wise quantile normalization

quantile_norm_feat

if(quantile_norm_feat == TRUE){
    final_data = apply(final_data, 2,stdnorm)
    final_data = data.frame(final_data)
    final_data[data_nas] = NA
    final_data$sample_id = rownames(final_data)
    data_long = melt(final_data)
    data_long$type = str_extract(data_long$variable, '.*_0_')
    data_long$type  = str_replace(data_long$type , '_0_', '')
    data_long$variable = str_replace(data_long$variable, '.*_0_', '')
    }

# Save Prepared Data

name_save

write.csv(data_long, paste0(result_path, '/Combined_Data_', name_save, '.csv'))

length(unique(data_long$variable))


