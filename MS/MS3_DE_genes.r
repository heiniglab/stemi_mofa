### Function to find DE genes between clusters/ cell-types for a selected comparison

de_genes_samples = function(preferred_grouping, de_comparisons, rna_seurat_data, formula = '~ group_id', adj_counts = FALSE){
    
    result_data = data.frame()
    error_data = data.frame()
    
    ## Iterate over all comparisons
    for (i in 1:nrow(de_comparisons)){
        
        ## Get relevant comparison data
        comparison = de_comparisons[i,]
        comparison_column = comparison$comparison_column
        comparison_values  = c(comparison$group_value, comparison$ref_value)
        print(comparison_values)
        
        ## Subset data for comparison
        rna_seurat_data_subset = rna_seurat_data[,rna_seurat_data[[]][,comparison_column] %in% comparison_values]
        
        print('')
        print(rna_seurat_data_subset) ## print subsettet object
        
        ## Define columns for aggregation
        sample_column = 'display_name' # to be sample-id   
        cluster_column =  preferred_grouping # to be cluster_id
        
        ## Generate new columns to add
        
        obs = rna_seurat_data_subset@meta.data
        
        obs$cell = rownames(obs)
        
        obs$library = as.character(obs$library)
        obs$group_id = as.character(obs[[comparison_column]])
        obs$sample_id   = paste0(as.character(obs[[sample_column]] ), '-', obs$library)
        obs$cluster_id = as.character(obs[[cluster_column]])
        
        obs = merge(obs, obs %>% group_by(library) %>% summarise(nCount_Lib = sum(nCount_RNA,na.rm = TRUE)))
        rownames(obs) = obs$cell
        
        ## Add new columns to seurat data
        
        rna_seurat_data_subset = AddMetaData(object = rna_seurat_data_subset, metadata = obs[,'library', drop = FALSE], col.name = 'library_char')
        rna_seurat_data_subset = AddMetaData(object = rna_seurat_data_subset, metadata = obs[,'group_id', drop = FALSE], col.name = 'group_id')
        rna_seurat_data_subset = AddMetaData(object = rna_seurat_data_subset, metadata = obs[,'cluster_id', drop = FALSE], col.name = 'cluster_id')
        rna_seurat_data_subset = AddMetaData(object = rna_seurat_data_subset, metadata = obs[,'sample_id', drop = FALSE], col.name = 'sample_id')
        rna_seurat_data_subset = AddMetaData(object = rna_seurat_data_subset, metadata = obs[,'nCount_Lib', drop = FALSE], col.name = 'nCount_Lib')
                                                        
        ## Convert format to sce
        
        rna_sce <- as.SingleCellExperiment(rna_seurat_data_subset)
        
        ## Optional: adjustment of count values
        
        if(adj_counts == TRUE){
            lib_count_constant = mean(unique(rna_seurat_data_subset$nCount_Lib))
            amount_genes = dim(rna_seurat_data_subset)[1]
            
            counts(rna_sce ) = counts(rna_sce)/ matrix(rep(rna_sce$nCount_Lib, amount_genes), nrow = amount_genes, byrow = TRUE) * lib_count_constant
            }
    
        ## Add cluster_id, sample_id and group_id columns
        rna_sce <- prepSCE(rna_sce, 
                           kid = 'cluster_id', # subpopulation assignments
                           gid = 'group_id',  # group IDs (ctrl/stim)
                           sid = 'sample_id',   # sample IDs (ctrl/stim.1234)
                           drop = FALSE)  # drop all other colData columns
        
        nk <- length(kids <- levels(rna_sce$cluster_id))
        ns <- length(sids <- levels(rna_sce$sample_id))
        names(kids) <- kids; names(sids) <- sids
        
        print('')
        print(paste0('Number of samples: ', ns))
        #print(unique(rna_sce$sample_id))
        
        ## Aggregate to Pseudobulk Data
        
        pb <- aggregateData(rna_sce,
                            assay = "counts", fun = "sum",
                            by = c("cluster_id", "sample_id")) # one sheet per subpopulation
        
        ## Define model to estimate
        
        # Data preparation: very specific might have to be adapted for further use-cases
        
        pb$sample_id2 = colnames(pb)
        #print(colData(pb))
        
        model_matrix_basis = unique(colData(pb)[,c('sample_id2', 'group_id', 'library')]) #TBD could be made only on sample basis -- joinable matrix needs to be handed over to function
        
        obs_matrix = as.data.frame(colData(rna_sce))
        
        counts_per_lib = obs_matrix %>% group_by(library) %>%  summarise(amount_counts = sum(nCount_RNA, na.rm = TRUE)) ## TBD maybe counts per library need to be calculated differently
        
        ## Apply some normalization steps to counts_per_lib to adjust scale
        total_counts = sum(counts_per_lib$amount_counts)
        counts_per_lib$amount_counts_norm = counts_per_lib$amount_counts/total_counts
        
        amount_genes = dim(pb)[1]
        counts_per_lib$amount_counts_gene_basis = counts_per_lib$amount_counts/amount_genes
        
        
        model_matrix_basis = merge(model_matrix_basis, counts_per_lib)
        
        rownames(model_matrix_basis) = model_matrix_basis$sample_id2
        
        ### Adjust levels for library column
        model_matrix_basis$library = as.character(model_matrix_basis$library)
        model_matrix_basis$library = as.factor(model_matrix_basis$library)
        
        ### Adjust levels for group_id column
        group_value = comparison$group_value
        ref_value = comparison$ref_value
        
        model_matrix_basis$group_id = as.character(model_matrix_basis$group_id)
        model_matrix_basis$group_id = factor(model_matrix_basis$group_id, levels = c(ref_value, group_value))
        
        ###  Generate model matrix
        
        formula = as.formula(formula)
        # alternatives: '~  amount_counts + group_id' // '~  library + group_id' // '~ 0 + library +  group_id' // '~  amount_counts_norm + group_id'
        mm <- model.matrix(formula, model_matrix_basis)
        
        print('')
        print(paste0('Dim of model matrix- cols', dim(mm)[2])) # # tp be printed out for checking
        print(paste0('Rank of model matrix', rankMatrix(mm))) # to be printed out for checking
        print(head(mm))  # to be printed out for checking
        
        ### Integrate checks on whether model can be estimated
        
        rank = rankMatrix(mm)
        dim_col = dim(mm)[2]
        
        ## Run DS based on model
        
        if(rank >= dim_col){
            
            result = pbDS(pb, verbose = TRUE, method = 'DESeq2' , design = mm, min_cells = 10)  # 
        
            ## Prepare output to return

            result_formatted = resDS(rna_sce, result, bind = "row")
            final_result = result_formatted[,c('gene', 'baseMean', 'contrast', 'lfcSE', 'stat', 'p_val', 'p_adj.loc',  'p_adj.glb')]
            final_result$cluster = result_formatted$cluster_id
            final_result$group = group_value
            final_result$reference = ref_value
            final_result$logfoldchanges = result_formatted$logFC
            final_result$pvals = result_formatted$p_val
            final_result$comparison_nr = comparison$nr_compariosn
            print(as.character(formula))
            final_result$formula = as.character(deparse(formula))

            ## Add information about sample sizes
            comparison_summary = as.data.frame(model_matrix_basis)
            comparison_summary = comparison_summary  %>% group_by(group_id) %>% count()

            comparison_summary$group = comparison_summary$group_id
            comparison_summary$group_n = comparison_summary$n

            final_result = merge(final_result, comparison_summary[,c('group_n', 'group')])

            comparison_summary$reference = comparison_summary$group_id
            comparison_summary$reference_n = comparison_summary$n
            final_result = merge(final_result, comparison_summary[,c('reference_n', 'reference')])

            print('')
            print(head(final_result))

            print('Design')
            head(result$args$design)
            print('Contrast')
            head(result$args$contrast)

            result_data = rbind(result_data, final_result)
        }
        
        else{
            error_data = rbind(error_data, comparison)
        }
        
        }
    result_return = list(result_data, error_data)
    return(result_return)
    ## TBD save all comparisons in one dataset
    }