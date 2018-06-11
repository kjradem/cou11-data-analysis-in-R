# PROVIDE ANNOTATION FOR GENES
AnnotateGenes <- function(genes){
  # Annotate genes.
  annotated_data <- cbind(genes, Annotation[rownames(genes),])
  return(annotated_data)
}

# CREATE GENE EXPRESSION FILE (.TXT) FOR GSEA
CreateGeneExpressionFile <- function(df_1, df_2, names){
  # Combine expression data of shared genes for two conditions, store in data frame and write to .txt file.
  shared_genes <- intersect(row.names(df_1), row.names(df_2))
  df_1 <- as.data.frame(df_1[row.names(df_1) %in% shared_genes,])
  df_2 <- as.data.frame(df_2[row.names(df_2) %in% shared_genes,])
  df <- cbind(shared_genes, NA, df_1, df_2, row.names=shared_genes)
  colnames(df) <- c('NAME', 'DESCRIPTION', names)
  write.table(df, file='Results/gsea_gene_expression.txt', row.names=F, quote=F, sep='\t')
}

# CREATE GENE SETS FILE (.GMX) FOR GSEA
CreateGeneSetsFile <- function(){
  max_nrow = 0
  # Load unique subclasses from Annotation, filter out duplicate/null values.
  subclasses <- as.character(unique(Annotation$subclass))
  subclasses <- subclasses[!(duplicated(tolower(subclasses)))]
  subclasses <- subclasses[subclasses != '']
  # Get corresponding genes per subclass, store in data frame and write to .gmx file.
  df <- data.frame(All=row.names(Annotation))
  for(index in 1:length(subclasses)){
    temp_subset <- subset(Annotation, Annotation$subclass==subclasses[index])['subclass']
    genes <- c(NA, row.names(temp_subset))
    df$col <- c(genes, rep('', nrow(df)-length(genes)))
    names(df)[index+1] = subclasses[index]
    if(nrow(temp_subset) > max_nrow){
      max_nrow = nrow(temp_subset)
    }
  }
  df <- subset(df, select = -c(All))
  df <- df[-c(max_nrow+2:nrow(df)), ]
  write.table(df, file='Results/gsea_gene_sets.gmx', row.names=F, quote=F, sep='\t')
}

# STORE EXPERIMENTAL CONDITIONS IN FACTOR NAMED GROUP
CreateGroup <- function(conditions){
  exp <- rep(conditions, each=2)
  group <- factor(exp)
  return(group)
}

# CREATE glmLRT MODEL
CreateModel <- function(strain, data, group){
  # Group samples by condition into design matrix.
  design <- model.matrix(~0+group, data=data$samples)
  colnames(design) <- levels(data$samples$group)
  # Create model for top differentially expressed genes.
  fit <- glmFit(data, design)
  if(strain == 'WCFS1'){
    mc <- makeContrasts(exp.r=WCFS1.glc-WCFS1.rib, levels=design)
  } else if(strain == 'NC8'){
    mc <- makeContrasts(exp.r=NC8.glc-NC8.rib, levels=design)
  }
  fit <- glmLRT(fit, contrast=mc)
  return(fit)
}

# PERFORM VARIOUS DATA PROCESSING TASKS
DataProcessing <- function(group, start, stop, cpm_filter){
  # Create DGEList object for storage of RNA-Seq data.
  y <- DGEList(counts=Counts[,start:stop], group=group)
  # Filter out genes below cpm_filter.
  keep.genes <- rowSums(cpm(y)>cpm_filter) >= 2
  y <- y[keep.genes,]
  y$samples$lib.size <- colSums(y$counts)
  # Determine scale factors using Trimmed Mean of M-values (TMM).
  y <- calcNormFactors(y, method='TMM')
  # Group samples by condition into design matrix.
  design <- model.matrix(~0+group, data=y$samples)
  colnames(design) <- levels(y$samples$group)
  # Estimate dispersions.
  y <- estimateGLMCommonDisp(y,design)
  y <- estimateGLMTrendedDisp(y,design, method='power')
  y <- estimateGLMTagwiseDisp(y,design)
  return(y)
}

# DETERMINE TOP DIFFERENTIALLY EXPRESSED GENES FOR A GIVEN LOGFC FILTER
DetermineDEGenes <- function(fit, filter_logFC){
  df <- topTags(fit, n=nrow(fit))
  df <- df$table
  df <- subset(df, abs(logFC) >= filter_logFC)
  return(df)
}

# DETERMINE OVERREPRESENTATION OF GENES IN KEGG PATHWAYS
DeterminePathwayOverrep <- function(fit, n_results){
  kegg_pathways <- kegga(fit, species.KEGG=KEGG_SPECIES)
  top_OR_pathways <- topKEGG(kegg_pathways, number=n_results)
  return(top_OR_pathways)
}

# GET KEGG PATHWAYS FOR GIVEN GENES
GetPathwaysForGenes <- function(genes){
  # Set up data frame.
  cols <- rownames(genes)
  n_pathways <- dbGetQuery(KEGG_dbconn(), 'SELECT COUNT(*) FROM pathway2name')[1,1]
  df <- data.frame(matrix(ncol=length(cols), nrow=n_pathways))
  colnames(df) <- cols
  # Store pathways per gene in data frame.
  max_n_pathways = 0
  for(index in 1:length(cols)){
    gene <- cols[index]
    try(query <- keggGet(c(paste('lpl:', gene, sep=''))), silent=T)
    if(exists('query')){
      pathways <- query[[1]]$PATHWAY
      if(!is.null(pathways)){
        for(index_2 in 1:length(pathways)){
          df[index_2, index] = pathways[index_2]
        }
        if(length(pathways) > max_n_pathways){
          max_n_pathways = length(pathways)
        }
      }
    }
  }
  df <- df[-c(max_n_pathways+1:nrow(df)), ]
  return(df)
}

# PLOT LOGFC VALUES OF GENES INTO A HEAT MAP
PlotHeatMap <- function(WCFS1_df, NC8_df, n_genes, ...){
  # Cuts data frames down to selected number of genes.
  WCFS1_df <- WCFS1_df[1:n_genes,]
  NC8_df <- NC8_df[1:n_genes,]
  shared_genes <- intersect(row.names(WCFS1_df), row.names(NC8_df))
  # Reduces data frames to only contain shared genes.
  WCFS1_df <- as.data.frame(WCFS1_df[row.names(WCFS1_df) %in% shared_genes,])
  NC8_df <- as.data.frame(NC8_df[row.names(NC8_df) %in% shared_genes,])
  # Combines data frames into single data frame and provides annotation.
  annotation_df <- as.data.frame(Annotation[row.names(Annotation) %in% shared_genes,])
  annotation_df <- subset(annotation_df, select=c('name', 'subclass'))
  df <- cbind(WCFS1_df['logFC'], NC8_df['logFC'], row.names=shared_genes)
  annotation_df <- cbind(df, annotation_df[rownames(df),], row.names=shared_genes)
  # Visualizes heat map and sets column and row names.
  rownames(df) <- paste(annotation_df[,'name'], annotation_df[,'subclass'], sep=' - ')
  colnames(df) <- c('WCFS1', 'NC8')
  color_palette <- colorRampPalette(c('green','blue'))(n = 64)
  heatmap.2(as.matrix(df),
            adjCol=c(NA, 0.5),
            cexCol=1.5,
            col=color_palette,
            Colv=F,
            dendrogram='none',
            density.info='none',
            margins=c(3,22),
            notecol='black',
            srtCol=0,
            trace='none')
}

# PLOT SAMPLE DISTANCES
PlotSampleDistances <- function(title, data, group, ...){
  # Set up colors and symbols for plot.
  if(length(levels(group)) == 4){
    colors <- rep(c('red', 'red', 'blue', 'blue'), 2)
  } else if(length(levels(group)) == 2){
    colors <- rep(c('red', 'blue'), 2)
  }
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  pch <- rep(c(PCH_1, PCH_2), length(levels(group)))
  pch_legend <- rep(c(16, 18), length(levels(group))/2)
  # Visualize plot.
  plotMDS(data, bg=colors[group], cex=2, col=1, pch=pch[group], xlab='Dimension 1', ylab='Dimension 2')
  title(title, line=0.5)
  legend('topright', col=colors, inset=c(-0.25,0), legend=levels(group), ncol=1, pch=pch_legend, title='Samples')
}

# PLOT A GIVEN KEGG PATHWAY WITH LOGFC VALUES OF GENES
PlotKEGGpathway <- function(df, pathway_id, pathway_name, fit, ...){
  # Select logFC values from genes in a specific pathway.
  matches <- vector('numeric')
  for(index in 1:ncol(df)){
    temp_vector <- df[,index]
    if(pathway_name %in% temp_vector){
      matches <- append(matches, colnames(df)[index])
    }
  }
  input_df <- subset(fit$table[row.names(fit$table) %in% matches,], select=logFC)
  input <- input_df[,length(input_df)]
  names(input) <- row.names(input_df)
  # Visualizes logFC values inside KEGG pathway.
  pathview(gene.data=input, species=KEGG_SPECIES, pathway=pathway_id, gene.idtype='KEGG')
}

# WRITE RESULTS TO EXCEL FILE
WriteResults <- function(file_name, annotated_results, sheet_name_1, or_pathways, sheet_name_2, pathways_de_genes, sheet_name_3){
  # Write results to sheets in Excel file.
  write.xlsx(annotated_results, file=file_name, sheetName=sheet_name_1, col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)
  write.xlsx(or_pathways, file=file_name, sheetName=sheet_name_2, col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)
  write.xlsx(pathways_de_genes, file=file_name, sheetName=sheet_name_3, col.names=TRUE, row.names=FALSE, append=TRUE, showNA=FALSE)
}