# INSTALL/LOAD REQUIRED LIBRARIES
source('http://bioconductor.org/biocLite.R')
if(is.element('biocLite', installed.packages()[,1])==FALSE){
  biocLite()
}
required_libraries <- c('DBI', 'edgeR', 'gplots', 'KEGGREST', 'KEGG.db', 'pathview', 'RColorBrewer', 'xlsx')
for(library in required_libraries){
  if(is.element(library, installed.packages()[,1])==FALSE){
    biocLite(pkgs=library)
  }
}
library('DBI')
library('edgeR')
library('gplots')
library('KEGGREST')
library('KEGG.db')
library('pathview')
library('RColorBrewer')
library('xlsx')

# LOAD FUNCTIONS
source('data-analysis-functions.r')

# LOAD DATA
Counts <- read.delim(file.choose(), header=TRUE, skip=1, row.names=1)
Annotation <- read.delim(file.choose(), header=TRUE, skip=1, row.names=1)

# SET GLOBAL VARIABLES
CPM = 10
PCH_1 = 21
PCH_2 = 23
LOGFC_FILTER = 2
KEGG_SPECIES = 'lpl'

# VISUALIZE AND CHECK SEPARATION OF SAMPLES ON GROWTH MEDIUM AND STRAIN
All_group <- CreateGroup(c('WCFS1.glc', 'WCFS1.rib', 'NC8.glc', 'NC8.rib'))
All_data <- DataProcessing(All_group, 1, 8, CPM)
PlotSampleDistances('Distances between RNA-Seq samples', All_data, All_group)

# PROCESS DATA FOR WCFS1
WCFS1_group <- CreateGroup(c('WCFS1.glc', 'WCFS1.rib'))
WCFS1_data <- DataProcessing(WCFS1_group, 1, 4, CPM)
WCFS1_fit <- CreateModel('WCFS1', WCFS1_data, WCFS1_group)

# PROCESS DATA FOR NC8
NC8_group <- CreateGroup(c('NC8.glc', 'NC8.rib'))
NC8_data <- DataProcessing(NC8_group, 5, 8, CPM)
NC8_fit <- CreateModel('NC8', NC8_data, NC8_group)

# DETERMINE DE GENES
WCFS1_de_genes <- DetermineDEGenes(WCFS1_fit, LOGFC_FILTER)
NC8_de_genes <- DetermineDEGenes(NC8_fit, LOGFC_FILTER)

# DETERMINE PATHWAYS FOR ALL GENES
WCFS1_pathways_genes <- GetPathwaysForGenes(WCFS1_fit$table)
NC8_pathways_genes <- GetPathwaysForGenes(NC8_fit$table)

# DETERMINE PATHWAYS FOR DE GENES
WCFS1_pathways_de_genes <- GetPathwaysForGenes(WCFS1_de_genes)
NC8_pathways_de_genes <- GetPathwaysForGenes(NC8_de_genes)

# VALIDATE RESULTS FOR WCFS1
WCFS1_overrep_pathways <- DeterminePathwayOverrep(WCFS1_fit, Inf)
WCFS1_annotated_results <- AnnotateGenes(WCFS1_fit$table)

# VALIDATE RESULTS FOR NC8
NC8_overrep_pathways <- DeterminePathwayOverrep(NC8_fit, Inf)
NC8_annotated_results <- AnnotateGenes(NC8_fit$table)

# PLOT DATA
PlotSampleDistances('Distances between WCFS1 RNA-Seq samples', WCFS1_data, WCFS1_group)
PlotSampleDistances('Distances between NC8 RNA-Seq samples', NC8_data, NC8_group)
PlotHeatMap(WCFS1_de_genes, NC8_de_genes, 50)
PlotKEGGpathway(WCFS1_pathways_de_genes, 'lpl00620', 'Pyruvate metabolism', WCFS1_fit)
PlotKEGGpathway(NC8_pathways_de_genes, 'lpl00620', 'Pyruvate metabolism', NC8_fit)

# EXPORT RESULTS
WriteResults('Results/RNA_Seq_analysis_results.xlsx', WCFS1_annotated_results, 'WCFS1 DE genes', WCFS1_overrep_pathways, 'WCFS1 Overrep pathways', WCFS1_pathways_de_genes, 'WCFS1 DE genes pathways')
WriteResults('Results/RNA_Seq_analysis_results.xlsx', NC8_annotated_results, 'NC8 DE genes', NC8_overrep_pathways, 'NC8 Overrep pathways', NC8_pathways_de_genes, 'NC8 DE genes pathways')

# CREATE GSEA INPUT FILES
CreateGeneSetsFile()
CreateGeneExpressionFile(WCFS1_de_genes['logFC'], NC8_de_genes['logFC'], c('WCFS1', 'NC8'))
