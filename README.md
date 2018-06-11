# cou11-data-analysis-in-R
A data analysis pipeline for RNA-Seq data in R

<h2>Description</h2>
This pipeline was developed to analyze RNA-Seq data for two Lactobacillus plantarum strains (WCFS1 & NC8) grown on different media (glucose & ribose). Methods include data processing, determination of differentially expressed (DE) genes, gene set enrichment analysis (GSEA) and visualization of results.

<h2>Installation</h2>

_Versions as of 10/Jun/2018_

Software:

* R version 3.5.0 (Joy in Playing)
* RStudio version 1.1.453

Libraries (installation automatically performed by script):

* Bioconductor version 3.7
* DBI version 1.0.0
* gplots version 3.0.1
* edgeR version 3.20.9
* KEGGREST version 1.18.1
* KEGG.db version 3.2.3
* pathview version 1.20.0
* RColorBrewer version 1.1-2
* xlsx version 0.5.7

Required libraries are loaded automatically, missing libraries are installed automatically. During installation the following console output may appear:

```Update all/some/none? [a/s/n]:```

reply to this with:

```a```

<h2>Workflow</h2>

![Workflow](https://raw.github.com/kjradem/cou11-data-analysis-in-R/master/Data/Legend.png)

<h2>Usage</h2>

This R script can be used to analyse RNA-seq data by:
 * Determining differentially expressed genes
 * Creating plots to interpret the difference between the experiments
 * Creating a heatmap of the most differentially expressed genes
 * Determining pathways of the genes
 * Determining overrepresented pathways
 * Creating the input for a GSEA analysis
 
<h3>Running the script</h3>

All analysis steps can be found in- and performed with 'data-analysis.r', functions are stored backstage in 'data-analysis-functions.r', which is automatically loaded.


It is advised to skip ```DETERMINE PATHWAYS FOR ALL GENES``` (line 55-57) unless necessary, as this may take up to 30 minutes to be performed for all genes. Running ```DETERMINE PATHWAYS FOR DE GENES``` (line 59-61) still yields biologically relevant results.

<h2>Functions</h2>

Brief descriptions of functions can be found below, more detailed and technical commentary on the functions can be found in the source code.

|              | AnnotateGenes |
| ------------ |:---------------------------------|
| Description  | Add annotation to gene data      |
| Returns      | Data frame containing annotation |

|               | CreateGeneExpressionFile |
| ------------- |:-------------------------------------|
| Description   | Create gene expression file for GSEA |
| Returns       | .txt file with gene expression data  |

|             | CreateGeneSetsFile |
| ----------- |:-------------------------------|
| Description | Create gene sets file for GSEA |
| Returns     | .gmx file with gene sets       |

|             | CreateGroup |
| ----------- |:------------------------------|
| Description | Store experimental conditions |
| Returns     | Factor 'group'                |

|             | CreateModel |
| ----------- |:--------------------|
| Description | Create glmLRT model |
| Returns     | glmLRT model        |

|             | DataProcessing |
| ----------- |:--------------------------------------|
| Description | Perform various data processing tasks |
| Returns     | RNA-Seq data as DGEList               |

|             | DetermineDEGenes |
| ----------- |:----------------------------------------------------------------------|
| Description | Determine top differentially expressed genes for a given logFC filter |
| Returns     | Data frame                                                            |

|             | DeterminePathwayOverrep |
| ----------- |:-------------------------------------------------------|
| Description | Determine overrepresentation of genes in KEGG pathways |
| Returns     | topKEGG object                                         |

|             | GetPathwaysForGenes |
| ----------- |:----------------------------------|
| Description | Get KEGG pathways for given genes |
| Returns     | Data frame                        |

|             | PlotHeatMap |
| ----------- |:-------------------------------------------|
| Description | Plot logFC values of genes into a heat map |
| Returns     | heatmap.2 heat map                         |

|             | PlotSampleDistances |
| ----------- |:----------------------|
| Description | Plot sample distances |
| Returns     | plotMDS plot          |

|             | PlotKEGGpathway |
| ----------- |:-----------------------------------------------------|
| Description | Plot a given KEGG pathway with logFC values of genes |
| Returns     | Pathview pathway plot                                |

|             | WriteResults |
| ----------- |:----------------------------|
| Description | Write results to Excel file |
| Returns     | Excel file                  |

<h2>Credits</h2>

Sjors Bongers, Daan Gilissen, Martijn Landman, Koen Rademaker, Ronald van den Hurk
