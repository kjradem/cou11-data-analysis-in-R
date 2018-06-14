# cou11-data-analysis-in-R
A data analysis pipeline for RNA-Seq data in R.

## I - Description
This pipeline was developed to analyse RNA-Seq data for two Lactobacillus plantarum strains (WCFS1 & NC8) grown on different media (glucose & ribose). Functionalities include data processing, determining differentially expressed (DE) genes and visualization of data.

## II - Installation 
_Versions as of 10/Jun/2018_

Software:
* R version 3.5.0 (Joy in Playing)
* RStudio version 1.1.453

Libraries (installation automatically performed):
* Bioconductor version 3.7
* DBI version 1.0.0
* gplots version 3.0.1
* edgeR version 3.20.9
* KEGGREST version 1.18.1
* KEGG.db version 3.2.3
* pathview version 1.20.0
* RColorBrewer version 1.1-2
* xlsx version 0.5.7

Make sure you run RStudio with administrator rights while running the installation.

During installation the following console output may appear:

```Update all/some/none? [a/s/n/]```

reply to this with:

```a```

## III - Usage
This R-script can be used to analyze RNA-Seq data by:
* Performing data processing to make RNA-Seq data suitable for analysis
* Creating plots to interpret the differences between experiments
* Determining differentially expressed genes
* Creating a heat map of the most differentially expressed genes
* Determining pathways for genes
* Visualizing logFC values in pathways
* Determining overrepresented pathways
* Creating input files for Gene Set Enrichment Analysis (GSEA)

### IIIa - Data
The files in the 'Data' folder were used while developing this pipeline, however other files may be used as well. Make sure the files follow these formats:

#### Counts data:
The first row above column names is reserved for comments.

| ID     | Condition 1 | Condition ...  |
| :----: | :---------: | :------------: |
| gene A | {value}     | ...            |
| gene B | {value}     | ...            |
| ...    | ...         | ...            |

#### Annotation data
The first row above column names is reserved for comments, annotation attributes can be added or removed freely depending on the situation. Keep in mind that some functions require specific annotation data to function properly.

| ID     | Gene name   | ...  |
| :----: | :---------: | :--: |
| gene A | {value}     | ...  |
| gene B | {value}     | ...  |
| ...    | ...         | ...  |

### IIIb - Running the script

#### Default
To run the entire pipeline as it was developed and tested, simply run the script 'data-analysis.r', with functions stores backstage in 'data-analysis-functions.r'.

#### Custom
To run specific functions of the pipeline, make sure the following segments (indivated by all-caps comments before the code) have been run previously:
* INSTALL/LOAD REQUIRED LIBRARIES
* LOAD FUNCTIONS
* LOAD DATA
* SET GLOBAL VARIABLES

### IIIc - Other
The function _GetPathwaysForGenes_ in lines 55 and 56 has been commented out on purpose. It is a time-intensive function (30+ minutes) to run for all genes and should only be ran if necessary.

## IV - Workflow

![Workflow](https://raw.github.com/kjradem/cou11-data-analysis-in-R/master/Data/Legend.png)

## V - Functions

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

## VI - Credits

Sjors Bongers, Daan Gilissen, Martijn Landman, Koen Rademaker, Ronald van den Hurk
