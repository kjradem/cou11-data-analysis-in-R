# cou11-data-analysis-in-R
A data analysis pipeline for RNA-Seq data in R

<h2>Description</h2>
This pipeline was developed to analyze RNA-Seq data for two Lactobacillus plantarum strains (WCFS1 & NC8) grown on different media (glucose & ribose). Methods include data processing, determination of differentially expressed (DE) genes, gene set enrichment analysis (GSEA) and visualization of results.

<h2>Installation</h2>

_Versions as of 04/Jun/2018_

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
 * Visualizing differentially expressed genes in KEGG pathways

<h2>Functions</h2>


<h2>Credits</h2>

Sjors Bongers, Daan Gilissen, Martijn Landman, Koen Rademaker, Ronald van den Hurk
