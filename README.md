# Scripts and Raw data for Wedge-tailed Shearwater Plastic analysis
The following files contain raw data for morphometric measurements, blood analyte measurements, location of sampling, and sampling date. 
  This file also contains which weight categories birds were assigned to and the 
  sample ID from RNA sequencing: 
  WTSH_data_set.xlsx 
  This file contains morphometric data:
  WTSH datasheet-Morphometrics.xlsx   

The following files contains the scripts for analyses and plots used for analyzing morphometric and chemical data:
  Script_complete_figures.R
  Barplots_ggplot_script.R
The following files contain code for processing and analyzing RNA sequencing data
  These files contain the code to align (STAR) and calculate expression (RSEM)
  RSEM-NicoleM.sh
  STAR-transcript.sh
  The following files contain code to quanitfy quality of alignment through mapped
  reads, uniquely mapped reads and multiple mapped reads
  WTSH-mapping-reads.R
  WTSH-uniq-map.txt
  WTSH-multiple-map.txt
  
  Files begining with "DE-" such as, "DE-nicole-SEX-NOT-plastic.R.", contain the 
  code for differential gene expression analyses. 

  These are the files for the outlier genes from the differential gene expression 
  analyses testing relationships between weight categories and presence of plastic
  weight-3factor-outliers-0.1.txt #these are the outlier genes under a p-value of 0.1
  weight-3factor-outliers.txt # these are the outlier genes under a p-value of 0.5
  These are the normalized counts from the two top differentially expressed genes 
  genes in the 3 weight category and presence of plastic analyses.
  WTSH_boxplot_HSPH1_counts - Sheet1.csv
  WTSH_boxplot_gene_ankrd - Sheet1.csv
  This scripts contain the code for creating the boxplots for the normalized counts  
  used in Figure 3C 
  boxplot_code.R  
  The following files contain code for the Gene enrichment and ontology analyses
  Nicole-extract-genes-Enrichment.sh
  GO-terms-table-nicol.R
