# Multisite phosphorylation of intrinsically disordered region of DVL facilitates Wnt signaling

## Authors:
Miroslav Micka, Jitender Kumar, Petra Paclíková, Zuzana Hayek, Kateřina Hanáková, Cherine Bechara, Hana Plešingerová, Sara Bologna, Ondrej Šedo, Elise Del Nero, Kristína Gömöryová, Tomáš Gybeľ, Marek Kravec, Zbyněk Zdráhal, Konstantinos Tripsianes**, Vítězslav Bryja*

*Corresponding author: bryja@sci.muni.cz

**Corresponding author: kostas.tripsianes@ceitec.muni.cz

## Reproducing Charge calculations
First run script "charge_functions.R". Here we edit functions from idpr R package (McFadden and Yanowitz, 2022), to enable charge calculation of phosphorylated proteins.
Next run script "charge_calculation.R". We either read analyzed protein sequences directly in the script or with read them from a text file. X stands for phosphoserine and Z stands for phosphothreonine in an analysed sequence.

## Reproducing the interactome analysis

### Deposition of raw data to PRIDE

Raw proteomic data can be accessed using the PRIDE with identifier [PXD067002](Interactome analysis of DVL3 S/T cluster mutants using proximity labeling). 

### Setting up the Docker container for running KNIME
In case of the analysis of protein complexes, raw data were searched using the [MaxQuant](https://www.maxquant.org/) software (v. 1.6.0.16). Resulting output, the proteinGroups.txt file, was further processed using the software container environment [OmicsWorkflows](https://github.com/OmicsWorkflows) (v. 3.7.2a): the workflow is stored within this repository as 6016_publication_template.knwf and can be inspected using the KNIME software.

To fully reproduce the analyses, run KNIME inside the Docker container using the 4.7.7 version of Docker image:
(Note: you need to have [Docker](https://docs.docker.com/get-docker/) installed on your computer)

1) Clone [this](https://github.com/OmicsWorkflows/KNIME_docker_vnc) repository locally to your computer
2) Adjust the start_container script for the folder which will contain your KNIME workspace (e.g. workspace-folder)
3) Run the start_container script to create a docker container as follows:
`cfprot/knime:4.7.7`, `5901`, `workspace-folder`

For more detailed instructions, please follow the tutorial [here](https://github.com/OmicsWorkflows/KNIME_docker_vnc)

### KNIME workflow import and description

Download the KNIME workflow (6016_publication_template.knwf) and import it into KNIME (`File` -> `Import KNIME workflow`). 

Run the particular nodes by right-clicking the node and `Execute` or directly pressing `F7`. 
Briefly, the workflow is as following:

* Data input (proteinGroups.txt file) DOHLEDAT A OPRAVIT
* Contaminants filtering (cRAP, Reverse, Only identified by Site)
* log2 transformation of protein intensities
* LoessF normalization
* Imputation of missing values by the global minimum
* Statistical testing using the limma test
* Exporting the results as .csv file 

The resulting output (ProteinGroups.csv) can be the directly used as an input for TurboID_volcano-plots.script and run within R (version 4.2.3). 


### Packages and R versions
Analysed at Windows 10, R version 4.2.2 or 4.2.3, dplyr_1.1.1, idpr_1.8.0, here_1.0.1, ggplot2_3.4.1, ggrepel_0.9.5, ggprism_1.0.5

