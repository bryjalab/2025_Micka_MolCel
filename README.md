# Multisite phosphorylation of intrinsically disordered region of DVL facilitates Wnt signaling

## Authors:
Miroslav Micka, Jitender Kumar, Petra Paclíková, Zuzana Hayek, Kateřina Hanáková, Cherine Bechara, Hana Plešingerová, Sara Bologna, Ondrej Šedo, Elise Del Nero, Kristína Gömöryová, Tomáš Gybeľ, Marek Kravec, Zbyněk Zdráhal, Konstantinos Tripsianes**, Vítězslav Bryja*

*Corresponding author: bryja@sci.muni.cz

**Corresponding author: kostas.tripsianes@ceitec.muni.cz

## Reproducing Charge calculations
First run script "charge_functions.R". Here we edit functions from idpr R package (McFadden and Yanowitz, 2022), to enable charge calculation of phosphorylated proteins.
Next run script "charge_calculation.R". We either read analyzed protein sequences directly in the script or with read them from a text file. X stands for phosphoserine and Z stands for phosphothreonine in an analysed sequence.

### Packages and R versions
Analysed at Windows 10, R version 4.2.2, dplyr_1.1.1, idpr_1.8.0, here_1.0.1, ggplot2_3.4.1

