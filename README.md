# Slovenia9YrFUP
Codes used in `Slovenia9YrFUP`

- `Codes/`: codes used in this analysis
- `Plots/`: plots used in this study
- `Results/`: output csv files

## Codes
- `data_preparation.R`
	- Data cleaning
	  - fixing dates 
		- column selection
		- wide -> long
		- inclusion/exclusion criteria
	- Creating subsets of data by age groups
   	- Creating subsets of data for co-test analysis
- `analysis.R`
  - main analysis 
  - clean output from the analysis for tables in the manuscript  
- `plots.R`
  - code used to generate the plots for the manuscript
- `ggvennXY.R`
  - modified `ggvenn()` functions to produce Venn diagrams with 2 lines of information (i.e., number test positive and number with CIN2/3+)

## Plots
- `Figure_1.png`: Cumulative risks of cervical intraepithelial neoplasia grade 2 or worse (CIN2+) and grade 3 or worse (CIN3+) by normal baseline cytology and negative baseline HPV result according to HPV assay used among women ≥30 years old (A, B) and in total study population (C, D). Risk bands are 95% confidence intervals. 
- `Figure_2.png`: Cumulative risks of cervical intraepithelial neoplasia grade 2 or worse (CIN2+) by normal baseline cytology (red), negative baseline HPV result (green), and co-testing (both normal baseline cytology and negative baseline HPV result; blue) according to HPV assay in women ≥30 years old. 
- `Figure_3.png`: Cumulative risks of cervical intraepithelial neoplasia grade 2 or worse (CIN2+) and grade 3 or worse (CIN3+) by HPV genotype during 9 years of follow-up among women ≥30 years old (A, B) and in total study population (C, D). 

## Results
- `Table_1.csv`: Baseline screening test positivity rates by four different clinically validated HPV assays used for HPV testing in women ≥30 years old and in total study population (women 20–64 years old).
- `Table_2.csv`: Cumulative risk of cervical intraepithelial neoplasia grade 2 or worse (CIN2+) and grade 3 or worse (CIN3+) in women with normal baseline cytology result or negative HPV baseline result at 3 years, 6 years, and 9 years of follow-up, stratified by age group.
- `Table_3.csv`: Concordance between four clinically validated HPV assays by number of concurrent baseline HPV-positive results and cervical intraepithelial neoplasia grade 2 or worse (CIN2+) identified over 9-year follow-up in women ≥30 years old and in total study population.
- `Table_S1.csv`: Longitudinal test characteristics using cervical intraepithelial neoplasia grade 2 or worse (CIN2+) and grade 3 or worse (CIN3+) as outcome in women by their baseline cytology and HPV result at 9 years of follow-up in women ≥ 30 years old and in total study population
