# Slovenia9YrFUP
Codes used in `Slovenia9YrFUP`

- `Codes/`: codes used in this analysis
- `Plots/`: plots used in this study
- `Results/`: output csv files

## Codes
- data_preparation.R
	- Data cleaning
	  - fixing dates 
		- column selection
		- wide -> long
		- inclusion/exclusion criteria
	- Creating subsets of data by age groups
   	- Creating subsets of data for co-test analysis
- analysis.R
  - main analysis 
  - clean output from the analysis for tables in the manuscript  
- plots.R
  - code used to generate the plots for the manuscript
- ggvennXY.R
  - modified `ggvenn()` functions to produce Venn diagrams with 2 lines of information (i.e., number test positive and number with CIN2/3+)
