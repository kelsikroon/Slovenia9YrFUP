# Slovenia9YrFUP
Codes used in `Slovenia9YrFUP`

- `Codes/`: codes used in this analysis
- `Plots/`: plots used in this study
- `Results/`: output csv files

## Codes
- data_preparation_Jun2024.R
	- Data cleaning
	  - fixing dates 
		- column selection
		- wide -> long
		- inclusion/exclusion criteria 
- analysis_Jun2024.R
  - main analysis 
  - clean output from the analysis for tables in the manuscript  
- plots_Jun2024.R
  - code used to generate the plots for the manuscript
- ggvennXY.R
  - modified `ggvenn()` functions to produce Venn diagrams with 2 lines of information (i.e., number test positive and number with CIN2/3+)
