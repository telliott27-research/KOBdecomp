# KOBdecomp

Package for running customized KOB decompositions.

This package was developed for the academic article [_Softer Policing or the Institutionalization of Protest? Decomposing Changes in Observed Protest Policing Over Time_](https://www.journals.uchicago.edu/doi/10.1086/719001) by Thomas Elliott, Jennifer Earl, Thomas V. Maher, and Heidi Reynolds-Stenson. Please refer to that paper, specifically the methodological appendix, for more details on how the model works.

The package provides a method for applying the modified Kitagawa-Oaxaca-Blinder (KOB) decomposition used in the article to binary and multinomial outcomes with continuous grouping variables (e.g. time). 

### Installing the package

The easiest way to install the package is to use the `install_github` function from the `remotes` package:

```R
remotes::install_github("telliott27-research/KOBdecomp")
```

For binary outcome models, the primary function for analysis is `gen_decomposed_results_probit`.

For multinomial outcome models, the primary function for analysis is `gen_decomposed_results_multinom`.