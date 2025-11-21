# Integration of relative proteomics data across multiple conditions with multiRENZO

This approach for the integration of relative proteomics data into constraint-based models of metabolism builds on the previously published approach [RENZO](https://github.com/pwendering/RENZO_paper). While RENZO was designed to integrate ratios of protein abundances for two condition, multiRENZO can integrate pairwise ratios across multiple conditions.



The approach was applied to relative proteomics data measured in _Arabidopsis thaliana_, after shifting the plants to different oxygen conditions:

PUBLICATION LOTRECK ET AL



The function `multiRENZO.m` implements the constraint-based approach. It requires the installation of the [Gurobi solver](https://www.gurobi.com).



The results of the study can be reproduced by running the functions `processProteomicsData.m`, `predictProteinAbundancesO2RENZO.m`, and`evaluatePredictionsMultiRenzo.m`.  Before, the code for the temperature-dependent ecAraCore metabolic model must be downloaded from [here](https://github.com/pwendering/AraTModel) and the Matlab code must be added to the Matlab path. Run the code from within the `multiRENZO/code/` directory or make sure to prepend the `multiRENZO/code/` folder to the Matlab path.
