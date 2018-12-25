# Tox21_MI


This project investigates using statistics from information theory, such as *Shannon entropy* and *mutual information*, to analyze dose-response data from drug screens. 

The results of this analysis is summarized in [this poster](Sankar_Poster.pdf)


User should only change the filepath and filenames as indicated in the script.

The following data is required as input (please organize the data in ONE folder):

	* Clump Results from the ANOVA process
	* Tox21 Assay chemical concentration-response data from qHTS screens (e.g. "AssayData.txt")
	* Corresponding metadata for the assay data (contians std dev. of negative controls

Pipeline will execute the following automatically and in this order:

	1. CAS IDs from Clump Results will be filtered to CAS IDs with homogeneous profiles (i.e. clump = 1)
		a. A list of these CAS IDs will be written out to "CASRN.txt"

	2. Concentration-response data for filtered CAS IDs will be subsetted and permuted 30 times.
		a. Filename "AssayData_0.txt" -> orignal data
		b. Filename "AssayData_1.txt", "AssayData_2.txt", ... , "AssayData_30.txt" -> permuted data

	3. MI, Pearson correlation, and Spearman correlation will be computed between all filtered CAS IDs¹
		a. Filename "AssayData_Spline_0.txt" -> data processing of original data ("AssayData_0.txt")
		b. Filename "AssayData_Spline_1.txt", "AssayData_Spline_2.txt", ... , "AssayData_Spline_30.txt" -> data processing of permuted data

	4. MI, Pearson correlation, and Spearman correlation estimates from permuted data are collated
		a. Collated data written out to "Permuted_Est.txt"

	5. Significance (p-values) of MI, Pearson, and Spearman measures obtained by permutation testing 
	
	6. Plot generation
		a. Scatter plots of Pearson vs. MI and Spearman vs. MI  using estimates from original data
			i. written out as "Pearson_MI_Scatter.jpg" and "Spearman_MI_Scatter.jpg", respectively

		b. Histograms of MI, Pearson correlation coefficient, and Spearman correlation coefficient distributions from permuted data
			i. written out as "Perm_MI.jpg", "Perm_Pearson.jpg", and "Perm_Spearman.jpg", respectively

		c. Concentration-response profiles of 10 potentially non-linear relationships between chemicals
			i. Significant MI value and insignificant Pearson correlation coefficient graphs written out to Pearson folder
			ii. Significant MI value and insignificant Spearman correlation coefficient graphs written out to Spearman folder

¹Interpolating splines are used to estimate biological response values at specific concentration intervals along each individual profile for every chemical. This strategy accounts for missing data and unequal concentration spacing in raw qHTS data.

NOTE: Input data remains in the home directory. Processed data are written out to Data folder. All plots are written out to Graphics folder.






