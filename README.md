# malaria-mosquito-analysis
Cleaning and preliminary analysis of mosquito and qPCR data.

Part of the [Spatial R21 project](https://sites.duke.edu/taylorlab/projects/#MolEpi), which investigates the molecular epidemiology of *P. falciparum* transmission in Western Kenya.

## Preliminary data management

### Cleaning data
- Standardize column names and sample IDs of all datasets.
- Reformat the descriptive *Anopheles* dataset such that each row corresponds to a single mosquito.
- Censor the qPCR dataset by Ct values of human β-tubulin (Hbtub) and *Plasmodium* (Pfr364) targets, as well as by parasitemia.

### Validating data
- Correct formatting errors in household and sample IDs.
- Correct inconsistencies in village, household, and sample IDs within each sample.
- Correct inconsistencies in descriptions of each mosquito.
- Remove duplicate entries of the same sample.

### Merging data
- Create a dataset pairing mosquito descriptions with corresponding qPCR results.

## Collaborators
- [Kelsey Sumner](https://github.com/kelseysumner/taylorlab)
