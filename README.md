# malaria-mosquito-analysis
Cleaning and preliminary analysis of mosquito and qPCR data.

Part of the [Spatial R21 project](https://sites.duke.edu/taylorlab/projects/#MolEpi), which investigates the molecular epidemiology of *P. falciparum* transmission in Western Kenya.

## Preliminary data management

### Cleaning data
- Standardize column names and sample IDs of all data sets.
- Reformat the descriptive *Anopheles* data set such that each row corresponds to a single mosquito.
- Censor the qPCR data set by Ct values of human β-tubulin (Hbtub) and *Plasmodium* (Pfr364) targets.

### Validating data
- Correct inconsistencies in village, household, and sample IDs within each sample.
- Correct inconsistencies in descriptions of each mosquito.

### Merging data
- Create a data set pairing mosquito descriptions with corresponding qPCR results.

## Collaborators
- [Kelsey Sumner](https://github.com/kelseysumner/taylorlab)
