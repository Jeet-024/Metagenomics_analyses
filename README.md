# Metagenomics_analyses
Microbial genomics analyses using DADA2 pipeline. Samples collected from different sources and taxonomical classification done. MiSeq amplicon sequencing done on V3-V4 region. 
Contains the R scripts necessary to run the pipeline as well as the results that were given as output for the analyses.
<img width="1480" height="975" alt="image" src="https://github.com/user-attachments/assets/20474c81-e328-4193-a04d-64acc861247d" />
The above image is a krona plot to show the classification of the bacteria population present across the 25 samples.


<img width="1511" height="965" alt="image" src="https://github.com/user-attachments/assets/71b87d61-f078-4260-a770-8ae38b0c5ce7" />
The above image is an interactive krona plot presenting the abundance of bacteria in all the samples individually.

Command to run the interactive krona plot generator:
```
ktImportText krona_nested_samples.txt -o krona_nested_plot.html
```
