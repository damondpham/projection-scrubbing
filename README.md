# projection-scrubbing-paper 

Public repository for code and selected results for the projection scrubbing paper. 

`/code`: R source code for the analysis.
* `/analysis`: Analysis scripts contributing to manuscript. Numbered in order of execution.
* `/exploratory`: Analysis scripts which did not directly contribute to the manuscript.

`/data`: Extra data beyond the HCP, which was stored elsewhere. Includes the parcellation and some intermediate analysis results. In this repo there's just the list of subjects used for prediction and the table of parcels.

`/analysis-results`: Saved results for a single scan: the one shown in Figure 1. All files except the FC results are exact copies; the FC results have been subsetted
to only show the selected scrubbing methods, since it's too large for GitHub otherwise.

Replicating the study would require adding missing data files to `/data`, pointing `/code/0_Secret-FilePaths.txt` to the HCP data and other file paths, and then running the scripts in `/code/analysis` in order. 