# EEG_NNK_GRAPH_CONSTRUCTION
NNK Graph Construction for EEG Signals


## Paper:
https://ieeexplore-ieee-org.libproxy1.usc.edu/document/9909594
E. Bozkurt and A. Ortega, "Non-Negative Kernel Graphs for Time-Varying Signals Using Visibility Graphs," 2022 30th European Signal Processing Conference (EUSIPCO), Belgrade, Serbia, 2022, pp. 1781-1785, doi: 10.23919/EUSIPCO55093.2022.9909594

## Description:
This work concerns the Non-negative kernel graph construction technique for multi-channel EEG signals.
1. It is compatible with any similarity/distance metric
2. It works in a sliding window manner
3. It is a time-domain approach
4. Time shifts can be pre-computed, and possible delays between the time windows are addressed

## Dataset:
This work is compatible with any multi-channel EEG datasets. Some of the public datasets that are explored in this study :
1. IID vs TDV in MUSIC/REST states:
   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7163311/
   Sareen E, Singh L, Varkey B, Achary K, Gupta A. EEG dataset of individuals with intellectual and developmental disorder and healthy controls under rest and music stimuli. Data Brief. 2020 Apr 7;30:105488. doi: 10.1016/j.dib.2020.105488. PMID: 32322626; PMCID: PMC7163311.
   
2. Neonatal Seizure:
   https://www.nature.com/articles/sdata201939#citeas
   Stevenson, N., Tapani, K., Lauronen, L. et al. A dataset of neonatal EEG recordings with seizure annotations. Sci Data 6, 190039 (2019). https://doi.org/10.1038/sdata.2019.39

## Dependencies:
Please download NNK graph construction: https://github.com/STAC-USC/NNK_graph_construction

