# 2020_Mayr
Reproducibility repo accompanying  Mayr et al. _"Integrative analysis of cell state changes in lung fibrosis with peripheral protein biomarkers"_ in EMBO Molecular Medicine.

<p align="center">
<img src="https://github.com/theislab/2020_Mayr/blob/master/figure%201/graph%20abstract.jpg" alt="drawing" width="750">
</p>

Repo for the preprint at medRxiv  accompanying  Mayr et al. _"Integrated single cell analysis of human lung fibrosis resolves cellular origins of predictive protein signatures in body fluids"_ can be found in the folder medRxiv.

Download all data from [here](https://drive.google.com/uc?export=download&id=13vf6Fcy6cCJUuGvbnj5sQDhayLRq7op1) and subsequently uncompress the files.  
The two anndata objects have following layers of information
- munich_dataset.5had:  
  raw counts in .layers["counts"], normalized/log-transformed counts in .raw.X and scaled/corrected counts in .X
- integrated_human_dataset.h5ad: contains raw counts of all integrated data sets in .X
