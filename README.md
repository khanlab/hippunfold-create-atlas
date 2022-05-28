# create probabilistic subfield atlas from **freesurfer** segmentations

This uses Freesurfer segmentations to generate a probabilistic (and maxprob) hippunfold atlas

Pre-requisites:
 0. bids dataset (e.g. ds002168 from openneuro) with T1w and hires T2w
 1. hippunfold run on this dataset
 2. freesurfer 7 `recon-all` run on this dataset
 3. freesurfer 7 `segmentHA_T2.sh` run on the T2s (e.g. `segmentHA_T2.sh sub-01 PATH_TO_T2 T2 0` )


This workflow will transform the freesurfer segmentations to corobl, sample the subfield labels, convert to cifti, generate a probablistic segmentation and a maximum-probability segmentation

Note: resampling to unfoldiso space requires the latest dev version of connectome-workbench with the `-bypass-sphere-check` flag 
