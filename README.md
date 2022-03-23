# create probabilistic subfield atlas from manual segmentations

This example uses the ASHS Magdeburg training data to generate a probabilistic segmentation.

You need to first run hippunfold on the training data to generate hippocampal surfaces, which will be used to sample the manual segmentations.

This workflow will transform the segmentations to corobl, sample the subfield labels, convert to cifti, generate a probablistic segmentation and a maximum-probability segmentation


ASHS Magdeburg training data can be downloaded on NITRC as `IKND Magdeburg Young Adult 7T Atlas` here: https://www.nitrc.org/frs/?group_id=370
