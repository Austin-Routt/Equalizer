# Equalizer

Equalizer is a synthetic gene circuit that compensates cell-to-cell variation due to plasmid copy number variations.

This repository contains code and models that can simulate and compare dosage compensation capacity of the Equalizer and other related circuits.

## MATLAB version

The code has been developed and tested on MATLAB R2020b. Required toolboxes are SimBiology, Statistics and Machine Learning Toolbox, Deep Learning Toolbox (`sumsqr`), and Econometrics Toolbox (`autocorr`).

## System requirements

Operating systems should allow running of MATLAB, and the system requirements can be found at: mathworks.com/support/requirements/matlab-system-requirements.html

## Installation

Download source code (should only take a minute):

`$ git clone https://github.com/stpierrelab/Equalizer`

## Running the code

The code reproduces the Equalizer paper simulation figures, with the provided models and experimental data in respective folders. Depending on the operating systems, some of the more computationally intensive code can take 10 - 30 minutes to finish running. The folders and code have been formatted so that the code can run directly once opened on MATLAB.

## Reference

[A synthetic circuit for buffering gene dosage variation between individual mammalian cells](https://doi.org/10.1038/s41467-021-23889-0)

St-Pierre Lab (stpierrelab.com)
