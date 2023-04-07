# Integrase-Excisionase Characterization in TX-TL
[![Paper](https://img.shields.io/badge/ACS--SynBio-2c00534-brightgreen)](https://pubs.acs.org/doi/10.1021/acssynbio.2c00534)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/license/mit/)

Data and models repository for paper on [Characterization of integrase and excisionase activity in cell-free protein expression system using a modeling and analysis pipeline](https://pubs.acs.org/doi/pdf/10.1021/acssynbio.2c00534). 

# Installation

1. Run `pip install -r requirements.txt` from your terminal.

Or, if you prefer to install everything manually, you can:
1. Install biocrnpyler: `pip install biocrnpyler>=1.1.1`
2. Install bioscrape: `pip install bioscrape>=1.2.0`
3. If Step 2 fails, you can clone the [Bioscrape](https://github.com/biocircuits/bioscrape/) repository and manually install the package. Run `python setup.py install` to install bioscrape once you meet all the [installation requirements](https://github.com/biocircuits/bioscrape/wiki/Installation).
4. Install autoreduce: `pip install autoreduce`
5. Install seaborn, corner, and tqdm for plotting: `pip install seaborn corner tqdm`

You can get started with the models using the [Getting Started](https://github.com/ayush9pandey/integrase-excisionase-characterization/blob/main/Getting%20started.ipynb) notebook. For all of the modeling details, you can run the individual notebooks for integrase and excisionase modeling or data analysis.

If you run into any issues, please feel free to raise an [issue](https://github.com/ayush9pandey/integrase-excisionase-characterization/issues).

# How to navigate this repository

## CRN models and their dimensionality reduction
We use BioCRNpyler to generate chemical reaction network models and AutoReduce to find reduced models. The generation of all models are in the `modeling/` directory in this repository.

## Experimental data
All experimental data used in the paper is available under the `experimental_data` directory. The plasmids used in this paper are available from [Addgene](https://www.addgene.org/browse/article/28233404/).

## Data analysis
We use Bioscrape to identify the most sensitive model parameters and identify them from the experimental data. The data analysis pipeline is available in the `data_analysis` directory.

## SBML models
All generated SBML models are available in the `sbml_files` directory. You can store other generated models here as well.

## Output files
All plots and visualizations are stored in the `outputs` directory. 

# Model parameters

All model parameters are available on the [Wiki page](https://github.com/ayush9pandey/integrase-excisionase-characterization/wiki/) of this Github repository.

# Contact

Contact [Ayush Pandey](https://www.its.caltech.edu/~apandey/) (apandey at caltech dot edu) for any other questions.
