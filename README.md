#   TAFPred
TAFPred: Torsion Angle Fluctuations Prediction from Protein Sequences

### Table of Content
- [TAFPred](#tafpred)
    - [Table of Content](#table-of-content)
- [Getting Started](#getting-started)
  - [Datasets](#datasets)
    - [Prerequisites](#prerequisites)
  - [Download and install code](#download-and-install-code)
  - [Authors](#authors)
  - [References](#references)

# Getting Started
 

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

 ## Datasets
The dataset can be found in the Dataset/completeDataset directory. The dataset is collected from [1].


### Prerequisites

We have tested Dispredict4.0 on Ubuntu 20.04. You would need to install the following software before replicating this framework in your local or server machine. 

1. pyenv latest version
    ```
    curl https://pyenv.run | bash
    exec $SHELL
    ```
    For more details, visit: https://github.com/pyenv/pyenv-installer

2. Python version 3.9.5

    ```
    pyenv install miniconda3-3.9-4.10.3
    pyenv local miniconda3-3.9-4.10.3 
    ```

    For more details, visit: https://github.com/pyenv/pyenv

3. Poetry version 1.3.2

    ```
    curl -sSL https://install.python-poetry.org | python3 - --version 1.3.2
    ```
    For more details, visit: https://python-poetry.org/docs/

4. Docker
    The feature extraction part depends on many other tools. So we created a docker image to easily extract features without any setup. It might take quit some time to get all the features.

  
## Download and install code

- Retrieve the code

```
github clone https://github.com/wasicse/TAFPred.git

```

To run the program, first install all required libraries by running the following command:

```
poetry install

```

Then execute the following command to run TAFPred from the script directory on the example dataset. You need to change the input of the Dataset/example directory to get prediction for new protein sequences.

```
cd script && poetry run python run_tafpred.py -f "taffeaturesv1" -o "./output/"

```

- Finally, check **output** folder for results. The output directory contains predicted lebels with probabilities for each residues.


## Authors

Md Wasi Ul Kabir, Md Tamjidul Hoque. For any issue please contact: Md Tamjidul Hoque, thoque@uno.edu 

## References

1. Zhang, Tuo, Eshel Faraggi, and Yaoqi Zhou. “Fluctuations of Backbone Torsion Angles Obtained from NMR-Determined Structures and Their Prediction: Fluctuations of Backbone Torsion Angles.” Proteins: Structure, Function, and Bioinformatics 78, no. 16 (December 2010): 3353–62. https://doi.org/10.1002/prot.22842.
