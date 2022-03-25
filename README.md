#   TAFPred


### Table of Content

- [Setup](#getting-started)
- [Dataset](#Dataset)
- [Prerequisites](#Prerequisites)
- [Download and install code](#download-and-install-code)
- [Demo](#demo)
- [Run with Docker](#Run-with-Docker)
- [References](#References) 

# Getting Started
 

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

 ## Dataset
The dataset can be found in the dataset directory. The dataset is collected from [1].


## Prerequisites

You would need to install the following software before replicating this framework in your local or server machine.

 ```
Python version 3.7.4

Poetry version 1.1.13

You can install poetry by running the following command:

curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -

To configure your current shell run `source $HOME/.poetry/env`


```
  
## Download and install code

- Retrieve the code

```
github clone 

```

## Demo

To run the program, first install all required libraries by running the following command:

```
poetry install

```

Then execute the following command to run TAFPred from the script directory.

```
cd tafpred
poetry run python TAFPred.py 

```

- Finally, check **output** folder for results. The output directory contains predicted lebels with probabilities for each residues.


## Run with Docker

- Build the docker image from Dockerfile.
```
export UID=$(id -u)
export GID=$(id -g)
docker build --build-arg USER=$USER \
             --build-arg UID=$UID \
             --build-arg GID=$GID \
             --build-arg PW=asdf \
             -t TAFPred\
             -f Dockerfile.txt\
             .
```

- Mount the Output direcotry in the Docker Container and run it.

```
docker run -ti  -v /$(pwd)/Output:/home/$USER/output TAFPred:latest
```

- Then, run following python command from the root directory.
```
source $HOME/.poetry/env
cd script
poetry run python TAFPred.py 
```

- Finally, check **output** folder for results. The output should be available in both host and docker. The output directory contains predicted lebels with probabilities for each residues.
 

## Authors

Md Wasi Ul Kabir, Md Tamjidul Hoque. For any issue please contact: Md Tamjidul Hoque, thoque@uno.edu 

## References

1. Zhang, Tuo, Eshel Faraggi, and Yaoqi Zhou. “Fluctuations of Backbone Torsion Angles Obtained from NMR-Determined Structures and Their Prediction: Fluctuations of Backbone Torsion Angles.” Proteins: Structure, Function, and Bioinformatics 78, no. 16 (December 2010): 3353–62. https://doi.org/10.1002/prot.22842.

