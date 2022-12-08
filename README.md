# scSHARP

Repository for "Consensus Label Propagation with Graph Convolutional Networks for Single-Cell RNA Sequencing Cell Type Annotation" extended abstract submission for Learning on Graphs Conference and full paper available at https://doi.org/10.1101/2022.11.23.517739.

## Installation (pypi)
### Create a new conda environment
```
conda create -n <env name> python=3.9
conda activate <env name>
```

### Install torch

Linux with GPU:
```
conda install pytorch pytorch-cuda=11.7 -c pytorch -c nvidia
``` 
Mac OS:
```
conda install pytorch -c pytorch
```

### Install torch geometric
```
conda install pyg -c pyg
```

### Install scSHARP

```
pip install scSHARP
```

Installation of torch and torch geometric required prior to pip install.

See demo.ipynb for example work flow.
Please input raw counts.
