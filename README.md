# S2F (Sequence to Function)

A protein function prediction tool that takes only the sequences as input

## Requirements

S2F relies on the following software:

* BLAST (`blastp` and `makeblastdb`)
* InterPro (`iprscan`)
* HMMer (`phmmer`)
* Python 3

S2F relies on several Python 3 packages. A script is provided that will install all the requirements simply by
running:

```bash
./install_requirements.sh
```

If you prefer to setup a virtualenv to use S2F, first create the virtual env, and then run the script passing the path
to the python interpreter within that virtualenv to the script:

```bash
./install_requirements.sh ~/.virtualenvs/s2f/bin/python
```

this will also generate a script named `S2F.sh` and configure it accordingly

## Installation

A one-command solution is available to setup the entire environment for S2F. 
You may simply run:
```bash
./S2F.sh install
```
and the installer will download and process the required databases. This, 
however, assumes that all the requirements are met  

## How to make a prediction

### Using a configuration file

Due to the number of parameters that can be set for S2F, it is often useful to
create a configuration file, which should then be provided to the `predict` 
command using the `--run-config` argument:

```bash
./S2F.sh predict --run-config my_organism.conf
```

### Using a command-line

```bash
./S2F.sh predict --run-config my_organism.conf
```


## Configuration Files

Fully explained options are available in the sample files

### Installation 

```ini
[section]
asd = asd
```

### Prediction

```ini
[section]
asd = asd
```