# S2F (Sequence to Function)

A protein function prediction tool that takes only the sequences as input

[ToC]

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

If you prefer to setup a `virtualenv` to use S2F, first create it, and then run the script passing the path to the python interpreter within that `virtualenv` to the script:

```bash
./install_requirements.sh ~/.virtualenvs/s2f/bin/python
```

this will also generate a script named `S2F.sh` and configure it accordingly

## Installation

A one-command solution is available to setup the entire environment for S2F (this will download HUNDREDS of GB of datasets, please look below to [configure the installer](#configure-the-installation-script) according to your needs, or the [command line options](#installer-options)). 
You may simply run:
```bash
python S2F.py install
```
You will be asked to confirm the configured options of the installation script. IF you're happy and reply with `y`, the installer will download and process the required databases. This, however, assumes that all requirements listed above are met. 

### Configure the installation script

The easiest way to configure the installation is through the configuration file. The default values are:

```ini
[directories]
installation_directory = ~/.S2F

[commands]
interpro = iprscan
hmmer = phmmer
blastp = blastp
makeblastdb = makeblastdb

[databases]
string_links = download
string_sequences = download
string_species = download
uniprot_sprot = download
uniprot_goa = download
filtered_goa = infer
filtered_sprot = infer

[options]
evidence_codes = experimental
```

| Option | Description |
| --- | --- |
| `installation_directory` | Path to the installation directory for S2F. |
| `interpro` | manually provide the path to the `iprscan` executable to avoid passing this parameter to the other commands every time. |
| `hmmer` | manually provide the path to the `phmmer` executable to avoid passing this parameter to the other commands every time. |
| `blastp` | manually provide the path to the `blastp` command in the system. If not provided, S2F will assume that the executable is available system-wide. |
| `makeblastdb` | manually provide the path to the `makeblastdb` command in the system. If not provided, S2F will assume that the executable is available system-wide. |
| `string_links` | 'manually provide the path to the STRING interactions database, it must be the full path to either `protein.links.full.vX.x.txt.gz` or `protein.links.detailed.vX.x.txt.gz`. If not provided, the installation script will attempt to download the full database using the `wget` command.' |
| `string_sequences` | manually provide the path to the STRING sequences database, it must be the full path to the `protein.sequences.vX.x.fa.gz` file. If not provided, the installation script will attempt to download it using the `wget` command.' |
| `string_species` | manually provide the path to the STRIN species list, it must be the full path to the `species.vX.x.txt` file. If not provided, the installation script will attempt to download it using the `wget` command. |
| `uniprot_sprot` | manually provide the path to the UniProt SwissProt sequences, it must be the full path to the "goa_uniprot_all.gaf.gz" file. If not provided, the installation script will attempt to download it using the `wget` command. |
| `uniprot_goa` | manually provide the path to the UniProt GOA, it must be the full path to the "goa_uniprot_all.gaf.gz" file. If not provided, the installation script will attempt to download it using the `wget` command. |
| `evidence_codes` | manually provide a list of evidence codes that will be used to filter the UniProt GOA. If not provided, S2F will be installed using only experimental evidence codes. |

### Installer options

The command line options for the installation are the following (but we highly recommend having a look at the [configuration file](#configure-the-installation-script) to avoid mistakes):

| Option | Description | Default Value |
| --- | --- | --- |
| `--installation-directory` | Path to the installation directory for S2F. | `~/.S2F` |
| `--config-file` | location of the configuration file that will be created. If not provided, the default configuration file will be loaded. | `s2f.conf` (found in the script's directory) |
| `--interpro` | manually provide the path to the `iprscan` executable to avoid passing this parameter to the other commands every time. | `iprscan` (assumes this is correctly configured in the `PATH` environment variable) |
| `--hmmer` | manually provide the path to the `phmmer` executable to avoid passing this parameter to the other commands every time. | `phmmer` (assumes this is correctly configured in the `PATH` environment variable) |
| `--blastp` | manually provide the path to the `blastp` command in the system. If not provided, S2F will assume that the executable is available system-wide. | `blastp` (assumes this is correctly configured in the `PATH` environment variable) |
| `--makeblastdb` | manually provide the path to the `makeblastdb` command in the system. If not provided, S2F will assume that the executable is available system-wide. | `makeblastdb` (assumes this is correctly configured in the `PATH` environment variable) |
| `--string-links` | 'manually provide the path to the STRING interactions database, it must be the full path to either `protein.links.full.vX.x.txt.gz` or `protein.links.detailed.vX.x.txt.gz`. If not provided, the installation script will attempt to download the full database using the `wget` command.' | `download` |
| `--string-sequences` | manually provide the path to the STRING sequences database, it must be the full path to the `protein.sequences.vX.x.fa.gz` file. If not provided, the installation script will attempt to download it using the `wget` command.' | `download` |
| `--string-species` | manually provide the path to the STRIN species list, it must be the full path to the `species.vX.x.txt` file. If not provided, the installation script will attempt to download it using the `wget` command. | `download` |
| `--uniprot-swissprot` | manually provide the path to the UniProt SwissProt sequences, it must be the full path to the "goa_uniprot_all.gaf.gz" file. If not provided, the installation script will attempt to download it using the `wget` command. | `download` |
| `--uniprot-goa` | manually provide the path to the UniProt GOA, it must be the full path to the "goa_uniprot_all.gaf.gz" file. If not provided, the installation script will attempt to download it using the `wget` command. | `download` |
| `--evidence-codes` | manually provide a list of evidence codes that will be used to filter the UniProt GOA. If not provided, S2F will be installed using only experimental evidence codes. | `experimental` |


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

S2F has 2 configuration files to simplify the installation and prediction processes. The installation file 



### Prediction

```ini
[section]
asd = asd
```