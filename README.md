# Blind_docking
Updating the scripts to utilize [GNINA](https://github.com/gnina/gnina) (both with and without CNN scoring functions).

Adds symmetry-corrected RMSD calculations (using `obmrs` from [Open-Babel](https://github.com/openbabel/openbabel))

## Prerequisites
- [GNINA](https://github.com/gnina/gnina) (if you don't want to do a full install, can use the binary)
- [Open-Babel](https://github.com/openbabel/openbabel)

## Installation

```bash
mamba create -n blind_docking python==3.10 numpy biopandas plumbum -c conda-forge
conda activate blind_docking
```

## Example

```bash
python equi_bind_diffdock_box_gnina.py
```
Then results including original output poses and RMSDs will be in the `results` folder.
