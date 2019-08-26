# ASE_calculators

This git rep. contains two ASE calculators from open source QM and ML packages:
 - QMLcode (http://www.qmlcode.org/)
 - PySCF (https://sunqm.github.io/pyscf/)

## How to use calculators
Put the calculator in the

```~/.local/lib/python3.7/site-packages/ase/calculators ```

directory and import it with:

```from ase.calculators.ml import ML_calculator ```

or

```from ase.calculators.pyscf import PySCF```


## TODO:

### ml.py:
 - still hardcoded variables which need to be adapted for every run

### pyscf.py:
 - the script was taken from the PySCF homepage (https://sunqm.github.io/pyscf/_modules/pyscf/pbc/tools/pyscf_ase.html) and adapted, reps. finished
 	- added forces
	- changed function ```ase_atoms_to_pyscf(ase_atoms)```
	- added an mp2 wrapper (needed for NEB runs)
- writing trajectories is disabled (JSON error while writing...)
