# ASE_calculators

This git rep. contains two ASE calculators from open source QM and ML packages:
 - QMLcode (http://www.qmlcode.org/)
 - PySCF (https://sunqm.github.io/pyscf/)

## How to use calculators
Put the calculator in the

```~/.local/lib/python3.7/site-packages/ase/calculators ```

directory and use it as follows:

```python
''' ML calculator'''
from ase.calculators.ml import ML_calculator

mol    = io.read('mol.xyz')
sigma  = 2.5
alphas = np.load('alphas.npy')
X      = np.load('X.npy')
Q      = np.load('Q.npy')

mol.set_calculator(ML_calculator(mol, sigma, alphas, X, Q))
```

Alphas, X, and Q are .npy files that were created from the training of the model (using FCHL19).
Sigma is a hyperparameter (also from the training of the model).

```python
''' PySCF calculator'''
from ase.calculators.pyscf import PySCF
mol = io.read('mol.xyz')
mol.set_calculator(PySCF(atoms=mol, molcell=gto.M(verbose=0), mf_class=scf.RHF, mf_dict={}))
```

For more details about PySCF (molcell, mf_class, and mf_dict) please see their documentation

## TODO:

### ml.py:
 - still hardcoded variables which need to be adapted for every run

### pyscf.py:
 - the script was taken from the PySCF homepage (https://sunqm.github.io/pyscf/_modules/pyscf/pbc/tools/pyscf_ase.html) and adapted, reps. finished
 	- added forces
	- changed function ```ase_atoms_to_pyscf(ase_atoms)```
	- added an mp2 wrapper (needed for NEB runs)
- writing trajectories is disabled (JSON error while writing...)
