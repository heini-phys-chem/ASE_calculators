# ASE_calculators

This git rep. contains three ASE calculators from open source QM and ML packages:
 - QMLcode (http://www.qmlcode.org/)
 - PySCF (https://sunqm.github.io/pyscf/)
 - PySCF_simple (HF, DFT, and MP2 for SP, GO and NEB)

## How to use calculators
Put the calculator in the

```~/.local/lib/python3.7/site-packages/ase/calculators ```

directory and use it as follows:

```python
''' ML calculator'''
from ase import io
from ase.atoms import Atoms
from ase.calculators.ml import ML_calculator

mol    = io.read('mol.xyz')
sigma  = 2.5
alphas = np.load('alphas.npy')
X      = np.load('X.npy')
Q      = np.load('Q.npy')

mol.set_calculator(ML_calculator(mol, sigma, alphas, X, Q))

print(mol.get_potential_energy())
```

Alphas, X, and Q are .npy files that were created from the training of the model (using FCHL19).
Sigma is a hyperparameter (also from the training of the model).

```python
''' PySCF calculator'''
from ase import io
from ase.atoms import Atoms
from ase.calculators.pyscf import PySCF

mol = io.read('mol.xyz')
mol.set_calculator(PySCF(atoms=mol, molcell=gto.M(verbose=0), mf_class=scf.RHF, mf_dict={}))

print(mol.get_potential_energy())
```

For more details about PySCF (molcell, mf_class, and mf_dict) please see their documentation

```python
''' PySCF_simple calculator'''
from ase import io
from ase.atoms import Atoms
from ase.calculators.pyscf_simple import PySCF_simple

mol = io.read('mol.xyz')
mol.set_calculator(PySCF_simple(atoms=mol, method='MP2', basis='6-31g*'))

print(mol.get_potential_energy())
print(mol.get_forces())
```

## NEB example
All in one python script:
- Read in reactant and product force field (FF) geometries with ASE
- Optimize both molecules
- Run an NEB calculation
- Plot results

![GitHub Logo](/images/irc.png)

## TODO:

### ml.py:
 - still hardcoded variables which need to be adapted for every run

### pyscf.py:
 - the script was taken from the PySCF homepage (https://sunqm.github.io/pyscf/_modules/pyscf/pbc/tools/pyscf_ase.html) and adapted, resp. finished
 	- added forces
	- changed function ```ase_atoms_to_pyscf(ase_atoms)```
	- add an mp2 wrapper (needed for NEB runs)
- writing trajectories is disabled (JSON error while writing...)

### pyscf_simple.py:
 - TODO: add remaining methods (semi empirical, CC, ...)
