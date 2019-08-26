# ASE_calculators
Put the calculator in the

```~/.local/lib/python3.7/site-packages/ase/calculators ```

directory and import it with:

```from ase.calculators.ml import ML_calculator ```

or
```from ase.calculators.pyscf import PySCF```



TODO:

ml.py:
 - 3-4 variables are hardcoded and need to be adapted for every run

pyscf.py:
- writing trajectories is disabled (JSON error while writing...)
