from ase import io
from ase.neb import NEB
from ase.neb import NEBTools
from ase.calculators.pyscf_simple import PySCF_simple
from ase.atoms import Atoms
from ase.optimize.lbfgs import LBFGS
from ase.optimize import FIRE

import numpy as np
import matplotlib.pyplot as plt

from pyscf import gto, scf, grad, mp

# Read initial and final states:
initial = io.read('react.xyz')
final = io.read('prod.xyz')

# set calculator
initial.set_calculator(PySCF_simple(initial, method='MP2', basis='6-31g*'))
final.set_calculator(PySCF_simple(final, method='MP2', basis='6-31g*'))

# Opt reactant & product
dyn_react = LBFGS(initial)
dyn_react.run(fmax=0.05)

dyn_prod = LBFGS(final)
dyn_prod.run(fmax=0.05)

# Make a band consisting of N + 2 images:
N = 17
images	= [initial]
images += [initial.copy() for i in range(N)]
images += [final]

# set calculators for images
for image in images:
    image.set_calculator(PySCF_simple(atoms=image, method='MP2', basis='6-31g*'))

# set NEB
neb = NEB(images, climb=True, k=0.6)
nebTools = NEBTools(images)

neb.interpolate('idpp')

# start NEB run
opt = FIRE(neb)
opt.run(fmax=0.05)

# get IRC data
Ef, dE = nebTools.get_barrier()
max_force = nebTools.get_fmax()
x, y, x_fit, y_fit, forces = nebTools.get_fit()

# save IRC data
np.save("x.npy", x)
np.save("y.npy", y)
np.save("x_fit.npy", x_fit)
np.save("y_fit.npy", y_fit)

# plot IRC data
y     *= 23.06
y_fit *= 23.06

plt.plot(x_fit, y_fit, color='C0', label='MP2')
plt.scatter(x, y, marker='x', color='k', lw=2)

leg = plt.legend(fontsize=30)
leg_lines = leg.get_lines()

plt.setp(leg_lines, linewidth=4)
plt.xlabel("IRC", fontsize=fs)
plt.ylabel("Energy [kcal/mol]", fontsize=fs)
plt.title(r"PySCF - S$_N$2", fontsize=fs)
plt.tick_params(labelsize=30)

plt.show()
