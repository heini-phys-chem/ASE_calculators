from ase import io
from ase.neb import NEB
from ase.neb import NEBTools
from ase.calculators.pyscf_simple import PySCF_simple
from ase.calculators.gaussian import Gaussian
from ase.atoms import Atoms
from ase.optimize.lbfgs import LBFGS
from ase.optimize import FIRE

from pyscf import gto, scf, grad, mp

import numpy as np

import matplotlib.pyplot as plt

# Read initial and final states:
print(" -> read xyz")
initial = io.read('react.xyz')
final = io.read('prod.xyz')

initial.set_calculator(PySCF_simple(initial, method='MP2', basis='6-31g*'))
final.set_calculator(PySCF_simple(final, method='MP2', basis='6-31g*'))

print(" -> opt react prod")
dyn_react = LBFGS(initial)
dyn_react.run(fmax=0.05)
#
dyn_prod = LBFGS(final)
dyn_prod.run(fmax=0.05)

# write output geometries of GO
io.write("initial.xyz", initial)
io.write("final.xyz", final)

# Make a band consisting of N images:
images	= [initial]
images += [initial.copy() for i in range(17)]
images += [final]

print(" -> set calculator")
for image in images:
    image.set_calculator(PySCF_simple(atoms=image, method='MP2', basis='6-31g*'))

neb = NEB(images,climb=True, k=0.6)
nebTools = NEBTools(images)

neb.interpolate('idpp')

print(" -> start neb run")
opt = FIRE(neb)
opt.run(fmax=0.05)

print(nebTools.get_barrier())

# get IRC data
Ef, dE = nebTools.get_barrier()
max_force = nebTools.get_fmax()
x, y, x_fit, y_fit, forces = nebTools.get_fit()

# save IRC data
np.save("x_claisen_mp2.npy", x)
np.save("y_claisen_mp2.npy", y)
np.save("x_claisen_fit_mp2.npy", x_fit)
np.save("y_claisen_fit_mp2.npy", y_fit)

# write NEB guess of interpolation
for i,image in enumerate(images):
    out_name = "claisen_{:03d}.xyz".format(i)
    io.write(out_name,	image)

# plot IRC
y2     *= 23.06
y_fit  *= 23.06

plt.plot(x_fit,y_fit, color='C0', label='MP2')
plt.scatter(x,y, marker='x', color='k', lw=2)

leg = plt.legend(fontsize=30)
leg_lines = leg.get_lines()

plt.setp(leg_lines, linewidth=4)
plt.xlabel("IRC", fontsize=fs)
plt.ylabel("Energy [kcal/mol]", fontsize=fs)
plt.title(r"PySCF - S$_N$2", fontsize=fs)
plt.tick_params(labelsize=30)

plt.show()
