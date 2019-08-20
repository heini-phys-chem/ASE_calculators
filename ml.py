import numpy as np
#sys.path.insert(0, "/home/heinen/qml_develop/build/lib.linux-x86_64-2.7")
import qml
from qml.math import cho_solve
from qml.representations import generate_fchl_acsf
from qml.kernels import get_atomic_local_gradient_kernel
from qml.kernels import get_atomic_local_kernel

from ase.calculators.general import Calculator
from ase.atoms import Atoms

convback    = 1.0
convback_E  = 1.0

class ML_calculator(Calculator):
	name = 'ML_Calculator'
	implemented_properties = ['energy', 'forces']

	def __init__(self, atoms, sigma, alphas, X, Q):
		self.alphas = alphas
		self.X			= X
		self.sigmas = sigma
		self.Q			= Q
		self.nAtoms = atoms.get_number_of_atoms()

	def get_potential_energy(self,atoms=None,force_consistent=False):
		x = []
		disp_x = []
		q = []

#		x1 = generate_fchl_acsf(atoms.get_atomic_numbers(), atoms.get_positions(), gradients=False, pad=9, elements=[1,6,7,9,17,35])
		x1 = generate_fchl_acsf(atoms.get_atomic_numbers(), atoms.get_positions(), gradients=False, pad=self.nAtoms)
		x.append(x1)
		q.append(atoms.get_atomic_numbers())

		Xs  = np.array(x)
		Qs  = q 

		Kse = get_atomic_local_kernel(self.X, Xs, self.Q, Qs, self.sigmas)
		energy = (float(np.dot(Kse, self.alphas)))*convback_E

		return energy

	def get_forces(self, atoms=None):
		x = []
		disp_x = []
		q = []

#		(x1, dx1) = generate_fchl_acsf(atoms.get_atomic_numbers(), atoms.get_positions(), gradients=True, pad=9, elements=[1,6,7,9,17,35])
		(x1, dx1) = generate_fchl_acsf(atoms.get_atomic_numbers(), atoms.get_positions(), gradients=True, pad = self.nAtoms)
		x.append(x1)
		disp_x.append(dx1)
		q.append(atoms.get_atomic_numbers())

		Xs  = np.array(x)
		dXs = np.array(disp_x)
		Qs  = q 

		Ks  = get_atomic_local_gradient_kernel(self.X, Xs, dXs, self.Q, Qs, self.sigmas)
		self.fYs = np.dot(Ks, self.alphas)
		Fss = self.fYs.reshape((self.nAtoms,3))*convback

		return Fss 
