#!/usr/bin/env python

import numpy as np
from ase.calculators.calculator import Calculator
from ase.atoms import Atoms

from pyscf import gto, scf, grad, mp

def ase_atoms_to_pyscf(ase_atoms):
  '''Convert ASE atoms to PySCF atom.

  Note: ASE atoms always use A.
  '''
#return [[ase_atoms.get_chemical_symbols(), ase_atoms.get_positions()] for i, atom in enumerate(ase_atoms)]
  return [ [ase_atoms.get_chemical_symbols()[i], ase_atoms.get_positions()[i]] for i in range(len(ase_atoms.get_positions()))]

class PySCF_simple(Calculator):
	name = 'PySCF_simple'
	implemented_properties = ['energies', 'forces']

	def __init__(self, atoms, method, basis):
		self.mol = gto.M(verbose=0)
		self.mol.atom = ase_atoms_to_pyscf(atoms)
		self.mol.basis = basis
		self.mol.build()

		self.method = method
		self.results = {}

	def get_potential_energy(self, atoms=None, force_consistent=False):
		if self.method != 'MP2':
			self.mf = scf.RHF(self.mol)
			energy = self.mf.kernel()

			self.results['energies'] = energy
			return energy

		else:
			self.mf  = scf.RHF(self.mol)
			hf_energy = self.mf.kernel()

			self.mp2 = mp.MP2(self.mf)
			corr_energy   = self.mp2.kernel()

			self.results['energies'] = corr_energy[0] + hf_energy
			return corr_energy[0] + hf_energy


	def get_forces(self, atoms=None):
		if self.method != 'MP2':
			forces  = -1*grad.RHF(self.mf).kernel()

			self.results['forces'] = forces
			return forces

		else:
			self.mf  = scf.RHF(self.mol)
			self.mf.kernel()

			self.mp2 = mp.MP2(self.mf)
			self.mp2.kernel()

			forces = -1*grad.mp2.Grad(self.mp2).kernel()

			self.results['forces'] = forces
			return forces
