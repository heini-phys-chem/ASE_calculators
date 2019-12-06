#!/usr/bin/env python3

import numpy as np
from ase.calculators.calculator import Calculator
from ase.atoms import Atoms

from pyscf import gto, scf, grad, mp, dft

convert_energy    =  27.2114
convert_forces    =  -27.2114 / 0.529177
convert_positions =  0.529177

def ase_atoms_to_pyscf(ase_atoms):
  return [ [ase_atoms.get_chemical_symbols()[i], ase_atoms.get_positions()[i]] for i in range(len(ase_atoms.get_positions()))]

class PySCF_simple(Calculator):
	name = 'PySCF_simple'
	implemented_properties = ['energies', 'forces']

	def __init__(self, atoms, method='HF', basis='sto-3g', xc='pbe'):
		self.basis  = basis
		self.method = method
		self.xc     = xc
		self.results = {}

		#self.hf_scanner  = gto.M().set(verbose=0).apply(scf.RHF).as_scanner()
		#self.mp2_scanner = gto.M().set(verbose=0).apply(scf.RHF).apply(mp.MP2).as_scanner()

	def get_potential_energy(self, atoms=None, force_consistent=False):
		if self.method == 'HF':
			#energy = self.hf_scanner(gto.M(atom=ase_atoms_to_pyscf(atoms), basis=self.basis))
			mf = scf.RHF(gto.M(atom=ase_atoms_to_pyscf(atoms), basis=self.basis, charge=0, verbose=0))

			energy  = mf.kernel()
			energy *= convert_energy

			self.results['energy'] = energy
			return energy

		elif self.method == "DFT":
			mf     = dft.RKS(gto.M(atom=ase_atoms_to_pyscf(atoms), basis=self.basis, charge=0, verbose=0))
			mf.xc = self.xc
			energy = mf.kernel()
			energy *= convert_energy

			self.results['energy'] = energy
			return energy

		elif self.method == "MP2":
			#energy = self.mp2_scanner(gto.M(atom=ase_atoms_to_pyscf(atoms), basis=self.basis))
			mf      = scf.RHF(gto.M(atom=ase_atoms_to_pyscf(atoms), basis=self.basis, charge=0, verbose=0))
			e_hf    = mf.kernel()

			mp2     = mp.MP2(mf)
			e_mp    = mp2.kernel()[0]

			energy  = e_hf + e_mp
			energy *= convert_energy

			self.results['energy'] = energy
			return energy

		#else:
		#	except Exception as e: print("wrong method chosen. available methdos: HF, DFT, and MP2")


	def get_forces(self, atoms=None):
		if self.method == 'HF':
			mf       = scf.RHF(gto.M(atom=ase_atoms_to_pyscf(atoms), basis=self.basis, charge=0, verbose=0)).run()
			gradient = grad.RHF(mf).kernel()
			forces   = gradient * convert_forces

			self.results['forces'] = forces
			return forces

		elif self.method == "DFT":
			mf       = dft.RKS(gto.M(atom=ase_atoms_to_pyscf(atoms), basis=self.basis, charge=0, verbose=0)).run()
			mf.xc = self.xc
			#gradient = grad.rks(mf).kernel()
			gradient = mf.nuc_grad_method().kernel()
			forces = gradient * convert_forces

			self.results['forces'] = forces
			return forces

		elif self.method == "MP2":
			mf       = scf.RHF(gto.M(atom=ase_atoms_to_pyscf(atoms), basis=self.basis, charge=0, verbose=0)).run()
			mp2      = mp.MP2(mf).run()
			gradient = mp2.nuc_grad_method().kernel()

			forces   = gradient * convert_forces

			self.results['forces'] = forces
			return forces

		#else:
		#	except Exception as e: print("wrong method chosen. available methdos: HF, DFT, and MP2")
