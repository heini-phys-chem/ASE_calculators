#!/usr/bin/env python
#
# Author: Garnet Chan <gkc1000@gmail.com>
# adapted by Stefan Heinen:
# - added forces
# - changed the ase_atoms_to_pyscf function
# - added an mp2 wrapper

'''
ASE package interface
'''

import numpy as np
from ase.calculators.calculator import Calculator
import ase.dft.kpoints
from ase.lattice import bulk

from pyscf import gto, scf, grad, mp2

def get_mp2_energy():
	test = 1

def ase_atoms_to_pyscf(ase_atoms):
		'''Convert ASE atoms to PySCF atom.

		Note: ASE atoms always use A.
		'''
#return [[ase_atoms.get_chemical_symbols(), ase_atoms.get_positions()] for i, atom in enumerate(ase_atoms)]
		return [ [ase_atoms.get_chemical_symbols()[i], ase_atoms.get_positions()[i]] for i in range(len(ase_atoms.get_positions()))]

atoms_from_ase = ase_atoms_to_pyscf

class PySCF(Calculator):
		implemented_properties = ['energy', 'forces']

		def __init__(self, restart=None, ignore_bad_restart_file=False,
								 label='PySCF', atoms=None, scratch=None, **kwargs):
				"""Construct PySCF-calculator object.

				Parameters
				==========
				label: str
						Prefix to use for filenames (label.in, label.txt, ...).
						Default is 'PySCF'.

				mfclass: PySCF mean-field class
				molcell: PySCF :Mole: or :Cell:
				"""
				Calculator.__init__(self, restart=None, ignore_bad_restart_file=False,
														label='PySCF', atoms=None, scratch=None, **kwargs)

				# TODO
				# This explicitly refers to "cell". How to refer
				# to both cell and mol together?

				self.mf=None
				self.initialize(**kwargs)

		def initialize(self, molcell, mf_class, mf_dict):
#				if not molcell.unit.startswith(('A','a')):
#						raise RuntimeError("PySCF unit must be A to work with ASE")

				self.molcell=molcell
				self.mf_class=mf_class
				self.mf_dict=mf_dict

		def set(self, **kwargs):
				changed_parameters = Calculator.set(self, **kwargs)
				if changed_parameters:
						self.reset()

		def calculate(self, atoms=None, properties=['energy', 'forces'],
									system_changes=['positions', 'numbers', 'cell',
																	'pbc', 'charges','magmoms']):

				Calculator.calculate(self, atoms)

				calc_molcell = self.molcell.copy()
				calc_molcell.atom = ase_atoms_to_pyscf(atoms)
#				calc_molcell.a = atoms.cell
				calc_molcell.build(None,None)
				self.mf = self.mf_class(calc_molcell)
				for key in self.mf_dict:
						self.mf.__dict__[key] = self.mf_dict[key]

				self.results['energy']=self.mf.scf(verbose=0)
				self.results['forces']=-1*grad.RHF(self.mf).kernel() # convert forces to gradient (*-1) !!!!! for the NEB run
				self.results['mf']=self.mf

def make_kpts(cell, nks):
		raise DeprecationWarning('Use cell.make_kpts(nks) instead.')
		return kpts

