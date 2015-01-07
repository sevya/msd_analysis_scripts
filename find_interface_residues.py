#!/usr/bin/env python2.7

from Bio.PDB import *
from Bio.PDB import Polypeptide
import warnings
import sys
from Bio.PDB import PDBExceptions
import math
from optparse import OptionParser

def distance_cutoff(structure_res, group1_chains, group2_chains, nearby_atom_cutoff, CB_dist_cutoff):
	interface_residues = set()
	for index, residue1 in enumerate(structure_res):
		for residue2 in structure_res[index:]:
			if residue1 == residue2:
				continue
			elif (residue1.get_parent().get_id() in group1_chains and residue2.get_parent().get_id() in group1_chains) or (residue1.get_parent().get_id() in group2_chains and residue2.get_parent().get_id() in group2_chains):
				continue
			else: 
				#these are interesting residues so let's assign atoms
				try:
					residue1_atom = residue1["CB"]
				except KeyError:
					residue1_atom = residue1["CA"]
				try:
					residue2_atom = residue2["CB"]
				except KeyError:
					residue2_atom = residue2["CA"]
				if residue1_atom - residue2_atom < float(CB_dist_cutoff):
					#check the distance of any atom to any atom
					for atom1 in residue1:
						if (atom1.get_id()[0] == "H") or (len(atom1.get_id()) > 1) and (atom1.get_id()[1] == "H"):
							continue
						for atom2 in residue2:
							if  (atom2.get_id()[0] == "H") or (len(atom2.get_id()) > 1) and (atom2.get_id()[1] == "H"):
								continue
							if atom1 - atom2 < float(nearby_atom_cutoff):
								interface_residues.add(residue1)
								interface_residues.add(residue2)
								continue
	return interface_residues

def calculate_vector(residue1, residue2, CB_dist_cutoff, vector_dist_cutoff, vector_angle_cutoff):
	try:
		if residue1["CB"] - residue2["CB"] > float(CB_dist_cutoff):
			return False
		if residue1["CB"] - residue2["CB"] > float(vector_dist_cutoff):
			return False
	except KeyError:
		try:
			if residue1["CB"] - residue2["CA"] > float(CB_dist_cutoff):
				return False
			if residue1["CB"] - residue2["CA"] > float(vector_dist_cutoff):
				return False
		except KeyError:
			return False

	ca1_cb1 = residue1["CA"] - residue1["CB"] # P12
	cb1_ca2 = residue1["CB"] - residue2["CA"] # P23
	ca1_ca2 = residue1["CA"] - residue2["CA"] # P13
	P12 = ca1_cb1
	P23 = cb1_ca2
	P13 = ca1_ca2
	#click cb1,ca1,cb2 in pymol angle
	angle = math.degrees(math.acos( ( (P12*P12) + (P13*P13) - (P23*P23) ) / (2 * P12 * P13) ) )

	if angle <= float(vector_angle_cutoff):
		return True

	return False

def vector_cutoff(structure_res, interface_residues, CB_dist_cutoff, vector_dist_cutoff, vector_angle_cutoff, group1_chains, group2_chains):
	for index, residue1 in enumerate(structure_res):
		for residue2 in structure_res[index:]:
			if residue1 == residue2:
				continue
			elif (residue1.get_parent().get_id() in group1_chains and residue2.get_parent().get_id() in group1_chains) or (residue1.get_parent().get_id() in group2_chains and residue2.get_parent().get_id() in group2_chains): # only residues on opposite sides of interface
				continue
						
			if calculate_vector(residue1, residue2, CB_dist_cutoff, vector_dist_cutoff, vector_angle_cutoff):
				interface_residues.add(residue1)
			if calculate_vector(residue2, residue1, CB_dist_cutoff, vector_dist_cutoff, vector_angle_cutoff):
				interface_residues.add(residue2)
	return interface_residues

