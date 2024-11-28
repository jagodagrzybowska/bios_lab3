from Bio import PDB
from Bio.PDB.vectors import calc_dihedral
import numpy
import matplotlib.pyplot as plt

pdb_parser = PDB.PDBParser()
structure = pdb_parser.get_structure('4ywo', '4ywo.pdb')
phi_angles = []
psi_angles = []

for model in structure:
    for chain in model:
        residues = list(chain.get_residues())
        for i in range(1, len(residues)-1):
            residue1 = residues[i-1]
            residue2 = residues[i]
            residue3 = residues[i+1]

            if ('N' in residue2 and 'CA' in residue2 and 'C' in residue2
                and 'N' in residue3 and 'C' in residue1): 
                    c1 = residue1['C'].get_coord()
                    c2 = residue2['C'].get_coord()
                    n2 = residue2['N'].get_coord()
                    n3 = residue3['N'].get_coord()
                    ca = residue2['CA'].get_coord()

                    phi = calc_dihedral(PDB.vectors.Vector(c1),
                                        PDB.vectors.Vector(n2),
                                        PDB.vectors.Vector(ca),
                                        PDB.vectors.Vector(c2))
                    psi = calc_dihedral(PDB.vectors.Vector(n2),
                                        PDB.vectors.Vector(ca),
                                        PDB.vectors.Vector(c2),
                                        PDB.vectors.Vector(n3))

                    phi_angles.append(numpy.degrees(phi))
                    psi_angles.append(numpy.degrees(psi))

plt.figure(figsize=(8, 8))
plt.grid()
plt.scatter(phi_angles, psi_angles, c='orange', edgecolors='black')
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.xticks([-180, -120, -60, 0, 60, 120, 180])
plt.yticks([-180, -120, -60, 0, 60, 120, 180])
plt.title('Wykres Ramachandrana dla struktury 4YWO')
plt.xlabel('Phi (°)')
plt.ylabel('Psi(°)')
# plt.show()
plt.savefig('wykres_ramachandran.png')
