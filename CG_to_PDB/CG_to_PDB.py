# -*- coding: utf-8 -*-
import argparse
from Bio import PDB
import os

argParser = argparse.ArgumentParser(description="Rebuild RNA structure from coarse-grain representation")
argParser.add_argument("-i", "--input", help="Input coarse-grain PDB file", required=True)
argParser.add_argument("-o", "--output", help="Output full-atom PDB file", required=True)
argParser.add_argument("-t", "--templates", help="Path to the directory containing template PDB files", required=True)

args = argParser.parse_args()

pdb_parser = PDB.PDBParser()

template_residues = {
    "A": pdb_parser.get_structure("A", os.path.join(args.templates, "A_template.pdb"))[0]["A"][26],
    "U": pdb_parser.get_structure("U", os.path.join(args.templates, "U_template.pdb"))[0]["A"][11],
    "G": pdb_parser.get_structure("G", os.path.join(args.templates, "G_template.pdb"))[0]["A"][24],
    "C": pdb_parser.get_structure("C", os.path.join(args.templates, "C_template.pdb"))[0]["A"][29],
    "R": pdb_parser.get_structure("R", os.path.join(args.templates, "ribose_template.pdb"))[0]["A"][1],  # Ryboza
}

structureId = os.path.splitext(os.path.basename(args.input))[0]
coarseGrainStructure = pdb_parser.get_structure(structureId, args.input)

fullAtomStructure = PDB.Structure.Structure(structureId + "_cg")
sup = PDB.Superimposer()

for model in coarseGrainStructure:
    newModel = PDB.Model.Model(model.id)
    fullAtomStructure.add(newModel)
    for chain in model:
        newChain = PDB.Chain.Chain(chain.id)
        newModel.add(newChain)
        for i, res in enumerate(chain):
            resname = res.get_resname()
            if resname not in ["A", "C", "G", "U", "R"]:
                continue

            fullAtomResidue = template_residues[resname].copy()
            fullAtomResidue.id = (" ", i + 1, " ")

            if resname in ["A", "G"]:  
                fixedAtomList = [res[atom] for atom in ["N9", "C2", "C6"]]
                movingAtomList = [fullAtomResidue[atom] for atom in ["N9", "C2", "C6"]]
            elif resname in ["C", "U"]: 
                fixedAtomList = [res[atom] for atom in ["N1", "C2", "C4"]]
                movingAtomList = [fullAtomResidue[atom] for atom in ["N1", "C2", "C4"]]
            elif resname == "R":  
                fixedAtomList = [res[atom] for atom in ["P", "C4'"]]
                movingAtomList = [fullAtomResidue[atom] for atom in ["P", "C4'"]]

            sup.set_atoms(fixedAtomList, movingAtomList)

            for atom in fullAtomResidue:
                atom.transform(sup.rotran[0], sup.rotran[1])

            newChain.add(fullAtomResidue.copy())

io = PDB.PDBIO()
io.set_structure(fullAtomStructure)
io.save(args.output)
