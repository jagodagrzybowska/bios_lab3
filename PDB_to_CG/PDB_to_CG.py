# -*- coding: utf-8 -*-

from Bio import PDB
import argparse
import os

def filter_atoms_to_cg(input_pdb, output_pdb):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", input_pdb)
    
    cg_structure = PDB.Structure.Structure("RNA_CG")
 
    for model in structure:
        new_model = PDB.Model.Model(model.id)
        cg_structure.add(new_model)
        for chain in model:
            new_chain = PDB.Chain.Chain(chain.id)
            new_model.add(new_chain)
            for i, residue in enumerate(chain):
         
                if residue.get_resname() not in {"A", "C", "G", "U"}:
                    continue
                
                new_residue = PDB.Residue.Residue((" ", i + 1, " "), residue.get_resname(), " ")
                new_chain.add(new_residue)

                if "P" in residue:
                    new_residue.add(residue["P"])
                new_residue.add(residue["C4'"])
                
                if residue.get_resname() in {"U", "C"}:  
                    for atom_name in ["N1", "C2", "C4"]:
                        if atom_name in residue:
                            new_residue.add(residue[atom_name])
                else:  
                    for atom_name in ["N9", "C2", "C6"]:
                        if atom_name in residue:
                            new_residue.add(residue[atom_name])
    

    io = PDB.PDBIO()
    io.set_structure(cg_structure)
    io.save(output_pdb)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Konwersja PDB do reprezentacji gruboziarnistej")
    parser.add_argument("--input", required=True, help="Plik wejsciowy PDB")
    parser.add_argument("--output", required=True, help="Plik wyjsciowy PDB")
    args = parser.parse_args()
    
    filter_atoms_to_cg(args.input, args.output)
