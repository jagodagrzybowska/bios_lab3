from Bio import PDB
import numpy
import matplotlib.pyplot as plt

filename = input('Podaj nazwę pliku PDB: ')
try:
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure('structure', filename)
    coords = []

    for model in structure:
        for chain in model:
            for res in chain:
                for atom in res.get_atoms():
                    if atom.get_name() == 'CA':
                        coords.append(atom.get_coord())

    atom_count = len(coords)
    contact_map = numpy.zeros((atom_count,atom_count))
    threshold = 8.0

    for i in range(atom_count):
        for j in range(atom_count):
            if i!=j:
                distance = numpy.linalg.norm(coords[i]-coords[j])
                if distance<threshold:
                    contact_map[i,j]=1

    plt.figure(figsize=(8, 8))
    plt.imshow(contact_map, cmap='Greys') 
    plt.title('Mapa kontaktow')
    plt.xlabel('Numer atomu')
    plt.ylabel('Numer atomu')
    # plt.show()
    plt.savefig('mapa.png')

except FileNotFoundError:
    print(f'Nie znaleziono pliku {filename}. Upewnij sie ze nazwa jest poprawna i plik znajduje sie w tym katalogu')
