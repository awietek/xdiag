#!/usr/bin/env python
import sys
import numpy as np
import toml

assert(len(sys.argv)==2)
filename = sys.argv[1]

header = ""
with open(filename, 'r') as fl:
    for line in fl.readlines():
        if line[0] == "#":
            clean_line = line[1:].strip()
            if (len(clean_line)>0):
                header += "# " + line[1:].strip() + "\n"
print(header)

coordinate_str = "Coordinates = [\n"
with open(filename, 'r') as fl:
    is_coordinate_line=False
    cline = 0
    for line in fl.readlines():
        if "[Interactions]" in line:
            coordinate_str = coordinate_str[:-2] + "\n]\n"
            break

        if is_coordinate_line:
            spt = line.split(" ")
            if cline == 0:
                coordinate_str += "  " + str([float(spt[0]), float(spt[1])]) + ",\n"
            else:
                coordinate_str += "  " + str([float(spt[0]), float(spt[1])]) + ",\n"                
            cline += 1
            
            
        if "[Sites]" in line:
            is_coordinate_line = True

print(coordinate_str)


interactions_str = "Interactions = [\n"
with open(filename, 'r') as fl:
    is_interaction_line=False
    iline = 0
    for line in fl.readlines():
        if "[SymmetryOps]" in line:
            interactions_str = interactions_str[:-2] + "\n]\n"
            break

        if is_interaction_line:
            spt = line.split(" ")
            bond = [str(spt[0]), str(spt[1])] + list(map(int, spt[2:]))
            if iline == 0:
                interactions_str += "  " + str(bond) + ",\n"
            else:
                interactions_str += "  " + str(bond) + ",\n"                
            iline += 1
            
        if "[Interactions]" in line:
            is_interaction_line = True


print(interactions_str)

symmetries_str = "Symmetries = [\n"
with open(filename, 'r') as fl:
    is_symmetry_line=False
    sline = 0
    for line in fl.readlines():
        if "[Irreps]" in line:
            symmetries_str = symmetries_str[:-2] + "\n]\n"
            break

        if is_symmetry_line:
            spt = line.split(" ")
            symmetries = list(map(int, spt[1:]))
            if sline == 0:
                symmetries_str += "  " + str(symmetries) + ",\n"
            else:
                symmetries_str += "  " + str(symmetries) + ",\n"                
            sline += 1
            
        if "[SymmetryOps]" in line:
            is_symmetry_line = True

print(symmetries_str)

irreps_str = "# Irreducible representations\n"
def irrep_str(name, allowed_ops, characters, momentum):
    s = "[" + name + "]\n"
    charstr = "[\n"
    for chi in characters:
        charstr += "  [{:.16f}, {:.16f}],\n".format(np.real(chi), np.imag(chi))
    charstr = charstr[:-2] + "\n]"
        
    s += "characters = " + charstr + "\n"
    s += "allowed_symmetries = " + str(allowed_ops) + "\n"
    s += "momentum = [{:.16f}, {:.16f}]".format(momentum[0], momentum[1]) 
    s += "\n\n"
    return s

with open(filename, 'r') as fl:
    is_irrep_line=False
    iline = 0
    inum = 0
    for line in fl.readlines():
        if "[Representation]" in line:
            if inum != 0:
                irreps_str += irrep_str(name, allowed_ops, characters, momentum)
            iline = 0
            is_irrep_line = True
            name = line.split("=")[1].strip()
            characters = []
            inum += 1
        if is_irrep_line:
            if iline == 1:
                k_spt = line.replace("[K=(", "").replace(")]", "").split(",")
                momentum = list(map(float, k_spt))
            if iline== 3:
                allowed_ops = list(map(int, line.split(" ")))
            if iline > 3:
                spt = line.split(" ")
                characters.append(complex(float(spt[0]), float(spt[1])))

            iline += 1

irreps_str += irrep_str(name, allowed_ops, characters, momentum)

print(irreps_str)

