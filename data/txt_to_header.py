import argparse

program_desc = "Generates a valid c++ header with arrays from mass-radius data stored in .txt file\n"
program_desc += ".txt file must have a 1 line header and 2 comma seperated columns (mass,radius), units should be in solar mass and solar radii"

parser = argparse.ArgumentParser(prog="C++ Header Generator", description=program_desc)
parser.add_argument("--in_file", "-i", required=True, help="filepath to mass radius relationship txt")
parser.add_argument("--out_file", "-o", required=True, help="filepath to generated header file (../include/mass_radius.hh)")
args = parser.parse_args()

mass_array_str = "static const double white_dwarf_mass[] = {"
radi_array_str = "static const double white_dwarf_radius[] = {"
n_entries = 0

with open(args.in_file) as in_file:
    next(in_file)
    for line in in_file:
        m,r = [float(e) for e in line.strip().split(',')]
        mass_array_str += f"{m*1.989100e+33},"
        radi_array_str += f"{r*6.9599e10},"
        n_entries += 1

mass_array_str = mass_array_str[:-1] + "};"
radi_array_str = radi_array_str[:-1] + "};"

with open(args.out_file, "w+") as out_file:
    out_file.write("#pragma once\n\n")
    out_file.write(f"static const int mass_radius_length = {n_entries};\n")
    out_file.write(mass_array_str+"\n")
    out_file.write(radi_array_str+"\n")
