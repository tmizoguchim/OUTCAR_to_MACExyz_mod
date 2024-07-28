import re
import numpy as np
import random
import argparse

def output_data_extract(input_data, interval=1, only_final=True):
    with open(input_data, "r") as f:
        p_found = False
        l_found = False
        step_counter = 0
        start = False
        energies = []
        all_atomic_positions = []
        all_forces = []
        all_atom = []
        all_lattice = []
        
        for line in f:
            if "POSCAR = " in line:
                s_line = line.replace("POSCAR = ", "") 
                elements = s_line.split()
                atoms_list = []
                for element in elements:
                    if element.startswith("-"):
                        print(f"Warning: Skipping element starting with '-': {element}")
                        continue  # Skip elements starting with "-"
                    match = re.match(r"([A-Za-z]+)(\d+)", element)
                    if match:
                        atom_type = match.group(1)
                        count = int(match.group(2))
                        atoms_list.extend([atom_type] * count)
                    else:
                        print(f"Warning: Skipping invalid element format: {element}")
                        
            if "aborting loop because EDIFF is reached" in line:
                start = True
                step_counter += 1
                if step_counter % interval == 0 or step_counter == 1:
                    all_atom.append(atoms_list)
                
            if "Ionic step " in line:
                start = False
                
            if start:
                if "free  energy" in line:
                    if step_counter % interval == 0 or step_counter == 1:
                        e_line = " ".join(line.split()).split()
                        energy = float(e_line[-2])
                        if energy > 0:
                            print(f"Warning: Skipping positive energy value: {energy}")
                            continue  # Skip positive energy values
                        energies.append(energy)

                if "POSITION" in line and "TOTAL-FORCE" in line:
                    p_found = True
                    l_p = 0
                    atomic_position = []
                    forces = []
                if p_found:
                    l_p += 1
                    if "------------" in line:
                        continue
                    if "total drift" in line:
                        p_found = False
                        if step_counter % interval == 0 or step_counter == 1:
                            all_atomic_positions.append(atomic_position)
                            all_forces.append(forces)
                        continue
                    if l_p >= 3:
                        p_line = " ".join(line.split()).split()
                        atomic_position.append(p_line[:3])
                        forces.append(p_line[3:])

                if "direct lattice vectors" in line:
                    l_found = True
                    l_l = 0
                    lattice = []
                if l_found:
                    l_l += 1
                    if line.isspace():
                        l_found = False
                        if step_counter % interval == 0 or step_counter == 1:
                            all_lattice.append(lattice)
                        continue
                    if l_l >= 2:
                        # Insert a space before any negative sign that follows a number without a space
                        line = re.sub(r'(\d)-', r'\1 -', line)
                        l_line = " ".join(line.split()).split()
                        lattice.extend(l_line[:3])
                           
    if not only_final:
        return all_atom, all_atomic_positions, all_forces, all_lattice, energies
    else:
        return [atoms_list], [atomic_position], [forces], [lattice], [energy]
    

def multi_extract(input_data_list, interval=1, only_final=False):
    all_atom = []
    all_atomic_positions = []
    all_forces = []
    all_lattice = []
    energies = []
    for data in input_data_list:
        try:
            atom, atomic_position, forces, lattice, energy = output_data_extract(data, interval=interval, only_final=only_final)
            all_atom.extend(atom)
            all_atomic_positions.extend(atomic_position)
            all_forces.extend(forces)
            all_lattice.extend(lattice)
            energies.extend(energy)
        except Exception as e:
            print(f"Error processing file {data}: {e}")
    return all_atom, all_atomic_positions, all_forces, all_lattice, energies


def convert_data_to_xyz_format(input_data_list, interval=1, only_final=False, train_name="train.xyz", test_name="test.xyz", train_test_split=False, test_ratio=0.1, random_seed=42):
    all_atom, all_atomic_positions, all_forces, all_lattice, energies = multi_extract(input_data_list, interval=interval, only_final=only_final)
    if not train_test_split:
        with open(train_name, "w") as train_xyz:
            for n in range(len(energies)):
                if len(all_atom[n]) == 0:
                    continue  # Skip entries with zero atoms
                train_xyz.write(str(len(all_atom[n])) + "\n")
                train_xyz.write("Lattice=" + '"' + " ".join(all_lattice[n]) + '" ')
                train_xyz.write("Properties=species:S:1:pos:R:3:forces:R:3 ")
                train_xyz.write("energy=" + '"' + str(energies[n]) + '"' + "\n")
                atom = all_atom[n]
                p = all_atomic_positions[n]
                f = all_forces[n]
                for a in range(len(atom)):
                    train_xyz.write(str(atom[a]) + "    " + "    ".join(p[a]) + "    " + "    ".join(f[a]) + "\n")
    else:
        idx = list(range(len(energies)))
        random.seed(random_seed)
        random.shuffle(idx)
        test_size = int(len(idx) * test_ratio)
        test_idx = idx[:test_size]
        train_idx = idx[test_size:]
        print(f"The total number of data is {len(idx)}")
        print(f"{len(train_idx)} data will be utilized in {train_name}, and {len(test_idx)} data will be utilized in {test_name}")

        with open(train_name, "w") as train_xyz, open(test_name, "w") as test_xyz:
            for n in range(len(energies)):
                if len(all_atom[n]) == 0:
                    continue  # Skip entries with zero atoms
                if n in test_idx:
                    test_xyz.write(str(len(all_atom[n])) + "\n")
                    test_xyz.write("Lattice=" + '"' + " ".join(all_lattice[n]) + '" ')
                    test_xyz.write("Properties=species:S:1:pos:R:3:forces:R:3 ")
                    test_xyz.write("energy=" + '"' + str(energies[n]) + '"' + "\n")
                    atom = all_atom[n]
                    p = all_atomic_positions[n]
                    f = all_forces[n]
                    for a in range(len(atom)):
                        test_xyz.write(str(atom[a]) + "    " + "    ".join(p[a]) + "    " + "    ".join(f[a]) + "\n")
                
                if n in train_idx:
                    train_xyz.write(str(len(all_atom[n])) + "\n")
                    train_xyz.write("Lattice=" + '"' + " ".join(all_lattice[n]) + '" ')
                    train_xyz.write("Properties=species:S:1:pos:R:3:forces:R:3 ")
                    train_xyz.write("energy=" + '"' + str(energies[n]) + '"' + "\n")
                    atom = all_atom[n]
                    p = all_atomic_positions[n]
                    f = all_forces[n]
                    for a in range(len(atom)):
                        train_xyz.write(str(atom[a]) + "    " + "    ".join(p[a]) + "    " + "    ".join(f[a]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract and convert data to XYZ format for ML training.")
    parser.add_argument("input_data_list", nargs='+', help="List of input data files.")
    parser.add_argument("--interval", type=int, default=1, help="Interval for extracting data.")
    parser.add_argument("--only_final", action='store_true', help="Extract only final step data.")
    parser.add_argument("--train_name", type=str, default="train.xyz", help="Name of the output training XYZ file.")
    parser.add_argument("--test_name", type=str, default="test.xyz", help="Name of the output testing XYZ file.")
    parser.add_argument("--train_test_split", action='store_true', help="Whether to split the data into training and testing sets.")
    parser.add_argument("--test_ratio", type=float, default=0.1, help="Ratio of data to be used for testing.")
    parser.add_argument("--random_seed", type=int, default=42, help="Random seed for data shuffling.")
    
    args = parser.parse_args()
    
    convert_data_to_xyz_format(
        args.input_data_list,
        interval=args.interval,
        only_final=args.only_final,
        train_name=args.train_name,
        test_name=args.test_name,
        train_test_split=args.train_test_split,
        test_ratio=args.test_ratio,
        random_seed=args.random_seed
    )
