# OUTCAR_to_MACExyz_mod
2024/7/28, TM
This program was originally generated by Poyen Chen, The University of Tokyo.
The original program was modified to be used for large supercells.

In particular around Line 74 as below
---
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
---
was changed for more than two digit with minus lattice. I was generated by Cryspy.

