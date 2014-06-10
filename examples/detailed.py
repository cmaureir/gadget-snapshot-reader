#!/usr/bin/env python

from gsr import *
import argparse as ap

class OptionParser:
    def __init__(self):
        desc = "Gadget Snapshot Reader (GSR)"
        self.parser = ap.ArgumentParser(prog="GSR",
                                        description=desc,
                                        formatter_class=ap.RawTextHelpFormatter)

        # General options
        self.general_group = self.parser.add_argument_group("General options")
        self.add_general_options()

    def add_general_options(self):
        file_help = "Input filename"
        self.general_group.add_argument("-f",
                                        dest="fname",
                                        type=str,
                                        nargs="*",
                                        help=file_help,
                                        required=True)

    def get_options(self):
        return self.parser.parse_args()



if __name__ == '__main__':
    op = OptionParser().get_options()


    # Each element of this lists will be the information of one snapshot.
    bh1_pos = []
    bh2_pos = []
    gas_pos = []


    for f in op.fname:
        print("Processing snapshot", f,"...")
        snap = Snapshot(f)
        bhs = snap.get_data_by_type(5)
        gas = snap.get_data_by_type(0)


        # Every element that we are appending will have 3 members (posx, posy, posz)
        bh1_pos.append(bhs[2][0])
        bh2_pos.append(bhs[2][1])

        # Every element that we are appending will have ALL the gas positions
        # for every snapshot
        gas_pos.append(gas[2])  # 0: Ids, 1: Masses, 2: Positions, 3: Velocities, ...

    # Now that we have the filled lists, we can iterate over them of the following
    # ways

    for element in bh1_pos:
        posx, posy, posz = element
        print("Example 1:", posx, posy, posz)

    # zero if you prefer index
    for i in range(len(bh1_pos)):
        posx, posy, posz = bh1_pos[i]
        print("Example 2:", posx, posy, posz)

    # To iterate the information of two BH, we can:
    for a, b in zip(bh1_pos, bh2_pos):
        a_posx, a_posy, a_posz = a
        b_posx, b_posy, b_posz = b
        print("Example 3a:", a_posx, a_posy, a_posz)
        print("Example 3b:", b_posx, b_posy, b_posz)

    # ...or using indexes
    for i in range(len(bh1_pos)): # Solo usamos 1, ambos tienen el mismo largo
        a_posx, a_posy, a_posz = bh1_pos[i]
        b_posx, b_posy, b_posz = bh2_pos[i]
        print("Example 3b:", a_posx, a_posy, a_posz)
        print("Example 3b:", b_posx, b_posy, b_posz)


    # For the gas particles, we have a deeper structure

    # The element will contain ALL the information of ALL the gas particles
    # for every snapshot
    for element in gas_pos:
        # The part variable will iterate over the gas particle array of one snapshot
        for part in element:
            posx, posy, posz = part
            # Particle's position inside a snapshot
            print("Gas:", posx, posy, posz)
