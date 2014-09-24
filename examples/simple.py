#!/usr/bin/env python

from gsr import *
import sys

if __name__ == '__main__':

    # Processing the Snapshot
    snap = Snapshot(sys.argv[1])

    # Printing header
    snap.print_header()

    # Getting all the information of the type-5 particle
    bhs = snap.get_data_by_type(5)

    # Getting type-1 particle mass
    gas_mass = snap.get_data_by_type(1)['mass'][0]

    print(bhs['pos'])
    print(gas_mass)

