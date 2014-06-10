#!/usr/bin/env python

from gsr import *

if __name__ == '__main__':

    snap = Snapshot('/home/cmaureir/repos/gadget/out0.001/ND_M001_acc_003')
    snap.print_data_by_type(5)  # Print information of the BlackHoles (Type 5)
