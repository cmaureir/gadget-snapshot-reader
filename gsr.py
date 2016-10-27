#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" GADGET Snapshot Reader  """

__author__ = "Cristián Maureira and Patrick Brem"
__copyright__ = "Copyright 2014, Cristián Maureira and Patrick Brem"
__credits__ = ["Cristián Maureira", "Patrick Brem"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Cristián Maureira"
__email__ = "cmaurei@aei.mpg.de"
__status__ = "Production"

import numpy as np

from os.path import isfile
from struct import unpack, calcsize
from sys import exit, stderr

class Snapshot:
    """Class in charge of reading the binary snapshot and process all the
    information, providing also a ser of methods to access and manipulate
    the data"""

    def __init__(self, filename,\
                 enable_potential = False, \
                 enable_accelerations = False, \
                 enable_entropy_production = False, \
                 enable_timesteps = False):

        # Check and open the snapshot file
        if not isfile(filename):
            stderr.write(filename, ": No such file")
            exit(1)
        else:
            self.fname = filename
            try:
                self.binfile = open(self.fname, "rb")
            except IOError:
                stderr.write(self.fname, ": Cannot open file")

        # Extra options for the snapshot format
        self.enable_potential = enable_potential
        self.enable_accelerations = enable_accelerations
        self.enable_entropy_production = enable_entropy_production
        self.enable_timesteps = enable_timesteps

        # Reading data
        self.raw_data = self.binfile.read()
        self.data = {}
        self.byte_count = 0
        self.missing_masses = 0

        # Starting the process
        self.process_data()

    def process_data(self):

        # Snapshot elements
        #
        # 01. Header                     [256 bytes]
        # 02. Positions                  [real*4 pos(3,N)]
        # 03. Velocities                 [real*4 vel(3,N)]
        # 04. Identifications            [int*4 id(N)]
        # 05. Masses                     [real*4 masses(Nm)]
        # 06. Internal energy            [real*4 u(Ngas)]
        # 07. Density                    [real*4 rho(Ngas)]
        # 08. Smoothing lenght           [real*4 hsml(Ngas)]
        # 09. Gravitational potential    [real*4 pot(N)]        (*)
        # 10. Accelerations              [real*4 acc(3,N)]      (*)
        # 11. Rate of entropy production [real*4 dAdt(Ngas)]    (*)
        # 12. Timesteps                  [real*4 dt(N)]         (*)
        #
        # (*) Optional: Must be activated on compilation of Gadget

        # Processing header
        self.header = self.process_header()

        # Useful variables
        self.npart  = self.header['Npart']
        self.mpart  = self.header['Mpart']
        self.time   = self.header['Time']
        self.ntotal = np.sum(self.header['Npart'])
        self.ngas   = np.int(self.header['Npart'][0])


        # Processing the rest of the Snapshot
        self.data['pos']  = self.process_positions()
        self.data['vel']  = self.process_velocities()
        self.data['id']   = self.process_identifications()
        self.data['mass'] = self.process_masses()
        self.data['u']    = self.process_internal_energy()
        self.data['rho']  = self.process_density()
        self.data['hsml'] = self.process_smoothing_lenght()

        if self.enable_potential:
            self.data['pot'] = self.process_potential()
        if self.enable_accelerations:
            self.data['acc'] = self.process_accelerations()
        if self.enable_entropy_production:
            self.data['dadt'] = self.process_entropy_production()
        if self.enable_timesteps:
            self.data['dt'] = self.process_timesteps()

    #@profile
    def process_header(self):
        self.skip_empty_bytes()

        # Header format
        nbytes = 196
        fmtstring = "6i8d10i4d9i".format(nbytes)

        # Print the amount of bytes
        #print(calcsize(fmtstring))

        # Header data section
        chunk_data = self.raw_data[self.byte_count:self.byte_count + nbytes]

        # Unpacking the raw_data
        fdata = unpack(fmtstring, chunk_data)

        # Move forward to continue reading the raw_data
        self.byte_count += nbytes

        # Creating header data structure
        header = {'Npart': None,
                  'Mpart': None,
                  'Time': None,
                  'Redshift': None,
                  'FlagSfr': None,
                  'FlagFeedback': None,
                  'Nall': None,
                  'FlagCooling': None,
                  'NumFiles': None,
                  'BoxSize': None,
                  'Omega0': None,
                  'OmegaLambda': None,
                  'HubbleParam': None,
                  'FlagAge': None,
                  'FlagMetals': None,
                  'NallHW': None,
                  'flag_entr_ics': None}


        # Filling header
        header['Npart']         = np.array(fdata[0:6], dtype = np.int)
        header['Mpart']         = np.array(fdata[6:12], dtype = np.float)
        header['Time']          = np.float(fdata[12:13][0])
        header['Redshift']      = np.float(fdata[13:14][0])
        header['FlagSfr']       = np.int(fdata[14:15][0])
        header['FlagFeedback']  = np.int(fdata[15:16][0])
        header['Nall']          = np.array(fdata[16:22], dtype = np.int)
        header['FlagCooling']   = np.int(fdata[22:23][0])
        header['NumFiles']      = np.int(fdata[23:24][0])
        header['BoxSize']       = np.float(fdata[24:25][0])
        header['Omega0']        = np.float(fdata[25:26][0])
        header['OmegaLambda']   = np.float(fdata[26:27][0])
        header['HubbleParam']   = np.float(fdata[27:28][0])
        header['FlagAge']       = np.int(fdata[28:29][0])
        header['FlagMetals']    = np.int(fdata[29:30][0])
        header['NallHW']        = np.array(fdata[30:36], dtype=np.int)
        header['flag_entr_ics'] = np.int(fdata[36:37][0])

        # Advancing bytes to reach the 256 bytes of the header
        # 256 - 196 = 60
        self.byte_count += (256 - nbytes)

        self.skip_empty_bytes()

        return header

    #@profile
    def unpack_data_section(self, elements, format="f"):

        # String format
        if format is "f":
            fmtstring = "{0:d}f".format(elements)
        elif format is "i":
            fmtstring = "{0:d}i".format(elements)
        else:
            stderr.write("Invalid unpack format")
            exit(1)

        # Amount of bytes
        #print(calcsize(fmtstring))

        # Data section:
        # The amount of elements are `float` or `int`,
        # both data types are 4 bytes long.
        chunk_data = self.raw_data[self.byte_count:self.byte_count + elements * 4]

        # Unpacking the raw_data
        everything = unpack(fmtstring, chunk_data)

        return everything

    #@profile
    def get_chunk_of_data(self, dim, data_type, everything, ptypes=range(6)):

        if dim == 1:
            #data = [np.zeros(i, dtype=data_type) for i in self.npart]
            data = [np.zeros(i, dtype=data_type) for i in self.npart]
        elif dim == 3:
            #data = [np.zeros((i, 3)) for i in self.npart]
            data = [np.zeros((i, 3)) for i in self.npart]
        else:
            print("[Error] Incorrect dimensions reading the snapshot")
            exit(1)

        offset = 0
        for i in ptypes:
            npart = self.npart[i]
            if npart:
                ini = offset * dim
                end = ini + npart * dim
                chunk = everything[ini:end]
                if dim == 1:
                    data[i] = np.array(chunk, dtype=data_type)
                elif dim == 3:
                    data[i] = np.reshape(chunk, (npart, 3))
                offset += npart

        return data

    def skip_empty_bytes(self):
        # Empty space between every snapshot element
        self.byte_count += 4

    def print_header(self):
        print("{")
        for i in self.header:
            print("\t",i,":", self.header[i])
        print("}")

    def process_positions(self):
        self.skip_empty_bytes()

        # Amount of elements to read
        elements = self.ntotal * 3
        # Getting the data
        everything = self.unpack_data_section(elements)
        pos = self.get_chunk_of_data(3, np.float, everything)

        # Move forward to continue reading the raw_data
        self.byte_count += elements*4

        self.skip_empty_bytes()

        return pos

    def process_velocities(self):
        self.skip_empty_bytes()

        # Amount of elements to read
        elements = self.ntotal * 3
        # Getting the data
        everything = self.unpack_data_section(elements)
        vel = self.get_chunk_of_data(3, np.float, everything)
        # Move forward to continue reading the raw_data
        self.byte_count += elements * 4

        self.skip_empty_bytes()

        return vel


    def process_identifications(self):
        self.skip_empty_bytes()

        # Amount of particles
        elements = self.ntotal
        # Getting the data
        everything = self.unpack_data_section(elements, "i")
        id = self.get_chunk_of_data(1, np.int, everything)

        self.byte_count += elements*4

        self.skip_empty_bytes()

        return id

    def process_masses(self):
        #mass = [np.zeros(i) for i in self.npart]
        mass = [np.zeros(i) for i in self.npart]
        missing_types = []

        for i in range(6):
            parts = self.npart[i]
            m  = self.mpart[i]

            if parts > 0:
                # Check how many particles have a different mass as an array
                # instead of using a value from the header
                if not m:
                    self.missing_masses += parts
                    missing_types.append(i)
                else:
                    mass[i].fill(m)

        if self.missing_masses > 0:
            self.skip_empty_bytes()

            elements = self.missing_masses
            everything = self.unpack_data_section(elements)

            offset = 0
            for i in missing_types:
                nparts = self.npart[i]
                chunk = everything[offset:offset+nparts]
                mass[i] = np.array(chunk)
                offset += nparts
            self.byte_count += elements * 4

            self.skip_empty_bytes()

        return mass

    def process_internal_energy(self):
        self.skip_empty_bytes()

        elements = np.int(self.ngas)
        everything = self.unpack_data_section(elements)
        u = np.array(everything)
        self.byte_count += elements * 4

        self.skip_empty_bytes()

        return u

    def process_density(self):

        rho = np.zeros(self.ngas)
        self.skip_empty_bytes()

        # Just Gas particles
        elements = np.int(self.ngas)
        everything = self.unpack_data_section(elements)
        rho = np.array(everything)
        self.byte_count += elements * 4

        self.skip_empty_bytes()

        return rho

    def process_smoothing_lenght(self):

        hsml = np.zeros(self.ngas)

        self.skip_empty_bytes()

        elements = np.int(self.ngas)
        everything = self.unpack_data_section(elements)
        hsml = np.array(everything)
        self.byte_count += elements * 4

        self.skip_empty_bytes()

        return hsml


    def process_potential(self):

        pot = [np.zeros(i, dtype=np.float) for i in self.npart]

        self.skip_empty_bytes()

        elements = np.sum(self.npart)
        everything = self.unpack_data_section(elements)
        pot = self.get_chunk_of_data(1, np.float, everything)

        self.skip_empty_bytes()

        return pot


    def process_accelerations(self):
        #acc = [np.zeros((i, 3)) for i in self.npart]
        acc = [np.zeros((i, 3)) for i in self.npart]

        self.skip_empty_bytes()

        # Amount of elements to read
        elements = np.ntotal * 3
        everything = self.unpack_data_section(elements)
        acc = self.get_chunk_of_data(3, np.float, everything)
        # Move forward to continue reading the raw_data
        self.byte_count += elements*4
        self.skip_empty_bytes()

        return acc

    def process_entropy_production(self):

        #dadt = np.zeros(self.ngas)
        dadt = np.zeros(self.ngas)

        self.skip_empty_bytes()

        elements = self.ngas
        everything = self.unpack_data_section(elements)
        dadt = np.array(everything)
        self.byte_count += elements * 4

        self.skip_empty_bytes()

        return dadt

    def process_timesteps(self):

        dt= [np.zeros(i, dtype=np.float) for i in self.npart]
        self.skip_empty_bytes()

        # Amount of elements to read
        assert(np.sum(self.npart) == self.ntotal)
        elements = self.ntotal
        # Getting the data
        everything = self.unpack_data_section(elements)
        dt = self.get_chunk_of_data(1, np.float, everything)

        # Move forward to continue reading the raw_data
        self.byte_count += elements*4

        self.skip_empty_bytes()

        return dt

    # Utilities

    def get_data_by_type(self, ptype):
        d = {}
        if ptype < 0 or ptype > 5:
            print(pytpe, "Invalid data type, use only 0, 1, 2, 3, 4 or 5")
            return d

        d['pos']  = self.data['pos'][ptype]
        d['vel']  = self.data['vel'][ptype]
        d['id']   = self.data['id'][ptype]
        d['mass'] = self.data['mass'][ptype]

        if not ptype:
            d['u']   = self.data['u']
            d['rho'] = self.data['rho']
            d['hsml'] = self.data['hsml']

        d['dt'] = self.data['dt'][ptype]

        return d

    def to_ascii(self):
        def get_tuple(key):
            if len(key) > 1:
                return tuple(i for i in self.data[key])

        ids  = np.concatenate(get_tuple('id'),   axis = 0)
        mass = np.concatenate(get_tuple('mass'), axis = 0)
        pos  = np.concatenate(get_tuple('pos'),  axis = 0)
        vel  = np.concatenate(get_tuple('vel'),  axis = 0)

        u    = self.data['u']
        rho  = self.data['rho']
        hsml = self.data['hsml']

        u.resize(self.ntotal, refcheck=False)
        rho.resize(self.ntotal, refcheck=False)
        hsml.resize(self.ntotal, refcheck=False)

        fmtstring = ['%8d', '% .5e',
                     '% .5e', '% .5e', '% .5e',
                     '% .5e', '% .5e', '% .5e',
                     '% .5e', '% .5e', '% .5e']

        assert(len(id) == len(mass) == len(pos) == len(vel))
        assert(len(id) == len(u) == len(rho) == len(hsml))

        data_to_print = np.c_[id, mass, pos, vel, u, rho, hsml]
        np.savetxt(self.fname+'.asc', data_to_print, fmt = fmtstring)

    def print_data_by_type(self, ptype):
        if ptype < 0 or ptype > 5:
            print(pytpe, "Invalid data type, use only 0, 1, 2, 3, 4 or 5")
            return None
        for i in range(self.npart[ptype]):
            pid = self.data['id'][ptype][i]
            m = self.data['mass'][ptype][i]
            rx, ry, rz = self.data['pos'][ptype][i]
            vx, vy, vz = self.data['vel'][ptype][i]

            if ptype == 0:
                u = self.data['u'][i]
                rho = self.data['rho'][i]
                hsml = self.data['hsml'][i]
            else:
                u = 0.0
                rho = 0.0
                hsml = 0.0

            fmts = '%8d %1.5e '
            fmts += '% .5e % .5e % .5e '
            fmts += '% .5e % .5e % .5e '
            fmts += '% .5e % .5e % .5e'
            print(fmts % (pid, m, rx, ry, rz, vx, vy, vz, u, rho, hsml))
