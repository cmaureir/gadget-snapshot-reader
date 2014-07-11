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
from struct import unpack
from sys import exit, stderr

class Snapshot:
    """Class in charge of the read-process of every snapshot"""

    def __init__(self, filename):
        if not isfile(filename):
            stderr.write(filename, ": No such file")
            exit(1)
        else:
            self.fname = filename
            try:
                self.binfile = open(self.fname, "rb")
            except IOError:
                stderr.write(self.fname, ": Cannot open file")

            # Information of the Snapshot
            self.sdata = {}

            # Process header
            self.sdata['header'] = self.ProcessHeader()

            # Process positions float[N][3]
            self.sdata['pos'] = self.ProcessParticles(self.unpackPositions)

            # Process velocities float[N][3]
            self.sdata['vel'] = self.ProcessParticles(self.unpackVelocities)

            # Process ID int[N]
            self.sdata['id'] = self.ProcessParticles(self.unpackID)

            # Process Masses float[Nm]
            # If the information is in the header, the snapshot
            # will not contain any array with the masses, so that is
            # why we need to skip the process
            if self.check_empty_masses():
                self.sdata['mass'] = self.ProcessParticles(self.unpackMasses)
            else:
                self.sdata['mass'] = self.m

            # Process Internal Energy float[Ngas]
            self.sdata['u'] = self.ProcessParticles(self.unpackEnergy)

            # Process density float[Ngas]
            self.sdata['rho'] = self.ProcessParticles(self.unpackRho)

            # TO DO
            # Process smoothing length float[Ngas]
            # Process gravitational potential float[N]
            # Process accelerations float[N][3]
            # Process Rate of entropy production float[Ngas]
            # Process timesteps of particles float[N]

    def __exit__(self):
        self.binfile.close()

    def check_empty_masses(self):
        self.m = [np.zeros(i) for i in self.Npart]
        self.missing_masses = 0
        check = False

        for i in range(6):
            total = self.Npart[i]
            mass = self.mpart[i]

            if total > 0:
                if mass == 0:
                    self.missing_masses += total
                    check = True
                elif mass > 0:
                    self.m[i].fill(mass)
        return check

    def getRecordLength(self, instring):
        "Takes a string of some length and interprets it as a series of bits"
        return unpack('i', instring)[0]

    # Header processing
    def ProcessHeader(self):
        # because at the end another field is reserved for the length again
        self.headerlength = self.getRecordLength(self.binfile.read(4)) + 4
        header = self.binfile.read(self.headerlength)
        return self.unpackHeader(header)

    def unpackHeader(self, instring):
        fmtstring = "6i8d9i{0:d}x".format(self.headerlength-124)
        everything = unpack(fmtstring, instring)
        # list of 6 items giving number of each particle
        self.Npart = np.array(everything[:6])
        self.mpart = np.array(everything[6:12])
        self.time = everything[12]
        self.Ntot = self.Npart.sum()
        self.Ngas = self.Npart[0]
        return {'Npart': self.Npart,
                'Mpart': self.mpart,
                'Time': self.time,
                'Ntot': self.Ntot}

    def ProcessParticles(self, unpack):
        nbytes = self.getRecordLength(self.binfile.read(4)) + 4
        body = self.binfile.read(nbytes)
        return unpack(body)

    # Positions processing
    def unpackPositions(self, instring):
        fmtstring = "{0:d}f4x".format(self.Ntot*3)
        everything = unpack(fmtstring, instring)
        self.pos = [np.zeros((i, 3)) for i in self.Npart]

        offset = 0
        for i in range(6):
            chunk = everything[offset*3:offset*3+self.Npart[i]*3]
            self.pos[i] = np.reshape(chunk, (self.Npart[i], 3))
            offset += self.Npart[i]
        return self.pos

    # Velocities processing
    def unpackVelocities(self, instring):
        fmtstring = "{0:d}f4x".format(self.Ntot*3)
        everything = unpack(fmtstring, instring)

        self.vel = [np.zeros((i, 3)) for i in self.Npart]

        offset = 0
        for i in range(6):
            chunk = everything[offset*3:offset*3 + self.Npart[i]*3]
            self.vel[i] = np.reshape(chunk, (self.Npart[i], 3))
            offset += self.Npart[i]
        return self.vel

    # Id processing
    def unpackID(self, instring):
        fmtstring = "{0:d}i4x".format(self.Ntot)
        everything = unpack(fmtstring, instring)

        self.ID = [np.zeros(i, dtype=np.int) for i in self.Npart]

        offset = 0
        for i in range(6):
            chunk = everything[offset:offset+self.Npart[i]]
            self.ID[i] = np.array(chunk, dtype=np.int)
            offset += self.Npart[i]
        return self.ID

    # Mass processing
    def unpackMasses(self, instring):
        if self.missing_masses > 0:
            fmtstring = "{0:d}f4x".format(self.missing_masses)
            everything = unpack(fmtstring, instring)
            offset = 0
            for i in range(6):
                if self.Npart[i] > 0 and self.mpart[i] == 0:
                    chunk = everything[offset:offset+self.Npart[i]]
                    self.m[i] = np.array(chunk)
                    offset += self.Npart[i]
        return self.m

    # Energy processing
    def unpackEnergy(self, instring):
        fmtstring = "{0:d}f4x".format(self.Ngas)
        everything = unpack(fmtstring, instring)

        self.Energy = np.array(everything)
        return self.Energy

    # Rho processing
    def unpackRho(self, instring):
        fmtstring = "{0:d}f4x".format(self.Ngas)
        everything = unpack(fmtstring, instring)

        chunk = everything[0:self.Ngas]
        self.Rho = np.array(chunk)
        return self.Rho


    # Utils
    def computeCOM(self, parts=range(6)):
        '''
        Computes center of mass for all the particle types
        given in the list parts,
        default all
        '''
        com = np.zeros(3)
        totM = 0.0
        for i in parts:
            for j in range(self.Npart[i]):
                com += self.pos[i][j, :]*self.m[i][j]
                totM += self.m[i][j]
        com = com/totM
        self.com = com

    def to_ascii(self):
        def get_tuple(key):
            if len(key) > 1:
                return tuple(i for i in self.sdata[key])

        id = np.concatenate(get_tuple('id'), axis = 0)
        mass = np.concatenate(get_tuple('mass'), axis = 0)
        pos = np.concatenate(get_tuple('pos'), axis = 0)
        vel = np.concatenate(get_tuple('vel'), axis = 0)
        u = self.sdata['u']
        rho = self.sdata['rho']
        u.resize(self.Ntot, refcheck=False)
        rho.resize(self.Ntot, refcheck=False)

        fmtstring = ['%8d', '%1.5e',
                     '% 1.5e', '% 1.5e', '% 1.5e',
                     '% 1.5e', '% 1.5e', '% 1.5e',
                     '% 1.5e', '% 1.5e']

        np.savetxt(self.fname+'.asc',
                   np.hstack([zip(id, mass), pos, vel, zip(u, rho)]),
                   fmt=fmtstring)

    def get_header(self):
        return self.sdata['header']

    def get_data_by_type(self, ptype):
        d = {}
        if ptype < 0 or ptype > 5:
            print(pytpe, "Invalid data type, use only 0, 1, 2, 3, 4 or 5")
            return d

        d['id'] = self.sdata['id'][ptype]
        d['mass'] = self.sdata['mass'][ptype]
        d['pos'] = self.sdata['pos'][ptype]
        d['vel'] = self.sdata['vel'][ptype]

        # Return Internal Energy and Density if the requested type is Gas.
        if ptype == 0:
            d['u'] = self.sdata['u']
            d['rho'] = self.sdata['rho']

        return d

    def print_data_by_type(self, ptype):
        if ptype < 0 or ptype > 5:
            print(pytpe, "Invalid data type, use only 0, 1, 2, 3, 4 or 5")
            return None
        for i in range(self.Npart[ptype]):
            pid = self.sdata['id'][ptype][i]
            mass = self.sdata['mass'][ptype][i]
            posx, posy, posz = self.sdata['pos'][ptype][i]
            velx, vely, velz = self.sdata['vel'][ptype][i]

            if ptype == 0:
                u = self.sdata['u'][i]
                rho = self.sdata['rho'][i]
                fmtstring = '%8d %1.5e '
                fmtstring += '% 1.5e % 1.5e % 1.5e '
                fmtstring += '% 1.5e % 1.5e % 1.5e '
                fmtstring += '% 1.5e % 1.5e'
                print(fmtstring % (pid, mass, posx, posy, posz, velx, vely, velz,\
                                   u, rho))

            else:
                fmtstring = '%8d %1.5e % 1.5e % 1.5e % 1.5e % 1.5e % 1.5e % 1.5e'
                print(fmtstring % (pid, mass, posx, posy, posz, velx, vely, velz))

## Print utils
#
#def print_header(snap):
#    for key, val in snap.sdata['header'].iteritems():
#        print(key, val)
#
#
#def print_pos(snap):
#    ptype = 0
#    for p in snap.sdata['pos']:
#        print("Type", ptype, p)
#        ptype += 1
#
#
#def print_vel(snap):
#    vtype = 0
#    for v in snap.sdata['vel']:
#        print("Type", vtype, v)
#        vtype += 1
#
#
#def print_id(snap):
#    itype = 0
#    for i in snap.sdata['id']:
#        print("Type", itype, i)
#        itype += 1
#
#
#def print_mass(snap):
#    mtype = 0
#    for m in snap.sdata['mass']:
#        print("Type", mtype, m)
#        mtype += 1

