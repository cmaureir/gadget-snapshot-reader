gadget-snapshot-reader (GSR)
============================

[GADGET](http://www.mpa-garching.mpg.de/gadget/) is a cosmological SPH and
N-body code, written by [Volker Springel](http://www.mpa-garching.mpg.de/~volker/).

Since this code is widely used by the cosmological community,
and considering the growth of Python use among scientist,
[Patrick Brem](https://github.com/skele/gadget-snapshot-reader) decided to write a
simple script to read the snapshot and be able to access the data.

I decided to improve his code, and add more features, to provide a nice/simple
module to deal with the GADGET binary snapshots.

Be aware, the project is still **under development**.

Installation
--------------

To use this module, you can install is using the `setup.py`.

```bash
    python setup.py --user
```

Examples
---------

You can visit the `examples` directory to see a couple of ways to use
GSR.

Usage
------

For a simple transformation from the gadget-binary file
to ascii, you can do the following:

```python
from gsr import *

snap = Snapshot(filename)
snap.to_ascii()
```

Every ascii file is generated using the following order:

| id | mass | posx | posy | posz | velx | vely | velz | u | rho |


> It is possible to process more parameters, just enabling them in the
> constructor call, for example:
>
> > * snap = Snapshot(filename, enable_potential = True)
>
> Extra options are:
>
> > *  enable_potential
> > *  enable_accelerations
> > *  enable_entropy_production
> > *  enable_timesteps
> > *  enable_density
> > *  enable_smoothing_lenght

If you want to print certain information from the snapshot
by type:

```python
particle_type = 5  # From GADGET, only 0, 1, 2, 3, 4, 5
snap.print_data_by_type(particle_type)
```

To use the data from the snapshot, it is possible to get
a list with all the information of a particle type:

```python
particle_type = 5
data = snap.get_data_by_type(particle_type)
```

which will return a dictionary with the following structure:

```python
{'id':[id0, id1, ...],'mass':[mass0, mass1, ...], 'pos':[[pos0x, pos0y, pos0z], [pos1x, pos1y, pos2z]], 'vel':[[vel0x, vel0y, vel0z], [vel1x, vel1y, vel1z], ...]]
```

if the requested type is `0` (Gas) the dictionary will also include the Internal energy `u` and the Density `rho`.

Feedback
---------

Since the project is still under development, a lot of feature are missing,
and will be really nice if you can collaborate, asking for a new functionallity
or even submitting a pull request to improve the code.

