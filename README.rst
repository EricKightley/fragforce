Fragforce
---------

Numerical simulation of fluid ellipsoid in simple shear. Computes the motion,
deformation, surface force, and fragmentation force on an ellipsoidal droplet.
Companion to Kightley et al. (in prep). See /examples/usage.py for usage.

Source Code and Local Modifications
___________________________________
Code is distributed into several independent files:

::

  /fragforce/fragforce/surface_triangulation.py
  /fragforce/fragforce/deformation.py
  /fragforce/fragforce/pywrappers.py

None of these files imports any of the others, which will I hope
reduce some headache in making any local modifications. The only
potential issue is that pywrappers.py imports two files. One
of them is a pickle file containing information on surface
triangulations that has been hard-coded to save on computation 
time (actually the one it imports is fast to generate but this
system is in place in case we at some point want more accurate
integration, at which point it will be desirable to import the
triangulation). The other is a .so (see the Compiled C code
section below). Both of these files are in the same directory
as pywrappers.py and pywrappers imports them using a relative
filepath constructing using the os package. I think that this 
means just importing pywrappers like any other local file, 
for example using sys:

::

    >>> import sys
    >>> sys.path.append('/home/cleveruser/path/to/fragforce/')
    >>> import deformation

will work on anyone's system. If it doesn't, look in pywrappers.py
to modify the filepath for the pickle load and the ctypes import,
both at the start of the file. 


Dependencies
____________

ctypes and pickle are required. I get an error when including this
in setup.py so they are (hopefully temporarily) excluded from it.


Compiled C Code
_______________

The code is in Python, but calls on some C code (in /fragforce/force.c).
A compiled executable (/fragforce/force.so) is included, which is used by
/fragforce/pywrappers.py. I don't know if this will work on all systems;
in case it doesn't the c code is included (/fragforce/force.c and
/fragforce/force.h) as well as a shell script to compile it 
(/fragforce/compile.sc) it using gcc. If the .so doesn't work you'll have to 
compile force.c, and if you're not using a nix machine with gcc then the
script won't work either so you'll have to manually compile force.c and
make whatever changes necessary to /fragforce/pywrappers.py so it imports
your compiled one instead of the .so bundled here. 
