Fragforce
---------

Numerical simulation of fluid ellipsoid in simple shear. Computes the motion,
deformation, surface force, and fragmentation force on an ellipsoidal droplet.
Companion to Kightley et al. (in prep). See /examples/usage.py for usage.


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
