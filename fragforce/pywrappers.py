
import numpy as np
import pickle as pickle
import ctypes
import os
from numpy.ctypeslib import ndpointer

# load the compiled c code
_LOCALDIR = os.path.dirname(__file__)
force=ctypes.CDLL(os.path.join(_LOCALDIR, "force.so"))

# set the surface triangulation as a global variable
srf_centers_sph, srf_normals_sph, srf_crosses_sph = \
    pickle.load(open(os.path.join(_LOCALDIR, "sphere_data6.p"), "rb"))


###############################################################################
###############################################################################
################                    WRAPPERS                     ##############
###############################################################################
###############################################################################

# The following are python wrappers for c functions defined in force.c. 
# See force.c for documentation.

# scale_triangulation
force.scale_triangulation.restype = None
force.scale_triangulation.argtypes = [ctypes.c_int,
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_scale_triangulations(a):
  srf_centers_scaled = np.zeros_like(srf_centers_sph)
  srf_normals_scaled = np.zeros_like(srf_normals_sph)
  srf_areas_scaled = np.zeros(len(srf_crosses_sph))
  SIZE = int(len(srf_centers_sph))
  force.scale_triangulation(SIZE, a, srf_centers_scaled, srf_areas_scaled, srf_normals_scaled, 
                        srf_centers_sph, srf_crosses_sph, srf_normals_sph)
  return [srf_centers_scaled, srf_areas_scaled, srf_normals_scaled]

# scale_edge
force.scale_edge.restype = None
force.scale_edge.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_scale_edge(axes, edge_normal_sph, edge_center_sph):
  edge_normal_scaled = np.zeros(3)
  edge_center_scaled = np.zeros(3)
  force.scale_edge(axes, edge_normal_sph, edge_center_sph, 
                   edge_normal_scaled, edge_center_scaled)
  return [edge_normal_scaled, edge_center_scaled]


# set_chi
force.set_chi.restype = None
force.set_chi.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), 
                          ctypes.c_double,
                          ctypes.c_double,
                          ctypes.c_double]
def py_set_chi(axes):
  chi = np.zeros(3)
  force.set_chi(chi, axes[0], axes[1], axes[2])
  return chi


# set_A
force.set_A.restype = None
force.set_A.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_set_A(w, a, L, X):
  A = np.zeros((3,3))
  force.set_A(A, w, a, L, X)
  return A


# set_L
force.set_L.restype = None
force.set_L.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                        ctypes.c_double,
                        ctypes.c_double,
                        ctypes.c_double]

def py_set_L(angles,gammadot):
    c = angles[0]
    s = angles[1]
    L = np.zeros([3,3])
    force.set_L(L,c,s,gammadot)
    return L


# set_farg
force.set_farg.restype = None
force.set_farg.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                           ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                           ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                           ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                           ctypes.c_double,
                           ctypes.c_double]

def py_set_farg(a, w, L, p0, mu):
  farg = np.zeros((3,3))
  force.set_farg(farg, a, w, L, p0, mu)
  return farg


# set_force_facets
force.set_force_facets.restype = None
force.set_force_facets.argtypes = [ctypes.c_int,
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_set_force_facets(farg, srf_normals_scaled, srf_areas_scaled):
  force_on_facets = np.zeros_like(srf_normals_scaled)
  SIZE = int(len(srf_normals_scaled))

  force.set_force_facets(SIZE, force_on_facets, 
                         farg, srf_normals_scaled,
                         srf_areas_scaled)
  return force_on_facets


# area_of_intersection
force.area_of_intersection.restype = ctypes.c_double
force.area_of_intersection.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_area_of_intersection(a, pn, px):
    aoi = force.area_of_intersection(a, pn, px)
    return aoi


# correct_pndotf
force.correct_pndotf.restype = ctypes.c_double
force.correct_pndotf.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_correct_pndotf(f, c, pn, px):
    nf = force.correct_pndotf(f, c, pn, px)
    return nf


# set_stress
force.set_stress.restype = ctypes.c_double
force.set_stress.argtypes = [ctypes.c_int, 
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_set_stress(a, fonfV, cV, pn, px):
    N = np.shape(fonfV)[0]
    stress = force.set_stress(N, a, fonfV, cV, pn, px)
    return stress


# set_force
force.set_force.restype = ctypes.c_double
force.set_force.argtypes = [ctypes.c_int, 
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_set_force(a, fonfV, cV, pn, px):
    N = np.shape(fonfV)[0]
    total_force = force.set_force(N, a, fonfV, cV, pn, px)
    return total_force



# set_force_Vectorized
force.set_force_Vectorized.restype = None
force.set_force_Vectorized.argtypes = [
                         ctypes.c_int,
                         ctypes.c_int,
                         ctypes.c_int,
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ctypes.c_double,
                         ctypes.c_double,
                         ctypes.c_double,
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS") ]

def py_set_force_Vectorized(aV, rotAngles, wV, pnV, pxV, gammadot, p0, mu):

    [NumTimes, NumPlanes, NumFacets] = [1,1,1]
    if aV.ndim > 1: 
        NumTimes = np.shape(aV)[0]
    if pnV.ndim > 1:
        NumPlanes = np.shape(pnV)[0]
    if srf_centers_sph.ndim > 1:
        NumFacets = np.shape(srf_centers_sph)[0]
    
    forces = np.zeros([NumTimes,NumPlanes])

    fonfV = np.zeros([NumFacets,3])

    force.set_force_Vectorized(
        NumTimes,
        NumPlanes,
        NumFacets,
        fonfV,
        forces,
        aV,
        rotAngles,
        wV,
        pnV,
        pxV,
        gammadot,
        p0,
        mu,
        srf_centers_sph,
        srf_crosses_sph,
        srf_normals_sph)

    return forces

###############################################################################
###############################################################################
################               EXTERNAL FUNCTIONS                ##############
###############################################################################
###############################################################################

 
def surface_forces(a, w, rotAngles, gammadot, p0, mu):
    """ Compute the force on the facets of the surface triangulation. 
    Wrapper for py_set_force_facets that takes more natural inputs and also
    returns the centers of the facets to be used in the fragmentation force 
    computations. Here fonfV can be thought of as the force vector heads 
    and srf_centers_scaled as the force vector tails.

    Note: computes the force (force_density * area of triangle) and not the
          force density. Since the areas are also returned, one can recover
          the force density by taking fonfV / srf_areas_scaled.

    Inputs:
        a          np.array(3)  ellipsoid axes lengths
        w          np.array(3)  angular velocity
        rotAngles  np.array(2)  cos theta, sin theta defining the rotation R. 
                   Note: only a single rotation, not an array of rotations
        gammadot   float  shear rate
        p0         float  external pressure
        mu         float  matrix viscosity

    Outputs:
        farg                  np.array([3,3])  farg (see force.c)
        fonfV                 np.array([N,3])  force vector on each facet
        srf_normals_scaled    np.array([N,3])  normals to facets
        srf_centers_scaled    np.array([N,3])  centers of facets
        srf_areas_scaled      np.array(N)      areas of facets

    """
    # set the velocity gradient in the body frame. see /tests/determine_R.py
    # for explanation of the direction of rotation.

    L = py_set_L(rotAngles,gammadot)
    
    # scale the surface triangulation quantities
    srf_centers_scaled, srf_areas_scaled, srf_normals_scaled = py_scale_triangulations(a)
    # compute farg
    farg = py_set_farg(a,w,L,p0,mu)
    # get the force on the facets
    fonfV = py_set_force_facets(farg, srf_normals_scaled, srf_areas_scaled)
    return [farg, fonfV, srf_centers_scaled, srf_areas_scaled]

