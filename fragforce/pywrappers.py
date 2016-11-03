
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

def py_scale_triangulation(a):
    srf_centers_scaled = np.zeros_like(srf_centers_sph)
    srf_normals_scaled = np.zeros_like(srf_normals_sph)
    srf_areas_scaled = np.zeros(len(srf_crosses_sph))
    SIZE = int(len(srf_centers_sph))
    force.scale_triangulation(SIZE,
                              np.ascontiguousarray(a),
                              np.ascontiguousarray(srf_centers_scaled),
                              np.ascontiguousarray(srf_areas_scaled),
                              np.ascontiguousarray(srf_normals_scaled), 
                              np.ascontiguousarray(srf_centers_sph),
                              np.ascontiguousarray(srf_crosses_sph),
                              np.ascontiguousarray(srf_normals_sph))
    return [srf_centers_scaled, srf_areas_scaled, srf_normals_scaled]

# scale_edge
force.scale_plane.restype = None
force.scale_plane.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_scale_plane(a_initial, a_current, pn_initial, px_initial):
    pn_scaled = np.zeros(3)
    px_scaled = np.zeros(3)
    force.scale_edge(np.ascontiguousarray(a_initial),
                     np.ascontiguousarray(a_current),
                     np.ascontiguousarray(pn_initial),
                     np.ascontiguousarray(px_initial), 
                     np.ascontiguousarray(pn_scaled),
                     np.ascontiguousarray(px_scaled))
    return [edge_normal_scaled, edge_center_scaled]


# set_chi
force.set_chi.restype = None
force.set_chi.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), 
                          ctypes.c_double,
                          ctypes.c_double,
                          ctypes.c_double]
def py_set_chi(a):
    chi = np.zeros(3)
    force.set_chi(np.ascontiguousarray(chi),
                  a[0],
                  a[1],
                  a[2])
    return chi


# set_L
force.set_L.restype = None
force.set_L.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                        ctypes.c_double,
                        ctypes.c_double,
                        ctypes.c_double]

def py_set_L(c,s,gammadot):
    L = np.zeros([3,3])
    force.set_L(np.ascontiguousarray(L),
                c,
                s,
                gammadot)
    return L


# set_A
force.set_A.restype = None
force.set_A.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_set_A(a, w, L, chi):
  A = np.zeros((3,3))
  force.set_A(np.ascontiguousarray(A),
              np.ascontiguousarray(a),
              np.ascontiguousarray(w),
              np.ascontiguousarray(L),
              np.ascontiguousarray(chi))
  return A


# set_farg
force.set_farg.restype = None
force.set_farg.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                           ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                           ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                           ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                           ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                           ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                           ctypes.c_double,
                           ctypes.c_double]

def py_set_farg(a, w, L, A, chi, p0, mu):
    farg = np.zeros((3,3))
    force.set_farg(np.ascontiguousarray(farg),
                   np.ascontiguousarray(a),
                   np.ascontiguousarray(w),
                   np.ascontiguousarray(L),
                   np.ascontiguousarray(A),
                   np.ascontiguousarray(chi),
                   p0,
                   mu)
    return farg


# set_force_density
force.set_force_density.restype = None
force.set_force_density.argtypes = [ctypes.c_int,
                                    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_set_force_density(farg, srf_normals_scaled):
    fdonfV = np.zeros_like(srf_normals_scaled)
    SIZE = int(len(srf_normals_scaled))

    force.set_force_density(SIZE,
                           np.ascontiguousarray(fdonfV), 
                           np.ascontiguousarray(farg),
                           np.ascontiguousarray(srf_normals_scaled))
    return fdonfV


# set_force_facets
force.set_force_facets.restype = None
force.set_force_facets.argtypes = [ctypes.c_int,
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_set_force_facets(fdonfV, srf_areas_scaled):
    SIZE = int(len(srf_areas_scaled))
    fonfV = np.zeros_like(fdonfV)

    force.set_force_facets(SIZE,
                           np.ascontiguousarray(fonfV), 
                           np.ascontiguousarray(fdonfV),
                           np.ascontiguousarray(srf_areas_scaled))
    return fonfV


# area_of_intersection
force.area_of_intersection.restype = ctypes.c_double
force.area_of_intersection.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_area_of_intersection(a, pn, px):
    aoi = force.area_of_intersection(np.ascontiguousarray(a), 
                                     np.ascontiguousarray(pn),
                                     np.ascontiguousarray(px))
    return aoi


# correct_pndotf
force.correct_pndotf.restype = ctypes.c_double
force.correct_pndotf.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_correct_pndotf(f, c, pn, px):
    nf = force.correct_pndotf(np.ascontiguousarray(f),
                              np.ascontiguousarray(c),
                              np.ascontiguousarray(pn),
                              np.ascontiguousarray(px))
    return nf

# sum_forces
force.sum_forces.restype = ctypes.c_double
force.sum_forces.argtypes = [ctypes.c_int, 
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def py_sum_forces(a, fonfV, cV, pn, px):
    N = np.shape(fonfV)[0]
    total_force = force.sum_forces(N, 
                                  np.ascontiguousarray(a),
                                  np.ascontiguousarray(fonfV),
                                  np.ascontiguousarray(cV),
                                  np.ascontiguousarray(pn),
                                  np.ascontiguousarray(px))
    return total_force

# frag_force
force.frag_force.restype = None
force.frag_force.argtypes = [ctypes.c_int,
                             ctypes.c_int,
                             ctypes.c_int,
                             ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                             ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
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
                             ctypes.c_int]

def py_frag_force(aV, RV, wV, pnV, pxV, gammadot, p0, mu, scale_planes_bool):
    [NTimes, NPlanes, NFacets] = [1,1,1]
    if aV.ndim > 1: 
        NTimes = np.shape(aV)[0]
    if pnV.ndim > 1:
        NPlanes = np.shape(pnV)[0]
    if srf_centers_sph.ndim > 1:
        NFacets = np.shape(srf_centers_sph)[0]
    forces = np.zeros([NTimes,NPlanes])
    force.frag_force(
        NTimes,
        NPlanes,
        NFacets,
        np.ascontiguousarray(forces),
        np.ascontiguousarray(aV),
        np.ascontiguousarray(RV),
        np.ascontiguousarray(wV),
        np.ascontiguousarray(pnV),
        np.ascontiguousarray(pxV),
        np.ascontiguousarray(srf_centers_sph),
        np.ascontiguousarray(srf_crosses_sph),
        np.ascontiguousarray(srf_normals_sph),
        gammadot,
        p0,
        mu,
        scale_planes_bool)
    return forces

