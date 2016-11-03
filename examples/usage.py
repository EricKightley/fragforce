import numpy as np
import fragforce as frag

# ----- Inputs -----
lam = 50                                 # viscosity ratio, unitless
mu_si = 8.953e-4                         # matrix viscosity,  Pa s=(N s)/(m ^2)
gammadot_si = 10.                        # shear rate, 1/s
Gamma_si = 4.1e-9                        # interfacial tension, N/m
p0_si = 0.0                              # external pressure, Pa = N / m^2

# ----- Unit Conversions -----
gammadot = gammadot_si                   # 1/s
mu = mu_si * 10**(-6)                    # microNewton s / micrometer^2
Gamma = Gamma_si                         # microNewton / micrometer
p0 = p0_si * 10**(-6)                    # microNewton / micrometer^2

# ----- Set Initial Axes -----
a0 = np.array([180., 160., 120.])        # micrometer

#######################
# Deformation functions

# compute the integration parameters
t0, t1, dt, tau, cap = frag.set_integration_params(a0, lam, mu, gammadot, Gamma)
# solve for the motion of the ellipsoid over time
aV, RV, wV, T = frag.evolve(t0, t1, dt, a0, lam, mu, gammadot, Gamma, JustAngles=True)
# solve for the motion of the ellipsoid treating it as a solid
aVsolid, RVsolid, wVsolid, Tsolid = frag.evolve_solid(t0, t1, dt, a0, gammadot, JustAngles=True)

#################
# Force functions

# choose a time index
k = 17
a = aV[k]
w = wV[k]
R = RV[k]

# obtain the scaled surface triangulations
srf_centers_scaled, srf_areas_scaled, srf_normals_scaled = frag.py_scale_triangulation(a)
# set L
L = frag.py_set_L(R[0],R[1],gammadot)
# set chi
chi = frag.py_set_chi(a)
# set A
A = frag.py_set_A(a, w, L, chi)
# set the force argument
farg = frag.py_set_farg(a, w, L, A, chi, p0, mu)
# compute force density
fdonfV = frag.py_set_force_density(farg, srf_normals_scaled)
# compute the force on the facets
fonfV = frag.py_set_force_facets(fdonfV, srf_areas_scaled)
# pick a plane and compute the fragmentation force against the plane
pn = np.array([1.0, 0.0, 0.0])
px = np.zeros(3)
# sum the forces
total_force = frag.py_sum_forces(a, fonfV, srf_centers_scaled, pn, px)

# wrapper for the above in C
sample_force = frag.py_frag_force(a, R, w, pn, px, gammadot, p0, mu, False)
# do it again but scale the edge (won't change anything here)
sample_force_scaled = frag.py_frag_force(a, R, w, pn, px, gammadot, p0, mu, True)
# now do it over all of the times and a bunch of edges
pnV = np.random.rand(150,3)
pxV = np.random.rand(150,3) * aV[0]
forcesV = frag.py_frag_force(aV, RV, wV, pnV, pxV, gammadot, p0, mu, True)
