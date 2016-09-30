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
axes, R, w, T = frag.evolve(t0, t1, dt, a0, lam, mu, gammadot, Gamma, JustAngles=True)
# solve for the motion of the ellipsoid treating it as a solid
axes_s, R_s, w_s, T_s = frag.evolve_solid(t0, t1, dt, a0, gammadot, JustAngles=True)

#################################
# Surface Triangulation functions

centers, normals, edge_crosses = frag.generate_triangulation(6)


#################
# Force functions
# pick an arbitrary timepoint and compute the surface forces at that time
k = 17
farg, fonf, srf_centers_scaled, srf_areas_scaled = frag.surface_forces(axes[k], w[k], R[k], gammadot, p0, mu)

