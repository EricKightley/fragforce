
import numpy as np
import matplotlib.pyplot as plt
import fragforce as frc
import pickle
import csv

# set the constants

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

# set the initial axes
a0 = np.array([180., 160., 140.])

# compute the time interval
t0,t1,dt,tau,cap = frc.set_integration_params(a0, lam, mu, gammadot, Gamma)
# modify it so we get a perfectly symmetric interval
scl = 1.22
t1 = scl * t1
dt = (scl * dt) / 2.0
# get the rotations, axes and angular velocity
aV, RV, wV, T = frc.evolve(t0,t1,dt,a0,lam,mu,gammadot,Gamma)

# use this to check the interval we picked
#plt.plot(wV[:,2])
#plt.show()

M = 150
scale = .9
xcoords = np.linspace(-scale*a0[0],scale*a0[0],M)
pxV = np.array( [np.array([x,0,0]) for x in xcoords ] )
pnV = np.array([[1.0,0.0,0.0],]*M)
#theta = np.linspace(0, .5*np.pi, M)

N = np.shape(RV)[0]
K = int(np.floor(N / M))
force = np.zeros([M,M])

RV = RV[0::K,:]
wV = wV[0::K,:]
aV = aV[0::K,:]

if ( np.shape(RV)[0] == M+1 ):
    RV = RV[:-1]
    wV = wV[:-1]
    aV = aV[:-1]


forces = frc.py_frag_force(aV, RV, wV, pnV, pxV, gammadot, p0, mu)
pickle.dump([forces, xcoords, T],open("force_data.p","wb"))


ofile = open("forces.csv","w")
writer = csv.writer(ofile)
for row in forces:
  writer.writerow( [val for val in row] )
ofile.close()
 
ofile = open("bounds.csv","w")
writer = csv.writer(ofile)
writer.writerow( [ xcoords[0], xcoords[-1], T[0], T[-1] ] )
ofile.close()


"""
forces_iter = np.zeros([M,M])
for i in range(M):
    angles = RV[i]
    a = aV[i]
    w = wV[i]

    L = frc.py_set_L(angles, gammadot)
    chi = frc.py_set_chi(a)
    A = frc.py_set_A(a, w, L, chi)
    farg = frc.py_set_farg(a, w, L, A, chi, p0, mu)
    srf_centers_scaled, srf_areas_scaled, srf_normals_scaled = frc.py_scale_triangulation(a)
    force_density = frc.py_set_force_density(farg, srf_normals_scaled)
    force_facets = frc.py_set_force_facets(force_density, srf_areas_scaled)

    for j in range(M):
        pn = pnV[j]
        px = pxV[j]
        forces_iter[i,j] = frc.py_sum_forces(a, force_facets, srf_centers_scaled, pn, px)




for i in range(M):
    angles = RV[i]
    a = aV[i]
    w = wV[i]
    L = frc.py_set_L(angles, gammadot)
    A = frc.py_set_A(w, a, L)
    farg = frc.py_set_farg(A,a, w, L, p0, mu)
    srf_centers_scaled, srf_areas_scaled, srf_normals_scaled = frc.py_scale_triangulations(a)
    force_density = frc.py_set_force_density(farg, srf_normals_scaled)
    fonf = frc.py_set_force_facets(force_density, srf_areas_scaled)
    for j in range(M):
        pn = pnV[j]
        px = pxV[j]
        forces_iter[i,j] = frc.py_sum_forces(a, fonf, srf_centers_scaled, pn, px)

for i in range(M):
  R = RVs[i]
  w = wVs[i]
  a = aVs[i]
  fonfV, cV = frc.set_fonfV(a, w, R, gammadot, p0, mu)
  size = np.shape(fonfV)[0]
  for j in range(M):
    px = np.array([xcoords[j], 0., 0.])
    force[i,j] = frag.frag_force(a, w, R, gammadot, p0, mu, pn, px)


ymin = T[0]
ymax = T[-1]
xmin = xcoords[0]
xmax = xcoords[-1]

aspectratio = (xmax-xmin)/(ymax-ymin)

pickle.dump([force,xcoords,T,[xmin,xmax,ymin,ymax]],open("force_data.p","wb"))

fig = plt.figure(figsize = (6,6))
ax = plt.subplot(111)
cutoff = 1.0
maxval = cutoff*np.max(force)
minval = cutoff*np.min(force)

force[force > maxval] = maxval
force[force < minval] = minval


im = ax.imshow(force,
               extent = [xmin,xmax,ymin,ymax], 
               cmap = 'PiYG',
               aspect = aspectratio,
               origin = 'lower',
               interpolation = 'nearest')

cbar = fig.colorbar(im)

ax.set_xlabel("x-coordinate of intersecting plane (microns)")
ax.set_ylabel("time (seconds)")
plt.show()
"""
