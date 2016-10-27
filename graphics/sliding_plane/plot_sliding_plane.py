import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cPickle as pickle

# ellipse for subplots
theta = np.linspace(0,2*np.pi,300)
a = 1.5
b = 1.0
scale = 1.5
coords = np.array([a*np.cos(theta), b*np.sin(theta)])
line = np.array([np.zeros(50), 1.2 * b * np.linspace(-1, 1, 50)])

# import force data
force, X, Y = pickle.load(open("force_data.p","rb"))

force *= 1/np.max(np.abs(force))

xmin = X[0]
xmax = X[-1]
ymin = Y[0]
ymax = Y[-1]

# set up rotations
ang2 = -np.pi/4
ang3 = -np.pi/2 
ang4 = -3*np.pi/4
ang5 = -np.pi

c2 = np.cos(ang2)
s2 = np.sin(ang2)
c3 = np.cos(ang3)
s3 = np.sin(ang3)
c4 = np.cos(ang4)
s4 = np.sin(ang4)
c5 = np.cos(ang5)
s5 = np.sin(ang5)

R2 = np.array([[c2, -s2],[s2, c2]])
R3 = np.array([[c3, -s3],[s3, c3]])
R4 = np.array([[c4, -s4],[s4, c4]])
R5 = np.array([[c5, -s5],[s5, c5]])

coords2 = np.dot(R2, coords)
coords3 = np.dot(R3, coords)
coords4 = np.dot(R4, coords)
coords5 = np.dot(R5, coords)
line2 = np.dot(R2, line)
line3 = np.dot(R3, line)
line4 = np.dot(R4, line)
line5 = np.dot(R5, line)

fig = plt.figure(figsize = (7,6))
gs = gridspec.GridSpec(6,7)
gs.update(hspace = 1.3, wspace = 1.3)
ax0 = plt.subplot(gs[0:-1,1:-1])

ax1 = plt.subplot(gs[4,0])
ax2 = plt.subplot(gs[3,0])
ax3 = plt.subplot(gs[2,0])
ax4 = plt.subplot(gs[1,0])
ax5 = plt.subplot(gs[0,0])

ax6 = plt.subplot(gs[5,1])
ax7 = plt.subplot(gs[5,2])
ax8 = plt.subplot(gs[5,3])
ax9 = plt.subplot(gs[5,4])
ax10 = plt.subplot(gs[5,5])

ax11 = plt.subplot(gs[5,0])

ax12 = plt.subplot(gs[:-1,6])

ax = ax1
ax.plot(coords[0,:],coords[1,:],linewidth = 1, color = 'k')
ax.plot(line[0,:],line[1,:],linewidth = 1, color = 'k')

ax = ax2
ax.plot(coords2[0,:],coords2[1,:],linewidth = 1, color = 'k')
ax.plot(line2[0,:],line2[1,:],linewidth = 1, color = 'k')

ax = ax3
ax.plot(coords3[0,:],coords3[1,:],linewidth = 1, color = 'k')
ax.plot(line3[0,:],line3[1,:],linewidth = 1, color = 'k')

ax = ax4
ax.plot(coords4[0,:],coords4[1,:],linewidth = 1, color = 'k')
ax.plot(line4[0,:],line4[1,:],linewidth = 1, color = 'k')

ax = ax5
ax.plot(coords5[0,:],coords5[1,:],linewidth = 1, color = 'k')
ax.plot(line5[0,:],line5[1,:],linewidth = 1, color = 'k')

ax = ax6
ax.plot(coords[0,:],coords[1,:],linewidth = 1, color = 'k')
ax.plot(line[0,:]-.8*a,line[1,:],linewidth = 1, color = 'k')

ax = ax7
ax.plot(coords[0,:],coords[1,:],linewidth = 1, color = 'k')
ax.plot(line[0,:]-.4*a,line[1,:],linewidth = 1, color = 'k')

ax = ax8
ax.plot(coords[0,:],coords[1,:],linewidth = 1, color = 'k')
ax.plot(line[0,:],line[1,:],linewidth = 1, color = 'k')

ax = ax9
ax.plot(coords[0,:],coords[1,:],linewidth = 1, color = 'k')
ax.plot(line[0,:]+.4*a,line[1,:],linewidth = 1, color = 'k')

ax = ax10
ax.plot(coords[0,:],coords[1,:],linewidth = 1, color = 'k')
ax.plot(line[0,:]+.8*a,line[1,:],linewidth = 1, color = 'k')

ax = ax0
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.contour(force)
ax.contour(force)
ax.contour(force)
aspectratio = (xmax-xmin)/(ymax-ymin)
Y = np.linspace(Y[0],Y[-1],np.shape(force)[0])
ax.contour(X,Y,force,colors='k')
im = ax.imshow(force,
               extent = [xmin,xmax,ymin,ymax], 
               #cmap = 'OrRd', 
               #cmap = 'bwr',
               cmap = 'PiYG',
               #cmap = 'PuRd',
               aspect = aspectratio,
               origin = 'lower',
               interpolation = 'lanczos')
               #interpolation = 'nearest')
               #interpolation = 'bilinear')
               #interpolation = 'none')
#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.8, -.2, 0.1,.9])
ax.set_xlabel("x-coordinate of intersecting plane (microns)")
ax.set_ylabel("time (seconds)")
#ax.plot(np.linspace(0,1,100), np.linspace(0,2,100),color = 'k')

ax = ax12
cbar = fig.colorbar(im, cax = ax)
ax.set_xticks([])
ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)


# hide the splines on all axes except the main one

all_axes = fig.get_axes()
all_axes = all_axes[1:-1]
for ax in all_axes:
    ax.set_xlim([-scale*a,scale*a])
    ax.set_ylim([-scale*a,scale*a])
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

plt.show()
#plt.savefig("stress_sliding.png",dpi=200,transparent=True, bbox_inches='tight', pad_inches=0.1)
