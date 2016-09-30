import numpy as np
import scipy.special as sp
from scipy.integrate import ode
from scipy.integrate import odeint

#----- Transcribed from TJW Code -----

def ellip1(xt,yt,zt):
    """Incomplete elliptic integral function.  This routine is 
    transliterated from the function "RD" from Eric Wetzel's
    calcesh.f file.  I don't know where Eric got it from."""

    C1=3.0/14.0
    C2=1.0/6.0
    C3=9.0/22.0 
    C4=3.0/26.0
    C5=.25*C3
    C6=1.5*C4
    sigma=0.
    fac=1.
    rtx=np.sqrt(xt)
    rty=np.sqrt(yt)
    rtz=np.sqrt(zt)
    alamb=rtx*(rty+rtz)+rty*rtz
    sigma=sigma+fac/(rtz*(zt+alamb))
    fac=.25*fac
    xt=.25*(xt+alamb)
    yt=.25*(yt+alamb)
    zt=.25*(zt+alamb)
    ave=.2*(xt+yt+3.*zt)
    delx=(ave-xt)/ave
    dely=(ave-yt)/ave
    delz=(ave-zt)/ave
    ea=delx*dely
    eb=delz*delz
    ec=ea-eb
    ed=ea-6.*eb
    ee=ed+ec+ec
    rd=3.*sigma+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee) + \
                     delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*np.sqrt(ave))
    return rd

def ellip2(xt,yt,zt):
    """  Incomplete elliptic integral function.  This routine is 
    transliterated from the function "RF" from Eric Wetzel's
    calcesh.f file.  I don't know where Eric got it from. """

    C1=1.0/24.0
    C2=0.1
    C3=3.0/44.0
    C4=1.0/14.0
    rtx=np.sqrt(xt)
    rty=np.sqrt(yt)
    rtz=np.sqrt(zt)
    alamb=rtx*(rty+rtz)+rty*rtz
    xt=.25*(xt+alamb)
    yt=.25*(yt+alamb)
    zt=.25*(zt+alamb)
    ave=1.0/3.0*(xt+yt+zt)
    delx=(ave-xt)/ave
    dely=(ave-yt)/ave
    delz=(ave-zt)/ave
    e2=delx*dely-delz**2
    e3=delx*dely*delz
    rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/np.sqrt(ave)
    return rf

def vec2tens(Vec):
    """Convert 6x1 contracted SYMMETRIC tensor to 3x3 matrix form"""
    tens = np.array([[Vec[0], Vec[5], Vec[4]],
                   [Vec[5], Vec[1], Vec[3]],
                   [Vec[4], Vec[3], Vec[2]]])
    return tens

def tens2vec(Tens):
    """Convert symmetric 2nd order tensor Tens (in 3x3 matrix) form to 6x1 
    column vector. EPKNOTE: no column vectors in numpy. Check multiplication 
    in subsequent functions. """
    vec = np.array([Tens[0,0], Tens[1,1], Tens[2,2], 
                  Tens[1,2], Tens[2,0], Tens[0,1]])
    return vec
  
def vec2skewTens(Vec):
    """Convert 6x1 contracted ANTISYMMETRIC tensor into 3x3 matrix form"""
    skewTens = np.array([[ Vec[0],  Vec[5], -Vec[4]],
                       [-Vec[5],  Vec[1],  Vec[3]],
                       [ Vec[4], -Vec[3],  Vec[2]]])
    return skewTens

def dropAxes(Gv):
    """Find the droplet semi-axes, daxes = [a b c], for the droplet shape 
    tensor Gv, which is stored in 6x1 column form.  
    Also find the angle theta (in degrees) between the major axis and the x1 axis
    (assuming that this axis lies in the 1-2 plane)"""
    D, V      = np.linalg.eigh(vec2tens(Gv))  # eigenvectors V and eigenvalues D
    daxesPrim = np.sqrt(1.0/np.abs(D))        # drop axes, as a row vector
    # sort axis lengths in descending order
    daxes     = np.sort(daxesPrim)[::-1]      # [::-1] reverses the order 
    sortOrder = np.argsort(daxesPrim)[::-1]
    # compute theta, in degrees I guess
    theta = np.arctan2(V[1,sortOrder[0]] , V[0,sortOrder[0]]) * 180/np.pi
    if theta < -90:
        theta = theta+180
    elif theta > 90:
        theta = theta-180

    return [daxes, theta]

def eshtens(G, lam):
    """Compute concentration tensors for droplet shape tensor G and viscosity 
    ratio lam. Bm and Cm are the concentration tensors for strain rate and 
    vorticity, and are stored in 6x6 matrix form. G is stored in 3x3 matrix form. 
    lam is a scalar. Equation numbers refer to Eric Wetzel's thesis"""

    #---- Principal values and axes of G ----
    evalsPrim, Evecs =  np.linalg.eigh(G)
    evals = np.abs(evalsPrim)  # no negative evals
    dropaxesPrim = 1./np.sqrt(evals)

    # sort axis lengths in ascending order
    dropaxes = np.sort(dropaxesPrim)
    sortOrder = np.argsort(dropaxesPrim)

    # define a rotation by the sorted eigenvectors 
    RPrim = np.array([Evecs[:,ind] for ind in sortOrder])
    Rt = np.reshape(RPrim, (3,3))
    R = np.transpose(Rt) # I guess matlab uses fortan-like indexing of arrays

    #---- Axis ratios (B.85) ----
    a3,a2,a1 = dropaxes
    C = a3/a1             # C = c/a
    D = a3/a2             # D = c/b
    TINY = 10**(-14)      # tol for using C or D = 0 formulae

    #---- Build Eshelby tensors Sm and Tm in the principal axis system ---
    Sm = np.zeros([6,6])
    Tm = np.zeros([6,6])

    C2 = C**2

    # independent components for various cases of C and D:
    if D < TINY:
        # disk shape (U corner)
        Sm[0,0] = 1.0
        Sm[0,2] = 1.0

    elif (1-(D-C)) < TINY:
        # circular cylinder (B corner)
        Sm[0,0] = 0.75
        Sm[1,1] = 0.75
        Sm[0,2] = 0.50
        Sm[1,0] = 0.25

    elif 1-C < TINY:
        # sphere (T corner)
        Sm[0,0] = 0.60
        Sm[1,1] = 0.60
        Sm[2,2] = 0.60
        Sm[2,1] = 0.20
        Sm[0,2] = 0.20
        Sm[1,0] = 0.20

    elif C < TINY:
        # ellipsoidal cylinder (UB edge)
        Sm[0,0] = (1+2*D) / (1+D)**2
        Sm[1,1] = D*(2+D) / (1+D)**2
        Sm[0,2] =   1 / (1+D)
        Sm[1,0] = D**2 / (1+D)**2

    elif (1-D) < TINY:
        # prolate spheriod (BT edge)
        exp1 = 0.5*np.log( (1+np.sqrt(1-C2)) / (1-np.sqrt(1-C2)) )
        it   =  (np.sqrt(1-C2) - C2*exp1) / (2*((1-C2)**(1.5)))
        Sm[0,0] = 0.75 * (1-3*it*C2) / (1-C2)
        Sm[1,1] = Sm[0,0]
        Sm[2,2] = (3-6*it-C2)   / (1-C2)
        Sm[2,1] = C2 * (3*it-1) / (1-C2)
        Sm[0,2] = (3*it-1) / (1-C2)
        Sm[1,0] = 0.25*(1-3*it*C2) / (1-C2)

    elif (D-C) < TINY:
        # oblate spheriod (TU edge)
        it = C * (np.arccos(C) - C*np.sqrt(1-C2)) / (2*((1-C2)**(1.5)))
        Sm[0,0] = (1+6*C2*it-3*C2) / (1-C2)
        Sm[1,1] = 0.75*(3*it-C2) / (1-C2)
        Sm[2,2] = Sm[1,1]
        Sm[2,1] = 0.25*(3*it-C2) / (1-C2)
        Sm[0,2] = (1-3*it) / (1-C2)
        Sm[1,0] = C2 * (1-3*it) / (1-C2)

    else:
        # general ellipsoid

        # Elliptic integral quantities (B.87, 88, 91, and Wetzel's "calcesh.f")
        RD = ellip1(C2, (C/D)**2, 1)
        RF = ellip2(C2, (C/D)**2, 1)
        F  = np.sqrt(1-C2)*RF
        E  = F - ((1-(C/D)**2)*np.sqrt(1-C2)*RD) / 3.0

        Ja = (C2*D*(F-E)) / ((D**2-C2)*np.sqrt(1-C2))
        Jc = (np.sqrt(1-C2) - D*E) / ((1-D**2) * np.sqrt(1-C2))
        Jb = 1 - Ja - Jc

        # Eqns (B.79-84)

        # Main components.  From Eqns. (B.79-84) with 1 and 3 indices swapped
        #   The 1-3 swapping occurs because Wetzel's appendix B associates S(1,1)
        #   with the longest semi-axis of the ellipse, which is the smallest
        #   eigenvalue/vector of G.
        Sm[0,0] = 1 + C2*(Ja-Jc) / (1-C2) + D**2*(Jb-Jc) / (1-D**2)
        Sm[1,1] = 1 + (Jb-Jc) / (1-D**2) + C2*(Ja-Jb) / (D**2-C2)
        Sm[2,2] = 1 + D**2*(Ja-Jb) / (D**2-C2) + (Ja-Jc) / (1 -C2)
        Sm[2,1] = C2*(Jb-Ja) / (D**2-C2)
        Sm[0,2] = (Jc-Ja) / (1-C2)
        Sm[1,0] = D**2*(Jc-Jb) / (1-D**2)

    # End of main principal-axis components for various cases of C and D


    #---- Fill in remaining principal-axis components using  (B.103-105) ----
    Sm[1,2] = 1.0 - Sm[0,2] - Sm[2,2]
    Sm[2,0] = 1.0 - Sm[1,0] - Sm[0,0]
    Sm[0,1] = 1.0 - Sm[2,1] - Sm[1,1]

    Sm[3,3] = 0.5*( Sm[1,2] + Sm[2,1] )
    Sm[4,4] = 0.5*( Sm[2,0] + Sm[0,2] )
    Sm[5,5] = 0.5*( Sm[0,1] + Sm[1,0] )

    Tm[3,3] = 0.5*( Sm[2,1] - Sm[1,2] )  # these eqns. match (B.105).  In
    Tm[4,4] = 0.5*( Sm[0,2] - Sm[2,0] )  # calcesh.f they differ by a sign,
    Tm[5,5] = 0.5*( Sm[1,0] - Sm[0,1] )  # which is later corrected in the C eqn.


  #---- 6x6 rotation matrices and associated quantities -----
    Qa=np.array(
         [[R[0,0]*R[0,0], R[0,1]*R[0,1], R[0,2]*R[0,2], 
           R[0,1]*R[0,2], R[0,2]*R[0,0], R[0,0]*R[0,1]],
          [R[1,0]*R[1,0], R[1,1]*R[1,1], R[1,2]*R[1,2], 
           R[1,1]*R[1,2], R[1,2]*R[1,0], R[1,0]*R[1,1]],
          [R[2,0]*R[2,0], R[2,1]*R[2,1], R[2,2]*R[2,2], 
           R[2,1]*R[2,2], R[2,2]*R[2,0], R[2,0]*R[2,1]],
          [R[1,0]*R[2,0], R[1,1]*R[2,1], R[1,2]*R[2,2], 
           R[1,1]*R[2,2], R[1,2]*R[2,0], R[1,0]*R[2,1]],
          [R[2,0]*R[0,0], R[2,1]*R[0,1], R[2,2]*R[0,2], 
           R[2,1]*R[0,2], R[2,2]*R[0,0], R[2,0]*R[0,1]],
          [R[0,0]*R[1,0], R[0,1]*R[1,1], R[0,2]*R[1,2], 
           R[0,1]*R[1,2], R[0,2]*R[1,0], R[0,0]*R[1,1]]])

    Qb=np.array(
         [[0, 0, 0, R[0,2]*R[0,1], R[0,0]*R[0,2], R[0,1]*R[0,0]],
          [0, 0, 0, R[1,2]*R[1,1], R[1,0]*R[1,2], R[1,1]*R[1,0]],
          [0, 0, 0, R[2,2]*R[2,1], R[2,0]*R[2,2], R[2,1]*R[2,0]],
          [0, 0, 0, R[1,2]*R[2,1], R[1,0]*R[2,2], R[1,1]*R[2,0]],
          [0, 0, 0, R[2,2]*R[0,1], R[2,0]*R[0,2], R[2,1]*R[0,0]],
          [0, 0, 0, R[0,2]*R[1,1], R[0,0]*R[1,2], R[0,1]*R[1,0]]])

    Q  = Qa+Qb  # rotation matrix for   symmetric tensors
    Qu = Qa-Qb  # rotation matrix for unsymmetric tensors

    #---- contracted 4th-order identity tensor and its inverse ----
    Id4 = np.diag(np.array([1, 1, 1, 0.5, 0.5, 0.5])) 
    R4  = np.diag(np.array([1, 1, 1, 2,   2,   2  ]))

    #---- Rotate Eshelby tensors into laboratory coordinates ----
    Qi = np.linalg.inv(Q)
    Sm = Q.dot(Sm.dot(R4.dot(Qi.dot(Id4))))  # Eshelby tensor in lab coords.
    Tm = Qu.dot(Tm.dot(R4.dot(Qi.dot(Id4)))) # alternate tensor in lab coords

    #---- Calculate the concentration tensors ----
    Smsi = np.linalg.inv(Id4 - (1-lam)*Sm)
    Bm   = Id4.dot(Smsi.dot(Id4))                    # eqn. (B.57) 
    # Cm = (1-lam)*Tm * inv(Id4-(1-lam)*Sm) * Id4    # eqn. (B.66)
    Cm = (1-lam)*Tm.dot(R4.dot(Bm))                  # eqn. (B.62)

    return [Bm, Cm, Sm, Tm]

#----- Begin EPK Code ----


def dgdt(Gv, t, L, lam, mu, Gamma):
    """ Evaluate dG/dt as per the Tucker Jackson Wetzel model. This function 
    is to be integrated numerically. Translated from the TJW matlab code.

    Input:
        Gv          6x1 np array; contracted vector form of G(t), the shape
                    tensor. 
        L           function, 3x3 array. The velocity gradient, which may have a
                    time dependence. Must therefore be defined as a function of 
                    time. See the testing section at the end of this document 
                    (after `if __name__ == "main": ') for an example on how to 
                    define L.
        lam         float. Viscosity ratio. 
        mu          float. "Matrix" (external fluid) viscosity. 
        Gamma       float. Interfacial tension.  
        vol         float. Initial volume of the ellipsoid. Used to preserve 
                    volume
    Output: 
        dgdt        6x1 np array. dG/dt, in contracted vector form. """

    # Compute axis lengths
    G = vec2tens(Gv)
    evals, V = np.linalg.eigh(G)
    a=1./np.sqrt(np.abs(evals))
 
    # Compute appropriate elliptic integrals of the first (ellipk) and second 
    # (ellipe) kinds, store them in F (3x1 array)
    F = np.zeros(3)			  
    if a[1] > a[2]:
        K,E = sp.ellipk(1-(a[2]/a[1])**2), sp.ellipe(1-(a[2]/a[1])**2)
        F[0] = E/a[2]
    else:
        K,E= sp.ellipk(1-(a[1]/a[2])**2), sp.ellipe(1-(a[1]/a[2])**2)
        F[0] = E/a[1]

    if a[0] > a[2]:
        K,E = sp.ellipk(1-(a[2]/a[0])**2), sp.ellipe(1-(a[2]/a[0])**2)
        F[1] = E/a[2]
    else:
        K,E = sp.ellipk(1-(a[0]/a[2])**2), sp.ellipe(1-(a[0]/a[2])**2)
        F[1] = E/a[0]

    if a[1] > a[0]:
        K,E = sp.ellipk(1-(a[0]/a[1])**2), sp.ellipe(1-(a[0]/a[1])**2)
        F[2] = E/a[0]
    else:
        K,E = sp.ellipk(1-(a[1]/a[0])**2), sp.ellipe(1-(a[1]/a[0])**2)
        F[2] = E/a[1]

    # Modify the elliptic integrals (not sure what c is)
    c = (40.0*(lam+1))/(19.0*lam+16)
    F = c * F * ((2*Gamma)/(np.pi*lam*mu))

    # Set an array to help with contracted-notation (not sure what this does)
    R4  = np.diag(np.array([1, 1, 1, 2, 2, 2]))

    # Compute the interfacial tension in the direction of the droplet axes
    interfacial = np.diag(-F + np.sum(F) * np.ones(3)/3)
                                          
    # "Go backwards from eig to get back to normal axes" (do not understand)
    I = V.dot(interfacial.dot(np.linalg.inv(V)))

    # Get interfacial tension as a vector
    Iv = tens2vec(I)

    # Compute Eshelby tensors 
    Bm, Cm, Sm, Tm = eshtens(G, lam)

    # Compute vorticity and deformation rate tensors from velocity gradient
    Llocal  = L    # replace L with L(t) if we want L as a function
    Ltlocal = np.transpose(Llocal)
    W = (Llocal-Ltlocal)/2.0
    D = (Llocal+Ltlocal)/2.0
    Dv = tens2vec(D)

    # "Inclusion velocity gradient" (do not understand)
    Lstar = W + vec2tens(Bm.dot(R4.dot(Dv))) + \
            vec2skewTens(Cm.dot(R4.dot(Dv))) + \
            lam * (vec2tens(Sm.dot(R4.dot(Bm.dot(R4.dot(Iv))))))
    Lstart = np.transpose(Lstar)

    # Get dgdt
    dgdt_tens = (-Lstart.dot(G) - G.dot(Lstar))
    # Ensure symmetry
    dgdt_tens = .5 * (dgdt_tens + dgdt_tens.T)
    # Get dgdt
    dgdt = tens2vec(dgdt_tens)
    return dgdt


def integrate_dgdt(t0,t1,dt,a0,lam,mu,gammadot,Gamma):
    """ Integrate the function dgdt using scipy.integrate.ode. Note that in
    general the velocity gradient can be anything; here we assume it is 0 
    everywhere except du/dy = gammadot. 

    Inputs:
        t0          float, start time
        t1          float, stop time (this may be changed)
        dt          float, time step
        a0          initial axes lengths
        lam         float. Viscosity ratio. 
        mu          float. "Matrix" (external fluid) viscosity. 
        gammadot    float. shear rate.
        Gamma       float. Interfacial tension. 

    Outputs:
        Y           6xM np array, where shapesV[:,i] is the shape tensor in 
                    contracted vector form at the ith time interval. M is
                    determined by the integration parameters:
                    M := floor( (stoptime-starttime)/dt )
        T           np.array(M) of times. 
    """
    # set up the velocity gradient L defined by du/dy=gammadot
    L = np.zeros([3,3])
    L[0,1] = gammadot

    # set up the initial shape tensor
    G0 = np.diag(1/a0**2)
    G0v = tens2vec(G0)

    T = np.arange(t0,t1,dt)
    args = (L, lam, mu, Gamma)
    Yvec = odeint(dgdt, G0v, T, args=args)
    Y = np.array([vec2tens(y) for y in Yvec])
    return [Y,T]


def eigensort(newevecs, newevals, oldevecs):
    """Corrects the order and the signs of a set of eigenvectors (and the order 
    of corresponding eigenvalues) to "match" (*) that of another, which is 
    taken to be the preceding set of eigenvectors in a time-series. Assuming
    a small enough timestep (hence a small enough rotation), a current evec 
    corresponds to the preceding evec with which it makes the largest 
    dot-product (in magnitude, to avoid sign errors); this is how we sort the 
    current evecs. After sorting, we correct the signs.  

    (*) When np.linalg.eigh is called on a time-series of shape tensors, two
    problems can arise: for a given set of eigenvectors, the order and/or
    the signs can be wrong. For example, suppose the solution S(t) has 
    eigenvectors [v1(t), v2(t), v3(t)]. We want the solution S(t+dt) to have 
    eigenvectors [v1(t+dt), v2(t+dt), v3(t+dt)]. Sometimes the solver returns 
    these in the wrong order (it has no memory of the last set it obtained). 
    Since linalg.eigh returns unit vectors, the only error that can arise from 
    the fact that evecs are closed under scalar multiplication is a sign error.
    This also happens sometimes. Sequential sets of eigenvectors V(t) and 
    V(t+dt) "match" if no such ordering or sign errors are present.   

    Inputs:
        newevecs          np.array([3,3]) of eigenvectors to be sorted and corrected
        newevals          np.array(3) of eigenvalues to be sorted
        oldevecs          np.array([3,3]) of eigenvectors against which to sort

    Outputs:
        outevecs          np.array[(3,3)] of sorted and corrected eigenvectors
        outevals          np.array(3) of sorted eigenvalues

    """
    dots = np.dot(np.transpose(newevecs), oldevecs)
    indices = np.absolute(dots).argmax(axis=0)
    outevecs = newevecs[:,indices]
    outevals = newevals[indices]
    # check that sign of evecs are correct
    for i in range(3):
        if np.dot(outevecs[:,i],oldevecs[:,i]) < 0:
          # this would mean the new one points in the wrong dir.
            outevecs[:,i] = -outevecs[:,i]
    return(outevecs, outevals)

def shapetensors_to_axes_rots(Y, X0 = np.identity(3)):
    """ Given a time-series of shape tensors and the basis vectors of 
    the body coordinate system at the initial time,  returns time-series of 
    eigenvectors and axes lenghts of these tensors. These are sorted and 
    corrected as per eigensort. 

    Inputs:
        Y               np.array([M,6]) of shape tensors in contracted tensor
                        notation,  from solving dgdt (time-series)
        X0              np.array[(3,3)] whose columns are the basis vectors of the
                        body frame at time 0 in order. 
                        (i.e., first column is the x axis). 
    Outputs:
        axes            np.array[(M,3)] of axes lengths
        R               np.array[(3,3,M)] of rotation matrices. 
    """


    #----- Initialize outputs -----
    M = Y.shape[0]
    evals = np.zeros([M,3])
    evecs = np.zeros([M,3,3])
    evals_sorted = np.zeros([M,3])
    evecs_sorted = np.zeros([M,3,3])

    #----- Obtain (unsorted) eigensystem
    for i in range(M):
        evals[i], evecs[i] = np.linalg.eigh(Y[i])

    #----- Sort eigensystem -----

    # set first entries (used to sort the rest)
    V = evecs[0]
    l = evals[0]
    # V is the set of eigenvectors of G0, the initial shape tensor. This is to
    # say that V.T G0 V = Gdiag, where Gdiag is diagonal; i.e., 
    # G0 = V Ggiad V.T, meaning that V.T is the rotation matrix sending an
    # ellipsoid aligned with the lab frame axes into the G0 frame; (final) i.e.,
    # V.T ought to be X0, the matrix of basis vectors of the body frame. 

    # sort initial ones
    indices = l.argsort()
    l_sorted = l[indices] # this is sorted in INCREASING order, so that axes 1 = 1/l are DECREASING order
    V_sorted = V[:,indices]

    evals_sorted[0] = l_sorted
    evecs_sorted[0] = V_sorted
  
    # sort each adjacent set of evecs iteratively
    for iii in range(M)[1:]:   # note we skip 0th entry
        newevecs = evecs[iii]
        newevals = evals[iii]
        oldevecs = evecs_sorted[iii-1]
        # sort newevecs and newevals by oldevecs
        evecs_sorted[iii], evals_sorted[iii] = eigensort(newevecs, 
                                                   newevals, oldevecs)
    axes = 1/np.sqrt(evals_sorted)
    R = evecs_sorted
    return(axes, R)

def angular_velocity(R, dt):
    """ Compute the angular velocity of a rotating object whose orientation at
    time t is given by R(t). 

    Inputs:
        R           np.array[(N,3,3)], each entry is rotation matrix at time t
        dt          float. time step.


    Outputs:
        w           np.array[(N,3)], angular velocity at each time. 
    """
    N = np.shape(R)[0]
    w = np.zeros([N,3])
    for i in range(N-1):
        Wx = np.dot( R[i+1] - R[i] , R[i].T ) / dt
        w[i] = np.array([Wx[2,1], Wx[0,2], Wx[1,0]])
    # approximate the last one so that w has the same length as R
    w[-1] = w[-2] + (w[-2] - w[-3])
    return w


###############################################################################
###############################################################################
##########                 Functions for External Use               ###########
###############################################################################
###############################################################################

def set_integration_params(a0, lam, mu, gammadot, Gamma):
    """ Set dimensionless parameters tau and capillary number and use them to
    determine an appropriate time-frame and step size to integrate dgdt.

    Inputs:
        a0          np.array(3), initial axes lengths
        lam         float, viscosity ratio
        mu          float, matrix viscosity
        gammadot    float, shear rate
        Gamma       float, interfacial tension

    Returns:
        t0          float, initial time (hard-coded to 0 right now)
        t1          float, end time
        dt          float, time-step
        tau         dimensionless "time"
        cap         capillary number
    
    """
    rad0 = np.prod(a0)**(1./3)                                                  
    tau = lam * mu * rad0 / Gamma                                               
    cap = lam * mu * gammadot * rad0 / Gamma                                    
                                                                              
    t0 = 0                                                                      
    t1 =  200 * lam / cap *  (4e-9 / Gamma)                                     
    dt = (t1 - t0) / 1000.
    return [t0,t1,dt,tau,cap]


def evolve(t0, t1, dt, a0, lam, mu, gammadot, Gamma, JustAngles = True):
    """ A simple wrapper for integrate_dgdt that gives rotations, axes,
    and the angular velocity. 

    Inputs:
        t0            float, start time
        t1            float, stop time (this may be changed)
        dt            float, time step
        a0            initial axes lengths
        lam           float. Viscosity ratio. 
        mu            float. "Matrix" (external fluid) viscosity. 
        gammadot      float. shear rate.
        Gamma         float. Interfacial tension. 
        JustAngles    bool, sets whether we return full rotations or just
                      the first two entries

    Outputs:
        axes          np.array[(M,3)] of axes lengths
        R             if JustAngles == True, np.array[(M,3,3)] of rotation matrices.
                      if JustAngles == False, np.array[(M,2)] where we have
                      replaced R with [R[0,0], R[1,0]]. This is to save space.
        w             np.array[(N,3)], angular velocity at each time. 
        T             np.array(M) of times. 
    """ 
    Y,T = integrate_dgdt(t0,t1,dt,a0,lam,mu,gammadot,Gamma)
    axes, R = shapetensors_to_axes_rots(Y)
    w = angular_velocity(R, dt)
    if JustAngles == True:
        rotAngles = np.zeros([np.shape(T)[0],2])
        rotAngles[:,0] = R[:,0,0]
        rotAngles[:,1] = R[:,1,0]
        R = rotAngles
    return([axes, R, w, T])

def evolve_solid(t0, t1, dt, a0, gammadot, JustAngles = True):
    """ Compute the motion of a solid ellipsoid (prescribed by Blaser instead
    of JWT model). Much easier!

    Inputs:
        t0            float, start time
        t1            float, stop time (this may be changed)
        dt            float, time step
        a0            initial axes lengths
        gammadot      float. shear rate.
        JustAngles    bool, sets whether we return full rotations or just
                      the first two entries

    Outputs:
        axes          np.array[(M,3)] of axes lengths
        R             if JustAngles == True, np.array[(M,3,3)] of rotation matrices.
                      if JustAngles == False, np.array[(M,2)] where we have
                      replaced R with [R[0,0], R[1,0]]. This is to save space.
        w             np.array[(N,3)], angular velocity at each time. 
        T             np.array(M) of times. 

    """
    [a1, a2, a3] = a0
    r = (a1**2 - a2**2) / (a1**2 + a2**2)
    alpha = np.sqrt((1+r)/(1-r))
    beta = np.sqrt(1-r**2) * gammadot / 2
    T = np.arange(t0,t1,dt)
    # compute the angle phi(t) over the grid
    # we want the rotation to be in the clockwise direction because of the 
    # orientation of the floc in the shear field. 
    Period = 2 * np.pi * (a1**2 + a2**2) / ( a1 * a2 * gammadot  )
    phi = -np.arctan( (a2/a1) * np.tan(2 * np.pi * T / Period)  )
    # arctan's range means we get period jumps from -pi/2 to pi/2. We don't 
    # want this because it causes huge spikes in the angular velocity.
    #  We go through and fix it so the angle is always decreasing.
    flag = True
    while flag == True:
        for i in range(len(phi)-1):
            flag = False
            if (phi[i+1]>phi[i]):
                phi[i+1:] += -np.pi
                flag = True
    # set up the rotations
    cp = np.cos(phi)
    sp = np.sin(phi)
    R = np.array( [np.array([[  cp[i], sp[i], 0 ],
                             [ -sp[i], cp[i], 0 ],
                             [  0,     0,     1 ]]).T 
                  for i in range(len(phi)) ] )
    w = angular_velocity(R, dt)
    if JustAngles == True:
        # overwrite the rotation matrix with just angles
        R = np.zeros([np.shape(T)[0],2])
        R[:,0] = cp
        R[:,1] = sp

    # make an array of the axes at each time step (note: axes are constant)
    axes_unsafe = np.reshape(np.repeat(a0,len(phi)), (3,len(phi))).T
    axes = np.ascontiguousarray(axes_unsafe)

    return([axes, R, w, T])

def angles_to_matrix(angles):
    """ Given the short form output of the angles from deform, reconstructs the matrix R. 
    Inputs: 
        angles        np.array(2), corresponding to [R[0,0], R[1,0]]

    Outputs:
        R             np.array([3,3)], reconstructed rotation
    """
    [c,s] = angles
    R = np.array([[ c , -s , 0 ], \
                  [ s ,  c , 0 ], \
                  [ 0 ,  0 , 1 ]])
    return R

