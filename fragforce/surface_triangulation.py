import numpy as np
import pickle as pickle

""" Functions to generate a triangulation of a sphere, and compute the revelant
quantities needed for the simulation. These functions are not optimized, and
are pretty slow. They are not meant to be imported for use into the main
program. Instead run this once to generate and write the data to disk as a 
pickle file. The force functions assume this file exists. 
"""

def triangulate_unit_sphere(iterMax):
    """ Iteratively triangulate a unit sphere, beginning with a tetrahedron. The procedure
    is to take the midpoint of each triangle, pull it out to the surface of the sphere,
    and then create 4 new triangles out of this one. At each step we do this to each triangle:


    . . . . . . . . . . . . . . . .      . . . . . . . . . . . . . . . .
     .                           .        .             .             .
      .                         .          .           . .           .
       .                       .            .         .   .         .
        .                     .              .       .     .       .
         .                   .                .     .       .     .
          .                 .                  .   .         .   .
           .               .                    . . . . . . . . .
            .             .        ===>          .             .
             .           .                        .           .
              .         .                          .         .
               .       .                            .       .
                .     .                              .     .
                 .   .                                .   .
                  . .                                  . .
                   .                                    .


    Inputs:
        iterMax         int, how many iterations to do.
                        CAUTION: there are 4^iterMax triangles generated.
                        using 9 creates about 130 MB of data.

    Outputs:
        tLnew           list of lists, each sublist is a list of 4 np.array(3),
                        which constitute the 4 vertices of a triangle.

    """

    # define the vertices of a tetrahedron centered at the origin
    rtwo = 1/np.sqrt(2.0)
    # we want the vertices to lie on the unit circle
    norm = np.sqrt(3/2.)
    v0 = np.array([ 1,  0, -rtwo])/norm
    v1 = np.array([-1,  0, -rtwo])/norm
    v2 = np.array([ 0,  1,  rtwo])/norm
    v3 = np.array([ 0, -1,  rtwo])/norm
    t0 = [v0, v1, v2]
    t1 = [v0, v1, v3]
    t2 = [v0, v2, v3]
    t3 = [v1, v2, v3]
    tL = [t0, t1, t2, t3]
    tLnew = tL
    iter = 0
    while iter < iterMax:
        tLold = tLnew
        tLnew = []
        for t in tLold:
            a = (t[0] + t[1]) / 2.0
            b = (t[1] + t[2]) / 2.0
            c = (t[2] + t[0]) / 2.0
            a = a / np.linalg.norm(a)
            b = b / np.linalg.norm(b)
            c = c / np.linalg.norm(c)
            tnew0 = [t[0], a, c]
            tnew1 = [a, t[1], b]
            tnew2 = [b, t[2], c]
            tnew3 = [c,a,b]
            tLnew.extend([tnew0,tnew1,tnew2,tnew3])
        iter += 1
    return tLnew

def sphere_triangle_comps(triangleList):
    """ Do some computations on the triangles in the surface triangulation. For computing
    the surface forces we need the (outward) normals and the centers of each triangle.
    Various other quantities of interest can be computed from the cross-product of two edges
    of the triangles (like surface area) and so we also save this.

    Inputs:
        triangleList         output of triangulate_unit_sphere. 
                             list of lists, each sublist is a list of 4 np.array(3),
                             which constitute the 4 vertices of a triangle.

    Outputs:
        centers              np.array([4^n, 3]), centers of the triangles
        normals              np.array([4^n, 3]), normals to the triangles
        edge_crosses         np.array([4^n, 3]), cross-product of the edges
    """

    N = len(triangleList)
    centers = np.zeros([N,3])
    normals = np.zeros([N,3])
    edge_crosses = np.zeros([N,3])
    for i in range(N):
        triangle = triangleList[i]
        v1,v2,v3 = triangle[0],triangle[1],triangle[2]
        centers[i] = 1/3. * (v1+v2+v3)
        e1 = v2-v1
        e2 = v3-v1
        edge_crosses[i] = np.cross(e1,e2)
        normal = np.sign(np.dot(edge_crosses[i],centers[i])) * edge_crosses[i]
        normals[i] = normal / np.linalg.norm(normal)
    return [centers, normals, edge_crosses]

def generate_triangulation(n_recursion, write=False):
    """ Generates and writes to disk a .p file containing the centers, normals,
    and edge cross products for the recursive triangulation of depth n_recursion.
    CAUTION: This can be way slow and take tons of disk space. Using n = 9,
    for example, takes up ~160 MB, since we generate 4^n triangles. 

    Inputs:
        n_recursion        int, how many interations to apply.
        write              bool, whether to write the triangulation to disk

    Outputs:
        centers            np.array([4^n, 3]), centers of the triangles
        normals            np.array([4^n, 3]), normals to the triangles
        edge_crosses       np.array([4^n, 3]), cross-product of the edges
    
        NOTE: If write==True, the outputs are saved to a file 
        "sphere_data%(n_recursion).p" in the directory from which the function is called.
    """

    triangleList = triangulate_unit_sphere(n_recursion)
    centers, normals, edge_crosses = sphere_triangle_comps(triangleList)
    filename = "sphere_data" + str(n_recursion) + ".p"
    if write == True:
        pickle.dump([centers, normals, edge_crosses], open(filename, "wb"))
        return
    else:
        return([centers, normals, edge_crosses])

