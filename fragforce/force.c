
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ------------------------------------------------------------------------*/
/*                      Elliptic Integral                                  */
/* ------------------------------------------------------------------------*/

double rd_recursion(double x,
                    double y,
                    double z,
                    int n,
                    int N)

    /* Used to compute the elliptic integral in Blaser. This can be done using a
    recursion relation which is implemented here. The recursion is of the form
    R_d(x_n) = R_d(x_{n+1}) + f(x_n), where R_d(x_n) goes to 0 as x_n goes to 
    infinity. Thus we compute the recusion some N times and then call the last 
    one 0. Note that the final integral needs to be scaled by 2/3. This is 
    recursive function. 

    Inputs:
        x,y,z                  parameters of integral
        n                      current step
        N                      maximum steps

    Outputs:
        R_d(x_N)               double, result after N iterations

    */
{
    double lam;
    if (n < N) {
        lam = sqrt(x*y) + sqrt(y*z) + sqrt(x*z);
        return ( 2.0 * rd_recursion( x + lam, y + lam, z + lam, n + 1, N) \
                 + 3. / (sqrt(z) * (z +lam)) );
    }
    else {
        lam = sqrt(x*y) + sqrt(y*z) + sqrt(x*z);
        return( 3.0 / ( sqrt( z ) * ( z + lam ) ) );
    }
}


double rd_converge(double x,
                   double y,
                   double z)
    /* Implements rd_recursion until we reach acceptable relative error.
    This means we compute rd_recursion (which is recursive N times) for some
    N and then again for 2N, and if the relative difference between these
    isn't small enough, we do it again for 4N, and so on, until it converges
    or until we reach the maximum permissible iterations. 

    Inputs:
        x,y,x                  parameters of integral

    Outputs:
        R_d(X_M)               result after M iterations

    */
{
    int N = 20;
    double tol = 10E-10;
    double error = tol + 1;   // Make sure we trigger the first while call.
    int MAX_IT = 8;           // Every iteration doubles the number of recursions.
    int it = 0;               // Count the iterations.

    double Rd1;
    double Rd2;

    while (error > tol && it <= MAX_IT) {
        Rd1 = rd_recursion(x,y,z,0,N);
        Rd2 = rd_recursion(x,y,z,0,2*N);
        error = fabs( (Rd1 - Rd2) / Rd2 );
        N *= 2;
        it += 1;
    }

    // Exit the program with error message if integral does not converge. 
    if (it > MAX_IT) {
        printf("Error: rd_converge failed to converge after %i iterations", MAX_IT * N);
        exit(0);
    }

    return Rd2;
}

void set_chi(double chi[3], 
             double a,
             double b,
             double c)
    /* Computes the three elliptic integrals using the recursion relation. 
    Note the scaling factor.

    Inputs:
        chi                    modified
        a,b,c                  parameters of the integral

    Modifies:
        chi                    three elliptic integrals                                  

    */
{
    double scale = 2.0 / 3.0;
    chi[0] = scale * rd_converge(c*c,b*b,a*a);
    chi[1] = scale * rd_converge(a*a,c*c,b*b);
    chi[2] = scale * rd_converge(b*b,a*a,c*c);
}


/* ------------------------------------------------------------------------*/
/*                     Scaling Functions                                   */
/* ------------------------------------------------------------------------*/


void scale_edge(double axes[3],
                double edge_normal_sph[3], 
                double edge_center_sph[3], 
                double edge_normal_scaled[3],
                double edge_center_scaled[3])
    /* Scale a single edge from the sphere to the axes. Scales the normals and 
    centers. Currently not used and needs proper documentation. */
{
    int i;
    double prenorm_edge_normal[3];

    for ( i = 0 ; i < 3 ; i++ ) {
        prenorm_edge_normal[i] = edge_normal_sph[i] * axes[i];
        edge_center_scaled[i] = edge_center_sph[i] * axes[i];      
    }

    // set normalization constant 
    double normalization = sqrt( 
        prenorm_edge_normal[0] * prenorm_edge_normal[0] + 
        prenorm_edge_normal[1] * prenorm_edge_normal[1] + 
        prenorm_edge_normal[2] * prenorm_edge_normal[2] );
    for ( i = 0 ; i < 3 ; i++ ) {
        edge_normal_scaled[i] = prenorm_edge_normal[i] / normalization;
    }
}


void scale_triangulation(int NFacets,
                         double a[3], 
                         double srf_centers_scaled[NFacets][3],
                         double srf_areas_scaled[NFacets],
                         double srf_normals_scaled[NFacets][3],
                         double srf_centers_sph[NFacets][3], 
                         double srf_crosses_sph[NFacets][3],
                         double srf_normals_sph[NFacets][3])

    /* The unit sphere's surface triangulation is loaded for each floc. At each
    time-step, the sphere triangulation needs to be transformed according to the
    ellipsoid's axes at that time. This scaling is hard-coded for effeciency.
    Note that the quantities ending in _sph are read from disk and passed
    through the python wrapper for this function.
    
    Note that the centers scale like the axes and the normals scale like inverse 
    axes. The areas are computed by taking half the cross-product of the adjusted
    edges. This computation is the most obscure - we take the norm of a dot that's
    all hard coded here. 

    Inputs: 
        NFacets                number of facets in the triangulation
        a                      axes lengths
        srf_centers_scaled     modified
        srf_areas_scaled       modified
        srf_normals_scaled     modified
        srf_centers_sph        facet centers, sphere
        srf_crosses_sph        cross-product of the facet edges, sphere
        srf_normals_sph        normals to facets, sphere 

    Modifies:
        srf_centers_scaled     facet centers scaled to ellipsoid
        srf_areas_scaled       facet areas scaled to ellipsoid
        srf_normals_scaled     facet normals scaled to ellipsoid
    
    */
{
  double d1 = a[0];
  double d2 = a[1];
  double d3 = a[2];
  int i;

  double pnr0, pnr1, pnr2, sc, ar0, ar1, ar2;
  for ( i = 0 ; i < NFacets ; i++ ) {
    srf_centers_scaled[i][0] = d1 * srf_centers_sph[i][0];
    srf_centers_scaled[i][1] = d2 * srf_centers_sph[i][1];
    srf_centers_scaled[i][2] = d3 * srf_centers_sph[i][2];

    ar0 = d2 * d3 * srf_crosses_sph[i][0];
    ar1 = d1 * d3 * srf_crosses_sph[i][1];
    ar2 = d1 * d2 * srf_crosses_sph[i][2];

    srf_areas_scaled[i]  = 0.5 * sqrt(ar0 * ar0 + ar1 * ar1 + ar2 * ar2);
    pnr0 = srf_normals_sph[i][0] / d1;
    pnr1 = srf_normals_sph[i][1] / d2;
    pnr2 = srf_normals_sph[i][2] / d3;

    // Normalize
    sc = 1.0 / sqrt(pnr0*pnr0 + pnr1*pnr1 + pnr2*pnr2);
    srf_normals_scaled[i][0] = sc * pnr0;
    srf_normals_scaled[i][1] = sc * pnr1;
    srf_normals_scaled[i][2] = sc * pnr2;
  }
}



/* ------------------------------------------------------------------------*/
/*                      Force Functions                                    */
/* ------------------------------------------------------------------------*/


void set_L(double L[3][3], double c, double s, double gammadot)
    /* Constructs the velocity gradient L in the ellipsoid frame out of
    the output of the deformation simulation. 

    Inputs:
        L                      modified
        c, s                   0,1 and 1,0 entry of rotation matrix
        gammadot               shear rate

    Modifies:
        L                      velocity gradient in the ellipsoid frame. 
                               Equivalent to R^T L0 R where R is the appropriate
                               rotation and L0 is the velocity gradient in the
                               lab frame. L0 must be of the form 

                               0 gammadot 0
                               0    0     0
                               0    0     0 
    */
{
    
    L[0][2] = 0.0;
    L[1][2] = 0.0;
    L[2][0] = 0.0;
    L[2][1] = 0.0;
    L[2][2] = 0.0;

    L[0][0] = -c * s * gammadot;
    L[0][1] =  c * c * gammadot;
    L[1][0] = -s * s * gammadot;
    L[1][1] = -L[0][0];
}


void set_A(double A[3][3], 
           double a[3],
           double w[3],
           double L[3][3],
           double chi[3])
    /* Sets the matrix A defined in equation 18 in Blaser. Entries are hard-
    coded. Notice that in each case the sums collapse into a single term; I 
    worked these out by hand and coded them in here. Tests check out with the
    explicit summation representation coded in Mathematica. 

    Inputs:
        A                      modified
        a                      axes lengths
        w                      angular velocity
        L                      velocity gradient, ellipsoid frame
        chi                    elliptic integrals

    Modifies:
        A                      matrix A in equation 18 of Blaser

  */

{

    // We need to change the sign on the angular velocity. I can't find where I
    // made the error; maybe there's a typo in Blaser. The reason I know we need
    // to change this is the force arrow plots - they look wrong if we don't and
    // right if we do. For transparency we use a new vector ww = - w

    int iter;
    double ww[3];
    for (iter = 0 ; iter < 3 ; iter++) {
        ww[iter] = -w[iter];
    }


    double aa[3] = { a[0]*a[0] , a[1]*a[1] , a[2]*a[2] };
    double Xp[3];
    double Xpp[3];

    Xp[0] = (chi[2] - chi[1]) / (aa[1] - aa[2]);
    Xp[1] = (chi[0] - chi[2]) / (aa[2] - aa[0]);
    Xp[2] = (chi[1] - chi[0]) / (aa[0] - aa[1]);

    Xpp[0] = ( aa[1] * chi[1] - aa[2] * chi[2] ) / (aa[1] - aa[2]);
    Xpp[1] = ( aa[2] * chi[2] - aa[0] * chi[0] ) / (aa[2] - aa[0]);
    Xpp[2] = ( aa[0] * chi[0] - aa[1] * chi[1] ) / (aa[0] - aa[1]);

    // So this is great: aa[0], aa[1], and aa[2] can be the same, which gives
    // a division by 0. When this happens you get a nan for the Xp and Xpp
    // above and then no error shows up later so you have no clue what
    // happened. I tried some limits a[2] -> a[1] to see what Xpp[0] is, and I
    // get that it should be X[1] which is equal to X[2]. I'm just going to
    // assume that when aa[1] == aa[2] then X[2] == X[1] and so I can just
    // set Xp[0] to 0 and Xpp[0] to X[2]. Note that due to the specifics of
    // my application of this code, I only ever get repeated axes in the 
    // first and second (and never zeroth) indices, so I've only fixed these
    // here. There's more of this hacky garbage later on when we set the
    // entries to A. 

    if ( aa[1] == aa[2] ) {
        Xp[0] = 0;
        Xpp[0] = chi[1];
    }

    double E[3][3], W[3][3];

    int i,j;
    for ( i = 0 ; i < 3 ; i++ ) {
        for (j = 0 ; j < 3 ; j++ ) {  
            E[i][j] = 0.5 * ( L[i][j] + L[j][i] );
            W[i][j] = 0.5 * ( L[i][j] - L[j][i] );
        }
    }


    double a01n, a02n, a10n, a12n, a20n, a21n;
    a01n = chi[1] * E[0][1] - aa[0] * Xp[2] * (ww[2] - W[0][1]);
    a02n = chi[2] * E[0][2] + aa[0] * Xp[1] * (ww[1] + W[0][2]);
    a10n = chi[0] * E[1][0] + aa[1] * Xp[2] * (ww[2] + W[1][0]);
    a12n = chi[2] * E[1][2] - aa[1] * Xp[0] * (ww[0] - W[1][2]);
    a20n = chi[0] * E[2][0] - aa[2] * Xp[1] * (ww[1] - W[2][0]);
    a21n = chi[1] * E[1][2] + aa[2] * Xp[0] * (ww[0] + W[2][1]);

    double a01d, a02d, a10d, a12d, a20d, a21d;
    a01d = 2 * (aa[0] * chi[0] + aa[1] * chi[1]) * Xp[2];
    a10d = a01d;
    a02d = 2 * (aa[0] * chi[0] + aa[2] * chi[2]) * Xp[1];
    a20d = a02d;
    a12d = 2 * (aa[1] * chi[1] + aa[2] * chi[2]) * Xp[0];
    a21d = a12d;

    double dd = 6*(Xpp[0]*Xpp[1] + Xpp[0]*Xpp[2] + Xpp[1]*Xpp[2]);
    double a00n = 2 * Xpp[0] * E[0][0] - Xpp[1] * E[1][1] - Xpp[2] * E[2][2];
    double a11n = 2 * Xpp[1] * E[1][1] - Xpp[0] * E[0][0] - Xpp[2] * E[2][2];
    double a22n = 2 * Xpp[2] * E[2][2] - Xpp[1] * E[1][1] - Xpp[0] * E[0][0];

    A[0][0] = a00n / dd;
    A[1][1] = a11n / dd;
    A[2][2] = a22n / dd;
    A[0][1] = a01n / a01d;
    A[0][2] = a02n / a02d;
    A[1][0] = a10n / a10d;
    A[1][2] = a12n / a12d;
    A[2][0] = a20n / a20d;
    A[2][1] = a21n / a21d;
     
    //if a[12] == a[1] then we STILL get some horrible divisions by 0 in here.
    //I think that when this happens the A entry should be 0. Maybe check this
    //later. 

    for ( i = 0 ; i < 3 ; i++ ) {
        for ( j = 0 ; j < 3 ; j++ ) {
            if ( isnan( A[i][j] ) ) {
                A[i][j] = 0;
            }
        }
    }

}


void set_farg(double farg[3][3],
              double a[3],
              double w[3],
              double L[3][3],
              double A[3][3],
              double chi[3],
              double p0,
              double mu)
    /* Computes the matrix in equation 17 in blaser (the result of adding 
    everything in those parentheses together).

    Inputs:
        farg                   modified
        a                      axes at the current time
        w                      angular velocity at the current time
        L                      velocity gradient, ellipsoid frame
        p0                     external pressure 
        mu                     matrix viscosity

    Modifies:
        farg                   argument to the force function

    */
{ 
    double c = ( 8.0 * mu ) / ( a[0] * a[1] * a[2] );
    double diag = - p0 - 4.0 * mu * ( chi[0] * A[0][0] + chi[1] * A[1][1] + \
                                      chi[2] * A[2][2]);
    int i,j;
    for ( i = 0 ; i < 3 ; i++ ) {
        for ( j = 0 ; j < 3 ; j++ ) {
            farg[i][j] = c * A[j][i];
        }
        farg[i][i] += diag;
    }
}


void set_force_density(int NFacets, 
                       double fdonfV[NFacets][3],
                       double farg[3][3],
                       double srf_normals_scaled[NFacets][3])
    /* Computes the force density on the triangulated surface. 
    All this does is take farg dot srf_normals_scaled. 

    Inputs:
        NFacets                number of facets in the triangulation
        fdonfV                 modified
        farg                   argument to the force function.
        srf_normals_scaled     facet normals scaled to ellipsoid

    Modifies:
        fdonfV                 force density on triangulation facets 
    */
{
    int i;
    for ( i = 0 ; i < NFacets ; i++ ) {    
    fdonfV[i][0] = (srf_normals_scaled[i][0] * farg[0][0]  + \
                    srf_normals_scaled[i][1] * farg[0][1]  + \
                    srf_normals_scaled[i][2] * farg[0][2]  ); 

    fdonfV[i][1] = (srf_normals_scaled[i][0] * farg[1][0]  + \
                    srf_normals_scaled[i][1] * farg[1][1]  + \
                    srf_normals_scaled[i][2] * farg[1][2]  ); 

    fdonfV[i][2] = (srf_normals_scaled[i][0] * farg[2][0]  + \
                    srf_normals_scaled[i][1] * farg[2][1]  + \
                    srf_normals_scaled[i][2] * farg[2][2]  ); 
    }       
}


void set_force_facets(int NFacets, 
                      double fonfV[NFacets][3], 
                      double fdonfV[NFacets][3], 
                      double srf_areas_scaled[NFacets])

    /* Computes the force on the triangulated surface. Just
    takes fdonfV and and scales it by the area of the facet.
    
    Inputs:
        NFacets                number of facets in the triangulation
        fonfV                  modified
        fdonfV                 force density on triangulation facets
        srf_areas_scaled       facet areas scaled to ellipsoid

    Modifies:
        fonfV                  force on triangulation facets
    */

{
    int i;
    for ( i = 0 ; i < NFacets ; i++ ) {    
        fonfV[i][0] = fdonfV[i][0] * srf_areas_scaled[i];
        fonfV[i][1] = fdonfV[i][1] * srf_areas_scaled[i];
        fonfV[i][2] = fdonfV[i][2] * srf_areas_scaled[i];
    }       
}


double area_of_intersection(double a[3],
                            double pn_scaled[3],
                            double px_scaled[3])
    /* Compute the area of intersection of a plane defined by normal pn
    and interior point px with an ellipsoid defined by axes lengths a.
    Note that the plane normal and interior point must be in the ellipsoid
    frame. 

    Inputs:
        a                      axes lengths
        pn_scaled              normal to intersecting plane
        px_scaled              point in intersecting plane

    Outpus:
        area                   area of intersection 
    */
{
    double k, kt, novera, area;
    k = fabs(pn_scaled[0] * px_scaled[0] + \
             pn_scaled[1] * px_scaled[1] + \
             pn_scaled[2] * px_scaled[2]);
    // kt = sqrt(dot(a^2,pn_scaled^2))
    kt = sqrt(a[0] * a[0] * pn_scaled[0] * pn_scaled[0] + \
              a[1] * a[1] * pn_scaled[1] * pn_scaled[1] + \
              a[2] * a[2] * pn_scaled[2] * pn_scaled[2]);
    novera = ( pn_scaled[0] / a[0] ) * ( pn_scaled[0] / a[0] ) + \
             ( pn_scaled[1] / a[1] ) * ( pn_scaled[1] / a[1] ) + \
             ( pn_scaled[2] / a[2] ) * ( pn_scaled[2] / a[2] );
    // check if the plane intersects the ellipsoid at all
    if (novera >= kt*kt*kt*kt / (kt * kt) ) {
        area = 0;
    }
    else {
        double pi = 3.14159265358979323846;
        area = pi * ( 1 - k*k / (kt*kt) ) * (a[0]*a[1]*a[2] ) / kt;
    }

    return(area);
}


double correct_pndotf(double fonf[3],
                      double srf_center_scaled[3],
                      double pn_scaled[3],
                      double px_scaled[3])
    /* Compute the magnitude of the component of the force vector fonf acting
    against the plane specified by pn_scaled, px_scaled. The plane must be
    specified in the ellipsoid frame. 

    Inputs:
        fonf                   force on facet
        srf_center_scaled      facet center scaled to ellipsoid
        pn_scaled              normal to intersecting plane
        px_scaled              point in intersecting plane

    Outputs:
        nf                     double, magnitude of fonf acting
                               against the plane  
    */
{
    //Compute the initial dot product
    double nf;
    int i;
    for ( i = 0 ; i < 3 ; i++ ) {
        nf = nf + pn_scaled[i] * fonf[i];
    }
    // if it's 0, exit and return 0
    if (nf == 0.0) return nf;
    // compute pn dot c-px
    double ncx;
    for ( i = 0 ; i < 3 ; i++ ) {
        ncx = ncx + pn_scaled[i] * ( srf_center_scaled[i] - px_scaled[i] );
    }
    double t = - ncx / nf;
    // treat the cases for t
    if ( t<0.0 || t>1.0) {
        // then f does not intersect the plane
        double quant = ( fabs(nf + ncx) - fabs(ncx) );
        double sign = 1;
        if (quant < 0) sign = -1.0;
        nf = fabs(nf) * sign;
        // right now ignore == 0 case. Not sure what to do with this
    }
    else if (t != 0.0) nf = -fabs(nf); // then f intersects the plane 
    else nf = fabs(nf); // then the force originates on the plane
    return nf;
}

double sum_forces(int NFacets,
                  double a[3], 
                  double fonfV[NFacets][3],
                  double srf_centers_scaled[NFacets][3],
                  double pn_scaled[3], 
                  double px_scaled[3])

    /* Sum the components of fonfV acting against the plane; in other words,
    compute the integral of the force density over the surface.

    Inputs
        NFacets                number of facets in the triangulation
        a                      axes lengths
        fonfV                  force on triangulation facets
        srf_centers_scaled     facet centers scaled to ellipsoid
        pn_scaled              normal to intersecting plane
        px_scaled              point in intersecting plane

    Outputs
        total_force            surface force (magnitude) acting against plane
    */
{
    double aoi = area_of_intersection(a, pn_scaled, px_scaled);
    if ( aoi <= 0) return 0; //the the plane does not intersect
    double fonf[3], c[3], nf, total_force;
    int i, j;
    for ( i = 0 ; i < NFacets ; i++ )
    {
        for ( j = 0 ; j < 3 ; j++ )
        {
            fonf[j] = fonfV[i][j];
            c[j] = srf_centers_scaled[i][j];
        }
        nf = correct_pndotf(fonf, c, pn_scaled, px_scaled);
        total_force = total_force + nf;
    }
    return total_force;
    // return (total_force * 2.0) / 2.0;
    // No idea why the hell it is happening but returning total_force gives
    // 0 - somehow doing something to it other than multiplying by 1 gets the
    // correct answer...
}


void frag_force(
        int NTimes, 
        int NPlanes, 
        int NFacets,
        double fragforceV[NTimes][NPlanes],
        double aV[NTimes][3], 
        double RV[NTimes][2],
        double wV[NTimes][3],
        double pnV[NPlanes][3], 
        double pxV[NPlanes][3],
        double srf_centers_sph[NFacets][3],
        double srf_crosses_sph[NFacets][3],
        double srf_normals_sph[NFacets][3],
        double gammadot,
        double p0,
        double mu)
{

    double a[3];
    double c,s;
    double w[3];

    double chi[3];
    double L[3][3];
    double A[3][3];
    double farg[3][3];

    double srf_centers_scaled[NFacets][3];
    double srf_areas_scaled[NFacets];
    double srf_normals_scaled[NFacets][3];
    double fdonfV[NFacets][3];
    double fonfV[NFacets][3];

    double pn[3];
    double px[3];


    int TimeStep, PlaneStep;

    for ( TimeStep = 0 ; TimeStep < NTimes ; TimeStep++ )
    {
        // Assign the time-dependent variables
        c = RV[TimeStep][0];
        s = RV[TimeStep][1]; 
        a[0] = aV[TimeStep][0];
        a[1] = aV[TimeStep][1];
        a[2] = aV[TimeStep][2];
        w[0] = wV[TimeStep][0];
        w[1] = wV[TimeStep][1];
        w[2] = wV[TimeStep][2];
        set_L(L, c, s, gammadot);
        set_chi(chi, a[0], a[1], a[2]);
        set_A(A, a, w, L, chi);
        set_farg(farg, a, w, L, A, chi, p0, mu);
        
        /*
        int nt;
        printf("farg (c) = \n");
        for ( nt = 0 ; nt < 3 ; nt++ )
        {
        printf("%e %e %e\n", farg[nt][0], farg[nt][1], farg[nt][2]);
        }
        */ 
        scale_triangulation(NFacets, 
                            a,
                            srf_centers_scaled,
                            srf_areas_scaled,
                            srf_normals_scaled,
                            srf_centers_sph, 
                            srf_crosses_sph,
                            srf_normals_sph);
        set_force_density(NFacets, fdonfV, farg, srf_normals_scaled);
        set_force_facets(NFacets, fonfV, fdonfV, srf_areas_scaled);
        //printf("%e %e %e \n", fonfV[17][0], fonfV[17][1], fonfV[17][2]);

        for (PlaneStep = 0 ; PlaneStep < NPlanes ; PlaneStep++ )
        {
            pn[0] = pnV[PlaneStep][0];
            pn[1] = pnV[PlaneStep][1];
            pn[2] = pnV[PlaneStep][2];
            px[0] = pxV[PlaneStep][0];
            px[1] = pxV[PlaneStep][1];
            px[2] = pxV[PlaneStep][2];
            //printf("%e %e %e \n", pn[0], pn[1], pn[2]);
            //printf("%e %e %e \n", px[0], px[1], px[2]);
            
            fragforceV[TimeStep][PlaneStep] = sum_forces(NFacets, a, fonfV, srf_centers_scaled, pn, px);
        }
    }
}


