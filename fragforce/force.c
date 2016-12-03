
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "force.h"


/* ------------------------------------------------------------------------*/
/*                      DE Transform Quadrature                            */
/* ------------------------------------------------------------------------*/

/*  intdeiini and intdei are modified implementations of the double-exponential 
    transform quadrature method from this page:

        http://www.kurims.kyoto-u.ac.jp/~ooura/intde.html
  
    The licensing on these functions, which can be found in the readme.txt of
    the package on the website above, is:
        
    "copyright
         Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
         You may use, copy, modify this code for any purpose and 
         without fee. You may distribute this ORIGINAL package."

*/

void intdeiini(int lenaw,
               double tiny,
               double eps,
               double *aw)
{
    /* ---- adjustable parameter ---- */
    double efs = 0.1, hoff = 11.0;
    /* ------------------------------ */
    int noff, nk, k, j;
    double pi4, tinyln, epsln, h0, ehp, ehm, h, t, ep, em, xp, xm, 
        wp, wm;
    
    pi4 = atan(1.0);
    tinyln = -log(tiny);
    epsln = 1 - log(efs * eps);
    h0 = hoff / epsln;
    ehp = exp(h0);
    ehm = 1 / ehp;
    aw[2] = eps;
    aw[3] = exp(-ehm * epsln);
    aw[4] = sqrt(efs * eps);
    noff = 5;
    aw[noff] = 1;
    aw[noff + 1] = 4 * h0;
    aw[noff + 2] = 2 * pi4 * h0;
    h = 2;
    nk = 0;
    k = noff + 6;
    do {
        t = h * 0.5;
        do {
            em = exp(h0 * t);
            ep = pi4 * em;
            em = pi4 / em;
            j = k;
            do {
                xp = exp(ep - em);
                xm = 1 / xp;
                wp = xp * ((ep + em) * h0);
                wm = xm * ((ep + em) * h0);
                aw[j] = xm;
                aw[j + 1] = xp;
                aw[j + 2] = xm * (4 * h0);
                aw[j + 3] = xp * (4 * h0);
                aw[j + 4] = wm;
                aw[j + 5] = wp;
                ep *= ehp;
                em *= ehm;
                j += 6;
            } while (ep < tinyln && j <= lenaw - 6);
            t += h;
            k += nk;
        } while (t < 1);
        h *= 0.5;
        if (nk == 0) {
            if (j > lenaw - 12) j -= 6;
            nk = j - noff;
            k += nk;
            aw[1] = nk;
        }
    } while (2 * k - noff - 6 <= lenaw);
    aw[0] = k - 6;
}


void intdei(double (*func)(double, double *),
            double *funcpars,
            double a,
            double *aw,
            double *integral, 
            double *err)
{
    int noff, lenawm, nk, k, j, jtmp, jm, m, klim;
    double epsh, ir, fp, fm, errt, errh, errd, h, iback, irback;
    
    noff = 5;
    lenawm = (int) (aw[0] + 0.5);
    nk = (int) (aw[1] + 0.5);
    epsh = aw[4];
    *integral = (*func)(a + aw[noff], funcpars);
    ir = *integral * aw[noff + 1];
    *integral *= aw[noff + 2];
    *err = fabs(*integral);
    k = nk + noff;
    j = noff;
    do {
        j += 6;
        fm = (*func)(a + aw[j], funcpars);
        fp = (*func)(a + aw[j + 1], funcpars);
        ir += fm * aw[j + 2] + fp * aw[j + 3];
        fm *= aw[j + 4];
        fp *= aw[j + 5];
        *integral += fm + fp;
        *err += fabs(fm) + fabs(fp);
    } while (aw[j] > epsh && j < k);
    errt = *err * aw[3];
    errh = *err * epsh;
    errd = 1 + 2 * errh;
    jtmp = j;
    while (fabs(fm) > errt && j < k) {
        j += 6;
        fm = (*func)(a + aw[j], funcpars);
        ir += fm * aw[j + 2];
        fm *= aw[j + 4];
        *integral += fm;
    }
    jm = j;
    j = jtmp;
    while (fabs(fp) > errt && j < k) {
        j += 6;
        fp = (*func)(a + aw[j + 1], funcpars);
        ir += fp * aw[j + 3];
        fp *= aw[j + 5];
        *integral += fp;
    }
    if (j < jm) jm = j;
    jm -= noff + 6;
    h = 1;
    m = 1;
    klim = k + nk;
    while (errd > errh && klim <= lenawm) {
        iback = *integral;
        irback = ir;
        do {
            jtmp = k + jm;
            for (j = k + 6; j <= jtmp; j += 6) {
                fm = (*func)(a + aw[j], funcpars);
                fp = (*func)(a + aw[j + 1], funcpars);
                ir += fm * aw[j + 2] + fp * aw[j + 3];
                *integral += fm * aw[j + 4] + fp * aw[j + 5];
            }
            k += nk;
            j = jtmp;
            do {
                j += 6;
                fm = (*func)(a + aw[j], funcpars);
                ir += fm * aw[j + 2];
                fm *= aw[j + 4];
                *integral += fm;
            } while (fabs(fm) > errt && j < k);
            j = jtmp;
            do {
                j += 6;
                fp = (*func)(a + aw[j + 1], funcpars);
                ir += fp * aw[j + 3];
                fp *= aw[j + 5];
                *integral += fp;
            } while (fabs(fp) > errt && j < k);
        } while (k < klim);
        errd = h * (fabs(*integral - 2 * iback) + fabs(ir - 2 * irback));
        h *= 0.5;
        m *= 2;
        klim = 2 * klim - noff;
    }
    *integral *= h;
    if (errd > errh) {
        *err = -errd * m;
    } else {
        *err *= aw[2] * m;
    }
}


double hypergeo3(double t,
                 double pV[6])
{
    /* The integrand to a hypergeometric integral, of the form

        f(t) = (t + z1)^(-b1) * (t + z2)^(-b2) * (t + z3)^(-b3).

    Used in the integration routine integrate_hypergeo3. 

    Inputs:
        t                      function variable
        pV                     zV and bV parameters (note the signs on the bi):
                               pV = [z1, z2, z3, b1, b2, b3]

    Returns:
        f(t)                   f as above

    */
    double out = pow(t+pV[0], -pV[3]) * pow(t+pV[1], -pV[4]) * pow(t+pV[2], -pV[5]);
    return out;
}


double integrate_hypergeo3(double pV[6])
{
    /* Integrates the hypergeometric function hypergeo3 using the DE Transform
    quadrature method implemented in intdei. Note that the quadrature function
    is a black box here. 

    Inputs:
        pV                     zV and bV parameters (note the signs on the bi):
                               pV = [z1, z2, z3, b1, b2, b3]

    Returns:
        integration_result     integral from 0 to infinity of hypergeo3 with
                               pars pV.

    */

    double tiny = 1.0e-307;
    int lenaw = 8000;
    double aw[lenaw], integration_result, err;

    intdeiini(lenaw, tiny, 1.0e-10, aw);
    intdei(hypergeo3, pV, 0.0, aw, &integration_result, &err);

    return integration_result;
}

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
    double lam = 0;
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

    double Rd1 = 0;
    double Rd2 = 0;

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


void scale_plane(double a_initial[3],
                 double a_current[3],
                 double pn_initial[3], 
                 double px_initial[3], 
                 double pn_scaled[3],
                 double px_scaled[3])
    /* Rescales the plane normal and interior point given the original axes
    lengths for which the quantities were defined and the current axes lengths.
    This is needed when we specify the plane normals and interior points using
    the minimum spanning tree on the bacterial centers of mass for our
    application, because we assume the bacteria move wrt eachother as the 
    ellipsoid changes shape.

    Inputs:
        a_initial           axes lengths at time = 0
        a_current           axes lengths at current time
        pn_initial             normal to intersecting plane, time = 0
        px_initial             point in intersecting plane, time = 0
        pn_scaled              modified 
        px_scaled              modified

    Modifies:
        pn_scaled              normal to intersecting plane
        px_scaled              point in intersecting plane
        
    */
{
    int i = 0;
    double scale;

    for ( i = 0 ; i < 3 ; i++ ) {
        scale = a_current[i] / a_initial[i];
        pn_scaled[i] = pn_initial[i] * scale;
        px_scaled[i] = px_initial[i] * scale;      
    }

    // set normalization constant 
    double normalization = sqrt( 
        pn_scaled[0] * pn_scaled[0] + 
        pn_scaled[1] * pn_scaled[1] + 
        pn_scaled[2] * pn_scaled[2] );
    for ( i = 0 ; i < 3 ; i++ ) {
        pn_scaled[i] = pn_scaled[i] / normalization;
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
        nfacets                number of facets in the triangulation
        a                      axes lengths
        srf_centers_scaled     modified
        srf_areas_scaled       modified
        srf_normals_scaled     modified
        srf_centers_sph        facet centers, sphere
        srf_crosses_sph        cross-product of the facet edges, sphere
        srf_normals_sph        normals to facets, sphere 

    modifies:
        srf_centers_scaled     facet centers scaled to ellipsoid
        srf_areas_scaled       facet areas scaled to ellipsoid
        srf_normals_scaled     facet normals scaled to ellipsoid
    
    */
{
    double d1 = a[0];
    double d2 = a[1];
    double d3 = a[2];
    int i = 0;

    double pnr0 = 0, pnr1 = 0, pnr2 = 0, sc = 0, ar0 = 0, ar1 = 0, ar2 = 0;
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

    int iter = 0;
    double ww[3] = { 0 };
    for (iter = 0 ; iter < 3 ; iter++) {
        ww[iter] = -w[iter];
    }


    double aa[3] = { a[0]*a[0] , a[1]*a[1] , a[2]*a[2] };
    double Xp[3] = { 0 };
    double Xpp[3] = { 0 };

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

    double E[3][3];
    memset(E, 0, sizeof E);
    double W[3][3];
    memset(W, 0, sizeof W);

    int i,j;
    for ( i = 0 ; i < 3 ; i++ ) {
        for (j = 0 ; j < 3 ; j++ ) {  
            E[i][j] = 0.5 * ( L[i][j] + L[j][i] );
            W[i][j] = 0.5 * ( L[i][j] - L[j][i] );
        }
    }


    double a01n = 0, a02n = 0, a10n = 0, a12n = 0, a20n = 0, a21n = 0;
    a01n = chi[1] * E[0][1] - aa[0] * Xp[2] * (ww[2] - W[0][1]);
    a02n = chi[2] * E[0][2] + aa[0] * Xp[1] * (ww[1] + W[0][2]);
    a10n = chi[0] * E[1][0] + aa[1] * Xp[2] * (ww[2] + W[1][0]);
    a12n = chi[2] * E[1][2] - aa[1] * Xp[0] * (ww[0] - W[1][2]);
    a20n = chi[0] * E[2][0] - aa[2] * Xp[1] * (ww[1] - W[2][0]);
    a21n = chi[1] * E[1][2] + aa[2] * Xp[0] * (ww[0] + W[2][1]);

    double a01d = 0, a02d = 0, a10d = 0, a12d = 0, a20d = 0, a21d = 0;
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
    double k = 0, kt = 0, novera = 0, area = 0;
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
    double nf = 0;
    int i;
    for ( i = 0 ; i < 3 ; i++ ) {
        nf += pn_scaled[i] * fonf[i];
    }
    // if it's 0, exit and return 0
    if (nf == 0.0) return nf;
    // compute pn dot c-px
    double ncx = 0;
    for ( i = 0 ; i < 3 ; i++ ) {
        ncx += pn_scaled[i] * ( srf_center_scaled[i] - px_scaled[i] );
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
    compute the integral of the force density over the surface. Normalizes
    the plane normal. 

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
    // make sure the plane normal is unit length
    double scale = 1 / sqrt( pn_scaled[0] * pn_scaled[0] +
                             pn_scaled[1] * pn_scaled[1] + 
                             pn_scaled[2] * pn_scaled[2] );
    pn_scaled[0] *= scale;
    pn_scaled[1] *= scale;
    pn_scaled[2] *= scale;
    
    // get the area of intersection to check if plane intersects
    double aoi = area_of_intersection(a, pn_scaled, px_scaled);
    if ( aoi <= 0) return 0; //the the plane does not intersect
    double fonf[3] = {0}, c[3] = {0}, nf = 0, total_force = 0;
    int i, j;
    for ( i = 0 ; i < NFacets ; i++ )
    {
        for ( j = 0 ; j < 3 ; j++ )
        {
            fonf[j] = fonfV[i][j];
            c[j] = srf_centers_scaled[i][j];
        }
        nf = correct_pndotf(fonf, c, pn_scaled, px_scaled);
        total_force += nf;
    }
    return total_force;
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
        double mu,
        int scale_planes_bool)

    /* Computes the fragmentation force over the evolution of a droplet, for
    a set of intersecting planes. The intersecting planes can be scaled
    according to the shape of the droplet if desired.  

    Inputs:
        NTimes                 number of time steps in the shape evolution
        NPlanes                number of intersecting planes
        NFacets                number of facets in the triangulation
        fragforceV             modified
        aV                     axes lengths at each time point
        RV                     0,1 and 1,0 entry of rotation matrix at each time
        wV                     angular velocity at each time point
        pnV                    normals to intersecting planes
        pxV                    interior points to intersecting planes
        srf_centers_sph        facet centers, sphere
        srf_crosses_sph        cross-product of the facet edges, sphere
        srf_normals_sph        normals to facets, sphere 
        gammadot               shear rate
        p0                     external pressure 
        mu                     matrix viscosity
        scale_planes_bool      whether or not to scale the planes

    Modifies:
        fragforceV             fragmentation force, [i,j]th entry is force at
                               ith time w.r.t. jth plane. 
    */
{

    double a[3] = {0};
    double c = 0 ,s = 0;
    double w[3] = {0};

    double chi[3] = {0};
    double L[3][3];
    memset(L, 0, sizeof L);
    double A[3][3];
    memset(A, 0, sizeof A);
    double farg[3][3];
    memset(farg, 0, sizeof farg);

    //VLA's can't be initialized without explicitly assigning the values.
    //Since initializing everything to 0 was kind of overkill anyway I 
    //am going to leave these ones for now. 
    double srf_centers_scaled[NFacets][3];
    double srf_areas_scaled[NFacets];
    double srf_normals_scaled[NFacets][3];
    //double fdonfV[NFacets][3];
    double fonfV[NFacets][3];

    double pn_scaled[3] = {0};
    double px_scaled[3] = {0};


    int TimeStep, PlaneStep;

    // initialize the quantities for scaling the planes
    double a_initial[3];
    double pn[3] = {0};
    double px[3] = {0};
    a_initial[0] = aV[0][0];
    a_initial[1] = aV[0][1];
    a_initial[2] = aV[0][2];

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
        
        scale_triangulation(NFacets, 
                            a,
                            srf_centers_scaled,
                            srf_areas_scaled,
                            srf_normals_scaled,
                            srf_centers_sph, 
                            srf_crosses_sph,
                            srf_normals_sph);
        //normal use would be:
        //    set_force_density(NFacets, fdonfV, farg, srf_normals_scaled);
        // but we are going to save space and use fonfV only, first setting
        // fdonfV and then overwriting it with fonfv.
        set_force_density(NFacets, fonfV, farg, srf_normals_scaled);
        set_force_facets(NFacets, fonfV, fonfV, srf_areas_scaled);

        for (PlaneStep = 0 ; PlaneStep < NPlanes ; PlaneStep++ )
        {
            if ( scale_planes_bool ) {
                // then scale the plane quantities
                pn[0] = pnV[PlaneStep][0];
                pn[1] = pnV[PlaneStep][1];
                pn[2] = pnV[PlaneStep][2];
                px[0] = pxV[PlaneStep][0];
                px[1] = pxV[PlaneStep][1];
                px[2] = pxV[PlaneStep][2];
                scale_plane(a_initial, a, pn, px, pn_scaled, px_scaled);
            }
            else {
                // then don't scale them
                pn_scaled[0] = pnV[PlaneStep][0];
                pn_scaled[1] = pnV[PlaneStep][1];
                pn_scaled[2] = pnV[PlaneStep][2];
                px_scaled[0] = pxV[PlaneStep][0];
                px_scaled[1] = pxV[PlaneStep][1];
                px_scaled[2] = pxV[PlaneStep][2];
            }


            fragforceV[TimeStep][PlaneStep] = sum_forces(NFacets, a, fonfV, srf_centers_scaled,
                                                         pn_scaled, px_scaled);
        }
    }
}


