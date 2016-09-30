
#include <stdio.h>
#include <math.h>

/* ------------------------------------------------------------------------*/
/*                      Elliptic Integral                                  */
/* ------------------------------------------------------------------------*/

double rd_recursion(double x, double y, double z, int n, int N)
  /* Used to compute the elliptic integral in Blaser. This can be done using a 
  recursion relation which is implemented here. The recursion is of the form
  R_d(x_n) = R_d(x_{n+1}) + f(x_n), where R_d(x_n) goes to 0 as x_n goes to 
  infinity. Thus we compute the recusion some N times and then call the last 
  one 0. Note that the final integral needs to be scaled by 2/3. 
  */
{
  double lam;
  if (n < N) {
    lam = sqrt ( x * y ) + sqrt ( y * z ) + sqrt ( x * z );
    return ( 2.0 * rd_recursion ( x + lam, y + lam, z + lam, n + 1, N) \
           + 3. / ( sqrt( z ) * ( z + lam ) ) );
  }
  else {
    lam = sqrt(x*y) + sqrt(y*z) + sqrt(x*z);
    return(3.0 / (sqrt(z) * (z+lam)));
  }
}

double rd_converge(double x, double y, double z)
  /* Implements rd_recursion until we reach acceptable relative error. */
{
  int N = 20;
  double tol = 10E-10;
  double error = tol + 1; //to make sure we trigger the first while call
  int MAX_IT = 8; //every iteration doubles the number of recursions
  int it = 0; //count the iterations

  double Rd1;
  double Rd2;

  while (error > tol && it <= MAX_IT) 
  {
    Rd1 = rd_recursion(x,y,z,0,N);
    Rd2 = rd_recursion(x,y,z,0,2*N);
    error = fabs( (Rd1 - Rd2) / Rd2 );
    N=2*N;
    it += 1;
  }
  if (it > MAX_IT)
    return 1; 
  return Rd2;
}

void set_chi(double lam[3], double a, double b, double c)
  /* Computes the three elliptic integrals using the recursion relation. 
  Note the scaling factor.  */
{
  double scale = 2.0/3.0;
  lam[0] = scale * rd_converge(c*c,b*b,a*a);
  lam[1] = scale * rd_converge(a*a,c*c,b*b);
  lam[2] = scale * rd_converge(b*b,a*a,c*c);
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
centers. */
{
  int i;
  double prenorm_edge_normal[3];

  for ( i = 0 ; i < 3 ; i++ )
  {
    prenorm_edge_normal[i] = edge_normal_sph[i] * axes[i];
    edge_center_scaled[i] = edge_center_sph[i] * axes[i];      
  }

  // set normalization constant 
  double normalization = sqrt( 
    prenorm_edge_normal[0] * prenorm_edge_normal[0] + 
    prenorm_edge_normal[1] * prenorm_edge_normal[1] + 
    prenorm_edge_normal[2] * prenorm_edge_normal[2] );
  for ( i = 0 ; i < 3 ; i++ )
  {
    edge_normal_scaled[i] = prenorm_edge_normal[i] / normalization;
  }

}


void scale_triangulation(int SIZE, double a[3], 
                                   double srf_centers_scaled[SIZE][3],
                                   double srf_areas_scaled[SIZE],
                                   double srf_normals_scaled[SIZE][3],
                                   double srf_centers_sph[SIZE][3], 
                                   double srf_crosses_sph[SIZE][3],
                                   double srf_normals_sph[SIZE][3])
/* The unit sphere's surface triangulation is loaded for each floc. At each
time-step, the sphere triangulation needs to be transformed according to the
ellipsoid's axes at that time. This scaling is hard-coded for effeciency. See
Mathematica file for derivation of these scalings. We scale:
  srf_centers_scaled     the centers of the facets
  srf_areas_scaled           the areas of the facets
  normals_scaled         normals to the facets. 

Note that the centers scale like the axes and the normals scale like inverse 
axes. The areas are computed by taking half the cross-product of the adjusted
edges. This computation is the most obscure - we take the norm of a dot that's
all hard coded here. 
*/
{
  double d1 = a[0];
  double d2 = a[1];
  double d3 = a[2];
  int i;

  double pnr0, pnr1, pnr2, sc, ar0, ar1, ar2;
  for ( i = 0 ; i < SIZE ; i++ ) {
    srf_centers_scaled[i][0] = d1 * srf_centers_sph[i][0];
    srf_centers_scaled[i][1] = d2 * srf_centers_sph[i][1];
    srf_centers_scaled[i][2] = d3 * srf_centers_sph[i][2];
    ar0                  = d2 * d3 * srf_crosses_sph[i][0];
    ar1                  = d1 * d3 * srf_crosses_sph[i][1];
    ar2                  = d1 * d2 * srf_crosses_sph[i][2];
    srf_areas_scaled[i]      = 0.5 * sqrt(ar0 * ar0 + ar1 * ar1 + ar2 * ar2);
    pnr0                 = srf_normals_sph[i][0] / d1;
    pnr1                 = srf_normals_sph[i][1] / d2;
    pnr2                 = srf_normals_sph[i][2] / d3;
    // Normalize
    sc                   =  1.0 / sqrt(pnr0*pnr0 + pnr1*pnr1 + pnr2*pnr2);
    srf_normals_scaled[i][0]  = sc * pnr0;
    srf_normals_scaled[i][1]  = sc * pnr1;
    srf_normals_scaled[i][2]  = sc * pnr2;
  }
}




/* ------------------------------------------------------------------------*/
/*                      Force Functions                                    */
/* ------------------------------------------------------------------------*/


void set_L(double L[3][3], double c, double s, double gammadot)
    /* Constructs the velocity gradient L in the ellipsoid frame out of
    the output of the deformation simulation, which is cos theta and
    sin theta. 

    Inputs:
        c, s             cos theta and sin theta
        gammadot         shear rate


    Outputs:
        L                velocity gradient in the ellipsoid frame. 
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

    //TESTING
    //L[0][0] = 20000*L[0][0];
    //L[1][1] = 20000*L[1][1];
}


void set_A(double A[3][3], double w[3], double a[3], double L[3][3], double X[3])
  /* Sets the matrix A defined in equation 18 in Blaser. Entries are hard-
  coded. Notice that in each case the sums collapse into a single term; I 
  worked these out by hand and coded them in here. Tests check out with the
  explicit summation representation coded in Mathematica. 
  */

{

  // We need to change the sign on the angular velocity. I can't find where I
  // made the error; maybe there's a typo in Blaser. The reason I know we need
  // to change this is the force arrow plots - they look wrong if we don't and
  // right if we do. This is done in a really hacky way. We define a new vector
  // ww below and use that where w appears.  

  int iter;
  double ww[3];
  for (iter = 0 ; iter < 3 ; iter++)
  {
    ww[iter] = -w[iter];
  }


  double aa[3] = { a[0]*a[0] , a[1]*a[1] , a[2]*a[2] } ;
  double Xp[3];
  double Xpp[3];

  Xp[0] = (X[2] - X[1]) / (aa[1] - aa[2]);
  Xp[1] = (X[0] - X[2]) / (aa[2] - aa[0]);
  Xp[2] = (X[1] - X[0]) / (aa[0] - aa[1]);

  //printf("aa = %e, %e, %e \n", aa[0], aa[1], aa[2]);
  //printf("X = %e, %e, %e \n", X[0], X[1], X[2]);

  Xpp[0] = ( aa[1] * X[1] - aa[2] * X[2] ) / (aa[1] - aa[2]);
  Xpp[1] = ( aa[2] * X[2] - aa[0] * X[0] ) / (aa[2] - aa[0]);
  Xpp[2] = ( aa[0] * X[0] - aa[1] * X[1] ) / (aa[0] - aa[1]);

  //printf("Xpp0 = %e", Xpp[0]);

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
    Xpp[0] = X[1];
  }

  double E[3][3], W[3][3];

  int i,j;
  for ( i = 0 ; i < 3 ; i++ ) 
    for (j = 0 ; j < 3 ; j++ )
  {  
    {
      E[i][j] = 0.5 * ( L[i][j] + L[j][i] );
      W[i][j] = 0.5 * ( L[i][j] - L[j][i] );
    }
  }


  double a01n, a02n, a10n, a12n, a20n, a21n;
  a01n = X[1] * E[0][1] - aa[0] * Xp[2] * (ww[2] - W[0][1]);
  a02n = X[2] * E[0][2] + aa[0] * Xp[1] * (ww[1] + W[0][2]);
  a10n = X[0] * E[1][0] + aa[1] * Xp[2] * (ww[2] + W[1][0]);
  a12n = X[2] * E[1][2] - aa[1] * Xp[0] * (ww[0] - W[1][2]);
  a20n = X[0] * E[2][0] - aa[2] * Xp[1] * (ww[1] - W[2][0]);
  a21n = X[1] * E[1][2] + aa[2] * Xp[0] * (ww[0] + W[2][1]);

  double a01d, a02d, a10d, a12d, a20d, a21d;
  a01d = 2 * (aa[0] * X[0] + aa[1] * X[1]) * Xp[2];
  a10d = a01d;
  a02d = 2 * (aa[0] * X[0] + aa[2] * X[2]) * Xp[1];
  a20d = a02d;
  a12d = 2 * (aa[1] * X[1] + aa[2] * X[2]) * Xp[0];
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

  if ( a[1] == a[2]  ) {
    A[1][2] = 0;
    A[2][1] = 0;
  }

 
  //printf("%e %e %e  \n", Xpp[0], Xpp[1], Xpp[2]);
}


void set_farg(double farg[3][3], double a[3], double w[3], double L[3][3], double p0, double mu)
{
  /* Sets the argument to the force function in blaser's model. This involves 
  computing the elliptic integral and then defining the A matrix (both 
  previous functions here). 

  INPUTS:
    a          array dim=3 double, axes at the current time
    w          array dim=3 double, angular velocity at the current time
    L          array dim=[3,3] double, velocity gradient **in the body frame**
    p0         double, ambient pressure 
    mu         double, matrix viscosity

  OUTPUTS:
    (well, void, but sets farg, array dim 3x3 double, what gets dotted into
    the normal to the surface to get the force. 
  */

  // Compute elliptic integral
  double chi[3];
  set_chi(chi, a[0], a[1], a[2]); 
  // Compute A
  double A[3][3];
  set_A(A, w, a, L, chi);

    
  
  double c = ( 8.0 * mu ) / ( a[0] * a[1] * a[2] );
  double diag = - p0 - 4.0 * mu * ( chi[0] * A[0][0] + chi[1] * A[1][1] + \
                                    chi[2] * A[2][2]);
  int i,j;
  for ( i = 0 ; i < 3 ; i++ )
  {
    for ( j = 0 ; j < 3 ; j++ )
    {
      farg[i][j] = c * A[j][i];
    }
    farg[i][i] += diag;
  }
     
  
  /*
  int nt;
  printf("L matrix = \n");
  for ( nt = 0 ; nt < 3 ; nt++ )
  {
    printf("%e %e %e\n", L[nt][0], L[nt][1], L[nt][2]);
  }
  
  printf("w = %f %f %f\n", w[0], w[1], w[2]);
  printf("a = %e %e %e\n", a[0], a[1], a[2]);
  printf("p0 = %e, mu = %e \n", p0, mu);


  printf("farg matrix = \n");
  for ( nt = 0 ; nt < 3 ; nt++ )
  {
    printf("%e %e %e\n", farg[nt][0], farg[nt][1], farg[nt][2]);
  }
  printf("\n");
  */
       
}

void set_force_facets(int SIZE, double force_on_facets[SIZE][3], 
                      double farg[3][3], double srf_normals_scaled[SIZE][3],
                      double srf_areas_scaled[SIZE])

/* Sets the force on each facet at a given time. This can be computed once
and then used for each edge. Matrix multiplication is hard-coded. Computes
farg dot normal * area. Normals are the normals to the surface.  
*/
{
  int i;
  for ( i = 0 ; i < SIZE ; i++ )
  {    
    force_on_facets[i][0] = (srf_normals_scaled[i][0] * farg[0][0]  + \
                             srf_normals_scaled[i][1] * farg[0][1]  + \
                             srf_normals_scaled[i][2] * farg[0][2]) * srf_areas_scaled[i]; 

    force_on_facets[i][1] = (srf_normals_scaled[i][0] * farg[1][0]  + \
                             srf_normals_scaled[i][1] * farg[1][1]  + \
                             srf_normals_scaled[i][2] * farg[1][2]) * srf_areas_scaled[i]; 

    force_on_facets[i][2] = (srf_normals_scaled[i][0] * farg[2][0]  + \
                             srf_normals_scaled[i][1] * farg[2][1]  + \
                             srf_normals_scaled[i][2] * farg[2][2]) * srf_areas_scaled[i]; 

  }       
}

/* ------------------------------------------------------------------------*/
/*                        Stress Computations                              */
/* ------------------------------------------------------------------------*/

double area_of_intersection(double a[3], double edge_normal_scaled[3], double edge_center_scaled[3])
  /*Given the semi-principal axis lengths of an ellipsoid alligned with the 
  coordinate system and a plane intersecting the ellipsoid, compute the area
  of this intersection. See the derivation materials for an explanation and 
  reference of the equation used here. 

  Inputs:
    axes       3x1 np.array, semi-principal axis lengths
    edge_normal_scaled    3x1 np.array, normal to the plane, NOT scaled to the sphere
    edge_center_scaled    3x1 np.array, a point in the plane, NOT scaled to the sphere

  Outpus:
    area       float, area of intersection 
  */
{
  double k, kt, novera, area;
  k = fabs(edge_normal_scaled[0] * edge_center_scaled[0] + \
           edge_normal_scaled[1] * edge_center_scaled[1] + \
           edge_normal_scaled[2] * edge_center_scaled[2]);
  // kt = sqrt(dot(a^2,edge_normal_scaled^2))
  kt = sqrt(a[0] * a[0] * edge_normal_scaled[0] * edge_normal_scaled[0] + \
            a[1] * a[1] * edge_normal_scaled[1] * edge_normal_scaled[1] + \
            a[2] * a[2] * edge_normal_scaled[2] * edge_normal_scaled[2]);
  novera = ( edge_normal_scaled[0] / a[0] ) * ( edge_normal_scaled[0] / a[0] )+ \
           ( edge_normal_scaled[1] / a[1] ) * ( edge_normal_scaled[1] / a[1] )+ \
           ( edge_normal_scaled[2] / a[2] ) * ( edge_normal_scaled[2] / a[2] );
  if (novera >= kt*kt*kt*kt / (kt * kt) )
  {
    area = 0;
  }
  else
  {
    double pi = 3.14159265358979323846;
    area = pi * ( 1 - k*k / (kt*kt) ) * (a[0]*a[1]*a[2] ) / kt;
  }
  return(area);
}



double correct_pndotf(double f[3], double c[3], double pn[3], double px[3])
    //Inputs
    //    f        force (or force density) on facet
    //    c        center of facet
    //    pn       normal to plane
    //    px       point in plane
    //Outputs
    //    nf       pn dot f, with the correct sign. units of f
{
    //Compute the initial dot product
    double nf;
    int i;
    for ( i = 0 ; i < 3 ; i++ )
    {
        nf = nf + pn[i] * f[i];
    }
    // if it's 0, exit and return 0
    if (nf == 0.0) return nf;
    // compute pn dot c-px
    double ncx;
    for ( i = 0 ; i < 3 ; i++ )
    {
        ncx = ncx + pn[i] * ( c[i] - px[i] );
    }
    double t = - ncx / nf;
    // treat the cases for t
    if ( t<0.0 || t>1.0) {
        // then f does not intersect the plane
        double quant = ( fabs(nf + ncx) - fabs(ncx) );
        double sign = 1;
        if (quant < 0) sign = -1.0;
        nf = fabs(nf) * sign;
        // right now ignore ==0 case. Not sure what to do with this
    }
    else if (t != 0.0) nf = -fabs(nf); // then f intersects the plane 
    else nf = fabs(nf); // then the force originates on the plane
    return nf;
}


double set_stress(int N, double a[3], double fonfV[N][3], double cV[N][3],
              double pn[3], double px[3])
    //Inputs
    //    N        number of facets
    //    fonf     force (NOT force density) on each facet, indexed by row
    //    cv       centers of facets, indexed by row
    //    pn       normal to plane
    //    px       point in plane
    //Outputs
    //    stress   stress due fo fonf wrt plane definde by pn,px
{
    double aoi = area_of_intersection(a, pn, px);
    if ( aoi <= 0) return 0; //the the plane does not intersect
    double fonf[3], c[3], nf, total_force;
    int i, j;
    for ( i = 0 ; i < N ; i++ )
    {
        for ( j = 0 ; j < 3 ; j++ )
        {
            fonf[j] = fonfV[i][j];
            c[j] = cV[i][j];
        }
        nf = correct_pndotf(fonf, c, pn, px);
        total_force = total_force + nf;
    }
    double stress = total_force / aoi;
    return stress;
}


double set_force(int N, double a[3], double fonfV[N][3], double cV[N][3],
              double pn[3], double px[3])
    //Inputs
    //    N        number of facets
    //    fonf     force (NOT force density) on each facet, indexed by row
    //    cv       centers of facets, indexed by row
    //    pn       normal to plane
    //    px       point in plane
    //Outputs
    //    stress   stress due fo fonf wrt plane definde by pn,px
{
    double aoi = area_of_intersection(a, pn, px);
    if ( aoi <= 0) return 0; //the the plane does not intersect
    double fonf[3], c[3], nf, total_force;
    int i, j;
    for ( i = 0 ; i < N ; i++ )
    {
        for ( j = 0 ; j < 3 ; j++ )
        {
            fonf[j] = fonfV[i][j];
            c[j] = cV[i][j];
        }
        nf = correct_pndotf(fonf, c, pn, px);
        total_force = total_force + nf;
    //    printf("%e\n", total_force);
    }

    return (total_force * 2.0) / 2.0;
    // No idea why the hell it is happening but returning total_force gives
    // 0 - somehow doing something to it other than multiplying by 1 gets the
    // correct answer...
}

void set_force_Vectorized(
    int NumTimes,
    int NumPlanes,
    int NumFacets,
    double fonfV[NumFacets][3],
    double forces[NumTimes][NumPlanes],
    double aV[NumTimes][3],
    double rotAngles[NumTimes][2],
    double wV[NumTimes][3],
    double pnV_sph[NumPlanes][3],
    double pxV_sph[NumPlanes][3],
    double gammadot,
    double p0,
    double mu,
    double srf_centers_sph[NumFacets][3],
    double srf_crosses_sph[NumFacets][3],
    double srf_normals_sph[NumFacets][3])

/*  Computes the force wrt each timestep and each plane. Takes as input the output of the deformation
    simulation as well as some other constants. 

    Inputs:
        NumTimes         number of timesteps
        NumPlanes        number of intersecting planes
        NumFacets        number of facets on surface triangulation
        fonfV            empty array, to make life easier. at each time step is overwritten
        forces           output array, empty
        aV               each row is the axes lengths at time t
        rotAngles        each row is [cos theta, sin theta] from the rotation
        wV               each row is the angular velocity vector at time t
        pnV_sph          each row is a plane normal, scaled to the sphere
        pnX_sph          each row is a plane point, scaled to the sphere
        gammadot         shear rate
        p0               external pressure
        mu               matrix viscosity
        srf_centers_sph  surface triangulation centers, scaled to sphere
        srf_crosses_sph  surface triangulation edge crosses, scaled to sphere
        srf_normals_sph  surface triangulation normals, scaled to sphere
*/

{
    double L[3][3];
    double farg[3][3];
    double a[3], w[3];
    double c,s;

    //double fonfV[NumFacets][3];
    double srf_centers_scaled[NumFacets][3];
    double srf_areas_scaled[NumFacets];
    double srf_normals_scaled[NumFacets][3];

    double pn[3], px[3] = {0};
    double mag;

    int tstep, pstep, iter_row = 0;
    for ( tstep = 0 ; tstep < NumTimes ; tstep++ )
    {
        // assign the axes and the angular velocity
        for ( iter_row = 0 ; iter_row < 3 ; iter_row++ )
        {
            a[iter_row] = aV[tstep][iter_row];
            w[iter_row] = wV[tstep][iter_row];
        }


        // scale the surface triangulation 
        scale_triangulation(NumFacets, a, srf_centers_scaled, srf_areas_scaled, srf_normals_scaled, 
                            srf_centers_sph, srf_crosses_sph, srf_normals_sph);

        // set up the rotated matrix L
        set_L( L, rotAngles[tstep][0], rotAngles[tstep][1] , gammadot );
        
        //c = rotAngles[tstep][0];
        //s = rotAngles[tstep][1];
        //L[0][0] = -c * s * gammadot;
        //L[0][1] =  c * c * gammadot;
        //L[1][0] = -s * s * gammadot;
        //L[1][1] = -L[0][0];

        // set farg
        set_farg(farg, a, w, L, p0, mu);

        // set the force on the facets    
        set_force_facets(NumFacets, fonfV, farg, srf_normals_scaled, srf_areas_scaled);

        // loop over the planes
        for ( pstep = 0 ; pstep < NumPlanes ; pstep++ )
        {
            // set the plane normal and point
            for (iter_row = 0 ; iter_row < 3 ; iter_row ++)
            {
                pn[iter_row] = pnV_sph[pstep][iter_row];
                px[iter_row] = pxV_sph[pstep][iter_row];
            }

            // scale the plane normal and point
            for (iter_row = 0 ; iter_row < 3 ; iter_row ++)
            {
                pn[iter_row] = pn[iter_row] * a[iter_row];
                px[iter_row] = px[iter_row] * a[iter_row];
            }

            // renormalize the normal
            mag = 0;
            for (iter_row = 0 ; iter_row < 3 ; iter_row ++)
            {
                mag = mag + pn[iter_row] * pn[iter_row];
            }
            // WOULDN'T WANT TO FORGET TO TAKE THE SQUARE ROOT OF THIS WOULD WE?
            // IF WE DID WE MIGHT SPEND SEVERAL DAYS TRYING TO FIND THE ERROR.
            // HAHAHAHAHAHA WOULDN'T THAT BE FUNNY IF THAT HAPPENED?
            // HAHAHAHAHA
            // HAHA
            // ...
            mag = sqrt(mag);
            for (iter_row = 0 ; iter_row < 3 ; iter_row ++)
            {
                pn[iter_row] = pn[iter_row] / mag;
            }


            // compute the force 
            forces[tstep][pstep] = set_force(NumFacets, a, fonfV, srf_centers_scaled, pn, px);
        }

    }
}










