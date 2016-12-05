// DE Transform Quadrature

void intdeiini(int lenaw, 
               double tiny,
               double eps,
               double *aw);

void intdei(double (*func)(double, double *),
            double *funcpars,
            double a,
            double *aw,
            double *integral, 
            double *err);


double hypergeo3(double t,
                 double pV[6]);

double integrate_hypergeo3(double z1,
                           double z2,
                           double z3,
                           double b1,
                           double b2,
                           double b3);


// Scaling functions

void scale_plane(double a_initial[3],
                 double a_current[3],
                 double pn_initial[3], 
                 double px_initial[3], 
                 double pn_scaled[3],
                 double px_scaled[3]);

void scale_triangulation(int NFacets, 
                         double a[3], 
                         double srf_centers_scaled[NFacets][3],
                         double srf_areas_scaled[NFacets],
                         double srf_normals_scaled[NFacets][3],
                         double srf_centers_sph[NFacets][3], 
                         double srf_crosses_sph[NFacets][3],
                         double srf_normals_sph[NFacets][3]);


// Force functions 

void set_farg(double farg[3][3],
              double a[3],
              double R[2],
              double w,
              double gammadot,
              double mu);

void set_force_density(int NFacets,
                      double fdonf[NFacets][3],
                      double farg[3][3], 
                      double srf_normals_scaled[NFacets][3]);

void set_force_facets(int NFacets, 
                      double fonfV[NFacets][3], 
                      double fdonf[NFacets][3], 
                      double srf_areas_scaled[NFacets]);

double area_of_intersection(double a[3],
                            double pn_scaled[3],
                            double px_scaled[3]);

double correct_pndotf(double fonf[3],
                      double srf_center_scaled[3],
                      double pn_scaled[3],
                      double px_scaled[3]);

double sum_forces(int NFacets,
                  double a[3],
                  double fonfV[NFacets][3],
                  double srf_centers_scaled[NFacets][3],
                  double pn_scaled[3], 
                  double px_scaled[3]);

void frag_force(
        int NTimes, 
        int NPlanes, 
        int NFacets,
        double fragforceV[NTimes][NPlanes],
        double aV[NTimes][3], 
        double RV[NTimes][2],
        double wV[NTimes],
        double pnV[NPlanes][3], 
        double pxV[NPlanes][3],
        double srf_centers_sph[NFacets][3],
        double srf_crosses_sph[NFacets][3],
        double srf_normals_sph[NFacets][3],
        double gammadot,
        double mu,
        int scale_planes_bool);




