

// Elliptic Integral functions

double rd_recursion(double x, double y, double z, int n, int N);

double rd_converge(double x, double y, double z);

void set_chi(double chi[3],
             double a,
             double b,
             double c);


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

void set_L(double L[3][3],
           double c,
           double s, 
           double gammadot);

void set_A(double A[3][3], 
           double a[3],
           double w[3],
           double L[3][3],
           double chi[3]);

void set_farg(double farg[3][3],
              double a[3],
              double w[3],
              double L[3][3],
              double A[3][3],
              double chi[3],
              double p0,
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
        double wV[NTimes][3],
        double pnV[NPlanes][3], 
        double pxV[NPlanes][3],
        double srf_centers_sph[NFacets][3],
        double srf_crosses_sph[NFacets][3],
        double srf_normals_sph[NFacets][3],
        double gammadot,
        double p0,
        double mu,
        int scale_planes_bool);




