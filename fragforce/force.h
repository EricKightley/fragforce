

// Elliptic Integral functions

double rd_recursion(double x, double y, double z, int n, int N);

double rd_converge(double x, double y, double z);

void set_chi(double lam[3], double a, double b, double c);


// Scaling functions

void scale_edge(double axes[3]
                 double edge_normal_sph[3], 
                 double edge_center_sph[3], 
                 double edge_normal_scaled[3],
                 double edge_center_scaled[3]);

void scale_triangulation(int SIZE, double a[3], 
                                   double srf_centers_scaled[SIZE][3],
                                   double srf_areas_scaled[SIZE],
                                   double srf_normals_scaled[SIZE][3],
                                   double srf_centers_sph[SIZE][3], 
                                   double srf_crosses_sph[SIZE][3],
                                   double srf_normals_sph[SIZE][3]);


// Force functions 

void set_L(double L[3][3], double c, double s, double gammadot);

void set_A(double A[3][3], double w[3], double a[3], double L[3][3], double X[3]);

void set_farg(double farg[3][3], double a[3], double w[3], double L[3][3], double p0, double mu);

void set_force_facets(int SIZE, double force_on_facets[SIZE][3], 
                      double farg[3][3], double srf_normals_scaled[SIZE][3],
                      double srf_areas_scaled[SIZE]);


// Stress functions

double area_of_intersection(double a[3], double edge_normal_scaled[3], double edge_center_scaled[3])

double correct_pndotf(double f[3], double c[3], double pn[3], double px[3]);

double set_stress(int N, double a[3], double fonfV[N][3], double cV[N][3], double pn[3], double px[3]);

double set_force(int N, double a[3], double fonfV[N][3], double cV[N][3], double pn[3], double px[3]);

void set_force_Vectorized(
    int NumTimes,
    int NumPlanes,
    int NumFacets,
    double forces[NumFacets][3],
    double forces[NumTimes][NumPlanes],
    double aV[NumTimes][3],
    double rotAngles[NumTimes][2],
    double wV[NumTimes][3],
    double pnV[NumPlanes][3],
    double pxV[NumPlanes][3],
    double gammadot,
    double p0,
    double mu,
    double srf_centers_sph[NumFacets][3],
    double srf_crosses_sph[NumFacets][3],
    double srf_normals_sph[NumFacets][3]);
