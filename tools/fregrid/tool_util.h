/***********************************************************************
                      tool_util.h
    This header file provide some utilities routine that will be used in many tools.
    
    contact: Zhi.Liang@noaa.gov
***********************************************************************/
#ifndef TOOL_UTIL_H_
#define TOOL_UTIL_H_
void get_file_path(const char *file, char *dir);
int get_int_entry(char *line, int *value); 
int get_double_entry(char *line, double *value);
double spherical_dist(double x1, double y1, double x2, double y2);
double spherical_area(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4 ); 
double bipolar_dist(double x1, double y1, double x2, double y2, double bpeq, double bpsp, double bpnp, double rp );
double bipolar_area(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4 );
void tp_trans(double *lon, double *lat, double lon_ref, double lon_start, 
              double lam0, double bpeq, double bpsp, double bpnp, double rp );
double* compute_grid_bound(int nb, const double *bnds, const int *npts, int *grid_size, const char *center);
#endif
