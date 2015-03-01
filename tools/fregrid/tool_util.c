

#include <stdlib.h> 
#include <string.h>
#include <math.h> 
#include "constant.h" 
#include "mosaic_util.h" 
#include "tool_util.h"
#include "interp.h"
#include "mpp.h"
#define  D2R (M_PI/180.)
#define  R2D (180./M_PI)

const double SMALL = 1.0e-4;
double distant(double a, double b, double met1, double met2);
double bp_lam(double x, double y, double bpeq, double rp);
double bp_phi(double x, double y, double bpsp, double bpnp);
double lon_in_range(double lon, double lon_strt);
void vtx_insert(double *x, double *y, int *n, int n_ins);
void vtx_delete(double *x, double *y, int *n, int n_del);
int lon_fix(double *x, double *y, int n_in, double tlon);

/***************************************************************************
  void get_file_path(const char *file, char *dir)
  get the directory where file is located. The dir will be the complate path
  before the last "/". If no "/" exist in file, the path will be current ".".
***************************************************************************/
void get_file_path(const char *file, char *dir)
{
  int len;
  char *strptr = NULL;

  /* get the diretory */
 
  strptr = strrchr(file, '/');
  if(strptr) {
    len = strptr - file;
    strncpy(dir, file, len);
  }
  else {
    len = 1;
    strcpy(dir, ".");
  }
  dir[len] = 0;

}; /* get_file_path */

int get_int_entry(char *line, int* value)
{
  char* pch;
  int num;
  
  pch = strtok(line, ", ");
  num = 0;
  while( pch != NULL) {
    value[num++] = atoi(pch);
    pch = strtok(NULL, ", ");
  }
  return num;
    
};

int get_double_entry(char *line, double *value)
{
  char* pch;
  int num;
  
  pch = strtok(line, ", ");
  num = 0;
  while( pch != NULL) {
    value[num++] = atof(pch);
    pch = strtok(NULL, ", ");
  }
  return num;
};

/*********************************************************************
  double spherical_dist(double x1, double y1, double x2, double y2)
  return distance between spherical grid on the earth
*********************************************************************/

double spherical_dist(double x1, double y1, double x2, double y2)
{
  double dist = 0.0;
  double h1, h2;
  
  if(x1 == x2) {
    h1 = RADIUS;
    h2 = RADIUS;
    dist = distant(y1,y2,h1,h2);
  }
  else if(y1 == y2) {
    h1 = RADIUS * cos(y1*D2R);
    h2 = RADIUS * cos(y2*D2R);
    dist = distant(x1,x2,h1,h2);
  }
  else 
    mpp_error("tool_till: This is not rectangular grid");

  return dist;
}; /* spherical_dist */
  

/*********************************************************************
  void double bipolar_dist(double x1, double y1, double x2, double y2)
  return distance of bipolar grids
*********************************************************************/
double bipolar_dist(double x1, double y1, double x2, double y2,
		    double bpeq, double bpsp, double bpnp, double rp )
{
  double dist, x[2],y[2], bp_lon[2], bp_lat[2], metric[2];
  double h1[2], h2[2], chic;
  int n;
  
  x[0] = x1;  x[1] = x2;
  y[0] = y1;  y[1] = y2;
  
  /*--- get the bipolar grid and metric term ----------------------------*/
  for(n=0; n<2; n++){
    bp_lon[n] = bp_lam(x[n],y[n],bpeq, rp);     /* longitude (degrees) in bipolar grid system */
    bp_lat[n] = bp_phi(x[n],y[n],bpsp, bpnp);  /* latitude (degrees) in bipolar grid system */
    h1[n]     = RADIUS*cos(bp_lat[n]*D2R);
    h2[n]     = RADIUS;
    metric[n] = 1.0;
    if (fabs(y[n]-90.0) < SMALL || fabs(bp_lon[n]*D2R) >= SMALL
	|| fabs(bp_lat[n]*D2R) >= SMALL) {
      chic = acos(cos(bp_lon[n]*D2R)*cos(bp_lat[n]*D2R));            /* eqn. 6 */
      metric[n] = rp*(1/pow(cos(chic/2),2))/(1+(pow(rp,2))*(pow(tan(chic/2),2)));/* eq 3 */
    }
  }

  /*--- then calculate the distance -------------------------------------*/
  if(x1 == x2) 
    dist = distant(bp_lon[0],bp_lon[1],metric[0]*h1[0],metric[1]*h1[1]);
  else if(y1 == y2) 
    dist = distant(bp_lat[0],bp_lat[1],metric[0]*h2[0],metric[1]*h2[1]);
  else
    mpp_error("tool_util: This tripolar grid not transformed from rectangular grid");    

  return dist;
  
}; /* bipolar_dist */

/*********************************************************************
  double distant(double a, double b, double met1, double met2)
  return distant on the earth
*********************************************************************/
double distant(double a, double b, double met1, double met2)
{
   return fabs(a-b)*D2R*(met1+met2)/2. ;
}; /* distant */

/*********************************************************************
   double spherical_area(double x1, double y1, double x2, double y2,
                   double x3, double y3, double x4, double y4 )            
   rectangular grid box area
 ********************************************************************/
double spherical_area(double x1, double y1, double x2, double y2,
		      double x3, double y3, double x4, double y4 )
{
  double area, dx, lat1, lat2, x[4],y[4];
  int i, ip;
  
  x[0] = x1; y[0] = y1;
  x[1] = x2; y[1] = y2;
  x[2] = x3; y[2] = y3;
  x[3] = x4; y[3] = y4;

  area = 0.0;

  for(i=0; i<4; i++) {
    ip = i+1;
    if(ip ==4) ip = 0;
    dx = (x[ip] - x[i])*D2R;
    lat1 = y[ip]*D2R;
    lat2 = y[i]*D2R;
    if(dx==0.0) continue;
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx < -M_PI) dx = dx + 2.0*M_PI;

    if (lat1 == lat2) /* cheap area calculation along latitude  */
      area = area - dx*sin(lat1);
    else 
      area = area - dx*(sin(lat1)+sin(lat2))/2;   /*  TRAPEZOID_RULE */
  }

  area = area * RADIUS * RADIUS;

  return area;
}; /* spherical_area */

/*********************************************************************
   double bipolar_area(double x1, double y1, double x2, double y2,
                       double x3, double y3, double x4, double y4 )            
   bipolar grid  area
 ********************************************************************/
double bipolar_area(double x1, double y1, double x2, double y2,
			  double x3, double y3, double x4, double y4 )
{
  double area, dx, lat1, lat2, x[8],y[8];
  int i, ip, n;
  
  x[0] = x1; y[0] = y1;
  x[1] = x2; y[1] = y2;
  x[2] = x3; y[2] = y3;
  x[3] = x4; y[3] = y4;


  /*--- first fix the longitude at the pole -----------------------------*/
  n = lon_fix(x, y, 4, 180.);

  /*--- calculate the area ----------------------------------------------  */
  area = 0.0;  
  for(i=0; i<n; i++){
    ip = i+1;
    if(ip == n) ip = 0;
    dx   = (x[ip] - x[i])*D2R;
    lat1 = y[ip]*D2R;
    lat2 = y[i]*D2R;
    if(dx==0.0) continue;
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx < -M_PI) dx = dx + 2.0*M_PI;

    if (lat1 == lat2)  /* cheap area calculation along latitude */
      area = area - dx*sin(lat1);
    else
      area = area - dx*(sin(lat1)+sin(lat2))/2;   /*  TRAPEZOID_RULE */
  }
  
  area = area * RADIUS * RADIUS;

  return area;
}; /* bipolar_area */

/*********************************************************************
  double lat_dist(double x1, double x2)
  distance (in degrees) between points on lat. circle
 ********************************************************************/
  double lat_dist(double x1, double x2)
{
  return min(fmod(x1-x2+720,360.),fmod(x2-x1+720,360.));
};


/*********************************************************************
  double bp_lam(double x, double y, double bpeq)
  find bipolar grid longitude given geo. coordinates
 ********************************************************************/
  double bp_lam(double x, double y, double bpeq, double rp)
{
  double bp_lam;

  /*  bp_lam = ((90-y)/(90-lat_join))*90 */
  /* invert eqn. 5 with phic=0 to place point at specified geo. lat */
  bp_lam = 2.*atan(tan((0.5*M_PI-y*D2R)/2)/rp)*R2D;
  if (lat_dist(x,bpeq)<90.) bp_lam = -bp_lam;
  return bp_lam;
}; /* bp_lam */

/*********************************************************************
   double bp_phi(double x, double y, double bpsp, double bpnp)
   find bipolar grid latitude given geo. coordinates
 ********************************************************************/
   double bp_phi(double x, double y, double bpsp, double bpnp)
{
  double bp_phi;

  if (lat_dist(x,bpsp)<90.)
    return (-90+lat_dist(x,bpsp));
  else
    return ( 90-lat_dist(x,bpnp));
}; /* bp_phi */


/*********************************************************************
  void tp_trans(double& lon, double& lat, double lon_ref)
  calculate tripolar grid
 ********************************************************************/
void tp_trans(double *lon, double *lat, double lon_ref, double lon_start, 
		    double lam0, double bpeq, double bpsp, double bpnp, double rp )
{
  double lamc, phic, lams, chic, phis;
  
  lamc = bp_lam(*lon, *lat, bpeq, rp )*D2R;
  phic = bp_phi(*lon, *lat, bpsp, bpnp)*D2R;

  if (fabs(*lat-90.) < SMALL) {
       if (phic > 0) 
	 *lon=lon_in_range(lon_start,lon_ref);
       else
	 *lon=lon_start+180.;
       chic = acos(cos(lamc)*cos(phic));                     /* eqn. 6 */
       phis = M_PI*0.5-2*atan(rp*tan(chic/2));                   /* eqn. 5 */
       *lat = phis*R2D;
       return;
  }

  if (fabs(lamc) < SMALL && fabs(phic) < SMALL) {
    *lat=90.;
    *lon=lon_ref;
  }
  else {
    lams = fmod(lam0+M_PI+M_PI/2-atan2(sin(lamc),tan(phic)),2*M_PI);  /* eqn. 5 */
    chic = acos(cos(lamc)*cos(phic));                          /* eqn. 6 */
    phis = M_PI*0.5-2*atan(rp*tan(chic/2));                        /* eqn. 5 */
    *lon = lams*R2D;
    *lon = lon_in_range(*lon,lon_ref); 
    *lat = phis*R2D;
  }
}; /* tp_trans */

/*********************************************************************
  double Lon_in_range(double lon, double lon_strt)
  Returns lon_strt <= longitude <= lon_strt+360
 ********************************************************************/
double lon_in_range(double lon, double lon_strt)
{
  double lon_in_range, lon_end;

  lon_in_range = lon;
  lon_end = lon_strt+360.;

  if (fabs(lon_in_range - lon_strt) < SMALL) 
    lon_in_range = lon_strt;
  else if (fabs(lon_in_range - lon_end) < SMALL)
    lon_in_range = lon_strt;
  else {
    while(1) {
      if (lon_in_range < lon_strt)          
	lon_in_range = lon_in_range +  360.;
      else if (lon_in_range  >  lon_end)
	lon_in_range  = lon_in_range - 360.;
      else
	break;
    }
  }
  return lon_in_range;
}; /* lon_in_range */


/*********************************************************************
   int lon_fix(double *x, double *y, int n_in, double tlon) 
   fix longitude at pole.
 ********************************************************************/
int lon_fix(double *x, double *y, int n_in, double tlon)
{
  int i, ip, im, n_out;
  double x_sum, dx;
  
  n_out = n_in;
  i     = 0;
  while( i < n_out) {
    if(fabs(y[i]) >= 90.-SMALL) {
      im = i - 1;
      if(im < 0) im = im + n_out;
      ip = i + 1;
      if(ip >= n_out) ip = ip - n_out;
      /*--- all pole points must be paired ---------------------------- */
      if(y[im] == y[i] && y[ip] == y[i] ) {
	vtx_delete(x,y, &n_out, i);
	i = i - 1;
      }
      else if(y[im] != y[i] && y[ip] != y[i] ) {
        vtx_insert(x,y,&n_out,i);
	i = i + 1;
      }
    }
    i = i + 1;
  }

  /*--- first of pole pair has longitude of previous vertex -------------
    --- second of pole pair has longitude of subsequent vertex ---------- */
  for(i=0;i<n_out;i++){
    if(fabs(y[i]) >= 90.-SMALL) {
      im= i - 1;
      if(im < 0) im = im + n_out;
      ip = i + 1;
      if(ip >= n_out) ip = ip - n_out;

      if(y[im] != y[i]) x[i] = x[im];
      if(y[ip] != y[i]) x[i] = x[ip];
    }
  }

  if(n_out == 0) return 0;

  x_sum = x[1];
  for(i=1;i< n_out;i++){
    dx = x[i] - x[i-1];
    if(dx < -180) 
      dx = dx + 360;
    else if (dx >  180)
      dx = dx - 360;

    x[i] = x[i-1] + dx;
    x_sum = x_sum + x[i];
  }

  dx = x_sum/(n_out) - tlon;
  if (dx < -180.) 
    for(i=0;i<n_out;i++) x[i] = x[i] + 360.;
  else if (dx > 180.)
    for(i=0;i<n_out;i++) x[i] = x[i] - 360.;

  return n_out;
  
}; /* lon_fix */


/*********************************************************************
   void vtx_delete(double *x, double *y, int *n, int n_del)
   delete vertex
 ********************************************************************/
void vtx_delete(double *x, double *y, int *n, int n_del)
{
  int i;

  for(i=n_del; i<=*n-2; i++)
    {
      x[i] = x[i+1];
      y[i] = y[i+1];
    }
  (*n)--;
}; /* vtx_delete */

/*********************************************************************
   void Vtx_insert(double *x, double *y, int *n, int n_del)
   insert vertex
 ********************************************************************/
void vtx_insert(double *x, double *y, int *n, int n_ins)
{
  int i;

  for(i=*n-1; i>=n_ins; i--){
    x[i+1] = x[i];
    y[i+1] = y[i];
  }
  (*n)++;

}; /* vtx_insert */


/*----------------------------------------------------------------------
    void vect_cross(e, p1, p2)
    Perform cross products of 3D vectors: e = P1 X P2
    -------------------------------------------------------------------*/
    
/********************************************************************************
  void compute_grid_bound(int nb, const couble *bnds, const int *npts, int *grid_size, const char *center_cell)
  compute the 1-D grid location.
********************************************************************************/
double* compute_grid_bound(int nb, const double *bnds, const int *npts, int *grid_size, const char *center)
{
  int    refine, i, n, np;
  double *grid=NULL, *tmp=NULL;
  double *grid1=NULL, *grid2=NULL;

  if(!strcmp(center, "none") )
    refine = 1;
  else if(!strcmp(center, "t_cell") || !strcmp(center, "c_cell") )
    refine = 2;
  else
    mpp_error("tool_util: center should be 'none', 'c_cell' or 't_cell' ");
	  
  grid1 = (double *)malloc(nb*sizeof(double));
  grid1[0] = 1;
  n = 0;
  for(i=1; i<nb; i++) {
    if(npts[i-1]%refine) mpp_error("tool_util: when center_cell is not 'none', npts should be divided by 2");
    n += npts[i-1]/refine;
    grid1[i] = n+1;
  }
  np = n + 1;
  *grid_size = n*refine;
  tmp   = (double *)malloc(np*sizeof(double));
  grid  = (double *)malloc((*grid_size+1)*sizeof(double));
  grid2 = (double *)malloc(np*sizeof(double));
  for(i=0;i<np;i++) grid2[i] = i + 1.0;

  cubic_spline( nb, np, grid1, grid2, bnds, tmp, 1e30, 1e30);
  if(!strcmp(center, "none")) {
    for(i=0; i<np; i++) grid[i] = tmp[i];
  }
  else if(!strcmp(center, "t_cell")) {
    for(i=0; i<np; i++) grid[2*i] = tmp[i];
    for(i=0; i<n;  i++) grid[2*i+1] = 0.5*(tmp[i]+tmp[i+1]);
  }
  else if( !strcmp(center, "c_cell")) {
    for(i=0; i<np; i++) grid[2*i] = tmp[i];
    grid[1] = 0.5*(tmp[0]+tmp[1]);
    for(i=1; i<n;  i++) grid[2*i+1] = 2*grid[2*i] - grid[2*i-1];
  }
    
  free(grid1);
  free(grid2);
  free(tmp);  

  return grid;
  
};/* compute_grid_bound */
