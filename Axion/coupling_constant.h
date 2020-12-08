#include <iostream>
#include <fstream>

using namespace std;

const double B0 = 1; 

double dfunction_alpha_TMTE_cos_p(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double z, double zs, double rho, double rhos, double phi, double phis)
{
  int n1 = (int) N1;
  int p1 = (int) P1;
  int q1 = (int) Q1;
  int n2 = (int) N2;
  int p2 = (int) P2;
  int q2 = (int) Q2; 
  int n3 = (int) N3;
  int p3 = (int) P3;
  int q3 = (int) Q3; 

  double R = sqrt(rs(z,zs,rho,rhos,phi,phis));
  double omega_p = omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);

  if (R == 0)
     return 0;
  else
  {
  
     if (omega_p > m)
        return (rho*B0*Detection_modeTM_Ez(n3,p3,q3,dist,R2,L2,z,rho,phi))*(rhos*F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*cos(R*sqrt(omega_p*omega_p - m*m))/R);
     else
        return (rho*B0*Detection_modeTM_Ez(n3,p3,q3,dist,R2,L2,z,rho,phi))*(rhos*F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*exp(-R*sqrt(m*m - omega_p*omega_p))/R);
 
  } 
}

int ifunction_alpha_TMTE_cos_p(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double *p = ((double *) fdata); 
    fval[0] = dfunction_alpha_TMTE_cos_p(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[16],x[0],x[1],x[2],x[3],x[4],x[5]);
    return 0; 
}


double alpha_TMTE_cos_p(int rank, int size, double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double eps)
{

   double dzs = L1/size;
   double zs_min = dzs*rank;
   double zs_max = dzs*(rank + 1);

   ofstream out("Messages.txt",std::ios::app); 
 
   double xmin[6] = {dist,zs_min,0,0,0,0}, xmax[6] = {dist + L2,zs_max,R2,R1,2*pi,2*pi}, p[17] = {N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m}, val, err;   
   hcubature(1, ifunction_alpha_TMTE_cos_p, &p, 6, xmin, xmax, 0, 0, eps, ERROR_INDIVIDUAL, &val, &err);
   out  << "alpha_TMTE_cos_p, m = " << m << ", rank = " << rank << ", val = " << val << ", err =  " << err << ", fabs(err/val) = " << fabs(err/val) << endl; 

   out.close();
 
   return val;


}

double dfunction_alpha_TMTE_sin_p(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double z, double zs, double rho, double rhos, double phi, double phis)
{
  int n1 = (int) N1;
  int p1 = (int) P1;
  int q1 = (int) Q1;
  int n2 = (int) N2;
  int p2 = (int) P2;
  int q2 = (int) Q2; 
  int n3 = (int) N3;
  int p3 = (int) P3;
  int q3 = (int) Q3; 

  double R = sqrt(rs(z,zs,rho,rhos,phi,phis));
  double omega_p = omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);

  if (R == 0)
     return 0;
  else
  {
  
     if (omega_p > m)
        return (rho*B0*Detection_modeTM_Ez(n3,p3,q3,dist,R2,L2,z,rho,phi))*(rhos*F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*sin(R*sqrt(omega_p*omega_p - m*m))/R);
     else
        return 0.;
 
  } 
}

int ifunction_alpha_TMTE_sin_p(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double *p = ((double *) fdata); 
    fval[0] = dfunction_alpha_TMTE_sin_p(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[16],x[0],x[1],x[2],x[3],x[4],x[5]);
    return 0; 
}


double alpha_TMTE_sin_p(int rank, int size, double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double eps)
{
   double dzs = L1/size;
   double zs_min = dzs*rank;
   double zs_max = dzs*(rank + 1);

   ofstream out("Messages.txt",std::ios::app); 
 
   double xmin[6] = {dist,zs_min,0,0,0,0}, xmax[6] = {dist + L2,zs_max,R2,R1,2*pi,2*pi}, p[17] = {N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m}, val, err;     
   hcubature(1, ifunction_alpha_TMTE_sin_p, &p, 6, xmin, xmax, 0, 0, eps, ERROR_INDIVIDUAL, &val, &err);
   out << "alpha_TMTE_sin_p, m = " << m << ", rank = " << rank <<", val = " << val << ", err =  " << err << ", fabs(err/val) = " << fabs(err/val) << endl; 

   out.close(); 

   return val;
  
}


double dfunction_beta_TMTE_cos_p(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double z, double zs, double rho, double rhos, double phi, double phis)
{
  int n1 = (int) N1;
  int p1 = (int) P1;
  int q1 = (int) Q1;
  int n2 = (int) N2;
  int p2 = (int) P2;
  int q2 = (int) Q2; 
  int n3 = (int) N3;
  int p3 = (int) P3;
  int q3 = (int) Q3; 

  double R = sqrt(rs(z,zs,rho,rhos,phi,phis));
  double omega_p = omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);

  if (R == 0)
     return 0;
  else
  {
  
     if (omega_p > m)
        return (rho*B0*Detection_modeTM_Ez(n3,p3,q3,dist,R2,L2,z,rho,phi))*(rhos*partial2_z_F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*cos(R*sqrt(omega_p*omega_p - m*m))/R);
     else
        return (rho*B0*Detection_modeTM_Ez(n3,p3,q3,dist,R2,L2,z,rho,phi))*(rhos*partial2_z_F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*exp(-R*sqrt(m*m - omega_p*omega_p))/R);
 
  } 
}

int ifunction_beta_TMTE_cos_p(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double *p = ((double *) fdata); 
    fval[0] = dfunction_beta_TMTE_cos_p(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[16],x[0],x[1],x[2],x[3],x[4],x[5]);
    return 0; 
}


double beta_TMTE_cos_p(int rank, int size, double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double eps)
{

   double dzs = L1/size;
   double zs_min = dzs*rank;
   double zs_max = dzs*(rank + 1);

   ofstream out("Messages.txt",std::ios::app); 
 
   double xmin[6] = {dist,zs_min,0,0,0,0}, xmax[6] = {dist + L2,zs_max,R2,R1,2*pi,2*pi}, p[17] = {N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m}, val, err;     
   hcubature(1, ifunction_beta_TMTE_cos_p, &p, 6, xmin, xmax, 0, 0, eps, ERROR_INDIVIDUAL, &val, &err);
   out  << "beta_TMTE_cos_p, m = " << m << ", rank = " << rank << ", val = " << val << ", err =  " << err << ", fabs(err/val) = " << fabs(err/val) << endl; 

   out.close();
 
   return val;

}

double dfunction_beta_TMTE_sin_p(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double z, double zs, double rho, double rhos, double phi, double phis)
{
  int n1 = (int) N1;
  int p1 = (int) P1;
  int q1 = (int) Q1;
  int n2 = (int) N2;
  int p2 = (int) P2;
  int q2 = (int) Q2; 
  int n3 = (int) N3;
  int p3 = (int) P3;
  int q3 = (int) Q3; 

  double R = sqrt(rs(z,zs,rho,rhos,phi,phis));
  double omega_p = omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);

  if (R == 0)
     return 0;
  else
  {
  
     if (omega_p > m)
        return (rho*B0*Detection_modeTM_Ez(n3,p3,q3,dist,R2,L2,z,rho,phi))*(rhos*partial2_z_F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*sin(R*sqrt(omega_p*omega_p - m*m))/R);
     else
        return 0.;
 
  } 
}

int ifunction_beta_TMTE_sin_p(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double *p = ((double *) fdata); 
    fval[0] = dfunction_beta_TMTE_sin_p(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[16],x[0],x[1],x[2],x[3],x[4],x[5]);
    return 0; 
}


double beta_TMTE_sin_p(int rank, int size, double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double eps)
{
   double dzs = L1/size;
   double zs_min = dzs*rank;
   double zs_max = dzs*(rank + 1);

   ofstream out("Messages.txt",std::ios::app); 

   double xmin[6] = {dist,zs_min,0,0,0,0}, xmax[6] = {dist + L2,zs_max,R2,R1,2*pi,2*pi}, p[17] = {N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m}, val, err;     
   hcubature(1, ifunction_beta_TMTE_sin_p, &p, 6, xmin, xmax, 0, 0, eps, ERROR_INDIVIDUAL, &val, &err);
   out << "beta_TMTE_sin_p, m = " << m << ", rank = " << rank <<", val = " << val << ", err =  " << err << ", fabs(err/val) = " << fabs(err/val) << endl; 
 
   out.close();

   return val;
  
}


double dfunction_gamma1_TMTE_cos_p(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double z, double rho, double rhos, double phi, double phis)
{
  int n1 = (int) N1;
  int p1 = (int) P1;
  int q1 = (int) Q1;
  int n2 = (int) N2;
  int p2 = (int) P2;
  int q2 = (int) Q2; 
  int n3 = (int) N3;
  int p3 = (int) P3;
  int q3 = (int) Q3; 

  double zs = L1;

  double R = sqrt(rs(z,zs,rho,rhos,phi,phis));
  double omega_p = omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);

  if (R == 0)
     return 0;
  else
  {
  
     if (omega_p > m)
        return (rho*B0*Detection_modeTM_Ez(n3,p3,q3,dist,R2,L2,z,rho,phi))*(rhos*partial_z_F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*cos(R*sqrt(omega_p*omega_p - m*m))/R);
     else
        return (rho*B0*Detection_modeTM_Ez(n3,p3,q3,dist,R2,L2,z,rho,phi))*(rhos*partial_z_F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*exp(-R*sqrt(m*m - omega_p*omega_p))/R);
 
  } 
}

int ifunction_gamma1_TMTE_cos_p(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double *p = ((double *) fdata); 
    fval[0] = dfunction_gamma1_TMTE_cos_p(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[16],x[0],x[1],x[2],x[3],x[4]);
    return 0; 
}


double gamma1_TMTE_cos_p(int rank, int size, double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double eps)
{

   double drhos = R1/size;
   double rhos_min = drhos*rank;
   double rhos_max = drhos*(rank + 1);

   ofstream out("Messages.txt",std::ios::app); 
 
   double xmin[5] = {dist,0,rhos_min,0,0}, xmax[5] = {dist + L2,R2,rhos_max,2*pi,2*pi}, p[17] = {N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m}, val, err;     
   hcubature(1, ifunction_gamma1_TMTE_cos_p, &p, 5, xmin, xmax, 0, 0, eps, ERROR_INDIVIDUAL, &val, &err);
   out  << "gamma1_TMTE_cos_p, m = " << m << ", rank = " << rank << ", val = " << val << ", err =  " << err << ", fabs(err/val) = " << fabs(err/val) << endl; 

   out.close();
 
   return val;

}



double dfunction_gamma1_TMTE_sin_p(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double z, double rho, double rhos, double phi, double phis)
{
  int n1 = (int) N1;
  int p1 = (int) P1;
  int q1 = (int) Q1;
  int n2 = (int) N2;
  int p2 = (int) P2;
  int q2 = (int) Q2; 
  int n3 = (int) N3;
  int p3 = (int) P3;
  int q3 = (int) Q3; 

  double zs = L1;

  double R = sqrt(rs(z,zs,rho,rhos,phi,phis));
  double omega_p = omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);

  if (R == 0)
     return 0;
  else
  {
  
     if (omega_p > m)
        return (rho*B0*Detection_modeTM_Ez(n3,p3,q3,dist,R2,L2,z,rho,phi))*(rhos*partial_z_F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*sin(R*sqrt(omega_p*omega_p - m*m))/R);
     else
        return 0.;
 
  } 
}

int ifunction_gamma1_TMTE_sin_p(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double *p = ((double *) fdata); 
    fval[0] = dfunction_gamma1_TMTE_sin_p(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[16],x[0],x[1],x[2],x[3],x[4]);
    return 0; 
}


double gamma1_TMTE_sin_p(int rank, int size, double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double eps)
{

   double drhos = R1/size;
   double rhos_min = drhos*rank;
   double rhos_max = drhos*(rank + 1);

   ofstream out("Messages.txt",std::ios::app); 
 
   double xmin[5] = {dist,0,rhos_min,0,0}, xmax[5] = {dist + L2,R2,rhos_max,2*pi,2*pi}, p[17] = {N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m}, val, err;     
   hcubature(1, ifunction_gamma1_TMTE_sin_p, &p, 5, xmin, xmax, 0, 0, eps, ERROR_INDIVIDUAL, &val, &err);
   out  << "gamma1_TMTE_sin_p, m = " << m << ", rank = " << rank << ", val = " << val << ", err =  " << err << ", fabs(err/val) = " << fabs(err/val) << endl; 

   out.close();
 
   return val;

}


double dfunction_gamma2_TMTE_cos_p(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double z, double rho, double rhos, double phi, double phis)
{
  int n1 = (int) N1;
  int p1 = (int) P1;
  int q1 = (int) Q1;
  int n2 = (int) N2;
  int p2 = (int) P2;
  int q2 = (int) Q2; 
  int n3 = (int) N3;
  int p3 = (int) P3;
  int q3 = (int) Q3; 

  double zs = 0.;

  double R = sqrt(rs(z,zs,rho,rhos,phi,phis));
  double omega_p = omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);

  if (R == 0)
     return 0;
  else
  {
  
     if (omega_p > m)
        return (rho*B0*Detection_modeTM_Ez(n3,p3,q3,dist,R2,L2,z,rho,phi))*(rhos*partial_z_F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*cos(R*sqrt(omega_p*omega_p - m*m))/R);
     else
        return (rho*B0*Detection_modeTM_Ez(n3,p3,q3,dist,R2,L2,z,rho,phi))*(rhos*partial_z_F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*exp(-R*sqrt(m*m - omega_p*omega_p))/R);
 
  } 
}

int ifunction_gamma2_TMTE_cos_p(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double *p = ((double *) fdata); 
    fval[0] = dfunction_gamma2_TMTE_cos_p(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[16],x[0],x[1],x[2],x[3],x[4]);
    return 0; 
}


double gamma2_TMTE_cos_p(int rank, int size, double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double eps)
{

   double drhos = R1/size;
   double rhos_min = drhos*rank;
   double rhos_max = drhos*(rank + 1);

   ofstream out("Messages.txt",std::ios::app); 
 
   double xmin[5] = {dist,0,rhos_min,0,0}, xmax[5] = {dist + L2,R2,rhos_max,2*pi,2*pi}, p[17] = {N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m}, val, err;     
   hcubature(1, ifunction_gamma2_TMTE_cos_p, &p, 5, xmin, xmax, 0, 0, eps, ERROR_INDIVIDUAL, &val, &err);
   out  << "gamma2_TMTE_cos_p, m = " << m << ", rank = " << rank << ", val = " << val << ", err =  " << err << ", fabs(err/val) = " << fabs(err/val) << endl; 

   out.close();
 
   return val;

}



double dfunction_gamma2_TMTE_sin_p(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double z, double rho, double rhos, double phi, double phis)
{
  int n1 = (int) N1;
  int p1 = (int) P1;
  int q1 = (int) Q1;
  int n2 = (int) N2;
  int p2 = (int) P2;
  int q2 = (int) Q2; 
  int n3 = (int) N3;
  int p3 = (int) P3;
  int q3 = (int) Q3; 

  double zs = 0.;

  double R = sqrt(rs(z,zs,rho,rhos,phi,phis));
  double omega_p = omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);

  if (R == 0)
     return 0;
  else
  {
  
     if (omega_p > m)
        return (rho*B0*Detection_modeTM_Ez(n3,p3,q3,dist,R2,L2,z,rho,phi))*(rhos*partial_z_F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*sin(R*sqrt(omega_p*omega_p - m*m))/R);
     else
        return 0.;
 
  } 
}

int ifunction_gamma2_TMTE_sin_p(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double *p = ((double *) fdata); 
    fval[0] = dfunction_gamma2_TMTE_sin_p(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[16],x[0],x[1],x[2],x[3],x[4]);
    return 0; 
}


double gamma2_TMTE_sin_p(int rank, int size, double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double dist, double N3, double P3, double Q3, double R2, double L2, double m, double eps)
{

   double drhos = R1/size;
   double rhos_min = drhos*rank;
   double rhos_max = drhos*(rank + 1);

   ofstream out("Messages.txt",std::ios::app); 
 
   double xmin[5] = {dist,0,rhos_min,0,0}, xmax[5] = {dist + L2,R2,rhos_max,2*pi,2*pi}, p[17] = {N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m}, val, err;     
   hcubature(1, ifunction_gamma2_TMTE_sin_p, &p, 5, xmin, xmax, 0, 0, eps, ERROR_INDIVIDUAL, &val, &err);
   out  << "gamma2_TMTE_sin_p, m = " << m << ", rank = " << rank << ", val = " << val << ", err =  " << err << ", fabs(err/val) = " << fabs(err/val) << endl; 

   out.close();
 
   return val;

}



void Coupling_Constant_TMTE_p(int argn, char** argv, int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double dist, int n3, int p3, int q3, double R2, double L2, double m_min, double m_max, double dm, double eps)
{

   double N1 = (double) n1;
   double P1 = (double) p1;
   double Q1 = (double) q1;
   double N2 = (double) n2;
   double P2 = (double) p2;
   double Q2 = (double) q2;
   double N3 = (double) n3;
   double P3 = (double) p3;
   double Q3 = (double) q3;

   char FN1 = n1 + '0';
   char FP1 = p1 + '0';                           
   char FQ1 = q1 + '0';

   char FN2 = n2 + '0';
   char FP2 = p2 + '0';                           
   char FQ2 = q2 + '0';

   char FileName[17] = {'T','M',FN1,FP1,FQ1,',',' ','T','E',FN2,FP2,FQ2,'.','t','x','t'};

   ofstream out1(FileName,std::ios::app);

   ofstream out2("Messages.txt",std::ios::app); 
 
   int start = MPI_Wtime();

   double Max_TM = Find_Max_B_TM(n1,p1,q1,R1,L1);
   double Max_TE = Find_Max_B_TE(n2,p2,q2,R1,L1);
 
   double omega_p = omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);

   double G_FACTOR = 2.7e-12*sqrt(omega_p/omega_TMTE_p(0,1,0,0,1,1,1.,1.))*sqrt((dist - L1)/0.2)/sqrt(pi*R1*R1*L1);

   MPI_Init(&argn, &argv); 
 
   MPI_Status status;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   for (double m = m_max; m >= m_min; m -= dm)
   {  
      double alpha_cos = alpha_TMTE_cos_p(rank,size,N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m,eps); 
      double alpha_sin = alpha_TMTE_sin_p(rank,size,N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m,eps); 
  
      double global_alpha_cos;
      double global_alpha_sin;

      MPI_Barrier(MPI_COMM_WORLD);
 
      MPI_Reduce(&alpha_cos,&global_alpha_cos,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&alpha_sin,&global_alpha_sin,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

   
      double beta_cos = beta_TMTE_cos_p(rank,size,N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m,eps); 
      double beta_sin = beta_TMTE_sin_p(rank,size,N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m,eps); 
  
      double global_beta_cos;
      double global_beta_sin;

      MPI_Barrier(MPI_COMM_WORLD);
 
      MPI_Reduce(&beta_cos,&global_beta_cos,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&beta_sin,&global_beta_sin,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      double gamma1_cos = gamma1_TMTE_cos_p(rank,size,N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m,eps); 
      double gamma1_sin = gamma1_TMTE_sin_p(rank,size,N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m,eps); 
  
      double global_gamma1_cos;
      double global_gamma1_sin;

      MPI_Barrier(MPI_COMM_WORLD);
 
      MPI_Reduce(&gamma1_cos,&global_gamma1_cos,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&gamma1_sin,&global_gamma1_sin,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      double gamma2_cos = gamma2_TMTE_cos_p(rank,size,N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m,eps); 
      double gamma2_sin = gamma2_TMTE_sin_p(rank,size,N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,dist,N3,P3,Q3,R2,L2,m,eps); 
  
      double global_gamma2_cos;
      double global_gamma2_sin;

      MPI_Barrier(MPI_COMM_WORLD);

      MPI_Reduce(&gamma2_cos,&global_gamma2_cos,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&gamma2_sin,&global_gamma2_sin,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      double Delta = dist - L1;
      double V1 = pi*R1*R1*L1;
      double V2 = pi*R2*R2*L2;  

   
      if (rank == 0)
      {  
         double gcos = global_alpha_cos + (global_beta_cos - global_gamma1_cos + global_gamma2_cos)/(omega_p*omega_p);
         double gsin = global_alpha_sin + (global_beta_sin - global_gamma1_sin + global_gamma2_sin)/(omega_p*omega_p);
         
         double kappa = sqrt(gcos*gcos + gsin*gsin)*Delta/(V1*V2);
         out1 << m*FACTOR_m << " " << kappa << " " << G_FACTOR/sqrt(kappa) << endl;
         out2 << "Coupling_Constant_TMTE_p, m = " << m << endl; 
      }  
   
      MPI_Barrier(MPI_COMM_WORLD);
 
   }

   int finish = MPI_Wtime();

   if (rank == 0)
      cout << "Time = " << (finish - start)/60 << " min " << (finish - start)%60 << " sec." << endl;  



   MPI_Finalize();

   out1.close();
   out2.close();

}







