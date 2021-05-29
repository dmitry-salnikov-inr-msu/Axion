#include <iostream>
#include <fstream>

using namespace std;


double dfunction_TMTE_cos_p(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double m, double z, double zs, double rho, double rhos, double phi, double phis)
{
 
  int n1 = (int) N1;
  int p1 = (int) P1;
  int q1 = (int) Q1;
  int n2 = (int) N2;
  int p2 = (int) P2;
  int q2 = (int) Q2;

  double R = sqrt(rs(z,zs,rho,rhos,phi,phis));
  double omega_p = omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);
     
  if (R == 0)
     return 0;
  else
  {
     if (omega_p > m)
        return rhos*F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*cos(R*sqrt(omega_p*omega_p - m*m))/R;
     else
        return rhos*F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*exp(-R*sqrt(m*m - omega_p*omega_p))/R;
  }
}


int ifunction_TMTE_cos_p(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double *p = ((double *) fdata); 
    fval[0] = dfunction_TMTE_cos_p(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],x[0],p[12],x[1],p[13],x[2]);
    return 0; 
}


double TMTE_cos_p(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double m, double z, double rho, double phi, double eps)
{
   ofstream out("Messages.txt",std::ios::app); 

   double xmin[3] = {0,0,0}, xmax[3] = {L1,R1,2*pi}, p[14] = {N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi}, val, err;
   hcubature(1, ifunction_TMTE_cos_p, &p, 3, xmin, xmax, 0, 0, eps, ERROR_INDIVIDUAL, &val, &err);  
   out << "m = " << m << ", z =  " << z << ", rho =  " << rho << ", phi = " << phi << ", TMTE_cos_p = " << val/(4*pi) << ", err =  " << err/(4*pi) << ", fabs(err/val) = " << fabs(err/val) << endl; 
   
   out.close(); 

   return val/(4*pi); 
}


double dfunction_TMTE_sin_p(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double m, double z, double zs, double rho, double rhos, double phi, double phis)
{
 
  int n1 = (int) N1;
  int p1 = (int) P1;
  int q1 = (int) Q1;
  int n2 = (int) N2;
  int p2 = (int) P2;
  int q2 = (int) Q2;

  double R = sqrt(rs(z,zs,rho,rhos,phi,phis));
  double omega_p = omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);
     
  if (R == 0)
     return 0;
  else
  {
     if (omega_p > m)
        return rhos*F_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*sin(R*sqrt(omega_p*omega_p - m*m))/R;
     else
        return 0.;
  }
}


int ifunction_TMTE_sin_p(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double *p = ((double *) fdata); 
    fval[0] = dfunction_TMTE_sin_p(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],x[0],p[12],x[1],p[13],x[2]);
    return 0; 
}


double TMTE_sin_p(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double m, double z, double rho, double phi, double eps)
{
   ofstream out("Messages.txt",std::ios::app); 

   double xmin[3] = {0,0,0}, xmax[3] = {L1,R1,2*pi}, p[14] = {N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi}, val, err;
   hcubature(1, ifunction_TMTE_sin_p, &p, 3, xmin, xmax, 0, 0, eps, ERROR_INDIVIDUAL, &val, &err);  
   out << "m = " << m << ", z =  " << z << ", rho =  " << rho << ", phi = " << phi << ", TMTE_sin_p = " << val/(4*pi) << ", err =  " << err/(4*pi) << ", fabs(err/val) = " << fabs(err/val) << endl; 

   out.close();
   
   return val/(4*pi);
}



double dfunction_TMTE_cos_m(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double m, double z, double zs, double rho, double rhos, double phi, double phis)
{
 
  int n1 = (int) N1;
  int p1 = (int) P1;
  int q1 = (int) Q1;
  int n2 = (int) N2;
  int p2 = (int) P2;
  int q2 = (int) Q2;

  double R = sqrt(rs(z,zs,rho,rhos,phi,phis));
  double omega_m = omega_TMTE_m(n1,p1,q1,n2,p2,q2,R1,L1);
     
  if (R == 0)
     return 0;
  else
  {
     if (omega_m > m)
        return rhos*F_TMTE_m(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*cos(R*sqrt(omega_m*omega_m - m*m))/R;
     else
        return rhos*F_TMTE_m(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*exp(-R*sqrt(m*m - omega_m*omega_m))/R;
  }
}


int ifunction_TMTE_cos_m(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double *p = ((double *) fdata); 
    fval[0] = dfunction_TMTE_cos_m(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],x[0],p[12],x[1],p[13],x[2]);
    return 0; 
}


double TMTE_cos_m(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double m, double z, double rho, double phi, double eps)
{
   ofstream out("Messages.txt",std::ios::app); 

   double xmin[3] = {0,0,0}, xmax[3] = {L1,R1,2*pi}, p[14] = {N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi}, val, err;
   hcubature(1, ifunction_TMTE_cos_m, &p, 3, xmin, xmax, 0, 0, eps, ERROR_INDIVIDUAL, &val, &err);  
   out << "m = " << m << ", z =  " << z << ", rho =  " << rho << ", phi = " << phi << ", TMTE_cos_m = " << val/(4*pi) << ", err =  " << err/(4*pi) << ", fabs(err/val) = " << fabs(err/val) << endl; 

   out.close();  

   return val/(4*pi);
}



double dfunction_TMTE_sin_m(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double m, double z, double zs, double rho, double rhos, double phi, double phis)
{
 
  int n1 = (int) N1;
  int p1 = (int) P1;
  int q1 = (int) Q1;
  int n2 = (int) N2;
  int p2 = (int) P2;
  int q2 = (int) Q2;

  double R = sqrt(rs(z,zs,rho,rhos,phi,phis));
  double omega_m = omega_TMTE_m(n1,p1,q1,n2,p2,q2,R1,L1);
     
  if (R == 0)
     return 0;
  else
  {
     if (omega_m > m)
        return rhos*F_TMTE_m(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,zs,rhos,phis)*sin(R*sqrt(omega_m*omega_m - m*m))/R;
     else
        return 0.;
  }
}


int ifunction_TMTE_sin_m(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{
    double *p = ((double *) fdata); 
    fval[0] = dfunction_TMTE_sin_m(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],x[0],p[12],x[1],p[13],x[2]);
    return 0; 
}


double TMTE_sin_m(double N1, double P1, double Q1, double N2, double P2, double Q2, double R1, double L1, double Max_TM, double Max_TE, double m, double z, double rho, double phi, double eps)
{

   ofstream out("Messages.txt",std::ios::app); 

   double xmin[3] = {0,0,0}, xmax[3] = {L1,R1,2*pi}, p[14] = {N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi}, val, err;
   hcubature(1, ifunction_TMTE_sin_m, &p, 3, xmin, xmax, 0, 0, eps, ERROR_INDIVIDUAL, &val, &err);  
   out << "m = " << m << ", z =  " << z << ", rho =  " << rho << ", phi = " << phi << ", TMTE_sin_m = " << val/(4*pi) << ", err =  " << err/(4*pi) << ", fabs(err/val) = " << fabs(err/val) << endl; 
 
   out.close();
   
   return val/(4*pi);
}


const double h = 1e-3;
const double FACTOR_E = 1.2e-2;
const double FACTOR_m = 2e-7;


double Axion_Energy_Density_TMTE_p(int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double m, double z, double rho, double phi, double eps)
{
   double N1 = (double) n1;
   double P1 = (double) p1;
   double Q1 = (double) q1;
   double N2 = (double) n2;
   double P2 = (double) p2;
   double Q2 = (double) q2;
 
   double omega_p = omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);

   double aCp = TMTE_cos_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi,eps);
   double aSp = TMTE_sin_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi,eps);

   double drho = h*R1;
   double dphi = h*(2*pi);
   double dz   = h*L1;

   double aCp_rho_p = TMTE_cos_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho + drho,phi,eps);
   double aCp_rho_m = TMTE_cos_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho - drho,phi,eps);
  
   double aCp_phi_p = TMTE_cos_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi + dphi,eps);
   double aCp_phi_m = TMTE_cos_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi - dphi,eps);
 
   double aCp_z_p   = TMTE_cos_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z + dz,rho,phi,eps);
   double aCp_z_m   = TMTE_cos_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z - dz,rho,phi,eps);
 
   double partial_rho_aCp = (aCp_rho_p - aCp_rho_m)/(2*drho);
   double partial_phi_aCp = (aCp_phi_p - aCp_phi_m)/(2*dphi);
   double partial_z_aCp   = (aCp_z_p - aCp_z_m)/(2*dz);
 
   double aSp_rho_p = TMTE_sin_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho + drho,phi,eps);
   double aSp_rho_m = TMTE_sin_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho - drho,phi,eps);
  
   double aSp_phi_p = TMTE_sin_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi + dphi,eps);
   double aSp_phi_m = TMTE_sin_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi - dphi,eps);
 
   double aSp_z_p   = TMTE_sin_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z + dz,rho,phi,eps);
   double aSp_z_m   = TMTE_sin_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z - dz,rho,phi,eps);
 
   double partial_rho_aSp = (aSp_rho_p - aSp_rho_m)/(2*drho);
   double partial_phi_aSp = (aSp_phi_p - aSp_phi_m)/(2*dphi);
   double partial_z_aSp   = (aSp_z_p - aSp_z_m)/(2*dz);
   
   double nabla_aCp_sqr = partial_rho_aCp*partial_rho_aCp + partial_phi_aCp*partial_phi_aCp/(rho*rho) + partial_z_aCp*partial_z_aCp;
   double nabla_aSp_sqr = partial_rho_aSp*partial_rho_aSp + partial_phi_aSp*partial_phi_aSp/(rho*rho) + partial_z_aSp*partial_z_aSp; 

   ofstream out("Messages.txt",std::ios::app); 
 
   out << "Axion_Energy_Density_TMTE_p, m = " << m << ", z = " << z << ", rho = " << rho << ", phi = " << phi << endl; 

   out.close();

   double rho_E = (m*m + omega_p*omega_p)*(aCp*aCp + aSp*aSp)/4  +  (nabla_aCp_sqr + nabla_aSp_sqr)/4;

   return rho_E*FACTOR_E;

}


double Axion_Energy_Density_TMTE_m(int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double m, double z, double rho, double phi, double eps)
{
   double N1 = (double) n1;
   double P1 = (double) p1;
   double Q1 = (double) q1;
   double N2 = (double) n2;
   double P2 = (double) p2;
   double Q2 = (double) q2;
 
   double omega_m = omega_TMTE_m(n1,p1,q1,n2,p2,q2,R1,L1);

   double aCm = TMTE_cos_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi,eps);
   double aSm = TMTE_sin_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi,eps);

   double drho = h*R1;
   double dphi = h*(2*pi);
   double dz   = h*L1;

   double aCm_rho_p = TMTE_cos_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho + drho,phi,eps);
   double aCm_rho_m = TMTE_cos_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho - drho,phi,eps);
  
   double aCm_phi_p = TMTE_cos_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi + dphi,eps);
   double aCm_phi_m = TMTE_cos_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi - dphi,eps);
 
   double aCm_z_p   = TMTE_cos_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z + dz,rho,phi,eps);
   double aCm_z_m   = TMTE_cos_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z - dz,rho,phi,eps);
 
   double partial_rho_aCm = (aCm_rho_p - aCm_rho_m)/(2*drho);
   double partial_phi_aCm = (aCm_phi_p - aCm_phi_m)/(2*dphi);
   double partial_z_aCm   = (aCm_z_p - aCm_z_m)/(2*dz);
 
   double aSm_rho_p = TMTE_sin_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho + drho,phi,eps);
   double aSm_rho_m = TMTE_sin_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho - drho,phi,eps);
  
   double aSm_phi_p = TMTE_sin_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi + dphi,eps);
   double aSm_phi_m = TMTE_sin_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z,rho,phi - dphi,eps);
 
   double aSm_z_p   = TMTE_sin_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z + dz,rho,phi,eps);
   double aSm_z_m   = TMTE_sin_m(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,z - dz,rho,phi,eps);
 
   double partial_rho_aSm = (aSm_rho_p - aSm_rho_m)/(2*drho);
   double partial_phi_aSm = (aSm_phi_p - aSm_phi_m)/(2*dphi);
   double partial_z_aSm   = (aSm_z_p - aSm_z_m)/(2*dz);
   
   double nabla_aCm_sqr = partial_rho_aCm*partial_rho_aCm + partial_phi_aCm*partial_phi_aCm/(rho*rho) + partial_z_aCm*partial_z_aCm;
   double nabla_aSm_sqr = partial_rho_aSm*partial_rho_aSm + partial_phi_aSm*partial_phi_aSm/(rho*rho) + partial_z_aSm*partial_z_aSm; 
 
   ofstream out("Messages.txt",std::ios::app); 
 
   out << "Axion_Energy_Density_TMTE_m, m = " << m << ", z = " << z << ", rho = " << rho << ", phi = " << phi << endl; 

   out.close();


   double rho_E = (m*m + omega_m*omega_m)*(aCm*aCm + aSm*aSm)/4  +  (nabla_aCm_sqr + nabla_aSm_sqr)/4;

   return rho_E*FACTOR_E;

}


void Distribution_of_Axion_Field_RhoZ_TMTE_cos_p(int argn, char** argv, int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double m, double z_min, double z_max, double rho_min, double rho_max, double phi, double eps)
{

    
   double Max_TM = Find_Max_B_TM(n1,p1,q1,R1,L1);
   double Max_TE = Find_Max_B_TE(n2,p2,q2,R1,L1);
 
   double rho;
   double z;
 
   ofstream out("Distribution of Axion field cos (rho,z) TMTE omega_p.txt",std::ios::app);
    
   MPI_Init(&argn, &argv); 

   int start = MPI_Wtime();
 
   MPI_Status status;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   double dz = (z_max - z_min)/(size - 1);  

   z = z_min + dz*((double) rank);
 
   for (rho = rho_min; rho <= (rho_max + 0.01); rho += 0.05)
   {
      out << rho << " " << z << " " << TMTE_cos_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,m,z + 0.005,rho + 0.005,phi,eps) << endl;
   }
 
   MPI_Barrier(MPI_COMM_WORLD);

   int finish = MPI_Wtime();
   
   if (rank == 0)
      cout << "Time = " << (finish - start)/60 << " min " << (finish - start)%60 << " sec." << endl;  


   MPI_Finalize();

   out.close();

}


void Distribution_of_Axion_Field_RhoZ_TMTE_sin_p(int argn, char** argv, int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double m, double z_min, double z_max, double rho_min, double rho_max, double phi, double eps)
{
 
   double Max_TM = Find_Max_B_TM(n1,p1,q1,R1,L1);
   double Max_TE = Find_Max_B_TE(n2,p2,q2,R1,L1);
 
   double rho;
   double z;

   ofstream out("Distribution of Axion field sin (rho,z) TMTE omega_p.txt",std::ios::app);
     
   MPI_Init(&argn, &argv); 
 
   int start = MPI_Wtime();
    
   MPI_Status status;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   double dz = (z_max - z_min)/(size - 1);  

   z = z_min + dz*((double) rank);
 
   for (rho = rho_min; rho <= (rho_max + 0.01); rho += 0.05)
   {
      out << rho << " " << z << " " << TMTE_sin_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,m,z + 0.005,rho + 0.005,phi,eps) << endl;
   }
 
   MPI_Barrier(MPI_COMM_WORLD);

   int finish = MPI_Wtime();
   
   if (rank == 0)
      cout << "Time = " << (finish - start)/60 << " min " << (finish - start)%60 << " sec." << endl;  


 
   MPI_Finalize();


   out.close();

}

void Distribution_of_Axion_Field_RhoZ_TMTE_cos_m(int argn, char** argv, int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double m, double z_min, double z_max, double rho_min, double rho_max, double phi, double eps)
{
    
   double Max_TM = Find_Max_B_TM(n1,p1,q1,R1,L1);
   double Max_TE = Find_Max_B_TE(n2,p2,q2,R1,L1);
 
   double rho;
   double z;
 
   ofstream out("Distribution of Axion field cos (rho,z) TMTE omega_m.txt",std::ios::app);
    
   MPI_Init(&argn, &argv); 

   int start = MPI_Wtime();
  
   MPI_Status status;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   double dz = (z_max - z_min)/(size - 1);  

   z = z_min + dz*((double) rank);
 
   for (rho = rho_min; rho <= (rho_max + 0.01); rho += 0.05)
   {
      out << rho << " " << z << " " << TMTE_cos_m(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,m,z + 0.005,rho + 0.005,phi,eps) << endl;
   }
 
   MPI_Barrier(MPI_COMM_WORLD);
 
   int finish = MPI_Wtime();
   
   if (rank == 0)
      cout << "Time = " << (finish - start)/60 << " min " << (finish - start)%60 << " sec." << endl;  


   MPI_Finalize();

   out.close();

}


void Distribution_of_Axion_Field_RhoZ_TMTE_sin_m(int argn, char** argv, int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double m, double z_min, double z_max, double rho_min, double rho_max, double phi, double eps)
{
 
   double Max_TM = Find_Max_B_TM(n1,p1,q1,R1,L1);
   double Max_TE = Find_Max_B_TE(n2,p2,q2,R1,L1);
 
   double rho;
   double z;

   ofstream out("Distribution of Axion field sin (rho,z) TMTE omega_m.txt",std::ios::app);
     
   MPI_Init(&argn, &argv); 

   int start = MPI_Wtime();
 
   MPI_Status status;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   double dz = (z_max - z_min)/(size - 1);  

   z = z_min + dz*((double) rank);
 
   for (rho = rho_min; rho <= (rho_max + 0.01); rho += 0.05)
   {
      out << rho << " " << z << " " << TMTE_sin_m(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,m,z + 0.005,rho + 0.005,phi,eps) << endl;
   }
 

   int finish = MPI_Wtime();
   
   if (rank == 0)
      cout << "Time = " << (finish - start)/60 << " min " << (finish - start)%60 << " sec." << endl;  


 
   MPI_Finalize();


   out.close();

}


void Distribution_of_Axion_Energy_Density_MZ_TMTE_p(int argn, char** argv, int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double m_min, double m_max, double z_min, double z_max, double rho, double phi, double eps)
{
   double Max_TM = Find_Max_B_TM(n1,p1,q1,R1,L1);
   double Max_TE = Find_Max_B_TE(n2,p2,q2,R1,L1);
 
 
   ofstream out("Distribution of Axion Energy Density (m,z) TMTE omega_p.txt",std::ios::app);
    
   MPI_Init(&argn, &argv); 
  
   int start = MPI_Wtime();

 
   MPI_Status status;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   double dm = (m_max - m_min)/(size - 1);

   double m = m_min + dm*((double) rank);

   for (double z = z_min; z <= (z_max + 0.01); z += 0.05)
      out << m*FACTOR_m << " " << z << " " << Axion_Energy_Density_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,m,z + 0.005,rho + 0.005,phi,eps) << endl;  


   MPI_Barrier(MPI_COMM_WORLD);
  
   int finish = MPI_Wtime();
 
   if (rank == 0)
      cout << "Time = " << (finish - start)/60 << " min " << (finish - start)%60 << " sec." << endl;  


   MPI_Finalize();


   out.close();

}


void Distribution_of_Axion_Energy_Density_RhoZ_TMTE_p(int argn, char** argv, int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double m, double z_min, double z_max, double rho_min, double rho_max, double phi, double eps)
{
   double Max_TM = Find_Max_B_TM(n1,p1,q1,R1,L1);
   double Max_TE = Find_Max_B_TE(n2,p2,q2,R1,L1);
 
   ofstream out("Distribution of Axion Energy Density (rho,z) TMTE omega_p.txt",std::ios::app);
    
   MPI_Init(&argn, &argv); 
 
   int start = MPI_Wtime();

   MPI_Status status;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   double dz = (z_max - z_min)/(size - 1);  

   double z = z_min + dz*((double) rank);
     
   for (double rho = rho_min; rho <= (rho_max + 0.01); rho += 0.05)
   {
      out << rho << " " << z << " " << Axion_Energy_Density_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,m,z + 0.005,rho + 0.005,phi,eps) << endl;
   }
 
   MPI_Barrier(MPI_COMM_WORLD);

   int finish = MPI_Wtime();
   
   if (rank == 0)
      cout << "Time = " << (finish - start)/60 << " min " << (finish - start)%60 << " sec." << endl;  

  
   MPI_Finalize();


   out.close();

}



void Distribution_of_Axion_Energy_Density_RhoPhi_TMTE_p(int argn, char** argv, int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double m, double z, double rho_min, double rho_max, double eps)
{
   double Max_TM = Find_Max_B_TM(n1,p1,q1,R1,L1);
   double Max_TE = Find_Max_B_TE(n2,p2,q2,R1,L1);

   ofstream out("Distribution of Axion Energy Density (rho, phi) TMTE omega_p.txt",std::ios::app);
    
   MPI_Init(&argn, &argv); 
 
   int start = MPI_Wtime();

   MPI_Status status;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   double phi = 2*pi*((double)rank)/((double) size);
  
   if (rank == 0)
   {
      for (double rho = rho_min; rho <= (rho_max + 0.01); rho += 0.1)
      {
         out << phi << " " << rho << " " <<  Axion_Energy_Density_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,m,z + 0.005,rho + 0.005,phi,eps)  << endl;
      }
   }
   else
   {
      for (double rho = (rho_min + 0.1); rho <= (rho_max + 0.01); rho += 0.1)
      {
         out << phi << " " << rho << " " <<  Axion_Energy_Density_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,m,z + 0.005,rho + 0.005,phi,eps)   << endl;
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);
  
   int finish = MPI_Wtime();
 
   if (rank == 0)
      cout << "Time = " << (finish - start)/60 << " min " << (finish - start)%60 << " sec." << endl;  


   MPI_Finalize();

   out.close();

}


void Distribution_of_Axion_Energy_Density_MZ_TMTE_m(int argn, char** argv, int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double m_min, double m_max, double z_min, double z_max, double rho, double phi, double eps)
{
   double Max_TM = Find_Max_B_TM(n1,p1,q1,R1,L1);
   double Max_TE = Find_Max_B_TE(n2,p2,q2,R1,L1);
 
 
   ofstream out("Distribution of Axion Energy Density (m,z) TMTE omega_m.txt",std::ios::app);
    
   MPI_Init(&argn, &argv); 
  
   int start = MPI_Wtime();

 
   MPI_Status status;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   double dm = (m_max - m_min)/(size - 1);

   double m = m_min + dm*((double) rank);

   for (double z = z_min; z <= (z_max + 0.01); z += 0.05)
      out << m*FACTOR_m << " " << z << " " << Axion_Energy_Density_TMTE_m(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,m,z + 0.005,rho + 0.005,phi,eps) << endl;  


   MPI_Barrier(MPI_COMM_WORLD);
  
   int finish = MPI_Wtime();
 
   if (rank == 0)
      cout << "Time = " << (finish - start)/60 << " min " << (finish - start)%60 << " sec." << endl;  



   MPI_Finalize();


   out.close();

}


void Distribution_of_Axion_Energy_Density_RhoZ_TMTE_m(int argn, char** argv, int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double m, double z_min, double z_max, double rho_min, double rho_max, double phi, double eps)
{
   double Max_TM = Find_Max_B_TM(n1,p1,q1,R1,L1);
   double Max_TE = Find_Max_B_TE(n2,p2,q2,R1,L1);
 
   ofstream out("Distribution of Axion Energy Density (rho,z) TMTE omega_m.txt",std::ios::app);
    
   MPI_Init(&argn, &argv); 
 
   int start = MPI_Wtime();

   MPI_Status status;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   double dz = (z_max - z_min)/(size - 1);  

   double z = z_min + dz*((double) rank);
     
   for (double rho = rho_min; rho <= (rho_max + 0.01); rho += 0.05)
   {
      out << rho << " " << z << " " << Axion_Energy_Density_TMTE_m(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,m,z + 0.005,rho + 0.005,phi,eps) << endl;
   }
 
   MPI_Barrier(MPI_COMM_WORLD);
  
   int finish = MPI_Wtime();
 
   if (rank == 0)
      cout << "Time = " << (finish - start)/60 << " min " << (finish - start)%60 << " sec." << endl;  

  
   MPI_Finalize();


   out.close();

}



void Distribution_of_Axion_Energy_Density_RhoPhi_TMTE_m(int argn, char** argv, int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double m, double z, double rho_min, double rho_max, double eps)
{
   double Max_TM = Find_Max_B_TM(n1,p1,q1,R1,L1);
   double Max_TE = Find_Max_B_TE(n2,p2,q2,R1,L1);

   ofstream out("Distribution of Axion Energy Density (rho, phi) TMTE omega_m.txt",std::ios::app);
    
   MPI_Init(&argn, &argv); 
 
   int start = MPI_Wtime();

   MPI_Status status;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   double phi = 2*pi*((double)rank)/((double) size);
  
   if (rank == 0)
   {
      for (double rho = rho_min; rho <= (rho_max + 0.01); rho += 0.1)
      {
         out << phi << " " << rho << " " <<  Axion_Energy_Density_TMTE_m(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,m,z + 0.005,rho + 0.005,phi,eps)  << endl;
      }
   }
   else
   {
      for (double rho = (rho_min + 0.1); rho <= (rho_max + 0.01); rho += 0.1)
      {
         out << phi << " " << rho << " " <<  Axion_Energy_Density_TMTE_m(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,m,z + 0.005,rho + 0.005,phi,eps)   << endl;
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);
  
   int finish = MPI_Wtime();
 
   if (rank == 0)
      cout << "Time = " << (finish - start)/60 << " min " << (finish - start)%60 << " sec." << endl;  


   MPI_Finalize();

   out.close();

}




 
