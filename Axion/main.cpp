#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include "cubature.h"
#include "TMnpq,TEnpq.h"
#include "energy_density.h"
#include "coupling_constant.h"
#include "radiation_pattern.h"

using namespace std;


int main (int argn, char** argv)
{
   int n1, p1, q1,
       n2, p2, q2,
       n3, p3, q3;

   double R1, L1, R2, L2; 
 
   double eps;

   double m, m_min, m_max, dm,
          z, z_min, z_max,
          rho, rho_min, rho_max,
          phi, dist;

   /*Distribution of Axion Field*/

   n1 = 0; p1 = 1; q1 = 0;
   n2 = 0; p2 = 1; q2 = 1;
 
   R1 = 1.; L1 = 1.;

   eps = 1e-3;

   m = 5.;

   z_min = 0.5*L1 + 0.;
   z_max = 0.5*L1 + 2.5;

   rho_min = 0.;
   rho_max = 5.;

   phi = 0.;

//   Distribution_of_Axion_Field_RhoZ_TMTE_cos_p(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z_min,z_max,rho_min,rho_max,phi,eps);
//   Distribution_of_Axion_Field_RhoZ_TMTE_sin_p(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z_min,z_max,rho_min,rho_max,phi,eps);
 
//   Distribution_of_Axion_Field_RhoZ_TMTE_cos_m(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z_min,z_max,rho_min,rho_max,phi,eps);
//   Distribution_of_Axion_Field_RhoZ_TMTE_sin_m(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z_min,z_max,rho_min,rho_max,phi,eps);


   /*Distribution of Axion Energy Density*/

   n1 = 0; p1 = 1; q1 = 0;
   n2 = 0; p2 = 1; q2 = 1;
 
   R1 = 1.; L1 = 1.;

   eps = 1e-5;

   m_min = 0.;
   m_max = 10.;

   z_min = 0.5*L1;
   z_max = 0.5*L1 + 1.;

   rho = 0.;
 
   phi = 0.;

//   Distribution_of_Axion_Energy_Density_MZ_TMTE_p(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m_min,m_max,z_min,z_max,rho,phi,eps);
 
   m = 7.;

   z_min = 0.5*L1;
   z_max = 0.5*L1 + 1.;

   rho_min = 0.;
   rho_max = 1.4; 

   phi = 0.;
  
//   Distribution_of_Axion_Energy_Density_RhoZ_TMTE_p(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z_min,z_max,rho_min,rho_max,phi,eps);

   m = 7.;

   z = 0.5*L1;

   rho_min = 0.;
   rho_max = 1.6; 
  
//   Distribution_of_Axion_Energy_Density_RhoPhi_TMTE_p(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z,rho_min,rho_max,eps);

   m_min = 0.;
   m_max = 10.;

   z_min = 0.5*L1;
   z_max = 0.5*L1 + 1.;

   rho = 0.;
 
   phi = 0.;

//   Distribution_of_Axion_Energy_Density_MZ_TMTE_m(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m_min,m_max,z_min,z_max,rho,phi,eps);
 
   m = 0.;

   z_min = 0.5*L1;
   z_max = 0.5*L1 + 1.;

   rho_min = 0.;
   rho_max = 1.4; 

   phi = 0.;
  
//   Distribution_of_Axion_Energy_Density_RhoZ_TMTE_m(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z_min,z_max,rho_min,rho_max,phi,eps);

   m = 0.;

   z = 0.5*L1;

   rho_min = 0.;
   rho_max = 1.4; 
  
//   Distribution_of_Axion_Energy_Density_RhoPhi_TMTE_m(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z,rho_min,rho_max,eps);

 
   /*Coupling Constant*/

   n1 = 0; p1 = 1; q1 = 0;
   n2 = 0; p2 = 1; q2 = 1;

   R1 = 1.; L1 = 1.;

   n3 = 0; p3 = 1; q3 = 0;
    
   eps = 1e-3;

   L2 = 0.2;
   R2 = x_np(n3,p3)/omega_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1);

   dist = L1 + 0.2;
   dm = 0.1;   
   m_min = 0.;
   m_max = 10.; 

//   Coupling_Constant_TMTE_p(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,dist,n3,p3,q3,R2,L2,m_min,m_max,dm,eps);


   /*Radiation Pattern*/

   n1 = 0; p1 = 1; q1 = 0;
   n2 = 0; p2 = 1; q2 = 1;

   phi = 0.;

   R1 = 1.;
   L1 = 1.;

   m = 5.;
   
   eps = 1e-3;   

//   Radiation_Pattern_Theta_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,m,phi,eps);

   return 0;
}



