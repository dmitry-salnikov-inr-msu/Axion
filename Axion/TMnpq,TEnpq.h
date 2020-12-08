#include <cmath>

double pi = 4*atan(1);

double x_np(int n, int p)
{
   double x[3][5] = {{0., 2.405, 5.520, 8.654 , 11.792},
                     {0., 3.832, 7.016, 10.174, 13.324},
                     {0., 5.135, 8.417, 11.620, 14.796}};

   return x[n][p];   
}

double xs_np(int n, int p)
{
   double x[3][5] = {{0., 3.832, 7.016, 10.174, 13.324},
                     {0., 1.841, 5.331, 8.536 , 11.706},
                     {0., 3.054, 6.706, 9.970 , 13.170}};

   return x[n][p];   
}

double k_sqr(int n, int p, int q, double R1, double L1)
{
   double Q = (double) q;
   double X = x_np(n,p);
   return (Q*pi/L1)*(Q*pi/L1) + (X/R1)*(X/R1);
}

double omega_TM(int n, int p, int q, double R1, double L1)
{
   return sqrt(k_sqr(n,p,q,R1,L1));
}

double ks_sqr(int n, int p, int q, double R1, double L1)
{
   double Q = (double) q;
   double X = xs_np(n,p);
   return (Q*pi/L1)*(Q*pi/L1) + (X/R1)*(X/R1);
}

double omega_TE(int n, int p, int q, double R1, double L1)
{
   return sqrt(ks_sqr(n,p,q,R1,L1));
}


double omega_TMTE_p(int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1)
{
   return fabs(omega_TM(n1,p1,q1,R1,L1) + omega_TE(n2,p2,q2,R1,L1));
}

double omega_TMTE_m(int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1)
{
   return fabs(omega_TM(n1,p1,q1,R1,L1) - omega_TE(n2,p2,q2,R1,L1));
}


/* differential of Bessel function  */
double jns(int n, double x)
{
   return 0.5*(jn(n - 1,x) - jn(n + 1,x));
}


/* TM_npq modes    */

double modeTM_Ez (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return jn(n,X*rho/R1)*cos(N*phi)*cos(Q*pi*z/L1);
}

double modeTM_Erho (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double K_sqr = k_sqr(n,p,q,R1,L1);
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return (-1./(K_sqr - (Q*pi/L1)*(Q*pi/L1)))*(Q*pi/L1)*(X/R1)*jns(n,X*rho/R1)*cos(N*phi)*sin(Q*pi*z/L1);
}

double modeTM_Ephi (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double K_sqr = k_sqr(n,p,q,R1,L1);
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   if (rho == 0)
      return 0;
   else
      return (-1./(K_sqr - (Q*pi/L1)*(Q*pi/L1)))*(1./rho)*(N*Q*pi/L1)*jn(n,X*rho/R1)*(-sin(N*phi))*sin(Q*pi*z/L1);
}

/* modeTM_Hz = 0   by definition  */

double modeTM_Hrho (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double W = omega_TM(n,p,q,R1,L1);
   double K_sqr = k_sqr(n,p,q,R1,L1);
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;


   if (rho == 0)
      return 0;
   else
      return (-W/(K_sqr - (Q*pi/L1)*(Q*pi/L1)))*(N/rho)*jn(n,X*rho/R1)*(-sin(N*phi))*cos(Q*pi*z/L1);
}

double modeTM_Hphi (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double W = omega_TM(n,p,q,R1,L1);
   double K_sqr = k_sqr(n,p,q,R1,L1);
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return (W/(K_sqr - (Q*pi/L1)*(Q*pi/L1)))*(X/R1)*jns(n,X*rho/R1)*cos(N*phi)*cos(Q*pi*z/L1);
}


/* TE_npq modes    */

double modeTE_Hz (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return jn(n,X*rho/R1)*cos(N*phi)*sin(Q*pi*z/L1);
}

double modeTE_Hrho (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double Ks_sqr = ks_sqr(n,p,q,R1,L1);
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return (1./(Ks_sqr-(Q*pi/L1)*(Q*pi/L1)))*(Q*pi/L1)*(X/R1)*jns(n,X*rho/R1)*cos(N*phi)*cos(Q*pi*z/L1);
}

double modeTE_Hphi (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double Ks_sqr = ks_sqr(n,p,q,R1,L1);
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   if (rho == 0)
      return 0;
   else
      return (1./(Ks_sqr - (Q*pi/L1)*(Q*pi/L1)))*(Q*pi/L1)*(N/rho)*jn(n,X*rho/R1)*(-sin(N*phi))*cos(Q*pi*z/L1);
}

/* modeTE_Ez = 0   by definition  */

double modeTE_Erho (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double W = omega_TE(n,p,q,R1,L1);
   double Ks_sqr = ks_sqr(n,p,q,R1,L1);
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;


   if (rho == 0)
      return 0;
   else
      return (W/(Ks_sqr - (Q*pi/L1)*(Q*pi/L1)))*(N/rho)*jn(n,X*rho/R1)*(-sin(N*phi))*sin(Q*pi*z/L1);
}

double modeTE_Ephi (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double W = omega_TE(n,p,q,R1,L1);
   double Ks_sqr = ks_sqr(n,p,q,R1,L1);
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;


   return -(W/(Ks_sqr -(Q*pi/L1)*(Q*pi/L1)))*(X/R1)*jns(n,X*rho/R1)*cos(N*phi)*sin(Q*pi*z/L1);
}


/* Detection TM_npq mode  */

double Detection_modeTM_Ez (int n, int p, int q, double dist, double R2, double L2, double z, double rho, double phi)
{
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return jn(n,X*rho/R2)*cos(N*phi)*cos(Q*pi*(z - dist)/L2);
}


/* Partial_z TM_npq modes    */

double partial_z_modeTM_Ez (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return jn(n,X*rho/R1)*cos(N*phi)*(Q*pi/L1)*(-sin(Q*pi*z/L1));
}

double partial_z_modeTM_Erho (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double K_sqr = k_sqr(n,p,q,R1,L1);
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return (-1./(K_sqr - (Q*pi/L1)*(Q*pi/L1)))*(Q*pi/L1)*(X/R1)*jns(n,X*rho/R1)*cos(N*phi)*(Q*pi/L1)*cos(Q*pi*z/L1);
}

double partial_z_modeTM_Ephi (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double K_sqr = k_sqr(n,p,q,R1,L1);
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   if (rho == 0)
      return 0;
   else
      return (-1./(K_sqr - (Q*pi/L1)*(Q*pi/L1)))*(1./rho)*(N*Q*pi/L1)*jn(n,X*rho/R1)*(-sin(N*phi))*(Q*pi/L1)*cos(Q*pi*z/L1);
}

/* modeTM_Hz = 0   by definition  */

double partial_z_modeTM_Hrho (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double W = omega_TM(n,p,q,R1,L1);
   double K_sqr = k_sqr(n,p,q,R1,L1);
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;


   if (rho == 0)
      return 0;
   else
      return (-W/(K_sqr - (Q*pi/L1)*(Q*pi/L1)))*(N/rho)*jn(n,X*rho/R1)*(-sin(N*phi))*(Q*pi/L1)*(-sin(Q*pi*z/L1));
}

double partial_z_modeTM_Hphi (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double W = omega_TM(n,p,q,R1,L1);
   double K_sqr = k_sqr(n,p,q,R1,L1);
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return (W/(K_sqr - (Q*pi/L1)*(Q*pi/L1)))*(X/R1)*jns(n,X*rho/R1)*cos(N*phi)*(Q*pi/L1)*(-sin(Q*pi*z/L1));
}


/* Partial_z TE_npq modes    */

double partial_z_modeTE_Hz (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return jn(n,X*rho/R1)*cos(N*phi)*(Q*pi/L1)*cos(Q*pi*z/L1);
}

double partial_z_modeTE_Hrho (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double Ks_sqr = ks_sqr(n,p,q,R1,L1);
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return (1./(Ks_sqr-(Q*pi/L1)*(Q*pi/L1)))*(Q*pi/L1)*(X/R1)*jns(n,X*rho/R1)*cos(N*phi)*(Q*pi/L1)*(-sin(Q*pi*z/L1));
}

double partial_z_modeTE_Hphi (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double Ks_sqr = ks_sqr(n,p,q,R1,L1);
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   if (rho == 0)
      return 0;
   else
      return (1./(Ks_sqr - (Q*pi/L1)*(Q*pi/L1)))*(Q*pi/L1)*(N/rho)*jn(n,X*rho/R1)*(-sin(N*phi))*(Q*pi/L1)*(-sin(Q*pi*z/L1));
}

/* modeTE_Ez = 0   by definition  */

double partial_z_modeTE_Erho (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double W = omega_TE(n,p,q,R1,L1);
   double Ks_sqr = ks_sqr(n,p,q,R1,L1);
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;


   if (rho == 0)
      return 0;
   else
      return (W/(Ks_sqr - (Q*pi/L1)*(Q*pi/L1)))*(N/rho)*jn(n,X*rho/R1)*(-sin(N*phi))*(Q*pi/L1)*cos(Q*pi*z/L1);
}

double partial_z_modeTE_Ephi (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double W = omega_TE(n,p,q,R1,L1);
   double Ks_sqr = ks_sqr(n,p,q,R1,L1);
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;


   return -(W/(Ks_sqr - (Q*pi/L1)*(Q*pi/L1)))*(X/R1)*jns(n,X*rho/R1)*cos(N*phi)*(Q*pi/L1)*cos(Q*pi*z/L1);
}



/* Partial2_z TM_npq modes    */

double partial2_z_modeTM_Ez (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return jn(n,X*rho/R1)*cos(N*phi)*(Q*pi/L1)*(Q*pi/L1)*(-cos(Q*pi*z/L1));
}

double partial2_z_modeTM_Erho (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double K_sqr = k_sqr(n,p,q,R1,L1);
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return (-1./(K_sqr - (Q*pi/L1)*(Q*pi/L1)))*(Q*pi/L1)*(X/R1)*jns(n,X*rho/R1)*cos(N*phi)*(Q*pi/L1)*(Q*pi/L1)*(-sin(Q*pi*z/L1));
}

double partial2_z_modeTM_Ephi (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double K_sqr = k_sqr(n,p,q,R1,L1);
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   if (rho == 0)
      return 0;
   else
      return (-1./(K_sqr - (Q*pi/L1)*(Q*pi/L1)))*(1./rho)*(N*Q*pi/L1)*jn(n,X*rho/R1)*(-sin(N*phi))*(Q*pi/L1)*(Q*pi/L1)*(-sin(Q*pi*z/L1));
}

/* modeTM_Hz = 0   by definition  */

double partial2_z_modeTM_Hrho (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double W = omega_TM(n,p,q,R1,L1);
   double K_sqr = k_sqr(n,p,q,R1,L1);
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;


   if (rho == 0)
      return 0;
   else
      return (-W/(K_sqr - (Q*pi/L1)*(Q*pi/L1)))*(N/rho)*jn(n,X*rho/R1)*(-sin(N*phi))*(Q*pi/L1)*(Q*pi/L1)*(-cos(Q*pi*z/L1));
}

double partial2_z_modeTM_Hphi (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double W = omega_TM(n,p,q,R1,L1);
   double K_sqr = k_sqr(n,p,q,R1,L1);
   double X = x_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return (W/(K_sqr - (Q*pi/L1)*(Q*pi/L1)))*(X/R1)*jns(n,X*rho/R1)*cos(N*phi)*(Q*pi/L1)*(Q*pi/L1)*(-cos(Q*pi*z/L1));
}


/* Partial2_z TE_npq modes    */

double partial2_z_modeTE_Hz (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return jn(n,X*rho/R1)*cos(N*phi)*(Q*pi/L1)*(Q*pi/L1)*(-sin(Q*pi*z/L1));
}

double partial2_z_modeTE_Hrho (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double Ks_sqr = ks_sqr(n,p,q,R1,L1);
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   return (1./(Ks_sqr-(Q*pi/L1)*(Q*pi/L1)))*(Q*pi/L1)*(X/R1)*jns(n,X*rho/R1)*cos(N*phi)*(Q*pi/L1)*(Q*pi/L1)*(-cos(Q*pi*z/L1));
}

double partial2_z_modeTE_Hphi (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double Ks_sqr = ks_sqr(n,p,q,R1,L1);
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;

   if (rho == 0)
      return 0;
   else
      return (1./(Ks_sqr - (Q*pi/L1)*(Q*pi/L1)))*(Q*pi/L1)*(N/rho)*jn(n,X*rho/R1)*(-sin(N*phi))*(Q*pi/L1)*(Q*pi/L1)*(-cos(Q*pi*z/L1));
}

/* modeTE_Ez = 0   by definition  */

double partial2_z_modeTE_Erho (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double W = omega_TE(n,p,q,R1,L1);
   double Ks_sqr = ks_sqr(n,p,q,R1,L1);
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;


   if (rho == 0)
      return 0;
   else
      return (W/(Ks_sqr - (Q*pi/L1)*(Q*pi/L1)))*(N/rho)*jn(n,X*rho/R1)*(-sin(N*phi))*(Q*pi/L1)*(Q*pi/L1)*(-sin(Q*pi*z/L1));
}

double partial2_z_modeTE_Ephi (int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double W = omega_TE(n,p,q,R1,L1);
   double Ks_sqr = ks_sqr(n,p,q,R1,L1);
   double X = xs_np(n,p);
   double N = (double) n;
   double Q = (double) q;


   return -(W/(Ks_sqr -(Q*pi/L1)*(Q*pi/L1)))*(X/R1)*jns(n,X*rho/R1)*cos(N*phi)*(Q*pi/L1)*(Q*pi/L1)*(-sin(Q*pi*z/L1));
}


int N = 100;

double B_TM(int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double B_TM_z, B_TM_rho, B_TM_phi;
   
   B_TM_z   = 0.;
   B_TM_rho = modeTM_Hrho(n,p,q,R1,L1,z,rho,phi);
   B_TM_phi = modeTM_Hphi(n,p,q,R1,L1,z,rho,phi);

   double B = sqrt(B_TM_z*B_TM_z + B_TM_rho*B_TM_rho + B_TM_phi*B_TM_phi);

   return B;
}


double Find_Max_B_TM(int n, int p, int q, double R1, double L1)
{
   double max = 0;
   
   double hz   = L1/N;
   double hrho = R1/N;
   double hphi = 2*pi/N;

   
   for (int irho = 0; irho <= N; irho++)
   {
      for (int iphi = 0; iphi <= N; iphi++)
      {
         double BB = B_TM(n,p,q,R1,L1,0,irho*hrho,iphi*hphi);
  
         if (BB > max)
            max = BB;   
      }  
   }


   for (int irho = 0; irho <= N; irho ++) 
   {
      for (int iphi = 0; iphi <= N; iphi++)
      {
         double BB = B_TM(n,p,q,R1,L1,L1,irho*hrho,iphi*hphi);
  
         if (BB > max)
            max = BB;   
      }  
   }


   for (int iz = 0; iz <= N; iz ++) 
   {
      for (int iphi = 0; iphi <= N; iphi++)
      {
         double BB = B_TM(n,p,q,R1,L1,iz*hz,R1,iphi*hphi);
  
         if (BB > max)
            max = BB;   
      }  
   }


   return max;
}


double B_TE(int n, int p, int q, double R1, double L1, double z, double rho, double phi)
{
   double  B_TE_z, B_TE_rho, B_TE_phi;
   
   B_TE_z   = modeTE_Hz(n,p,q,R1,L1,z,rho,phi);
   B_TE_rho = modeTE_Hrho(n,p,q,R1,L1,z,rho,phi);
   B_TE_phi = modeTE_Hphi(n,p,q,R1,L1,z,rho,phi);

   double B = sqrt(B_TE_z*B_TE_z + B_TE_rho*B_TE_rho + B_TE_phi*B_TE_phi);

   return B;
}


double Find_Max_B_TE(int n, int p, int q, double R1, double L1)
{
   double max = 0;
   
   double hz   = L1/N;
   double hrho = R1/N;
   double hphi = 2*pi/N;

   
   for (int irho = 0; irho <= N; irho++)
   {
      for (int iphi = 0; iphi <= N; iphi++)
      {
         double BB = B_TE(n,p,q,R1,L1,0,irho*hrho,iphi*hphi);
  
         if (BB > max)
            max = BB;   
      }  
   }


   for (int irho = 0; irho <= N; irho++) 
   {
      for (int iphi = 0; iphi <= N; iphi++)
      {
         double BB = B_TE(n,p,q,R1,L1,L1,irho*hrho,iphi*hphi);
  
         if (BB > max)
            max = BB;   
      }  
   }


   for (int iz = 0; iz <= N; iz ++) 
   {
      for (int iphi = 0; iphi <= N; iphi++)
      {
         double BB = B_TE(n,p,q,R1,L1,iz*hz,R1,iphi*hphi);
  
         if (BB > max)
            max = BB;   
      }  
   }


   return max;
}



double F_TM_TE (int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double z, double rho, double phi)
{

   double Ez,Erho,Ephi;
   double Hz,Hrho,Hphi;
   double S;
  
   Ez   = 0.5*modeTM_Ez(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   Erho = 0.5*modeTM_Erho(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   Ephi = 0.5*modeTM_Ephi(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));

   Hz   = 0.5*modeTE_Hz(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   Hrho = 0.5*modeTE_Hrho(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   Hphi = 0.5*modeTE_Hphi(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));


   S = Ez*Hz + Erho*Hrho + Ephi*Hphi;

   return S;
}


double partial_z_F_TM_TE (int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double z, double rho, double phi)
{

   double Ez,Erho,Ephi;
   double Hz,Hrho,Hphi;
  
   double partial_z_Ez, partial_z_Erho, partial_z_Ephi;
   double partial_z_Hz, partial_z_Hrho, partial_z_Hphi; 

   double S;
 
  
   Ez   = 0.5*modeTM_Ez(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   Erho = 0.5*modeTM_Erho(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   Ephi = 0.5*modeTM_Ephi(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));

   Hz   = 0.5*modeTE_Hz(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   Hrho = 0.5*modeTE_Hrho(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   Hphi = 0.5*modeTE_Hphi(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));

   partial_z_Ez   = 0.5*partial_z_modeTM_Ez(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   partial_z_Erho = 0.5*partial_z_modeTM_Erho(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   partial_z_Ephi = 0.5*partial_z_modeTM_Ephi(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));

   partial_z_Hz   = 0.5*partial_z_modeTE_Hz(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   partial_z_Hrho = 0.5*partial_z_modeTE_Hrho(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   partial_z_Hphi = 0.5*partial_z_modeTE_Hphi(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));

   S = Ez*partial_z_Hz + partial_z_Ez*Hz + Erho*partial_z_Hrho + partial_z_Erho*Hrho + Ephi*partial_z_Hphi + partial_z_Ephi*Hphi;

   return S;
}


double partial2_z_F_TM_TE (int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double z, double rho, double phi)
{

   double Ez,Erho,Ephi;
   double Hz,Hrho,Hphi;
  
   double partial_z_Ez, partial_z_Erho, partial_z_Ephi;
   double partial_z_Hz, partial_z_Hrho, partial_z_Hphi; 

   double partial2_z_Ez, partial2_z_Erho, partial2_z_Ephi;
   double partial2_z_Hz, partial2_z_Hrho, partial2_z_Hphi; 

   double S;
 
  
   Ez   = 0.5*modeTM_Ez(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   Erho = 0.5*modeTM_Erho(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   Ephi = 0.5*modeTM_Ephi(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));

   Hz   = 0.5*modeTE_Hz(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   Hrho = 0.5*modeTE_Hrho(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   Hphi = 0.5*modeTE_Hphi(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));

   partial_z_Ez   = 0.5*partial_z_modeTM_Ez(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   partial_z_Erho = 0.5*partial_z_modeTM_Erho(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   partial_z_Ephi = 0.5*partial_z_modeTM_Ephi(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));

   partial_z_Hz   = 0.5*partial_z_modeTE_Hz(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   partial_z_Hrho = 0.5*partial_z_modeTE_Hrho(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   partial_z_Hphi = 0.5*partial_z_modeTE_Hphi(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));

   partial2_z_Ez   = 0.5*partial2_z_modeTM_Ez(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   partial2_z_Erho = 0.5*partial2_z_modeTM_Erho(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   partial2_z_Ephi = 0.5*partial2_z_modeTM_Ephi(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));

   partial2_z_Hz   = 0.5*partial2_z_modeTE_Hz(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   partial2_z_Hrho = 0.5*partial2_z_modeTE_Hrho(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   partial2_z_Hphi = 0.5*partial2_z_modeTE_Hphi(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));


   S = Ez*partial2_z_Hz + 2*partial_z_Ez*partial_z_Hz + partial2_z_Ez*Hz + Erho*partial2_z_Hrho + 2*partial_z_Erho*partial_z_Hrho + partial2_z_Erho*Hrho + Ephi*partial2_z_Hphi + 2*partial_z_Ephi*partial_z_Hphi + partial2_z_Ephi*Hphi;

   return S;
}



double F_TE_TM (int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double z, double rho, double phi)
{

   double Ez,Erho,Ephi;
   double Hz,Hrho,Hphi;
   double S;

   Ez   = 0.;
   Erho = 0.5*modeTE_Erho(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   Ephi = 0.5*modeTE_Ephi(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));


   Hz   = 0.;
   Hrho = 0.5*modeTM_Hrho(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   Hphi = 0.5*modeTM_Hphi(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));

   S = Ez*Hz + Erho*Hrho + Ephi*Hphi;

   return S;
}


double partial_z_F_TE_TM (int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double z, double rho, double phi)
{

   double Ez,Erho,Ephi;
   double Hz,Hrho,Hphi;
 
   double partial_z_Ez, partial_z_Erho, partial_z_Ephi;
   double partial_z_Hz, partial_z_Hrho, partial_z_Hphi;
  
   double S;

   Ez   = 0.;
   Erho = 0.5*modeTE_Erho(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   Ephi = 0.5*modeTE_Ephi(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));


   Hz   = 0.;
   Hrho = 0.5*modeTM_Hrho(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   Hphi = 0.5*modeTM_Hphi(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));

   partial_z_Ez   = 0.;
   partial_z_Erho = 0.5*partial_z_modeTE_Erho(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   partial_z_Ephi = 0.5*partial_z_modeTE_Ephi(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));


   partial_z_Hz   = 0.;
   partial_z_Hrho = 0.5*partial_z_modeTM_Hrho(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   partial_z_Hphi = 0.5*partial_z_modeTM_Hphi(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));


   S = Ez*partial_z_Hz + partial_z_Ez*Hz + Erho*partial_z_Hrho + partial_z_Erho*Hrho + Ephi*partial_z_Hphi + partial_z_Ephi*Hphi;

   return S;
}


double partial2_z_F_TE_TM (int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double z, double rho, double phi)
{

   double Ez,Erho,Ephi;
   double Hz,Hrho,Hphi;
 
   double partial_z_Ez, partial_z_Erho, partial_z_Ephi;
   double partial_z_Hz, partial_z_Hrho, partial_z_Hphi;

   double partial2_z_Ez, partial2_z_Erho, partial2_z_Ephi;
   double partial2_z_Hz, partial2_z_Hrho, partial2_z_Hphi;
   
   double S;

   Ez   = 0.;
   Erho = 0.5*modeTE_Erho(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   Ephi = 0.5*modeTE_Ephi(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));


   Hz   = 0.;
   Hrho = 0.5*modeTM_Hrho(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   Hphi = 0.5*modeTM_Hphi(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));


   partial_z_Ez   = 0.;
   partial_z_Erho = 0.5*partial_z_modeTE_Erho(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   partial_z_Ephi = 0.5*partial_z_modeTE_Ephi(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
 
   partial_z_Hz   = 0.;
   partial_z_Hrho = 0.5*partial_z_modeTM_Hrho(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   partial_z_Hphi = 0.5*partial_z_modeTM_Hphi(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));


   partial2_z_Ez   = 0.;
   partial2_z_Erho = 0.5*partial2_z_modeTE_Erho(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));
   partial2_z_Ephi = 0.5*partial2_z_modeTE_Ephi(n2,p2,q2,R1,L1,z,rho,phi)/(Max_TE*sqrt(2));

   partial2_z_Hz   = 0.;
   partial2_z_Hrho = 0.5*partial2_z_modeTM_Hrho(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));
   partial2_z_Hphi = 0.5*partial2_z_modeTM_Hphi(n1,p1,q1,R1,L1,z,rho,phi)/(Max_TM*sqrt(2));


   S = Ez*partial2_z_Hz + 2*partial_z_Ez*partial_z_Hz + partial2_z_Ez*Hz + Erho*partial2_z_Hrho + 2*partial_z_Erho*partial_z_Hrho + partial2_z_Erho*Hrho + Ephi*partial2_z_Hphi + 2*partial_z_Ephi*partial_z_Hphi + partial2_z_Ephi*Hphi;

  
   return S;
}


double F_TMTE_p (int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double z, double rho, double phi)
{
    return (F_TE_TM(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,z,rho,phi) - F_TM_TE(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,z,rho,phi));
}

double partial_z_F_TMTE_p (int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double z, double rho, double phi)
{
    return (partial_z_F_TE_TM(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,z,rho,phi) - partial_z_F_TM_TE(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,z,rho,phi));
}

double partial2_z_F_TMTE_p (int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double z, double rho, double phi)
{
    return (partial2_z_F_TE_TM(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,z,rho,phi) - partial2_z_F_TM_TE(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,z,rho,phi));
}


double F_TMTE_m (int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double z, double rho, double phi)
{
    return (F_TE_TM(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,z,rho,phi) + F_TM_TE(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,z,rho,phi));
}

double partial_z_F_TMTE_m (int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double z, double rho, double phi)
{
    return (partial_z_F_TE_TM(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,z,rho,phi) + partial_z_F_TM_TE(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,z,rho,phi));
}

double partial2_z_F_TMTE_m (int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double Max_TM, double Max_TE, double z, double rho, double phi)
{
    return (partial2_z_F_TE_TM(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,z,rho,phi) + partial2_z_F_TM_TE(n1,p1,q1,n2,p2,q2,R1,L1,Max_TM,Max_TE,z,rho,phi));
}


double rs(double z, double zs, double rho, double rhos, double phi, double phis)
{
   return (z - zs)*(z - zs) + rho*rho + rhos*rhos - 2*rho*rhos*cos(phi - phis);
}


