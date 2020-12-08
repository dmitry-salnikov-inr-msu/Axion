#include <iostream>
#include <fstream>

using namespace std;


void Radiation_Pattern_Theta_TMTE_p(int n1, int p1, int q1, int n2, int p2, int q2, double R1, double L1, double m, double phi, double eps)
{

   double N1 = (double) n1;
   double P1 = (double) p1;
   double Q1 = (double) q1;
   double N2 = (double) n2;
   double P2 = (double) p2;
   double Q2 = (double) q2;
 
   double Max_TM = Find_Max_B_TM(n1,p1,q1,R1,L1);
   double Max_TE = Find_Max_B_TE(n2,p2,q2,R1,L1);
 
   ofstream out("Radiation_Pattern (theta) TMTE omega_p.txt",std::ios::app);
   
   int start = MPI_Wtime(); 
    
   double R = 10.;

  
   for (double theta = 0.; theta <= (pi + 0.001); theta += 0.01*pi)
   {
      double aCp = TMTE_cos_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,R*cos(theta) + 0.5*L1,R*sin(theta),phi,eps);
      double aSp = TMTE_sin_p(N1,P1,Q1,N2,P2,Q2,R1,L1,Max_TM,Max_TE,m,R*cos(theta) + 0.5*L1,R*sin(theta),phi,eps);

      out << theta << " " << 0.5*R*R*(aCp*aCp + aSp*aSp) << endl;
   }
 
  
   int finish = MPI_Wtime();

   cout << "Time = " << (finish - start)/60 << " min " << (finish - start)%60 << " sec." << endl;  

   out.close();

}

  
