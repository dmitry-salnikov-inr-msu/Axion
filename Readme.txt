This is a program for calculating spacial distribution of axion field, generated by two electromagnetic modes located in a cylindrical cavity, 
written by Dmitry Salnikov (salnikov.dv16@physics.msu.ru).

In addition to spacial distribution of axion field, the program shows energy distribution of axion field, 
and a predicted constarint on the axion-photon coupling in an experimental setup with two coaxial cavities. 
The physical results are presented in [1].

The numerical integrations in the program are based on package for adaptive multidimensional integration (cubature) 
of vector-valued integrands over hypercubes, written by Steven G. Johnson (https://github.com/stevengj/cubature). MPI is used for parallelization. 

The program includes five files: "TMnpq,TEnpq.h", "coupling_constant.h", "energy_density.h", "radiation_pattern.h" and "main.cpp".
The file "main.cpp" consists of the following functions:

1) Distribution_of_Axion_Field_RhoZ_TMTE_cos_p(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z_min,z_max,rho_min,rho_max,phi,eps);
2) Distribution_of_Axion_Field_RhoZ_TMTE_sin_p(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z_min,z_max,rho_min,rho_max,phi,eps);

These functions calculate formulae (11) from [1] for TMn1p1q1 and TEn2p2q2 modes for axion wave with frequency  omega_p [m^-1] 
for cavity dimensions R = R1[m], L = L1[m] for fixed axion mass m_a [m^-1] (m_a[eV] = 2*10^-7*m_a[m^-1]) and fixed angle phi
for different rho and z. The area of calculation  [(rho,z)| z_min [m] < z < z_max [m], rho_min [m] < rho < rho_max[m]] should be defined. 
Cylinder's borders are z = 0, z = L1, rho = R1. Step by coordinate rho equals 0.05 [m], step by coordinate z is defined by the relation 
(z_max - z_min)/(size - 1), where size is the number of processes. eps is relative precision of calculations. 
The data in .txt file is presented as

rho  z  a_cos(sin)_p(rho,z)
... ...        ... 

3) Distribution_of_Axion_Field_RhoZ_TMTE_cos_m(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z_min,z_max,rho_min,rho_max,phi,eps);
4) Distribution_of_Axion_Field_RhoZ_TMTE_sin_m(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z_min,z_max,rho_min,rho_max,phi,eps);

These functions are similar to the previous ones but they calculate spatial distribution of axion field for axion wave 
with frequency equals omega_m [m^-1].

The following functions calculate average axion energy density for axion wave with frequency equals omega_p [m^-1]
as a function of axion mass and z for fixed rho and phi:
5) Distribution_of_Axion_Energy_Density_MZ_TMTE_p(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m_min,m_max,z_min,z_max,rho,phi,eps);

as a function of rho and z for fixed axion mass and phi:
6) Distribution_of_Axion_Energy_Density_RhoZ_TMTE_p(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z_min,z_max,rho_min,rho_max,phi,eps);

as a function of rho and phi (0 < phi < 2*pi) for fixed axion mass and z:
7) Distribution_of_Axion_Energy_Density_RhoPhi_TMTE_p(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z,rho_min,rho_max,eps);

The following functions are similar to the previous ones but they calculate average axion energy density for axion wave 
with frequency equals omega_m [m^-1].

8) Distribution_of_Axion_Energy_Density_MZ_TMTE_m(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m_min,m_max,z_min,z_max,rho,phi,eps);
9) Distribution_of_Axion_Energy_Density_RhoZ_TMTE_m(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z_min,z_max,rho_min,rho_max,phi,eps);
10) Distribution_of_Axion_Energy_Density_RhoPhi_TMTE_m(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,m,z,rho_min,rho_max,eps);

In all these functions step by the first coordinate (or mass) is fixed, step by the second coordinate is defined by 
the relation (second_coordinate_max - second_coordinate_min)/(size - 1), where size is the number of processes.

The following function compute the constarint on the coupling constant g_{a\gamma\gamma} for production modes 
TMn1p1q1, TEn2p2q2 and detection mode TMn3p3q3 (TM010 mode in [1]), for production cavity dimensions R = R1, L = L1 and detection 
cavity dimensions R = R2, L = L2 in axion mass range m_min < m_a < m_max. Step by axion mass is dm.
dist is z coordinate of border of detection cavity (Delta = dist - L1). Parallelization is used to
calculate integrals:

11)Coupling_Constant_TMTE_p(argn,argv,n1,p1,q1,n2,p2,q2,R1,L1,dist,n3,p3,q3,R2,L2,m_min,m_max,dm,eps);

The following function calculate radiation pattern (0 < theta < pi) for fixed axion mass and angle phi.

12)Radiation_Pattern_Theta_TMTE_p(n1,p1,q1,n2,p2,q2,R1,L1,m,phi,eps);

!!!You can use only one function in a single compilation.

!Compilation and launch!

module load openmpi-x86_64
mpic++ -o Axion main.cpp hcubature.c -lm
mpirun -np number_of_processes ./Axion &

If you want convert data for 2d contour plots in origin matrix format you can use converter.cpp.
You should name input file "input.txt". You will receive output data in file "output.txt".
nx is number of columns, ny is number of lines in output file.  
Note that the converter also outputs data displayed relative to the x-axis along with the source data.

Also you can make animation of propagation of axion wave using program in directory "Axion animation (python)".
You should calculate formulae (11) in ranges [0 < rho < 5, 0.5*R1 < z < 0.5*R1 + 2.5] and use converter 
(number of columns is nx = 101, number of lines is ny = 51).

-----------------------------
[1] Dmitry Salnikov, Petr Satunin, D. V. Kirpichnikov, Maxim Fitkevich. 
"Examining axion-like particles with superconducting radio-frequency cavity."
arXiv: https://arxiv.org/abs/2011.12871 
doi: https://link.springer.com/article/10.1007/JHEP03(2021)143
