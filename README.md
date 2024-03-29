# Axion
This is a program for calculating spacial distribution of axion field, generated by two electromagnetic modes located in a cylindrical cavity, written by Dmitry Salnikov (salnikov.dv16@physics.msu.ru).

In addition to spacial distribution of axion field, the program shows energy distribution of axion field, and a predicted constarint on the axion-photon coupling in an experimental setup with two coaxial cavities. The physical results are presented in [1].

The numerical integrations in the program are based on package for adaptive multidimensional integration (cubature) of vector-valued integrands over hypercubes, written by Steven G. Johnson (https://github.com/stevengj/cubature). MPI is used for parallelization.

Description of program see in Readme.txt.

-------------------------------
[1] Dmitry Salnikov, Petr Satunin, D. V. Kirpichnikov, Maxim Fitkevich. "Examining axion-like particles with superconducting radio-frequency cavity."

arXiv: https://arxiv.org/abs/2011.12871

doi: https://link.springer.com/article/10.1007/JHEP03(2021)143
