% A script to add paths.
clear;

currentpath = pwd;

addpath (genpath ('geopdes-periodic'))
addpath (genpath ('geopdes'))
addpath (genpath ('nurbs'))
addpath (genpath ('geopdes_hierarchical'))
addpath (genpath ('data'))
addpath (genpath ('quasi_interpolants'))
addpath (genpath ('solvers'))
addpath (genpath ('New_qi'))

ex_laplace_iso_ring
close all
clear
clc