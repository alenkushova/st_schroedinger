% A script to add paths.
clear;

currentpath = pwd;

addpath (genpath ('geopdes'))
addpath (genpath ('nurbs'))
addpath (genpath ('geopdes_hierarchical'))
addpath (genpath ('data'))
addpath (genpath ('solvers'))
addpath (genpath ('examples'))
addpath (genpath ('inf_sup_test'))

close 
clear
clc