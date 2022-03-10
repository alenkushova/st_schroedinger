% A script to add paths.
clear;

currentpath = pwd;

addpath (genpath ('geopdes'))
addpath (genpath ('nurbs'))
addpath (genpath ('geopdes_hierarchical'))
addpath (genpath ('data'))
addpath (genpath ('solvers'))

close 
clear
clc