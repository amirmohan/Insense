clc
clear all
close all

%% dimension of the data
D = 100;
N = 100;

%% generate a random dataset
Phi = unifrnd(0,1,[N,D]);

%% number of sensors 
M = 10;

%% set parameters    
param.tol = 1e-5;
param.maxit = 20;
param.backtrack = 1;
param.gamma = 1;
param.init = sqrt(1/D)*ones(D,1);

%% run Insense
[z] = Insense(Phi,M,param);
[zval,zind] = sort(z,'descend'); 
zind(1:M)'












