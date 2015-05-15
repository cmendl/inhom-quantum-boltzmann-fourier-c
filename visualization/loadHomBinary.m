function loadHomBinary()
%LOADHOMBINARY - Load and visualize simulation data of the spatially homogeneous quantum Boltzmann equation
%
%    LOADHOMBINARY()

N = 32;
R = 7.5;
% periodic volume [-L,L]^2
L = ceil(0.5*(3+sqrt(2))/sqrt(2) * R);
J = 32;
M = N;

% time steps
dt = 0.001;
tmax = 0.1;
t = 0:dt:tmax;

fprintf('N: %d\n',N);
fprintf('L: %g\n',L);
fprintf('J: %d\n',J);
fprintf('M: %d\n',M);
fprintf('R: %g\n',R);
fprintf('dt: %g\n',dt);

% load simulation data from disk
% Wigner function represented in physical velocity space
datafile = 'Wevolv_hom.dat';
W = loadData(['../data/',datafile],'real',[N,N,4,length(t)]);

visualizeHomBoltzmann(N,L,J,M,R,dt,W);
