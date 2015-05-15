function [Un1P,Un1D] = slope_lim_ref()
%SLOPE_LIM_REF - Slope limiter step in 1D with periodic and Dirichlet boundary conditions, reference for comparison with C implementation
%
%    [Un1P,Un1D] = SLOPE_LIM_REF()

% spatial mesh width
h = 0.08;
% time step
k = 0.02;

% number of finite volumes
N = 128;

% A matrix, here 1 x 1
A = -1.2;

% Courant number, Eq. (10.56), must be <= 1
courant = abs(A)*k/h;
fprintf('Courant number: %g (should be <= 1)\n',courant);

% load finite volume state from disk
Un = loadData('data/Un.dat','real',[1,N]);

% slope limiter step with periodic boundary conditions
Un1P = slopeLimStep(h,k,A,Un);
Un1D = slopeLimStepDirichlet(h,k,A,Un);

% save to disk
fprintf('Saving results to disk...\n');
writeData('data/Un1_periodic_ref.dat', Un1P);
writeData('data/Un1_Dirichlet_ref.dat',Un1D);
