function Un1M = slope_lim_ref2()
%SLOPE_LIM_REF2 - Slope limiter step with Maxwell boundary condition, reference for comparison with C implementation
%
%    Un1M = SLOPE_LIM_REF2()

% spatial mesh width
h = 0.01;
% time step
dt = 0.005;

% number of finite volumes
N = 101;

% A matrix
A = diag([1,0.5,-0.5,-1]);

% accommodation coefficient for Maxwell reflection operator
alpha = 0.4;

% pre-determined fluxes
fluxL = 0.7;
fluxR = 0.45;

% incoming state from the left
UmaxwL = [0.2;-pi/4;0;0];		% non-zero for velocities > 0
% incoming state from the right
UmaxwR = [0;0;exp(-1);-0.8];	% non-zero for velocities < 0

% Courant number, Eq. (10.56), must be <= 1
courant = max(abs(eig(A)))*dt/h;
fprintf('Courant number: %g (should be <= 1)\n',courant);

% load finite volume state from disk
Un = loadData('data/UnM.dat','real',[4,N]);

% slope limiter step with Maxwell boundary conditions
Un1M = slopeLimStepMaxwell(h,dt,A,alpha,fluxL,fluxR,UmaxwL,UmaxwR,Un);

% save to disk
fprintf('Saving results to disk...\n');
writeData('data/Un1_Maxwell_ref.dat',Un1M);
