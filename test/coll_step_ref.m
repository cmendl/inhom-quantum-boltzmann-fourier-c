function Wnext = coll_step_ref()
%COLL_STEP_REF - Perform a collision time step at fixed parameters, reference for comparison with C implementation
%
%    Wnext = COLL_STEP_REF()

% parameters
N = 32;
R = 7.5;
L = 12;
J = 32;
M = 32;

% time step
dt = 0.001;

% calculate quadrature weights
quadwI1 = fourierI1(N,L,J,R);
quadwI2 = fourierI2(N,L,J,R);
quadwI3 = fourierI3(N,L,J,M,R);
quadwI4 = fourierI4(N,L,J,R);

% load matrix entries from disk
W = cell(1,4);
for j=1:4
	W{j} = loadData(sprintf('data/W%d.dat',j-1),'real',[N,N]);
end

% represent in Fourier space
for j=1:4
	% divide by N^2 due to normalization convention
	W{1,j} = fft2(W{1,j})/N^2;
end

% one collision time step

Wnext = cell(1,4);
Wast  = cell(1,4);

fprintf('Collision time step...\n');

tic;

% midpoint rule

Cc = CcInt(W,quadwI1);
Cd = CdInt(W,quadwI2,quadwI3,quadwI4);
for j=1:4
	Wast{j} = W{1,j} + dt/2 * (Cc{j} + Cd{j});
end

Cc = CcInt(Wast,quadwI1);
Cd = CdInt(Wast,quadwI2,quadwI3,quadwI4);
for j=1:4
	Wnext{1,j} = W{1,j} + dt * (Cc{j} + Cd{j});
end

toc

% save to disk
fprintf('Saving results to disk...\n');
for j=1:4
	writeData(sprintf('data/WFnext%d_ref.dat',j-1),Wnext{j});
end
