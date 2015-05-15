function Cc = Cc_ref()
%CC_REF - Calculate Cc integral at fixed parameters, reference for comparison with C implementation
%
%    Cc = CC_REF()

% parameters
N = 32;
R = 7.5;
L = 12;
J = 48;

% calculate quadrature weights
quadwI1 = fourierI1(N,L,J,R);

% load matrix entries from disk
W = cell(1,4);
for j=1:4
	W{j} = loadData(sprintf('data/W%d.dat',j-1),'real',[N,N]);
end

lambdaT = eig(pauliToMatrix([W{1}(5),W{2}(5),W{3}(5),W{4}(5)]));
fprintf('example eigenvalues (should be in [0,1]):\n');
disp(lambdaT);

% represent in Fourier space
for j=1:4
	% divide by N^2 due to normalization convention
	W{1,j} = fft2(W{1,j})/N^2;
end

% calculate Cc integral
fprintf('Calculating Cc integral...\n');
tic;
Cc = CcInt(W,quadwI1);
toc
% store as complex numbers
Cc{1} = complex(Cc{1});

% save to disk
fprintf('Saving results to disk...\n');
for j=1:4
	writeData(sprintf('data/Cc%d_ref.dat',j-1),Cc{j});
end
