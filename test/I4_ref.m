function I4 = I4_ref()
%I4_REF - Evaluate the I4 integral at fixed parameters, reference for comparison with C implementation
%
%    I4 = I4_REF()

N = 32;
L = 12;
J = 32;
R = 7.5;

% load matrix entries from disk
f = loadData('data/f.dat','complex',[N,N]);
g = loadData('data/g.dat','complex',[N,N]);
h = loadData('data/h.dat','complex',[N,N]);

quadwI4 = fourierI4(N,L,J,R);

I4 = I4integral(f,g,h,quadwI4.psiR1,quadwI4.psiR2,quadwI4.weight);

% save to disk
writeData('data/I4_ref.dat',I4);
