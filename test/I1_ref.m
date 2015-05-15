function I1 = I1_ref()
%I1_REF - Evaluate the I1 integral at fixed parameters, reference for comparison with C implementation
%
%    I1 = I1_REF()

N = 32;
L = 12;
J = 32;
R = 7.5;

% load matrix entries from disk
f = loadData('data/f.dat','complex',[N,N]);
g = loadData('data/g.dat','complex',[N,N]);
h = loadData('data/h.dat','complex',[N,N]);

quadwI1 = fourierI1(N,L,J,R);

I1 = I4integral(f,g,h,quadwI1.psiR1,quadwI1.psiR2,quadwI1.weight);

% save to disk
writeData('data/I1_ref.dat',I1);
