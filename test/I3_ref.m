function I3 = I3_ref()
%I3_REF - Evaluate the I3 integral at fixed parameters, reference for comparison with C implementation
%
%    I3 = I3_REF()

N = 32;
L = 12;
J = 32;
R = 7.5;
M = 7;

% load matrix entries from disk
f = loadData('data/f.dat','complex',[N,N]);
g = loadData('data/g.dat','complex',[N,N]);
h = loadData('data/h.dat','complex',[N,N]);

quadwI3 = fourierI3(N,L,J,M,R);

I3 = I3integral(f,g,h,quadwI3.psiR1,quadwI3.psiR2,quadwI3.weight);

% save to disk
writeData('data/I3_ref.dat',I3);
