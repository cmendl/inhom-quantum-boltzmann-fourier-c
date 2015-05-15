function I2 = I2_ref2()
%I2_REF2 - Evaluate the I2 integral at fixed parameters, reference for comparison with C implementation
%
%    I2 = I2_REF2()

N = 32;
L = 12;
J = 32;
R = 7.5;

% load matrix entries from disk
f = loadData('data/f.dat','complex',[N,N]);
g = loadData('data/g.dat','complex',[N,N]);
h = loadData('data/h.dat','complex',[N,N]);

quadwI2 = fourierI2(N,L,J,R);

I2 = I2integral(f,g,h,quadwI2.psiR1,quadwI2.psiR2,quadwI2.weight);

% save to disk
writeData('data/I2_ref2.dat',I2);
