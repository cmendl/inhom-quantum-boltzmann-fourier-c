function I2 = I2_ref()
%I2_REF - Evaluate the I2 integral at fixed parameters, reference for comparison with C implementation
%
%    I2 = I2_REF()

N = 32;
J = 7;

% load matrix entries from disk
f = loadData('data/f.dat','complex',[N,N]);
g = loadData('data/g.dat','complex',[N,N]);
h = loadData('data/h.dat','complex',[N,N]);

weight = loadData('data/weight.dat','real',[1,J]);
psiR1 = cell(1,J);
psiR2 = cell(1,J);
for j=1:J
	psiR1{j} = loadData(['data/psiR1_',num2str(j),'.dat'],'real',[N,N]);
	psiR2{j} = loadData(['data/psiR2_',num2str(j),'.dat'],'real',[N,N]);
end

I2 = I2integral(f,g,h,psiR1,psiR2,weight);

% save to disk
writeData('data/I2_ref.dat',I2);
