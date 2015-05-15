function loadInhomBinary(paramfile)
%LOADINHOMBINARY - Load and visualize simulation data of the spatially inhomogeneous quantum Boltzmann equation
%
%    LOADINHOMBINARY(paramfile)

params = parseParameterFile(paramfile);

% N is fixed in C simulation
N = 32;

fprintf('N: %d\n',N);
fprintf('L: %g\n',params.L);
fprintf('J: %d\n',params.J);
fprintf('M: %d\n',params.M);
fprintf('R: %g\n',params.R);
fprintf('h: %g\n',params.h);	% spatial mesh width
fprintf('dt: %g\n',params.dt);	% time step
fprintf('number of time steps: %d\n',params.numsteps);
fprintf('number of finite volumes: %d\n',params.numVol);

% load simulation data from disk
% Wigner function represented in physical velocity space
W = loadData(params.filenameWevolv,'real',[N,N,4,params.numVol,params.numsteps]);

% extract file name without path
[~,name,ext] = fileparts(params.filenameWevolv);
datafile = [name,ext];

visualizeInhomBoltzmann(N,params.L,params.J,params.M,params.R,params.h,params.dt,params.boundaryType,...
	W,['C simulation data in "',strrep(datafile,'_','-'),'"']);
