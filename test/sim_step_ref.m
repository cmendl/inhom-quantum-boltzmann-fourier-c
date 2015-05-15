function Wnext = sim_step_ref()
%SIM_STEP_REF - Perform a simulation time step, reference for comparison with C implementation
%
%    Wnext = SIM_STEP_REF()

% parameters
N = 32;
R = 7.5;
L = 12;
J = 16;
M = 32;

% spatial mesh
h = 0.08;
x = h * (0:7);

% time step
dt = 0.002;

% Courant number must be <= 1 according to CFL condition
courant = L*dt/h;
fprintf('Courant number: %g (should be <= 1)\n',courant);

fprintf('N: %d\n',N);
fprintf('L: %g\n',L);
fprintf('J: %d\n',J);
fprintf('M: %d\n',M);
fprintf('R: %g\n',R);
fprintf('h: %g\n',h);
fprintf('dt: %g\n',dt);
fprintf('number of finite volumes: %d\n',length(x));

% velocity grid
vgrid = [0:N/2-1,-N/2:-1]*(2*L)/N;
[vx,vy] = meshgrid(vgrid,vgrid);

% calculate quadrature weights
quadwI1 = fourierI1(N,L,J,R);
quadwI2 = fourierI2(N,L,J,R);
quadwI3 = fourierI3(N,L,J,M,R);
quadwI4 = fourierI4(N,L,J,R);

% load matrix entries from disk
W = loadData('data/Wx.dat','real',[N,N,4,length(x)]);

% one simulation time step

Wnext = zeros(N,N,4,length(x));

% temporary variable
Wtmp = cell(length(x),4);
for ix=1:length(x)
	for j=1:4
		Wtmp{ix,j} = zeros(N);
	end
end

fprintf('Simulation time step...\n');

tic;

% half step of transport term (Trotter splitting)
for k1=1:size(vx,1)	% for each velocity component...
	for k2=1:size(vx,2)
		for j=1:4
			% copy x values into array
			Wx = zeros(size(x));
			for ix=1:length(x)
				Wx(ix) = W(k1,k2,j,ix);
			end
			% transport with velocity vx(k) for half a time step
			Wx = slopeLimStep(h,0.5*dt,vx(k1,k2),Wx);
			% copy values to temporary variable
			for ix=1:length(x)
				Wtmp{ix,j}(k1,k2) = Wx(ix);
			end
		end
	end
end

% represent in Fourier space
for ix=1:length(x)
	for j=1:4
		% divide by N^2 due to normalization convention
		Wtmp{ix,j} = fft2(Wtmp{ix,j})/N^2;
	end
end

% local collisions at each x
for ix=1:length(x)
	% midpoint rule
	Wloc = Wtmp(ix,:);
	Cc = CcInt(Wloc,quadwI1);
	Cd = CdInt(Wloc,quadwI2,quadwI3,quadwI4);
	% CB = magneticRot(Wloc,Bext);
	% temporary variable for midpoint rule
	Wast = cell(1,4);
	for j=1:4
		Wast{j} = Wloc{j} + dt/2 * (Cc{j} + Cd{j});	% + CB{j}
	end

	Cc = CcInt(Wast,quadwI1);
	Cd = CdInt(Wast,quadwI2,quadwI3,quadwI4);
	% CB = magneticRot(Wast,Bext);
	for j=1:4
		Wtmp{ix,j} = Wtmp{ix,j} + dt * (Cc{j} + Cd{j});	% + CB{j}
	end
end

% transform to physical velocity space
for ix=1:length(x)
	for j=1:4
		% multiply by N^2 due to normalization convention
		Wtmp{ix,j} = N^2*real(ifft2(Wtmp{ix,j}));
	end
end

% half step of transport term (Trotter splitting)
for k1=1:size(vx,1)	% for each velocity component...
	for k2=1:size(vx,2)
		for j=1:4
			% copy x values into array
			Wx = zeros(size(x));
			for ix=1:length(x)
				Wx(ix) = Wtmp{ix,j}(k1,k2);
			end
			% transport with velocity vx(k) for half a time step
			Wx = slopeLimStep(h,0.5*dt,vx(k1,k2),Wx);
			% copy values back
			for ix=1:length(x)
				Wnext(k1,k2,j,ix) = Wx(ix);
			end
		end
	end
end

toc

% save to disk
fprintf('Saving results to disk...\n');
writeData('data/Wx_next_ref.dat',Wnext);
