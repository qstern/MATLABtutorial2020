%% Matlab tutorial - Session 3 - Basics of NMR simulation
% This document contains the elements presented in a tutorial on MATLAB's
% basics on 23.4.20 for the CRMN members by Quentin Chappuis

clear all
%% Defining Pauli matrices
% Uses the scalar ai as the imaginary number i
ai=sqrt(-1);

% Creates the three angular momentum bases on the Pauli matrices
Ix=1/2*[0 1;1 0];
Iy=1/2*[0 -ai;ai 0];
Iz=1/2*[1 0;0 -1];

%% Single spin simulation
% Simulation of the spectrum for a single spins 1/2 under the effect of 
% field B0

% Initial density matrix
rho0=[1 0; 0 0];

% Expectation values of the angular momentum of the initial state
disp('Before pulsing')
disp(['x: ' num2str(trace(rho0*Ix))])
disp(['y: ' num2str(trace(rho0*Iy))])
disp(['z: ' num2str(trace(rho0*Iz))])

% Free evolution Hamiltonian
w0=30*2*pi;     % Larmor frequency in rad/s
H0=w0*Iz;       % Hamiltonian in natural units

% Pulse Hamiltonian
w1=30e3*2*pi;   % Nutation frequency in rad/s
H1=w1*Iy;       % Hamiltonian in natural units
taup=pi/2/w1;   % Pulse length for pi/2

% Propagator for the pulse Hamiltonian during the length of the pulse
U1=expm(-ai*H1*taup);

% Propagation of the density matrix under the effect of the pulse
% Hamiltonian using the sandwitch formula
rho1=U1'*rho0*U1;

% Expectation values of the angular momentum of the system after pulsing
disp('After pulsing')
disp(['x: ' num2str(trace(rho1*Ix))])
disp(['y: ' num2str(trace(rho1*Iy))])
disp(['z: ' num2str(round(trace(rho1*Iz),5))])

% Parameters for the Free Induction Decay (FID)
aq=5;           % Acquisition times in s
SW=100*2*pi;    % Spectral window in rad/s

% Number of points of the FID in the next power of 2
N=2^round(log2(SW*aq/2/pi));

% Time axis of the FID
t=linspace(0,aq,N);

% Empty vectors for the expectation values along time/propagation of the
% free evolution Hamiltonian. They correspond to the signal along x and y
% assuming quadrature detection.
Sxt=zeros(N,1);
Syt=zeros(N,1);

% Computation of the FID
for i=1:N
    % Propagator taking the density matrix from t=0 to t=t(i), that is,
    % from after the pulse to time t.
    U0t=expm(-ai*H0*t(i));
    
    % Propagation of the denisty matrix using the sandwitch formula
    rhot=U0t'*rho1*U0t;
    
    % Expectation value of the angular momentum along x and y.
    Sxt(i)=real(trace(rhot*Ix));
    Syt(i)=real(trace(rhot*Iy));
end

% Function performing the Fourier transform
% The first two arguments are mendatory: time and signal of the FID. The
% signal must be expressed as complex number (quadratur detection)
% In this example, I used some options: 'ph' for the phase, 'lb' for line
% broadening (uses a momo exponential window), 'zero' for zero filling
% 'spectrum' and 'fid' to plot the spectrum and FID.
Signal1D(t, Sxt+ai*Syt, 'ph', pi, 'lb', 1,'zero',16*N, 'spectrum', true, 'fid', true);

%% Single spin simulation -  in the laboratory frame
% Simulation of the spectrum for a single spins 1/2 under the effect of 
% field B0 expressed in the laboratory frame

% Initial density matrix
rho0=[1 0; 0 0];

% Expectation values of the angular momentum of the initial state
disp('Before pulsing')
disp(['x: ' num2str(trace(rho0*Ix))])
disp(['y: ' num2str(trace(rho0*Iy))])
disp(['z: ' num2str(trace(rho0*Iz))])

% Free evolution Hamiltonian
% For simplicity, the carrier frequency and the reference frequency are
% assumed to be the same.
O1=400e6*2*pi;  % Carrier frequency
d=1.5e-6;       % Chemical shift
w0=O1*(1+d);    % Larmor frequency in rad/s
H0=w0*Iz;       % Hamiltonian in natural units

% Pulse Hamiltonian
w1=30e3*2*pi;   % Nutation frequency in rad/s
H1=w1*Iy;       % Hamiltonian in natural units
taup=pi/2/w1;   % Pulse length for pi/2

% Propagator for the pulse Hamiltonian during the length of the pulse
U1=expm(-ai*H1*taup);

% Propagation of the density matrix under the effect of the pulse
% Hamiltonian using the sandwitch formula
rho1=U1'*rho0*U1;

% Expectation values of the angular momentum of the system after pulsing
disp('After pulsing')
disp(['x: ' num2str(trace(rho1*Ix))])
disp(['y: ' num2str(trace(rho1*Iy))])
disp(['z: ' num2str(round(trace(rho1*Iz),5))])

% Parameters for the Free Induction Decay (FID)
aq=5;           % Acquisition times in s
SW=2000*2*pi;    % Spectral window in rad/s

% Number of points of the FID in the next power of 2
N=2^round(log2(SW*aq/2/pi));

% Time axis of the FID
t=linspace(0,aq,N);

% Empty vectors for the expectation values along time/propagation of the
% free evolution Hamiltonian. They correspond to the signal along x and y
% assuming quadrature detection.
Sxt=zeros(N,1);
Syt=zeros(N,1);

% Computation of the FID
for i=1:N
    % Propagator taking the density matrix from t=0 to t=t(i), that is,
    % from after the pulse to time t.
    U0t=expm(-ai*H0*t(i));
    
    % Propagation of the denisty matrix using the sandwitch formula
    rhot=U0t'*rho1*U0t;
    
    % Expectation value of the angular momentum along x and y.
    Sxt(i)=real(trace(rhot*Ix));
    Syt(i)=real(trace(rhot*Iy));
end

% Fourier transform
% The spectrum can be plotted in ppm unit (option 'unit', value 'ppm')
% provided the carrier frequency is defined (option 'O1').
% The carrier is mixed with the FID in a smilar way that is performed by a
% spectrometer.
Signal1D(t, Sxt+ai*Syt,'O1', O1,'unit','ppm', 'ph', pi,'zero',16*N, 'lb', 1, 'spectrum', true);

% REMARK: the error that appeared during the tutorial session (cf youtube
% video) was a bug in Signal1D when zero filling was used and the frequency
% was expresed in ppm. It has been fixed.

%% Two spin system heteronuclear
% Simulation of the spectrum for a 1H spin coupled to a 13C spin
% Spin 1: 1H
% Spin 2: 13C

% Angular momentum operators in the Hilbert space of a single spin
Ix=1/2*[0 1;1 0];
Iy=1/2*[0 -ai;ai 0];
Iz=1/2*[1 0;0 -1];

% Angular momentum operators in the Hilbert space of a two-spin system^
% kron computes the Kronecker product.
I1z=kron(Iz, eye(2));
I2z=kron(eye(2),Iz);
I1x=kron(Ix, eye(2));
I2x=kron(eye(2),Ix);
I1y=kron(Iy, eye(2));
I2y=kron(eye(2),Iy);
Ix=I1x+I2x;
Iy=I1y+I2y;
Iz=I1z+I2z;

% Function computing the density matrix for several spins as a function of
% field B0 in T and temperature in K.
% [1 13] refers to the nuclei using their atomic mass
% The field is 9.8 T and the temperature is 298 K.
rho0=ThermalEqDensityOp([1 13], 9.8, 298);

O1=400e6*2*pi;  % Carrier frequency for the 1H channel
d=1e-6;         % 1H chemical shift
wI1=O1*(1+d);   % 1H Larmor frequency
wI2=O1/4;       % 13C Larmor frequency
J12=20;         % J-coupling

% Free precession Hamiltonian
H0=wI1*I1z + wI2*I2z + 2*pi*J12*(I1x*I2x+I1y*I2y+I1z*I2z);

% Pulse Hamiltonian
% Operator I1y is used rather than Iy because the pulse only acts on the 1H
% spin
w1=30e3*2*pi;   % Nutation frequency in rad/s
H1=w1*I1y;      % Hamiltonian in natural units
taup=pi/2/w1;   % Pulse length for pi/2

% Propagator for the pulse Hamiltonian during the length of the pulse
U1=expm(-ai*H1*taup);

% Propagation of the density matrix under the effect of the pulse
% Hamiltonian using the sandwitch formula
rho1=U1'*rho0*U1;

% Parameters for the Free Induction Decay (FID)
aq=3;           % Acquisition times in s
SW=100*2*pi;    % Spectral window in rad/s

% Number of points of the FID in the next power of 2
N=2^round(log2(SW*aq*2*pi));

% Time axis for the FID
t=linspace(0,aq,N);

% Empty vectors for the expectation values along time/propagation of the
% free evolution Hamiltonian. They correspond to the signal along x and y
% assuming quadrature detection.
Ixt=zeros(N,1);
Iyt=zeros(N,1);

% Computation of the FID
for i=1:N
    U=expm(-ai*t(i)*H0);
    rhot=U'*rho1*U;
    
    % Computation of the expectation values. Operators I1x and I1y are used
    % rather than Ix and Iy since only spin 1 (1H) is detected.
    Ixt(i)=real(trace(I1x*rhot));
    Iyt(i)=real(trace(I1y*rhot));
end

Signal1D(t, Ixt+ai*Iyt,'lb',1,'spectrum', true,'fid',true,'O1', O1, 'unit', 'ppm', 'ph', pi);

%% Two spin system homonuclear case
% Simulation of the spectrum for a pair of J-coupled 1H spins using the
% full J Hamiltonian. The parameters chosen here imply that the two 1H
% spins are stronhly coupled, resulting in the 'roof effect'

% Angular momentum operators in the Hilbert space of a single spin
Ix=1/2*[0 1;1 0];
Iy=1/2*[0 -ai;ai 0];
Iz=1/2*[1 0;0 -1];

% Angular momentum operators in the Hilbert space of a two-spin system^
% kron computes the Kronecker product.
I1z=kron(Iz, eye(2));
I2z=kron(eye(2),Iz);
I1x=kron(Ix, eye(2));
I2x=kron(eye(2),Ix);
I1y=kron(Iy, eye(2));
I2y=kron(eye(2),Iy);
Ix=I1x+I2x;
Iy=I1y+I2y;
Iz=I1z+I2z;

% Function computing the density matrix for several spins as a function of
% field B0 in T and temperature in K.
% [1 1] refers to the nuclei using their atomic mass
% The field is 9.8 T and the temperature is 298 K.
rho0=ThermalEqDensityOp([1 1], 9.8, 298);

O1=400e6*2*pi;  % Carrier frequency
d1=1e-6;        % 1H chemical shift - spin 1
d2=-1.5e-6;     % 1H chemical shift - spin 2
wI1=O1*(1+d1);  % Larmor frequency - spin 1
wI2=O1*(1+d2);  % Larmor frequency - spin 2
J12=20;         % J-coupling

% Free precession Hamiltonian
H0=wI1*I1z + wI2*I2z + 2*pi*J12*(I1x*I2x+I1y*I2y+I1z*I2z);

% Pulse Hamiltonian
% Operator Iy=I1y+I2y is used because the pulse acts on both spins
w1=30e3*2*pi;   % Nutation frequency in rad/s
H1=w1*Iy;       % Hamiltonian in natural units
taup=pi/2/w1;   % Pulse length for pi/2

% Propagator for the pulse Hamiltonian during the length of the pulse
U1=expm(-ai*H1*taup);

% Propagation of the density matrix under the effect of the pulse
% Hamiltonian using the sandwitch formula
rho1=U1'*rho0*U1;

% Parameters for the Free Induction Decay (FID)
aq=3;           % Acquisition times in s
SW=50*2*pi;    % Spectral window in rad/s

% Number of points of the FID in the next power of 2
N=2^round(log2(SW*aq*2*pi));

% Time axis for the FID
t=linspace(0,aq,N);

% Empty vectors for the expectation values along time/propagation of the
% free evolution Hamiltonian. They correspond to the signal along x and y
% assuming quadrature detection.
Ixt=zeros(N,1);
Iyt=zeros(N,1);

% Computation of the FID
for i=1:N
    U=expm(-ai*t(i)*H0);
    rhot=U'*rho1*U;
    
    % Computation of the expectation values
    Ixt(i)=real(trace(Ix*rhot));
    Iyt(i)=real(trace(Iy*rhot));
end

Signal1D(t, Ixt+ai*Iyt,'lb',1,'spectrum', true,'fid',true,'O1', O1, 'unit', 'ppm', 'ph', pi);

%% Two spin system homonuclear case - truncated J Hamiltonian
% Simulation of the spectrum for a pair of J-coupled 1H spins using the
% truncated J Hamiltonian, that is, with only z terms. The 'roof effect' is
% lost

% Angular momentum operators in the Hilbert space of a single spin
Ix=1/2*[0 1;1 0];
Iy=1/2*[0 -ai;ai 0];
Iz=1/2*[1 0;0 -1];

% Angular momentum operators in the Hilbert space of a two-spin system^
% kron computes the Kronecker product.
I1z=kron(Iz, eye(2));
I2z=kron(eye(2),Iz);
I1x=kron(Ix, eye(2));
I2x=kron(eye(2),Ix);
I1y=kron(Iy, eye(2));
I2y=kron(eye(2),Iy);
Ix=I1x+I2x;
Iy=I1y+I2y;
Iz=I1z+I2z;

% Function computing the density matrix for several spins as a function of
% field B0 in T and temperature in K.
% [1 1] refers to the nuclei using their atomic mass
% The field is 9.8 T and the temperature is 298 K.
rho0=ThermalEqDensityOp([1 1], 9.8, 298);

O1=400e6*2*pi;  % Carrier frequency
d1=1e-6;        % 1H chemical shift - spin 1
d2=-1.5e-6;     % 1H chemical shift - spin 2
wI1=O1*(1+d1);  % Larmor frequency - spin 1
wI2=O1*(1+d2);  % Larmor frequency - spin 2
J12=20;         % J-coupling

% Free precession Hamiltonian TRUNCATED
H0=wI1*I1z + wI2*I2z + 2*pi*J12*I1z*I2z;

% Pulse Hamiltonian
% Operator Iy=I1y+I2y is used because the pulse acts on both spins
w1=30e3*2*pi;   % Nutation frequency in rad/s
H1=w1*Iy;       % Hamiltonian in natural units
taup=pi/2/w1;   % Pulse length for pi/2

% Propagator for the pulse Hamiltonian during the length of the pulse
U1=expm(-ai*H1*taup);

% Propagation of the density matrix under the effect of the pulse
% Hamiltonian using the sandwitch formula
rho1=U1'*rho0*U1;

% Parameters for the Free Induction Decay (FID)
aq=3;           % Acquisition times in s
SW=50*2*pi;    % Spectral window in rad/s

% Number of points of the FID in the next power of 2
N=2^round(log2(SW*aq*2*pi));

% Time axis for the FID
t=linspace(0,aq,N);

% Empty vectors for the expectation values along time/propagation of the
% free evolution Hamiltonian. They correspond to the signal along x and y
% assuming quadrature detection.
Ixt=zeros(N,1);
Iyt=zeros(N,1);

% Computation of the FID
for i=1:N
    U=expm(-ai*t(i)*H0);
    rhot=U'*rho1*U;
    
    % Computation of the expectation values
    Ixt(i)=real(trace(Ix*rhot));
    Iyt(i)=real(trace(Iy*rhot));
end

Signal1D(t, Ixt+ai*Iyt,'lb',1,'spectrum', true,'fid',true,'O1', O1, 'unit', 'ppm', 'ph', pi);
