%% Matlab tutorial - Session 4 - Numerical propagation and INEPT
% This document contains the elements presented in a tutorial on MATLAB's
% basics on 30.4.20 for the CRMN members by Quentin Chappuis
%
% Subjects:
%   - Methods for numerical propagation
%   - INEPT pulse sequence, phase cycling and monitoring of multiple spin
%     order

%% Comparison of two propagation methods
% The simulation of the 1H spectrum of ethanol is taken as example
% The pair of 1H spins of the CH2 moiety are labelled a
% The three 1H spins of the CH3 moiety are labelled b.
clear all
ai=sqrt(-1);

% Initial density matrix of 5 1H spins at 9.6T and 298K.
rhoini=ThermalEqDensityOp([1 1 1 1 1], 9.6, 298);
% For more detail, type 'help ThermalEqDensityOp' in the command window

O1=400e6*2*pi;  % Carrier frequency
da=1e-6;        % Chemical shift of spins 1Ha
db=-1.5e-6;     % Chemical shift of spins 1Hb
wa=O1*(1+da);   % Larmor frequecy of spins 1Ha
wb=O1*(1+db);   % Larmor frequecy of spins 1Hb
Jab=6;          % J-coupling between spins 1Ha and 1Hb

% Angular momentum operator in a 5-spin Hilbert space
Ix=I(5, 'x');
Iy=I(5, 'y');
Iz=I(5, 'z');

I1z=I(5,1,'z');
I2z=I(5,2,'z');
I3z=I(5,3,'z');
I4z=I(5,4,'z');
I5z=I(5,5,'z');

% Zeeman Hamiltonian
HZ=wa*(I1z+I2z) + wb*(I3z+I4z+I5z);

% J Hamiltonian
HJab1=2*pi*Jab*(I1z*I3z + I1z*I4z + I1z*I5z);
HJab2=2*pi*Jab*(I2z*I3z + I2z*I4z + I2z*I5z);

% Total Hamiltonian (of free precession)
H0=HZ+HJab1+HJab2;

% Pulse along y on all spins
H1taup=pi/2*Iy;
U1=expm(-ai*H1taup);
rhoini=U1'*rhoini*U1;

% Acquisition parameters
aq=20;
SW=2000*2*pi;
N=2^round(log2(SW*aq/2/pi));
t=((1:N)-1)/(N-1)*aq;
St=zeros(N,1);

% Propagation method 1
tic
for k=1:N
    Ut=expm(-ai*H0*t(k));
    rhot=Ut'*rhoini*Ut;
    St(k)=real(trace(rhot*Ix))+ai*real(trace(rhot*Iy));
end
Tprop_M1=toc;

% Propagation method 2
tic
dt=aq/(N-1);
Udt=expm(-ai*H0*dt);
Udt_inv=Udt';
rhok=rhoini;
for k=1:N
    St(k)=real(trace(rhok*Ix))+ai*real(trace(rhok*Iy));
    rhok=Udt_inv*rhok*Udt;  
end
Tprop_M2=toc;

disp(['Method 2 is ' num2str(Tprop_M1/Tprop_M2) ' times faster than method 1.'])

% Fourier transform after method 2
Signal1D(t,St,'lb',1,'ph',pi, 'spectrum', true, 'o1', O1, 'unit', 'ppm');

%% INEPT
% Simulation of the transfer of polarization from spin I (1H) to spin S
% (13C).
clear all
ai=sqrt(-1);

O1=100e6*2*pi;  % Carrier frequency
dS=1e-6;        % Chemical shift of spin S
wI=O1*4;        % Larmor frequency of spin I
wS=O1*(1+dS);   % Larmor frequency of spin S
JIS=20;         % J coupling
tau=1/4/JIS;    % Free evolution delay for maximum transfer

% Angular momentum operator in a 2-spin Hilbert space
Ix=I(2,1,'x');
Sx=I(2,2,'x');
Iy=I(2,1,'y');
Sy=I(2,2,'y');
Iz=I(2,1,'z');
Sz=I(2,2,'z');

% Initial density matrix of 1H spin and a 13C spin at 9.6T and 298K.
rho0=ThermalEqDensityOp([1 13], 9.6, 298);

% Hamiltonian for free evolution
H0=wI*Iz+wS*Sz + 2*pi*JIS*(Ix*Sx+Iy*Sy+Iz*Sz);

% Propagator of free evolution during tau
U0tau=expm(-ai*H0*tau);

% Pi/2 pulse along x on spin I
H1=Ix*pi/2;
rho1=expm(+ai*H1)*rho0*expm(-ai*H1);

% Free evolution during tau
rho2=U0tau'*rho1*U0tau;

% Pi pulse along x on both spins
H3=(Ix+Sx)*pi;
rho3=expm(+ai*H3)*rho2*expm(-ai*H3);

% Free evolution during tau
rho4=U0tau'*rho3*U0tau;

% Pi/2 pulse along y on spin I and along x on spin S
H5=(Iy+Sx)*pi/2;
rho5=expm(+ai*H5)*rho4*expm(-ai*H5);

% Acquistion parameters
aq=5;
SW=500*2*pi;
N=2^round(log2(SW*aq/2/pi));

% Propation during the FID at the end of the INEPT sequence
t=((1:N)-1)/(N-1)*aq;
dt=aq/(N-1);
St=zeros(N,1);
Udt=expm(-ai*H0*dt);
Udt_inv=Udt';
rhok=rho5;
for k=1:N
    St(k)=real(trace(rhok*Sx))+ai*real(trace(rhok*Sy));
    rhok=Udt_inv*rhok*Udt;  
end

% Propation during the FID without INEPT (simple zg)
St_noINEPT=zeros(N,1);
Hp=pi/2*Sx;
rhok=expm(+ai*Hp)*rho0*expm(-ai*Hp);
for k=1:N
    St_noINEPT(k)=real(trace(rhok*Sx))+ai*real(trace(rhok*Sy));
    rhok=Udt_inv*rhok*Udt;  
end

% Fourier transform for each spectrum
[w, St_INEPT_]=Signal1D(t, St,'lb',1,'ph',-pi/2, 'o1', O1, 'unit', 'ppm');
[w, St_]=Signal1D(t, St_noINEPT,'lb',1,'ph',-pi/2, 'o1', O1, 'unit', 'ppm');

% Comparison of the two signals
plot(w, real(St_INEPT_), w, real(St_))
xlim([0.5 1.5])
xlabel 'Offset / ppm'

%% Refocused INEPT
% Simulation of refocused INEPT with same conditions as above
% Including a phase cycling
clear all
ai=sqrt(-1);

O1=100e6*2*pi;  % Carrier frequency
dS=1e-6;        % Chemical shift of spin S
wI=O1*4;        % Larmor frequency of spin I
wS=O1*(1+dS);   % Larmor frequency of spin S
JIS=20;         % J coupling
tau=1/4/JIS;    % Free evolution delay for maximum transfer

% Angular momentum operator in a 2-spin Hilbert space
Ix=I(2,1,'x');
Sx=I(2,2,'x');
Iy=I(2,1,'y');
Sy=I(2,2,'y');
Iz=I(2,1,'z');
Sz=I(2,2,'z');

% Initial density matrix of 1H spin and a 13C spin at 9.6T and 298K.
rho0=ThermalEqDensityOp([1 13], 9.6, 298);

% Hamiltonian for free evolution
H0=wI*Iz+wS*Sz + 2*pi*JIS*(Ix*Sx+Iy*Sy+Iz*Sz);

% Propagator of free evolution during tau
U0tau=expm(-ai*H0*tau);

% Acquisition parameters
aq=5;
SW=500*2*pi;
N=2^round(log2(SW*aq/2/pi));

t=((1:N)-1)/(N-1)*aq;
dt=aq/(N-1);
St=zeros(N,2);
Udt=expm(-ai*H0*dt);
Udt_inv=Udt';

% Phase lists
ph1_list=[0 2]; % First pulse
ph31_list=[0 2];% Acquisition
% Here I follow Topspin numbering:
% 0 -> x
% 1 -> y
% 2 -> -x
% 3 -> -y

% The pulse sequence is repeated twice including the phase cycling
% j is the increment of the phase cycling
for j=1:2
    % Pi/2 pulse along x and -x on spin I (phase cycling)
    H1=(Ix*cos(pi/2*ph1_list(j)) + Iy*sin(pi/2*ph1_list(j))) *pi/2;
    rho1=expm(+ai*H1)*rho0*expm(-ai*H1);

    % Free evolution during tau
    rho2=U0tau'*rho1*U0tau;

    % Pi pulse along x on both spins
    H3=(Ix+Sx)*pi;
    rho3=expm(+ai*H3)*rho2*expm(-ai*H3);

    % Free evolution during tau
    rho4=U0tau'*rho3*U0tau;

    % Pi/2 pulse along y on spin I and along x on spin S
    H5=(Iy+Sx)*pi/2;
    rho5=expm(+ai*H5)*rho4*expm(-ai*H5);

    % Free evolution during tau
    rho6=U0tau'*rho5*U0tau;

    % Pi pulse along x on both spins
    H7=(Ix+Sx)*pi;
    rho7=expm(+ai*H7)*rho6*expm(-ai*H7);

    % Free evolution during tau
    rho8=U0tau'*rho7*U0tau; 
    
    % Operators for detection incluing phase cycling
    Dx=(Sx*cos(pi/2*ph31_list(j)) + Sy*sin(pi/2*ph31_list(j)));
    Dy=(Sx*sin(pi/2*ph31_list(j)) + Sy*cos(pi/2*ph31_list(j))) ;
    
    % Propation during the FID at the end of the INEPT sequence
    rhok=rho8;
    for k=1:N
        St(k,j)=real(trace(rhok*Dx))+ai*real(trace(rhok*Dy));
        rhok=Udt_inv*rhok*Udt;  
    end  
end

% Propation during the FID without INEPT (simple zg)
S_noINEPT=zeros(N,1);
Hp=pi/2*Sx;
rhok=expm(+ai*Hp)*rho0*expm(-ai*Hp);
for k=1:N
    S_noINEPT(k)=real(trace(rhok*Sx))+ai*real(trace(rhok*Sy));
    rhok=Udt_inv*rhok*Udt;  
end

% Fourier transform for each spectrum
[w, Sph1]=Signal1D(t, St(:,1),'lb',1,'ph',pi, 'o1', O1, 'unit', 'ppm', 'zero', 16*N);
[w, Sph2]=Signal1D(t, St(:,2),'lb',1,'ph',pi, 'o1', O1, 'unit', 'ppm', 'zero', 16*N);
[w, S_noINEPT]=Signal1D(t, S_noINEPT,'lb',1,'ph',-pi/2, 'o1', O1, 'unit', 'ppm', 'zero', 16*N);

% Normalisation factor
m=max(real(S_noINEPT));

figure
subplot(1,2,1)
plot(w, real(Sph1)/m, w, real(Sph2)/m)
xlim([0.75 1.25])
ylim([-0.02 5.2])
xlabel 'Offset / ppm'
ylabel 'Relative intensity / a.u.'
legend({"Phase 1","Phase 2"})

subplot(1,2,2)
plot(w, real((Sph1+Sph2)/2)/m, w, real(S_noINEPT)/m)
xlim([0.75 1.25])
ylim([-0.02 5.2])
xlabel 'Offset / ppm'
legend({"INEPT","No INEPT"})
%% INEPT - Followin two-spin order
clear all
ai=sqrt(-1);

O1=100e6*2*pi;  % Carrier frequency
dS=1e-6;        % Chemical shift of spin S
wI=O1*4;        % Larmor frequency of spin I
wS=O1*(1+dS);   % Larmor frequency of spin S
JIS=20;         % J coupling
tau=1/4/JIS;    % Free evolution delay for maximum transfer

% Angular momentum operator in a 2-spin Hilbert space
Ix=I(2,1,'x');
Sx=I(2,2,'x');
Iy=I(2,1,'y');
Sy=I(2,2,'y');
Iz=I(2,1,'z');
Sz=I(2,2,'z');

% Two spin order operator 2IxSz and 2IySz
twoIxSz=2*Ix*Sz;
twoIySz=2*Iy*Sz;

% Initial density matrix of 1H spin and a 13C spin at 9.6T and 298K.
rho0=ThermalEqDensityOp([1 13], 9.6, 298);
rho0=PolarizedDensityOp(2, [1 0.25]);

% Hamiltonian for free evolution
H0=wI*Iz+wS*Sz + 2*pi*JIS*(Ix*Sx+Iy*Sy+Iz*Sz);

% Propagator of free evolution during tau
U0tau=expm(-ai*H0*tau);

% Pi/2 pulse along x on spin I
H1=Ix*pi/2;
rho1=expm(+ai*H1)*rho0*expm(-ai*H1);

% Discretization of the free evolution during INEPT transfer
M=1000;
IxSzt=zeros(2*M,1);
IySzt=zeros(2*M,1);
Ixt=zeros(2*M,1);
Iyt=zeros(2*M,1);

dt=tau/(M-1);
Udt=expm(ai*H0*dt);
rhok=rho1;
t=linspace(0, 2*tau, 2*M)';

% Propagation
for i=1:2*M
    IxSzt(i)=real(trace(twoIxSz*rhok)/trace(twoIxSz^2));
    IySzt(i)=real(trace(twoIySz*rhok)/trace(twoIySz^2));
    Ixt(i)=real(trace(Ix*rhok)/trace(Ix^2));
    Iyt(i)=real(trace(Iy*rhok)/trace(Iy^2));
    rhok=Udt'*rhok*Udt;
    if i==M
        Up=expm(-ai*pi*(Ix+Sx));
        rhok=Up'*rhok*Up;
    end
end

% Monitoring of the transfer of from single spin order to two spin order
figure
plot(t, sqrt(Ixt.^2+Iyt.^2),'r', 'linewidth', 2)
hold on
plot(t, sqrt(IxSzt.^2+IySzt.^2),'k', 'linewidth', 2)
plot(t, sqrt(IxSzt.^2+IySzt.^2+Ixt.^2+Iyt.^2),'b', 'linewidth', 2)
plot(t,[Ixt Iyt],'r', 'linewidth', 0.5)
plot(t,[IxSzt IySzt],'k', 'linewidth', 0.5)
hold off
title 'Conversion to two-spin order during INEPT'
legend({"Two-spin order","Single spin order","Sum"},'location','southeast')
xlabel 'Time/s'
ylabel 'Operator projections'