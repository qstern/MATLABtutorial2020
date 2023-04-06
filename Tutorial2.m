%% Matlab tutorial - Session 2 - Fitting models to experimental data
% This document contains the elements presented in a tutorial on MATLAB's
% basics on 16.4.20 for the CRMN members by Quentin Chappuis

% The experimental data used in this fit are available on this link:
% https://mycore.core-cloud.net/index.php/s/X6abY5o4JQTrCfD
% Download the file named 150.zip and unzip it in the locatio of your
% choice. Modify the the variable Path (in the next section 'Importing NMR
% data'

clear all

%% Importing NMR data
% Set the path to the data
% See Session 1 of the tutorial for more details
Path='C:\Bruker\TopSpin4.0.7\examdata\20200107_Dissolution\150';
NMRdata=rbnmr(Path);
NMRdata=NMRdata{1};

%% Playing with NMR data
% See Session 1 of the tutorial for more details
w=NMRdata.XAxis;    % Takes the frequency axis created by rbnmr and 
                    % stores in a vector called w
S2D=NMRdata.Data';  % Saves the 2D signal saved in Data

%% Region selection and integration
% See Session 1 of the tutorial for more details
% Here I define two regions (in point number) where I want to integrate my
% signal.
SB=[3100 3500];     % Signal boundaries
NB=[13000 15000];   % Noise boundaries 

figure
subplot(1, 3, 1)
plot(S2D(:,1))
title('Whole spectrum')
xlabel('Freq/ppm')
ylabel('Relative intensity/a.u.')

subplot(1, 3, 2)
plot(w(SB(1):SB(2)), S2D(SB(1):SB(2),1))
title('Signal')
xlabel('freq/ppm')
set (gca,'xdir','reverse')

subplot(1, 3, 3)
plot(w(NB(1):NB(2)), S2D(NB(1):NB(2),1))
title('Noise')
xlabel('freq/ppm')

%% Remove empty spectrum
% See Session 1 of the tutorial for more details
N=size(S2D,2);  
while sum(S2D(:,N))==0
     N=N-1;
end
S2D=S2D(:,1:N);

%% Integration
% See Session 1 of the tutorial for more details
D1=NMRdata.Acqus.D(2);

% Creates a time vector
t=linspace(0,D1*(size(S2D,2)-1),size(S2D,2));

% Sums the signal over the range defined earlier
I=sum(S2D(SB(1):SB(2),:),1)*abs(w(2)-w(1));

% Creates a figure with the signal integral
figure
plot(t, I/max(I))
grid on
xlim([min(t) max(t)])
ylim([-0.05 1.05])
xlabel 'time/s'
ylabel 'Relative signal integral/a.u.'

%% Plot a single spectrum
figure
plot(w(SB(1):SB(2)), S2D(SB(1):SB(2),2))
w_=w(SB(1):SB(2));
S_=S2D(SB(1):SB(2),2);

%% Fit each spectrum
% The first N spectra of the decay are fitted. For each of them, a
% quadruplet if fitted with 4 independent parameters:
% I0 - Signal integral
% J  - J-coupling in ppm
% g  - FWHM of the peaks
% w0 - Centre of the quadruplet in ppm
N=100;

% These four vectors will contain the fit results for each independent
% parameter. It has N rows (for the N spectra to be fitted) and one column.
I0=zeros(N,1);
J=zeros(N,1);
g=zeros(N,1);
w0=zeros(N,1);

% These four vectors will contain the errors on the fit for each 
% independent parameter. It has N rows (for the N spectra to be fitted) and
% one column.
eI0=zeros(N,2);
eJ=zeros(N,2);
eg=zeros(N,2);
ew0=zeros(N,2);

% The fit algorithm needs a good guess of the initial parameters. Here we
% put parameters that we read from the spectrum. Note that we make the
% process much faster by redefining the StartPoint after each fit. Indeed,
% we use the results of the previous fit as the StartPoint of the next. 
StartPoint=[1e10 0.3 0.025 180.4];

tic
for i=1:N
    % Calls the function that performs the Lorentzian deconvution (defined
    % below)
    [fitresult, gof] = QuadrupletFit(w_, S2D(SB(1):SB(2),i),StartPoint, false);
    
    % Stores the result of the fit
    I0(i)=fitresult.I0;
    J(i) =fitresult.J;
    g(i) =fitresult.g;
    w0(i)=fitresult.w0;
    
    % Gets the confidence interval of each parameter (with 95 % confidence)
    error=confint(fitresult);
    
    % Stores the confidence interval (that I loosely called 'errors')
    eI0(i,:)=error(:,1);
    eJ(i,:)=error(:,2);
    eg(i,:)=error(:,3);
    ew0(i,:)=error(:,4);
    
    % Shows the progress of the process in the command window. num2str
    % converts the integers i and N into characters. The the [] concatenate
    % the characters and disp displays them in the command window.
    disp([num2str(i) '/' num2str(N)])
    
    % Define the starting point of the next fit by the result of the
    % previous one.
    StartPoint=[fitresult.I0 fitresult.J fitresult.g fitresult.w0];
end
toc

%% Plot fit results
figure
subplot(2,2,1)
plot(t(1:N),I0)
hold on
plot(t(1:N), I(1:N),'x')
plot(t(1:N), eI0(:,2), 'r:')
plot(t(1:N), eI0(:,1), 'r:')
hold off
xlabel 't', ylabel 'I0'
grid on

subplot(2,2,2)
plot(t(1:N),J)
hold on
plot(t(1:N), eJ(:,2), 'r:')
plot(t(1:N), eJ(:,1), 'r:')
hold off
xlabel 't', ylabel 'J/ppm'
grid on

subplot(2,2,3)
plot(t(1:N),g)
hold on
plot(t(1:N), eg(:,2), 'r:')
plot(t(1:N), eg(:,1), 'r:')
hold off
xlabel 't', ylabel 'FWHM/ppm'
grid on

subplot(2,2,4)
plot(t(1:N),w0)
hold on
plot(t(1:N), ew0(:,2), 'r:')
plot(t(1:N), ew0(:,1), 'r:')
hold off
xlabel 't', ylabel 'w0/ppm'
grid on

%% Comparison of fit
% These two fits are the monoexponential decay fits of the same data. In
% the first case, we fit the intensity obtained by Lorentzian
% deconvolution. In the second, we fit the numerical integral of the
% signal.
[FitRes1, gof1]=MonoexponentialFit(t(1:N), I0')
[FitRes2, gof2]=MonoexponentialFit(t(1:N), I(1:N))
% Using the fitted intensity slightly improves the monoexponential fit. In
% case where the signal to noise ratio is low, fitting the fitted intensity
% can bring a tremendous improvement.

%% Creating 2D spectra to be simulated
% Creates a Lorentzian function with N points
N=100;
w=linspace(0, 100, N);  % linspace(a,b,N) creates a row vector with N
                        % elements equality spread over a and b
w0=50;                  % Centre of the Lorentzian
FWHM=10;                % Full width at half maximum of the Lorentzian
S=1/pi*(FWHM/2)./((FWHM/2)^2+(w-w0).^2);
plot(w, S)

% Computes a monoexponential decay with time constant T1.
T1=5;                   % decay time constant
t=linspace(0,15,12);    % row vector corresponsing to time
S2D=S'*exp(-t/T1);

% This codes simulates a 2D spectrum with gaussian noise.
S2DwithNoise=S2D+normrnd(0,0.001,size(S2D));
S2DwithNoise=S2DwithNoise';
figure, waterfall(w,t,S2DwithNoise)

% Fits the simulated spectrum.
FitRes2D = Fit2D(w, t, S2DwithNoise)

%% functions
function [fitresult, gof] = QuadrupletFit(w_, S_, StartPoint, PlotFig)
    %CREATEFIT(W_,S_)
    %  Create a fit.
    %
    %  Data for 'Quadruplet deconvolution' fit:
    %      X Input : w_
    %      Y Output: S_
    %  Output:
    %      fitresult : a fit object representing the fit.
    %      gof : structure with goodness-of fit info.
    %
    %  See also FIT, CFIT, SFIT.

    %  Auto-generated by MATLAB on 16-Apr-2020 15:34:54


    %% Fit: 'Quadruplet deconvolution'.
    [xData, yData] = prepareCurveData( w_, S_ );

    % Set up fittype and options.
    ft = fittype( '3/8*I0/pi*((g/2)/((g/2)^2+(w0-w+J*1/2)^2))+  3/8*I0/pi*((g/2)/((g/2)^2+(w0-w-J*1/2)^2))+  1/8*I0/pi*((g/2)/((g/2)^2+(w0-w+J*3/2)^2))+  1/8*I0/pi*((g/2)/((g/2)^2+(w0-w-J*3/2)^2))', 'independent', 'w', 'dependent', 'S' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.MaxIter = 4000;
                        %I0   J   g    w0
    opts.Lower =        [0    0.2 0     180.3];
    opts.StartPoint =   StartPoint;
    opts.Upper =        [Inf  0.4 1     180.6];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    % Plot fit with data.
    if PlotFig==true
        figure( 'Name', 'Quadruplet deconvolution' );
        h = plot( fitresult, xData, yData );
        legend( h, 'S vs. w', 'Quadruplet deconvolution', 'Location', 'NorthEast' );
        % Label axes
        xlabel w
        ylabel S
        grid on       
    end
end

function [fitresult, gof] = MonoexponentialFit(t, I)
    %CREATEFIT(T,I)
    %  Create a fit.
    %
    %  Data for 'Monoexponential fit' fit:
    %      X Input : t
    %      Y Output: I
    %  Output:
    %      fitresult : a fit object representing the fit.
    %      gof : structure with goodness-of fit info.
    %
    %  See also FIT, CFIT, SFIT.

    %  Auto-generated by MATLAB on 16-Apr-2020 16:20:49


    %% Fit: 'Monoexponential fit'.
    [xData, yData] = prepareCurveData( t, I );

    % Set up fittype and options.
    ft = fittype( 'I0*exp(-t/T1)', 'independent', 't', 'dependent', 'I' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0 0];
    opts.MaxIter = 4000;
    opts.StartPoint = [5000000000 100];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    % Plot fit with data.
    figure( 'Name', 'Monoexponential fit' );
    h = plot( fitresult, xData, yData );
    legend( h, 'I vs. t', 'Monoexponential fit', 'Location', 'NorthEast' );
    % Label axes
    xlabel t
    ylabel I
    grid on

end

function [fitresult, gof] = Fit2D(x, y, Z)
    %CREATEFIT(x,y,Z)
    %  Create a fit.
    %
    %  Data for '2D fit' fit:
    %      X Input : frequency
    %      Y Input : time
    %      Z Output: Signal
    %  Output:
    %      fitresult : a fit object representing the fit.
    %      gof : structure with goodness-of fit info.
    %
    %  See also FIT, CFIT, SFIT.

    %  Auto-generated by MATLAB on 21-Apr-2020 12:40:07


    %% Fit: 'untitled fit 1'.
    [xData, yData, zData] = prepareSurfaceData( x, y, Z );

    % Set up fittype and options.
    ft = fittype( 'exp(-t/T1)*I0/pi*(g/2)/((g/2)^2+(w-w0)^2)', 'independent', {'w', 't'}, 'dependent', 'z' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
                     % I0  T1  g   w0 
    opts.Lower =      [0   0   0   0];
    opts.StartPoint = [1   1   1   1];
    opts.Upper =      [Inf Inf Inf Inf];

    % Fit model to data.
    [fitresult, gof] = fit( [xData, yData], zData, ft, opts );

    % Plot fit with data.
    figure( 'Name', '2D fit' );
    h = plot( fitresult, [xData, yData], zData );
    legend( h, '2D fit', 'Signal vs. w, t', 'Location', 'NorthEast' );
    % Label axes
    xlabel w
    ylabel t
    zlabel S2DwithNoise
    grid on

end
