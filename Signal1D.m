function [w, S]=Signal1D(varargin)
    % Quentin Chappuis 29.4.2020 - HMRlab	
    %
    % Returns the Fourier transform from a free induction decay (FID)
    %
    % ARGUMENT 1: Time axis in s
    % ARGUMENT 2: Complex FID
    %
    % OPTIONS: (not case sensitive)
    %
    % 'SPECTRUM' and 'FID' - set to true to plot the spectrum and FID,
    % respectively (use boolean variable)
    %
    % 'ZERO' - zero filling. The value must be an integer. If it is smaller
    % than the size of the FID, the FID is truncated.
    %
    % 'PH' - add a phase to the complex signal. The value must be a scalar.
    %
    % 'O1' - set a carrrier frequency that will be mixed with the complex
    % FID
    %
    % 'UNIT' - change the unit of the frequency axis. Default: Hz. Options
    % are 'ppm' and 'rad' (rad for rad/s). Ppm is only available if O1 is
    % defined.
    %
    % 'LB' - add line broadening to the FIDin Hz. Adds an apodization 
    % window using a monoexponential decay function.
    %
    % Example 1: Signal1D(t, S), where S is a complex signal
    %
    % Example 2: Signal1D(t, Sx +sqrt(-1)*Sy), where Sx and Sy are real
    %
    % Example 3: Signal1D(t, S, 'spectrum', true, 'o1', 500e6*2*pi,
    % 'unit','ppm', 'lb', 1), plot the spectrum with O1 is 500 MHz, with
    % units of ppm and a line broadening of 1 Hz.
    
    k=size(varargin,2);
    t=varargin{1};
    I=varargin{2};
    aq=max(t);          % Acquisition lenght [s]
    N=max(size(t));     % Number of points
    N0=N;
    ShowSpectrum=false;
    ShowFID=false;
    Ph0=0;
    O1=0;
    unit='hz';
    lb=0;
    Xlim=0;
    
    if k==3
        opt=varargin{3};
        if isfield(opt,'zero')
            N0=opt.zero;
        end
        if isfield(opt,'spectrum')
            ShowSpectrum=opt.spectrum;
        end
        if isfield(opt,'fid')
            ShowFID=opt.fid;
        end
        if isfield(opt,'ph')
            Ph0=opt.ph;
        end
        if isfield(opt,'o1')
            O1=opt.o1;
        end
        if isfield(opt,'unit')
            unit=opt.unit;
        end
        if isfield(opt,'lb')
            lb=opt.lb;
        end
        if isfield(opt,'xlim')
            Xlim=opt.xlim;
        end
    elseif k>2
        for i=3:k-1
            if round(i/2)~=i/2
               if strcmpi(varargin{i}, 'zero') 
                   N0=varargin{i+1};
               elseif strcmpi(varargin{i}, 'spectrum')
                   ShowSpectrum=varargin{i+1};
               elseif strcmpi(varargin{i}, 'fid')
                   ShowFID=varargin{i+1};
               elseif strcmpi(varargin{i}, 'ph')
                   Ph0=varargin{i+1};
               elseif strcmpi(varargin{i}, 'o1') 
                   O1=varargin{i+1};    
               elseif strcmpi(varargin{i}, 'unit')
                   unit=varargin{i+1}; 
               elseif strcmpi(varargin{i}, 'lb')
                   lb=varargin{i+1}; 
               elseif strcmpi(varargin{i}, 'xlim')
                   Xlim=varargin{i+1}; 
               else 
                   disp(['The parameter "' varargin{i} '" is not known'])
               end
            end
        end
    end
    
    % Phase correction, offset mixing and line broadening
    t=linspace(0,aq,N)'; 
    I=I.*exp(sqrt(-1)*(Ph0+t*O1)-t*pi*lb);
    
    % Creating frequency axis
    SW = N/aq/2;    % Sectral width before zero filling [Hz]
    w  = linspace(+SW,-SW,N0); % frequency axis with zero filling
    
    if strcmpi(unit,'rad')
        w=w*2*pi;
    elseif strcmpi(unit,'ppm')
        w=w*2*pi/O1*1e6;
    end
    
    % Fourier transform
    S = fftshift(fft(I,N0))/N0;
    
    % Plots the FID
    if ShowFID==true
        figure
        t=linspace(0,aq,N)'; 
        plot(t, real(I), t,imag(I))
        title('FID')
        xlabel('time / s'), xlim([0 aq])
        ylabel('Signal intensity')
        grid on
    end

    % Plots the spectrum
    if ShowSpectrum==true
        figure
        plot(w, real(S))
        if strcmpi(unit,'hz')
            xlabel('freq / Hz')
        elseif strcmpi(unit,'ppm')
            xlabel('offset / ppm')
        elseif strcmpi(unit,'rad/s')
            xlabel('freq / rad.s_{-1}')
        else
            disp('Unknown frequency unit')
        end
        title('Spectrum')
        if Xlim == 0
            xlim([min(w) max(w)])
        else
            xlim(Xlim)
        end
        minS=min(real(S));
        maxS=max(real(S));
        diff=maxS-minS;
        ylim([minS-0.05*diff maxS+0.05*diff])
        set(gca, 'xdir', 'reverse')
    end
end