classdef motionprofile < handle
    %MOTIONPROFILE represents a motion-load profile to be used in a
    %geareval instance.
    
    properties
        description                         % text description of load profile
        time                                % time vector in seconds
        angle                               % position vector in radians
        anglevel                            % velocity vector in rad/s            
        angleaccel                          % acceleration vector in rad/s²
        load                                % load vector in Nm
        loadvel                             % load velocity vector in Nm/s
        loadaccel                           % load acceleration vector in Nm/s²
        rmsload
        peakload
        period                              % period in seconds
        Npoints                             % number of profile points
        poslp                               % position vector lowpass filter frequency [Hz]
        loadlp                              % load vector lowpass filter frequency [Hz]
        psdplot (1,1) matlab.ui.Figure      % Power Spectral density plot of input signal
        
    end
    
    methods
        function obj = motionprofile(description, timevec, positionvec, loadvec, varargin)
            %MOTIONPROFILE Construct an instance of this class
            %   Detailed explanation goes here
            arguments (Input)                
                description {mustBeTextScalar}
                timevec double
                positionvec double {mustBeVector}
                loadvec double {mustBeVector}
            end
            arguments (Repeating)
                varargin
            end
            
            assert(length(positionvec) == length(loadvec),'positionvec input must be of equal length as loadvec')

            p = inputParser;

            checkDoubleScalar = @(x) isscalar(x)&&isa(x,'double') ;
            addParameter(p, 'period', [], checkDoubleScalar);
            addParameter(p, 'poslp', 6, checkDoubleScalar);
            addParameter(p, 'loadlp', 3, checkDoubleScalar);
            
            parse(p,varargin{:});

            if ~isempty(p.Results.period)
                obj.period = p.Results.period;
                N = length(positionvec);
                temptimevec = linspace(0,obj.period,N);
            else
                assert(length(positionvec) == length(timevec), 'timevec input must be of equal length as positionvec and loadvec')
                temptimevec = timevec(:) - timevec(1);
                obj.period = temptimevec(end) - temptimevec(1);
            end
            obj.time = temptimevec(:);

            obj.poslp = p.Results.poslp;
            obj.loadlp = p.Results.loadlp;

            obj.angle = positionvec(:);
            obj.load = loadvec(:);
            obj.description = description;
            obj.Npoints = length(obj.time);
            
            timestep = obj.period/(obj.Npoints-1);
            fs = 1/timestep;
            
            obj.angle = positionvec(:);
            obj.anglevel = lowpass(gradient(obj.angle,timestep),obj.poslp,fs);                
            obj.angleaccel = lowpass(gradient(obj.anglevel,timestep),obj.poslp,fs);

            obj.load = lowpass(loadvec(:),obj.loadlp,fs);
            obj.loadvel = lowpass(gradient(obj.load, timestep),obj.loadlp,fs);
            obj.loadaccel = lowpass(gradient(obj.loadvel, timestep),obj.loadlp,fs);

            obj.rmsload = rms(obj.load);
            if abs(max(obj.load)) >= abs(min(obj.load))
                obj.peakload = max(obj.load);
            else
                obj.peakload = min(obj.load);
            end

        end
        function f=plot(obj)
            f=figure('Name', char(obj.description));
            
            subplot(4,1,1);
            plot(obj.time,obj.angle);
            grid on
            ylabel({'Angle', '(rad)'})
            title(['Motion-load Profile: ' char(obj.description)])
            
            subplot(4,1,2);
            plot(obj.time,obj.anglevel);
            grid on
            ylabel({'Velocity', '(rad/s)'})
            
            subplot(4,1,3);
            plot(obj.time, obj.angleaccel);
            grid on
            ylabel({'Acceleration' '(rad/s^2)'})
            
            subplot(4,1,4);
            plot(obj.time,obj.load);
            grid on
            ylabel({'Load moment' '(Nm)'})
            xlabel('Time (s)')
        
            % sgtitle()            
            shg
        end

        function showSpectrum(obj)
            obj.psdplot = figure();
            [pxx,f] = periodogram(obj.angle, hamming(obj.Npoints), 1024, obj.Npoints/(obj.time(end) - obj.time(1)),"onesided","psd");
            plot(gca,f,10*log10(pxx));
            xlabel("Frequency [Hz]");
            ylabel("Power spectral density [dB/Hz]")
            title(string(obj.description) + ": " + "Power spectral density of Excitation", Interpreter="none")
            obj.psdplot.Visible = true;
        end
    end
end

