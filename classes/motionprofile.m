classdef motionprofile < handle
    %MOTIONPROFILE represents a motion-load profile to be used in a
    %geareval instance.
    
    properties
        description                         % text description of load profile
        time                                % time vector in seconds
        freq                                % frequency vector in radians
        angle                               % position vector in radians
        anglefft                            % Single-side band FFT of angle
        anglevel_unfiltered                 % UNFILTERED velocity vector in rad/s 
        anglevel                            % velocity vector in rad/s   
        angleaccel_unfiltered               % UNFILTERED acceleration vector in rad/s²
        angleaccel                          % acceleration vector in rad/s²
        load                                % load vector in Nm
        load_unfiltered                     % UNFILTERED load vector in Nm
        loadvel                             % load velocity vector in Nm/s
        loadvel_unfiltered                  % UNFILTERED load velocity vector in Nm/s
        loadaccel                           % load acceleration vector in Nm/s²
        loadaccel_unfiltered                % UNFILTERED load acceleration vector in Nm/s²
        rmsload
        peakload
        period                              % period in seconds
        Npoints                             % number of profile points
        poslp                               % position vector lowpass filter frequency [Hz]
        loadlp                              % load vector lowpass filter frequency [Hz]
               
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
            addParameter(p, 'poslp', [], checkDoubleScalar);
            addParameter(p, 'loadlp', [], checkDoubleScalar);
            
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

            if ~isempty(p.Results.poslp)
                obj.poslp = p.Results.poslp;
            else
                obj.poslp = 1/obj.time(end)*7; % 7x the fundamental frequency
            end

            if ~isempty(p.Results.loadlp)
                obj.loadlp = p.Results.loadlp;
            else
                obj.loadlp = 1/obj.time(end)*4; % 4x the fundamental frequency
            end

            obj.angle = positionvec(:);
            obj.load = loadvec(:);
            obj.description = description;
            obj.Npoints = length(obj.time);
            
            timestep = obj.period/(obj.Npoints-1);
            fs = 1/timestep;
            
            % [b1,a1] = butter(min(obj.Npoints/3-1,20),obj.poslp/(fs/2)); % high order BW filter for filtfilt
            % [b2,a2] = butter(min(obj.Npoints/3-1,20),obj.loadlp/(fs/2)); % high order BW filter for filtfilt

            obj.angle = positionvec(:);
            [obj.freq, obj.anglefft] = obj.FFT(obj.angle);

            obj.anglevel_unfiltered = gradient(obj.angle,timestep);  
            obj.anglevel = lowpass(obj.anglevel_unfiltered,obj.poslp,fs);
            %obj.anglevel = filtfilt(b1,a1, obj.anglevel_unfiltered);

            obj.angleaccel_unfiltered = gradient(obj.anglevel_unfiltered,timestep);
            obj.angleaccel = lowpass(gradient(obj.anglevel,timestep),obj.poslp,fs);
            %obj.angleaccel = filtfilt(b1,a1,gradient(obj.anglevel,timestep));
            
            obj.load_unfiltered = loadvec(:);
            obj.load = lowpass(obj.load_unfiltered,obj.loadlp,fs);
            %obj.load = filtfilt(b2,a2,obj.load_unfiltered);

            obj.loadvel_unfiltered = gradient(obj.load_unfiltered, timestep);
            obj.loadvel = lowpass(gradient(obj.load, timestep),obj.loadlp,fs);
            %obj.loadvel = filtfilt(b2,a2,gradient(obj.load, timestep));

            obj.loadaccel_unfiltered = gradient(obj.loadvel_unfiltered, timestep);
            obj.loadaccel = lowpass(gradient(obj.loadvel, timestep),obj.loadlp,fs);
            %obj.loadaccel = filtfilt(b2,a2,gradient(obj.loadvel, timestep));

            obj.rmsload = rms(obj.load);
            if abs(max(obj.load)) >= abs(min(obj.load))
                obj.peakload = max(obj.load);
            else
                obj.peakload = min(obj.load);
            end

        end
        function f=plot(obj)
            f=figure('Name', char(obj.description));
            
            subplot(6,1,1);
            plot(obj.time,obj.angle);
            grid on
            ylabel({'Angle', '(rad)'})
            %title(['Motion-load Profile: ' char(obj.description)])
            
            subplot(6,1,2);
            plot(obj.time,obj.anglevel_unfiltered);
            hold on
            plot(obj.time,obj.anglevel);
            hold off
            grid on
            ylim(1.5*[min(obj.anglevel) max(obj.anglevel)])
            ylabel({'Velocity', '(rad/s)'})
            
            subplot(6,1,3);
            plot(obj.time, obj.angleaccel_unfiltered);
            hold on
            plot(obj.time,obj.angleaccel);
            hold off
            grid on
            ylim(1.5*[min(obj.angleaccel) max(obj.angleaccel)])
            ylabel({'Acceleration' '(rad/s²)'})
            
            subplot(6,1,4);
            plot(obj.time,obj.load_unfiltered);
            hold on
            plot(obj.time,obj.load);
            hold off
            grid on
            %ylim(2*[min(obj.load) max(obj.load)])
            ylabel({'Load moment' '(Nm)'})

            subplot(6,1,5);
            plot(obj.time,obj.loadvel_unfiltered);
            hold on
            plot(obj.time,obj.loadvel);
            hold off
            grid on
            ylim(1.5*[min(obj.loadvel) max(obj.loadvel)])
            ylabel({'Load moment' 'velocity' '(Nm/s)'})

            subplot(6,1,6);
            plot(obj.time,obj.loadaccel_unfiltered);
            hold on
            plot(obj.time,obj.loadaccel);
            hold off
            grid on
            ylim(1.5*[min(obj.loadaccel) max(obj.loadaccel)])            
            ylabel({'Load moment' 'acceleration' '(Nm/s²)'})
            xlabel('Time (s)')
        
            % sgtitle()            
            shg

            legend({'Unfiltered', 'Filtered'},Orientation="horizontal")
            legend("Position", [0.59852,0.050059,0.29472,0.025316])
        end

        function fftplot=analyseSpectrum(obj)
            
            %angle
            [f,~, am,~, ~, cpow] = obj.FFT(obj.angle);
            
            figure;
            fftplot=gcf;
            subplot(3,2,1);
            yyaxis left
            bar(f,am, EdgeColor="none", FaceColor=	"#0072BD") 
            title("Load motion")
            xlabel("f (Hz)")
            ylabel("Angle (rad)")            
            hold on
            %bar(f,amfil,'cyan')
            hold off

            yyaxis right
            plot(f,cpow/cpow(end),LineWidth=1)
            hold on
            %plot(f,cpowfil/cpow(end))
            %ylabel("Relative cumulative power")
            hold off
                       
            %angular velocity
            [f,~, am,~, ~, cpow] = obj.FFT(obj.anglevel_unfiltered);
            [f,~, amfil,~, ~, cpowfil] = obj.FFT(obj.anglevel);
            
            subplot(3,2,3);
            yyaxis left
            bar(f,am, EdgeColor="none", FaceColor=	"#0072BD") 
            %title("Angular Velocity")
            xlabel("f (Hz)")
            ylabel("Velocity (rad/s)")            
            hold on
            bar(f,amfil,'cyan', EdgeColor="none")
            hold off

            yyaxis right
            plot(f,cpow/cpow(end),LineWidth=1)
            hold on
            plot(f,cpowfil/cpow(end),LineWidth=1)
            %ylabel("Relative cumulative power")
            hold off
            
            legend({'Unfilt. ampl.', 'Filt. ampl.','Unfilt. cum. power', 'Filt. cum. power'},Orientation="horizontal")
            legend("Position", [0.044494,0.014913,0.94813,0.032755]);

            %angular acceleration
            [f,~, am,~, ~, cpow] = obj.FFT(obj.angleaccel_unfiltered);
            [f,~, amfil,~, ~, cpowfil] = obj.FFT(obj.angleaccel);
            
            subplot(3,2,5);
            yyaxis left
            bar(f,am, EdgeColor="none", FaceColor=	"#0072BD") 
            %title("Angular Acceleration")
            xlabel("f (Hz)")
            ylabel("Acceleration (rad/s²)")            
            hold on
            bar(f,amfil,'cyan', EdgeColor="none")
            hold off

            yyaxis right
            plot(f,cpow/cpow(end),LineWidth=1)
            hold on
            plot(f,cpowfil/cpow(end),LineWidth=1)
            %ylabel("Relative cumulative power")
            hold off

            %load
            [f,~, am,~, ~, cpow] = obj.FFT(obj.load_unfiltered);
            [f,~, amfil,~, ~, cpowfil] = obj.FFT(obj.load);
            
            subplot(3,2,2);
            yyaxis left
            bar(f,am, EdgeColor="none", FaceColor=	"#0072BD") 
            title("Load moment")
            xlabel("f (Hz)")
            ylabel("Load moment (Nm)")            
            hold on
            bar(f,amfil,'cyan', EdgeColor="none")
            hold off

            yyaxis right
            plot(f,cpow/cpow(end),LineWidth=1)
            hold on
            plot(f,cpowfil/cpow(end),LineWidth=1)
            ylabel("Relative cumulative power")
            hold off

            %load velocity
            [f,~, am,~, ~, cpow] = obj.FFT(obj.loadvel_unfiltered);
            [f,~, amfil,~, ~, cpowfil] = obj.FFT(obj.loadvel);
            
            subplot(3,2,4);
            yyaxis left
            bar(f,am, EdgeColor="none", FaceColor=	"#0072BD") 
            %title("Load velocity")
            xlabel("f (Hz)")
            ylabel({'Load moment' 'velocity (Nm/s)'})            
            hold on
            bar(f,amfil,'cyan', EdgeColor="none")
            hold off

            yyaxis right
            plot(f,cpow/cpow(end),LineWidth=1)
            hold on
            plot(f,cpowfil/cpow(end),LineWidth=1)
            ylabel("Relative cumulative power")
            hold off;

            %load acceleration
            [f,~, am,~, ~, cpow] = obj.FFT(obj.loadaccel_unfiltered);
            [f,~, amfil,~, ~, cpowfil] = obj.FFT(obj.loadaccel);
            
            subplot(3,2,6);
            yyaxis left
            bar(f,am, EdgeColor="none", FaceColor=	"#0072BD") 
            %title("Load acceleration")
            xlabel("f (Hz)")
            ylabel({'Load moment' 'acceleration (Nm/s²)'})            
            hold on
            bar(f,amfil,'cyan', EdgeColor="none")
            hold off

            yyaxis right
            plot(f,cpow/cpow(end),LineWidth=1)
            hold on
            plot(f,cpowfil/cpow(end),LineWidth=1)
            ylabel("Relative cumulative power")
            hold off

         end
    end
    methods (Access=private)
        function [f,X, Mag, Phase,Pow,CumPow] = FFT(obj, signal)
            L=length(signal);
            signal = signal(:);
            Fs = (L-1)/obj.period ;           % Sampling frequency
            Ts=1/Fs;
            t = 0:Ts:(L-1)*Ts;    %time vector

            N = 2^nextpow2(L); % FFT points
            
            % upsample to N using linear interpollation
            Fs_up = (N-1)/obj.period;
            Ts_up = 1/Fs_up;
            t_up = 0:Ts_up:(N-1)*Ts_up;   
            signal_upsampled = interp1(t,signal', t_up)';

            % figure
            % plot(t,signal)
            % hold on
            % plot(t_up, signal_upsampled)
            % hold off

            % FFT
            X = fft(signal_upsampled,N);
            SSB = X(1:N/2);
            SSB(2:end) = 2*SSB(2:end);
            f = (0:N/2-1)*(Fs_up/N);

            Mag = abs(SSB/N);
            Phase = angle(SSB);
            Pow =Mag.^2;
            CumPow = cumsum(Pow);

        end

    end
end

