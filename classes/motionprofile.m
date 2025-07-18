classdef motionprofile < handle
    %MOTIONPROFILE represents a motion-load profile to be used in a
    %geareval instance.
    
    properties
        description                         % text description of load profile
        time                                % time vector in seconds
        freqwrap                            % wrapped frequency vector in Hz (incl negative frequencies at end)
        freqssb                             % SSB frequency vector in Hz
        angle                               % position vector in radians
        angle_unfiltered                    % UNFILTERED position vector in radians
        anglefft                            % FFT of angle
        anglevel_unfiltered                 % UNFILTERED velocity vector in rad/s 
        anglevel                            % velocity vector in rad/s   
        angleaccel_unfiltered               % UNFILTERED acceleration vector in rad/s²
        angleaccel                          % acceleration vector in rad/s²
        load                                % load vector in Nm
        load_unfiltered                     % UNFILTERED load vector in Nm
        loadfft                             % FFT of load
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
        fourier_ord                         % Order of Fourier series approx., if applicable. 2-vector with first element position and 2nd load                            
               
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
            checkintergerNonnegScalar = @(x) isvector(x) && all(x>=0);
            addParameter(p, 'period', [], checkDoubleScalar);
            addParameter(p, 'poslp', [], checkDoubleScalar);
            addParameter(p, 'loadlp', [], checkDoubleScalar);
            addParameter(p, 'fourier', [], checkintergerNonnegScalar); 
            
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

            if ~isempty(p.Results.fourier)
                obj.fourier_ord = floor(p.Results.fourier);
            else
                obj.fourier_ord = 0; % do not use Fourier series
            end


            obj.description = description;
            obj.Npoints = length(obj.time);
            
            timestep = obj.period/(obj.Npoints-1);
            fs = 1/timestep;
            
            % [b1,a1] = butter(min(obj.Npoints/3-1,20),obj.poslp/(fs/2)); % high order BW filter for filtfilt
            % [b2,a2] = butter(min(obj.Npoints/3-1,20),obj.loadlp/(fs/2)); % high order BW filter for filtfilt

            obj.angle_unfiltered = positionvec(:);
            obj.anglevel_unfiltered = obj.first_derivative_5pt(obj.angle_unfiltered,timestep); 
            obj.angleaccel_unfiltered = obj.second_derivative_5pt(obj.angle_unfiltered,timestep);

            obj.load_unfiltered = loadvec(:);
            obj.loadvel_unfiltered = obj.first_derivative_5pt(obj.load_unfiltered, timestep);
            obj.loadaccel_unfiltered = obj.second_derivative_5pt(obj.load_unfiltered, timestep);
            
            % do FFT
            [obj.freqssb, obj.anglefft, ~,~,~,~, obj.freqwrap] = obj.FFT(obj.angle_unfiltered);
            [obj.freqssb, obj.loadfft] = obj.FFT(obj.load_unfiltered);

            if obj.fourier_ord == 0
                obj.angle = obj.angle_unfiltered;
                 
                obj.anglevel = lowpass(obj.anglevel_unfiltered,obj.poslp,fs);
                %obj.anglevel = filtfilt(b1,a1, obj.anglevel_unfiltered);
    
                obj.angleaccel = lowpass(obj.angleaccel_unfiltered,obj.poslp,fs);
                %obj.angleaccel = filtfilt(b1,a1,gradient(obj.anglevel,timestep));
                
                obj.load = lowpass(obj.load_unfiltered,obj.loadlp,fs);
                %obj.load = filtfilt(b2,a2,obj.load_unfiltered);    
                
                obj.loadvel = lowpass(obj.first_derivative_5pt(obj.load,timestep),obj.loadlp,fs);
                %obj.loadvel = filtfilt(b2,a2,gradient(obj.load, timestep));
    
                obj.loadaccel = lowpass(obj.second_derivative_5pt(obj.load,timestep),obj.loadlp,fs);
                %obj.loadaccel = filtfilt(b2,a2,gradient(obj.loadvel, timestep));

            else % use fourier
                % position and derivatives
                obj.angle = obj.eval_dfourier(obj.freqwrap, obj.anglefft,  obj.fourier_ord(1),0, obj.time);
                obj.anglevel = obj.eval_dfourier(obj.freqwrap, obj.anglefft, obj.fourier_ord(1),1, obj.time);
                obj.angleaccel = obj.eval_dfourier(obj.freqwrap, obj.anglefft, obj.fourier_ord(1),2, obj.time);

                obj.load = obj.eval_dfourier(obj.freqwrap, obj.loadfft,  obj.fourier_ord(1),0, obj.time);
                obj.loadvel = obj.eval_dfourier(obj.freqwrap, obj.loadfft, obj.fourier_ord(1),1, obj.time);
                obj.loadaccel = obj.eval_dfourier(obj.freqwrap, obj.loadfft, obj.fourier_ord(1),2, obj.time);

            end

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
            plot(obj.time,obj.angle_unfiltered);
            hold on
            plot(obj.time,obj.angle);
            hold off
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
        function [fssb,X, Mag, Phase,Pow,CumPow, fwrap] = FFT(obj, signal)
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
            
            % frequency vectors
            fssb = (0:N/2-1)*(Fs_up/N);
            fwrap = [(0:Fs_up/N:Fs_up/2) (-Fs_up/2+Fs_up/N:Fs_up/N:-Fs_up/N)];

            % FFT
            X = fft(signal_upsampled,N);
            SSB = X(1:N/2);
            SSB(2:end) = 2*SSB(2:end);

            Mag = abs(SSB/N);
            Phase = angle(SSB);
            Pow =Mag.^2;
            CumPow = cumsum(Pow);

        end
        
        % Note: ifft not needed 
        % function y = eval_fourier_ifft(obj,Y,order)  
        %     % this function returns the truncated Fourier series up to
        %     % order 'order', of the fft output vector Y, and donwsamples it to the
        %     % samplerate of the original signal
        %     Nfft = length(Y);
        %     bins = order+1; % amount of bins to retain
        %     Y_trunc = Y; 
        %     Y_trunc(bins+1:end-bins+1) = 0; % remove unwanted bins
        % 
        %     y = interp1(1:Nfft, ifft(Y_trunc,"symmetric"), linspace(1,Nfft,obj.Npoints),"linear");

        % end
    end % Methods

     methods (Static)
        function df = first_derivative_5pt(f, delta)
            % Computes the first derivative using the fourth-order (five-point) finite difference.
            % Input:
            %   f     - vector of function values [f1, f2, ..., fn]
            %   delta - uniform spacing between points
            % Output:
            %   df    - vector of first derivative approximations

            assert(length(f) >=5 ,'Less than 5 points for finite difference derivatives!')
            
            n = length(f);
            df = zeros(size(f)); % preallocate
            
            % Left boundary: f'_1
            df(1) = (1/(12*delta)) * (-25*f(1) + 48*f(2) - 36*f(3) + 16*f(4) - 3*f(5));
            
            % Near-left boundary: f'_2
            df(2) = (1/(12*delta)) * (-3*f(1) - 10*f(2) + 18*f(3) - 6*f(4) + f(5));
            
            % Interior points: f'_k for k = 3 to n-2
            for k = 3:n-2
                df(k) = (1/(12*delta)) * (-f(k+2) + 8*f(k+1) - 8*f(k-1) + f(k-2));
            end
            
            % Near-right boundary: f'_{n-1}
            df(n-1) = (1/(12*delta)) * (-f(n-4) + 6*f(n-3) - 18*f(n-2) + 10*f(n-1) + 3*f(n));
            
            % Right boundary: f'_n
            df(n) = (1/(12*delta)) * (3*f(n-4) - 16*f(n-3) + 36*f(n-2) - 48*f(n-1) + 25*f(n));
            
        end

        function d2f = second_derivative_5pt(f, delta)
            % Computes the second derivative using the fourth-order (five-point) finite difference.
            % Input:
            %   f     - vector of function values [f1, f2, ..., fn]
            %   delta - uniform spacing between points
            % Output:
            %   d2f   - vector of second derivative approximations

            assert(length(f) >=5 ,'Less than 5 points for finite difference derivatives!')
            
            n = length(f);
            d2f = zeros(size(f)); % preallocate
            
            % Left boundary: f''_1
            d2f(1) = (1/(12*delta^2)) * ( ...
                45*f(1) - 154*f(2) + 214*f(3) - 156*f(4) + 61*f(5) - 10*f(6) );
            
            % Near-left boundary: f''_2
            d2f(2) = (1/(12*delta^2)) * ( ...
                10*f(1) - 15*f(2) - 4*f(3) + 14*f(4) - 6*f(5) + f(6) );
            
            % Interior points: f''_k for k = 3 to n-2
            for k = 3:n-2
                d2f(k) = (1/(12*delta^2)) * ( ...
                    -f(k-2) + 16*f(k-1) - 30*f(k) + 16*f(k+1) - f(k+2) );
            end
            
            % Near-right boundary: f''_{n-1}
            d2f(n-1) = (1/(12*delta^2)) * ( ...
                f(n-5) - 6*f(n-4) + 14*f(n-3) - 4*f(n-2) - 15*f(n-1) + 10*f(n) );
            
            % Right boundary: f''_n
            d2f(n) = (1/(12*delta^2)) * ( ...
                -10*f(n-5) + 61*f(n-4) - 156*f(n-3) + 214*f(n-2) - 154*f(n-1) + 45*f(n) );
            
        end


        function f_val = eval_dfourier(fwrap, fft_vals, k, m, t_eval)
        % EVAL_DFOURIER Evaluates a derivative of the the Fourier series at a specific time
        %   f_val = evaluate_fourier_series(fft_vals, T, t_eval, k, m)
        %
        %   fwrap    : frequency bin center values corresponding with fft_vals
        %   fft_vals : The FFT output of the original signal (complex values)        
        %   k        : The number of harmonics (positive and negative) to include
        %   m        : Order of the derivative, 0 means no derivative
        %   t_eval   : The time (or array of times) at which to evaluate the series
        
            N = length(fft_vals);          % Number of FFT points
            n = -k:k;                      % Harmonic indices
            t_eval = t_eval(:)';           % Ensure row vector for broadcasting
            fwrap = fwrap(:)';
        
            % FFT indices wrap around from 0 to N-1 (1-based index for MATLAB)
            % Negative frequencies are in fft_vals(N-k+1:N), positive in fft_vals(2:k+1)
            X_k = zeros(size(n));
            fwrap_k = X_k;
        
            for k = 1:length(n)
                idx = mod(n(k), N) + 1;    % Wrap index to 1-based MATLAB indexing
                fwrap_k(k) = fwrap(idx);
                X_k(k) = fft_vals(idx) / N * (1i * 2 * pi * fwrap_k(k))^m; %get fft val and modify for m-th derivative
            end
        
            % Evaluate the Fourier series
            exponents = exp(1i * 2 * pi * fwrap_k' * t_eval);
            f_val = real(X_k * exponents);  % Take real part, as signal is real

        end
    

     end % methods (Static)
end

