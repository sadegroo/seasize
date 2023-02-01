classdef geareval < handle
    %GEAREVAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mot
        profile
        efficiency
        assistfactor
        issea
        results
        fig
        NKfigMrmsContour
        NKfigPmaxContour
        NKfigPavgContour

    end

    properties (Constant)
        zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
    end
    
    methods
        function obj = geareval(mot, profile, efficiency, assistfactor, varargin)
            %GEAREVAL Construct an instance of this class
            %   Detailed explanation goes here
            arguments (Input)
                mot (1,1) motor
                profile (1,1) motionprofile
                efficiency (1,1) double
                assistfactor (1,1) double
            end
            arguments (Repeating)
                varargin
            end

             p = inputParser;

            checkLogicalScalar = @(x) islogical(x)&&(isequal(size(x), [1 1]));
            addParameter(p, 'sea', false, checkLogicalScalar);

            parse(p,varargin{:});

            obj.issea = p.Results.sea;

            obj.mot = mot;
            obj.profile = profile;
            obj.efficiency = efficiency;
            obj.assistfactor = assistfactor;
                        
        end
        
        function eval(obj,Nmax, nN, varargin)
            arguments (Input)
                obj geareval
                Nmax (1,1) double        % maximum transmission ratio to evaluate
                nN (1,1) double          % amount of points to generate between 1 and Nmax to find feasible gearbox range         
            end
            arguments (Repeating)
                varargin
            end

            p = inputParser;

            check =  @(x) isa(x,'double')&&isequal(size(x),[1 2]) ;
            addParameter(p, 'Krange', [], check);            % A range of K to be evaluated for each evaluated gearbox ratio
            addParameter(p, 'nNK', [], check);                                          % size of N-K space feasibility space
            
            parse(p,varargin{:});

            Krange = p.Results.Krange;
            nNK = p.Results.nNK;
            
            %normalisation factors
            f1 = sqrt(obj.mot.inertia/obj.efficiency);
            f2 = sqrt(obj.mot.inertia*obj.efficiency);
            wmot = obj.mot.nomspeed * f2;
            Tmotmax = obj.mot.peaktorque/f1;
            Tmotrms = obj.mot.rmstorque/f1;    
            
            N = linspace(1,Nmax, nN);
            N=N(:) * f2;             % make column vector and normalize
            omeganorm = zeros(nN,1);
            momentnormRms = zeros(nN,1);
            momentnormMax = zeros(nN,1);

            for i = 1:nN
                omeganorm(i) = max(abs(obj.profile.anglevel)).* N(i);
                momentnormRms(i) =  rms(obj.assistfactor*obj.profile.load ./ N(i) + N(i).*obj.profile.angleaccel);
                momentnormMax(i) = max(abs(obj.assistfactor*obj.profile.load./N(i) + N(i).*obj.profile.angleaccel));
            end
           
            cond1 = Tmotrms > momentnormRms;
            cond2 = Tmotmax > momentnormMax;
            cond3 = wmot > omeganorm;
            [~,idmin] = min(momentnormRms);
            %cond4 = [true(idmin,1);false(Npoints-idmin,1)];
            Nfeas = N(cond1&cond2&cond3); %all conditions
            %Nfeas = N(cond1&cond3); %ignore peak condition
            %Nfeasrange = [min(Nfeas),min(Nfeas)/f2; max(Nfeas),max(Nfeas)/f2];
            if isempty(Nfeas)
                disp("No feasible reduction range found, evaluation terminated for this combination: motor="...
                    + string(obj.mot.name) + ", profile=" + string(obj.profile.description)...
                    + ", efficiency=" + string(obj.efficiency) + ", assistfactor=" + obj.assistfactor);
                return;
            end
            obj.results.minN = min(Nfeas)/f2;
            obj.results.maxN = max(Nfeas)/f2;
            obj.results.optN = min(max(Nfeas),N(idmin))/f2;
            
            % figure
            % -------
            obj.fig = figure();
            % rms curve
            plot(omeganorm,momentnormRms,'b',omeganorm,cummin(momentnormRms),'b:',wmot,Tmotrms,'xb')
            hold on
            % get intersects
            idxminmrms = obj.zci(momentnormRms-Tmotrms);
            idxmax = obj.zci(omeganorm-wmot);
            its1 = [omeganorm(idxminmrms(1)) Tmotrms];
            its2 = [wmot, momentnormRms(idxmax(1))];
            % max curve
            plot(omeganorm,momentnormMax,'r',omeganorm,cummin(momentnormMax),'r:',wmot,Tmotmax,'xr')
            % get intersects
            idxminpk = obj.zci(momentnormMax-Tmotmax);
            its3 = [omeganorm(idxminpk(1)) Tmotmax];
            its4 = [wmot, momentnormMax(idxmax(1))];
            % plot intersects
            plot([its1(1) its2(1)], [its1(2) its2(2)], '.b');
            line([its1(1) wmot its2(1)],[its1(2), Tmotrms, its2(2)],'LineStyle','--');
            plot([its3(1) its4(1)], [its3(2) its4(2)], '.r');
            line([its3(1) wmot its4(1)],[its3(2), Tmotmax its4(2)],'LineStyle','--', 'Color', 'r');
            
            hold off
            legend('rms norm','MLB rms', ['rms ' obj.mot.name],'max norm','MLB max', ['max ' obj.mot.name]);
            xlim([0 2*wmot]);
            xlabel('$\omega^*\,[\sqrt{kg\,m^2}\,rad/s]$', "Interpreter","latex", "FontSize",14);
            ylim([0 2*Tmotmax]);
            ylabel('$\tau^*\,[\frac{Nm}{\sqrt{kg\,m^2}}]$', "Interpreter","latex", "FontSize",14);
            
            title(obj.profile.description);
            set(gcf,'Visible','on')

            if obj.issea
                % additional SEA evaluations. Target metrics are RMS torque
                % reduction and average power reduction (= energy
                % consumption reduction)

                Nrange_eval = linspace(obj.results.minN,obj.results.maxN,nNK(1));
                Krange_eval = linspace(Krange(1), Krange(2), nNK(2));

                mspring = zeros(nNK(1),nNK(2),obj.profile.Npoints);
                omegaspring = zeros(nNK(1),nNK(2),obj.profile.Npoints);
                pspring = zeros(nNK(1),nNK(2),obj.profile.Npoints);
                
                mspring_rms = zeros(nNK(1),nNK(2));
                mspring_rms_rel = zeros(nNK(1),nNK(2));
                pspring_peak = zeros(nNK(1),nNK(2));
                pspring_peak_rel = zeros(nNK(1),nNK(2));
                pspring_avg = zeros(nNK(1),nNK(2));
                pspring_avg_rel = zeros(nNK(1),nNK(2));

                mnospring = zeros(nNK(1),obj.profile.Npoints);
                pnospring = zeros(nNK(1),obj.profile.Npoints);

                idxminmrms = zeros(nNK(1),1);
                idxminpmax = zeros(nNK(1),1);
                idxminpavg = zeros(nNK(1),1);                

                for n = 1:nNK(1)
                    Neval = Nrange_eval(n);

                    mnospring_temp = obj.profile.load*obj.assistfactor + obj.mot.inertia*Neval^2*(obj.profile.angleaccel);                                %torque no spring
                    mnospring(n,:) = mnospring_temp';
                    pnospring(n,:) = (mnospring_temp.*(obj.profile.anglevel))';                                                                  %power no spring
                    
                    % no spring metrics as basis for relative "with spring" metrics
                    rms_mnospring = rms(mnospring_temp);
                    peak_pnospring = max(abs(pnospring(n,:)));
                    avg_pnospring = trapz(obj.profile.time, pnospring(n,:))/obj.profile.period;
                    

                    for k = 1:nNK(2)
                        Keval = Krange_eval(k);

                        mspring_temp = (obj.profile.load*obj.assistfactor + obj.mot.inertia*Neval^2*(obj.profile.angleaccel+obj.profile.loadaccel*obj.assistfactor/Keval));
                        mspring(n,k,:) = mspring_temp';

                        omegaspring_temp = (obj.profile.anglevel + obj.profile.loadvel*obj.assistfactor/Keval);
                        omegaspring(n,k,:) = omegaspring_temp';
                        
                        pspring(n,k,:) = (mspring_temp .* omegaspring_temp)';       % power profile

                        mspring_rms(n,k) = rms(mspring(n,k,:));                                         % rms torque
                        mspring_rms_rel(n,k) =  mspring_rms(n,k)/rms_mnospring;                         % relative rms torque
                        pspring_peak(n,k) = max(abs(pspring(n,k,:)));                                   % peak power
                        pspring_peak_rel(n,k) = pspring_peak(n,k)/peak_pnospring;                       % relative peak power
                        pspring_avg(n,k) = trapz(obj.profile.time, pspring(n,k,:))/obj.profile.period;  % average power using trapezoidal integration
                        pspring_avg_rel(n,k) = pspring_avg(n,k)/avg_pnospring;                          % relative average power                        

                                                
                    end
                    [~, idxminmrms(n)] = min(mspring_rms_rel(n,:));  % find K which results in minimum rms torque for this reduction
                    [~, idxminpmax(n)] = min(pspring_peak_rel(n,:));  % find K which results in minimum peak power for this reduction
                    [~, idxminpavg(n)] = min(pspring_avg_rel(n,:));  % find K which results in minimum average power for this reduction
                    
                end

                % write to results output
                obj.results.Nrange_evel = Nrange_eval;
                obj.results.Krange_eval = Krange_eval;
                
                obj.results.mnospring = mnospring;
                obj.results.pnospring = pnospring;

                obj.results.mspring = mspring;
                obj.results.omegaspring = omegaspring;
                obj.results.pspring = pspring;

                obj.results.mspring_rms = mspring_rms;
                obj.results.mspring_rms_rel = mspring_rms_rel;
                obj.results.pspring_peak = pspring_peak;
                obj.results.pspring_peak_rel = pspring_peak_rel;
                obj.results.pspring_avg = pspring_avg;
                obj.results.pspring_avg_rel = pspring_avg_rel;

                obj.results.minpmaxcurve = Krange_eval(idxminpmax(:))'; % K value which results in minimum peak power for each value in Nrange_eval.
                obj.results.minpavgcurve = Krange_eval(idxminpavg(:))'; % K value which results in minimum average power for each value in Nrange_eval.
                obj.results.minmrmscurve = Krange_eval(idxminmrms(:))'; % K value which results in minimum rms torque for each value in Nrange_eval.
                
                % figures
                %--------

                % common meshgrid
                [X,Y] = meshgrid(Krange_eval, Nrange_eval);

                obj.NKfigMrmsContour = figure;
                [C1,h1] = contourf(X,Y,mspring_rms_rel,0.7:0.1:1.1);
                clabel(C1,h1)
                ylabel('Reduction Ratio', "Interpreter","none");                
                xlabel('Spring rate [Nm/rad]', "Interpreter","none");
                hold on
                plot(Krange_eval(idxminmrms(:)),Nrange_eval, 'r');
                hold off
                title('relative rms torque');
                
                obj.NKfigPmaxContour = figure;
                [C2,h2] = contourf(X,Y,pspring_peak_rel, 0:0.1:1.5);
                clabel(C2,h2)
                ylabel('Reduction Ratio', "Interpreter","none");                
                xlabel('Spring rate [Nm/rad]', "Interpreter","none");
                hold on
                plot(Krange_eval(idxminmrms(:)),Nrange_eval, 'r');
                plot(Krange_eval(idxminpmax(:)),Nrange_eval, 'g');
                hold off
                title('relative peak absolute power');

                obj.NKfigPavgContour = figure;
                [C3,h3] = contourf(X,Y,pspring_avg_rel, 0:0.1:1.1);
                clabel(C3,h3)
                ylabel('Reduction Ratio', "Interpreter","none");                
                xlabel('Spring rate [Nm/rad]', "Interpreter","none");
                hold on
                plot(Krange_eval(idxminmrms(:)),Nrange_eval, 'r');
                plot(Krange_eval(idxminpavg(:)),Nrange_eval, 'g');
                hold off
                title('relative average power');


            end

        end
            
        function [results, fig] = margins(obj,Nsel,varargin)
            % Calculates how much margin the motor has in terms of nominal
            % speed, nominal torque (overloaded) and maximum torque
            % (usually controller limited). A plot to determine if
            % the profile remains within the operating range of the motor
            % is returned in fig. A check to see if the profile crosses a
            % motor limit must be performed (not implemented
            % programmatically)
            arguments (Input)
                obj geareval
                Nsel (1,1) double
            end
            arguments (Repeating)
                varargin
            end

            p = inputParser;

            check =  @(x) isa(x,'double')&&isequal(size(x),[1 1]) ;
            addParameter(p, 'K', [], check);                     % A range of K to be evaluated for each evaluated gearbox ratio
                        
            parse(p,varargin{:});

            Ksel = p.Results.K;

            m_temp = obj.assistfactor*obj.profile.load ./ (Nsel*obj.efficiency) + obj.mot.inertia * Nsel * obj.profile.angleaccel;
            omega_temp = obj.profile.anglevel * Nsel;

            m_temp_nospring = m_temp;
            omega_temp_nospring = omega_temp;

            if obj.issea && ~isempty(Ksel)
                m_temp = (obj.profile.load*obj.assistfactor/obj.efficiency + obj.mot.inertia*Nsel^2*(obj.profile.angleaccel+obj.profile.loadaccel*obj.assistfactor/Ksel))/Nsel;
                omega_temp = (obj.profile.anglevel + obj.profile.loadvel*obj.assistfactor/Ksel)*Nsel;
            end

            omega_load_max_atmot = max(omega_temp);
            moment_load_rms_atmot =  rms(m_temp);
            moment_load_max_atmot = max(m_temp);
            
            results.speed_margin = (obj.mot.nomspeed - omega_load_max_atmot)/obj.mot.nomspeed;
            results.moment_margin_rms = (obj.mot.rmstorque - moment_load_rms_atmot)/obj.mot.rmstorque;
            results.moment_margin_max = (obj.mot.peaktorque - moment_load_max_atmot)/obj.mot.peaktorque;

            fig = figure;
            
            if obj.issea && ~isempty(Ksel)
                % Verify that profiles are withing motor physical limits due to back emf
                plot(omega_temp,m_temp)
                hold on
                plot(omega_temp_nospring,m_temp_nospring)
                vec_w= [0 obj.mot.nomspeed obj.mot.noloadspeed];
                vec_T= [obj.mot.stalltorque obj.mot.nomtorque 0];
                plot(vec_w,vec_T,'g',-1*vec_w, -1*vec_T, 'g')
                yline(obj.mot.peaktorque,'g')
                yline(-obj.mot.peaktorque,'g')
                xline(obj.mot.noloadspeed,'g')
                xline(-obj.mot.noloadspeed,'g')
                xlim([-obj.mot.noloadspeed*1.1, obj.mot.noloadspeed*1.1]);
                ylim([obj.mot.peaktorque*(-1.1) obj.mot.peaktorque*(1.1)]);
                legend({'motor torque - velocity profile with spring', 'motor torque - velocity profile without spring', 'motor limits'}, "Location","best")
                xlabel('Angular velocity at motor shaft (rad/s)')
                ylabel('Torque at motor shaft (Nm)')
                title("Motor physical limit check, profile=" + obj.profile.description + ", N=" + string(Nsel) + ", K=" +string(Ksel) + " Nm/rad")
            else
                 % Verify that profiles are withing motor physical limits due to back emf
                plot(omega_temp_nospring,m_temp_nospring)
                hold on
                vec_w= [0 obj.mot.nomspeed obj.mot.noloadspeed];
                vec_T= [obj.mot.stalltorque obj.mot.nomtorque 0];
                plot(vec_w,vec_T,'g',-1*vec_w, -1*vec_T, 'g')
                yline(obj.mot.peaktorque,'g')
                yline(-obj.mot.peaktorque,'g')
                xline(obj.mot.noloadspeed,'g')
                xline(-obj.mot.noloadspeed,'g')
                xlim([-obj.mot.noloadspeed*1.1, obj.mot.noloadspeed*1.1]);
                ylim([obj.mot.peaktorque*(-1.1) obj.mot.peaktorque*(1.1)]);
                legend({'motor torque - velocity profile without spring', 'motor limits'}, "Location","best")
                xlabel('Angular velocity at motor shaft (rad/s)')
                ylabel('Torque at motor shaft (Nm)')
                title("Motor physical limit check, profile=" + obj.profile.description + ", N=" + string(Nsel))
            end
            hold off

        end

        function [results_selected, figs] = plotSelected(obj, varargin)
            arguments (Input)
                obj geareval
            end
            arguments (Repeating)
                varargin
            end

            p = inputParser;

            check =  @(x) isa(x,'double')&&isequal(size(x),[1 1]) ;
            addParameter(p, 'N', [], check);                        % A range of K to be evaluated for each evaluated gearbox ratio
            addParameter(p, 'K', [], check);                                          % size of N-K space feasibility space
            
            parse(p,varargin{:});
            
            Nsel = p.Results.N;
            Ksel = p.Results.K;

            exception_sea = MException('geareval:plotSelected:InvalidKinput', ...
                                       'K can only be inputted when the issea property is true and eval() was performed');
            not_impl = MException('geareval:plotSelected:NotImplemented', ...
                                       'not implmented combination of inputs');
            if ~isempty(Ksel) && ~obj.issea
                throw(exception_sea);
            end               

            if ~isempty(Ksel) && ~isempty(Nsel)
                % K and N provided
                % profiles without spring
                results_selected.mnospring = (obj.profile.load*obj.assistfactor + obj.mot.inertia*Nsel^2*(obj.profile.angleaccel))/Nsel;
                results_selected.omeganospring = Nsel*obj.profile.anglevel;
                results_selected.accelnospring = Nsel*obj.profile.angleaccel;
                results_selected.pnospring = results_selected.mnospring.* results_selected.omeganospring;

                % Altered gait profiles due to spring
                results_selected.mspring = (obj.profile.load*obj.assistfactor + obj.mot.inertia*Nsel^2*(obj.profile.angleaccel+obj.profile.loadaccel*obj.assistfactor/Ksel))/Nsel;
                results_selected.omegaspring = Nsel*(obj.profile.anglevel + obj.profile.loadvel*obj.assistfactor/Ksel);
                results_selected.accelspring = Nsel*(obj.profile.angleaccel+obj.profile.loadaccel*obj.assistfactor/Ksel);
                results_selected.pspring = results_selected.mspring.*results_selected.omegaspring;

                % rms torque
                results_selected.mnospring_rms = rms(results_selected.mnospring);
                results_selected.mspring_rms = rms(results_selected.mspring);
                results_selected.mspring_rms_rel = results_selected.mspring_rms/results_selected.mnospring_rms;

                % peak power
                results_selected.pnospring_max = max(abs(results_selected.pnospring));
                results_selected.pspring_max = max(abs(results_selected.pspring));
                results_selected.pspring_max_rel =  results_selected.pspring_max/results_selected.pnospring_max;

                % mean power
                results_selected.pnospring_avg = trapz(obj.profile.time, results_selected.pnospring)/obj.profile.period;
                results_selected.pspring_avg = trapz(obj.profile.time, results_selected.pspring)/obj.profile.period;
                results_selected.pspring_avg_rel =  results_selected.pspring_avg/results_selected.pnospring_avg;                
                
                % energy consumption
                time = obj.profile.time;                
                results_selected.enospring = trapz(time, results_selected.pnospring);
                results_selected.espring = trapz(time, results_selected.pspring);
                results_selected.espring_rel =  results_selected.espring/results_selected.enospring;

                
                % Motor power plot, with and without spring
                figs(1) = figure;
                plot(obj.profile.time, results_selected.pspring)
                hold on
                plot(obj.profile.time, results_selected.pnospring)
                legend({'motor power with spring','motor power without spring'},"Location","best");
                xlabel('Time (s)')
                ylabel('Power (W)')
                title("Motor power, profile = " + obj.profile.description +", N="  + string(Nsel) +", K=" +string(Ksel) + " Nm/rad")
                hold off
                
                % Motor torque plot, with and without spring
                figs(2) = figure;
                plot(obj.profile.time, results_selected.mspring);
                hold on
                plot(obj.profile.time, results_selected.mnospring)
                legend({'motor torque with spring','motor torque without spring'},"Location","best");
                xlabel('Time (s)')
                ylabel('Torque (Nm)')
                title("Motor torque, profile = " + obj.profile.description +", N="  + string(Nsel) +", K=" +string(Ksel) + " Nm/rad")
                hold off
                
                % Motor velocity plot, with and without spring
                figs(3) = figure;
                plot(obj.profile.time, results_selected.omegaspring)
                hold on
                plot(obj.profile.time, results_selected.omeganospring)
                legend({'motor speed with spring','motor speed without spring'},"Location","best");
                xlabel('Time (s)')
                ylabel('Velocity (rad/s)')
                title("Motor velocity, profile = " + obj.profile.description +", N="  + string(Nsel) +", K=" +string(Ksel) + " Nm/rad")
                hold off

                % Motor velocity plot, with and without spring
                figs(4) = figure;
                plot(obj.profile.time, results_selected.accelspring)
                hold on
                plot(obj.profile.time, results_selected.accelnospring)
                legend({'motor acceleration with spring','motor acceleration without spring'},"Location","best");
                xlabel('Time (s)')
                ylabel('Acceleration (rad/s²)')
                title("Motor acceleration, profile = " + obj.profile.description +", N="  + string(Nsel) +", K=" +string(Ksel) + " Nm/rad")
                hold off

            elseif ~isempty(Nsel)
                % only N provided
                % not implemented yet
                throw(not_impl);
            elseif ~isempty(Ksel)
                % only K provided
                % not implemented yet
                throw(not_impl);
            else
                % invalid input
                exception = MException('geareval:plotSelected:InvalidInput', ...
                                       'K and/or N name value pair must be provided');
                % Throw the exception
                throw(exception)
            end

        end

    end

    methods (Static)
    end
end

