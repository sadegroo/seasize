classdef geareval < handle
    %GEAREVAL evaluation of a motor vs. motion-load profile
    % to find a feasible reduction ratio, and simulate a series spring
    
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
        NKfigPstdContour

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
        
        function [result, fig] = eval(obj,Nmax, nN, varargin)
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
            Nfeas1 = N(cond1); %rms torque conditions
            Nfeas2 = N(cond2); %peak torque conditions
            Nfeas = N(cond1&cond2&cond3); %all conditions
            %Nfeas = N(cond1&cond3); %ignore peak condition
            %Nfeasrange = [min(Nfeas),min(Nfeas)/f2; max(Nfeas),max(Nfeas)/f2];
            if isempty(Nfeas)
                warning("No feasible reduction range found, evaluation terminated for this combination: motor="...
                    + string(obj.mot.name) + ", profile=" + string(obj.profile.description)...
                    + ", efficiency=" + string(obj.efficiency) + ", assistfactor=" + obj.assistfactor);
                return;
            end
            obj.results.minNTrms = min(Nfeas1)/f2;
            obj.results.minNTmax = min(Nfeas2)/f2;
            obj.results.minN = min(Nfeas)/f2;
            obj.results.maxN = max(Nfeas)/f2;
            obj.results.optN = min(max(Nfeas),N(idmin))/f2;
            
            % figure
            % -------
            obj.fig = figure();
            % rms curve
            plot(omeganorm,momentnormRms,'b:',omeganorm,cummin(momentnormRms),'b',wmot,Tmotrms,'.b', MarkerSize=10)
            hold on
            % max curve
            plot(omeganorm,momentnormMax,'r:',omeganorm,cummin(momentnormMax),'r',wmot,Tmotmax,'.r', MarkerSize=10)
            try
                % get intersects
                idxminmrms = obj.zci(momentnormRms-Tmotrms);
                idxmax = obj.zci(omeganorm-wmot);
                its1 = [omeganorm(idxminmrms(1)) Tmotrms];
                its2 = [wmot, momentnormRms(idxmax(1))];
                % plot intersects
                plot([its1(1) its2(1)], [its1(2) its2(2)], 'xb', MarkerSize=10);
                line([its1(1) wmot its2(1)],[its1(2), Tmotrms, its2(2)],'LineStyle','--');
            catch ME
                warning("RMS curve intersects could not be found: motor="...
                    + string(obj.mot.name) + ", profile=" + string(obj.profile.description)...
                    + ", efficiency=" + string(obj.efficiency) + ", assistfactor=" + obj.assistfactor);
            end
            try
                % get intersects
                idxminpk = obj.zci(momentnormMax-Tmotmax);
                its3 = [omeganorm(idxminpk(1)) Tmotmax];
                its4 = [wmot, momentnormMax(idxmax(1))];
                % plot intersects max curve
                plot([its3(1) its4(1)], [its3(2) its4(2)], 'xr', MarkerSize=10);
                line([its3(1) wmot its4(1)],[its3(2), Tmotmax its4(2)],'LineStyle','--', 'Color', 'r');
            catch ME
                warning("MAX curve intersects could not be found: motor="...
                    + string(obj.mot.name) + ", profile=" + string(obj.profile.description)...
                    + ", efficiency=" + string(obj.efficiency) + ", assistfactor=" + obj.assistfactor);
            end
            
            hold off
            legend({'$\omega^*-\tau_2^*$','MLB $\omega^*-\tau_2^*$', '($\omega_{mot,max}^*,\tau_{mot,2}^*$)','$\omega^*-\tau_\infty^*$','MLB $\omega^*-\tau_\infty^*$', '($\omega_{mot,max}^*,\tau_{mot,\infty}^*$)'},Interpreter="latex");
            %xlim([0 3*wmot]);
            xlabel('$\omega^*\,\left (\mathrm{\sqrt{kg\,m^2}\,rad/s}\right )$', "Interpreter","latex", "FontSize",12);
            ylim([0 2*Tmotmax]);
            ylabel('$\tau^*\,\left (\mathrm{\frac{Nm}{\sqrt{kg\,m^2}}}\right )$', "Interpreter","latex", "FontSize",12);
            grid on
            
            title("Motor " +obj.mot.name+ " vs. motion-load profile " +obj.profile.description);
            %set(gcf,'Visible','on')

            if obj.issea
                % additional SEA evaluations. Target metrics are RMS torque
                % reduction and standard deviation of power reduction

                Nrange_eval = linspace(obj.results.minN,obj.results.maxN,nNK(1));
                Krange_eval = linspace(Krange(1), Krange(2), nNK(2));

                mspring = zeros(nNK(1),nNK(2),obj.profile.Npoints);
                omegaspring = zeros(nNK(1),nNK(2),obj.profile.Npoints);
                pspring = zeros(nNK(1),nNK(2),obj.profile.Npoints);
                
                mspring_rms = zeros(nNK(1),nNK(2));
                mspring_rms_rel = zeros(nNK(1),nNK(2));
                pspring_peak = zeros(nNK(1),nNK(2));
                pspring_peak_rel = zeros(nNK(1),nNK(2));
                pspring_std = zeros(nNK(1),nNK(2));
                pspring_std_rel = zeros(nNK(1),nNK(2));

                mnospring = zeros(nNK(1),obj.profile.Npoints);
                pnospring = zeros(nNK(1),obj.profile.Npoints);

                idxminmrms = zeros(nNK(1),1);
                idxminpmax = zeros(nNK(1),1);
                idxminpstd = zeros(nNK(1),1);                

                for n = 1:nNK(1)
                    Neval = Nrange_eval(n);

                    mnospring_temp = obj.profile.load*obj.assistfactor + obj.mot.inertia*Neval^2*(obj.profile.angleaccel);                                %torque no spring
                    mnospring(n,:) = mnospring_temp';
                    pnospring(n,:) = (mnospring_temp.*(obj.profile.anglevel))';                                                                  %power no spring
                    
                    % no spring metrics as basis for relative "with spring" metrics
                    rms_mnospring = rms(mnospring_temp);
                    peak_pnospring = max(abs(pnospring(n,:)));
                    std_pnospring = std(pnospring(n,:));
                    

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
                        pspring_std(n,k) = std(pspring(n,k,:));  % standard deviation of power 
                        pspring_std_rel(n,k) = pspring_std(n,k)/std_pnospring;                          % relative std power                        

                                                
                    end
                    [~, idxminmrms(n)] = min(mspring_rms_rel(n,:));  % find K which results in minimum rms torque for this reduction
                    [~, idxminpmax(n)] = min(pspring_peak_rel(n,:));  % find K which results in minimum peak power for this reduction
                    [~, idxminpstd(n)] = min(pspring_std_rel(n,:));  % find K which results in minimum std of power for this reduction
                    
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
                obj.results.pspring_std = pspring_std;
                obj.results.pspring_std_rel = pspring_std_rel;

                obj.results.minpmaxcurve = Krange_eval(idxminpmax(:))'; % K value which results in minimum peak power for each value in Nrange_eval.
                obj.results.minpstdcurve = Krange_eval(idxminpstd(:))'; % K value which results in minimum power standard deviation for each value in Nrange_eval.
                obj.results.minmrmscurve = Krange_eval(idxminmrms(:))'; % K value which results in minimum rms torque for each value in Nrange_eval.
                
                % figures
                %--------

                % common meshgrid
                [X,Y] = meshgrid(Krange_eval, Nrange_eval);

                obj.NKfigMrmsContour = figure;
                [C1,h1] = contourf(X,Y,mspring_rms_rel,0.7:0.1:1.1);
                clabel(C1,h1)
                ylabel('Reduction Ratio N', "Interpreter","none");                
                xlabel('Spring rate K (Nm/rad)', "Interpreter","none");
                hold on
                plot(Krange_eval(idxminmrms(:)),Nrange_eval, 'r',LineWidth=1);
                plot(Krange_eval(idxminpstd(:)),Nrange_eval, 'g',LineWidth=1);
                hold off
                title('relative rms torque');
                colorbar;
                
                obj.NKfigPmaxContour = figure;
                [C2,h2] = contourf(X,Y,pspring_peak_rel, 0:0.1:1.5);
                clabel(C2,h2)
                ylabel('Reduction Ratio N', "Interpreter","none");                
                xlabel('Spring rate K (Nm/rad)', "Interpreter","none");
                hold on
                plot(Krange_eval(idxminmrms(:)),Nrange_eval, 'r',LineWidth=1);
                plot(Krange_eval(idxminpstd(:)),Nrange_eval, 'g',LineWidth=1);
                hold off
                title('relative peak absolute power');
                colorbar;

                obj.NKfigPstdContour = figure;
                [C3,h3] = contourf(X,Y,pspring_std_rel, 0:0.1:1.1);
                clabel(C3,h3)
                ylabel('Reduction Ratio N', "Interpreter","none");                
                xlabel('Spring rate K (Nm/rad)', "Interpreter","none");
                hold on
                plot(Krange_eval(idxminmrms(:)),Nrange_eval, 'r', LineWidth=1);
                plot(Krange_eval(idxminpstd(:)),Nrange_eval, 'g',LineWidth=1);
                hold off
                title('relative power standard deviation');
                colorbar;


            end
            
            result = obj.results;
            fig(1) = obj.fig;
            fig(2) = obj.NKfigMrmsContour;
            fig(3) = obj.NKfigPmaxContour;
            fig(4) = obj.NKfigPstdContour;
        end
            
        function [results, fig] = margins(obj,Nsel,varargin)
            % Calculates how much margin the motor has in terms of nominal
            % speed, nominal torque (overloaded) and maximum torque
            % (usually controller limited). A plot to determine if
            % the profile remains within the operating range of the motor
            % is returned in fig. A check to see if the profile crosses a
            % motor limit must be performed
            % results.limitcheck: 1x3 logical, true = limit breached
            % results.limitcheck(1): peak torque limit
            % results.limitcheck(2): no load speed limit
            % results.limitcheck(3): torque-speed characteristic limit

            arguments (Input)
                obj geareval
                Nsel (1,1) double
            end
            arguments (Repeating)
                varargin
            end

            obj=obj(:)';          

            p = inputParser;

            check =  @(x) isa(x,'double')&&isequal(size(x),[1 1]) ;
            addParameter(p, 'K', [], check);                     % A range of K to be evaluated for each evaluated gearbox ratio
                        
            parse(p,varargin{:});

            Ksel = p.Results.K;

            for j = 1:numel(obj)

            m_temp = obj(j).assistfactor*obj(j).profile.load ./ (Nsel*obj(j).efficiency) + obj(j).mot.inertia * Nsel * obj(j).profile.angleaccel;
            omega_temp = obj(j).profile.anglevel * Nsel;

            m_temp_nospring = m_temp;
            omega_temp_nospring = omega_temp;

            if obj(j).issea && ~isempty(Ksel)
                m_temp = (obj(j).profile.load*obj(j).assistfactor/obj(j).efficiency + obj(j).mot.inertia*Nsel^2*(obj(j).profile.angleaccel+obj(j).profile.loadaccel*obj(j).assistfactor/Ksel))/Nsel;
                omega_temp = (obj(j).profile.anglevel + obj(j).profile.loadvel*obj(j).assistfactor/Ksel)*Nsel;
            end
            
            % motor limit vals
            motpktor = obj(j).mot.peaktorque; % in Nm
            motnomtor = obj(j).mot.nomtorque;  % in Nm
            motstalltor = obj(j).mot.stalltorque; % in Nm
            motnomspd = obj(j).mot.nomspeed;    % in rad/s
            motnlspd = obj(j).mot.noloadspeed;  % in rad/s

           % motor speed-load char
           spdtorcharpos = @(x)  max((-motnomtor/(motnlspd-motnomspd))*(x-motnlspd), ...
               (motnomtor-motstalltor)/(motnomspd)*x + motstalltor);
           spdtorcharneg = @(x) -spdtorcharpos(-x);

            % check motor limits here: true= limit breached
           results(j).limitcheck(1) = any(m_temp>motpktor | m_temp<-motpktor);    % peaktorque
           results(j).limitcheck(2) = any(omega_temp>motnlspd | omega_temp<-motnlspd);         % max speed
           results(j).limitcheck(3) = any(m_temp > spdtorcharpos(omega_temp) | m_temp < spdtorcharneg(omega_temp)); %speed-torque char
           
            omega_load_max_atmot = max(abs(omega_temp));
            moment_load_rms_atmot =  rms(m_temp);
            moment_load_max_atmot = max(abs(m_temp));
           
            results(j).speed_margin = (obj(j).mot.nomspeed - omega_load_max_atmot)/obj(j).mot.nomspeed;
            results(j).moment_margin_rms = (obj(j).mot.rmstorque - moment_load_rms_atmot)/obj(j).mot.rmstorque;
            results(j).moment_margin_max = (obj(j).mot.peaktorque - moment_load_max_atmot)/obj(j).mot.peaktorque;

            fig(j) = figure;

            vec_w= 0:0.1:motnlspd;
            vec_T= spdtorcharpos(vec_w);
            
            if obj(j).issea && ~isempty(Ksel)
                % Verify that profiles are withing motor physical limits due to back emf
                plot(omega_temp,m_temp)
                hold on
                plot(omega_temp_nospring,m_temp_nospring)
                plot(vec_w,vec_T,'k--',-1*vec_w, -1*vec_T, 'k--')
                yline(obj(j).mot.peaktorque,'k--')
                yline(-obj(j).mot.peaktorque,'k--')
                xline(obj(j).mot.noloadspeed,'k--')
                xline(-obj(j).mot.noloadspeed,'k--')
                xlim([-obj(j).mot.noloadspeed*1.1, obj(j).mot.noloadspeed*1.1]);
                ylim([obj(j).mot.peaktorque*(-1.1) obj(j).mot.peaktorque*(1.1)]);
                legend({'With spring', 'No spring', 'Motor limits'}, "Location","best")
                xlabel('$\omega_{mot}$ (rad/s)', Interpreter='latex')
                ylabel('$\tau_{mot}$ (Nm)', Interpreter='latex')
                title("Motor physical limit check, profile=" + obj(j).profile.description + ", N=" + string(Nsel) + ", K=" +string(Ksel) + " Nm/rad")
            else
                 % Verify that profiles are withing motor physical limits due to back emf
                plot(omega_temp,m_temp)
                hold on
                plot(vec_w,vec_T,'k--',-1*vec_w, -1*vec_T, 'k--')
                yline(obj(j).mot.peaktorque,'k--')
                yline(-obj(j).mot.peaktorque,'k--')
                xline(obj(j).mot.noloadspeed,'k--')
                xline(-obj(j).mot.noloadspeed,'k--')
                xlim([-obj(j).mot.noloadspeed*1.1, obj(j).mot.noloadspeed*1.1]);
                ylim([obj(j).mot.peaktorque*(-1.1) obj(j).mot.peaktorque*(1.1)]);
                legend({'No spring', 'motor limits'}, "Location","best")
                xlabel('Angular velocity at motor shaft (rad/s)')
                ylabel('Torque at motor shaft (Nm)')
                title("Motor physical limit check, profile=" + obj(j).profile.description + ", N=" + string(Nsel))
            end
            hold off
            end

        end

        function [r, figs] = plotSelected(obj, varargin)
            arguments (Input)
                obj (1,1) geareval 
            end
            arguments (Repeating)
                varargin
            end

            p = inputParser;

            check =  @(x) isa(x,'double')&&(isequal(size(x),[1 1]) || isequal(size(x),[1 3])) ;
            addParameter(p, 'N', [1 100 100], check);                        % A range of K to be evaluated for each evaluated gearbox ratio
            addParameter(p, 'K', [10 1000 100], check);                                          % size of N-K space feasibility space

            % if K/N is a scalar, the exact value is selected. if it is a
            % 1x3 vector, the first 2 values are min/max limits and the
            % third is the amount of points
            
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

            if isscalar(Ksel) && isscalar(Nsel)
                % K and N provided
                % profiles without spring
                r.mnospring = obj.profile.load*obj.assistfactor/Nsel + obj.mot.inertia*Nsel*(obj.profile.angleaccel);
                r.omeganospring = Nsel*obj.profile.anglevel;
                r.accelnospring = Nsel*obj.profile.angleaccel;
                r.pnospring = r.mnospring.* r.omeganospring;

                % Altered gait profiles due to spring
                r.omegaspring = Nsel*(obj.profile.anglevel + obj.profile.loadvel*obj.assistfactor/Ksel);
                r.accelspring = Nsel*(obj.profile.angleaccel+obj.profile.loadaccel*obj.assistfactor/Ksel);
                r.mspring = obj.profile.load*obj.assistfactor/Nsel + obj.mot.inertia*r.accelspring;

                r.pspring = r.mspring.*r.omegaspring;

                % rms torque
                r.mnospring_rms = rms(r.mnospring);
                r.mspring_rms = rms(r.mspring);
                r.mspring_rms_rel = r.mspring_rms/r.mnospring_rms;

                % peak power
                r.pnospring_max = max(abs(r.pnospring));
                r.pspring_max = max(abs(r.pspring));
                r.pspring_max_rel =  r.pspring_max/r.pnospring_max;

                % std power
                r.pnospring_std = std(r.pnospring);
                r.pspring_std = std(r.pspring);
                r.pspring_std_rel =  std(r.pspring_std);                
                
                % energy consumption
                % time = obj.profile.time;                
                % r.enospring = trapz(time, r.pnospring);
                % r.espring = trapz(time, r.pspring);
                % r.espring_rel =  r.espring/r.enospring;

                
                % Motor power plot, with and without spring
                figs(1) = figure;
                plot(obj.profile.time, r.pspring)
                hold on
                plot(obj.profile.time, r.pnospring)
                legend({'motor power with spring','motor power without spring'},"Location","best");
                xlabel('Time (s)')
                ylabel('Power (W)')
                title("Motor power, profile = " + obj.profile.description +", N="  + string(Nsel) +", K=" +string(Ksel) + " Nm/rad")
                hold off
                
                % Motor torque plot, with and without spring
                figs(2) = figure;
                plot(obj.profile.time, r.mspring);
                hold on
                plot(obj.profile.time, r.mnospring)
                legend({'motor torque with spring','motor torque without spring'},"Location","best");
                xlabel('Time (s)')
                ylabel('Torque (Nm)')
                title("Motor torque, profile = " + obj.profile.description +", N="  + string(Nsel) +", K=" +string(Ksel) + " Nm/rad")
                hold off
                
                % Motor velocity plot, with and without spring
                figs(3) = figure;
                plot(obj.profile.time, r.omegaspring)
                hold on
                plot(obj.profile.time, r.omeganospring)
                legend({'motor speed with spring','motor speed without spring'},"Location","best");
                xlabel('Time (s)')
                ylabel('Velocity (rad/s)')
                title("Motor velocity, profile = " + obj.profile.description +", N="  + string(Nsel) +", K=" +string(Ksel) + " Nm/rad")
                hold off

                % Motor velocity plot, with and without spring
                figs(4) = figure;
                plot(obj.profile.time, r.accelspring)
                hold on
                plot(obj.profile.time, r.accelnospring)
                legend({'motor acceleration with spring','motor acceleration without spring'},"Location","best");
                xlabel('Time (s)')
                ylabel('Acceleration (rad/s²)')
                title("Motor acceleration, profile = " + obj.profile.description +", N="  + string(Nsel) +", K=" +string(Ksel) + " Nm/rad")
                hold off

            elseif isscalar(Nsel)
                % only N provided

                % profiles without spring
                r.mnospring = obj.profile.load*obj.assistfactor/Nsel + obj.mot.inertia*Nsel*(obj.profile.angleaccel);
                r.omeganospring = Nsel*obj.profile.anglevel;
                r.accelnospring = Nsel*obj.profile.angleaccel;
                r.pnospring = r.mnospring.* r.omeganospring;

               % rms torque
                r.mnospring_rms = rms(r.mnospring);
                % peak power
                r.pnospring_max = max(abs(r.pnospring));
                % std power
                r.pnospring_std = std(r.pnospring);
                % energy consumption
                % time = obj.profile.time;                
                % r.enospring = trapz(time, r.pnospring);

                Kspace = linspace(Ksel(1), Ksel(2), Ksel(3));

                for i=1:numel(Kspace)
                    Kloop = Kspace(i);
    
                    % Altered gait profiles due to spring
                    omegaspring = Nsel*(obj.profile.anglevel + obj.profile.loadvel*obj.assistfactor/Kloop);
                    accelspring = Nsel*(obj.profile.angleaccel+obj.profile.loadaccel*obj.assistfactor/Kloop);
                    mspring = obj.profile.load*obj.assistfactor/Nsel + obj.mot.inertia*accelspring;
                    pspring = mspring.*omegaspring;
    
                    % rms torque
                    r.mspring_rms(i) = rms(mspring);
                    r.mspring_rms_rel(i) = r.mspring_rms(i)/r.mnospring_rms;
    
                    % peak power
                    r.pspring_max(i) = max(abs(pspring));
                    r.pspring_max_rel(i) =  r.pspring_max(i)/r.pnospring_max;
    
                    % std power
                    r.pspring_std(i) =std(pspring);
                    r.pspring_std_rel(i) =  r.pspring_std(i)/r.pnospring_std;                
                    
                    % energy consumption
                    % time = obj.profile.time;                
                    % r.espring(i) = trapz(time, pspring);
                    % r.espring_rel(i) =  r.espring(i)/r.enospring;
              
                end

                % relative metrics figure, varied K
                figs(1) = figure;
                lgd = cell(0);
                hold on
                plot(Kspace, r.mspring_rms_rel, 'r')   
                lgd = horzcat(lgd, {'r_{\tau,rms}'});
                xlabel('K (Nm/rad)')
                ylabel('r')
                title("Relative metrics, N="  + string(Nsel))
                
                Kzci = obj.zci(r.mspring_rms_rel-1);
                try
                    r.K_be_rms_tor = Kspace(Kzci(1));
                    [min_rms_tor,idxmin] = min(r.mspring_rms_rel);
                    r.K_min_rms_tor = Kspace(idxmin);
                    % plot break-even and min markers
                    plot(r.K_be_rms_tor,1,'xr','MarkerSize',10, LineWidth=1.5);
                    plot(r.K_min_rms_tor,min_rms_tor,'or','MarkerSize',10,  LineWidth=1.5);
                    lgd = horzcat(lgd, {'K_{eq,\tau, rms}', 'K_{min,\tau, rms}'});
                catch ME
                     r.K_be_rms_tor = NaN;
                     r.K_min_rms_tor = NaN;

                end

                xlim([Ksel(1)/2 Ksel(2)])
                
                % peak power, varied K
                plot(Kspace, r.pspring_max_rel,'b')
                lgd = horzcat(lgd, {'r_{P,max}'});
                Kzci = obj.zci(r.pspring_max_rel-1);
                try
                    r.K_be_pk_p = Kspace(Kzci(1));
                    [min_pk_p,idxmin] = min(r.pspring_max_rel);
                    r.K_min_pk_p = Kspace(idxmin);
                    % plot break-even and min markers
                    plot(r.K_be_pk_p,1,'xb','MarkerSize',10, LineWidth=1.5);
                    plot(r.K_min_pk_p,min_pk_p,'ob','MarkerSize',10,  LineWidth=1.5);
                    lgd = horzcat(lgd, {'K_{eq,P, max}', 'K_{min,P, max}'});
                catch ME
                    r.K_be_pk_p = NaN;
                    r.K_min_pk_p = NaN;
                end
                
                % std of power, varied K
                plot(Kspace, r.pspring_std_rel,'g')  
                lgd = horzcat(lgd, {'r_{P,std}'});
                Kzci = obj.zci(r.pspring_std_rel-1);
                try
                    r.K_be_std_p = Kspace(Kzci(1));
                    [min_std_p,idxmin] = min(r.pspring_std_rel);
                    r.K_min_std_p = Kspace(idxmin);
                    % plot break-even and min markers
                    plot(r.K_be_std_p,1,'xg','MarkerSize',10, LineWidth=1.5);
                    plot(r.K_min_std_p,min_std_p,'og','MarkerSize',10,  LineWidth=1.5);
                    lgd = horzcat(lgd, {'K_{eq,P, std}', 'K_{min,P, std}'});
                catch ME
                    r.K_be_std_p = NaN;
                    r.K_min_std_p = NaN;
                end
                
                
                yline(1,'k--')
                lgd = horzcat(lgd, {'r = 1 (Break-even)'});
                hold off

                legend(lgd);


                % % energy consumption, varied K
                % figs(4) = figure;
                % plot(Kspace, r.espring_rel)
                % hold on
                % yline(1,'k')
                % hold off
                % xlabel('K (Nm/rad)')
                % ylabel('relative energy consumption')
                % title("Relative energy consumption, N="  + string(Nsel))
                % hold off

                % Kzci = obj.zci(r.espring_rel-1);
                % try
                %     r.K_be_e = Kspace(Kzci(1));
                % catch ME
                %     r.K_be_e = NaN;
                % end
                % 
                % [~,idxmin] = min(r.espring_rel);
                % r.K_min_e = Kspace(idxmin);
                                
            elseif isscalar(Ksel)
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

        function [result, fig] = feasableBoxplot(obj, varargin)
            arguments (Input)
                obj geareval
            end
            arguments (Repeating)
                varargin
            end

            obj=obj(:)';    % make row vector of objects

            for i = 1:numel(obj)
                try 
                    mins(i) = obj(i).results.minN;
                    maxs(i) = obj(i).results.maxN;
                catch ME
                    mins(i) = NaN;
                    maxs(i) = NaN;
                end
            end
            
            [result.maxmin, result.maxminidx] = max(mins);
            [result.minmax, result.minmaxidx] = min(maxs);

            fig=figure;
            boxplot([mins;maxs],varargin)
            xtickangle(90)
            set(gca,'XGrid','off','YGrid','on');
            ylim([min(mins)-10, max(maxs)+10])
            title("Feasible Reduction")
            ylabel('Reduction ratio N')
            %fig.Visible = "on";
        end

        function [fig] = feasableBoxplotBounds(obj)
            arguments (Input)
                obj geareval
            end

            obj=obj(:)';    % make row vector of objects

            for i = 1:numel(obj)
                try 
                    mins(i) = obj(i).results.minN;
                    maxs(i) = obj(i).results.maxN;
                    optTrms(i) = obj(i).results.optN;
                    center(i) = mean([mins(i) maxs(i)]);
                catch ME
                    mins(i)= NaN;
                    maxs(i) = NaN;
                    optTrms(i) = NaN;
                    center(i) = NaN;
                end
            end
            
            fig=figure;
            boxplot([mins' center' maxs'],{'Lower bound', 'Center', 'Upper bound'})
            xtickangle(90)
            set(gca,'XGrid','off','YGrid','on');
            ylim([min(mins)-10, max(maxs)+10])
            title("Feasible Reduction")
            ylabel('Reduction ratio N')
            %fig.Visible = "on";
        end

        function [fig] = feasableStaircase(obj,Npoints)
            arguments (Input)
                obj geareval
                Npoints (1,1) {mustBeInteger, mustBePositive, mustBeScalarOrEmpty, mustBeNonempty}
            end

            obj=obj(:)';    % make row vector of objects

            for i = 1:numel(obj)
                try 
                    mins(i) = obj(i).results.minN;
                    maxs(i) = obj(i).results.maxN;
                    optTrms(i) = obj(i).results.optN;
                    center(i) = mean([mins(i) maxs(i)]);
                catch ME
                    mins(i)= NaN;
                    maxs(i) = NaN;
                    optTrms(i) = NaN;
                    center(i) = NaN;
                end
            end
            
            fig=figure;
            startN = floor(min(mins));
            stopN = ceil(max(maxs));
            Nrange = linspace(startN,stopN,Npoints);
            for i = 1:Npoints
                Nthis = Nrange(i);
                feascount(i) = sum(Nthis > mins & Nthis<maxs,"all");
            end
            plot(Nrange,feascount)
            xlabel('Reduction ratio N')
            ylabel('Feasible profiles')
            ylim([0, numel(obj)]);
        end

    end

    methods (Static)
    end
end

