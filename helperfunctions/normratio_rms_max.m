function [N,omeganorm,momentnormRms,momentnormMax] = normratio_rms_max(JNT,Nrange,n,invertload,assistfactor)
    % This function returns, for a range of normalised transmission ratios [1],
    % the respective minimum normalised angular velocity [2] and minimum
    % normalised torque , evaluated by the 2-norm (RMS)[3] and infinity-norm (maximum)[4]
    % Inputs are: (1)angular velocity (rad/s) over 1 period, (2)angular acceleration (rad/s²)
    % over 1 period, load moment (3) (Nm) over 1 period, (4) range for
    % normalised transmission ratios 'N' [min max], (5) amount of intermediate
    % points for N    

    thetad = JNT.ad;
    thetadd = JNT.add;  
    ldmoment = JNT.m_abs;

    if invertload 
        ldmoment = -1*ldmoment;
    end
    
    omeganorm = zeros(n,1);
    momentnormRms = zeros(n,1);
    momentnormMax = zeros(n,1);
    N=linspace(Nrange(1),Nrange(2),n)';
    for i = 1:n
        omeganorm(i) = max(abs(thetad)).*N(i);
        momentnormRms(i) =  rms(assistfactor*ldmoment./N(i) + N(i).*thetadd);
        momentnormMax(i) = max(abs(assistfactor*ldmoment./N(i) + N(i).*thetadd));
    end
end