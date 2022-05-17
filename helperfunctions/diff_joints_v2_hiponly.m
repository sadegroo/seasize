function H = diff_joints_v2_hiponly(TBL,Tgait,flp,flp_load)
%DIFF_JOINTS Derivative and filtering of joint angles and torques

    % extract joint angles in degrees
    H.a = TBL.aHsag;
    H.name = 'Hip';
    
    % data vector size, timestep and sample frequency
    [~,pnts] = size(TBL.aHsag);
    timestep = Tgait/(pnts-1);
    fs = 1/timestep;
    
    % angular velocity (rad/s), angular acceleration (rad/s²)
    % through numerical differentiation
    [H.ad, H.add] = diff_vel_accel_v2(H.a*pi/180,flp);
    
    % extract load moments in Nm
    H.m_abs = lowpass(TBL.mHsag_abs',flp_load,fs)';
    
    %differentiatie torques
    [H.md, H.mdd] = diff_vel_accel_v2(H.m_abs,flp_load);
    
    H.m_rel = TBL.mHsag_rel;
    
    % extract powers in W
    H.p_abs = TBL.pHsag_abs;
    
    H.p_rel = TBL.pHsag_rel;

    function [vel,accel] = diff_vel_accel_v2(pos,flp)
        %DIFF_VEL_ACCEL numerical first and second pseudo-derivative for a row vector or 2D array of positions
        % <pos> must be a array or row vector. each column is a position
        % <period> represents the time needed to pass all positions
        % <flp> is a lowpass filter constant
        
        assert(isnumeric(pos),'First input <pos> must be a numeric array.')
        assert(numel(size(pos))<= 2,'First input <pos> must be 1-dimensional.')
        assert(size(pos,2)>=3, 'First input <pos> must contain at least 3 elements')
                
        % angular velocity (rad/s), , lowpass filtered
        vel = lowpass(gradient(pos,timestep),flp,fs);
        
        % angular acceleration (rad/s), lowpass filtered
        accel = lowpass(gradient(vel,timestep),flp,fs); 
        
    end

end

