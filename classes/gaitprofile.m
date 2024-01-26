classdef gaitprofile < motionprofile
    %GAITPROFILE represents one motion-load profile from gait data
    
    properties
    end
    
    methods
        function obj = gaitprofile(description, gaitdata, stridetime, assistload, assisttimevec, invertload,varargin)
            %GAITPROFILE 
            arguments (Input)
                description {mustBeTextScalar}
                gaitdata (1,:) table
                stridetime (1,1) double
                assistload double
                assisttimevec double 
                invertload (1,1) logical
            end
            arguments (Repeating)
                varargin
            end

            assert(isscalar(assistload)||isvector(assistload),'assistload must be scalar or a vector');
            N = length(gaitdata.aHsag(:));
            timevec = linspace(0,stridetime,N);                 % make time vector based on stride time

            if isscalar(assistload)
                loadvec = gaitdata.mHsag_abs(:)*assistload;
            elseif isvector(assistload)
                assert(length(assistload)==length(assisttimevec), 'assistload and assisttimevec must be of equal length')
                assisttimevec_norm = (assisttimevec - assisttimevec(1))/assisttimevec(end)*stridetime;       % make vector start at 0 and end at stridetime
                %interpollate assistload according to timevec
                loadvec = interp1(assisttimevec_norm,assistload,timevec,'spline');
            end

            anglevec = gaitdata.aHsag(:)*pi/180; % table in deg, convert to rad
            if invertload
                loadvec = (-1)*loadvec;
            end
            obj = obj@motionprofile(description, timevec, anglevec, loadvec, varargin{:});
        end
        
    end
end

