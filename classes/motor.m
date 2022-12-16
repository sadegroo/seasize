classdef motor < handle
    %MOTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        overloadfactor
        speedfactor
        inertia
        nomtorque
        rmstorque
        peaktorque
        stalltorque
        nomspeed
        noloadspeed
        km
        tc
    end
    
    methods
        function obj = motor(motordata, overloadfactor, usedvoltage)
            %MOTOR Construct an instance of this class
            %   Detailed explanation goes here
            arguments (Input)
                motordata (1,:) table
                overloadfactor double {mustBeScalarOrEmpty}
                usedvoltage double {mustBeScalarOrEmpty}
            end

            if isempty(overloadfactor)
                obj.overloadfactor = 1;
            else
                obj.overloadfactor = overloadfactor;
            end
            
            if isempty(usedvoltage)
                obj.speedfactor = 1;
            else
                obj.speedfactor = usedvoltage/motordata.NOMINALVOLTAGEV;
            end

            obj.name = char(string(motordata.NAME));
            obj.inertia = motordata.INERTIAkgm ;                                               % gearbox inertia is ignored
            obj.nomtorque = motordata.NOMINALTORQUENm;                                           % no overload continuous torque
            obj.rmstorque = obj.nomtorque * obj.overloadfactor;
            obj.peaktorque = motordata.PEAKTORQUENm;                                              % allowable peak torque (controller limited)                                      
            obj.stalltorque= motordata.STALLTORQUENm;                                            % not achievable in most setups
            obj.nomspeed = 	motordata.NOMINALSPEEDrpm * 2 * pi /60 * obj.speedfactor;                     % RPM to rad/s
            obj.noloadspeed = motordata.NOLOADSPEEDrpm * 2 * pi /60 * obj.speedfactor;
            obj.km = motordata.TORQUECFNmA;                                                     % motor torque constant
            obj.tc = motordata.WINDTCs;                                                                 % winding time constant chapter "3.11.2 Output Current Limitation according to I2t Method" of the "EPOS4 Firmware Specification.pdf"

        end
        
    end
end

