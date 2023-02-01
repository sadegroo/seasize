clear
clc
close all
load gaitdata.mat;
%keep only sagittal hip data
TBL = removevars(gaitdata,{'aKsag', 'aAsag', 'mKsag_abs',...
    'mKsag_max', 'mKsag_min', 'mKsag_rel', 'mKsag_rms', 'pKsag_abs',...
    'pKsag_max', 'pKsag_min', 'pKsag_rel', 'mAsag_abs', 'mAsag_max',...
    'mAsag_min', 'mAsag_rel', 'mAsag_rms', 'pAsag_abs', 'pAsag_max',...
    'pAsag_min', 'pAsag_rel'});

%% motion profile data
id=31;
stridetime = 1.2;
% construct gaitprofile instance: (name 'pat 31', 31st row of gaitdata
% table, 1.2s stride time, 100% of internal moment as assistance, [no time
% vector], invert load true)
gaitprofile_inst=gaitprofile(['pat ' char(string(id))],TBL(id,:),...
    stridetime,1,[],true);
gaitprofile_inst.plot;

%% motor data
load actuatordata.mat
% pick 6th motor from table
mot_id = 6;                                                 
% construct motor instance: 
% (38% RMS overload, 106% peak overload, [nominal voltage])
motor_inst = motor(actuatordata(mot_id,:),1.38,2.06,[]);  

%% gear evaluation
% construct geareval instance: 
% (motor_inst, gaitfprofile_inst, 100% efficiency, 100% assistance, 
% 'sea'=true
geareval_inst = geareval(motor_inst, gaitprofile_inst,1,1, 'sea',true); 
% Evaluate with: (Nmax=200, Npoints=500, Krange=[1,300], NK grid size =
% [150, 500])
geareval_inst.eval(200, 500,'Krange',[1,300], 'nNK',[150, 500])
% validate margins for N=81 and K=200(Nm/rad)
[margs, margfig] = geareval_inst.margins(81,'K',200);
% plot changed motor velocity, acceleration, torque and power profiles for N=81 and K=200(Nm/rad)
[results_sel, figsselected]=geareval_inst.plotSelected('N', 81, 'K', 200);
