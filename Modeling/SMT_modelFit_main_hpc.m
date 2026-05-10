% SMT_modelFit_main_hpc
% write by XR @ March 4th 2024
% The main script to run the model fitting procedure

function SMT_modelFit_main_hpc(hpc_platform, StudyWord_ii, iGrp, iSub)
% Adapted for SLURM parallel execution (each job fits one subject)

if nargin < 4
    error('Usage: SeqMemTask_modelFit_main_hpc(hpc_platform, StudyWord_ii, iGrp, iSub)');
end
%%
clearvars -except hpc_platform StudyWord_ii iGrp iSub
clc

if isequal(hpc_platform, 'tardis')
    addpath(genpath('/home/mpib/ren/rxj-neurocode/AgingStudy/'));
    folder = '/home/mpib/ren/rxj-neurocode/AgingStudy';

elseif isequal(hpc_platform, 'hummel')
    addpath(genpath('/beegfs/u/bbc7806/AgingStudy/'));
    folder = '/beegfs/u/bbc7806/AgingStudy';

elseif isequal(hpc_platform, 'local')
    addpath(genpath('Aging-SeqMemTask/'));
    folder = '/Users/ren/Projects-NeuroCode/MyExperiment/Aging-SeqMemTask';
end
model_folder = [folder, '/SeqMemTask-modeling'];

fprintf('=== Running model fitting for %s, hpc %s, group %d, subject %d ===\n', hpc_platform, StudyWord_ii, iGrp, iSub);

%% Initialize parameters and subject lists
run(fullfile(model_folder, 'SMT_modelFit_init_subjInfo.m')); 
% This file should contain all the long setup code you posted (nSes, subjList_YA, subjList_OA, etc.),
% but without the subject or model loops.
% You can just cut the entire initialization section from your main script into that file.

%% Choose experiment
if isequal(StudyWord_ii, 'Experiment1')
    bhvDataDir = [folder, '/AgingReplay-OnlineData'];
    expList    = {'interleaved-Exp1'};
    suffixWord = expList{1};

elseif isequal(StudyWord_ii, 'Experiment2Interleaved')
    bhvDataDir = [folder, '/AgingReplay-OnlineData'];
    expList    = {'interleaved'};
    suffixWord = expList{1};

else
    error('Unknown StudyWord_ii');
end

%% Select group and subject
if iGrp == 1
    groupName = 'younger';
    subj_list  = subjList_YA;
    if exist('subj_conds_young', 'var')
        subj_conds = subj_conds_young;
    else
        subj_conds = zeros(size(subj_list,1),1);
    end
elseif iGrp == 2
    groupName = 'older';
    subj_list  = subjList_OA;
    if exist('subj_conds_old', 'var')
        subj_conds = subj_conds_old;
    else
        subj_conds = zeros(size(subj_list,1),1);
    end
else
    error('Invalid group index: must be 1 (YA) or 2 (OA)');
end

if iSub > size(subj_list,1)
    error('Subject index exceeds subject list length');
end

%% Call fitting code for this single subject
fprintf('>>> Group: %s | Subject: %d/%d\n', groupName, iSub, size(subj_list,1));
fprintf('>>> Running model fitting...\n');

SMT_fit_wrapper_hpc(folder, bhvDataDir, StudyWord_ii, groupName, subj_list, subj_conds, iGrp, iSub, expList, params);

fprintf('=== Done with %s | %s | sub-%d ===\n', StudyWord_ii, groupName, iSub);
exit; % exit MATLAB after finishing


