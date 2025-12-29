% seqMemTask_BehvPaper_Figure2.m
% revision by XR @ Dec 21 2025
% ****** simplified based on seqMemTask_v2_anal_summary.m ******
% Plotting script for Figure 2 in the Behavioral Manuscript 1

clear
clc

%%
addpath('tight_subplot/');
addpath(genpath('HierarchicalCluster/'));
addpath(genpath('Aging-SeqMemTask/'));
% addpath(genpath('Violinplot-Matlab-master/'));
addpath('fdr_bh');

%% parameters
folder      = '/Users/ren/Projects-NeuroCode/MyExperiment/Aging-SeqMemTask';
bhvDataDir  = [folder, '/AgingReplay-OnlineData'];
statSaveDir = [folder, '/AgingReplay-StatAnal'];

RChunk = 0.35;
nSes    = 8*2; % each block contain 4 sequences, the 4 unique sequences will be repeated 8 times in 8 blocks
nImgSeq = 2; % 2 unique content sequence
nPosSeq = 2; % 2 unique position sequence
nPos    = 8; % 8 positions uniformly distributed on a circle
nImg    = 8; % 8 unique images
nTrans  = 5; % each sequence contains 5 transitions; 5 transitions = 5 categories
nEpi    = nImgSeq * nPosSeq * nSes;
nDtr    = 1; % number of distractor
postTn  = 8; % trials of post-testing
trlBlc  = 8; % each block has 8 trials
uniTrl  = 4; % 4 unique trials
nBlock  = nEpi / trlBlc; % 8 blocks
nComb   = 4; % 4 unique combinations between 2 content sequence and 2 position sequence
angCir = 0 : pi/50 : 2 * pi;
centerX = 0;
centerY = 0;
xCir   = RChunk * cos(angCir) + centerX;
yCir   = RChunk * sin(angCir) + centerY;

expList = {'interleaved', 'contentBlocked', 'positionBlocked'};
expId   = expList{1};
if isequal(expId, 'interleaved')
    subjList_young = {'5ad63c167f70c10001904bc5', '2023-08-30_17h17.39.428'; '5bdb51e1ba9b510001052364', '2023-08-30_15h12.00.151'; '5c4b06903566570001309394', '2023-08-30_16h55.13.543'; ...
                      '5d024a1fb58b6f001a58f74d', '2023-08-30_15h11.44.361'; '5d43404f1e6eef00011dec22', '2023-08-30_15h12.02.990'; '5ef25afb8ebcdf0b2b95d9cd', '2023-08-30_15h09.37.394'; ...
                      '5f15f96e54587538da27d452', '2023-08-30_15h40.43.668'; '5fd0c81fc79aef1882cbee94', '2023-08-30_16h25.12.136'; '60fecc838b1c231b1732cbb0', '2023-08-30_15h07.34.541'; ...
                      '601f93758d79b24eabff2e44', '2023-08-30_15h11.55.950'; '602fc5844525b3d343303a2a', '2023-08-30_14h05.56.283'; '604be8ac8e0c517878fd1d9f', '2023-08-30_14h05.21.886'; ...
                      '612ecc90331b627f7aaac5dc', '2023-08-30_16h23.19.001'; '614fca831894ddce32c1a342', '2023-08-30_15h20.32.432'; '615b5902e51bcad574d81203', '2023-08-30_15h18.25.454'; ...
                      '6016c8e7ea3f2387ae8b47d5', '2023-08-30_16h11.53.139'; '6103c08d411c6be73d9d78a7', '2023-08-30_15h20.28.394'; '6159f6b637bab134ea9bb92e', '2023-08-30_15h13.34.934'; ...
                      '61070b50a022d7360e46e985', '2023-08-30_16h08.26.937'; '61353c933f32fef782432cc7', '2023-08-30_15h10.45.788'; '605272be8568b6160f582f2e', '2023-08-30_14h38.14.393'; ...
                      '6107292e60892e4246db7425', '2023-08-30_15h12.11.729'; '61685478a9bd5239a9438f66', '2023-08-30_15h12.25.867'; '614831813dc412ccc8e2f563', '2023-08-30_15h31.25.179'};

    subjList_old   = {'5abb8dcb7ccedb0001b7f0d7', '2023-05-23_15h56.51.831'; '5be064114c6bd000013368f3', '2023-05-22_18h00.20.474'; '5c5df0475b87820001c4f21c', '2023-05-23_16h15.47.308'; ...
                      '5e9f0bc126557006ea49d1f4', '2023-05-23_16h58.16.506'; '5ea20fd571038c119083a8df', '2023-05-23_15h56.44.482'; '63b2d04ed0f53f75de4ba38e', '2023-05-23_16h30.06.501'; ...
                      '609a503448860549084c43ce', '2023-05-22_17h13.56.435'; '60534c39d754d351333bdd7c', '2023-05-23_16h15.50.361'; '597519f8262c480001bbaf8b', '2023-05-23_17h49.49.000'; ...
                      '61539b3fa541b182c0fadde1', '2023-05-26_11h37.17.532'; '574ce0a57fd0ec000db73aa6', '2023-05-26_12h33.38.572'; '55900dcffdf99b3f7aada3f5', '2023-05-26_10h08.56.439'; ...
                      '55e9aa1c735c45001043fbb6', '2023-05-26_17h56.53.674'; '64456ad3d3e7651a1dad232c', '2023-05-26_11h53.36.716'; '62aa26dd93252c8d69f7fc45', '2023-05-26_17h51.56.081'; ...
                      '5f53b958c8cfea6e2104c5b6', '2023-05-26_17h23.34.078'; '5f48e3d7f998433ac6356ad4', '2023-05-26_11h52.29.081'; '62f0f033178f89dd6f416590', '2023-05-26_17h21.29.494'; ...
                      '5c79a584670f87001646cef6', '2023-05-26_17h42.42.514'; '630be3605287a0f49b87c709', '2023-05-26_16h37.54.676'; '6121190671d1042b24d8d67b', '2023-05-26_16h16.28.267'; ...
                      '5c4cdcb14cb4630001ec4955', '2023-05-26_16h17.13.674'; '5f6e83419dd5cb3c85325fc6', '05-26-2023_16h37.34.939'}; % '63beebaa4c5884797ff00a98', '2023-05-26_17h39.52.151': attend contentBlocked

    ageList_younger    = [24, 22, 23, ...
                          25, 24, 21, ...
                          22, 23, 28, ...
                          21, 22, 23, ...
                          24, 30, 28, ...
                          27, 22, 24, ...
                          21, 20, 26, ...
                          26, 33, 28];
    genderList_younger = {'F', 'M', 'F', ...
                          'M', 'F', 'M', ... % the 4th participant: male(trans)
                          'M', 'M', 'F', ...
                          'M', 'F', 'M', ...
                          'M', 'M', 'F', ...
                          'F', 'F', 'F', ...
                          'F', 'M', 'F', ...
                          'M', 'F', 'M'}; 

    ageList_older    = [67, 68, 67, ...
                        69, 65, 66, ... % 5th participant: 64???
                        73, 73, 70, ...
                        65, 72, 72, ...
                        65, 69, 67, ...
                        67, 65, 68, ...
                        70, 72, 65, ...
                        67, 67];
    genderList_older = {'F', 'F', 'M', ...
                        'M', 'F', 'F', ...
                        'F', 'M', 'F', ...
                        'F', 'F', 'M', ...
                        'F', 'F', 'M', ...
                        'M', 'F', 'F', ...
                        'F', 'F', 'M', ...
                        'M', 'F'};

end
nGroup = 2; %% younger and older adults

%% quantify the age and gender
% ------ YAs ------
disp('---------- Age of YAs ----------');
[mean_age, sem_age] = Mean_and_Se(ageList_younger') 

disp('---------- Females in YAs----------');
sum(strcmp(genderList_younger, 'F'))
length(genderList_younger)

% ------ OAs ------
disp('---------- Age of OAs ----------');
[mean_age, sem_age] = Mean_and_Se(ageList_older') 

disp('---------- Females in OAs----------');
sum(strcmp(genderList_older, 'F'))
length(genderList_older)

%% data analysis
% overall accuracy & RTs across all blocks
acc_group     = cell(1, nGroup);
rt_group      = cell(1, nGroup);
acc_dim_group = cell(1, nGroup);
% overall accuracy for content and position in marginal report and
% reconstruction report (the latter in the third and first order)
acc_marginalDim_inJointRep_group = cell(1, nGroup);
% accuracy for each sub-block
acc_blc_group          = cell(1, nGroup);
% accuracy for content, position and reconstruction for two displaying orders
acc_subj_order_group   = cell(1, nGroup);
acc_subj_orderUP_group = cell(1, nGroup);
% percentage of fully correct trials
trialPerc_group = cell(1, nGroup);
% accuracy & RTs in each trial
acc_trial_group = cell(1, nGroup);
% ------accuracy for subsequent LMM analysis------
acc_lmm_group          = [];
choice_lmm_group       = [];
acc_change_conds_group = [];

% ------Binding related calculation------
binds_conPctr_group = cell(2, nGroup); % 2: marginal and reconstruction reports
% ------Transition evidence------
transAcc_count_group            = cell(2, nGroup); % 2: marginal and joint
transAcc_count_firstHalf_group  = cell(2, nGroup);
transAcc_count_secondHalf_group = cell(2, nGroup);

% ----------Integrate the following measures into one variable and use it for the subsequent correlation analysis with model parameters----------
allMeas_inOne_group = cell(1, nGroup);

% ----------False alarm from the lure stimuli----------
FA_lure_group = cell(1, nGroup);

suffixWord = expId;
for iGrp = 1 : nGroup %% younger and older adults
    if iGrp == 1
        groupName = 'younger';
        subj_list = subjList_young;
    elseif iGrp == 2
        groupName = 'older';
        subj_list = subjList_old;
    end
    subLen   = length(subj_list);
    subjPath = [bhvDataDir, '/AgingReplay-v2-', suffixWord, '/', suffixWord, '-', groupName, '/'];

    % ------accuracy & RTs in each trial------
    acc_trial_subj = nan(subLen, 5, (nEpi + postTn)); % 5: (1) content report; (2) position report; (3) both report; (4) content in both report; (5) position in both report;
    rt_trial_subj  = nan(subLen, 5, (nEpi + postTn)); 

    % ------overall accuracy & RTs across all blocks------
    acc_subj       = nan(subLen, 14); % post-test: (6) both; (7) content in both report; (8) position in both report;
    rt_subj        = nan(subLen, 8);

    % ------accuracy in each sub-block------
    acc_blc_subj   = nan(subLen, (nBlock+1)*2, 3); % 3: (1) content report; (2) position report; (3) both report;

    % ------accuracy in reconstruction report------
    acc_marginalDim_inJointRep_sub = nan(subLen, 2, 3); % 2: content and position reports; 3: all trials, recons report in the 3rd order, only-recons trials

    % ------accuracy for content, position and reconstruction separately for two dispalying order------
    acc_subj_order   = nan(subLen, 2, 3); % 2nd dimension: first and second report condition; 3rd dimension: content, position and reconstruction
    acc_subj_orderUp = nan(subLen, 3, 3); % 2nd dimension: first, second, and third report order; 3rd dimension: content, position and reconstruction

    % ------percentage of the fully correct trials------
    trialPerc_subj   = nan(subLen, 4); % 4: content, position, reconstruction, post-test

    % ------accuracy for subsequent LMM analysis------
    acc_lmm_subj          = [];
    choice_lmm_subj       = []; % choice in each slot: 1-correct, 0-incorrect
    acc_change_conds_subj = [];

    % ------overall accuracy in 11 measurements------
    acc_dim_subj          = nan(subLen, 11); % post-test: (6) both; (7) content in both report; (8) position in both report; (9) both in recons-only; (10) content in recons-only; (11) position in recons-only

    % ------Binidng related measures------
    binds_conPctr_marg_subj = nan(subLen, 4, 4); % the 2nd 2: (1-2) proportion: item on position and position on item; (3-4) detected response numbers
    binds_conPctr_join_subj = nan(subLen, 4, 4); 
    % ------Transition evidence------
    transAcc_count_marg_subj = nan(subLen, 4, 2); % 4: 4 different counts; 2: item and location
    transAcc_count_join_subj = nan(subLen, 4, 2);
    % ------Transition evidence for the 1st half and 2nd half per retrieval test------
    transAcc_count_marg_firstHalf_subj  = nan(subLen, 4, 2);
    transAcc_count_join_firstHalf_subj  = nan(subLen, 4, 2);
    transAcc_count_marg_secondHalf_subj = nan(subLen, 4, 2);
    transAcc_count_join_secondHalf_subj = nan(subLen, 4, 2);

    % ------False alarm from the lure stimuli------
    FA_lure_subj   = nan(subLen, 2, 2); % first 2: item and location; second 2: partial and full retrieval

    %%
    for iSub = 1 : subLen
        subjBv = subj_list{iSub, 1};
        subjTm = subj_list{iSub, 2};
        seqMem_subj = readtable([subjPath, subjBv, '_EpisodicMemoryTask-', suffixWord, '_', subjTm, '.csv']);
        %%
        %%% trial test order: 0-content firstly; 1-position firstly
        testOrd = seqMem_subj.trlTestOrd;
        testOrd = testOrd(~isnan(testOrd));
        %%% unique combination between two content and two position sequence
        uniCombSeq_tmp = seqMem_subj.trlComb;
        uniCombSeq_tmp = uniCombSeq_tmp(~cellfun('isempty', uniCombSeq_tmp));
        uniCombSeq = nan(nEpi, 2);
        for i =  1 : nEpi
            uniCombSeq(i, :) = str2num(uniCombSeq_tmp{i});
        end
        %%% shortWTI index
        shortWTI = seqMem_subj.sWTImark;
        shortWTI = shortWTI(~isnan(shortWTI));
        %%% reconstruction index
        reconsOnly = seqMem_subj.reconsMark;
        reconsOnly = reconsOnly(~isnan(reconsOnly));
        %%% shortWTI+recons
        shortWTI_recons = seqMem_subj.sWTIrecons;
        shortWTI_recons = shortWTI_recons(~isnan(shortWTI_recons));

        %%% longWTI+recons
        longWTI_recons = seqMem_subj.lWTIrecons;
        longWTI_recons = longWTI_recons(~isnan(longWTI_recons));

        %% positions of each slot
        posX_tmp = seqMem_subj.locSeqXTrl;
        posX_tmp = posX_tmp(~cellfun('isempty', posX_tmp));
        posY_tmp = seqMem_subj.locSeqYTrl;
        posY_tmp = posY_tmp(~cellfun('isempty', posY_tmp));
        posX_col = cell(nEpi, 1);
        posY_col = cell(nEpi, 1);
        for i =  1 : nEpi
            posX_col{i} = str2num(posX_tmp{i});
            posY_col{i} = str2num(posY_tmp{i});
        end

        %% ---------- content report ----------
        conTrue_tmp_raw = seqMem_subj.conReportTrue;
        conTrue_tmp_raw = conTrue_tmp_raw(~cellfun('isempty', conTrue_tmp_raw));
        conTrue_tmp = cell(nEpi, 1);
        conTrue_tmp(reconsOnly == 0) = conTrue_tmp_raw;
        conRep_tmp_raw  = seqMem_subj.conReportOrd;
        conRep_tmp_raw  = conRep_tmp_raw(~cellfun('isempty', conRep_tmp_raw));
        conRep_tmp = cell(nEpi, 1);
        conRep_tmp(reconsOnly == 0) = conRep_tmp_raw ;
        conRT_tmp_raw   = seqMem_subj.conRTs;
        conRT_tmp_raw   = conRT_tmp_raw(~cellfun('isempty', conRT_tmp_raw));
        conRT_tmp = cell(nEpi, 1);
        conRT_tmp(reconsOnly == 0) = conRT_tmp_raw;
        conTrue_col = cell(nEpi, 1);
        conRep_col  = cell(nEpi, 1);
        conRT_col   = cell(nEpi, 1);
        for i =  1 : nEpi
            if reconsOnly(i) == 0 %% non reconstruction only trial
                conTrue_col{i} = str2num(conTrue_tmp{i});
                conRep_noRef   = str2num(conRep_tmp{i});
                conRep_col{i}  = conRep_noRef;
                %%% get RT for each item by subtracting the RT of its
                %%% predecessor
                conRT_noRef = str2num(conRT_tmp{i});
                conRT_Ref   = nan(1, (nTrans + nDtr));
                for j = 1 : nTrans
                    if j == 1
                        conRT_j = conRT_noRef(conRep_noRef == j);
                    else
                        conRT_j = conRT_noRef(conRep_noRef == j) - conRT_noRef(conRep_noRef == (j - 1));
                    end
                    conRT_Ref(conRep_noRef == j) = conRT_j;
                end
                conRT_col{i} = conRT_Ref;
            end
        end

        %% ---------- position report ----------
        locTrue_tmp_raw = seqMem_subj.locReportTrue;
        locTrue_tmp_raw = locTrue_tmp_raw(~cellfun('isempty', locTrue_tmp_raw));
        locTrue_tmp = cell(nEpi, 1);
        locTrue_tmp(reconsOnly == 0) = locTrue_tmp_raw;
        locRep_tmp_raw  = seqMem_subj.locReportOrd;
        locRep_tmp_raw  = locRep_tmp_raw(~cellfun('isempty', locRep_tmp_raw));
        locRep_tmp  = cell(nEpi, 1);
        locRep_tmp(reconsOnly == 0) = locRep_tmp_raw;
        locRT_tmp_raw   = seqMem_subj.locRTs;
        locRT_tmp_raw   = locRT_tmp_raw(~cellfun('isempty', locRT_tmp_raw));
        locRT_tmp   = cell(nEpi, 1);
        locRT_tmp(reconsOnly == 0) = locRT_tmp_raw;
        locTrue_col = cell(nEpi, 1);
        locRep_col  = cell(nEpi, 1);
        locRT_col   = cell(nEpi, 1);
        for i =  1 : nEpi
            if reconsOnly(i) == 0 %% non reconstruction only trial
                locTrue_col{i} = str2num(locTrue_tmp{i});
                locRep_noRef   = str2num(locRep_tmp{i});
                locRep_col{i}  = locRep_noRef;
                %%% get RT for each item
                locRT_noRef = str2num(locRT_tmp{i});
                locRT_Ref   = nan(1, (nTrans + nDtr));
                for j = 1 : nTrans
                    if j == 1
                        locRT_j = locRT_noRef(locRep_noRef == j);
                    else
                        locRT_j = locRT_noRef(locRep_noRef == j) - locRT_noRef(locRep_noRef == (j - 1));
                    end
                    locRT_Ref(locRep_noRef == j) = locRT_j;
                end
                locRT_col{i} = locRT_Ref;
            end
        end

        %% ---------- reconstruction report ----------
        bothTrue_tmp = seqMem_subj.bothReportTrue;
        bothTrue_tmp = bothTrue_tmp(~cellfun('isempty', bothTrue_tmp));
        bothRep_tmp  = seqMem_subj.bothReportOrd;
        bothRep_tmp  = bothRep_tmp(~cellfun('isempty', bothRep_tmp));
        bothRT_tmp   = seqMem_subj.bothRTs;
        bothRT_tmp   = bothRT_tmp(~cellfun('isempty', bothRT_tmp));
        %%% accuracy and RT calculation based on integration of content and
        %%% position
        bothTrue_col = cell(nEpi + postTn, 1);
        bothRep_col  = cell(nEpi + postTn, 1);
        bothRT_col   = cell(nEpi + postTn, 1);
        for i =  1 : (nEpi + postTn)
            bothTrue_col{i} = str2num(bothTrue_tmp{i});
            bothCol = str2num(bothRep_tmp{i});
            if i <= nEpi
                bothCol = reshape(bothCol, 2, (nTrans+nDtr));
            else
                bothCol = reshape(bothCol, 2, nTrans);
            end
            bothRep_col{i} = bothCol;
            bothRep_con = bothCol(1, :); %% report order of content
            bothRep_loc = bothCol(2, :); %% report order of position

            %%% In the raw data, the RTs is aligned based on the position
            %%% report (the 2nd row of bothCol)
            bothRT_noRef = str2num(bothRT_tmp{i});
            if i <= nEpi
                bothRT_Ref = nan(2, (nTrans + nDtr)); % 1st row: content; 2nd row: position;
            else
                bothRT_Ref = nan(2, nTrans); % 1st row: content; 2nd row: position;
            end
            for j = 1 : nTrans
                if j == 1
                    bothRT_j = bothRT_noRef(bothRep_loc == j);
                else
                    bothRT_j = bothRT_noRef(bothRep_loc == j) - bothRT_noRef(bothRep_loc == (j - 1));
                end
                if ~isempty(bothRT_j)
                    bothRT_Ref(1, bothRep_con == j) = bothRT_j;
                    bothRT_Ref(2, bothRep_loc == j) = bothRT_j;
                end
            end
            bothRT_col{i} = bothRT_Ref;
        end

        %% !!!!!!!!!! Accuracy & RT calculation!!!!!!!!!!
        %% trial-by-trial accuracy & RT
        % single content and location report
        % content and report
        choice_con_iSub = nan(nEpi, 5); % 48 marginal report trials * 5 transitions
        choice_pos_iSub = nan(nEpi, 5);
        for i =  1 : nEpi
            if reconsOnly(i) == 0 %% non reconstruction only trial
                % content report
                conTrue_i = conTrue_col{i};
                conRep_i  = conRep_col{i};
                conRT_i   = conRT_col{i};
                conRep_i  = conRep_i(conTrue_i ~= 6);
                conRT_i   = conRT_i(conTrue_i ~= 6);
                conTrue_i = conTrue_i(conTrue_i ~= 6);
                acc_trial_subj(iSub, 1, i) = (sum(conRep_i == conTrue_i)) / nTrans;
                choice_con_iSub(i, :) = (conRep_i == conTrue_i); % 1-correct; 0-incorrect;
                %%% RT calculation based on single correct item
                if sum(conRep_i == conTrue_i) ~= 0
                    rt_trial_subj(iSub, 1, i) = nanmean(conRT_i(conRep_i == conTrue_i));
                end

                % position report
                locTrue_i = locTrue_col{i};
                locRep_i  = locRep_col{i};
                locRT_i   = locRT_col{i};
                locRep_i  = locRep_i(locTrue_i ~= 6);
                locRT_i   = locRT_i(locTrue_i ~= 6);
                locTrue_i = locTrue_i(locTrue_i ~= 6);
                acc_trial_subj(iSub, 2, i) = (sum(locRep_i == locTrue_i)) / nTrans;
                choice_pos_iSub(i, :) = (locRep_i == locTrue_i); % 1-correct; 0-incorrect;
                %%% RT calculation based on single correct item
                if sum(locRep_i == locTrue_i) ~= 0
                    rt_trial_subj(iSub, 2, i) = nanmean(locRT_i(locRep_i == locTrue_i));
                end
            end
        end
        choice_con_iSub(reconsOnly == 1, :) = [];
        choice_pos_iSub(reconsOnly == 1, :) = [];

        % reconstruction
        choice_both_iSub      = nan(nEpi, 5);
        choice_both_item_iSub = nan(nEpi, 5);
        choice_both_loc_iSub  = nan(nEpi, 5);
        for i =  1 : (nEpi + postTn)
            bothTrue_i = bothTrue_col{i};
            bothRep_i  = bothRep_col{i}; % 2 row: content and position
            bothRT_i   = bothRT_col{i};  % 2 row: content and position
            %%% accuracy integrate both content and position
            acc_j = zeros(nTrans, 1);
            for j = 1 : nTrans
                if bothRep_i(1, j) == j && bothRep_i(2, j) == j
                    acc_j(j) = 1;
                end
            end
            acc_trial_subj(iSub, 3, i) = sum(acc_j) / nTrans;
            if i <= nEpi
                choice_both_iSub(i, :) = acc_j;
                % ------ Marginal reports in the reconstruction report
                % ------
                choice_both_item_iSub(i, :) = (bothRep_i(1, 1 : nTrans) == (1 : 1 : nTrans));
                choice_both_loc_iSub(i, :)  = (bothRep_i(2, 1 : nTrans) == (1 : 1 : nTrans));

            end
            %%% RT calculation based on single correct item
            if sum(acc_j) ~= 0
                bothRT_ij = bothRT_i(1, 1 : nTrans); % only the items consistent between content and position are valid
                rt_trial_subj(iSub, 3, i) = nanmean(bothRT_ij(acc_j == 1));
            end

            %%% accuracy for content and position separately
            true_Tmp = 1 : 1 : nTrans;
            bothRT_ij = bothRT_i(:, 1 : nTrans);
            acc_trial_subj(iSub, 4, i) = (sum(bothRep_i(1, 1 : nTrans) == true_Tmp)) / nTrans;
            acc_trial_subj(iSub, 5, i) = (sum(bothRep_i(2, 1 : nTrans) == true_Tmp)) / nTrans;
            %%% RT calculation based on single correct item
            rt_trial_subj(iSub, 4, i) = nanmean(bothRT_ij(1, (bothRep_i(1, 1 : nTrans) == true_Tmp)));
            rt_trial_subj(iSub, 5, i) = nanmean(bothRT_ij(2, (bothRep_i(2, 1 : nTrans) == true_Tmp)));
        end

        %% overall accuracy & RT
        % 64 episodes before numJdgTask
        % acc_trial_iSub = squeeze(acc_trial_subj(iSub, :, 1 : nEpi));
        % acc_trial_iSub = acc_trial_iSub(:, reconsOnly == 0);
        % acc_subj(iSub, 1 : 5) = nanmean(acc_trial_iSub, 2);
        acc_subj(iSub, 1 : 5) = squeeze(nanmean(acc_trial_subj(iSub, :, 1 : nEpi), 3));
        rt_subj(iSub, 1 : 5)  = squeeze(nanmean(rt_trial_subj(iSub, :, 1 : nEpi), 3));
        % 8 episodes after numJdgTask
        acc_subj(iSub, 6 : 8) = squeeze(nanmean(acc_trial_subj(iSub, 3 : end, (nEpi + 1) : end), 3));
        rt_subj(iSub, 6 : 8)  = squeeze(nanmean(rt_trial_subj(iSub, 3 : end, (nEpi + 1) : end), 3));

        %% calculate the predicted joint accuracy from either marginal reports or the marginal dimension in the reconstruction reports
        % added by rxj @ Sep 14 2024
        % only non-reconstruction only trials
        acc_trial_con_iSub = squeeze(acc_trial_subj(iSub, 1, 1 : nEpi));
        acc_trial_pos_iSub = squeeze(acc_trial_subj(iSub, 2, 1 : nEpi));
        acc_trial_rec_iSub = squeeze(acc_trial_subj(iSub, 3, 1 : nEpi));
        acc_trial_con_iSub = acc_trial_con_iSub(reconsOnly == 0);
        acc_trial_pos_iSub = acc_trial_pos_iSub(reconsOnly == 0);
        acc_trial_rec_iSub = acc_trial_rec_iSub(reconsOnly == 0);
        %%% ------The non-reconsOnly trials: 3 consecutive reports (same as Experiment 1)------
        % ----content----
        acc_subj(iSub, 9)  = nanmean(acc_trial_con_iSub);
        % ----position----
        acc_subj(iSub, 10) = nanmean(acc_trial_pos_iSub);
        % ----reconstruction----
        acc_subj(iSub, 11) = nanmean(acc_trial_rec_iSub);
        % ----independent prediction----
        acc_subj(iSub, 12) = nanmean(acc_trial_con_iSub .* acc_trial_pos_iSub);
        acc_subj(iSub, 13) = nanmean(squeeze(acc_trial_subj(iSub, 4, 1 : nEpi)) .* squeeze(acc_trial_subj(iSub, 5, 1 : nEpi)));

        %% overall accuracy: but separate the single dimension report in reconstruction only trials
        % added by rxj @ 06/06/2023
        %acc_dim_subj = nan(subLen, 11); % post-test: (6) both; (7) content in both report; (8) position in both report; (9) both in recons-only; (10) content in recons-only; (11) position in recons-only
        acc_trial_iSub  = squeeze(acc_trial_subj(iSub, :, 1 : nEpi)); %% 5measures * 64nEpi
        acc_dim_subj(iSub, 1 : 2) = nanmean(acc_trial_iSub(1 : 2, :), 2);
        acc_dim_subj(iSub, [3, 4, 5])   = nanmean(acc_trial_iSub([3, 4, 5], reconsOnly == 0), 2);
        acc_dim_subj(iSub, [6, 7, 8])   = nanmean(acc_trial_iSub([3, 4, 5], reconsOnly == 1), 2); %% reconstruction only trials
        acc_dim_subj(iSub, [9, 10, 11]) = squeeze(nanmean(acc_trial_subj(iSub, 3 : end, (nEpi + 1) : end), 3));

        %% accuracy for marginal dimension in the joint reports
        % 3 cases: (1) all trials; (2) recons in the "3-reports" trials; (3) reconstruction only trial;
        % added by rxj @ 02/04/2024
        % acc_marginalDim_inJointRep_sub = nan(subLen, 2, 3);
        % --------marginal dimension in joint report across all trials--------
        acc_marginalDim_inJointRep_sub(iSub, 1 : 2, 1) = nanmean(acc_trial_iSub([4, 5], :), 2);
        % --------marginal dimension in joint report across the "3-reports trials"--------
        acc_marginalDim_inJointRep_sub(iSub, 1 : 2, 2) = nanmean(acc_trial_iSub([4, 5], reconsOnly == 0), 2);
        % --------marginal dimension in joint report across the "reconstruction only trials"--------
        acc_marginalDim_inJointRep_sub(iSub, 1 : 2, 3) = nanmean(acc_trial_iSub([4, 5], reconsOnly == 1), 2);

        %% block-wise pattern
        % acc_blc_subj = nan(subLen, (nBlock+1)*2, 3); % 3: (1) content report; (2) position report; (3) both report;
        if isequal(suffixWord, 'interleaved')
            for iBlc = 1 : ((nBlock + 1) * 2) % 1: post-test block
                iTrl = (iBlc - 1) * uniTrl + 1 : iBlc * uniTrl;
                if iBlc <= (nBlock * 2)
                    acc_blc_subj(iSub, iBlc, :) = squeeze(nanmean(acc_trial_subj(iSub, 1 : 3, iTrl), 3));
                else
                    acc_blc_subj(iSub, iBlc, 3) = squeeze(nanmean(acc_trial_subj(iSub, 3, iTrl), 3));
                end
            end
        end

        %% order effect: compare the accuacy within a dimension (content or position) between the 1st and 2nd report after the encoding stage
        % added by rxj @ 01/22/2024
        % acc_trial_subj = nan(subLen, 5, (nEpi + postTn)); % 5: (1) content report; (2) position report; (3) both report; (4) content in both report; (5) position in both report;
        acc_trialReport = (squeeze(acc_trial_subj(iSub, 1 : 3, 1 : nEpi)))'; % 64 trials * 3 columns
        trial_marginal  = testOrd(reconsOnly == 0); % delete the reconstruction only trials
        acc_marginal    = acc_trialReport((reconsOnly == 0), :); % 48 trials * 3 columns
        acc_reconsOnly  = acc_trialReport((reconsOnly == 1), :); % 16 trials * 3 columns (the first 2 columns are NAN)

        % acc_subj_order = nan(subLen, 2, 3); % 2nd dimension: first and second report order; 3rd dimension: content, position and reconstruction
        % content report
        acc_subj_order(iSub, 1, 1) = nanmean(acc_marginal(trial_marginal == 0, 1));
        acc_subj_order(iSub, 2, 1) = nanmean(acc_marginal(trial_marginal == 1, 1));
        % position report
        acc_subj_order(iSub, 1, 2) = nanmean(acc_marginal(trial_marginal == 1, 2));
        acc_subj_order(iSub, 2, 2) = nanmean(acc_marginal(trial_marginal == 0, 2));
        % reconstruction report
        acc_subj_order(iSub, 1, 3) = nanmean(acc_marginal(trial_marginal == 0, 3)); % content-position-reconstruction
        acc_subj_order(iSub, 2, 3) = nanmean(acc_marginal(trial_marginal == 1, 3)); % position-content-reconstruction

        %% accuracy for content, position and reconstruction for two displaying orders: 2nd version
        % added by rxj @ 01/24/2024
        % plotting according to the reporting order after the encoding stage
        % 3 kinds of reporting order:
        % (1) first order: content in content-first trials, posiiton in position-first trials,
        % recon report in reconstruction only trials;
        % (2) second order: pos in content-first trials, con in position-first trials;
        % (3) third order: recon report in non recon-only trials
        %acc_subj_orderUp = nan(subLen, 3, 3); % 2nd dimension: first, second, and third report order; 3rd dimension: content, position and reconstruction
        % ------first order------
        % content-first, position-first, and recons-only trials
        acc_subj_orderUp(iSub, 1, 1) = nanmean(acc_marginal(trial_marginal == 0, 1));
        acc_subj_orderUp(iSub, 1, 2) = nanmean(acc_marginal(trial_marginal == 1, 2));
        acc_subj_orderUp(iSub, 1, 3) = nanmean(acc_reconsOnly(:, 3));
        % ------second order------
        acc_subj_orderUp(iSub, 2, 1) = nanmean(acc_marginal(trial_marginal == 1, 1));
        acc_subj_orderUp(iSub, 2, 2) = nanmean(acc_marginal(trial_marginal == 0, 2));
        % ------third order------
        acc_subj_orderUp(iSub, 3, 3) = nanmean(acc_marginal(:, 3));
        
        % ------Independent hypothesis prediction in the first reporting
        % window from the marginal reports------
        acc_subj(iSub, 14) = acc_subj_orderUp(iSub, 1, 1) * acc_subj_orderUp(iSub, 1, 2);

        %% percentage of fully correct trials
        con_trials = squeeze(acc_trial_subj(iSub, 1, 1 : nEpi));
        pos_trials = squeeze(acc_trial_subj(iSub, 2, 1 : nEpi));
        rec_trials = squeeze(acc_trial_subj(iSub, 3, :));
        trialPerc_subj(iSub, 1) = length(find(con_trials == 1)) / (length(find(reconsOnly == 0)));
        trialPerc_subj(iSub, 2) = length(find(pos_trials == 1)) / (length(find(reconsOnly == 0)));
        trialPerc_subj(iSub, 3) = length(find(rec_trials(1 : nEpi) == 1)) / nEpi; % length(find(rec_trials == 1)) / (nEpi + postTn);
        trialPerc_subj(iSub, 4) = length(find(rec_trials(nEpi + 1 : end) == 1)) / postTn;

        %% save accuracy and choice for each subject for subsequent LMM and GLMM analysis
        acc_trial_col = (squeeze(acc_trial_subj(iSub, [1, 2, 3], 1 : nEpi)))'; % nEpi * 3 meas
        acc_col    = []; %% accuracy
        trlCor_col = []; %% fully correct: 1-yes, 0-no
        meas_col   = []; %% 3 measures
        trlNo_col  = [];
        repOrd_col = []; %% report orders after the encoding stage: 0-content+position+recon, 1-position+content+recon, 2-recon only
        for iM = 1 : 3 %% 3 measures: content, position and reconstruction
            acc_iM = acc_trial_col(:, iM);
            trlNo_iM = (1 : 1 : nEpi)';
            trlNo_iM(isnan(acc_iM)) = [];
            repOrd_iM = testOrd;
            repOrd_iM(isnan(acc_iM)) = []; % for the reconstruction only trials, the value in content and position reports will be NAN
            if iM == 1    % content report
                repOrd_iM(repOrd_iM == 1) = 2; % second report after the encoding stage
                repOrd_iM(repOrd_iM == 0) = 1; % first report after the encoding stage
            elseif iM == 2 % position report
                repOrd_iM(repOrd_iM == 1) = 1; % first report after the encoding stage
                repOrd_iM(repOrd_iM == 0) = 2; % second report after the encoding stage
            elseif iM == 3 % reconstruction report
                repOrd_iM(reconsOnly == 0) = 3; % third report after the encoding stage: (1) con+pos+recons; (2) pos+con+recons;
                repOrd_iM(reconsOnly == 1) = 1; % first report after the encoding stage
            end
            acc_iM(isnan(acc_iM)) = [];
            acc_col   = [acc_col; acc_iM];
            trlCor_iM = acc_iM;
            trlCor_iM(trlCor_iM ~= 1) = 0;
            trlCor_col = [trlCor_col; trlCor_iM];
            meas_col   = [meas_col; repmat(iM - 1, length(acc_iM), 1)];
            trlNo_col  = [trlNo_col; trlNo_iM];
            repOrd_col = [repOrd_col; repOrd_iM];
        end
        group_Col = repmat(iGrp - 1, length(acc_col), 1);
        subj_Col  = repmat(iSub - 1, length(acc_col), 1);
        acc_lmm_subj = [acc_lmm_subj; group_Col, subj_Col, trlNo_col, meas_col, acc_col, trlCor_col, repOrd_col];

        %%  1. overall choice for Figure 3B: save choice in each slot for every subject for subsequent GLMM analysis
        % added by rxj @ April 20 2025: also labeling the reporting window
        %%% 2. choice according to chronological order (reporting window after the encoding stage): save choice in each slot for every subject for subsequent GLMM analysis
        % added by rxj @ 01/30/2024
        % correct (1) or incorrect (0) in each slot, each report and each trial
        %          choice_con_iSub = nan(nEpi, 5); % 48 marginal report trials * 5 transitions
        %          choice_pos_iSub = nan(nEpi, 5);
        %          choice_both_iSub = nan(nEpi, 5);
        choice_slot_con  = reshape(choice_con_iSub', [size(choice_con_iSub, 1) * nTrans, 1]);
        choice_slot_pos  = reshape(choice_pos_iSub', [size(choice_pos_iSub, 1) * nTrans, 1]);
        choice_slot_both = reshape(choice_both_iSub', [nEpi * nTrans, 1]);
        choice_slot_both_con = reshape(choice_both_item_iSub', [size(choice_both_item_iSub, 1) * nTrans, 1]);
        choice_slot_both_pos = reshape(choice_both_loc_iSub', [size(choice_both_loc_iSub, 1) * nTrans, 1]);

        % label the reporting window for each report: 0-first, 1-second,
        % 2-third (only for reconstruction in the 3-report trials)
        testOrd_threeReports_con = (testOrd(reconsOnly == 0))';  % one row; 0-item report in 1st reporting window; 1-item in 2nd reporting window;
        testOrd_threeReports_pos = testOrd_threeReports_con - 1;
        testOrd_threeReports_pos(testOrd_threeReports_pos == -1) = 1; % 0-location report in 1st reporting window; 1-location in 2nd reporting window;

        testOrd_both = (reconsOnly - 1)';     % 0: the 1st reporting window
        testOrd_both(testOrd_both == -1) = 2; % 2: the 3rd reporting window

        repWin_slot_con  = reshape(repmat(testOrd_threeReports_con, nTrans, 1), [length(testOrd_threeReports_con) * nTrans, 1]);
        repWin_slot_pos  = reshape(repmat(testOrd_threeReports_pos, nTrans, 1), [length(testOrd_threeReports_pos) * nTrans, 1]);
        repWin_slot_both = reshape(repmat(testOrd_both, nTrans, 1), [length(testOrd_both) * nTrans, 1]);

        choice_col = [];
        repWin_col = [];
        meas_col   = []; %% 3 measures
        trlNo_col  = [];
        slotNo_col = []; % each sequence contains 5 slots
        % added by rxj @ May 3rd 2025
        slow_fast_label = []; % 0-fast; 1-slow
        for iM = 1 : 3 %% 3 measures: content, position and reconstruction
            trlNo_iM = repmat((1 : 1 : nEpi)', 1, nTrans);
            if iM == 1     % content report
                choice_slot_iM = choice_slot_con;
                repWin_slot_iM = repWin_slot_con;
                trlNo_iM(reconsOnly == 1, :) = [];
            elseif iM == 2 % position report
                choice_slot_iM = choice_slot_pos;
                repWin_slot_iM = repWin_slot_pos;
                trlNo_iM(reconsOnly == 1, :) = [];
            elseif iM == 3 % reconstruction report
                choice_slot_iM = choice_slot_both;
                repWin_slot_iM = repWin_slot_both;
            end
            choice_col = [choice_col; choice_slot_iM];
            repWin_col = [repWin_col; repWin_slot_iM];
            meas_col   = [meas_col; repmat(iM - 1, length(choice_slot_iM), 1)];
            trlNo_col  = [trlNo_col; reshape(trlNo_iM', [size(trlNo_iM, 1) * size(trlNo_iM, 2), 1])];
            slotNo_iM  = repmat((1 : 1 : nTrans) - 1, size(trlNo_iM, 1), 1);
            slotNo_col = [slotNo_col; reshape(slotNo_iM', [size(slotNo_iM, 1) * size(slotNo_iM, 2), 1])];

        end
        group_Col = repmat(iGrp - 1, length(choice_col), 1);
        subj_Col  = repmat(iSub - 1, length(choice_col), 1);
        choice_lmm_subj = [choice_lmm_subj; group_Col, subj_Col, trlNo_col, meas_col, slotNo_col, choice_col, repWin_col];

        %% calculate the decline percentage for each report and import to LMM to see if blocked design decline reduces
        % added by rxj @ 02/01/2024
        % !!!!!!!!!! To make this consistent with Exp1, we should consider
        % using the three-retrieval-trials !!!!!!!!!! 
        if iGrp == 2 % older group
            acc_trial_col = (squeeze(acc_trial_subj(iSub, [1, 2, 3], 1 : nEpi)))'; % nEpi * 3 meas
            % ----------calculate accuracy change relative to the younger----------
            acc_trial_col(:, 1) = (acc_trial_col(:, 1) -  mean_young(1)) ./ mean_young(1);
            acc_trial_col(:, 2) = (acc_trial_col(:, 2) -  mean_young(2)) ./ mean_young(2);
            acc_trial_col(:, 3) = (acc_trial_col(:, 3) -  mean_young(3)) ./ mean_young(3);
            acc_col    = []; %% accuracy
            trlCor_col = []; %% fully correct: 1-yes, 0-no
            meas_col   = []; %% 3 measures
            trlNo_col  = [];
            repOrd_col = []; %% report orders after the encoding stage: 0-content+position+recon, 1-position+content+recon, 2-recon only
            for iM = 1 : 3 %% 3 measures: content, position and reconstruction
                acc_iM = acc_trial_col(:, iM);
                trlNo_iM = (1 : 1 : nEpi)';
                trlNo_iM(isnan(acc_iM)) = [];
                repOrd_iM = testOrd;
                repOrd_iM(isnan(acc_iM)) = []; % for the reconstruction only trials, the value in content and position reports will be NAN
                if iM == 1    % content report
                    repOrd_iM(repOrd_iM == 1) = 2; % second report after the encoding stage
                    repOrd_iM(repOrd_iM == 0) = 1; % first report after the encoding stage
                elseif iM == 2 % position report
                    repOrd_iM(repOrd_iM == 1) = 1; % first report after the encoding stage
                    repOrd_iM(repOrd_iM == 0) = 2; % second report after the encoding stage
                elseif iM == 3 % reconstruction report
                    repOrd_iM(reconsOnly == 0) = 3; % third report after the encoding stage: (1) con+pos+recons; (2) pos+con+recons;
                    repOrd_iM(reconsOnly == 1) = 1; % first report after the encoding stage
                end
                acc_iM(isnan(acc_iM)) = [];
                acc_col   = [acc_col; acc_iM];
                trlCor_iM = acc_iM;
                trlCor_iM(trlCor_iM ~= 1) = 0;
                trlCor_col = [trlCor_col; trlCor_iM];
                meas_col   = [meas_col; repmat(iM - 1, length(acc_iM), 1)];
                trlNo_col  = [trlNo_col; trlNo_iM];
                repOrd_col = [repOrd_col; repOrd_iM];
            end
            group_Col = repmat(iGrp - 1, length(acc_col), 1);
            subj_Col  = repmat(iSub - 1, length(acc_col), 1);
            acc_change_conds_subj = [acc_change_conds_subj; group_Col, subj_Col, trlNo_col, meas_col, acc_col, trlCor_col, repOrd_col];
        end

        %% ----------Transition accuracy and binding accuracy----------
        conTrue_threeRep  = conTrue_col(reconsOnly == 0);
        conRep_threeRep   = conRep_col(reconsOnly == 0);
        locTrue_threeRep  = locTrue_col(reconsOnly == 0);
        locRep_threeRep   = locRep_col(reconsOnly == 0);
        bothRep_reconsRep = bothRep_col(1 : nEpi); % 2 row: content and position; only the learning-reports trials
        % organize the content report according to the ground truth
        conRep_threeRep_sort = cell(size(conTrue_threeRep, 1), 1);
        for iEpi = 1 : size(conTrue_threeRep, 1)
            conTrue_iEpi = conTrue_threeRep{iEpi};
            [~, I] = sort(conTrue_iEpi, 'ascend');
            conRep_iEpi  = conRep_threeRep{iEpi};
            conRep_iEpi_sort = conRep_iEpi(I);
            conRep_threeRep_sort{iEpi} = conRep_iEpi_sort(1 : nTrans);
        end
        % convert the cell to matrix
        conRep_threeRep_mat = cell2mat(conRep_threeRep_sort); % 5 columns
        locRep_threeRep_mat = cell2mat(locRep_threeRep);
        locRep_threeRep_mat = locRep_threeRep_mat(:, 1 : nTrans); % 5 columns
        conRep_reconRep_mat = nan(nEpi, nTrans); % 5 columns
        locRep_reconRep_mat = nan(nEpi, nTrans);
        for iEpi = 1 : nEpi
            conRep_reconRep_mat(iEpi, :) = bothRep_reconsRep{iEpi}(1, 1 : nTrans);
            locRep_reconRep_mat(iEpi, :) = bothRep_reconsRep{iEpi}(2, 1 : nTrans);
        end
        slotCorr_marg_con = nan(sum(~reconsOnly), nTrans); % 1: the response is correct; 0: incorrect response;
        slotCorr_marg_pos = nan(sum(~reconsOnly), nTrans);
        slotCorr_join_con = nan(nEpi, nTrans);
        slotCorr_join_pos = nan(nEpi, nTrans);
        for iTrans = 1 : (nTrans)
            % ----Content reports in the marginal reports trials----
            slotCorr_marg_con(:, iTrans) = (conRep_threeRep_mat(:, iTrans) == iTrans);
            % ----Position reports in the marginal reports trials----
            slotCorr_marg_pos(:, iTrans) = (locRep_threeRep_mat(:, iTrans) == iTrans);

            % ----Content reports in the reconstruction reports trials----
            slotCorr_join_con(:, iTrans) = (conRep_reconRep_mat(:, iTrans) == iTrans);
            % ----Position reports in the reconstruction reports
            % trials----
            slotCorr_join_pos(:, iTrans) = (locRep_reconRep_mat(:, iTrans) == iTrans);
        end

        %% ----------Transition evidence----------
        % ------Transitions in the marginal reports of the 3-reports trials------
        % ------as well as in the reconstruction reports------
        % categorize the four responses, given the transition is from X to
        % Y: (1) correct-X to Y; (2)incorrect-X to nonY; (3) incorrect-nonX
        % to Y; (4) incorrect-nonX to nonY
        for ij = 1 : 2 % marginal and reconstruction
            if ij == 1
                slotCorr_con = slotCorr_marg_con;
                slotCorr_pos = slotCorr_marg_pos;
                transError_con = nan(sum(~reconsOnly), (nTrans - 1));
                transError_pos = nan(sum(~reconsOnly), (nTrans - 1));
                trlLen_temp = sum(~reconsOnly);
            elseif ij == 2
                slotCorr_con = slotCorr_join_con;
                slotCorr_pos = slotCorr_join_pos;
                transError_con = nan(nEpi, (nTrans - 1));
                transError_pos = nan(nEpi, (nTrans - 1));
                trlLen_temp = nEpi;
            end
            for iTrans = 1 : (nTrans - 1)
                % ----Content----
                iTrans_from_con = slotCorr_con(:, iTrans);
                iTrans_to_con   = slotCorr_con(:, iTrans + 1);
                transError_con(find(iTrans_from_con == 1 & iTrans_to_con == 1), iTrans) = 1; % Pr(toY=1 | fromX=1); fully correct
                transError_con(find(iTrans_from_con == 1 & iTrans_to_con == 0), iTrans) = 2; % Pr(toY=0 | fromX=1)
                transError_con(find(iTrans_from_con == 0 & iTrans_to_con == 1), iTrans) = 3; % Pr(toY=1 | fromX=0)
                transError_con(find(iTrans_from_con == 0 & iTrans_to_con == 0), iTrans) = 4; % Pr(toY=0 | fromX=0); fully incorrect

                % ----Position----
                iTrans_from_pos = slotCorr_pos(:, iTrans);
                iTrans_to_pos   = slotCorr_pos(:, iTrans + 1);
                transError_pos(find(iTrans_from_pos == 1 & iTrans_to_pos == 1), iTrans) = 1;
                transError_pos(find(iTrans_from_pos == 1 & iTrans_to_pos == 0), iTrans) = 2;
                transError_pos(find(iTrans_from_pos == 0 & iTrans_to_pos == 1), iTrans) = 3;
                transError_pos(find(iTrans_from_pos == 0 & iTrans_to_pos == 0), iTrans) = 4;

            end
            %%% ------ Quantify the transition accuracy 1) previous report was correct, 2) previous repprt was incorrect------
            transError_con_col = reshape(transError_con, [trlLen_temp * (nTrans - 1), 1]);
            transError_loc_col = reshape(transError_pos, [trlLen_temp * (nTrans - 1), 1]);
            error_con_len = [length(find(transError_con_col == 1)), length(find(transError_con_col == 2)), ...
                             length(find(transError_con_col == 3)), length(find(transError_con_col == 4))];

            error_pos_len = [length(find(transError_loc_col == 1)), length(find(transError_loc_col == 2)), ...
                             length(find(transError_loc_col == 3)), length(find(transError_loc_col == 4))];
            %%% ------ First and second half per retrieval test ------
            transError_con_firstHalf_col = reshape(transError_con(:, 1:2), [trlLen_temp * 2, 1]);
            transError_loc_firstHalf_col = reshape(transError_pos(:, 1:2), [trlLen_temp * 2, 1]);
            error_con_firstThree_len = [length(find(transError_con_firstHalf_col == 1)), length(find(transError_con_firstHalf_col == 2)), ...
                                        length(find(transError_con_firstHalf_col == 3)), length(find(transError_con_firstHalf_col == 4))];
            error_pos_firstThree_len = [length(find(transError_loc_firstHalf_col == 1)), length(find(transError_loc_firstHalf_col == 2)), ...
                                        length(find(transError_loc_firstHalf_col == 3)), length(find(transError_loc_firstHalf_col == 4))];
            transError_con_secondHalf_col = reshape(transError_con(:, 3:4), [trlLen_temp * 2, 1]);
            transError_loc_secondHalf_col = reshape(transError_pos(:, 3:4), [trlLen_temp * 2, 1]);
            error_con_secondThree_len = [length(find(transError_con_secondHalf_col == 1)), length(find(transError_con_secondHalf_col == 2)), ...
                                         length(find(transError_con_secondHalf_col == 3)), length(find(transError_con_secondHalf_col == 4))];
            error_pos_secondThree_len = [length(find(transError_loc_secondHalf_col == 1)), length(find(transError_loc_secondHalf_col == 2)), ...
                                         length(find(transError_loc_secondHalf_col == 3)), length(find(transError_loc_secondHalf_col == 4))];

            % ---- proportion of correct responses if previous item
            % is 1) correctly or 2) incorrectly reported ----
            % transAcc_count_marg_subj = nan(subLen, 4, 2); % 4: 4 different counts; 2: item and location
            % transAcc_count_join_subj = nan(subLen, 4, 2);
            if ij == 1
                % ---- Item ----
                transAcc_count_marg_subj(iSub, :, 1) = error_con_len; % 4: 4 different counts; 2: item and location
                transAcc_count_marg_firstHalf_subj(iSub, :, 1)  = error_con_firstThree_len;
                transAcc_count_marg_secondHalf_subj(iSub, :, 1) = error_con_secondThree_len;

                % ---- Location ----
                transAcc_count_marg_subj(iSub, :, 2) = error_pos_len;
                transAcc_count_marg_firstHalf_subj(iSub, :, 2)  = error_pos_firstThree_len;
                transAcc_count_marg_secondHalf_subj(iSub, :, 2) = error_pos_secondThree_len;

            elseif ij == 2
                % ---- Item ----
                transAcc_count_join_subj(iSub, :, 1) = error_con_len;
                transAcc_count_join_firstHalf_subj(iSub, :, 1)  = error_con_firstThree_len;
                transAcc_count_join_secondHalf_subj(iSub, :, 1) = error_con_secondThree_len;

                % ---- Location ----
                transAcc_count_join_subj(iSub, :, 2) = error_pos_len;
                transAcc_count_join_firstHalf_subj(iSub, :, 2)  = error_pos_firstThree_len;
                transAcc_count_join_secondHalf_subj(iSub, :, 2) = error_pos_secondThree_len;
            end
        end

        %% ----------Binding evidence----------
        % Comparison between 1) the previous transition is correct, and
        % current other dimension is correct vs. 2) the previous transition
        % is correct, and current other dimension is incorrect
        % Hypothesis: if participants real use the binding information,
        % then there should be significant difference between the two kinds
        % of trials
        for ij = 1 : 2 % marginal and reconstruction
            if ij == 1
                slotCorr_con = slotCorr_marg_con;
                slotCorr_pos = slotCorr_marg_pos;
            elseif ij == 2
                slotCorr_con = slotCorr_join_con;
                slotCorr_pos = slotCorr_join_pos;
            end
            binds_conPctr_counts = cell(2, 4); % 2: (Item|Pos) and (Pos|Item)
            for ijT = 1 : size(slotCorr_pos, 1)
                slotCorr_con_ij = slotCorr_con(ijT, :);
                slotCorr_pos_ij = slotCorr_pos(ijT, :);
                % +++++++++++++++++ (Item|Pos) +++++++++++++++++
                for iTr = 2 : nTrans
                    if slotCorr_con_ij(iTr - 1) == 1 && slotCorr_pos_ij(iTr) == 1
                        binds_conPctr_counts{1, 1} = [binds_conPctr_counts{1, 1}; slotCorr_con_ij(iTr)];

                    elseif slotCorr_con_ij(iTr - 1) == 1 && slotCorr_pos_ij(iTr) == 0
                        binds_conPctr_counts{1, 2} = [binds_conPctr_counts{1, 2}; slotCorr_con_ij(iTr)];

                    elseif slotCorr_con_ij(iTr - 1) == 0 && slotCorr_pos_ij(iTr) == 1
                        binds_conPctr_counts{1, 3} = [binds_conPctr_counts{1, 3}; slotCorr_con_ij(iTr)];

                    elseif slotCorr_con_ij(iTr - 1) == 0 && slotCorr_pos_ij(iTr) == 0
                        binds_conPctr_counts{1, 4} = [binds_conPctr_counts{1, 4}; slotCorr_con_ij(iTr)];

                    end
                end
                % +++++++++++++++++ (Pos|Item) +++++++++++++++++
                for iTr = 2 : nTrans
                    if slotCorr_pos_ij(iTr - 1) == 1 && slotCorr_con_ij(iTr) == 1
                        binds_conPctr_counts{2, 1} = [binds_conPctr_counts{2, 1}; slotCorr_pos_ij(iTr)];

                    elseif slotCorr_pos_ij(iTr - 1) == 1 && slotCorr_con_ij(iTr) == 0
                        binds_conPctr_counts{2, 2} = [binds_conPctr_counts{2, 2}; slotCorr_pos_ij(iTr)];

                    elseif slotCorr_pos_ij(iTr - 1) == 0 && slotCorr_con_ij(iTr) == 1
                        binds_conPctr_counts{2, 3} = [binds_conPctr_counts{2, 3}; slotCorr_pos_ij(iTr)];

                    elseif slotCorr_pos_ij(iTr - 1) == 0 && slotCorr_con_ij(iTr) == 0
                        binds_conPctr_counts{2, 4} = [binds_conPctr_counts{2, 4}; slotCorr_pos_ij(iTr)];
                    end
                end
            end
            % ------Calculate the proportion of correctness------
            if ij == 1
                for iC = 1 : 4
                    % ------Proportion------
                    % ****** the denominator is the total data points in each catetory ******
                    binds_conPctr_marg_subj(iSub, iC, 1) = length(find(binds_conPctr_counts{1, iC} == 1)) ./ length(binds_conPctr_counts{1, iC});
                    binds_conPctr_marg_subj(iSub, iC, 2) = length(find(binds_conPctr_counts{2, iC} == 1)) ./ length(binds_conPctr_counts{2, iC});
                    % ------Trial numbers------
                    binds_conPctr_marg_subj(iSub, iC, 3) = length(find(binds_conPctr_counts{1, iC} == 1)) ./ length(cell2mat((binds_conPctr_counts(1, :))'));
                    binds_conPctr_marg_subj(iSub, iC, 4) = length(find(binds_conPctr_counts{2, iC} == 1)) ./ length(cell2mat((binds_conPctr_counts(2, :))'));
                end

            elseif ij == 2
                for iC = 1 : 4
                    % ------Proportion------
                    % ****** the denominator is the total data points in each catetory ******
                    binds_conPctr_join_subj(iSub, iC, 1) = length(find(binds_conPctr_counts{1, iC} == 1)) ./ length(binds_conPctr_counts{1, iC});
                    binds_conPctr_join_subj(iSub, iC, 2) = length(find(binds_conPctr_counts{2, iC} == 1)) ./ length(binds_conPctr_counts{2, iC});
                    % ------Trial numbers------
                    binds_conPctr_join_subj(iSub, iC, 3) = length(find(binds_conPctr_counts{1, iC} == 1)) ./ length(cell2mat((binds_conPctr_counts(1, :))'));
                    binds_conPctr_join_subj(iSub, iC, 4) = length(find(binds_conPctr_counts{2, iC} == 1)) ./ length(cell2mat((binds_conPctr_counts(2, :))'));
                end
            end
        end

        %% Calculate the lure effect from the distractor
        % item and location in partial and full retrieval separately
        conLure  = nan(nEpi, 1);
        locLure  = nan(nEpi, 1);
        bothLure = zeros(nEpi, 2); % 2: item and location dimension in the full retrieval
        % ---------- note down the trial index where the lure stimulus was picked ----------
        lure_trialIdx = zeros(nEpi, 3); % 3: item, location and full retrieval
        for i = 1 : nEpi
            if reconsOnly(i) == 0 %% non reconstruction only trial
                % item retrieval
                conTrue_i   = conTrue_col{i};
                conRep_i    = conRep_col{i};
                conLure_i   = find(conTrue_i == 6);
                conRep_lure = conRep_i(conLure_i);
                conLure(i)  = (conRep_lure ~= 0 | sum(conRep_i == 6) ~= 0); % check if the distractor image was picked

                % location retrieval
                locTrue_i   = locTrue_col{i};
                locRep_i    = locRep_col{i};
                locLure_i   = find(locTrue_i == 6); % always the last one
                locRep_lure = locRep_i(locLure_i);
                locLure(i)  = (locRep_lure ~= 0 | sum(locRep_i ==6) ~= 0); % check if the distractor location was picked: though sum(locRep_i ==6) ~= 0 is a redundant check in location retrieval
            end
            % full retrieval
            bothRep_i      = bothRep_col{i}; % 2 row: content and position
            bothLure(i, 1) = (sum(bothRep_i(1, :) == 6) | bothRep_i(1, end) ~= 0);
            bothLure(i, 2) = (sum(bothRep_i(2, :) == 6) | bothRep_i(2, end) ~= 0);

        end
        % ---------- partial retrieval ----------
        FA_lure_subj(iSub, 1, 1) = nansum(conLure) / sum(~reconsOnly);
        FA_lure_subj(iSub, 2, 1) = nansum(locLure) / sum(~reconsOnly);

        % ---------- full retrieval ----------
        FA_lure_subj(iSub, 1, 2) = sum(bothLure(:, 1)) / nEpi;
        FA_lure_subj(iSub, 2, 2) = sum(bothLure(:, 2)) / nEpi;

        lure_trialIdx(:, 1) = conLure;
        lure_trialIdx(:, 2) = locLure;
        lure_trialIdx(:, 3) = (bothLure(:, 1) == 1 | bothLure(:, 2) == 1);

    end
    %%
    acc_trial_group{iGrp} = acc_trial_subj;
    acc_group{iGrp}       = acc_subj;
    rt_group{iGrp}        = rt_subj;
    acc_blc_group{iGrp}    = acc_blc_subj;
    acc_marginalDim_inJointRep_group{iGrp} = acc_marginalDim_inJointRep_sub;
    acc_subj_order_group{iGrp}   = acc_subj_order;
    acc_subj_orderUP_group{iGrp} = acc_subj_orderUp;
    trialPerc_group{iGrp}        = trialPerc_subj;
    acc_dim_group{iGrp}          = acc_dim_subj;
    acc_lmm_group          = [acc_lmm_group; acc_lmm_subj];
    choice_lmm_group       = [choice_lmm_group; choice_lmm_subj];
    acc_change_conds_group = [acc_change_conds_group; acc_change_conds_subj];

    % ----binding evidence----
    binds_conPctr_group{1, iGrp} = binds_conPctr_marg_subj;
    binds_conPctr_group{2, iGrp} = binds_conPctr_join_subj;

    % ----transition evidence----
    transAcc_count_group{1, iGrp} = transAcc_count_marg_subj;
    transAcc_count_group{2, iGrp} = transAcc_count_join_subj;
    transAcc_count_firstHalf_group{1, iGrp} = transAcc_count_marg_firstHalf_subj;
    transAcc_count_firstHalf_group{2, iGrp} = transAcc_count_join_firstHalf_subj;
    transAcc_count_secondHalf_group{1, iGrp} = transAcc_count_marg_secondHalf_subj;
    transAcc_count_secondHalf_group{2, iGrp} = transAcc_count_join_secondHalf_subj;

    % ----------False alarm from the lure stimuli----------
    FA_lure_group{iGrp} = FA_lure_subj;

    %% calculate the mean across subjects, for each report
    if iGrp == 1 % younger group
        mean_young = nan(1, 3);
        mean_young(1) = nanmean(acc_group{1}(:, 1)); % content
        mean_young(2) = nanmean(acc_group{1}(:, 2)); % position
        mean_young(3) = nanmean(acc_group{1}(:, 3)); % reconstruction
    end

    %% ----------Integrate the following measures into one variable and use it for the subsequent correlation analysis with model parameters----------
    allMeas_inOne = nan(subLen, 8);
    % (1) accuracy for itemReport
    allMeas_inOne(:, 1) = acc_subj(:, 1);
    % (2) accuracy for posReport
    allMeas_inOne(:, 2) = acc_subj(:, 2);
    % (3) accuracy for the difference between itemReport and posReport (itemReport - posReport)
    allMeas_inOne(:, 3) = acc_subj(:, 1) - acc_subj(:, 2);
    % (4) accuracy for reconstruction report
    allMeas_inOne(:, 4) = acc_subj(:, 3);
    % (5) accuracy for prediction from marginal reports
    allMeas_inOne(:, 5) = acc_subj(:, 12);
    % (6) difference between reconstruction report and prediction from
    % marginal reports
    allMeas_inOne(:, 6) = acc_subj(:, 3) - acc_subj(:, 12);
    % (7) binding score for item | pos in marginal reports
    % ------calculate the cost (reversed) and benefit and then average the two------
    allMeas_inOne(:, 7) = nanmean([binds_conPctr_marg_subj(:, 1, 1) - binds_conPctr_marg_subj(:, 2, 1), ...
                                   binds_conPctr_marg_subj(:, 3, 1) - binds_conPctr_marg_subj(:, 4, 1)], 2);
    % (8) binding score for pos | item in marginal reports
    % ------calculate the cost (reversed) and benefit and then average the two------
    allMeas_inOne(:, 8) = nanmean([binds_conPctr_marg_subj(:, 1, 2) - binds_conPctr_marg_subj(:, 2, 2), ...
                                   binds_conPctr_marg_subj(:, 3, 2) - binds_conPctr_marg_subj(:, 4, 2)], 2);
    allMeas_inOne_group{iGrp} = allMeas_inOne;
end
%% save data for Figure 2C in the FlexibleBinding paper
FBdata_folder = [bhvDataDir, '/FlexibleBindingPaper-Data/']; % flexible binding data folder
save([FBdata_folder, 'Fig2C_acc_group.mat'], 'acc_group');
save([FBdata_folder, 'Fig2C_acc_subj_orderUP_group.mat'], 'acc_subj_orderUP_group');
save([FBdata_folder, 'Fig2D_Exp2_FA_lure_group.mat'], 'FA_lure_group');
save([FBdata_folder, 'Fig2GHI_Exp2_binds_conPctr_group.mat'], 'binds_conPctr_group');

%% Data for Figure 2E: plotting the evidence for transition learning –– proportion of correct responses when the previous item/loc correct versus incorrect
% added by XR @ Sep 14 2025
% ------Transition evidence------
% transAcc_count_group = cell(2, nGroup); % 2: marginal and joint
% transAcc_count_firstHalf_group  = cell(2, nGroup);
% transAcc_count_secondHalf_group = cell(2, nGroup);
% transAcc_count_marg_subj = nan(subLen, 4, 2); % 4: 4 different counts; 2: item and location
% transAcc_count_join_subj = nan(subLen, 4, 2);
transAcc_flg = 0;
if transAcc_flg == 0    % all responses per retrieval test
    transAcc_count_plot = transAcc_count_group;
elseif transAcc_flg == 1 % first half
    transAcc_count_plot = transAcc_count_firstHalf_group;
elseif transAcc_flg == 2 % second half
    transAcc_count_plot = transAcc_count_secondHalf_group;
end
transAcc_plot_group       = cell(2, nGroup); % 2: marginal and joint
transAcc_plot_mjAvg_group = cell(1, nGroup); % average across marginal and joint
transAcc_drop_mjAvg_group = cell(1, nGroup);
for iGrp = 1 : nGroup % YA and OA
    %% marginal and reconstruction report separately
    for ij = 1 : 2 % marginal and reconstruction
        transAcc_count_ij = transAcc_count_plot{ij, iGrp};
        transAcc_plot = nan(size(transAcc_count_ij, 1), 4); % 4: 1-2, accuracy for item transition; 3-4, accuracy for location transition
        for iSubj = 1 : size(transAcc_count_ij, 1)
            countLen_item_ii = squeeze(transAcc_count_ij(iSubj, :, 1));
            countLen_loc_ii  = squeeze(transAcc_count_ij(iSubj, :, 2));
            % ------ Item ------
            if (countLen_item_ii(1) + countLen_item_ii(2)) ~= 0
                transAcc_plot(iSubj, 1) = countLen_item_ii(1) / (countLen_item_ii(1) + countLen_item_ii(2));
            end
            if (countLen_item_ii(3) + countLen_item_ii(4)) ~= 0
                transAcc_plot(iSubj, 2) = countLen_item_ii(3) / (countLen_item_ii(3) + countLen_item_ii(4));
            end

            % ------ Location ------
            if (countLen_loc_ii(1) + countLen_loc_ii(2)) ~= 0
                transAcc_plot(iSubj, 3) = countLen_loc_ii(1) / (countLen_loc_ii(1) + countLen_loc_ii(2));
            end
            if (countLen_loc_ii(3) + countLen_loc_ii(4)) ~= 0
                transAcc_plot(iSubj, 4) = countLen_loc_ii(3) / (countLen_loc_ii(3) + countLen_loc_ii(4));
            end
        end
        transAcc_plot_group{ij, iGrp} = transAcc_plot;
    end

    %% average across marginal and reconstruction report
    transAcc_plot_iGrp_temp = nan(size(transAcc_plot_group{1, iGrp}, 1), 4, 2); % the last 2: marginal and reconstruction
    transAcc_plot_iGrp_temp(:, :, 1) = transAcc_plot_group{1, iGrp};
    transAcc_plot_iGrp_temp(:, :, 2) = transAcc_plot_group{2, iGrp};
    
    transAcc_plot_iGrp_temp_avg = nanmean(transAcc_plot_iGrp_temp, 3);
    transAcc_plot_mjAvg_group{1, iGrp} = transAcc_plot_iGrp_temp_avg;
    %%% ------ calculate the drop percentage for both YA and OA ------
    transAcc_drop_mjAvg_group{1, iGrp} = [
        (transAcc_plot_iGrp_temp_avg(:, 1) -  transAcc_plot_iGrp_temp_avg(:, 2)) ./ transAcc_plot_iGrp_temp_avg(:, 1), ...
        (transAcc_plot_iGrp_temp_avg(:, 3) -  transAcc_plot_iGrp_temp_avg(:, 4)) ./ transAcc_plot_iGrp_temp_avg(:, 3)
        ];
end
%% save data for Figure 2E in the FlexibleBinding paper
FBdata_folder = [bhvDataDir, '/FlexibleBindingPaper-Data/']; % flexible binding data folder
save([FBdata_folder, 'Fig2E_Exp2_transAcc_plot_mjAvg_group.mat'], 'transAcc_plot_mjAvg_group');

%% color settings
colorSets = [0.98, 0.72, 0.69; ...
             0.97, 0.85, 0.67; ...
             0.33, 0.73, 0.83; ...
             0.72, 0.80, 0.88; ...
             0.54, 0.67, 0.20; ...
             0.82, 0.92, 0.78; ...
             0.78, 0.50, 0.75; ...
             0.86, 0.80, 0.89; ...
             0.75, 0.56, 0; ...
             0.40, 0.40, 0.40];

iGrp = 1;
color_Grp = colorSets([iGrp+1, iGrp+3, iGrp+5], :);

colorGrp = [230, 85, 13; ...
            253, 141, 60; ...
            253, 190, 133; ...%% content report: 3 groups
            49, 130, 189; ...
            107, 174, 214; ...
            189, 215, 231; ...%% position report: 3 groups
            117, 107, 177; ...
            158, 154, 200; ...
            203, 201, 226; ...%% reconstruction report: 3 gruops
            49, 163, 84; ...
            116, 196, 118; ...
            186, 228, 179] ./ 255; %% recons post-test: 3 groups

%% Figure 2C: Interleaved condiiton in Experiment 2
%%% seqMemTask_v2_anal_early_late_window_Binding.m
%%% ---------Figure 2C in the behavioral manuscript: Memory gradient only for the Interleaved condition---------
nConds_comp = 2; % 1) three-retrieval trials; 2) only the trials in the first reporting window
acc_group_comp = cell(2, nConds_comp);
acc_recon_pred = cell(2, nConds_comp);
for iA = 1 : 2 % YA and OA group
    if iA == 1
        groupName = 'younger';
    elseif iA == 2
        groupName = 'older';
    end
    %%% ------Three-retrieval-trials------
    acc_group_comp{iA, 1} = acc_group{iA}(:, [9,10,11]); % subLen * 3 reports
    %%% ------Data from the 1st reporting window from the Interleaved condition------
    acc_group_comp{iA, 2} = squeeze(acc_subj_orderUP_group{iA}(:, 1, :)); % subLen * 3 reports (content-first, position-first and recons-only)
end

%% plotting
figKey = 1;
if figKey == 0
    barLineWid = 2;
    errLineWid = 3;
    refLineWid = 1;
elseif figKey == 1
    barLineWid = 1;
    errLineWid = 2;
    refLineWid = 0.5;
end
figure('Position', [100 100 300 150]), clf;
barPos = [1, 1.5, 2;...  % YA
          1.15, 1.65, 2.15; ...
          3, 3.5, 4; ...
          3.15, 3.65, 4.15]; % OA
color_gray = [0, 0, 0; ...
              82, 82, 82] ./ 255;
LineStyle_conds = {'-', ':'};
color_iGrp = 0.5 * colorGrp([1, 4, 7], :) + 0.5 * [1, 1, 1];

p_recons_nGroup  = nan(1, 2);
for iA = 1 : 2 % YA and OA group
    if iA == 1
        disp('==========YA==========')
    elseif iA == 2
        disp('==========OA==========')
    end
    barPos_iA = barPos((iA - 1) * 2 + 1 : iA * 2, :);
    acc_avg_group = zeros(nConds_comp, 3);
    for iConds = 1 : nConds_comp % % 1) three reports trials; 2) only the trials in the first reporting window
        acc_iGrp = acc_group_comp{iA, iConds}(:, 1 : 3);
        [acc_avg, acc_sem] = Mean_and_Se(acc_iGrp, 1);
        acc_avg_group(iConds, :) = acc_avg;
        if iConds == 1     % plot it as barplot
            for iC = 1 : 3 % content, position and both-pre
                barPos_i = barPos_iA(iConds, iC);
                b = bar(barPos_i, acc_avg(iC), 0.45, 'LineStyle', '-', 'LineWidth', barLineWid); hold on;
                b.FaceColor = color_iGrp(iC, :); %[1, 1, 1];
                %b.EdgeColor = color_iGrp(iC, :);
                errorbar(barPos_i, acc_avg(iC), acc_sem(iC), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
                plot(barPos_i, acc_avg(iC), 'Marker', '.', 'MarkerSize', 15, 'Color', 'k'); hold on;
            end

        elseif iConds == 2 % plot it as lineplot
            for iC = 3 % only plot the both-pre; 1 : 3 % content, position and both-pre
                errorbar(barPos_iA(iConds, iC), acc_avg(iC), acc_sem(iC), 'Color', color_gray(iConds, :), 'LineStyle', LineStyle_conds{iGrp}, 'LineWidth', errLineWid); hold on;
                plot(barPos_iA(iConds, iC), acc_avg_group(iConds, iC), 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', [1, 1, 1], 'MarkerEdgeColor', color_iGrp(iC, :), 'LineStyle', '-', 'LineWidth', 1); hold on;
            end
        end
    end

    % ------Reconstruction report from the 3-reports trials vs. recons-only
    % trials------
    disp('------Recons in the three-reports trials vs. recons-only trials------')
    [h, p, ci, stats] = ttest(acc_group_comp{iA, 1}(:, 3), acc_group_comp{iA, 2}(:, 3), 'tail', 'both')
    p_recons_nGroup(iA) = p;
end

disp('^^^^^^^^^^ Bonferroni-Holm corrected p: recons in the three-reports trials vs. recons-only trials- ^^^^^^^^^^ ')
[cor_p, h] = bonf_holm(p_recons_nGroup, 0.05)

xlim([0.7, 4.7]);
ylim([0, 1]);
if figKey == 0
    % ------For presentation------
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', '', 'XTickLabel', '');
    set(gca, 'YTick', [0.5, 1], 'YTickLabel', [0.5, 1]); 
elseif figKey == 1
    % ------For Adobe Illustrator------
    set(gca, 'LineWidth', 0.8);
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', '', 'XTickLabel', '');
    set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', {'', '', ''});
end
box off;

%% ---------- Figure 2D in Experiment 2 (replication in Experiment 1); lure effect ----------
figKey = 1;  % 0: figure for presentation; 1: figure for AI.
if figKey == 0
    barLineWid = 2;
    errLineWid = 3;
    refLineWid = 1;
elseif figKey == 1
    barLineWid = 1;
    errLineWid = 1.5; %2;
    refLineWid = 0.5;
end
barPos = [1, 1.7, 2.0, 2.7];
figure('Position', [100 100 160 120]), clf;
for iGrp = 1 : nGroup % YA and OA
    FA_lure_iGrp = FA_lure_group{iGrp}; % subLen * 2(item/loc) * 2(partial/full)
    FA_lure_iGrp = nanmean(FA_lure_iGrp, 3);
    if iGrp == 1
        disp('------YA------')
    elseif iGrp == 2
        disp('------OA------')
    end
    iDm_idx = (iGrp - 1) * 2 + 1 : iGrp * 2;
    [tAcc_avg, tAcc_sem] = Mean_and_Se(FA_lure_iGrp, 1);
    barPos_i = barPos(iDm_idx);
    plot(barPos_i, FA_lure_iGrp, 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', 0.4); hold on;
    plot(barPos_i, tAcc_avg, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
 
    for iDm = 1 : 2 % item or location dimension
        errorbar(barPos_i(iDm), tAcc_avg(iDm), tAcc_sem(iDm), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
        plot(barPos_i(iDm), tAcc_avg(iDm), 'Marker', 'o', 'MarkerSize', 4.5, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', color_Grp(iDm, :), 'LineStyle', '-'); hold on;
    end
    % ------ Statistical tests ------
    disp('======== FA: item vs. location ========');
    [h, p, ci, stats] = ttest(FA_lure_iGrp(:, 1), FA_lure_iGrp(:, 2))
end
xlim([0.6, 3.1]);
ylim([0, 1]);
if figKey == 0
    % ------For presentation------
    plot(xlim, [1/5, 1/5], 'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 1); hold on;
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', '', 'XTickLabel', '');
    set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', 0 : 0.5 : 1);
elseif figKey == 1
    % ------For Adobe Illustrator------
    %plot(xlim, [1/5, 1/5], 'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', refLineWid); hold on;
    set(gca, 'LineWidth', 0.6); % 0.8
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', '', 'XTickLabel', '');
    set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', {'', '', ''});
end
box off;

%% ---------- Figure 2E in Experiment 2 (replication in Experiment 1); transition accurcy ----------
figKey = 1;  % 0: figure for presentation; 1: figure for AI.
if figKey == 0
    barLineWid = 2;
    errLineWid = 3;
    refLineWid = 1;
elseif figKey == 1
    barLineWid = 1;
    errLineWid = 1.5; %2;
    refLineWid = 0.5;
end
barPos = [1, 1.7, 2.0, 2.7];
condsWord = 'marginal+full';
for iGrp = 1 : nGroup % YA and OA
    if iGrp == 1
        disp('------YA------')
    elseif iGrp == 2
        disp('------OA------')
    end
    figure('Position', [100 100 150 120]), clf;
    for iDm = 1 : 2 % item or location dimension
        if iDm == 1
            color_iCp = [color_Grp(1, :); 1, 1, 1]; % facecolor
            color_iEd = [0, 0, 0; 0, 0, 0]; % edgecolor
            dimWords  = 'item';
        elseif iDm == 2
            color_iCp = [color_Grp(2, :); 1, 1, 1];
            color_iEd = [0, 0, 0; 0, 0, 0];
            dimWords  = 'location';
        end
        iDm_idx = (iDm - 1) * 2 + 1 : iDm * 2;
        transAcc_plot_ij     = transAcc_plot_mjAvg_group{1, iGrp}(:, iDm_idx); % subj * 4
        [tAcc_avg, tAcc_sem] = Mean_and_Se(transAcc_plot_ij, 1);
        barPos_i = barPos(iDm_idx);
        plot(barPos_i, transAcc_plot_ij, 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', 0.4); hold on;
        plot(barPos_i, tAcc_avg, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
        for jDm = 1 : 2 % 2: previous report is correct or incorrect
            if jDm == 1
                preWords = 'preCorrect';
            elseif jDm == 2
                preWords = 'preIncorrect';
            end
            errorbar(barPos_i(jDm), tAcc_avg(jDm), tAcc_sem(jDm), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
            plot(barPos_i(jDm), tAcc_avg(jDm), 'Marker', 'o', 'MarkerSize', 4.5, 'MarkerEdgeColor', color_iEd(jDm, :), 'MarkerFaceColor', color_iCp(jDm, :), 'LineStyle', '-'); hold on;
            % ------ Statistical tests ------
            disp(['======== ', condsWord, '-', dimWords, '-', preWords, '========']);
            [h, p, ci, stats] = ttest(transAcc_plot_ij(:, jDm))
        end
        for jDm = 1 : 2 % 2: previous report is correct or incorrect
            errorbar(barPos_i(jDm), tAcc_avg(jDm), tAcc_sem(jDm), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
            plot(barPos_i(jDm), tAcc_avg(jDm), 'Marker', 'o', 'MarkerSize', 4.5, 'MarkerEdgeColor', color_iEd(jDm, :), 'MarkerFaceColor', color_iCp(jDm, :), 'LineStyle', '-'); hold on;
        end
    end
    xlim([0.6, 3.1]);
    ylim([0, 1]);
    if figKey == 0
        % ------For presentation------
        plot(xlim, [1/5, 1/5], 'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 1); hold on;
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', '', 'XTickLabel', '');
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', 0 : 0.5 : 1);
    elseif figKey == 1
        % ------For Adobe Illustrator------
        %plot(xlim, [1/5, 1/5], 'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', refLineWid); hold on;
        set(gca, 'LineWidth', 0.6); % 0.8
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', '', 'XTickLabel', '');
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', {'', '', ''});
    end
    box off;
end

%% ---------- Figure 2G in Experiment 2 (replication in Experiment 1): grand average of the binding score for the other dimension correct vs. incorrect in YA and OA ----------
figKey = 1;  % 0: figure for presentation; 1: figure for AI.
if figKey == 0
    barLineWid = 2;
    errLineWid = 3;
    refLineWid = 1;
elseif figKey == 1
    barLineWid = 1;
    errLineWid = 1.5; %2;
    refLineWid = 0.5;
end
barPos = [1, 1.7];
binds_dataAnal = binds_conPctr_group;
binds_dataAnal_grandAvg = cell(1, nGroup);
bindScore_grandAvg = cell(1, nGroup);
for iGrp = 1 : nGroup
    if iGrp == 1
        disp('------YA------')
    elseif iGrp == 2
        disp('------OA------')
    end
    %%% ------ Average across different dimension ------
    %%% ----Partial retrieval----
    binds_dataAnal_partial = nanmean(binds_dataAnal{1, iGrp}(:, :, 1 : 2), 3); % subj * 4 * 2
    binds_dataAnal_partial = [nanmean(binds_dataAnal_partial(:, [1, 3]), 2), nanmean(binds_dataAnal_partial(:, [2, 4]), 2)];
    %%% ----Full retrieval----
    binds_dataAnal_full = nanmean(binds_dataAnal{2, iGrp}(:, :, 1 : 2), 3);
    binds_dataAnal_full = [nanmean(binds_dataAnal_full(:, [1, 3]), 2), nanmean(binds_dataAnal_full(:, [2, 4]), 2)];
    %%% ----average across partial and full finally----
    binds_dataAnal_Tmp = nan(size(binds_dataAnal_partial, 1), 2, 2);
    binds_dataAnal_Tmp(:, :, 1) = binds_dataAnal_partial;
    binds_dataAnal_Tmp(:, :, 2) = binds_dataAnal_full;
    binds_dataAnal_grandAvg{1, iGrp} = nanmean(binds_dataAnal_Tmp, 3);
    bindScore_grandAvg{1, iGrp} = binds_dataAnal_grandAvg{1, iGrp}(:, 1) - binds_dataAnal_grandAvg{1, iGrp}(:, 2);

    figure('Position', [100 100 80 120]), clf;
    binds_grandAvg_iGrp = binds_dataAnal_grandAvg{1, iGrp};
    [tAcc_avg, tAcc_sem] = Mean_and_Se(binds_grandAvg_iGrp, 1);
    barPos_i = barPos;
    plot(barPos_i, binds_grandAvg_iGrp, 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', 0.4); hold on;
    plot(barPos_i, tAcc_avg, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
    for jDm = 1 : 2 % 2: previous report is correct or incorrect
        if jDm == 1
            colorFace_jDm = [0, 0, 0];
        elseif jDm == 2
            colorFace_jDm = [1, 1, 1];
        end
        errorbar(barPos_i(jDm), tAcc_avg(jDm), tAcc_sem(jDm), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
        plot(barPos_i(jDm), tAcc_avg(jDm), 'Marker', 'o', 'MarkerSize', 4.5, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorFace_jDm, 'LineStyle', '-'); hold on;
    end
    xlim([0.6, 2.1]);
    ylim([0, 1]);
    if figKey == 0
        % ------For presentation------
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', '', 'XTickLabel', '');
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', 0 : 0.5 : 1);
    elseif figKey == 1
        % ------For Adobe Illustrator------
        set(gca, 'LineWidth', 0.6); % 0.8
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', '', 'XTickLabel', '');
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', {'', '', ''});
    end
    box off;

end

%% ------ Check the binding score ------
for iGrp = 1 : nGroup % YA and OA
    if iGrp == 1
        disp('========== YA ==========');
    elseif iGrp == 2
        disp('========== OA ==========');
    end
    bindScore_grandAvg_iGrp = bindScore_grandAvg{1, iGrp};
    [bindScore_avg, bindScore_sem] = Mean_and_Se(bindScore_grandAvg_iGrp)
end

%% ---------- Figure 2H and 2I in Experiment 2 (replication in Experiment 1): simplying the effect for Partial retrieval vs. Full retrieval; Transition correct vs. incorrect ----------
% (1) Figure 2H: partial vs. full in YA and OA
% (2) Figure 2I: preceding report in the current dimension was correct vs.
% incorrect
figKey = 1;  % 0: figure for presentation; 1: figure for AI.
if figKey == 0
    barLineWid = 2;
    errLineWid = 3;
    refLineWid = 1;
elseif figKey == 1
    barLineWid = 1;
    errLineWid = 1.5; %2;
    refLineWid = 0.5;
end
binds_dataAnal = binds_conPctr_group;
bindsScore_PartialFull    = cell(1, nGroup);
bindsScore_PreceedingCorr = cell(1, nGroup);
for iGrp = 1 : nGroup % YA and OA
    %%% ------ Partial vs. Full ------
    %%% (1) Partial
    binds_dataAnal_partial = binds_dataAnal{1, iGrp}(:, :, 1 : 2);
    binds_dataAnal_partial_diff = nan(size(binds_dataAnal_partial , 1), 2, 2); % 1st 2: cost and benefit; 2nd 2: item | pos and pos | item
    binds_dataAnal_partial_diff(:, :, 1) = [binds_dataAnal_partial(:, 1, 1) - binds_dataAnal_partial(:, 2, 1), ...
                                            binds_dataAnal_partial(:, 3, 1) - binds_dataAnal_partial(:, 4, 1)]; % item | pos
    binds_dataAnal_partial_diff(:, :, 2) = [binds_dataAnal_partial(:, 1, 2) - binds_dataAnal_partial(:, 2, 2), ...
                                            binds_dataAnal_partial(:, 3, 2) - binds_dataAnal_partial(:, 4, 2)]; % pos | item
    binds_dataAnal_partial_diffAvg = nanmean(nanmean(binds_dataAnal_partial_diff, 3), 2);
    %%% (2) Full
    binds_dataAnal_full = binds_dataAnal{2, iGrp}(:, :, 1 : 2);
    binds_dataAnal_full_diff = nan(size(binds_dataAnal_full , 1), 2, 2);
    binds_dataAnal_full_diff(:, :, 1) = [binds_dataAnal_full(:, 1, 1) - binds_dataAnal_full(:, 2, 1), ...
                                         binds_dataAnal_full(:, 3, 1) - binds_dataAnal_full(:, 4, 1)];
    binds_dataAnal_full_diff(:, :, 2) = [binds_dataAnal_full(:, 1, 2) - binds_dataAnal_full(:, 2, 2), ...
                                         binds_dataAnal_full(:, 3, 2) - binds_dataAnal_full(:, 4, 2)];
    binds_dataAnal_full_diffAvg = nanmean(nanmean(binds_dataAnal_full_diff, 3), 2);
    bindsScore_PartialFull{1, iGrp} = [binds_dataAnal_partial_diffAvg, binds_dataAnal_full_diffAvg];

    %%% ------ Preceding report in the current dimension was correct vs.
    %%% incorrect ------
    binds_dataAnal_partial_CostBenefit = nanmean(binds_dataAnal_partial_diff, 3);
    binds_dataAnal_full_CostBenefit    = nanmean(binds_dataAnal_full_diff, 3);
    binds_dataAnal_CostBenefitAll = nan(size(binds_dataAnal_partial_CostBenefit, 1), 2, 2); % 1st 2: cost and benefit; 2nd 2: partial and full report
    binds_dataAnal_CostBenefitAll(:, :, 1) = binds_dataAnal_partial_CostBenefit;
    binds_dataAnal_CostBenefitAll(:, :, 2) = binds_dataAnal_full_CostBenefit;
    binds_dataAnal_CostBenefitAvg = nanmean(binds_dataAnal_CostBenefitAll, 3);
    bindsScore_PreceedingCorr{1, iGrp} = binds_dataAnal_CostBenefitAvg;
end
barPos = [1, 1.7];
colorFace_jDm = [0, 0, 0];
p_values_allConds = [];
for iPlot = 1 : 2 % 1-partial vs. full; 2-preceeding correct (cost) vs. incorrect (benefit)
    if iPlot == 1
        dataPlot = bindsScore_PartialFull;
    elseif iPlot == 2
        dataPlot = bindsScore_PreceedingCorr;
    end
    for iGrp = 1 : nGroup % YA and OA
        if iGrp == 1
            disp('------YA------')
        elseif iGrp == 2
            disp('------OA------')
        end
        figure('Position', [100 100 80 120]), clf;
        binds_grandAvg_iGrp = dataPlot{1, iGrp};
        [tAcc_avg, tAcc_sem] = Mean_and_Se(binds_grandAvg_iGrp, 1);
        barPos_i = barPos;
        plot(barPos_i, binds_grandAvg_iGrp, 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', 0.4); hold on;
        plot(barPos_i, tAcc_avg, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
        for jDm = 1 : 2 
            errorbar(barPos_i(jDm), tAcc_avg(jDm), tAcc_sem(jDm), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
            plot(barPos_i(jDm), tAcc_avg(jDm), 'Marker', 'o', 'MarkerSize', 4.5, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorFace_jDm, 'LineStyle', '-'); hold on;
            % -------- Statistical tests --------
            if jDm == 1
                disp('------Partial or Correct------')
            elseif jDm == 2
                disp('------Full or Incorrect------')
            end
            data_jDm = binds_grandAvg_iGrp(:, jDm);
            data_jDm(isnan(data_jDm)) = [];
            [h, p, ci, stats] = ttest(data_jDm)
            p_values_allConds = [p_values_allConds; p];
        end
        xlim([0.6, 2.1]);
        ylim([-0.2, 1]);
        if figKey == 0
            % ------For presentation------
            set(gca, 'LineWidth', 2);
            set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
            set(gca, 'XTick', '', 'XTickLabel', '');
            set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', 0 : 0.5 : 1);
        elseif figKey == 1
            % ------For Adobe Illustrator------
            set(gca, 'LineWidth', 0.6); % 0.8
            set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
            set(gca, 'XTick', '', 'XTickLabel', '');
            set(gca, 'YTick', [-0.2, 0 : 0.5 : 1], 'YTickLabel', {'', '', '', ''});
        end
        box off;
    end
end
%%
disp('^^^^^^^^^^ Bonferroni-Holm corrected p: only the Figure 2H ^^^^^^^^^^ ')
[cor_p, h] = bonf_holm(p_values_allConds(5 : end), 0.05)
