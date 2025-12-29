% seqMemTask_BehvPaper_Figure2.m
% revision by XR @ Dec 21 2025
% simplified based on seqMemTask_v1_anal_summary.m
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
folder      = '/Users/ren/Projects-NeuroCode/MyExperiment/Aging-SeqMemTask'; %'/Volumes/HierarchicalCluster';
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

expList = {'pathA', 'pathB', 'pathC', 'all'};
exp_id  = 4;
exp_tgt = expList{exp_id};
subjList_young = {'5ed8f768fd2b9c0705700c30', '03-08-2023_17h40.53.123'; '5f3e52355f3b5d02fe4dcf61', '2023-03-08_23h21.38.250'; '5f761e5106b786071f45b4aa', '03-09-2023_11h26.52.637'; ...
                  '5f915119067a16025ca4327c', '03-08-2023_18h33.19.422'; '601ee8728c09203d5489c743', '03-08-2023_11h22.14.281'; '605b52d6903ac0311ac2f941', '03-08-2023_18h33.14.900'; ...
                  '606f25c909b8c8c57c904dcc', '03-08-2023_15h26.09.357'; '6103872a37b1dd5b8f350237', '03-08-2023_11h32.33.617'; '6154b3a7e4bc93e1ce93c6bc', '2023-03-20_12h04.53.471'; ...
                  '5f444f41d9255617e1ba10ff', '2023-03-20_19h13.30.426'; '602e50ffe527c397a3914569', '2023-03-20_19h02.37.911'; '5f971b07a3ef980924b8702c', '2023-03-20_20h06.27.929'; ...
                  '614c797f4975041c463df13e', '2023-03-20_20h14.17.616'; '60d1a4da9e47f5d3b3e85a70', '2023-03-20_20h30.32.298'; '60d3867b6e7529fc77032b54', '2023-03-20_21h29.29.301'; ...
                  '5e9bf93992a7fb1111847755', '2023-03-08_18h33.20.911'; '5e9d9dd04ca6690944caf284', '03-08-2023_19h41.57.572'; '5f591a1be30bf81170acf5b9', '03-08-2023_17h34.21.699'; ...
                  '5fd663b8b06c1b570569dc2b', '03-08-2023_17h56.53.082'; '60fda80455b1001936d34ce8', '03-08-2023_19h38.28.682'; '616d6215bbc05a1429fcfc8f', '03-08-2023_19h35.11.383'; ...
                  '60185218bc28d50a56653c58', '03-08-2023_18h33.09.255'; '5e92267dcc3d6c553eb8c817', '2023-03-20_19h04.52.279'; '5ed7a7a467a98224295459ff', '03-20-2023_18h01.21.913'; ...
                  '5ed33b5cea44492cf6a28cdf', '03-20-2023_12h04.49.847'; '613ceea33c37d70e87c5232e', '2023-03-20_18h05.10.662'; '5fb45818d84de3123a4304f8', '2023-03-20_21h02.42.051'; ...
                  '5c0333c84e1b7b00016a83fd', '03-08-2023_17h45.41.644'; '5d575f44182f9f001900a380', '03-08-2023_17h32.57.282'; '5e760a9c865a94167c2f14ca', '03-08-2023_17h40.33.797'; ...
                  '5eb11bea8fb5a51eb2e73ecc', '03-08-2023_19h53.50.063'; '60b730d71e6cc3224b7ad99e', '03-08-2023_19h34.08.058'; '611e0d5e51e659a7bab55418', '03-08-2023_19h35.20.645'; ...
                  '61412d8b80e3b72f1612c408', '03-08-2023_11h35.41.649'; '6068354d363a73e3d5851e19', '2023-03-08_18h56.37.730'; '60ca6d332749493cf6326ab6', '2023-03-20_18h06.55.905'; ...
                  '61076ea4d11798351e01d677', '2023-03-20_12h04.56.561'; '63ced1948d33275a0d76bbdd', '2023-03-20_18h08.56.340'; '60ca1c436511b9fc8ab35615', '2023-03-20_21h06.49.702'};

subjList_old   = {'5a74c608527fd00001462f87', '03-09-2023_08h23.30.024'; '5d52e15363bdb100010883f1', '03-09-2023_05h06.54.849'; '5e2872fffcdb020467d12799', '03-08-2023_15h05.15.629'; ...
                  '5ea0b332440c6a0a2b78b2df', '03-09-2023_07h49.32.036'; '5f392fae9168cc472f0468f6', '03-08-2023_14h47.27.508'; '5fc3849ea684f631e16d2b64', '03-09-2023_04h17.58.718'; ...
                  '55abb1a5fdf99b501fab62e3', '2023-03-08_15h12.10.842'; '612f5beeadebe8b6969dc5a9', '03-09-2023_08h25.52.208'; '6165754a7219521b599335a0', '03-21-2023_06h54.45.638'; ...
                  '5d0ca909a37d8b0018913ca4', '03-21-2023_09h14.44.099'; '571164253d47e80009635a80', '2023-03-24_07h11.55.290'; ...
                  '5c1150e31e7d6900018207ef', '2023-03-08_10h35.55.529'; '63dfc2c8d2a0118a9552a5ac', '03-08-2023_18h28.20.295'; '599c1130638529000144ffdd', '03-08-2023_18h27.21.808'; ...
                  '6050acf587bc4d19589f8002', '03-08-2023_18h10.35.717'; '60375d67dde85cad87f3c99d', '03-08-2023_18h30.00.243'; '5c58244f62290e00012fa64c', '2023-03-20_18h43.49.184'; ...
                  '5ea026ebbc5a78011ea73ffb', '2023-03-20_19h38.49.000'; '5ea073cfb8640c040e830651', '2023-03-21_12h02.28.164'; '5ea1950819664f04542a50bf', '2023-03-21_12h34.44.630'; ...
                  '59bda4313c45a10001ccca80', '2023-03-21_15h33.15.927'; ...
                  '5e51c83155087a204d8e7f1f', '03-08-2023_13h58.57.484'; '59bb9e3a178f1b0001828b12', '03-08-2023_23h04.11.074'; '632651a9f1a00140dc2b9487', '03-08-2023_13h58.45.001';...
                  '6348801ddd19d0c7bfddac73', '03-08-2023_14h40.30.799'; '6394056ffb0438ad5373b6ec', '03-09-2023_14h49.28.749'; '63e3ae9a69fe37cb356d7867', '2023-03-21_13h09.12.908'};
suffixWord_young_list = {'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', ...
                         'pathB', 'pathB', 'pathB', 'pathB', 'pathB', 'pathB', 'pathB', 'pathB', 'pathB', 'pathB', 'pathB', 'pathB', ...
                         'pathC', 'pathC', 'pathC', 'pathC', 'pathC', 'pathC', 'pathC', 'pathC', 'pathC', 'pathC', 'pathC', 'pathC'};

suffixWord_old_list   = {'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA', 'pathA',...
                         'pathB', 'pathB', 'pathB', 'pathB', 'pathB', 'pathB', 'pathB', 'pathB', 'pathB', 'pathB', ...
                         'pathC', 'pathC', 'pathC', 'pathC', 'pathC', 'pathC'};

ageList_younger = [22, 21, 25, ...
                   21, 23, 20, ...
                   23, 23, 21, ...
                   29, 20, 20, ...
                   29, 23, 23, ...
                   21, 23, 23, ...
                   21, 23, 22, ...
                   22, 21, 25, ...
                   23, 25, 22, ...
                   24, 22, 21, ...
                   22, 22, 22, ...
                   30, 21, 20, ....
                   23, 21, 25];

genderList_younger = {'M', 'M', 'M', ...
                      'F', 'M', 'F', ...
                      'M', 'M', 'M', ...
                      'M', 'M', 'F', ...
                      'M', 'M', 'F', ...
                      'M', 'M', 'F', ...
                      'M', 'F', 'F', ...
                      'M', 'M', 'M', ...
                      'F', 'F', 'M', ...
                      'M', 'F', 'M', ...
                      'F', 'F', 'F', ...
                      'F', 'M', 'F', ...
                      'F', 'A', 'M'}; % agender

ageList_older    = [81, 78, 76, ...
                    76, 75, 77, ...
                    79, 76, 78, ...
                    85, 76, ...
                    77, 78, 75, ...
                    79, 75, 75, ... % the last participant this row: age unknown
                    78, 82, 82, ...
                    76, ...
                    75, 79, 78, ...
                    76, 75, 75];
genderList_older = {'F', 'M', 'M', ...
                    'M', 'M', 'M', ...
                    'M', 'M', 'M', ...
                    'M', 'F', ...
                    'F', 'F', 'F', ...
                    'F', 'F', '', ... % gender unknown
                    'M', 'F', 'F', ...
                    'M', ...
                    'F', 'M', 'F', ...
                    'F', 'F', 'M'};
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


%% data analysis: RT & accuracy
% overall accuracy & RTs across all blocks
acc_group = cell(1, nGroup);
rt_group  = cell(1, nGroup);
% time resistence analysis
acc_subj_orderUP_group = cell(1, nGroup);
% accuracy for each sub-block
acc_blc_group = cell(1, nGroup);
% accuracy & RTs in each trial
acc_trial_group = cell(1, nGroup);
rt_trial_group  = cell(1, nGroup);
% percentage of fully correct trials
trialPerc_group = cell(1, nGroup);

% ------accuracy for subsequent LMM analysis------
acc_lmm_group = [];
% ------ response (0/1) for subsequent GLMM analysis -----
% added by rxj @ April 12 2025
resp_glmm_group = []; % response (0/1) in each slot

% ------Binding related calculation------
binds_conPctr_group         = cell(2, nGroup); % 2: marginal and reconstruction reports
% ------Transition evidence------
transAcc_count_group            = cell(2, nGroup); % 2: marginal and joint
transAcc_count_firstHalf_group  = cell(2, nGroup);
transAcc_count_secondHalf_group = cell(2, nGroup);
% ----------Integrate the following measures into one variable and use it for the subsequent correlation analysis with model parameters----------
allMeas_inOne_group = cell(1, nGroup);
% ----------False alarm from the lure stimuli----------
FA_lure_group = cell(1, nGroup);

for iGrp = 1 : nGroup %% younger and older adults
    %%
    if iGrp == 1
        subj_list = subjList_young;
    elseif iGrp == 2
        subj_list = subjList_old;
    end
    subLen = length(subj_list);
    % ------accuracy & RTs in each trial------
    acc_trial_subj     = nan(subLen, 5, (nEpi + postTn)); % 5: (1) content report; (2) position report; (3) both report; (4) content in both report; (5) position in both report;
    rt_trial_subj      = nan(subLen, 5, (nEpi + postTn));
    % ------overall accuracy & RTs across all blocks------
    acc_subj           = nan(subLen, 10); % post-test: (6) both; (7) content in both report; (8) position in both report; (9) and (10) multiplication of content and position accuracy for marginal and reconstruction reports separately
    rt_subj            = nan(subLen, 8);
    % ------accuracy in each sub-block------
    acc_blc_subj       = nan(subLen, (nBlock+1)*2, 3); % 3: (1) content report; (2) position report; (3) both report;
    % ------percentage of the fully correct trials------
    trialPerc_subj     = nan(subLen, 3); % 3: content, position, content + position
    % ------accuracy for subsequent LMM analysis------
    acc_lmm_subj        = [];
    acc_change_lmm_subj = [];
    % ------binary response (0/1) in each slot for subvsequent GLMM analysis-------
    resp_glmm_subj      = [];
    % ------time consistency analysis------
    acc_subj_orderUp    = nan(subLen, 3, 2); % 2nd dimension: first, second, and third report order; 3rd dimension: content, position and reconstruction (for reconstruction, check if 1con-2pos differs from 1pos-1con)
    % ======There is another control condition, where we control the
    % correctness of the other dimension, but change the
    % correctness/incorrectness of the previous transition for one dimension
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
        %seqMemDir  = [bhvDataDir, subjBv, '_EpisodicMemoryTask_fourBlocks'];
        if isequal(exp_tgt, 'all')
            if iGrp == 1
                suffixWord = suffixWord_young_list{iSub};
            elseif iGrp == 2
                suffixWord = suffixWord_old_list{iSub};
            end
            subjPath   = [bhvDataDir, '/AgingReplay-v1-', suffixWord, '/'];
        end
        if isequal(suffixWord, 'pathA')
            seqMem_subj = readtable([subjPath, subjBv, '_EpisodicMemoryTask', '_', subjTm, '.csv']);
        else
            seqMem_subj = readtable([subjPath, subjBv, '_EpisodicMemoryTask-', suffixWord, '_', subjTm, '.csv']);
        end

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
        conTrue_tmp = seqMem_subj.conReportTrue;
        conTrue_tmp = conTrue_tmp(~cellfun('isempty', conTrue_tmp));
        conRep_tmp  = seqMem_subj.conReportOrd;
        conRep_tmp  = conRep_tmp(~cellfun('isempty', conRep_tmp));
        conRT_tmp   = seqMem_subj.conRTs;
        conRT_tmp   = conRT_tmp(~cellfun('isempty', conRT_tmp));
        conTrue_col = cell(nEpi, 1);
        conRep_col  = cell(nEpi, 1);
        conRT_col   = cell(nEpi, 1);
        for i =  1 : nEpi
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

        %% ---------- position report ----------
        locTrue_tmp = seqMem_subj.locReportTrue;
        locTrue_tmp = locTrue_tmp(~cellfun('isempty', locTrue_tmp));
        locRep_tmp  = seqMem_subj.locReportOrd;
        locRep_tmp  = locRep_tmp(~cellfun('isempty', locRep_tmp));
        locRT_tmp   = seqMem_subj.locRTs;
        locRT_tmp   = locRT_tmp(~cellfun('isempty', locRT_tmp));
        locTrue_col = cell(nEpi, 1);
        locRep_col  = cell(nEpi, 1);
        locRT_col   = cell(nEpi, 1);
        for i =  1 : nEpi
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
            bothRep_col{i}  = bothCol;
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
        % ------ quantify the correctness (1-correct, 0-incorrect) in each response ------
        conResp_binary = nan(nEpi * nTrans, 1);
        locResp_binary = nan(nEpi * nTrans, 1);
        for i =  1 : nEpi
            resp_idx = (i - 1) * nTrans + 1 : i * nTrans;
            % content report
            conTrue_i = conTrue_col{i};
            conRep_i  = conRep_col{i};
            conRT_i   = conRT_col{i};
            conRep_i  = conRep_i(conTrue_i ~= 6);
            conRT_i   = conRT_i(conTrue_i ~= 6);
            conTrue_i = conTrue_i(conTrue_i ~= 6);
            acc_trial_subj(iSub, 1, i)  = (sum(conRep_i == conTrue_i)) / nTrans;
            conResp_binary(resp_idx, 1) = (conRep_i == conTrue_i);
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
            acc_trial_subj(iSub, 2, i)  = (sum(locRep_i == locTrue_i)) / nTrans;
            locResp_binary(resp_idx, 1) = (locRep_i == locTrue_i);
            %%% RT calculation based on single correct item
            if sum(locRep_i == locTrue_i) ~= 0
                rt_trial_subj(iSub, 2, i) = nanmean(locRT_i(locRep_i == locTrue_i));
            end
        end

        % reconstruction
        % ------ quantify the correctness (1-correct, 0-incorrect) in each response ------
        bothResp_binary = nan(nEpi * nTrans, 1);
        for i =  1 : (nEpi + postTn)
            resp_idx = (i - 1) * nTrans + 1 : i * nTrans;
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
            if i <= nEpi
                bothResp_binary(resp_idx, 1) = acc_j;
            end
            acc_trial_subj(iSub, 3, i) = sum(acc_j) / nTrans;
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

        %% -------- correct or incorrect responses in each response of the report --------
        % added by rxj @ April 12 2025
        trialNo = 1 : 1 : nEpi;
        trialNo_mat = repmat(trialNo, nTrans, 1); % 5 responses in each report
        trialNo_mat = reshape(trialNo_mat, [length(trialNo) * nTrans, 1]);
        trialNo_mat = repmat(trialNo_mat, 1, 3); % 3 reporting types

        % trial number for each slot in each trial
        trialNo_col  = reshape(trialNo_mat, [size(trialNo_mat, 1) * size(trialNo_mat, 2), 1]);
        % transition for each report in each trial
        transTrial_col = repmat((1 : 1 : nTrans)', [nEpi * 3, 1]);
        % response types: 0-item report, 1-location report,
        % 2-reconstruction report
        respType_col = [zeros(length(trialNo) * nTrans, 1); ones(length(trialNo) * nTrans, 1); repmat(2, [length(trialNo) * nTrans, 1])];
        % label if the response is correct or incorrect per slot for each
        % response
        responseCol_binary = [conResp_binary; locResp_binary; bothResp_binary];
        % subject
        subj_col  = repmat(iSub - 1, [size(responseCol_binary, 1), 1]);
        % group
        group_col = repmat(iGrp - 1, [size(responseCol_binary, 1), 1]);
        % put all variable together
        resp_glmm_subj = [resp_glmm_subj; group_col, subj_col, trialNo_col, transTrial_col, respType_col, responseCol_binary];

        %% overall accuracy & RT
        % 64 episodes before numJdgTask
        acc_subj(iSub, 1 : 5) = squeeze(nanmean(acc_trial_subj(iSub, :, 1 : nEpi), 3));
        rt_subj(iSub, 1 : 5)  = squeeze(nanmean(rt_trial_subj(iSub, :, 1 : nEpi), 3));
        % 8 episodes after numJdgTask
        acc_subj(iSub, 6 : 8) = squeeze(nanmean(acc_trial_subj(iSub, 3 : end, (nEpi + 1) : end), 3));
        rt_subj(iSub, 6 : 8)  = squeeze(nanmean(rt_trial_subj(iSub, 3 : end, (nEpi + 1) : end), 3));

        %% calculate the predicted joint accuracy from either marginal reports or the marginal dimension in the reconstruction reports
        % added by rxj @ Sep 14 2024
        acc_subj(iSub, 9)  = nanmean(squeeze(acc_trial_subj(iSub, 1, 1 : nEpi)) .* squeeze(acc_trial_subj(iSub, 2, 1 : nEpi)));
        acc_subj(iSub, 10) = nanmean(squeeze(acc_trial_subj(iSub, 4, 1 : nEpi)) .* squeeze(acc_trial_subj(iSub, 5, 1 : nEpi)));
        margPred_iSub = squeeze(acc_trial_subj(iSub, 1, 1 : nEpi)) .* squeeze(acc_trial_subj(iSub, 2, 1 : nEpi));
        accDiff_real_pred_iSub = margPred_iSub - squeeze(acc_trial_subj(iSub, 3, 1 : nEpi)); % per subject: one vector (64 values)

        %% accuracy for content, position and reconstruction for two displaying orders
        % added by rxj @ 01/24/2024
        % plotting according to the reporting order after the encoding stage
        % 3 kinds of reporting order:
        % (1) first order: content in content-first trials, posiiton in position-first trials,
        % recon report in reconstruction only trials;
        % (2) second order: pos in content-first trials, con in position-first trials;
        % (3) third order: recon report in non recon-only trials
        %acc_subj_orderUp = nan(subLen, 3, 3); % 2nd dimension: first, second, and third report order; 3rd dimension: content, position and reconstruction
        acc_trialReport = (squeeze(acc_trial_subj(iSub, 1 : 3, 1 : nEpi)))'; % 64 trials * 3 columns
        % ------first order------
        % acc_subj_orderUp = nan(subLen, 3, 2);
        % content-first, position-first
        acc_subj_orderUp(iSub, 1, 1) = nanmean(acc_trialReport(testOrd == 0, 1));
        acc_subj_orderUp(iSub, 1, 2) = nanmean(acc_trialReport(testOrd == 1, 2));

        % ------second order------
        acc_subj_orderUp(iSub, 2, 1) = nanmean(acc_trialReport(testOrd == 1, 1));
        acc_subj_orderUp(iSub, 2, 2) = nanmean(acc_trialReport(testOrd == 0, 2));
        % ------third order: only reconstruction reports------
        acc_subj_orderUp(iSub, 3, 1) = nanmean(acc_trialReport(testOrd == 0, 3));
        acc_subj_orderUp(iSub, 3, 2) = nanmean(acc_trialReport(testOrd == 1, 3));

        %% block-wise pattern
        % acc_blc_subj = nan(subLen, (nBlock+1)*2, 3); % 3: (1) content report; (2) position report; (3) both report;
        for iBlc = 1 : ((nBlock + 1) * 2) % 1: post-test block
            iTrl = (iBlc - 1) * uniTrl + 1 : iBlc * uniTrl;
            if iBlc <= (nBlock * 2)
                acc_blc_subj(iSub, iBlc, :) = squeeze(nanmean(acc_trial_subj(iSub, 1 : 3, iTrl), 3));
            else
                acc_blc_subj(iSub, iBlc, 3) = squeeze(nanmean(acc_trial_subj(iSub, 3, iTrl), 3));
            end
        end
        %% percentage of fully correct trials
        con_trials = squeeze(acc_trial_subj(iSub, 1, 1 : nEpi));
        pos_trials = squeeze(acc_trial_subj(iSub, 2, 1 : nEpi));
        trialPerc_subj(iSub, 1) = length(find(con_trials == 1)) / nEpi;
        trialPerc_subj(iSub, 2) = length(find(pos_trials == 1)) / nEpi;
        trialPerc_subj(iSub, 3) = length(find(con_trials == 1 & pos_trials == 1)) / nEpi;

        %% save accuracy and choice for each subject for subsequent LMM and GLMM analysis
        acc_trial_col = (squeeze(acc_trial_subj(iSub, [1, 2, 3], 1 : nEpi)))'; % nEpi * 3 meas
        acc_col    = []; %% accuracy
        trlCor_col = []; %% fully correct: 1-yes, 0-no
        meas_col   = []; %% 3 measures
        trlNo_col  = [];
        for iM = 1 : 3 %% 3 measures: content, position and reconstruction
            acc_iM = acc_trial_col(:, iM);
            trlNo_iM = (1 : 1 : nEpi)';
            trlNo_iM(isnan(acc_iM)) = [];
            acc_iM(isnan(acc_iM)) = [];
            acc_col   = [acc_col; acc_iM];
            trlCor_iM = acc_iM;
            trlCor_iM(trlCor_iM ~= 1) = 0;
            trlCor_col = [trlCor_col; trlCor_iM];
            meas_col   = [meas_col; repmat(iM - 1, length(acc_iM), 1)];
            trlNo_col  = [trlNo_col; trlNo_iM];
        end
        group_Col = repmat(iGrp - 1, length(acc_col), 1);
        subj_Col  = repmat(iSub - 1, length(acc_col), 1);
        acc_lmm_subj = [acc_lmm_subj; group_Col, subj_Col, trlNo_col, meas_col, acc_col, trlCor_col];

        %% save accuracy change between older and younger adults, with the mean in the latter as references
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
            for iM = 1 : 3 %% 3 measures: content, position and reconstruction
                acc_iM = acc_trial_col(:, iM);
                trlNo_iM = (1 : 1 : nEpi)';
                trlNo_iM(isnan(acc_iM)) = [];
                acc_iM(isnan(acc_iM)) = [];
                acc_col   = [acc_col; acc_iM];
                trlCor_iM = acc_iM;
                trlCor_iM(trlCor_iM ~= 1) = 0;
                trlCor_col = [trlCor_col; trlCor_iM];
                meas_col   = [meas_col; repmat(iM - 1, length(acc_iM), 1)];
                trlNo_col  = [trlNo_col; trlNo_iM];
            end
            group_Col = repmat(iGrp - 1, length(acc_col), 1);
            subj_Col  = repmat(iSub - 1, length(acc_col), 1);
            acc_change_lmm_subj = [acc_change_lmm_subj; group_Col, subj_Col, trlNo_col, meas_col, acc_col, trlCor_col];
        end
        %     if iGrp == 2 % older group: save the data
        %         save([bhvDataDir, '/acc_change_lmm_subj_exp1.mat'], 'acc_change_lmm_subj');
        %     end

        %% Binding evidence: check if the correct/incorrect in the different dimension influence the accuracy in the dimension
        % by controling the previous transition in the same dimension as correct or incorrect.
        % added by rxj @ Oct 20 2024
        conTrue_threeRep = conTrue_col; % no reconsOnly trials in Experiment 1
        conRep_threeRep  = conRep_col;
        locTrue_threeRep = locTrue_col;
        locRep_threeRep  = locRep_col;
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
        slotCorr_marg_con = nan(nEpi, nTrans); % 1: the response is correct; 0: incorrect response;
        slotCorr_marg_pos = nan(nEpi, nTrans);
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
                transError_con = nan(nEpi, (nTrans - 1));
                transError_pos = nan(nEpi, (nTrans - 1));
            elseif ij == 2
                slotCorr_con = slotCorr_join_con;
                slotCorr_pos = slotCorr_join_pos;
                transError_con = nan(nEpi, (nTrans - 1));
                transError_pos = nan(nEpi, (nTrans - 1));
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
            transError_con_col = reshape(transError_con, [nEpi * (nTrans - 1), 1]);
            transError_loc_col = reshape(transError_pos, [nEpi * (nTrans - 1), 1]);
            error_con_len = [length(find(transError_con_col == 1)), length(find(transError_con_col == 2)), ...
                             length(find(transError_con_col == 3)), length(find(transError_con_col == 4))];

            error_pos_len = [length(find(transError_loc_col == 1)), length(find(transError_loc_col == 2)), ...
                             length(find(transError_loc_col == 3)), length(find(transError_loc_col == 4))];
            %%% ------ First and second half per retrieval test ------
            transError_con_firstHalf_col = reshape(transError_con(:, 1:2), [nEpi * 2, 1]);
            transError_loc_firstHalf_col = reshape(transError_pos(:, 1:2), [nEpi * 2, 1]);
            error_con_firstThree_len = [length(find(transError_con_firstHalf_col == 1)), length(find(transError_con_firstHalf_col == 2)), ...
                                        length(find(transError_con_firstHalf_col == 3)), length(find(transError_con_firstHalf_col == 4))];
            error_pos_firstThree_len = [length(find(transError_loc_firstHalf_col == 1)), length(find(transError_loc_firstHalf_col == 2)), ...
                                        length(find(transError_loc_firstHalf_col == 3)), length(find(transError_loc_firstHalf_col == 4))];
            transError_con_secondHalf_col = reshape(transError_con(:, 3:4), [nEpi * 2, 1]);
            transError_loc_secondHalf_col = reshape(transError_pos(:, 3:4), [nEpi * 2, 1]);
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
        conLure  = zeros(nEpi, 1);
        locLure  = zeros(nEpi, 1);
        bothLure = zeros(nEpi, 2); % 2: item and location dimension in the full retrieval
        % ---------- note down the trial index where the lure stimulus was picked ----------
        lure_trialIdx = zeros(nEpi, 3); % 3: item, location and full retrieval
        for i = 1 : nEpi
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

            % full retrieval
            bothRep_i      = bothRep_col{i}; % 2 row: content and position
            bothLure(i, 1) = (sum(bothRep_i(1, :) == 6) | bothRep_i(1, end) ~= 0);
            bothLure(i, 2) = (sum(bothRep_i(2, :) == 6) | bothRep_i(2, end) ~= 0);

        end
        % ---------- partial retrieval ----------
        FA_lure_subj(iSub, 1, 1) = sum(conLure) / nEpi;
        FA_lure_subj(iSub, 2, 1) = sum(locLure) / nEpi;

        % ---------- full retrieval ----------
        FA_lure_subj(iSub, 1, 2) = sum(bothLure(:, 1)) / nEpi;
        FA_lure_subj(iSub, 2, 2) = sum(bothLure(:, 2)) / nEpi;

        lure_trialIdx(:, 1) = conLure;
        lure_trialIdx(:, 2) = locLure;
        lure_trialIdx(:, 3) = (bothLure(:, 1) == 1 | bothLure(:, 2) == 1);

    end
    acc_trial_group{iGrp} = acc_trial_subj;
    acc_group{iGrp} = acc_subj;
    rt_group{iGrp} = rt_subj;
    acc_subj_orderUP_group{iGrp} = acc_subj_orderUp;
    acc_blc_group{iGrp}    = acc_blc_subj;
    trialPerc_group{iGrp} = trialPerc_subj;
    %%% lmm
    acc_lmm_group = [acc_lmm_group; acc_lmm_subj];
    %%% glmm
    resp_glmm_group = [resp_glmm_group; resp_glmm_subj];

    % ----proportions and total counts----
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
    allMeas_inOne(:, 5) = acc_subj(:, 9);
    % (6) difference between reconstruction report and prediction from
    % marginal reports
    allMeas_inOne(:, 6) = acc_subj(:, 3) - acc_subj(:, 9);
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
%% save data for Figure 2A and 2B in the FlexibleBinding paper
FBdata_folder = [bhvDataDir, '/FlexibleBindingPaper-Data/']; % flexible binding data folder
save([FBdata_folder, 'Fig2A_acc_blc_group.mat'], 'acc_blc_group');
save([FBdata_folder, 'Fig2B_acc_group.mat'], 'acc_group');
save([FBdata_folder, 'Fig2D_FA_lure_group.mat'], 'FA_lure_group');
save([FBdata_folder, 'Fig2GHI_binds_conPctr_group.mat'], 'binds_conPctr_group');

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
save([FBdata_folder, 'Fig2E_transAcc_plot_mjAvg_group.mat'], 'transAcc_plot_mjAvg_group');


%% color settings
colorSets = [0.98, 0.72, 0.69; ...
             0.97, 0.85, 0.67; ...
             0.33, 0.73, 0.83; ...
             0.72, 0.80, 0.88; ...
             0.78, 0.50, 0.75; ...
             0.86, 0.80, 0.89; ...
             0.54, 0.67, 0.20; ...
             0.82, 0.92, 0.78; ...
             0.75, 0.56, 0; ...
             0.40, 0.40, 0.40];

iGrp = 1;
color_Grp = colorSets([iGrp+1, iGrp+3, iGrp+5], :);

%% Figure 2A: Experiment 1
%%% seqMemTask_v1_anal_summary.m
%%% ---------Figure 2A. Block-wise learning curves---------
acc_group_comp = cell(1, nGroup);
for iGrp = 1 : nGroup
    acc_group_comp{iGrp} = acc_group{iGrp}(:, [1,2,3]);
end
barPos = [1, 1.5, 2; % younger: content, position, reconstruction
          3, 3.5, 4]; % older
x_rand_group = cell(nGroup, 3);
for iGrp = 1 : nGroup 
    for iC = 1 : 3 % content, position and both-pre
        barPos_i = barPos(iGrp, iC);
        acc_iGrp = acc_group_comp{iGrp}(:, iC);
        x_rand_group{iGrp, iC} = unifrnd(barPos_i - 0.1, barPos_i + 0.1, length(acc_iGrp), 1);
    end
end

color_Grp_curve = color_Grp;
blockLines = [2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5]; % using vertical line to label the block
for iGrp = 1 : nGroup %% younger and older adults
    acc_blc_iGrp = acc_blc_group{iGrp}(:, 1 : 2*8, :);
    [acc_avg, acc_sem] = Mean_and_Se(acc_blc_iGrp, 1);
    acc_avg = squeeze(acc_avg); % (nBlock)*2 * 3; % 3: (1) content report; (2) position report; (3) both report;
    acc_sem = squeeze(acc_sem);
    figure('Position', [100 100 300 150]), clf;
    for iP = 1 : 3 % the first 3 rows: content, position and both
        errorbar(1 : 1 : (nBlock*2), acc_avg(1 : 1 : (nBlock*2), iP), acc_sem(1 : 1 : (nBlock*2), iP), 'Color', color_Grp_curve(iP, :), 'LineStyle', '-', 'LineWidth', 1); hold on;
        plot(1 : 1 : (nBlock*2), acc_avg(1 : 1 : (nBlock*2), iP), 'Marker', '.', 'MarkerSize', 15, 'Color', color_Grp_curve(iP, :), 'LineStyle', 'none'); hold on;
    end
    ylim([0, 1]);
    xlim([0, ((nBlock)*2)+1]);
    for iBlc = 1 : length(blockLines)
        plot([blockLines(iBlc), blockLines(iBlc)], ylim, 'Color', [0.6, 0.6, 0.6], 'LineStyle', ':', 'LineWidth', 0.4); hold on;
    end
    %plot(xlim, [1/(nTrans+nDtr), 1/(nTrans+nDtr)], 'k--', 'LineWidth', 0.6); hold on;
    set(gca, 'LineWidth', 0.8);
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', blockLines, 'XTickLabel', '');
    set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', '');
    box off;
    %%% statistical test
    statMat = nan(nBlock*2, 2);
    for iB = 1 : (nBlock*2)
        [h, p, ci, stats] = ttest(acc_blc_iGrp(:, iB, 1), acc_blc_iGrp(:, iB, 2));
        statMat(iB, :) = [p, stats.tstat];
    end
end

%% Figure 2B: Experiment 1
%%% seqMemTask_v1_anal_summary.m
%% ---------Behavioral paper: figure 2B-left panel---------
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
for iGrp = 1 : nGroup 
    %color_Grp = colorSets([iGrp, iGrp+2, iGrp+4], :);
    for iC = 1 : 3 % content, position and both-pre
        barPos_i = barPos(iGrp, iC);
        acc_iGrp = acc_group_comp{iGrp}(:, iC);
        [acc_avg, acc_sem] = Mean_and_Se(acc_iGrp, 1);
        b = bar(barPos(iGrp, iC), acc_avg, 0.45, 'LineStyle', '-', 'LineWidth', barLineWid); hold on;
        b.FaceColor = color_Grp(iC, :);
        x_rand = x_rand_group{iGrp, iC};
        plot(x_rand, acc_iGrp, 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0.9, 0.9, 0.9], 'LineStyle', 'none'); hold on;
        errorbar(barPos(iGrp, iC), acc_avg, acc_sem, 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
    end
end
xlim([0.5, 4.5]);
ylim([0, 1]);
if figKey == 0
    % ------For presentation------
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', '', 'XTickLabel', '');
    set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', 0 : 0.5 : 1); 
elseif figKey == 1
    % ------For Adobe Illustrator------
    set(gca, 'LineWidth', 0.8);
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', '', 'XTickLabel', '');
    set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', {'', '', ''});
end
box off;

%% ---------Behavioral paper: Figure 2B-right panel---------
acc_change_comp = nan(length(subjList_old), 3); % content, position and reconstruction
for iC = 1 : 3 % content, position and both-pre
    acc_change_comp(:, iC) = (acc_group_comp{2}(:, iC) - nanmean(acc_group_comp{1}(:, iC))) ./ nanmean(acc_group_comp{1}(:, iC));
end
figKey = 1;
if figKey == 0
    errLineWid = 3;
    refLineWid = 1;
elseif figKey == 1
    errLineWid = 2;
    refLineWid = 0.5;
end
figure('Position', [100 100 80 150]), clf;
[acc_avg, acc_sem] = Mean_and_Se(acc_change_comp, 1);
for iC = 1 : 3 % content, position and both-pre
    errorbar(iC, acc_avg(iC), acc_sem(iC), 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
    plot(iC, acc_avg(iC), 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', color_Grp(iC, :), 'MarkerEdgeColor', [0, 0, 0], 'LineStyle', '-', 'LineWidth', refLineWid); hold on;
end
xlim([0.5, 3.5]);
ylim([-0.6, 0]);
if figKey == 0
    % ------For presentation------
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', '', 'XTickLabel', '');
    %set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', 0 : 0.5 : 1); 
elseif figKey == 1
    % ------For Adobe Illustrator------
    set(gca, 'LineWidth', 0.8);
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', '', 'XTickLabel', '');
    set(gca, 'YTick', -0.6 : 0.2 : 0, 'YTickLabel', {'', '', ''}); 
end
box off;

%% Figure 2D: Experiment 1
%%% seqMemTask_v1_lureEffect.m
%%% ------ Figure 2D in the behavioral manuscript: plotting the fasle alarm rate for both YA and OA: merge partial and full retrieval ------
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

%% Figure 2E: Experiment 1
%%% SeqMemTask_v1_anal_summary.m
%%% --------- Figure 2E: Average across partial (or margial) and full (or reconstruction) reports ----------
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
            %b = bar(barPos_i(jDm), tAcc_avg(jDm), 0.45, 'LineStyle', '-', 'LineWidth', barLineWid); hold on;
            %b.FaceColor = color_iCp(jDm, :);
            %b.EdgeColor = color_iEd(jDm, :);
            %b.LineWidth = barLineWid;
            %x_rand = unifrnd(barPos_i(jDm) - 0.1, barPos_i(jDm) + 0.1, size(transAcc_plot_ij, 1), 1);
            %plot(x_rand, transAcc_plot_ij(:, jDm), 'Marker', 'o', 'MarkerSize', 4, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0.9, 0.9, 0.9], 'LineStyle', 'none'); hold on;
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

%% Figure 2G: Experiment 1
%%% SeqMemTask_v1_anal_summary.m
%%% ---------- Figure 2G: grand average of the binding score for the other dimension correct vs. incorrect in YA and OA ----------
% added by XR @ Sep 16 2025
% average across the condition that the preceding reports in the correct
% dimension is correct or not; bidirectional binding bdtween item and
% loction; partial (or marginal) and full (or reconstruction) reports
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


%% Figure 2H and 2I: Experiment 1
%%% SeqMemTask_v1_anal_summary.m
%%% ---------- Figure 2H and 2I in behavioral manuscript: simplying the effect in the above figures ----------
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


