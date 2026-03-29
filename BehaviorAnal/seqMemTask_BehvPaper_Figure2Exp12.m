% seqMemTask_BehvPaper_Figure2Exp12.m
% added by XR @ Jan 5 2026
% ------ Integrating the data in Experiment 1 and Experiment 2 to plot Figures 2D, 2E, 2G, 2H and 2I ------

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
FBdata_folder = [bhvDataDir, '/FlexibleBindingPaper-Data/']; % flexible binding data folder

%% Read the corresponding data from Exp1 and Exp2, saved in 
% seqMemTask_BehvPaper_Figure2Exp1.m
% seqMemTask_BehvPaper_Figure2Exp2.m

N_Exp12 = [39, 27; ...
           24, 23];

%% Experiment 1: data for Figure 2D, 2E and 2GHI
% load([FBdata_folder, 'Fig2A_acc_blc_group_Exp1.mat'], 'acc_blc_group');
% load([FBdata_folder, 'Fig2B_acc_group_Exp1.mat'], 'acc_group');
FA_lure_group_Exp1 = load([FBdata_folder, 'Fig2D_FA_lure_group_Exp1.mat'], 'FA_lure_group');
FA_lure_group_Exp1 = FA_lure_group_Exp1.FA_lure_group; % cell(1, 2)

transAcc_plot_group_Exp1 = load([FBdata_folder, 'Fig2E_transAcc_plot_group_Exp1.mat'], 'transAcc_plot_group'); % transAcc_plot_group = cell(2, nGroup); % 2: marginal and joint
transAcc_plot_group_Exp1 = transAcc_plot_group_Exp1.transAcc_plot_group;

transAcc_plot_mjAvg_group_Exp1 = load([FBdata_folder, 'Fig2E_transAcc_plot_mjAvg_group_Exp1.mat'], 'transAcc_plot_mjAvg_group');
transAcc_plot_mjAvg_group_Exp1 = transAcc_plot_mjAvg_group_Exp1.transAcc_plot_mjAvg_group;

binds_conPctr_group_Exp1 = load([FBdata_folder, 'Fig2GHI_binds_conPctr_group_Exp1.mat'], 'binds_conPctr_group');
binds_conPctr_group_Exp1 = binds_conPctr_group_Exp1.binds_conPctr_group;

%% Experiment 2: data for Figure 2D, 2E and 2GHI
% load([FBdata_folder, 'Fig2C_acc_group_Exp2.mat'], 'acc_group');
% load([FBdata_folder, 'Fig2C_acc_subj_orderUP_group_Exp2.mat'], 'acc_subj_orderUP_group');
FA_lure_group_Exp2 = load([FBdata_folder, 'Fig2D_FA_lure_group_Exp2.mat'], 'FA_lure_group');
FA_lure_group_Exp2 = FA_lure_group_Exp2.FA_lure_group;

transAcc_plot_group_Exp2 = load([FBdata_folder, 'Fig2E_transAcc_plot_group_Exp2.mat'], 'transAcc_plot_group'); % transAcc_plot_group = cell(2, nGroup); % 2: marginal and joint
transAcc_plot_group_Exp2 = transAcc_plot_group_Exp2.transAcc_plot_group;

transAcc_plot_mjAvg_group_Exp2 = load([FBdata_folder, 'Fig2E_transAcc_plot_mjAvg_group_Exp2.mat'], 'transAcc_plot_mjAvg_group');
transAcc_plot_mjAvg_group_Exp2 = transAcc_plot_mjAvg_group_Exp2.transAcc_plot_mjAvg_group;

binds_conPctr_group_Exp2 = load([FBdata_folder, 'Fig2GHI_binds_conPctr_group_Exp2.mat'], 'binds_conPctr_group');
binds_conPctr_group_Exp2 = binds_conPctr_group_Exp2.binds_conPctr_group;


%% Integrating data from Exp1 and Exp2
nGroup = 2;
% ------ Figure 2D: False alarm rate ------
FA_lure_group_Exp12       = cell(1, nGroup);
FA_lure_group_Exp12{1, 1} = [FA_lure_group_Exp1{1, 1}; FA_lure_group_Exp2{1, 1}];
FA_lure_group_Exp12{1, 2} = [FA_lure_group_Exp1{1, 2}; FA_lure_group_Exp2{1, 2}];

% ------ Figure 2E: Transition accuracy ------
transAcc_plot_mjAvg_group_Exp12       = cell(1, nGroup);
transAcc_plot_mjAvg_group_Exp12{1, 1} = [transAcc_plot_mjAvg_group_Exp1{1, 1}; transAcc_plot_mjAvg_group_Exp2{1, 1}];
transAcc_plot_mjAvg_group_Exp12{1, 2} = [transAcc_plot_mjAvg_group_Exp1{1, 2}; transAcc_plot_mjAvg_group_Exp2{1, 2}];

% ------ Figure 2E: transition accuracy (statistics) ------
transAcc_plot_group_Exp12 = cell(2, nGroup); % 2: marginal and joint
% ------ YA ------
transAcc_plot_group_Exp12{1, 1} = [transAcc_plot_group_Exp1{1, 1}; transAcc_plot_group_Exp2{1, 1}];
transAcc_plot_group_Exp12{2, 1} = [transAcc_plot_group_Exp1{2, 1}; transAcc_plot_group_Exp2{2, 1}];

% ------ OA ------
transAcc_plot_group_Exp12{1, 2} = [transAcc_plot_group_Exp1{1, 2}; transAcc_plot_group_Exp2{1, 2}];
transAcc_plot_group_Exp12{2, 2} = [transAcc_plot_group_Exp1{2, 2}; transAcc_plot_group_Exp2{2, 2}];

% ------ Figure 2GHI: binding accuracy related ------
binds_conPctr_group_Exp12       = cell(2, nGroup);
% ------ YA ------
binds_conPctr_group_Exp12{1, 1} = [binds_conPctr_group_Exp1{1, 1}; binds_conPctr_group_Exp2{1, 1}];
binds_conPctr_group_Exp12{2, 1} = [binds_conPctr_group_Exp1{2, 1}; binds_conPctr_group_Exp2{2, 1}];
% ------ OA ------
binds_conPctr_group_Exp12{1, 2} = [binds_conPctr_group_Exp1{1, 2}; binds_conPctr_group_Exp2{1, 2}];
binds_conPctr_group_Exp12{2, 2} = [binds_conPctr_group_Exp1{2, 2}; binds_conPctr_group_Exp2{2, 2}];


%% Figure 2D statistics
% merge the data across Exp1 and Exp2
% based on the script in seqMemTask_v1_lureEffect.m or seqMemTask_v2_lureEffect.m
% added by XR @ Dec 4th 2025
% variables: item vs. location, partial vs. full retrieval, YA vs. OA
% FA_lure_subj   = nan(subLen, 2, 2); % first 2: item and location; second 2: partial and full retrieval
% ---------- partial retrieval ----------
% FA_lure_subj(iSub, 1, 1) = sum(conLure) / nEpi;
% FA_lure_subj(iSub, 2, 1) = sum(locLure) / nEpi;
% 
% % ---------- full retrieval ----------
% FA_lure_subj(iSub, 1, 2) = sum(bothLure(:, 1)) / nEpi;
% FA_lure_subj(iSub, 2, 2) = sum(bothLure(:, 2)) / nEpi;

FA_lure_overall_YAOA = []; 
for iGrp = 1 : nGroup % YA and OA
    if iGrp == 1
        disp('------YA------')
    elseif iGrp == 2
        disp('------OA------')
    end

    FA_lure_iGrp = FA_lure_group_Exp12{iGrp};

    for ij = 1 : 2 % marginal and reconstruction
        FA_lure_ij = FA_lure_iGrp(:, :, ij);

        %% ------ data for the LMM test ------
        for iCp = 1 : 2 % item and location
            %%% ----------FA rate----------
            FAlure_ij_iCp = FA_lure_ij(:, iCp);

            subLen_tmp = size(FAlure_ij_iCp, 1); % subLen

            %%% ----------Label age groups: 0-YA, 1-OA----------
            groupType_col = repmat((iGrp - 1), [subLen_tmp, 1]);

            %%% ----------Label report types: 0-marginal report, 1-reconstruction report----------
            reportType_col = repmat((ij - 1), [subLen_tmp, 1]);

            %%% ----------Item or Location sequence or dimension----------
            dimType_col = repmat((iCp - 1), [subLen_tmp, 1]);

            %%% ----------Label participant index----------
            subj_col = (1 : 1 : subLen_tmp)' - 1;

            %%% ----------Concatenate the transition calculation across loops----------
            FA_lure_temp = [];
            FA_lure_temp = [groupType_col, subj_col, reportType_col, dimType_col, FAlure_ij_iCp];
            FA_lure_overall_YAOA = [FA_lure_overall_YAOA; FA_lure_temp];
        end
    end
end
%%% save the FA_lure_group.mat for subsequent LMM analysis
%save([FBdata_folder, 'Fig2D_FA_lure_overall_YAOA_stats_Exp12.mat'], 'FA_lure_overall_YAOA');


%% Figure 2E statistics
%%% ------Reorganize the transAcc_plot_group.mat in order to run the statistical analysis------
% ----transition evidence----
% transAcc_plot_group       = cell(2, nGroup); % 2: marginal and joint
% transAcc_plot = nan(size(transAcc_count_ij, 1), 4); % 4: 1-2, accuracy for item transition; 3-4, accuracy for location transition
transAcc_overall_group = []; % ****** For Figure 2E and its SI ******
for iGrp = 1 : nGroup % YA and OA
    if iGrp == 1
        disp('------YA------')
    elseif iGrp == 2
        disp('------OA------')
    end
    for ij = 1 : 2 % marginal and reconstruction
        transAcc_ij = transAcc_plot_group_Exp12{ij, iGrp};

        %% ------ data for the LMM test ------
        for iCp = 1 : 2 % item and location
            iCp_idx = (iCp - 1) * 2 + 1 : iCp * 2;
            transAcc_ij_iCp = transAcc_ij(:, iCp_idx);
            subLen_tmp  = size(transAcc_ij_iCp, 1); % subLen
            condsNs_tmp = size(transAcc_ij_iCp, 2); % 2 columns: correct vs. incorrect

            %%% ----------Label age groups: 0-YA, 1-OA----------
            groupType_col = repmat((iGrp - 1), [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Label report types: 0-marginal report, 1-reconstruction report----------
            reportType_col = repmat((ij - 1), [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Item or Location sequence or dimension----------
            transReportType_col = repmat((iCp - 1), [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Label participant index----------
            subj_col = repmat((1 : 1 : subLen_tmp)' - 1, 1, condsNs_tmp);
            subj_col = reshape(subj_col, [subLen_tmp * condsNs_tmp, 1]);
            
            %%% ----------Label each column and then concatenate the 2 columns into one column ----------
            % 1: (item | item = R) or (loc | loc = R) [R: correct response; W: incorrect response]
            % 0: (item | item = W) or (loc | loc = W) 
            transConds = nan(subLen_tmp, condsNs_tmp);
            transConds(:, 1) = ones(subLen_tmp, 1);
            transConds(:, 2) = zeros(subLen_tmp, 1);
            transConds_col = reshape(transConds, [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Transition accuracy----------
            transAcc_ij_iCp_col = reshape(transAcc_ij_iCp, [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Concatenate the transition calculation across loops----------
            transAcc_cal_temp = [];
            transAcc_cal_temp = [groupType_col, subj_col, reportType_col, transReportType_col, transConds_col, transAcc_ij_iCp_col];
            transAcc_overall_group = [transAcc_overall_group; transAcc_cal_temp];

        end
    end
end
%%% ------save the above data for subsequent statistical tests------
%save([FBdata_folder, 'Fig2E_transAcc_overall_group_stats_Exp12.mat'], 'transAcc_overall_group');

%% Figure2G and Figure 2HI statistics
% Original script in seqMemTask_v1_anal_summary.m (Lines 2301-2390)
% added by rxj @ Dec 30 2024
binds_overall_group      = []; % ****** For Figure 2G ****** 
binds_overall_diff_group = []; % ****** For Figure 2H and 2I****** 
% ------- add a column of age as continuous variable ------
ageList_YA_Exp1 = [22, 21, 25, ...
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
                   23, 21, 25]; % 39 YA
ageList_YA_Exp2 = [24, 22, 23, ...
                   25, 24, 21, ...
                   22, 23, 28, ...
                   21, 22, 23, ...
                   24, 30, 28, ...
                   27, 22, 24, ...
                   21, 20, 26, ...
                   26, 33, 28]; % 24 YA
ageList_OA_Exp1 = [81, 78, 76, ...
                   76, 75, 77, ...
                   79, 76, 78, ...
                   85, 76, ...
                   77, 78, 75, ...
                   79, 75, 75, ... % the last participant this row: age unknown
                   78, 82, 82, ...
                   76, ...
                   75, 79, 78, ...
                   76, 75, 75]; % 27 OA
ageList_OA_Exp2 = [67, 68, 67, ...
                   69, 65, 66, ... % 5th participant: 64???
                   73, 73, 70, ...
                   65, 72, 72, ...
                   65, 69, 67, ...
                   67, 65, 68, ...
                   70, 72, 65, ...
                   67, 67];
ageList_YA_Exp12 = [ageList_YA_Exp1, ageList_YA_Exp2];
ageList_OA_Exp12 = [ageList_OA_Exp1, ageList_OA_Exp2];
for iGrp = 1 : nGroup % YA and OA
    if iGrp == 1
        disp('------YA------')
        ageList_iGrp = ageList_YA_Exp12;
    elseif iGrp == 2
        disp('------OA------')
        ageList_iGrp = ageList_OA_Exp12;
    end
    for ij = 1 : 2 % marginal and reconstruction
        binds_conPctr_group_ij = binds_conPctr_group_Exp12{ij, iGrp}(:, :, 1 : 2);
        for iCp = 1 : 2 % (Con|Pos) and (Pos|Con)
            %% ------ Data for Figure 2G ------
            binds_Cp_subj = binds_conPctr_group_ij(:, :, iCp); % subLen * 4
            subLen_tmp  = size(binds_Cp_subj, 1); % subLen
            condsNs_tmp = size(binds_Cp_subj, 2); % 4 columns (or sub-conditions)

            %%% ----------Label age groups: 0-YA, 1-OA----------
            groupType_col = repmat((iGrp - 1), [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Label report types: 0-marginal report, 1-reconstruction report----------
            reportType_col = repmat((ij - 1), [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Label binding directions: 0-(item | pos), 1-(pos | item)----------
            bindsDirection_col = repmat((iCp - 1), [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Label participant index----------
            subj_col = repmat((1 : 1 : subLen_tmp)' - 1, 1, condsNs_tmp);
            subj_col = reshape(subj_col, [subLen_tmp * condsNs_tmp, 1]);
            
            %%% ----------Label each column and then concatenate the 4 columns into one column ----------
            % ------(item | pos)------
            % 1: (item | item = R, pos = R), (item | item = W, pos = R) [R: correct response; W: incorrect response]
            % 0: (item | item = R, pos = W), (item | item = W, pos = W)
            % ------(pos | item)------
            % 1: (pos | pos = R, item = R), (pos | pos = W, item = R)
            % 0: (pos | pos = R, item = W), (pos | pos = W, item = W)
            bindsCond = nan(subLen_tmp, condsNs_tmp);
            bindsCond(:, [1, 3]) = ones(subLen_tmp, 2);
            bindsCond(:, [2, 4]) = zeros(subLen_tmp, 2);
            bindsCond_col = reshape(bindsCond, [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Binds calculations----------
            binds_Cp_subj_col = reshape(binds_Cp_subj, [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Age----------
            age_subj_col = reshape(repmat(ageList_iGrp, condsNs_tmp, 1), [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Concatenate the binds calculation across loops----------
            binds_cal_temp = [];
            binds_cal_temp = [groupType_col, subj_col, reportType_col, bindsDirection_col, bindsCond_col, binds_Cp_subj_col, age_subj_col];
            binds_overall_group = [binds_overall_group; binds_cal_temp];

            %% ------ Data for original Figure 2H and 2I ------
            %%% ----------Calculate the differences----------
            % !!!!!! make the cost and benefit the same direction !!!!!!
            binds_Cp_subj_diff = [binds_Cp_subj(:, 1) - binds_Cp_subj(:, 2), binds_Cp_subj(:, 3) - binds_Cp_subj(:, 4)];
            subLen_tmp  = size(binds_Cp_subj_diff, 1); % subLen
            condsNs_tmp = size(binds_Cp_subj_diff, 2); % 2 columns: differences

            %%% ----------Label age groups: 0-YA, 1-OA----------
            groupType_col = repmat((iGrp - 1), [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Label report types: 0-marginal report, 1-reconstruction report----------
            reportType_col = repmat((ij - 1), [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Label binding directions: 0-(item | pos), 1-(pos | item)----------
            bindsDirection_col = repmat((iCp - 1), [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Label participant index----------
            subj_col = repmat((1 : 1 : subLen_tmp)' - 1, 1, condsNs_tmp);
            subj_col = reshape(subj_col, [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Label the cost and benefit condition----------
            bindsCond = [zeros(subLen_tmp, 1), ones(subLen_tmp, 1)];
            bindsCond_col = reshape(bindsCond, [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Binds difference calculations----------
            binds_Cp_subj_col = reshape(binds_Cp_subj_diff, [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Age----------
            age_subj_col = reshape(repmat(ageList_iGrp, condsNs_tmp, 1), [subLen_tmp * condsNs_tmp, 1]);

            %%% ----------Concatenate the binds calculation across loops----------
            binds_cal_diff_temp = [];
            binds_cal_diff_temp = [groupType_col, subj_col, reportType_col, bindsDirection_col, bindsCond_col, binds_Cp_subj_col, age_subj_col];
            binds_overall_diff_group = [binds_overall_diff_group; binds_cal_diff_temp];

        end
    end
end
%%% ------save the above data for subsequent statistical tests------
% save([FBdata_folder, 'Fig2G_binds_overall_group_stats_Exp12.mat'], 'binds_overall_group');
% save([FBdata_folder, 'Fig2HI_binds_overall_diff_group_stats_Exp12.mat'], 'binds_overall_diff_group');

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

%% ------ Figure 2D in the behavioral manuscript: plotting the fasle alarm rate for both YA and OA: merge partial and full retrieval ------
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
    FA_lure_iGrp = FA_lure_group_Exp12{iGrp}; % subLen * 2(item/loc) * 2(partial/full)
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

%% --------- Figure 2E: Average across partial (or margial) and full (or reconstruction) reports ----------
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

transAcc_drop_mjAvg_group_Exp12 = cell(1, nGroup);
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
        transAcc_plot_ij     = transAcc_plot_mjAvg_group_Exp12{1, iGrp}(:, iDm_idx); % subj * 4
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

    %%% ------ calculate the drop percentage for both YA and OA ------
    transAcc_plot_iGrp_temp_avg = transAcc_plot_mjAvg_group_Exp12{1, iGrp};
    transAcc_drop_mjAvg_group_Exp12{1, iGrp} = [
        (transAcc_plot_iGrp_temp_avg(:, 1) -  transAcc_plot_iGrp_temp_avg(:, 2)) ./ transAcc_plot_iGrp_temp_avg(:, 1), ...
        (transAcc_plot_iGrp_temp_avg(:, 3) -  transAcc_plot_iGrp_temp_avg(:, 4)) ./ transAcc_plot_iGrp_temp_avg(:, 3)
    ];
end
%% ------ Check the accuracy drop percentage ------
for iGrp = 1 : nGroup % YA and OA
    if iGrp == 1
        disp('========== YA ==========');
    elseif iGrp == 2
        disp('========== OA ==========');
    end
    transAcc_drop_mjAvg_iGrp = transAcc_drop_mjAvg_group_Exp12{1, iGrp};
    [transAcc_avg, transAcc_sem] = Mean_and_Se(transAcc_drop_mjAvg_iGrp)
end


%% ---------- Figure 2G: grand average of the binding score for the other dimension correct vs. incorrect in YA and OA ----------
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
binds_dataAnal = binds_conPctr_group_Exp12;
binds_dataAnal_grandAvg = cell(1, nGroup);
bindScore_grandAvg = cell(1, nGroup);
for iGrp = 1 : nGroup
    if iGrp == 1
        disp('------YA------')
        ageList_iGrp = ageList_YA_Exp12;
    elseif iGrp == 2
        disp('------OA------')
        ageList_iGrp = ageList_OA_Exp12;
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
        errorbar(barPos_i(jDm), tAcc_avg(jDm), tAcc_sem(jDm), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 1); hold on; % errLineWid
        plot(barPos_i(jDm), tAcc_avg(jDm), 'Marker', 'o', 'MarkerSize', 4, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorFace_jDm, 'LineStyle', '-'); hold on;
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
%%% ------ correlation between age and binding score ------
ageCol = [ageList_YA_Exp12, ageList_OA_Exp12]';
binds_score_YAOA = [bindScore_grandAvg{1, 1}; bindScore_grandAvg{1, 2}];
val_idx = ~isnan(ageCol) & ~isnan(binds_score_YAOA);
[r, p] = corr(ageCol, binds_score_YAOA);

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

%% Correlation between age and binding score


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
binds_dataAnal = binds_conPctr_group_Exp12;
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
            plot(barPos_i(jDm), tAcc_avg(jDm), 'Marker', 'o', 'MarkerSize', 4, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorFace_jDm, 'LineStyle', '-'); hold on;
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








