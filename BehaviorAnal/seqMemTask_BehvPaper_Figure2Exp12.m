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
% load([FBdata_folder, 'Fig2A_acc_blc_group.mat'], 'acc_blc_group');
% load([FBdata_folder, 'Fig2B_acc_group.mat'], 'acc_group');
FA_lure_group_Exp1 = load([FBdata_folder, 'Fig2D_FA_lure_group.mat'], 'FA_lure_group');
FA_lure_group_Exp1 = FA_lure_group_Exp1.FA_lure_group; % cell(1, 2)

transAcc_plot_mjAvg_group_Exp1 = load([FBdata_folder, 'Fig2E_transAcc_plot_mjAvg_group.mat'], 'transAcc_plot_mjAvg_group');
transAcc_plot_mjAvg_group_Exp1 = transAcc_plot_mjAvg_group_Exp1.transAcc_plot_mjAvg_group;

binds_conPctr_group_Exp1 = load([FBdata_folder, 'Fig2GHI_binds_conPctr_group.mat'], 'binds_conPctr_group');
binds_conPctr_group_Exp1 = binds_conPctr_group_Exp1.binds_conPctr_group;

%% Experiment 2: data for Figure 2D, 2E and 2GHI
% save([FBdata_folder, 'Fig2C_acc_group.mat'], 'acc_group');
% save([FBdata_folder, 'Fig2C_acc_subj_orderUP_group.mat'], 'acc_subj_orderUP_group');
FA_lure_group_Exp2 = load([FBdata_folder, 'Fig2D_Exp2_FA_lure_group.mat'], 'FA_lure_group');
FA_lure_group_Exp2 = FA_lure_group_Exp2.FA_lure_group;

transAcc_plot_mjAvg_group_Exp2 = load([FBdata_folder, 'Fig2E_Exp2_transAcc_plot_mjAvg_group.mat'], 'transAcc_plot_mjAvg_group');
transAcc_plot_mjAvg_group_Exp2 = transAcc_plot_mjAvg_group_Exp2.transAcc_plot_mjAvg_group;

binds_conPctr_group_Exp2 = load([FBdata_folder, 'Fig2GHI_Exp2_binds_conPctr_group.mat'], 'binds_conPctr_group');
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

% ------ Figure 2GHI: binding accuracy related ------
binds_conPctr_group_Exp12       = cell(2, nGroup);
% ------ OA ------
binds_conPctr_group_Exp12{1, 1} = [binds_conPctr_group_Exp1{1, 1}; binds_conPctr_group_Exp2{1, 1}];
binds_conPctr_group_Exp12{2, 1} = [binds_conPctr_group_Exp1{2, 1}; binds_conPctr_group_Exp2{2, 1}];
% ------ YA ------
binds_conPctr_group_Exp12{1, 2} = [binds_conPctr_group_Exp1{1, 2}; binds_conPctr_group_Exp2{1, 2}];
binds_conPctr_group_Exp12{2, 2} = [binds_conPctr_group_Exp1{2, 2}; binds_conPctr_group_Exp2{2, 2}];

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








