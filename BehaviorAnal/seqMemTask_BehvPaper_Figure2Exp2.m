% seqMemTask_BehvPaper_Figure2.m
% revision by XR @ Dec 21 2025
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
folder     = '/Users/ren/Projects-NeuroCode/MyExperiment/Aging-SeqMemTask';
bhvDataDir = [folder, '/AgingReplay-OnlineData'];
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

expList = {'interleaved'};
nGroup  = 2; %% younger and older adults

%% ---------- Figure 2E in Experiment 2 (replication in Experiment 1); transition accurcy ----------

%% The binding relevant data was originally saved in seqMemTask_v2_anal_early_late_window_Binding.m
%%% --------- put younger and older adults together for the Interleaved curriculum ----------
binds_conPctr_group_YAOA = cell(2, nGroup); % 2: marginal and reconstruction reports
for iA = 1 : 2 % YA and OA group
    if iA == 1
        groupName = 'younger';
    elseif iA == 2
        groupName = 'older';
    end
    load([bhvDataDir, '/binds_conPctr_group_', groupName, '_exp2.mat'], 'binds_conPctr_group');
    binds_conPctr_group_YAOA(:, iA) = binds_conPctr_group(:, 3);
end

%% ---------- Figure 2G in Experiment 2 (replication in Experiment 1): grand average of the binding score for the other dimension correct vs. incorrect in YA and OA ----------
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

binds_dataAnal = binds_conPctr_group_YAOA;
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
binds_dataAnal = binds_conPctr_group_YAOA;
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

