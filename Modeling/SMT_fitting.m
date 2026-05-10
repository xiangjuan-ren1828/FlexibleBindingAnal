function [lliEst, paramsEst, lli_all, fn, fnDir] = SMT_fitting(folder, groupID, expID, Midx, subID, fitWord, refit, nFit, nImg, nPos, nTrans, trialLen, stim_encode, disp_con, disp_pos, disp_rec, resp_con, resp_pos, resp_rec, ...
    context_arr, testOrd_arr, reconsOnly, Minit, angList, angSeq_encode, bindingDirec, optimzerIdx, dataFittingSource)
% write by XR @ Mar 4th 2024
% Fitting the different models
% ------------------Descriptions------------------
% groupID: groups – younger vs. older
% expID: experimental conditions  – interleaved, contentBlocked &
% positionBlocked
% Midx: model name
% fitWord: only marginal report trials, only recons trial, or total

if optimzerIdx == 0 || ~exist("optimzerIdx", 'var') % fminsearchbnd
    fnDir = [folder, '/ModelFitting_Results/', groupID, '/', expID, '/', Midx, '-ModelFits-', fitWord, '/'];
elseif optimzerIdx == 1 % bads
    if dataFittingSource == 0     % responses from the marginal reports
        fnDir = [folder, '/ModelFitting_BadsResults/', groupID, '/', expID, '/', Midx, '-ModelFits-', fitWord, '/'];
    end
end
if ~exist(fnDir, 'dir')
    mkdir(fnDir);
end
fn = [fnDir, groupID, '-', expID, '-', subID, '-rep', num2str(nFit), '.mat'];
if exist(fn, 'file') && refit == 0
    load(fn, 'lliEst', 'paramsEst', 'lli_all');
else
    %MList = {'featureIndependent', 'featureIndependentVM', 'featureBindingWeightRW-cOp', 'featureBindingWeightRW', 'featureCompetitionRW'};
    if isequal(Midx, 'featureIndependent')
        theModel = @SMT_featureIndependent;
        LB = [0, 0, 0, 0, 0, 0, 0]; 
        UB = [1, 1, 1, 1, 1, Inf, Inf];

    elseif isequal(Midx, 'featureIndependentVM')
        theModel = @SMT_featureIndependentVM;
        LB = [0, 0, 0, 0, 0, 0, 0, 0]; 
        UB = [1, 1, 1, Inf, 1, 1, Inf, Inf];

    elseif isequal(Midx, 'featureBindingWeightRW') || isequal(Midx, 'featureBindingWeightRW-cOp')
        theModel = @SMT_featureBinding;
        LB = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        UB = [1, 1, 1, 1, 1, Inf, Inf, 1, 1, Inf];

    elseif isequal(Midx, 'featureCompetitionRW') % exactly the same as 'featureCompetition'
        theModel = @SMT_featureCompetition;
        % --------Boundary setting for the original fminsearchbnd--------
        LB = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        UB = [1, 1, 1, 1, 1, Inf, Inf, 1, 1, 1, Inf];
    end
    % ----------options for the fminsearchbnd----------
    options = optimset('MaxFunEvals', length(LB)*1000, 'MaxIter', length(LB)*1000);
    % ----------options for the bads----------
    % options = bads('defaults'); % Get a default OPTIONS struct

    iFit    = 1;
    lli_all    = zeros(1, nFit);
    params_all = zeros(length(LB), nFit);
    while iFit <= nFit
        display([groupID, '-', expID, '-', subID, '-', Midx, 'model-', fitWord, '-Fit', num2str(iFit)]);
        x0 = rand(1, length(LB));
        tic
        if optimzerIdx == 0 || ~exist("optimzerIdx", 'var')
            % ----------fminsearchbnd----------
            [paramsEst, minuslli, exitflag] = fminsearchbnd(@(params)theModel(params, nImg, nPos, nTrans, trialLen, stim_encode, disp_con, disp_pos, disp_rec, resp_con, resp_pos, resp_rec, context_arr, testOrd_arr, reconsOnly, Minit, angList, angSeq_encode, bindingDirec), ...
                x0, LB, UB, options)
        elseif optimzerIdx == 1
            % ----------bads----------
            [paramsEst, minuslli, exitflag, output] = bads(@(params)theModel(params, nImg, nPos, nTrans, trialLen, stim_encode, disp_con, disp_pos, disp_rec, resp_con, resp_pos, resp_rec, context_arr, testOrd_arr, reconsOnly, Minit, angList, angSeq_encode, bindingDirec), ...
                x0, LB, UB, pLB, pUB, [], options)
        end
        toc
        
        if exitflag == 1
            params_all(:, iFit) = paramsEst;
            lli_all(iFit)       = - minuslli;
            iFit = iFit + 1;
        end
    end
    [bestlli, llibest_id] = max(lli_all);
    
    paramsEst = params_all(:, llibest_id);
    lliEst    = bestlli;
    save(fn, 'lliEst', 'paramsEst', 'lli_all');
end