function [minuslli, Mc, Mp, Mfc_c, Mfc_p, Mb] = SMT_featureBinding(params, nImg, nPos, nTrans, trialLen, stim_encode, disp_con, disp_pos, disp_rec, resp_con, resp_pos, resp_rec, context_arr, testOrd_arr, reconsOnly, Minit, angList, angSeq_encode, bindingDirec)
% modified by XR @ Nov 7 2025
% This is modified based on SeqMem_featureBindingRW.m
% In the $SeqMem_featureBindingRW.m$, if we consider the binding
% association learning for one dimension, I simply remove the transition
% learning component, which means that participant makes decisions based on
% the binding component exclusively
% In this version, for the binding component, adding the weighted parameter 
% to consider the contributions from the transition learning and binding association
% learning
% =========================================================================
% revision by XR @ 03/21/2024
% Content-position binding hypothesis
% =================================================
% write by XR @ 01/20/2024
% =================================================
% nImg: 8 unique images
% nPos: 8 positions uniformly distributed on a circle
% nTrans: each sequence contains 5 transitions; 5 transitions = 5 categories
% stim: [trialLen * trans, 2], 2 denotes content and position transitions;
%       the transitions are denoted by the index of the corresponding
%       stimulus
% disp_con, disp_pos, disp_rec: 64 * 1 cells, the displayed stimuli during
%                               reporting stage
% resp_con, resp_pos, resp_rec: 64 * 1 cells
% context_arr: [trialLen, 1], the context of each trial; could be 1,2,3,4
% testOrd_arr: 0-content+position+recons; 1-position+content+recons;
% reconsOnly_arr: 0-3 reports trials; 1-reconstruction only trials;
% angList: the 8 unique angles used for a specific participant (ranges from
% -pi to pi)
% angSeq_encode: 5 angles in each trial

if isequal(bindingDirec, 'pOc') % 'position on content': position retrieval depends on content encoding
    % --------Independent transition encoding--------
    gamma   = params(1); % discount factor
    alpha_c = params(2); % learning rate for content transition
    alpha_p = params(3); % learning rate for position transition
    % --------Item-to-context binding--------
    eta0_c  = params(4); % scaling factor of the item-to-context association, content
    eta0_p  = params(5); % position
    eta     = params(6); % exponential parameter of the item-to-context association
    % --------Content-to-positon binding strength--------
    kappa   = params(7); % precision in the VM function
    phi     = params(8); % decay parameter for the binding strength
    % --------Weight for encoding strategy--------
    omega_p = params(9); % used in the location retrieval
    % --------Decision noise (inverse temperature)--------
    beta    = params(10); % inverse temperature

elseif isequal(bindingDirec, 'cOp') % 'content on position': content retrieval depends on position encoding
    % --------Independent transition encoding--------
    gamma   = params(1); % discount factor
    alpha_c = params(2); % learning rate for content transition
    alpha_p = params(3); % learning rate for position transition
    % --------Item-to-context binding--------
    eta0_c  = params(4); % scaling factor of the item-to-context association, content
    eta0_p  = params(5); % position
    eta     = params(6); % exponential parameter of the item-to-context association
    % --------Content-to-positon binding strength--------
    kappa   = params(7); % precision in the VM function
    phi     = params(8); % decay parameter for the binding strength
    % --------Weight for encoding strategy--------
    omega_c = params(9); % used in the item retrieval
    % --------Decision noise (inverse temperature)--------
    beta    = params(10); % inverse temperature
end

%%
states = nImg;
Mc  = Minit; % eye(states, states); % 8 * 8; eye(states, states);
Mp  = Minit; % eye(states, states);
I   = eye(states, states);
Mb  = zeros(nImg, nPos); % the initial matrix of the content-position binding
Mfc_c = zeros(nImg, 4); % 8 * 4; feature-to-context association matrix
Mfc_p = zeros(nPos, 4);

minuslli = 0;
% lliTrial_all = [];
for iT = 1 : trialLen % only those non-reconstruction-only trials, 48 trials
    % Encoding stage: transition sequence from the in the current trial 
    stim_idx = (iT - 1) * nTrans + 1 : iT * nTrans;
    content_i  = stim_encode{iT, 1}; % label with numbers
    position_i = stim_encode{iT, 2}; % label with numbers
    trans_ang  = angSeq_encode{iT};
    % transition sequence history
    context_history = context_arr(1 : iT); % history of context until current trial
    con_history_temp = (cell2mat(stim_encode(1 : iT, 1)))';
    pos_history_temp = (cell2mat(stim_encode(1 : iT, 2)))';
    % history of content stimulus
    con_history = reshape(con_history_temp, [size(con_history_temp, 1)*size(con_history_temp, 2), 1]); 
    % history of position stimulus
    pos_history = reshape(pos_history_temp, [size(pos_history_temp, 1)*size(pos_history_temp, 2), 1]); 
    % Reporting stage: stimuli displayed during the reporting stage: 5 targets + 1 distractors
    dispImg = (disp_con{iT})'; % one column, the index of the displayed stimulus: 5 targets + 1 distractor
    dispPos = (disp_pos{iT})'; 
    %dispRec = (disp_rec{iT})'; 
    % participants' report
    conRep_iT = resp_con{iT};
    posRep_iT = resp_pos{iT};
    %recRep_iT = resp_rec{iT};

    % report type
    reportType = testOrd_arr(iT); % 0-content+position+reconstruction; 1-position+content+reconstruction; 2-only reconstruction.
    % combination label: which context (or current context)
    comb_iT = context_history(iT);
    context_vector = zeros(1, 4); % 4 elements: if the context occur, the element will be 1.
    context_vector(comb_iT) = 1;
    %% --------Learning stage: update all associations--------
    content_occur_iT = zeros(nImg, 1);
    for ii = 1 : nTrans
        s_i_con = content_i(ii);
        s_i_pos = position_i(ii);
        content_occur_iT(s_i_con) = 1; % if the content occurs in the current trial
        if ii < nTrans % only 4 transitions: e.g., A-B, B-C, C-D, D-E
            s_Tgt_con = content_i(ii + 1);
            s_Tgt_pos = position_i(ii + 1);
            % --------Transition relations update--------
            Mc(s_i_con, :) = (1 - alpha_c) .* Mc(s_i_con, :) + alpha_c .* (I(s_Tgt_con, :) + gamma .* Mc(s_Tgt_con, :));
            Mp(s_i_pos, :) = (1 - alpha_p) .* Mp(s_i_pos, :) + alpha_p .* (I(s_Tgt_pos, :) + gamma .* Mp(s_Tgt_pos, :));
        end
        % --------Content-position binding update--------
        mu = trans_ang(ii); % assuming all 8 positions are known
        samePart = 1 / (2 * pi * besseli(0, kappa));
        delta_Mb = samePart * exp(kappa .* cos(angList - mu)); % angList: the eight angle used for each participant
        Mb(s_i_con, :) = phi * Mb(s_i_con, :) + (1 - phi) * delta_Mb'; % the binding strength from the target content to all possible positions

        % --------Feature-to-context update--------
        % the maximal item-to-context strength for content and position is
        % different
        % 1-8: content; 9-16: position
        for iCt = 1 : 4 % 4 contexts
            iCt_find = find(context_history == iCt);
            if ~isempty(iCt_find)
                iCt_find = iCt_find(end);
                iCt_idx = iT - iCt_find; % if iCt_idx = 0, the context occurs in the current trial; if iCt_idx > 0, the context occurs in previous trials
                Mfc_c(s_i_con, iCt) = eta0_c * exp(-eta * (iCt_idx * nTrans + ii - 1));
                Mfc_p(s_i_pos, iCt) = eta0_p * exp(-eta * (iCt_idx * nTrans + ii - 1));
            end
        end
    end
    % --------For the remaining items don't occur in the current trial,
    % their binding strength will decay linearly--------
    Mb(content_occur_iT ~= 1, :) = phi * Mb(content_occur_iT ~= 1, :);

    % ========Update the item-to-context strength for the remaining items (that
    % don't occur in the current trial) to the current context========
    % distance between the items and contexts (similar decay for positive
    % and negative distance)
    itemArr = 1 : 1 : nImg;
    % ------content-------
    resArr_con = setdiff(itemArr, content_i); % find out the items not occur in the current sequence
    for ir = 1 : length(resArr_con)
        res_ir_find = find(con_history == resArr_con(ir));
        if ~isempty(res_ir_find)
            res_ir_find = res_ir_find(end);
            distance = length(con_history) - res_ir_find;
            Mfc_c(resArr_con(ir), context_vector == 1) = eta0_c * exp(-eta * (distance));
        end
    end
    % ------position-------
    resArr_pos = setdiff(itemArr, position_i);
    for ir = 1 : length(resArr_pos)
        res_ir_find = find(pos_history == resArr_pos(ir));
        if ~isempty(res_ir_find)
            res_ir_find = res_ir_find(end);
            distance = length(pos_history) - res_ir_find;
            Mfc_p(resArr_pos(ir), context_vector == 1) = eta0_p * exp(-eta * (distance));
        end
    end

    %% --------Retrieval and reporting stage--------
    % Three kinds reporting trials: (1) content-position-reconstruction;
    % (2) position-content-reconstruction; (3) reconstruction
    % The displayed stimuli on the screen will determine the context
    % Five sequential retrieval
    % --------content report--------
    if reconsOnly(iT) == 0
        if isequal(bindingDirec, 'pOc') % 'position on content': position retrieval depends on content encoding
            if ~isempty(conRep_iT)
                dispImg_choose = zeros(nTrans + 1, 1); % one distractor
                for ic = 1 : nTrans
                    if ic == 1 % Initial item retrieval is based on the item-to-context strength
                        context_signal = sum(Mfc_c(dispImg, :), 1);
                        [argvalue, argmax] = max(context_signal);
                        % --------retrieve the starting state through sigmoid function--------
                        Vc_signal = Mfc_c(dispImg, argmax);
                    else
                        % after the 1st item retrieval, the following items are
                        % retrieved by the transition matrix
                        % ------the alternative values are selected from the remaining unselected items------
                        dispImg_ic = dispImg(dispImg_choose ~= 1);
                        Vc_signal = (Mc(dispImg(tgt_choice), dispImg_ic))';
                    end
                    p_denominator = exp(beta * Vc_signal);
                    p_denominator(isinf(p_denominator)) = exp(700);
                    pCorr = p_denominator ./ sum(p_denominator);
                    pCorr(pCorr < 1e-16)   = 1e-16;

                    %% calculate the likelihood
                    tgt_choice = find(conRep_iT == ic); % indicating participant made a decision
                    if ~isempty(tgt_choice)
                        dtrChoiceID = zeros(length(dispImg), 1);
                        dtrChoiceID(tgt_choice) = 1; % from data
                        dtrChoiceID(dispImg_choose == 1) = []; % delete already selected items

                        minuslli = minuslli + sum(dtrChoiceID .* log(pCorr));
                        dispImg_choose(tgt_choice) = 1; % label the already selected items
                        %lliTrial_all = [lliTrial_all; sum(dtrChoiceID .* log(pCorr))];
                    end
                end
            end

            %% --------position report--------
            if ~isempty(posRep_iT)
                dispPos_choose = zeros(nTrans + 1, 1); % one distractor
                for ip = 1 : nTrans
                    % ----------When reporting the posiiton sequence, the brain
                    % must maintain the information of the content sequence
                    % (which is a irrelevant information)----------
                    if ip == 1 % Initial item retrieval is based on the item-to-context strength
                        dispPos_ic = dispPos;
                        context_signal = sum(Mfc_p(dispPos, :), 1);
                        [argvalue, argmax] = max(context_signal);
                        % --------retrieve the starting state through sigmoid function--------
                        Vp_signal_ind = Mfc_p(dispPos, argmax);
                    else
                        % after the 1st item retrieval, the following items are
                        % retrieved by the transition matrix
                        % ------the alternative values are selected from the remaining unselected items------
                        dispPos_ic = dispPos(dispPos_choose ~= 1);
                        if isempty(tgt_choice) % ------This will only occur in the reconstruction reports------
                            Vp_signal_ind = zeros(length(dispPos_ic), 1);
                        else
                            Vp_signal_ind = (Mp(dispPos(tgt_choice), dispPos_ic))';
                        end
                    end
                    % ------choice probability from the independent learning------
                    p_denominator_ind = exp(beta * Vp_signal_ind);
                    p_denominator_ind(isinf(p_denominator_ind)) = exp(700);
                    pCorr_ind = p_denominator_ind ./ sum(p_denominator_ind);

                    % ------Using the true content transitions, but not participants' report------
                    Vp_signal_bid = (Mb(content_i(ip), dispPos_ic))'; % !!!!!!!!I use the true content to index the binding strength across 8 posiitons rather than the reported content from participants due to the reported value includes the response noise!!!!!!!!
                    p_denominator_bid = exp(beta * Vp_signal_bid);
                    p_denominator_bid(isinf(p_denominator_bid)) = exp(700);
                    pCorr_bid = p_denominator_bid ./ sum(p_denominator_bid);

                    % ------weighted average------
                    pCorr = omega_p * pCorr_ind + (1 - omega_p) * pCorr_bid;
                    pCorr(pCorr < 1e-16)   = 1e-16;

                    %% calculate the likelihood
                    tgt_choice = find(posRep_iT == ip);
                    if ~isempty(tgt_choice)
                        dtrChoiceID = zeros(length(dispPos), 1);
                        dtrChoiceID(tgt_choice) = 1; % from data
                        dtrChoiceID(dispPos_choose == 1) = []; % delete already selected items

                        minuslli = minuslli + sum(dtrChoiceID .* log(pCorr));
                        dispPos_choose(tgt_choice) = 1; % label the already selected items
                        %lliTrial_all = [lliTrial_all; sum(dtrChoiceID .* log(pCorr))];
                    end
                end
            end

        elseif isequal(bindingDirec, 'cOp') % 'content on position': content retrieval depends on position encoding
            if ~isempty(conRep_iT)
                dispImg_choose = zeros(nTrans + 1, 1); % one distractor
                for ic = 1 : nTrans
                    if ic == 1 % Initial item retrieval is based on the item-to-context strength
                        dispImg_ic = dispImg;
                        context_signal = sum(Mfc_c(dispImg, :), 1);
                        [argvalue, argmax] = max(context_signal);
                        % --------retrieve the starting state through sigmoid function--------
                        Vc_signal_ind = Mfc_c(dispImg, argmax);
                    else
                        % after the 1st item retrieval, the following items are
                        % retrieved by the transition matrix
                        % ------the alternative values are selected from the remaining unselected items------
                        dispImg_ic = dispImg(dispImg_choose ~= 1);
                        if isempty(tgt_choice) % ------This will only occur in the reconstruction reports------
                            % ----------The tgt_choice might be [] between if we use the
                            % dara from reconstruction reports----------
                            % added by rxj @ Sep 8 2024
                            % ------------------------------------------------------------
                            Vc_signal_ind = zeros(length(dispImg_ic), 1);
                        else
                            Vc_signal_ind = (Mc(dispImg(tgt_choice), dispImg_ic))';
                        end
                    end
                    % ------choice probability from the independent learning------
                    p_denominator_ind = exp(beta * Vc_signal_ind);
                    p_denominator_ind(isinf(p_denominator_ind)) = exp(700);
                    pCorr_ind = p_denominator_ind ./ sum(p_denominator_ind);

                    % ------choice probability from the binding learning------
                    Vc_signal_bid = Mb(dispImg_ic, position_i(ic));
                    p_denominator_bid = exp(beta * Vc_signal_bid);
                    p_denominator_bid(isinf(p_denominator_bid)) = exp(700);
                    pCorr_bid = p_denominator_bid ./ sum(p_denominator_bid);

                    % ------weighted average------
                    pCorr = omega_c * pCorr_ind + (1 - omega_c) * pCorr_bid;
                    pCorr(pCorr < 1e-16) = 1e-16;

                    %% calculate the likelihood
                    tgt_choice = find(conRep_iT == ic); % indicating participant made a decision
                    if ~isempty(tgt_choice)
                        dtrChoiceID = zeros(length(dispImg), 1);
                        dtrChoiceID(tgt_choice) = 1; % from data
                        dtrChoiceID(dispImg_choose == 1) = []; % delete already selected items

                        minuslli = minuslli + sum(dtrChoiceID .* log(pCorr));
                        dispImg_choose(tgt_choice) = 1; % label the already selected items
                        %lliTrial_all = [lliTrial_all; sum(dtrChoiceID .* log(pCorr))];
                    end
                end
            end
            
            %% --------position report--------
            if ~isempty(posRep_iT)
                dispPos_choose = zeros(nTrans + 1, 1); % one distractor
                for ip = 1 : nTrans
                    if ip == 1 % Initial item retrieval is based on the item-to-context strength
                        context_signal = sum(Mfc_p(dispPos, :), 1);
                        [argvalue, argmax] = max(context_signal);
                        % --------retrieve the starting state through sigmoid function--------
                        Vp_signal = Mfc_p(dispPos, argmax);
                    else
                        % after the 1st item retrieval, the following items are
                        % retrieved by the transition matrix
                        % ------the alternative values are selected from the remaining unselected items------
                        dispPos_ic = dispPos(dispPos_choose ~= 1);
                        Vp_signal = (Mp(dispPos(tgt_choice), dispPos_ic))';
                    end
                    p_denominator = exp(beta * Vp_signal);
                    p_denominator(isinf(p_denominator)) = exp(700);
                    pCorr = p_denominator ./ sum(p_denominator);
                    pCorr(pCorr < 1e-16)   = 1e-16;

                    %% calculate the likelihood
                    tgt_choice = find(posRep_iT == ip);
                    if ~isempty(tgt_choice)
                        dtrChoiceID = zeros(length(dispPos), 1);
                        dtrChoiceID(tgt_choice) = 1; % from data
                        dtrChoiceID(dispPos_choose == 1) = []; % delete already selected items

                        minuslli = minuslli + sum(dtrChoiceID .* log(pCorr));
                        dispPos_choose(tgt_choice) = 1; % label the already selected items
                        %lliTrial_all = [lliTrial_all; sum(dtrChoiceID .* log(pCorr))];
                    end
                end
            end
        end
    end
end
minuslli = -minuslli;

% Here is a hidden hypothesis that there would be no error in the other
% dimension when participants employ the binding association strategy

