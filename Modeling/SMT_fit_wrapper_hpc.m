function SMT_fit_wrapper_hpc(...
    folder, bhvDataDir, StudyWord_ii, groupName, subj_list, subj_conds, iGrp, iSub, expList, params)

% Contains the inner fitting loop (single subject)
nImg   = params.nImg;
nPos   = params.nPos;
nTrans = params.nTrans;
nEpi   = params.nEpi;
postTn = params.postTn;
nDtr   = params.nDtr;
suffixWord_young_list = params.suffixWord_young_list;
suffixWord_old_list   = params.suffixWord_old_list;

%% Preparing the data and fitting the models
if isequal(StudyWord_ii, 'Experiment1')
    %% model fitting for Experiment 1 (behavior study)
    % added by XR @ Oct 9th 2025
    StudyName = 'Exp1';
    MList       = {'featureIndependent', 'featureIndependentVM', 'featureBindingWeightRW-cOp', 'featureBindingWeightRW', 'featureCompetitionRW'};
    fitWordList = {'allLearning-marginalRep', 'allLearning-marginalRep', 'allLearning-marginalRep', 'allLearning-marginalRep', 'allLearning-marginalRep'};
    bindingDirecList = {'', '', 'cOp', 'pOc', ''};

    optimzerIdx = 0; % 0-fminsearchbnd; 1-bads
    dataFittingSource = 0; % 0-data from the marginal reports; 1-data from the reconstruction reports.
    states = nImg;
    refit  = 1;
    nFit = 100;

    suffixWord = expList{1};
    postTn_all = postTn;
    groupID  = groupName;
    subLen   = length(subj_list);
    expID    = suffixWord;
    %%
    disp(['======', groupName, '-', expID, '-sub', num2str(iSub), '======']);
    subjBv = subj_list{iSub, 1};
    subjTm = subj_list{iSub, 2};
    if iGrp == 1
        pathWord = suffixWord_young_list{iSub};
    elseif iGrp == 2
        pathWord = suffixWord_old_list{iSub};
    end
    subjPath = [bhvDataDir, '/AgingStudy-FlexibleBinding/', StudyName, '-raw/'];
    if isequal(pathWord, 'pathA')
        seqMem_subj = readtable([subjPath, subjBv, '_EpisodicMemoryTask_', subjTm, '.csv']);
    else
        seqMem_subj = readtable([subjPath, subjBv, '_EpisodicMemoryTask-', pathWord, '_', subjTm, '.csv']);
    end

    subID  = ['sub', num2str(iSub)];

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

    %%% no reconstruction-only trials
    reconsOnly = zeros(nEpi, 1);

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

    %% ---------- reconstruction report ----------
    bothTrue_tmp = seqMem_subj.bothReportTrue;
    bothTrue_tmp = bothTrue_tmp(~cellfun('isempty', bothTrue_tmp));
    bothRep_tmp  = seqMem_subj.bothReportOrd;
    bothRep_tmp  = bothRep_tmp(~cellfun('isempty', bothRep_tmp));
    %%% accuracy and RT calculation based on integration of content and
    %%% position
    bothTrue_col = cell(nEpi + postTn, 1);
    bothRep_col  = cell(nEpi + postTn, 1);
    for i =  1 : (nEpi + postTn)
        bothTrue_col{i} = str2num(bothTrue_tmp{i});
        bothCol = str2num(bothRep_tmp{i});
        if i <= nEpi
            bothCol = reshape(bothCol, 2, (nTrans+nDtr));
        else
            bothCol = reshape(bothCol, 2, nTrans);
        end
        bothRep_col{i}  = bothCol;
    end

    %% ---------- content report ----------
    conTrue_tmp = seqMem_subj.conReportTrue;
    conTrue_tmp = conTrue_tmp(~cellfun('isempty', conTrue_tmp));
    conRep_tmp  = seqMem_subj.conReportOrd;
    conRep_tmp  = conRep_tmp(~cellfun('isempty', conRep_tmp));
    conTrue_col = cell(nEpi, 1);
    conRep_col  = cell(nEpi, 1);
    for i =  1 : nEpi
        conTrue_col{i} = str2num(conTrue_tmp{i});
        conRep_noRef   = str2num(conRep_tmp{i});
        conRep_col{i}  = conRep_noRef;
    end

    %% ---------- position report ----------
    locTrue_tmp = seqMem_subj.locReportTrue;
    locTrue_tmp = locTrue_tmp(~cellfun('isempty', locTrue_tmp));
    locRep_tmp  = seqMem_subj.locReportOrd;
    locRep_tmp  = locRep_tmp(~cellfun('isempty', locRep_tmp));
    locTrue_col = cell(nEpi, 1);
    locRep_col  = cell(nEpi, 1);
    for i =  1 : nEpi
        locTrue_col{i} = str2num(locTrue_tmp{i});
        locRep_noRef   = str2num(locRep_tmp{i});
        locRep_col{i}  = locRep_noRef;
    end

    %% ========Prepararing data for model fitting========
    % index to indicate transitions in each sequence
    % 8 indexes ranging from 0 to 8 for each dimension
    seqLab = nan(nEpi, 1);
    for i =  1 : nEpi
        uniC_i = uniCombSeq(i, :);
        if uniC_i(1) == 0 && uniC_i(2) == 0     % first number: content; second number: position
            seqLab(i) = 1;
        elseif uniC_i(1) == 0 && uniC_i(2) == 1
            seqLab(i) = 2;
        elseif uniC_i(1) == 1 && uniC_i(2) == 0
            seqLab(i) = 3;
        elseif uniC_i(1) == 1 && uniC_i(2) == 1
            seqLab(i) = 4;
        end
    end

    % ================the two content transition sequences================
    % extracting the content string from the displaying period of
    % reporting
    conSeqAll = seqMem_subj.conSeqTrl;
    conSeqAll = conSeqAll(~cellfun('isempty', conSeqAll));
    contentString = cell(nEpi, nTrans + 1);
    for i =  1 : nEpi
        input_str = conSeqAll{i};
        file_paths = strsplit(input_str, ',');
        % Initialize an empty cell array to store the extracted names
        names = cell(size(file_paths));

        % Iterate over each file path
        for iS = 1:length(file_paths)
            % Remove leading and trailing whitespaces and quotation marks
            file_path = strtrim(file_paths{iS});
            file_path = erase(file_path, '"');

            % Split the file path into path, name, and extension
            [~, file_name, ~] = fileparts(file_path);

            % Split the file name at '/'
            parts = strsplit(file_name, '/');

            % Extract the desired name
            name = parts{end};

            % Store the name in the cell array
            names{iS} = name;
        end
        contentString(i, :) = names;
    end
    % trial index of the two content sequences
    conA_id = find(uniCombSeq(:, 1) == 0); % one of the content transition sequence
    conB_id = find(uniCombSeq(:, 1) == 1); % the other content transition sequence
    conAB_string = cell(nTrans, 2);

    % write down the number corresponding to each string
    itemIdxA = conTrue_col{conA_id(1)};
    [~, itemIdxA_sort] = sort(itemIdxA);
    conStrA = contentString(conA_id(1), itemIdxA_sort);
    conStrA(end) = []; % the last item is the distractor
    conAB_string(:, 1) = conStrA;

    itemIdxB = conTrue_col{conB_id(1)};
    [~, itemIdxB_sort] = sort(itemIdxB);
    conStrB = contentString(conB_id(1), itemIdxB_sort);
    conStrB(end) = [];
    conAB_string(:, 2) = conStrB;

    % the complete image list, corresponding to index from 1 to 8
    imgList = [conAB_string(:, 1); conAB_string([2, 4, 5], 2)];
    conA_seq = (1 : 1 : nTrans);
    conB_seq = [3, 6, 1, 7, 8];

    % --------Assign the correct index for the stimuli displaying
    % during the report--------
    conStimIdx_report = cell(nEpi, 1); % 5 targets + 1 distractor; the index of the corresponding stimilus in the whole image list (related to the update of the transition matrix)
    % the corresponding correct responding order is save in
    % !!!!!! conTrue_col.mat !!!!!!
    conSeq_encode = cell(nEpi, 1); % 5 targets; sequence displayed during encoding stage
    for i = 1 : nEpi
        conRep_str = contentString(i, :);
        % find the index of those displayed images in the whole image
        % list
        conStimIdx_i = nan(1, nTrans + 1);
        for id = 1 : (nTrans + 1)
            index = find(ismember(imgList, conRep_str{id}));
            conStimIdx_i(id) = index;
        end
        conStimIdx_report{i} = conStimIdx_i;
        % the corresponding displayed sequence during encoding stage
        itemIdx_i = conTrue_col{i};
        conStimIdx_i(itemIdx_i > nTrans) = [];
        if sum(sort(conStimIdx_i) == conA_seq) == nTrans
            conSeq_encode{i} = conA_seq;
        else
            conSeq_encode{i} = conB_seq;
        end
    end

    % ================the two position transition sequences================
    positionTrans = zeros(nTrans, 2);
    posA_id = find(uniCombSeq(:, 2) == 0); % one of the position transition sequence
    posB_id = find(uniCombSeq(:, 2) == 1); % the other position transition sequence
    % positions of one of the position transition sequence
    posA_seq = zeros(nTrans, 2);
    posA_seq(:, 1) = posX_col{posA_id(1)}(1 : nTrans);
    posA_seq(:, 2) = posY_col{posA_id(1)}(1 : nTrans);
    % positions of the other position transition sequence
    posB_seq = zeros(nTrans, 2);
    posB_seq(:, 1) = posX_col{posB_id(1)}(1 : nTrans);
    posB_seq(:, 2) = posY_col{posB_id(1)}(1 : nTrans);
    % ------assign index to the two position transition sequences------
    if isequal(pathWord, 'pathA')
        % [1, 3, 8, 6, 4] or [6, 2, 5, 7, 3]
        if (posA_seq(2, 1) == posB_seq(nTrans, 1) && posA_seq(2, 2) == posB_seq(nTrans, 2)) &&...
                (posA_seq(4, 1) == posB_seq(1, 1) && posA_seq(4, 2) == posB_seq(1, 2))
            positionTrans = [1, 3, 8, 6, 4; ...
                             6, 2, 5, 7, 3];

        elseif (posB_seq(2, 1) == posA_seq(nTrans, 1) && posB_seq(2, 2) == posA_seq(nTrans, 2)) &&...
                (posB_seq(4, 1) == posA_seq(1, 1) && posB_seq(4, 2) == posA_seq(1, 2))
            positionTrans = [6, 2, 5, 7, 3; ...
                             1, 3, 8, 6, 4];
        end

    elseif isequal(pathWord, 'pathB')
        % [1, 3, 8, 6, 4] or [4, 7, 5, 8, 2]
        if (posA_seq(3, 1) == posB_seq(4, 1) && posA_seq(3, 2) == posB_seq(4, 2)) &&...
                (posA_seq(nTrans, 1) == posB_seq(1, 1) && posA_seq(nTrans, 2) == posB_seq(1, 2))
            positionTrans = [1, 3, 8, 6, 4; ...
                             4, 7, 5, 8, 2];

        elseif (posB_seq(3, 1) == posA_seq(4, 1) && posB_seq(3, 2) == posA_seq(4, 2)) &&...
                (posB_seq(nTrans, 1) == posA_seq(1, 1) && posB_seq(nTrans, 2) == posA_seq(1, 2))
            positionTrans = [4, 7, 5, 8, 2; ...
                             1, 3, 8, 6, 4];
        end

    elseif isequal(pathWord, 'pathC')
        % [1, 3, 8, 6, 4] or [4, 7, 2, 5, 8]
        if (posA_seq(3, 1) == posB_seq(nTrans, 1) && posA_seq(3, 2) == posB_seq(nTrans, 2)) &&...
                (posA_seq(nTrans, 1) == posB_seq(1, 1) && posA_seq(nTrans, 2) == posB_seq(1, 2))
            positionTrans = [1, 3, 8, 6, 4; ...
                             4, 7, 2, 5, 8];

        elseif (posB_seq(3, 1) == posA_seq(nTrans, 1) && posB_seq(3, 2) == posA_seq(nTrans, 2)) &&...
                (posB_seq(nTrans, 1) == posA_seq(1, 1) && posB_seq(nTrans, 2) == posA_seq(1, 2))
            positionTrans = [4, 7, 2, 5, 8; ...
                             1, 3, 8, 6, 4];
        end

    end

    posList = nan(nPos, 2); % 2: x- and y-axis
    for ip = 1 : nPos
        ip_find = find(positionTrans(1, :) == ip);
        if ~isempty(ip_find)
            posList(ip, :) = posA_seq(ip_find, :);
        else
            ip_find = find(positionTrans(2, :) == ip);
            posList(ip, :) = posB_seq(ip_find, :);
        end
    end
    % --------convert posList to angles--------
    angList = nan(nPos, 1);
    for ip = 1 : nPos
        angList(ip) = atan2(posList(ip, 2), posList(ip, 1)); % the angle in radians between -π and π
    end
    % angle transitions in each sequence
    angList_seq = zeros(nTrans, 2);
    angList_seq(:, 1) = angList((positionTrans(1, :))');
    angList_seq(:, 2) = angList((positionTrans(2, :))');

    % --------Assign the correct index for the stimuli displaying
    % during the report--------
    posStimIdx_report = cell(nEpi, 1); % the index of the corresponding stimilus in the whole angle or position list (related to the update of the transition matrix)
    % the corresponding correct responding order is save in
    % !!!!!! locTrue_col.mat !!!!!!
    posSeq_encode = cell(nEpi, 1); % 5 targets; sequence displayed during encoding stage
    angSeq_encode = cell(nEpi, 1); % angles for the five targets during encoding stage
    for i = 1 : nEpi
        posX_col_i = posX_col{i};
        posY_col_i = posY_col{i};
        % find the index of those displayed images in the whole image
        % list
        posStimIdx_i = nan(1, nTrans + 1);
        for id = 1 : (nTrans + 1)
            index_X = find(ismember(posList(:, 1), posX_col_i(id))); % x-pos
            index_Y = find(ismember(posList(:, 2), posY_col_i(id))); % y-pos
            index = intersect(index_X, index_Y);
            posStimIdx_i(id) = index;
        end
        posStimIdx_report{i} = posStimIdx_i;

        % the corresponding displayed sequence during encoding stage
        itemIdx_i = locTrue_col{i};
        posStimIdx_i(itemIdx_i > nTrans) = []; % always deleting the last item
        if sum(sort(posStimIdx_i) == sort(positionTrans(1, :))) == nTrans
            posSeq_encode{i} = positionTrans(1, :);
        else
            posSeq_encode{i} = positionTrans(2, :);
        end
        angSeq_encode{i} = angList((posSeq_encode{i})');
    end

    %% ----------Reorganizing stimuli for subsequent fitting----------
    % Only the non-reconstruction-only trials
    stim_encode = cell(nEpi, 2); % 2: content and position sequence, 64-16=48 trials
    stim_encode(:, 1) = conSeq_encode; % conSeq_encode(reconsOnly == 0);
    stim_encode(:, 2) = posSeq_encode;

    disp_con = conStimIdx_report;
    disp_pos = posStimIdx_report;
    disp_rec = [];

    if dataFittingSource == 0     % responses from the marginal reports
        resp_con = conRep_col;
        resp_pos = locRep_col;
    elseif dataFittingSource == 1 % responses from the reconstruction reports
        resp_con = cell(nEpi, 1);
        resp_pos = cell(nEpi, 1);
        for i = 1 : nEpi
            resp_con{i} = bothRep_col{i}(1, :);
            resp_pos{i} = bothRep_col{i}(2, :);
        end
        reconsOnly = zeros(nEpi, 1);
    end
    resp_rec = [];

    context_arr = seqLab;
    testOrd_arr = testOrd;

    trialLen = length(context_arr);

    %% fitting
    for iM = 1 : length(MList)
        Midx = MList{iM};
        fitWord = fitWordList{iM};
        bindingDirec = bindingDirecList{iM};

        pattern = '-(\w+)$'; % This pattern captures the word at the end of the string after the last hyphen
        matches = regexp(fitWord, pattern, 'tokens', 'once');
        if ~isempty(matches)
            uniPrior = matches{1};
        else
            uniPrior = '';
        end
        if isequal(uniPrior, 'uniPrior')
            priorWord = 'uniform'; % {'diagonal', 'uniform'}
        else
            priorWord = 'diagonal'; % {'diagonal', 'uniform'}
        end

        % initial matrix for the learning stage
        if isequal(priorWord, 'diagonal')
            % ------the initial matrix------
            Minit = eye(states, states); % 8 * 8; eye(states, states);
        elseif isequal(priorWord, 'uniform')
            Minit = zeros(states, states);
            Minit(:) = 1 / states;
        end
        [lliEst, paramsEst, lli_all] = SMT_fitting(folder, groupID, expID, Midx, subID, fitWord, refit, nFit, nImg, nPos, nTrans, trialLen, ...
            stim_encode, disp_con, disp_pos, disp_rec, resp_con, resp_pos, resp_rec, context_arr, testOrd_arr, reconsOnly, Minit, angList, angSeq_encode, ...
            bindingDirec, optimzerIdx, dataFittingSource);

    end

elseif isequal(StudyWord_ii, 'Experiment2Interleaved') 
    StudyName   = 'Exp2';
    MList       = {'featureIndependent', 'featureIndependentVM', 'featureBindingWeightRW-cOp', 'featureBindingWeightRW', 'featureCompetitionRW'};
    fitWordList = {'allLearning-marginalRep', 'allLearning-marginalRep', 'allLearning-marginalRep', 'allLearning-marginalRep', 'allLearning-marginalRep'};
    bindingDirecList = {'', '', 'cOp', 'pOc', ''};

    optimzerIdx = 0; % 0-fminsearchbnd; 1-bads
    dataFittingSource = 0; % 0-data from the marginal reports; 1-data from the reconstruction reports.
    states = nImg;
    refit  = 1;
    nFit = 100;

    suffixWord = expList{1};
    postTn_all = postTn;
    groupID  = groupName;
    subLen   = length(subj_list);
    expID    = suffixWord;

    %%
    disp(['======', groupName, '-', expID, '-sub', num2str(iSub), '======']);
    subjBv = subj_list{iSub, 1};
    subjTm = subj_list{iSub, 2};
    subjPath    = [bhvDataDir, '/AgingStudy-FlexibleBinding/', StudyName, '-raw/'];
    seqMem_subj = readtable([subjPath, subjBv, '_EpisodicMemoryTask-', suffixWord, '_', subjTm, '.csv']);
    subID  = ['sub', num2str(iSub)];

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

    %% ---------- reconstruction report ----------
    bothTrue_tmp = seqMem_subj.bothReportTrue;
    bothTrue_tmp = bothTrue_tmp(~cellfun('isempty', bothTrue_tmp));
    bothRep_tmp  = seqMem_subj.bothReportOrd;
    bothRep_tmp  = bothRep_tmp(~cellfun('isempty', bothRep_tmp));
    %%% accuracy and RT calculation based on integration of content and
    %%% position
    bothTrue_col = cell(nEpi + postTn, 1);
    bothRep_col  = cell(nEpi + postTn, 1);
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
    end

    %% ---------- content report ----------
    conTrue_tmp_raw = seqMem_subj.conReportTrue;
    conTrue_tmp_raw = conTrue_tmp_raw(~cellfun('isempty', conTrue_tmp_raw));
    conTrue_tmp = cell(nEpi, 1);
    conTrue_tmp(reconsOnly == 0) = conTrue_tmp_raw;
    conRep_tmp_raw  = seqMem_subj.conReportOrd;
    conRep_tmp_raw  = conRep_tmp_raw(~cellfun('isempty', conRep_tmp_raw));
    conRep_tmp = cell(nEpi, 1);
    conRep_tmp(reconsOnly == 0) = conRep_tmp_raw;
    conTrue_col = cell(nEpi, 1);
    conRep_col  = cell(nEpi, 1);
    for i =  1 : nEpi
        if reconsOnly(i) == 0 %% non reconstruction only trial
            conTrue_col{i} = str2num(conTrue_tmp{i});
            conRep_noRef   = str2num(conRep_tmp{i});
            conRep_col{i}  = conRep_noRef;
        end
    end
    conTrue_col(reconsOnly == 1) = bothTrue_col(reconsOnly == 1);

    %% ---------- position report ----------
    locTrue_tmp_raw = seqMem_subj.locReportTrue;
    locTrue_tmp_raw = locTrue_tmp_raw(~cellfun('isempty', locTrue_tmp_raw));
    locTrue_tmp = cell(nEpi, 1);
    locTrue_tmp(reconsOnly == 0) = locTrue_tmp_raw;
    locRep_tmp_raw  = seqMem_subj.locReportOrd;
    locRep_tmp_raw  = locRep_tmp_raw(~cellfun('isempty', locRep_tmp_raw));
    locRep_tmp  = cell(nEpi, 1);
    locRep_tmp(reconsOnly == 0) = locRep_tmp_raw;
    locTrue_col = cell(nEpi, 1);
    locRep_col  = cell(nEpi, 1);
    for i =  1 : nEpi
        if reconsOnly(i) == 0 %% non reconstruction only trial
            locTrue_col{i} = str2num(locTrue_tmp{i});
            locRep_noRef   = str2num(locRep_tmp{i});
            locRep_col{i}  = locRep_noRef;
        end
    end
    locTrue_recon_trial = cell(sum(reconsOnly), 1);
    for i = 1 : sum(reconsOnly)
        locTrue_recon_trial{i} = 1 : 1 : (nTrans + 1);
    end
    locTrue_col(reconsOnly == 1) = locTrue_recon_trial;

    %% ---------- reconstruction report ----------
    bothTrue_tmp = seqMem_subj.bothReportTrue;
    bothTrue_tmp = bothTrue_tmp(~cellfun('isempty', bothTrue_tmp));
    bothRep_tmp  = seqMem_subj.bothReportOrd;
    bothRep_tmp  = bothRep_tmp(~cellfun('isempty', bothRep_tmp));
    %%% accuracy and RT calculation based on integration of content and
    %%% position
    bothTrue_col = cell(nEpi, 1); % only the non-postTest trials
    bothRep_col  = cell(nEpi, 1);
    for i =  1 : nEpi
        bothTrue_col{i} = str2num(bothTrue_tmp{i});
        bothCol = str2num(bothRep_tmp{i});
        bothCol = reshape(bothCol, 2, (nTrans+nDtr));
        bothRep_col{i}  = bothCol;
    end

    %% ========Prepararing data for model fitting========
    % index to indicate transitions in each sequence
    % 8 indexes ranging from 0 to 8 for each dimension
    seqLab = nan(nEpi, 1);
    for i =  1 : nEpi
        uniC_i = uniCombSeq(i, :);
        if uniC_i(1) == 0 && uniC_i(2) == 0     % first number: content; second number: position
            seqLab(i) = 1;
        elseif uniC_i(1) == 0 && uniC_i(2) == 1
            seqLab(i) = 2;
        elseif uniC_i(1) == 1 && uniC_i(2) == 0
            seqLab(i) = 3;
        elseif uniC_i(1) == 1 && uniC_i(2) == 1
            seqLab(i) = 4;
        end
    end

    % ================the two content transition sequences================
    % extracting the content string from the displaying period of
    % reporting
    conSeqAll = seqMem_subj.conSeqTrl;
    conSeqAll = conSeqAll(~cellfun('isempty', conSeqAll));
    contentString = cell(nEpi, nTrans + 1);
    for i =  1 : nEpi
        input_str = conSeqAll{i};
        file_paths = strsplit(input_str, ',');
        % Initialize an empty cell array to store the extracted names
        names = cell(size(file_paths));

        % Iterate over each file path
        for iS = 1:length(file_paths)
            % Remove leading and trailing whitespaces and quotation marks
            file_path = strtrim(file_paths{iS});
            file_path = erase(file_path, '"');

            % Split the file path into path, name, and extension
            [~, file_name, ~] = fileparts(file_path);

            % Split the file name at '/'
            parts = strsplit(file_name, '/');

            % Extract the desired name
            name = parts{end};

            % Store the name in the cell array
            names{iS} = name;
        end
        contentString(i, :) = names;
    end
    % trial index of the two content sequences
    conA_id = find(uniCombSeq(:, 1) == 0 & reconsOnly == 0); % one of the content transition sequence
    conB_id = find(uniCombSeq(:, 1) == 1 & reconsOnly == 0); % the other content transition sequence
    conAB_string = cell(nTrans, 2);

    % write down the number corresponding to each string
    itemIdxA = conTrue_col{conA_id(1)};
    [~, itemIdxA_sort] = sort(itemIdxA);
    conStrA = contentString(conA_id(1), itemIdxA_sort);
    conStrA(end) = []; % the last item is the distractor
    conAB_string(:, 1) = conStrA;

    itemIdxB = conTrue_col{conB_id(1)};
    [~, itemIdxB_sort] = sort(itemIdxB);
    conStrB = contentString(conB_id(1), itemIdxB_sort);
    conStrB(end) = [];
    conAB_string(:, 2) = conStrB;

    % the complete image list, corresponding to index from 1 to 8
    imgList = [conAB_string(:, 1); conAB_string([2, 4, 5], 2)];
    conA_seq = (1 : 1 : nTrans);
    conB_seq = [3, 6, 1, 7, 8];

    % --------Assign the correct index for the stimuli displaying
    % during the report--------
    conStimIdx_report = cell(nEpi, 1); % 5 targets + 1 distractor; the index of the corresponding stimilus in the whole image list (related to the update of the transition matrix)
    % the corresponding correct responding order is save in
    % !!!!!! conTrue_col.mat !!!!!!
    conSeq_encode = cell(nEpi, 1); % 5 targets; sequence displayed during encoding stage
    for i = 1 : nEpi
        conRep_str = contentString(i, :);
        % find the index of those displayed images in the whole image
        % list
        conStimIdx_i = nan(1, nTrans + 1);
        for id = 1 : (nTrans + 1)
            index = find(ismember(imgList, conRep_str{id}));
            conStimIdx_i(id) = index;
        end
        conStimIdx_report{i} = conStimIdx_i;
        % the corresponding displayed sequence during encoding stage
        itemIdx_i = conTrue_col{i};
        conStimIdx_i(itemIdx_i > nTrans) = [];
        if sum(sort(conStimIdx_i) == conA_seq) == nTrans
            conSeq_encode{i} = conA_seq;
        else
            conSeq_encode{i} = conB_seq;
        end
    end

    % ================the two position transition sequences================
    positionTrans = zeros(nTrans, 2);
    posA_id = find(uniCombSeq(:, 2) == 0 & reconsOnly == 0); % one of the position transition sequence
    posB_id = find(uniCombSeq(:, 2) == 1 & reconsOnly == 0); % the other position transition sequence
    % positions of one of the position transition sequence
    posA_seq = zeros(nTrans, 2);
    posA_seq(:, 1) = posX_col{posA_id(1)}(1 : nTrans);
    posA_seq(:, 2) = posY_col{posA_id(1)}(1 : nTrans);
    % positions of the other position transition sequence
    posB_seq = zeros(nTrans, 2);
    posB_seq(:, 1) = posX_col{posB_id(1)}(1 : nTrans);
    posB_seq(:, 2) = posY_col{posB_id(1)}(1 : nTrans);
    % ------assign index to the two position transition sequences------
    if (posA_seq(2, 1) == posB_seq(nTrans, 1) && posA_seq(2, 2) == posB_seq(nTrans, 2)) &&...
            (posA_seq(4, 1) == posB_seq(1, 1) && posA_seq(4, 2) == posB_seq(1, 2))
        positionTrans = [1, 3, 8, 6, 4; ...
            6, 2, 5, 7, 3];

    elseif (posB_seq(2, 1) == posA_seq(nTrans, 1) && posB_seq(2, 2) == posA_seq(nTrans, 2)) &&...
            (posB_seq(4, 1) == posA_seq(1, 1) && posB_seq(4, 2) == posA_seq(1, 2))
        positionTrans = [6, 2, 5, 7, 3; ...
            1, 3, 8, 6, 4];
    end
    posList = nan(nPos, 2); % 2: x- and y-axis
    for ip = 1 : nPos
        ip_find = find(positionTrans(1, :) == ip);
        if ~isempty(ip_find)
            posList(ip, :) = posA_seq(ip_find, :);
        else
            ip_find = find(positionTrans(2, :) == ip);
            posList(ip, :) = posB_seq(ip_find, :);
        end
    end
    % --------convert posList to angles--------
    angList = nan(nPos, 1);
    for ip = 1 : nPos
        angList(ip) = atan2(posList(ip, 2), posList(ip, 1)); % the angle in radians between -π and π
    end
    % angle transitions in each sequence
    angList_seq = zeros(nTrans, 2);
    angList_seq(:, 1) = angList((positionTrans(1, :))');
    angList_seq(:, 2) = angList((positionTrans(2, :))');

    % --------Assign the correct index for the stimuli displaying
    % during the report--------
    posStimIdx_report = cell(nEpi, 1); % the index of the corresponding stimilus in the whole angle or position list (related to the update of the transition matrix)
    % the corresponding correct responding order is save in
    % !!!!!! locTrue_col.mat !!!!!!
    posSeq_encode = cell(nEpi, 1); % 5 targets; sequence displayed during encoding stage
    angSeq_encode = cell(nEpi, 1); % angles for the five targets during encoding stage
    for i = 1 : nEpi
        posX_col_i = posX_col{i};
        posY_col_i = posY_col{i};
        % find the index of those displayed images in the whole image
        % list
        posStimIdx_i = nan(1, nTrans + 1);
        for id = 1 : (nTrans + 1)
            index_X = find(ismember(posList(:, 1), posX_col_i(id))); % x-pos
            index_Y = find(ismember(posList(:, 2), posY_col_i(id))); % y-pos
            index = intersect(index_X, index_Y);
            posStimIdx_i(id) = index;
        end
        posStimIdx_report{i} = posStimIdx_i;

        % the corresponding displayed sequence during encoding stage
        itemIdx_i = locTrue_col{i};
        posStimIdx_i(itemIdx_i > nTrans) = []; % always deleting the last item
        if sum(sort(posStimIdx_i) == sort(positionTrans(1, :))) == nTrans
            posSeq_encode{i} = positionTrans(1, :);
        else
            posSeq_encode{i} = positionTrans(2, :);
        end
        angSeq_encode{i} = angList((posSeq_encode{i})');
    end

    %% ----------Reorganizing stimuli for subsequent fitting----------
    % Only the non-reconstruction-only trials
    stim_encode = cell(nEpi, 2); % 2: content and position sequence, 64-16=48 trials
    stim_encode(:, 1) = conSeq_encode; % conSeq_encode(reconsOnly == 0);
    stim_encode(:, 2) = posSeq_encode;

    disp_con = conStimIdx_report;
    disp_pos = posStimIdx_report;
    disp_rec = [];

    if dataFittingSource == 0     % responses from the marginal reports
        resp_con = conRep_col;
        resp_pos = locRep_col;
    elseif dataFittingSource == 1 % responses from the reconstruction reports
        resp_con = cell(nEpi, 1);
        resp_pos = cell(nEpi, 1);
        for i = 1 : nEpi
            resp_con{i} = bothRep_col{i}(1, :);
            resp_pos{i} = bothRep_col{i}(2, :);
        end
        reconsOnly = zeros(nEpi, 1);
    end
    resp_rec = [];

    context_arr = seqLab;
    testOrd_arr = testOrd;

    trialLen = length(context_arr);

    %% fitting
    for iM = 1 : length(MList)
        Midx = MList{iM};
        fitWord = fitWordList{iM};
        bindingDirec = bindingDirecList{iM};

        pattern = '-(\w+)$'; % This pattern captures the word at the end of the string after the last hyphen
        matches = regexp(fitWord, pattern, 'tokens', 'once');
        if ~isempty(matches)
            uniPrior = matches{1};
        else
            uniPrior = '';
        end
        if isequal(uniPrior, 'uniPrior')
            priorWord = 'uniform'; % {'diagonal', 'uniform'}
        else
            priorWord = 'diagonal'; % {'diagonal', 'uniform'}
        end

        % initial matrix for the learning stage
        if isequal(priorWord, 'diagonal')
            % ------the initial matrix------
            Minit = eye(states, states); % 8 * 8; eye(states, states);
        elseif isequal(priorWord, 'uniform')
            Minit = zeros(states, states);
            Minit(:) = 1 / states;
        end
        [lliEst, paramsEst, lli_all] = SMT_fitting(folder, groupID, expID, Midx, subID, fitWord, refit, nFit, nImg, nPos, nTrans, trialLen, ...
            stim_encode, disp_con, disp_pos, disp_rec, resp_con, resp_pos, resp_rec, context_arr, testOrd_arr, reconsOnly, Minit, angList, angSeq_encode, ...
            bindingDirec, optimzerIdx, dataFittingSource);

    end

end