% example.m
% basic code for fitting a learning model
clear all
close all
cd 'T:\RLmodel\rlfit-master'

%% Set period of interest for pupil diamter
% Note: might want to recheck literature to find principled duration of
% interest.  I think some papers even suggest that slightly before the
% stimulus is more telling than after stimulus.
pars.velThresh = 15;
pars.pkAdjust = 50; % amount to shift from detected pks to find beginning/end of peak
pars.antRespLength = 500; % for measure of pupil dilation, period to ignore prior to choice (right before choice there's typically an anticipatory response)
pars.grabAnticResp=false;
condition = 0;
plotPupilRespOverTrials = false; % toggle exploratory visualization
for ppInd = 1:45
    try
        load(['T:\RLmodel\dataset-eyeTracker-choices\two_stage_task_sub' num2str(ppInd) '_bars' num2str(condition) '.mat'])
    catch
        disp(['No data for pp' num2str(ppInd)])
    end
    
    %% Code pupil resps and plot measures across trials
    if plotPupilRespOverTrials
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    numRows = 3;
    %% Measure and plot post-stage onset pupil resp over trials
    for gatherStimOnset = [0]%,1]
        pars.gatherStimOnset = gatherStimOnset;
        if gatherStimOnset
            measOpts = [2000];%500,1000,2000,-500,-1000,-2000];
        else
            measOpts = [-2500];%,1000,2000,-500,-1000,-2000];
        end
        numCols = length(measOpts); % column for selected measurement period

        for measOptInd = 1:length(measOpts)
            pars.stimDurForEyeMeans = measOpts(measOptInd); % one sec
            pars.gatherStimOnset = false; % false = gather data leading to choice instead
            pupilDiams = codePupilResps(two_stage_task_data,condition,pars);
            if plotPupilRespOverTrials
                for state = 1:3
                    row = state;
                    col = measOptInd;
                    plotInd = (numCols*row)-(numCols-col);
                    subtightplot(3,length(measOpts),plotInd)
                    hold on
                    scatter([1:size(two_stage_task_data.eton,1)],pupilDiams(:,state))
                    xlim([0,size(two_stage_task_data.eton,1)])
                    yMin = 0;
                    yMax = 5000;
                    if min(pupilDiams(:,state)) < yMin
                        keyboard
                    end
                    if max(pupilDiams(:,state)) > yMax
                        keyboard
                    end
                end
            end
        end
    end
    ppData.pupilDiams{ppInd} = pupilDiams;
    %% convert to 'choice' var and 'outcome' var
    % outcome is easy
    for row = 1:size(two_stage_task_data.payoff,1)
        rowData = two_stage_task_data.payoff(row,:);
        rowData = rowData(~isnan(rowData));
        outcome(row) = sum(rowData);
        if isnan(outcome(row))
            keyboard
        end
        % choices have to be coded
        choiceOpts = two_stage_task_data.action(row,:);
        if choiceOpts(2) == 0
            action = 1;
        elseif choiceOpts(2) == 1
            action = 2;
        elseif choiceOpts(3) == 0
            action = 3;
        elseif choiceOpts(3) == 1
            action = 4;
        else
            keyboard
        end    
        choice(row) = action;
    end
    blueStateOnlyIndices = find(choice<3);
    outcome = outcome(blueStateOnlyIndices);
    choice = choice(blueStateOnlyIndices);
    pupilDiamsForSingleTrialAttentionWeights = pupilDiams(blueStateOnlyIndices,2);
    attentionWeights = zeros(length(pupilDiamsForSingleTrialAttentionWeights),1);
    %% BIC gets down to 82 with below exponentially weighted system (limits based on rough distribution estimates) 
    attentionWeights(pupilDiamsForSingleTrialAttentionWeights < 2020) = .1;
    attentionWeights(pupilDiamsForSingleTrialAttentionWeights < 2070 & attentionWeights == 0) = .2;
    attentionWeights(pupilDiamsForSingleTrialAttentionWeights < 2120 & attentionWeights == 0) = .4;
    attentionWeights(pupilDiamsForSingleTrialAttentionWeights > 2120 & attentionWeights == 0) = .8;
%     % 2020 -> 2070 -> 2120 (pre-choice)
    %% BIC gets down to 84.92 (normalized.^.5); 
%     pupilDiamsForSingleTrialAttentionWeights = (pupilDiamsForSingleTrialAttentionWeights-min(pupilDiamsForSingleTrialAttentionWeights))/(max(pupilDiamsForSingleTrialAttentionWeights)-min(pupilDiamsForSingleTrialAttentionWeights));
%     pupilDiamsForSingleTrialAttentionWeights = pupilDiamsForSingleTrialAttentionWeights.^.5;
    %% Add four params for setting limits on attention weights
%     attentionWeights = pupilDiamsForSingleTrialAttentionWeights;

    %% Begin model fitting section
    % model has only one parameter, the learning rate, between 0 and 1
    % in general, there will be one of these for each parameter (excluding the
    % softmax parameter); softmax parameter is the inverse temp if I remember
    % correctly
%     lb = [0,1900,1900,1900,1900]; %lower bounds
%     ub = [1,3000,3000,3000,3000]; %upper bounds
%     lb = [0,1000,0];%,1000,0,0]; %lower bounds
%     ub = [1,3500,1];%,3500,1,1]; %upper bounds
    lb = [0];%,1000,0];%,1000,0,0]; %lower bounds
    ub = [1];%,3500,1];%,3500,1,1]; %upper bounds

    
    % % however, we can also use a decorator to add perseveration behavior 
    % % to the model:
    % % the following are limits for a perseveration bonus
    % lb = [-1, lb];
    % ub = [1, ub];
    % Qfun = add_perseveration(@Q_model);
    Qfun = (@Q_model);

    % now optmize to fit model
    numiter = 15;
    [beta, LL, Q] = rlfit(Qfun, choice, outcome,attentionWeights,lb, ub, numiter);
    ppData.params{ppInd} = beta;
    ppData.LL{ppInd} = LL;
    ppData.Q{ppInd} = Q;

%     keyboard
%     % plot results
%     figure
%     plot(Q)
    match = 0;
    exploratoryTrialsBoolean_diffToMax = zeros(length(Q),1);
    avgQacrossTrials = mean(Q,1);
    for act = 1:size(Q,1)
        qRow = Q(act,:);
        % only compare available options once in final state (decode
        % options)
        % Recall: 
        % action 1 leftImage of leftChoices (blue)
        % action 2 rightImage of leftChoices (blue)
        % action 3 leftImage of rightChoices (purple)
        % action 4 rightImage of leftChoices (purple)
        actionData = two_stage_task_data.action(blueStateOnlyIndices(act),:);
%         if isnan(actionData(2))
%             % purple state
%             actionOpts = [1,2];
%         else
%             % blue state
%             actionOpts = [3,4];
%         end
%         % update qRow to contain relevat options
%         qRow = qRow(actionOpts);
%         if isempty(find(choice(act) == qRow))
%             disp('somethings broken')
%             disp(actionOpts)
%             disp(num2str(choice(act)))
%         end
        [maxQ,maxChoice] = max(qRow);
        if maxChoice == choice(act)
            match = match + 1;
            exploratoryTrialsBoolean_diffToMax(2,act) = nan;
        else
            diffToMax = maxQ - qRow(choice(act));
            exploratoryTrialsBoolean_diffToMax(1,act) = 1;
            exploratoryTrialsBoolean_diffToMax(2,act) = diffToMax;
        end
    end
    acc = match/length(choice)
    ppData.acc{ppInd} = acc;
    ppData.exploratoryTrialsBoolean_diffToMax{ppInd} = exploratoryTrialsBoolean_diffToMax;
end

diffToMaxes=[];
bin1Thresh=.1;
bin2Thresh=.3;
bin3Thresh=inf;
for ppInd = 1:length(ppData.pupilDiams)
    pupilData = ppData.pupilDiams{ppInd};
    relRows = find(~isnan(pupilData(:,2)));
    pupilData = pupilData(relRows,2);
    %     bar([1,2],[exploratoryPupilSize,nonExplorPupilSize])
    % combine all diffToMax's
    diffToMax=ppData.exploratoryTrialsBoolean_diffToMax{ppInd};
    diffToMax=diffToMax(1:2,diffToMax(1,:)==1);
    bin1Inds = find(diffToMax(2,:) < bin1Thresh & ~isnan(diffToMax(2,:)));
    bin2Inds = find(diffToMax(2,:) < bin2Thresh & diffToMax(2,:) >= bin1Thresh);
    bin3Inds = find(diffToMax(2,:) >= bin2Thresh);
    explorTrials = ppData.exploratoryTrialsBoolean_diffToMax{ppInd};
    explorTrials=explorTrials(1,:);
    exploratoryPupilSize_bin1 = mean(pupilData(bin1Inds));
    exploratoryPupilSize_bin2 = mean(pupilData(bin2Inds));
    exploratoryPupilSize_bin3 = mean(pupilData(bin3Inds));
%     nonExplorPupilSize = mean(pupilData(explorTrials==0));

    dataAcrossPps(ppInd,1) = exploratoryPupilSize_bin1;
    dataAcrossPps(ppInd,2) = exploratoryPupilSize_bin2;
    dataAcrossPps(ppInd,3) = exploratoryPupilSize_bin3;

    if isempty(diffToMaxes)
        diffToMaxes = diffToMax(2,:);
    else
        diffToMaxes=[diffToMaxes,diffToMax(2,:)];
    end
end
%% Plot pupil size according to bins/data distribution
figure
avgPupilSize_bin1 = mean(dataAcrossPps(:,1),'omitnan');
avgPupilSize_bin2 = mean(dataAcrossPps(:,2),'omitnan');
avgPupilSize_bin3 = mean(dataAcrossPps(:,3),'omitnan');
bar([1,2,3],[avgPupilSize_bin1,avgPupilSize_bin2,avgPupilSize_bin3])
ylim([1900,2200])
title('Pupil size Increases Towards More Exploratory (maxQ - chosenQ) Choices')
% 1995 -> 2040 -> 2090 (post-stim)
% 2020 -> 2070 -> 2120 (pre-choice)

%% Plot exploratory param val vs. avg pupil diamter across all trials
for paramInd = 1:length(ppData.params{1})
    figure
    exploreVsDiam = [];
    for ppInd = 1:length(ppData.pupilDiams)
        pupilDiams = ppData.pupilDiams{ppInd};
        pupilDiams = pupilDiams(:,2); % blue state only
        avgDiam = mean(pupilDiams,'omitnan');
        exploreVal = ppData.params{ppInd};
        exploreVal = exploreVal(paramInd);
        exploreVsDiam(ppInd,1) = exploreVal;
        exploreVsDiam(ppInd,2) = avgDiam;
    end
    scatter(exploreVsDiam(:,1),exploreVsDiam(:,2))
    title(['param' num2str(paramInd)])
    mdl = fitlm(exploreVsDiam(:,1),exploreVsDiam(:,2))
end
    
%% Plot distribution of Q-diffs
figure
hist(diffToMaxes,100)

%% Calc Sum BIC values
sumBIC = 0;
for ppInd = 1:length(ppData.acc)
    BIC = (2*log(length(Q))*length(ppData.params{ppInd})) - (2*ppData.LL{ppInd});
    sumBIC = sumBIC + BIC;
end
disp(['sumBIC = ' num2str(sumBIC)])

