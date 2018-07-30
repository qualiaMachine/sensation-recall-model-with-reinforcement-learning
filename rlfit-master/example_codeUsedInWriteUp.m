% basic code for fitting a learning model
% % % clear all
clear all
close all

pars.dirs.baseDataDir='/Volumes/NEWUSB 1/data/sensation_recall_model_RL';
addpath('/Volumes/NEWUSB 1/code/3rdParty/subtightplot/subtightplot')
pars.dirs.figDir = '/Volumes/NEWUSB 1/analysisResults/rlModel/figs/';
pars.dirs.saveDir = '/Volumes/NEWUSB 1/analysisResults/rlModel/';

%% Set params 
plotPupilRespOverTrials = false; % toggle exploratory visualization
pars.grabAnticResp=0;
pars.velThresh = 15;
pars.pkAdjust = 50; % amount to shift from detected pks to find beginning/end of peak
pars.antRespLength = 500; % for measure of pupil dilation, period to ignore prior to choice (right before choice there's typically an anticipatory response)
condition = 0;

% skipCreateData = false;
% loadData = false;
% plotStuff = false;
% 
loadData = true;
plotStuff = true;

for trinaryPupilMeasSetting=[2];%,1,3]; % 3=subtract postDecision-preDecision; 2=preDecision, 1=postDecision
    %% Iterate through model options and save results
    for modelNum = [3]%6:9,2,3]
        disp(['model used = ' num2str(modelNum)])
        if ~loadData
            for ppInd = 1:45
                try
                    load([pars.dirs.baseDataDir filesep 'dataset-eyeTracker-choices' filesep 'two_stage_task_sub' num2str(ppInd) '_bars' num2str(condition) '.mat'])
                catch
                    disp(['No data for pp' num2str(ppInd)])
                    continue
                end

                %% Code pupil resps and plot measures across trials
                numRows = 3;
                %% Measure and plot post-stage onset pupil resp over trials
                if trinaryPupilMeasSetting == 1
                    pars.gatherStimOnset=1;
                    pars.stimDurForEyeMeans = 2500;%500,1000,2000,-500,-1000,-2000];
                    pars.grabAnticResp = 0;
                elseif trinaryPupilMeasSetting == 2
                    pars.gatherStimOnset=0;
                    if pars.grabAnticResp
                        pars.stimDurForEyeMeans = -500; % look specifically at anticipatory response
                    else
                        pars.stimDurForEyeMeans = -3000;%,1000,2000,-500,-1000,-2000];
                    end
                else
                    % measure pupilSize diff between postDecision and preDecison
                end

                pupil1=[];
                allDrops1=[];
                pupil2=[];
                allDrops2=[];
                if trinaryPupilMeasSetting == 3
                    pars.gatherStimOnset=0;
                    pars.stimDurForEyeMeans=-3000;
                    [pupilDiams,allPupDrops] = codePupilResps(two_stage_task_data,condition,pars);
                    pupil1=pupilDiams;
                    allDrops1=allPupDrops;

                    pars.gatherStimOnset=1;
                    pars.stimDurForEyeMeans=2500;
                    [pupilDiams,allPupDrops] = codePupilResps(two_stage_task_data,condition,pars);
                    pupil2=pupilDiams;
                    allDrops2=allPupDrops;
                    pupilDiams=pupil2-pupil1;
                    try
                        allPupDrops=allDrops2-allDrops1;
                    catch why
                        keyboard
                    end
                else
                    [pupilDiams,allPupDrops] = codePupilResps(two_stage_task_data,condition,pars);
                end

                if sum(sum(pupilDiams == 0)) > 3
                    keyboard
                end

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
                ppData.pupilDiams{ppInd} = pupilDiams;
                ppData.pupDrops{ppInd} = allPupDrops;
                if sum(sum(pupilDiams == 0)) > 3
                    keyboard
                end
                %% convert to 'choice' var and 'outcome' var
                action1 = two_stage_task_data.action(:,1);
                action2 = sum(two_stage_task_data.action(:,2:3),2,'omitnan');
                states = two_stage_task_data.state;

                % outcome is easy
                for row = 1:size(two_stage_task_data.payoff,1)
                    rowData = two_stage_task_data.payoff(row,:);
                    rowData = rowData(~isnan(rowData));
                    payoff(row) = sum(rowData,'omitnan');
                    if isnan(payoff(row))
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
            %         choice(row) = action;
                end
            %     blueStateOnlyIndices = find(choice<3);
            %     outcome = outcome(blueStateOnlyIndices);
            %     choice = choice(blueStateOnlyIndices);
                pupilDiamsForSingleTrialAttentionWeights = pupilDiams;%(blueStateOnlyIndices,2);
            %     figure;
            %     bar(pupilDiamsForSingleTrialAttentionWeights)
                attentionWeights = zeros(length(pupilDiamsForSingleTrialAttentionWeights),1);
                %% BIC gets down to 82 with below exponentially weighted system (limits based on rough distribution estimates) 
            %     attentionWeights(pupilDiamsForSingleTrialAttentionWeights < 2020) = .1;
            %     attentionWeights(pupilDiamsForSingleTrialAttentionWeights < 2070 & attentionWeights == 0) = .2;
            %     attentionWeights(pupilDiamsForSingleTrialAttentionWeights < 2120 & attentionWeights == 0) = .4;
            %     attentionWeights(pupilDiamsForSingleTrialAttentionWeights > 2120 & attentionWeights == 0) = .8;
            %     % 2020 -> 2070 -> 2120 (pre-choice)
                %% BIC gets down to 84.92 (normalized.^.5); 
            %     pupilDiamsForSingleTrialAttentionWeights = (pupilDiamsForSingleTrialAttentionWeights-min(pupilDiamsForSingleTrialAttentionWeights))/(max(pupilDiamsForSingleTrialAttentionWeights)-min(pupilDiamsForSingleTrialAttentionWeights));
            %     pupilDiamsForSingleTrialAttentionWeights = pupilDiamsForSingleTrialAttentionWeights.^.5;
                %% Add four params for setting limits on attention weights
                attentionWeights = pupilDiamsForSingleTrialAttentionWeights(:,1);

                %% Begin model fitting section
                % model has only one parameter, the learning rate, between 0 and 1
                % in general, there will be one of these for each parameter (excluding the
                % softmax parameter); softmax parameter is the inverse temp if I remember
                % correctly
            %     lb = [0,1900,1900,1900,1900]; %lower bounds
            %     ub = [1,3000,3000,3000,3000]; %upper bounds
            %     lb = [0,1000,0];%,1000,0,0]; %lower bounds
            %     ub = [1,3500,1];%,3500,1,1]; %upper bounds
                %% model2
                if modelNum == 2
                    lb = [0,0,0,0,1000];%,0,0];%,0,1000];%,1000,0,0]; %lower bounds
                    ub = [1,1,1,1,5000];%,1,1];%,.01,3000];%,3500,1,1]; %upper bounds
                elseif modelNum == 3
                    %% model3
                    lb = [0,0,0];%,0,0];%,0,1000];%,1000,0,0]; %lower bounds
                    ub = [1,1,1];%,1,1];%,.01,3000];%,3500,1,1]; %upper bounds
                elseif modelNum <10
                    %% model6,7,8,9
                    lb = [0,0,0,0,0];%,0,0];%,0,1000];%,1000,0,0]; %lower bounds
                    ub = [1,1,1,1,1];%,1,1];%,.01,3000];%,3500,1,1]; %upper bounds
                else
                    %% model-11
                    lb = [0,0,0,0,0,0];%,0,0];%,0,1000];%,1000,0,0]; %lower bounds
                    ub = [1,1,1,1,1,1];%,1
                end
                %% Beta limits set in rlfit function
                % beta = modelParams(1); %inverse temperature (softmax choice parameter)
                % alpha = modelParams(2);
                % lambda = modelParams(3); %eligibility trace decay parameter
            %     % p = modelParams{4}
                % w = modelParams(4);


                % % however, we can also use a decorator to add perseveration behavior 
                % % to the model:
                % % the following are limits for a perseveration bonus
                % lb = [-1, lb];
                % ub = [1, ub];
                % Qfun = add_perseveration(@Q_model);
            %     Qfun = (@Q_model6);
                eval(['Qfun = (@Q_model' num2str(modelNum) ');']);
                % now optmize to fit model
                numiter = 50;
                [beta, LL, Q] = rlfit3(Qfun, action1,action2,states,payoff,attentionWeights,lb, ub, numiter);
                action = [action1;action2];
                ppData.params{ppInd} = beta;
                ppData.LL{ppInd} = LL;
                ppData.Q{ppInd} = Q;
                if ppInd == 15
            %         keyboard
                end
                ppData.actions{ppInd}=action;

            %     keyboard
            %     % plot results
            %     figure
            %     plot(Q)
                useChoicePreds = false;
                if useChoicePreds
                    sumError = 0;
                end
                match = 0;
                exploratoryTrialsBoolean_diffToMaxMin = zeros(length(Q),2);
            %     avgQacrossTrials = mean(Q,1);

                for act = 1:length(Q)-1
                    if useChoicePreds
                        qVal = Q(act);
                    else
                        qVal = Q(act,:);
                    end
                    % only compare available options once in final state (decode
                    % options)
                    % Recall: 
                    % action 1 leftImage of leftChoices (blue)
                    % action 2 rightImage of leftChoices (blue)
                    % action 3 leftImage of rightChoices (purple)
                    % action 4 rightImage of leftChoices (purple)
            %         actionData = action(act);%two_stage_task_data.action(blueStateOnlyIndices(act),:);
                    if useChoicePreds
                        actionData=action(act);
                    else
                        actionData = action(act);
                        if actionData == 0
                            actionData = 1;
                        else
                            actionData = 2;
                        end
                    end

            %         predictedAct = round(qVal);
            %         if predictedAct ~= 1 && predictedAct ~=0
            %             keyboard
            %         end
                    if useChoicePreds
                        % use choice error
                        error = abs(qVal-actionData);
                        sumError = sumError + error;
                        if actionData == round(qVal) 
                            match = match + 1;
                            exploratoryTrialsBoolean_diffToMaxMin(act,2) = error;
                        else
                            exploratoryTrialsBoolean_diffToMaxMin(act,1) = 1;
                            exploratoryTrialsBoolean_diffToMaxMin(act,2) = error;
                        end
                    else
                        [maxQ,maxChoice] = max(qVal);
                        if maxChoice == actionData
                            match = match + 1;
                            inds = [1,2];
                            inds = inds(inds~=actionData);
                            diffToMin = maxQ - qVal(inds);

                            exploratoryTrialsBoolean_diffToMaxMin(act,2) = diffToMin;
                        else
                            diffToMax = qVal(actionData)-maxQ;
                            exploratoryTrialsBoolean_diffToMaxMin(act,1) = 1;
                            exploratoryTrialsBoolean_diffToMaxMin(act,2) = diffToMax;
                        end
                    end
                end
                acc = match/length(Q)
                if ~useChoicePreds
                    ppData.acc{ppInd} = acc;
                else
                    sumError
                    ppData.acc{ppInd} = acc;
                end
                ppData.exploratoryTrialsBoolean_diffToMaxMin{ppInd} = exploratoryTrialsBoolean_diffToMaxMin;
            end
        end
        if loadData
            load([pars.dirs.baseDataDir filesep 'ppData_model' num2str(modelNum) '_antic' num2str(pars.grabAnticResp) '_cond' num2str(condition) '_trinaryPupilMeasSetting' num2str(trinaryPupilMeasSetting) '.mat'],'ppData')
%             load([pars.dirs.saveDir filesep 'qDiffRel_ppData_model' num2str(modelNum) '_antic' num2str(pars.grabAnticResp) '_cond' num2str(condition) '_trinaryPupilMeasSetting' num2str(trinaryPupilMeasSetting) '.mat'],'ppData')
        else
            save([pars.dirs.baseDataDir filesep 'qDiffRel_ppData_model' num2str(modelNum) '_antic' num2str(pars.grabAnticResp) '_cond' num2str(condition) '_trinaryPupilMeasSetting' num2str(trinaryPupilMeasSetting) '.mat'],'ppData')
        end
        if ~plotStuff
            continue
        end

        %% Begin plot section
        condition
        trinaryPupilMeasSetting
        modelNum

        diffToMaxes=[];
        qDiffs_explore=[];
        qDiffs_exploit=[];

        pupilSizes_explore=[];
        pupilSizes_exploit=[];

        diffToMins=[];
        binThreshes=[.33,.66,1];

        dataAcrossPps=NaN(45,length(binThreshes),300,2);
        avgAcrossTrials = false;

        %% Compile data across participants and bins
        for ppInd = 1:length(ppData.pupilDiams)
            pupilData = ppData.pupilDiams{ppInd};
            if sum(sum(pupilData == 0)) > 3
                keyboard
            end
            % combine all diffToMax's
            diffToMax=ppData.exploratoryTrialsBoolean_diffToMaxMin{ppInd};
            if isempty(diffToMax)
                continue
            end
            exploreRows = find(diffToMax(:,1)==1);
            missingPupData = find(isnan(pupilData(:,2)) & isnan(pupilData(:,3)));
            pupilData(missingPupData,2) = mean(pupilData(:,2),'omitnan'); % arbitrarily replace nan value of 2rd/3rd col
            ppData.numExploreTrials{ppInd} = length(exploreRows);
            pupilDataAcrossConds = pupilData(:,1);
            pupilDataAcrossConds = [pupilDataAcrossConds;pupilData(~isnan(pupilData(:,2)),2)];
            pupilDataAcrossConds = [pupilDataAcrossConds;pupilData(~isnan(pupilData(:,3)),3)];
            pupilData= pupilDataAcrossConds;
            pupilData_explore = pupilData(exploreRows,1);
            if isempty(qDiffs_explore)
                qDiffs_explore = diffToMax(exploreRows,2);
                pupilSizes_explore = pupilData_explore;
            else
                qDiffs_explore = [qDiffs_explore;diffToMax(exploreRows,2)];
                pupilSizes_explore = [pupilSizes_explore;pupilData_explore];
            end
            diffToMax = diffToMax(exploreRows,:);

            for binThreshInd = 1:length(binThreshes)
                binThresh = binThreshes(binThreshInd);
                if binThreshInd == 1
                    binInds = find(diffToMax(:,2) < binThresh & ~isnan(diffToMax(:,2)));
                elseif binThreshInd == length(binThreshes)
                    binInds = find(diffToMax(:,2) >= binThreshes(binThreshInd-1));
                else
                    binInds = find(diffToMax(:,2) < binThresh & diffToMax(:,2) >= binThreshes(binThreshInd-1));
                end
                if avgAcrossTrials
                    dataAcrossPps(ppInd,binThreshInd,1,1) = mean(pupilData_explore(binInds),'omitnan');
                else
                    dataAcrossPps(ppInd,binThreshInd,1:length(pupilData_explore(binInds)),1) = pupilData_explore(binInds);
                end
            end

            if isempty(diffToMaxes)
                diffToMaxes = diffToMax(:,2);
            else
                diffToMaxes=[diffToMaxes;diffToMax(:,2)];
            end

            %% combine all diffToMax's
            diffToMin=ppData.exploratoryTrialsBoolean_diffToMaxMin{ppInd};
            exploitRows = find(diffToMin(:,1)==0);
            ppData.numExploitTrials{ppInd} = length(exploitRows);

            pupilData_exploit = pupilData(exploitRows,1);

            if isempty(qDiffs_exploit)
                qDiffs_exploit=diffToMin(exploitRows,2);
                pupilSizes_exploit = pupilData_exploit;
            else
                qDiffs_exploit = [qDiffs_exploit;diffToMin(exploitRows,2)];
                pupilSizes_exploit = [pupilSizes_exploit;pupilData_exploit];
            end

            try
                diffToMin = diffToMin(exploitRows,:);
            catch why
                keyboard
            end
            for binThreshInd = 1:length(binThreshes)
                binThresh = binThreshes(binThreshInd);
                if binThreshInd == 1
                    binInds = find(diffToMin(:,2) < binThresh & ~isnan(diffToMin(:,2)));
                elseif binThreshInd == length(binThreshes)
                    binInds = find(diffToMin(:,2) >= binThreshes(binThreshInd-1));
                else
                    binInds = find(diffToMin(:,2) < binThresh & diffToMin(:,2) >= binThreshes(binThreshInd-1));
                end
                if avgAcrossTrials
                    dataAcrossPps(ppInd,binThreshInd,1,2) = mean(pupilData_exploit(binInds),'omitnan');
                else
                    dataAcrossPps(ppInd,binThreshInd,1:length(pupilData_exploit(binInds)),2) = pupilData_exploit(binInds);
                end
            end

            if isempty(diffToMins)
                diffToMins = diffToMin(:,2);
            else
                diffToMins=[diffToMins;diffToMin(:,2)];
            end
            qDiffs_exploit=[];
            qDiffs_explore=[];
            pupilSizes_exploit=[];
            pupilSizes_explore=[];
        end

        %% Plot pupil size according to bins/data distribution
        avgDilationAcrossExploitinessFig=figure;
        globalPlotInd = 1;
        plotInd=1;
        dataForMeanComparisons={};
        dataForMeanComparisons_qDiffs={};
        dataForMeanComparisons_ppNum={};
        dataForMeanComparisons_cumSums={};
        %% Plot means for each bin for data across ALL participants
        for exploitVsExplore = 2:-1:1
            subplot(3,2,plotInd)
            allBins=[];
            for binThresh = 1:size(dataAcrossPps,2)
        %         allPps = squeeze(dataAcrossPps(:,binThresh,:,exploitVsExplore));
                data = squeeze(dataAcrossPps(:,binThresh,:,exploitVsExplore)); % first colon is all pps, second colon is all trials
                dataForMeanComparisons{exploitVsExplore,binThresh,1}=data;
                avgPupilSize_bin = nanmean(data(:));
                if exploitVsExplore == 2
                    allBins=[avgPupilSize_bin,allBins];
                else
                    allBins=[allBins,avgPupilSize_bin];
                end
            end

            bar([1:length(allBins)],smooth(allBins,1))
            if plotInd == 2
                set(gca,'YTickLabel',{'','',''})
                set(gca,'xtick',[]);
                title('explore trials: bestQ-chosenQ','FontSize',14)
            else
                set(gca,'YTickLabel',{'','',''})
                set(gca,'xtick',[]);
                title('exploit trials: chosenQ-worstQ','FontSize',14)
            end

            ylim([1750,2500])
            plotInd = plotInd + 1;
            globalPlotInd=globalPlotInd+1;
        end

        %% Stats for comparing bins (all participants)
        [h,p,ci,stats] = ttest2(dataForMeanComparisons{1,1,1}(:),dataForMeanComparisons{1,2,1}(:));
        [h,p,ci,stats] = ttest2(dataForMeanComparisons{1,1,1}(:),dataForMeanComparisons{1,3,1}(:));

        [h,p,ci,stats] = ttest2(dataForMeanComparisons{2,1,1}(:),dataForMeanComparisons{2,2,1}(:));
        [h,p,ci,stats] = ttest2(dataForMeanComparisons{2,1,1}(:),dataForMeanComparisons{2,3,1}(:));


        %% Plot exploratory param val vs. avg pupil diamter across all trials
        correlations=figure;
        subplot(1,10,1:2)
        for paramInd = 1:1
            paramValsVsDiam = [];
            for ppInd = 1:length(ppData.pupilDiams)
                pupilDiams = ppData.pupilDiams{ppInd};
                if isempty(pupilDiams)
                    paramValsVsDiam(ppInd,1) = nan;
                    paramValsVsDiam(ppInd,2) = nan;
                    continue
                end
                pupilDataAcrossConds = pupilDiams(:,1); % blue state only
                exploreBooleans = ppData.exploratoryTrialsBoolean_diffToMaxMin{ppInd};

        % %         pupilDataAcrossConds = pupilData(:,1);
        %         pupilDataAcrossConds = [pupilDataAcrossConds;pupilDiams(~isnan(pupilDiams(:,2)),2)];
        %         pupilDataAcrossConds = [pupilDataAcrossConds;pupilDiams(~isnan(pupilDiams(:,3)),3)];
                pupilDiams= pupilDataAcrossConds;
        %         pupilDiams = pupilDiams(exploreBooleans(:,1)==1);

                avgDiam = mean(pupilDiams,'omitnan');
                paramVal = ppData.params{ppInd};
                paramVal = paramVal(paramInd);
                paramValsVsDiam(ppInd,1) = paramVal;
                paramValsVsDiam(ppInd,2) = avgDiam;
            end
            figure(correlations)
            scatter(paramValsVsDiam(:,1),paramValsVsDiam(:,2))
            lsline
        %     title(['Beta correlates with avg pre-choice pupil diam.'])
            ylabel('Pupil Diameter','FontSize',12)
            xlabel('Beta Value','FontSize',12)
        %     saveas(gcf,[pars.dirs.figDir filesep 'avgPupilDiamVsBEta_stimOnsetBool' num2str(gatherStimOnset) '_meas' num2str(measOpts) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.png'])
            mdl = fitlm(paramValsVsDiam(:,1),paramValsVsDiam(:,2))
        end
        %% Next specifically look at learning rate
        subplot(1,10,3:4)
        for paramInd = 2:2
        %     figure
            paramValsVsDiam = [];
            for ppInd = 1:length(ppData.pupilDiams)
                pupilDiams = ppData.pupilDiams{ppInd};
                if isempty(pupilDiams)
                    paramValsVsDiam(ppInd,1) = nan;
                    paramValsVsDiam(ppInd,2) = nan;
                    continue
                end
                pupilDataAcrossConds = pupilDiams(:,1); % blue state only
                exploreBooleans = ppData.exploratoryTrialsBoolean_diffToMaxMin{ppInd};
                pupilDiams= pupilDataAcrossConds;

                avgDiam = mean(pupilDiams,'omitnan');
                paramVal = ppData.params{ppInd};
                paramVal = paramVal(paramInd);
                paramValsVsDiam(ppInd,1) = paramVal;
                paramValsVsDiam(ppInd,2) = avgDiam;
            end
            figure(correlations)
            scatter(paramValsVsDiam(:,1),paramValsVsDiam(:,2))
            lsline
            xlabel('Alpha Value','FontSize',12)
            set(gca,'YTickLabel',{'','','','','','',''})
            mdl = fitlm(paramValsVsDiam(:,1),paramValsVsDiam(:,2))
        end
        
        %% Plot #exploreTrials Vs Avg Pupil Diam
        paramValsVsDiam = [];
        for ppInd = 1:length(ppData.pupilDiams)
            pupilDiams = ppData.pupilDiams{ppInd};
            if isempty(pupilDiams)
                paramValsVsDiam(ppInd,1) = nan;
                paramValsVsDiam(ppInd,2) = nan;
                continue
            end
            pupilDiams = pupilDiams(:,1); % blue state only
            pupilDiams = (pupilDiams-min(pupilDiams))/(max(pupilDiams)-min(pupilDiams));
            avgDiam = mean(pupilDiams,'omitnan');
            paramVal = ppData.params{ppInd};
            beta = paramVal(1); % 4 is w (1 is completely model-based)

            numExploitTrials = ppData.numExploitTrials{ppInd};
            %% outlier
            if beta > 5 && numExploitTrials < 160
                paramValsVsDiam(ppInd,1) = nan;
                paramValsVsDiam(ppInd,2) = nan;
            else
                paramValsVsDiam(ppInd,1) = beta;
                paramValsVsDiam(ppInd,2) = numExploitTrials;
            end
        end
        figure(correlations)
        subplot(1,10,7:8)
        scatter(paramValsVsDiam(:,1),paramValsVsDiam(:,2))
        lsline
        % title(['Beta correlates with total exploit actions'])
        ylabel('Total Exploit Actions','FontSize',12)
        xlabel('Beta Value','FontSize',12)
        mdl = fitlm(paramValsVsDiam(:,1),paramValsVsDiam(:,2))

        %% Learning rate correlates w/ #exploit trials
        paramValsVsDiam = [];
        for ppInd = 1:length(ppData.pupilDiams)
            paramVal = ppData.params{ppInd};
            try
                alpha = paramVal(2); % 4 is w (1 is completely model-based)
            catch 
                paramValsVsDiam(ppInd,1) = nan; 
                paramValsVsDiam(ppInd,2) = nan;
                continue
            end
            numExploitTrials = ppData.numExploitTrials{ppInd};
            paramValsVsDiam(ppInd,1) = alpha; 
            paramValsVsDiam(ppInd,2) = numExploitTrials;
        end
        figure(correlations)
        subplot(1,10,9:10)

        scatter(paramValsVsDiam(:,1),paramValsVsDiam(:,2))
        lsline
        mdl = fitlm(paramValsVsDiam(:,1),paramValsVsDiam(:,2))
        xlabel('Alpha Value','FontSize',12)
        set(gca,'YTickLabel',{'','','','','','',''})
%         saveas(gcf,[pars.dirs.figDir filesep 'paramCorrelations_beta_alpha_exploitCt_pupilAvg_trinaryPupilMeasSetting' num2str(trinaryPupilMeasSetting) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.png'])

        %% Compile data for model predictions
        paramInd=4;
        thresh = .33;
        allWvals=[];
        allPupilDiams=[];
        allQdiffs=[];
        allExploreBooleans=[];
        allPps=[];
        allAlphas=[];
        allBetas=[];
        allLambdas=[];
        allCumSums=[];
        allStates=[];
        allOutcomes=[];
        paramValsVsDiam = [];

        for ppInd = 1:length(ppData.pupilDiams)
            try
                load([pars.dirs.baseDataDir filesep 'dataset-eyeTracker-choices' filesep 'two_stage_task_sub' num2str(ppInd) '_bars' num2str(condition) '.mat'])
            catch why
                continue
            end
            outcomes=two_stage_task_data.payoff;
            col1=find(~isnan(outcomes(:,1)));
            col2=find(~isnan(outcomes(:,2)));
            newOutcomes=nan(150,1);
            newOutcomes(col1)=outcomes(col1,1);
            newOutcomes(col2)=outcomes(col2,2);
            outcomes=[newOutcomes;newOutcomes];
            allOutcomes=[allOutcomes,outcomes];
            states = two_stage_task_data.state;
            allStates=[allStates;ones(150,1);states];
            cumSums=zeros(2,150,3);

            actions=ppData.actions{ppInd};
            action1=actions(1:150);
            action2=actions(151:300);

            cumSums(1,:,1)=cumsum(action1==0);
            cumSums(2,:,1)=cumsum(action1==1);

            blueState = NaN(150,1);
            blueState(states==2)=action2(states==2);
            purpleState = NaN(150,1);
            purpleState(states==3)=action2(states==3);

            cumSums(1,:,2)=cumsum(blueState==0);
            cumSums(2,:,2)=cumsum(blueState==1);

            cumSums(1,:,3)=cumsum(purpleState==0);
            cumSums(2,:,3)=cumsum(purpleState==1);
            cumSumsCollapsed=zeros(150,1);
            for stageOneInd = 1:150
                cumSumsCollapsed(stageOneInd)=cumSums(action1(stageOneInd)+1,stageOneInd,1);
            end
            for stageOneInd = 151:300

                cumSumsCollapsed(stageOneInd)=cumSums(action2(stageOneInd-150)+1,stageOneInd-150,states(stageOneInd-150));
            end
%             cumSumsCollapsed=cumSums([action1+1],:,1);
%             ppData.aCumSum{ppInd}=aCumSum;
%             ppData.bCumSum{ppInd}=bCumSum;
            ppData.cumSumsCollapsed{ppInd}=cumSumsCollapsed;            
            
            exploreBooleans = ppData.exploratoryTrialsBoolean_diffToMaxMin{ppInd};
            if isempty(exploreBooleans)
                paramValsVsDiam(ppInd,1) = nan;
                paramValsVsDiam(ppInd,2) = nan;
                paramValsVsDiam(ppInd,3) = nan;
                paramValsVsDiam(ppInd,4) = nan;
                paramValsVsDiam(ppInd,5) = nan;
                paramValsVsDiam(ppInd,6) = nan; 
                paramValsVsDiam(ppInd,7) = nan;
                paramValsVsDiam(ppInd,8) = nan;
                paramValsVsDiam(ppInd,9) = nan;
                paramValsVsDiam(ppInd,10) = nan;
                continue
            end
            
            paramVals = ppData.params{ppInd};
            wVal = paramVals(paramInd);

            qDiffs=(exploreBooleans(:,2));
            exploreBooleans=(exploreBooleans(:,1));
            allExploreBooleans=[allExploreBooleans;exploreBooleans];
            pupilDiams = ppData.pupilDiams{ppInd};
            missingPupData = find(isnan(pupilDiams(:,2)) & isnan(pupilDiams(:,3)));
            pupilDiams(missingPupData,2) = mean(pupilDiams(:,2),'omitnan');
            pupilDataAcrossConds = pupilDiams(:,1);
            pupilDataAcrossConds = [pupilDataAcrossConds;pupilDiams(~isnan(pupilDiams(:,2)),2)];
            pupilDataAcrossConds = [pupilDataAcrossConds;pupilDiams(~isnan(pupilDiams(:,3)),3)];
            pupilDiams= pupilDataAcrossConds;
            exploitPupilDiams = pupilDiams(find(exploreBooleans==0));
            
            %% find instances where exploit came after explore
            exploitAfterExploreIndices = find(diff(exploreBooleans)==-1)+1;
            exploitAfterExploit = find([1;diff(exploreBooleans)==0] & exploreBooleans==0);
            exploitAfterExplorePupilDiams=pupilDataAcrossConds(exploitAfterExploreIndices);
            exploitAfterExploitPupilDiams=pupilDataAcrossConds(exploitAfterExploit);

%             figure
%             scatter(qDiffs,pupilDiams)
            exploitOrExplore=0;
            allQdiffs=[allQdiffs;qDiffs];
%             allQdiffs=[allQdiffs;qDiffs(find(exploreBooleans==exploitOrExplore))];
%             pupilDiams = pupilDiams(find(exploreBooleans==exploitOrExplore));
%             cumSumsCollapsed=cumSumsCollapsed(find(exploreBooleans==exploitOrExplore));
            allPupilDiams=[allPupilDiams;pupilDiams];
            if length(allPupilDiams) ~= length(allQdiffs)
                keyboard
            end
            
            allCumSums=[allCumSums;cumSumsCollapsed];
            wVals = ones(length(pupilDiams),1)*wVal;
            allWvals=[allWvals;wVals];
            
            alphaVal = paramVals(2);
            alphas = ones(length(pupilDiams),1)*alphaVal;
            allAlphas=[allAlphas;alphas];

            beta = paramVals(1);
            betas = ones(length(pupilDiams),1)*beta;
            allBetas=[allBetas;betas];

            lambda = paramVals(3);
            lambdas = ones(length(pupilDiams),1)*lambda;
            allLambdas=[allLambdas;lambdas];

            pps = ones(length(pupilDiams),1)*ppInd;
            allPps=[allPps;pps];
            if length(allPupilDiams) ~= length(allPps)
                keyboard
            end
            avgDiam = mean(pupilDiams,'omitnan');
            paramValsVsDiam(ppInd,1) = wVal;
            paramValsVsDiam(ppInd,2) = avgDiam;
            meanQdiff = mean(qDiffs);
            paramValsVsDiam(ppInd,3) = meanQdiff;
            paramValsVsDiam(ppInd,4) = nanmean(exploreBooleans);
            paramValsVsDiam(ppInd,5) = nanmean(cumSumsCollapsed);
            paramValsVsDiam(ppInd,6) = sum(exploreBooleans); % total exploratory trials
            paramValsVsDiam(ppInd,7) = nanmean(exploitPupilDiams);
            paramValsVsDiam(ppInd,8) = beta;
            paramValsVsDiam(ppInd,9) = lambda;
            paramValsVsDiam(ppInd,10) = alphaVal;
            paramValsVsDiam(ppInd,11) = nanmean(exploitAfterExplorePupilDiams);
            paramValsVsDiam(ppInd,12) = nanmean(exploitAfterExploitPupilDiams);
            
        end
        
        %% Create model and plot cumulative model predictions
        allWvals=allWvals(:);
        %
%         thresh=.33;
%         allOutcomes=allOutcomes(allWvals>thresh);
%         allStates=allStates(allWvals>thresh);
%         allExploreBooleans=allExploreBooleans(allWvals>thresh);
%         allCumSums=allCumSums(allWvals>thresh);
%         allAlphas=allAlphas(allWvals>thresh);
%         allBetas=allBetas(allWvals>thresh);
%         allLambdas=allLambdas(allWvals>thresh);
%         allPupilDiams=allPupilDiams(allWvals>thresh);
%         allQdiffs=allQdiffs(allWvals>thresh);
%         allPps=allPps(allWvals>thresh);
%         allWvals=allWvals(allWvals>thresh);
        % model-free
%         allOutcomes=allOutcomes(allWvals<thresh);
%         allStates=allStates(allWvals<thresh);
%         allExploreBooleans=allExploreBooleans(allWvals<thresh);
%         allCumSums=allCumSums(allWvals<thresh);
%         allAlphas=allAlphas(allWvals<thresh);
%         allBetas=allBetas(allWvals<thresh);
%         allLambdas=allLambdas(allWvals<thresh);
%         allPupilDiams=allPupilDiams(allWvals<thresh);
%         allQdiffs=allQdiffs(allWvals<thresh);
%         allPps=allPps(allWvals<thresh);
%         allWvals=allWvals(allWvals<thresh);

        tbl3 = table(allOutcomes(:),allStates(:),allExploreBooleans(:),allCumSums(:),allAlphas(:),allBetas(:),allLambdas(:),allPupilDiams(:),allQdiffs(:),allWvals(:),allPps(:),'VariableNames',{...
            'allOutcomes','allStates','allExploreBooleans','allCumSums','allAlphas','allBetas','allLambdas','allPupilDiams','allQdiffs','allWvals','IDs'});
        tbl3.allExploreBooleans = nominal(tbl3.allExploreBooleans);
        tbl3.allStates = nominal(tbl3.allStates);
        lm2 = fitlm(tbl3,['allPupilDiams~((allOutcomes*allStates+allOutcomes*allWvals+allQdiffs+allCumSums+allAlphas*allWvals+allBetas*allWvals+allLambdas*allWvals)*allExploreBooleans)'])% ...
        
%         lm2 = fitlm(tbl3,['allPupilDiams~((allOutcomes*allStates+allOutcomes*allWvals+allCumSums+allAlphas*allWvals+allBetas*allWvals+allLambdas*allWvals)*allExploreBooleans)'])% ...

        % plot main effects
        figure
        plotEffects(lm2)
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Times','fontsize',14)
%         saveas(gcf,[pars.dirs.figDir filesep 'mainEffects' '_model_' num2str(modelNum) '_trinaryPupilMeasSetting' num2str(trinaryPupilMeasSetting) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.png'])
        
        figure
        plotAdded(lm2,'allExploreBooleans')
        
        figure
        plot(lm2)
        ylabel('Pre-choice pupil diameter')
        title('')
%         saveas(gcf,[pars.dirs.figDir filesep 'linearModel' '_model_' num2str(modelNum) '_trinaryPupilMeasSetting' num2str(trinaryPupilMeasSetting) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.png'])
        lmx = fitlm(tbl3,['allWvals~allBetas'])% ...
        lmx = fitlm(tbl3,['allPupilDiams~allBetas'])% ...
        lmx = fitlm(tbl3,['allPupilDiams~allAlphas'])% ...
        lmx = fitlm(tbl3,['allPupilDiams~allWvals'])% ...
        lmx = fitlm(tbl3,['allAlphas~allBetas'])% ...
        lmx = fitlm(tbl3,['allAlphas~allWvals'])% ...

        lmx = fitlm(tbl3,['allBetas~allExploreBooleans'])% ...

        figure
        plot(lmx)
        

        %% INTERACTIONS
        numCols = 3;
        numRows = 1;
%         linearInd = sub2ind(matrixSize, rowSub, colSub)
        spacing =.01;
        yMin=1500;
        yMax=2400;
        %% wval interacts with outcomes
        figure
        plotInteraction(lm2,'allWvals','allOutcomes','predictions')
        ylabel('Pre-choice pupil diameter')
        title('')
        ylim([yMin,yMax])
%         saveas(gcf,[pars.dirs.figDir filesep 'w-outcome-interaction' '_model_' num2str(modelNum) '_trinaryPupilMeasSetting' num2str(trinaryPupilMeasSetting) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.png'])

        figure
        %% Increasing alphas has steeper effect on exploitation trials than exploration trials in terms of increasing pupil diameter.  
        rowSub=1;
        colSub=1;
        linearInd = (numCols.*rowSub)-(numCols-colSub);
        subtightplot(numRows,numCols,linearInd,spacing)
        plotInteraction(lm2,'allExploreBooleans','allAlphas','predictions')
        ylabel('Pre-choice pupil diameter')
        title('')
        ylim([yMin,yMax])
        legend(gca,'off');
        %% Increasing beta has steeper affect on exploration trials than exploitation trials in terms of decreasing pupil diameter.        
        rowSub=1;
        colSub=2;
        linearInd = (numCols.*rowSub)-(numCols-colSub);
        subtightplot(numRows,numCols,linearInd,spacing)
        plotInteraction(lm2,'allExploreBooleans','allBetas','predictions')
        ylabel('')
        title('')
        ylim([yMin,yMax])
        set(gca,'yticklabel',{'','','','','','',''})
        legend(gca,'off');
        %% Increasing wVal has steeper affect on exploration trials than exploitation trials in terms of decreasing pupil diameter.        
        rowSub=1;
        colSub=3;
        linearInd = (numCols.*rowSub)-(numCols-colSub);
        subtightplot(numRows,numCols,linearInd,spacing)
        plotInteraction(lm2,'allExploreBooleans','allWvals','predictions')
        ylabel('')
        title('')
        ylim([yMin,yMax])
        set(gca,'yticklabel',{'','','','','','',''})
%         legend(gca,'off');
%         saveas(gcf,[pars.dirs.figDir filesep 'linearModel_exploreExploitInteraction' '_model_' num2str(modelNum) '_trinaryPupilMeasSetting' num2str(trinaryPupilMeasSetting) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.png'])

        figure        
        %% wVal vs alpha
        rowSub=1;
        colSub=1;
        linearInd = (numCols.*rowSub)-(numCols-colSub);
        subtightplot(numRows,numCols,linearInd,spacing) 
        plotInteraction(lm2,'allWvals','allAlphas','predictions')
        ylabel('')
        title('')
        ylim([yMin,yMax])
%         set(gca,'yticklabel',{'','','','','','',''})
        ylabel('Pre-choice pupil diameter')
        legend(gca,'off');
        %% wVal vs beta
        rowSub=1;
        colSub=2;
        linearInd = (numCols.*rowSub)-(numCols-colSub);
        subtightplot(numRows,numCols,linearInd,spacing) 
        plotInteraction(lm2,'allWvals','allBetas','predictions')
        ylabel('')
        title('')
        ylim([yMin,yMax])
        set(gca,'yticklabel',{'','','','','','',''})
        legend(gca,'off');
        %% wVal vs lambdas
        rowSub=1;
        colSub=3;
        linearInd = (numCols.*rowSub)-(numCols-colSub);
        subtightplot(numRows,numCols,linearInd,spacing)
        x=plotInteraction(lm2,'allWvals','allLambdas','predictions')
        ylabel('')
        title('')
        ylim([yMin,yMax])
        set(gca,'yticklabel',{'','','','','','',''})
%         saveas(gcf,[pars.dirs.figDir filesep 'linearModel_wInteraction' '_model_' num2str(modelNum) '_trinaryPupilMeasSetting' num2str(trinaryPupilMeasSetting) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.png'])


        %% Plot param val vs. avgPupDrops
        % for paramInd = 1:length(ppData.params{1})
        %     figure
        %     paramValsVsDiam = [];
        %     for ppInd = 1:length(ppData.pupilDiams)
        %         pupDropCts = ppData.pupDrops{ppInd};
        %         pupDropCts = mean(mean(pupDropCts));
        % %         avgDiam = mean(pupDropCts,'omitnan');
        %         paramVal = ppData.params{ppInd};
        %         paramVal = paramVal(paramInd);
        %         paramValsVsDiam(ppInd,1) = paramVal;
        %         paramValsVsDiam(ppInd,2) = pupDropCts;
        %     end
        %     scatter(paramValsVsDiam(:,1),paramValsVsDiam(:,2))
        %     title(['avgPupDrogs vs. param val' num2str(paramInd)])
        %     mdl = fitlm(paramValsVsDiam(:,1),paramValsVsDiam(:,2))
        % end

        %% Plot distribution of Q-diffs
        figure
        hist(diffToMaxes,100)
        title('Diff to Max Distribution')
        figure
        hist(diffToMins,100)
        title('Diff to Min Distribution')

        %% Calc BIC values
        %% Calc Sum BIC values
        sumBIC = 0;
        for ppInd = 1:length(ppData.acc)
            if isempty(ppData.LL{ppInd})
                continue
            end
            BIC = (2*log(150)*length(ppData.params{ppInd})) - (2*ppData.LL{ppInd});
            sumBIC = sumBIC + BIC;
        end
        disp(['sumBIC = ' num2str(sumBIC)])
        % keyboard

        for pp = 1:length(ppData.pupilDiams)
            params = ppData.params{pp};
            try
                w = params(4);
            catch
                continue
            end
            allWs(pp) = w;
        end
        figure
        hist(allWs)
        % title('Distribution of w Weights (Model-free Vs. Model-based Learning)')
        % w=1 is pure model based
        modelBasedPps = find(allWs>.66);
        modelFreePps = find(allWs<.33);
        % saveas(gcf,[pars.dirs.figDir filesep 'MF-MB-pp-distributione_stimOnsetBool' num2str(gatherStimOnset) '_meas' num2str(measOpts) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.png'])

        %% Plot exploratory param val vs. avg pupil diamter across all trials

        % bin1Thresh=.33;
        % bin3Thresh=.33;
        % bin4Thresh=.66;
        % bin1Thresh=.05;
        % bin2Thresh=.1;
        % bin3Thresh=.2;
        % bin4Thresh=.4;
        % bin5Thresh=bin4Thresh;
        % dataAcrossPps=NaN(45,length(binThreshes),150,2);
        for modelBasedModelFreeTests = 1:2
            if modelBasedModelFreeTests ==1
                pps = modelBasedPps;
                titleAppend = 'model-based';
            else
                pps = modelFreePps;
                titleAppend = 'model-free';
            end
            dataAcrossPps=NaN(45,length(binThreshes),300,2);
            cumSumDataAcrossPps=NaN(45,length(binThreshes),300,2);
            qDiffs=NaN(45,length(binThreshes),300,2);
            ppMatrix = NaN(45,300);
            ppCt2=1;
            for ppNum = 1:45
                if sum(find(pps==ppNum)) > 0
                    ppMatrix(ppCt2,:) = ppNum;
                    ppCt2=ppCt2+1;
                end
            end
            ppCt = 1;
            diffToMaxes=[];
            diffToMins=[];
            paramValsVsDiam = [];
            for ppInd = pps
                pupilDiams = ppData.pupilDiams{ppInd};
                if isempty(pupilDiams)
                    continue
                end
                pupilDiams = pupilDiams(:,1);
                avgDiam = mean(pupilDiams,'omitnan');
                for paramInd = 1:length(ppData.params{1})
                    paramVal = ppData.params{ppInd};
                    paramVal = paramVal(paramInd);
                    paramValsVsDiam(ppCt,paramInd,1) = paramVal;
                    paramValsVsDiam(ppCt,paramInd,2) = avgDiam;
                end
                %% Measure exploriness vs puipl size
                pupilData = ppData.pupilDiams{ppInd};

                missingPupData = find(isnan(pupilData(:,2)) & isnan(pupilData(:,3)));
                pupilData(missingPupData,2) = mean(pupilData(:,2),'omitnan');

                pupilDataAcrossConds = pupilData(:,1);
                pupilDataAcrossConds = [pupilDataAcrossConds;pupilData(~isnan(pupilData(:,2)),2)];
                pupilDataAcrossConds = [pupilDataAcrossConds;pupilData(~isnan(pupilData(:,3)),3)];
                pupilData= pupilDataAcrossConds;
                % combine all diffToMax's
                diffToMax=ppData.exploratoryTrialsBoolean_diffToMaxMin{ppInd};
                exploreRows = find(diffToMax(:,1)==1);
                pupilData_explore = pupilData(exploreRows,1);
                cumSumsCollapsed_explore = cumSumsCollapsed(exploreRows);
                if sum(size(pupilData_explore) == size(cumSumsCollapsed_explore)) ~= 2
                    keyboard
                end

                diffToMax = diffToMax(exploreRows,:);

                for binThreshInd = 1:length(binThreshes)
                    binThresh = binThreshes(binThreshInd);
                    if binThreshInd == 1
                        binInds = find(diffToMax(:,2) < binThresh & ~isnan(diffToMax(:,2)));
                    elseif binThreshInd == length(pupilSizes_explore)
                        binInds = find(diffToMax(:,2) >= binThreshes(binThreshInd-1));
                    else
                        binInds = find(diffToMax(:,2) < binThresh & diffToMax(:,2) >= binThreshes(binThreshInd-1));
                    end
                    if avgAcrossTrials
                        dataAcrossPps(ppCt,binThreshInd,1,1) = mean(pupilData_explore(binInds),'omitnan');
                    else
                        dataAcrossPps(ppCt,binThreshInd,1:length(pupilData_explore(binInds)),1) = pupilData_explore(binInds);
                        cumSumDataAcrossPps(ppCt,binThreshInd,1:length(cumSumsCollapsed_explore(binInds)),1) = cumSumsCollapsed_explore(binInds);
                        qDiffs(ppCt,binThreshInd,1:length(cumSumsCollapsed_explore(binInds)),1) = diffToMax(binInds,2);
                    end
                end

                if isempty(diffToMaxes)
                    diffToMaxes = diffToMax(:,2);
                else
                    diffToMaxes=[diffToMaxes;diffToMax(:,2)];
                end
                diffToMin=ppData.exploratoryTrialsBoolean_diffToMaxMin{ppInd};
                exploitRows = find(diffToMin(:,1)==0);
                cumSumsCollapsed_exploit = cumSumsCollapsed(exploitRows);
                if isempty(exploitRows)
                    keyboard
                end

                try
                    pupilData_exploit = pupilData(exploitRows,1);
                    diffToMin = diffToMin(exploitRows,:);
                catch why
                    keyboard
                end
                for binThreshInd = 1:length(binThreshes)
                    binThresh = binThreshes(binThreshInd);
                    if binThreshInd == 1
                        binInds = find(diffToMin(:,2) < binThresh & ~isnan(diffToMin(:,2)));
                    elseif binThreshInd == length(pupilSizes_explore)
                        binInds = find(diffToMin(:,2) >= binThreshes(binThreshInd-1));
                    else
                        binInds = find(diffToMin(:,2) < binThresh & diffToMin(:,2) >= binThreshes(binThreshInd-1));
                    end

                    if avgAcrossTrials
                        dataAcrossPps(ppCt,binThreshInd,1,2) = mean(pupilData_exploit(binInds),'omitnan');
                    else
                        dataAcrossPps(ppCt,binThreshInd,1:length(pupilData_exploit(binInds)),2) = pupilData_exploit(binInds);
                        cumSumDataAcrossPps(ppCt,binThreshInd,1:length(cumSumsCollapsed_exploit(binInds)),2) = cumSumsCollapsed_exploit(binInds);
                        qDiffs(ppCt,binThreshInd,1:length(cumSumsCollapsed_exploit(binInds)),2) = diffToMin(binInds,2);
                    end
                end

                ppCt = ppCt + 1;

                if isempty(diffToMins)
                    diffToMins = diffToMin(:,2);
                else
                    diffToMins=[diffToMins;diffToMin(:,2)];
                end
            end
            %% Plot parameter correlations with pppDiam
            for paramInd = 1:length(ppData.params{1})
        %         pupilDiams = ppData.pupilDiams{ppInd};
        %         pupilDiams = pupilDiams(:,1); % blue state only
        %         avgDiam = mean(pupilDiams,'omitnan');
        %         paramVal = ppData.params{ppInd};
        %         paramVal = paramVal(paramInd);
        %         paramValsVsDiam(:,paramInd,1);
        %         paramValsVsDiam(:,paramInd,2);
                figure
                scatter(paramValsVsDiam(:,paramInd,1),paramValsVsDiam(:,paramInd,2))
                title(['param' num2str(paramInd) ', ' titleAppend])
                disp([titleAppend ', param' num2str(paramInd)])
                mdl = fitlm(paramValsVsDiam(:,paramInd,1),paramValsVsDiam(:,paramInd,2))
            end



            %% Plot pupil size according to bins/data distribution
        %     exploitVsExploreFig_binned=figure;
        %     exploitVsExploreFig_avg=figure;

            ind=1;
            for exploitVsExplore = 2:-1:1
        %         figure(exploitVsExploreFig_binned)
                figure(avgDilationAcrossExploitinessFig)
                subplot(3,2,globalPlotInd)


                allBins=[];
                for binThresh = 1:size(dataAcrossPps,2)

                    qDiffData = squeeze(qDiffs(:,binThresh,:,exploitVsExplore));
                    data = squeeze(dataAcrossPps(:,binThresh,:,exploitVsExplore));
                    dataCumSum = squeeze(cumSumDataAcrossPps(:,binThresh,:,exploitVsExplore));

                    dataForMeanComparisons_qDiffs{exploitVsExplore,binThresh,1+modelBasedModelFreeTests} = qDiffData;
                    dataForMeanComparisons_ppNum{exploitVsExplore,binThresh,1+modelBasedModelFreeTests} = ppMatrix;
                    dataForMeanComparisons{exploitVsExplore,binThresh,1+modelBasedModelFreeTests} = data;
                    dataForMeanComparisons_cumSums{exploitVsExplore,binThresh,1+modelBasedModelFreeTests} = dataCumSum;
                    avgPupilSize_bin = nanmean(data(:));
                    if exploitVsExplore == 2
                        allBins=[avgPupilSize_bin,allBins];
                    else
                        allBins=[allBins,avgPupilSize_bin];
                    end
                end

                bar([1:length(allBins)],smooth(allBins,1))
                if ind == 2
                    set(gca,'ytick',[])
                    if globalPlotInd==6
                        set(gca,'XTickLabel',{'x<.33','.33<=x<=.66','x>.66'},'FontSize',16)
                    else
                        set(gca,'xtick',[])
                    end
        %             title('explore trials: bestQ-chosenQ','FontSize',12)
                else
                    if globalPlotInd==5
                        ylabel('Pupil Diameter','FontSize',30)
        %                 ylabs=get(gca,'YTickLabel')
        %                 set(gca,'YTickLabel',ylabs,'FontSize',10)
                        set(gca,'XTickLabel',{'x>.66','.33<=x<=.66','x<.33'},'FontSize',16)

                    else
                        set(gca,'xtick',[])
        %                 set(gca,'ytick',[])
                        set(gca,'YTickLabel',{'','',''})

                    end
        %             title('exploit trials: chosenQ-worstQ','FontSize',12)
                end

        %         ylim([0,1])
                if condition == 0
            %         ylim([1750,2000])
                end

        %         bar([1:length(allBins)],allBins)
        %         x=smooth(allBins(~isnan(allBins)),20);
        %         plot(x)
        %         ylim([0.1,.9])
        %         if condition == 0
                ylim([1750,2500])
        %         end
        % %         figure
                ind=ind+1;
                globalPlotInd=globalPlotInd+1;
            end
        end
        %% Plot pupilSize ~ exploitBoolean * qDiff * familiarity
        for modelBasedOrFree = 2:3
            for exploitVsExplore = 1:2
                qValPvals{modelBasedOrFree,exploitVsExplore}=[];
                cumSumPvals{modelBasedOrFree,exploitVsExplore}=[];
                qEsts{modelBasedOrFree,exploitVsExplore}=[];
                cumSumEsts{modelBasedOrFree,exploitVsExplore}=[];
            end
        end
        for pp = 1:45
    %         statsTests=figure;
    %         colorCoord=figure;
            figure
            for modelBasedOrFree = 2:3 % 2=model-based (row1); 3 is model-free
                if modelBasedOrFree == 2
                    figTitle = 'MB';
                else
                    figTitle='MF';
                end
                for exploitVsExplore = 1:2 % 1=explore
                    pupilData = dataForMeanComparisons{exploitVsExplore,1,modelBasedOrFree}(:);
                    pupilData=[pupilData;dataForMeanComparisons{exploitVsExplore,2,modelBasedOrFree}(:)];
                    pupilData=[pupilData;dataForMeanComparisons{exploitVsExplore,3,modelBasedOrFree}(:)];

                    qVals = dataForMeanComparisons_qDiffs{exploitVsExplore,1,modelBasedOrFree}(:);
                    qVals = [qVals;dataForMeanComparisons_qDiffs{exploitVsExplore,2,modelBasedOrFree}(:)];
                    qVals = [qVals;dataForMeanComparisons_qDiffs{exploitVsExplore,3,modelBasedOrFree}(:)];

                    pps = dataForMeanComparisons_ppNum{exploitVsExplore,1,modelBasedOrFree}(:);
                    pps = [pps;dataForMeanComparisons_ppNum{exploitVsExplore,2,modelBasedOrFree}(:)];
                    pps = [pps;dataForMeanComparisons_ppNum{exploitVsExplore,3,modelBasedOrFree}(:)];
                    pps = pps(~isnan(pupilData));

                    length(unique(pps(~isnan(pps))));

                    cumSums = dataForMeanComparisons_cumSums{exploitVsExplore,1,modelBasedOrFree}(:);
                    cumSums = [cumSums;dataForMeanComparisons_cumSums{exploitVsExplore,2,modelBasedOrFree}(:)];
                    cumSums = [cumSums;dataForMeanComparisons_cumSums{exploitVsExplore,3,modelBasedOrFree}(:)];
                    cumSums = cumSums(~isnan(pupilData));

                    qVals = qVals(~isnan(pupilData));
                    pupilData = pupilData(~isnan(pupilData));

                    uniquePps = unique(pps);
                    if isempty(pupilData(pps==pp))
                        continue
                    end
    %                 tbl = table(pupilData,qVals,cumSums,pps,'VariableNames',{'pupilData','qVals','cumSums','pps'});

                    tbl = table(pupilData(pps==pp),qVals(pps==pp),cumSums(pps==pp),pps(pps==pp),'VariableNames',{'pupilData','qVals','cumSums','pps'});
                    % Fit a linear regression model for miles per gallon (MPG) with weight and acceleration as the predictor variables.
        %             lm = fitlm(tbl,'pupilData~qVals+pps')
                    tbl.pps = nominal(tbl.pps);
    %                 try
                    lm = fitlm(tbl,'pupilData~qVals+cumSums')%+pps+cumSums')

    %                     lm = fitlm(tbl,'pupilData~qVals')%+pps+cumSums')
    %                     qValPval=lm.anova.pValue(1);
    %                     cumSumPval=lm.anova.pValue(2);
    %                     qEst=lm.Coefficients.Estimate(2);
    %                     cumSumEst=lm.Coefficients.Estimate(3);
    %                     qValPvals=[qValPvals,qValPval];
    %                     cumSumPvals=[cumSumPvals,cumSumPval];
    %                     
    %                     qEsts_temp=qEsts{modelBasedOrFree,exploitVsExplore};
    %                     qEsts_temp=[qEsts_temp,qEst];
    %                     qEsts{modelBasedOrFree,exploitVsExplore}=qEsts_temp;
    %                     cumSumEsts{modelBasedOrFree,exploitVsExplore}=cumSumEsts_temp;
    %                     
    %                     cumSumEsts_temp=cumSumEsts{modelBasedOrFree,exploitVsExplore};
    %                     cumSumEsts_temp=[cumSumEsts_temp,cumSumEst];
    %                     cumSumEsts{modelBasedOrFree,exploitVsExplore}=cumSumEsts_temp;
    %                     
    %                     
    %                     figTitle
    %                 catch why
    %                     keyboard
    %                 end
                    plotInd = (2*(modelBasedOrFree-1))-(2-exploitVsExplore);
        %             plotInd = (totalCols*plotRow)-(totalCols-plotCol);
        %                 x(plotInd)=mean(pupilData)
        %                 figure(statsTests)
                    subplot(2,2,plotInd)
                    plot(lm)
                    title(figTitle)
                    ylim([0,3500])


        %             figure(colorCoord)
        %             subplot(2,2,plotInd)
        %             gscatter(qVals,pupilData,pps)
        %             ylim([0,3500])
    %                 figTitle
    %                 pause(2)
                end
            end
        end
    %     keyboard


    % figure()
    % gscatter(pupilData,qVals,pps,'bgr','x.o')


    %     lm = fitlm(tbl,'pupilData~qVals')
    %     figure
    %     plot(lm)
        % titleAppend='model-free';
        % saveas(gcf,[pars.dirs.figDir filesep titleAppend '_stimOnsetBool' num2str(gatherStimOnset) '_meas' num2str(measOpts) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.png'])
        % titleAppend='model-based';
        % saveas(gcf,[pars.dirs.figDir filesep titleAppend '_stimOnsetBool' num2str(gatherStimOnset) '_meas' num2str(measOpts) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.png'])
        saveas(gcf,[pars.dirs.figDir filesep 'MF_BM_Qdifs_v_pupilSizes_stimOnsetBool' num2str(gatherStimOnset) '_meas' num2str(measOpts) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.png'])
        % saveas(gcf,[pars.dirs.figDir filesep 'MF_BM_Qdifs_v_pupilSizes_stimOnsetBool' num2str(gatherStimOnset) '_meas' num2str(measOpts) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.tif'])

        %% Average and compare pupil sizes between learning styles
        useMeanAcrossTrials = false;
        for ppType = 1:2
            if ppType == 1
                pps = modelFreePps;
            else
                pps = modelBasedPps;
            end
            allPupDiams = [];
            for modelPp = pps
                if useMeanAcrossTrials
                    if isempty(allPupDiams)
                        allPupDiams = mean(ppData.pupilDiams{modelPp},'omitnan');
                    else
                        allPupDiams = [allPupDiams;mean(ppData.pupilDiams{modelPp},'omitnan')];
                    end
                else
                    if isempty(allPupDiams)
                        allPupDiams = (ppData.pupilDiams{modelPp});
                    else
                        allPupDiams = [allPupDiams;(ppData.pupilDiams{modelPp})];
                    end
                end
            end
        %     avgs = [];
        %     for state = 2:2
        %         avgDiam = mean(allPupDiams(:,state),'omitnan');
        %         avgs=[avgs,avgDiam];
        %     end
            if ppType == 1
                avgs_mf = (allPupDiams);
            else
                avgs_mb = (allPupDiams);
            end
        end
    %     dataForMeanComparisons
        %% Test means of bins against each other
        cols=[2,1];
        for group = 2:3
            for exploitExplore = 2:-1:1
                for bin1 = [1,2];
                    for bin2=[2,3];
                        if bin1 == 2 && bin2 == 2
                            continue
                        end
                        col=cols(exploitExplore);
                        row=group;
                        disp(['row' num2str(row) ', col' num2str(col)])
                        disp(['bin1:' num2str(bin1) ', bin2:' num2str(bin2)])

                        [h,p,ci,stats] = ttest2(dataForMeanComparisons{exploitExplore,bin1,group}(:),dataForMeanComparisons{exploitExplore,bin2,group}(:));
                        p
                        stats.df
                        stats.tstat
                    end
                end
            end
        end
        datas={};
        avgs=[];
        groupInd=1;
        for group = 2:3
            ind=1;
            for exploitExplore = 2:-1:1
                data=[];
                for bin = 1:3
                    if isempty(data)
                        data=dataForMeanComparisons{exploitExplore,bin,group};
                    else
                        data=[data;dataForMeanComparisons{exploitExplore,bin,group}];
                    end

                end
                datas{ind,groupInd}=data;
                avgs(groupInd,ind) = nanmean(data(:));
                ind=ind+1;
            end
            groupInd=groupInd+1;
        end
        % avgs_mf=mean(avgs_mf,2);
        % avgs_mb=mean(avgs_mb,2);
        figure
        subplot(2,1,1)
        bar([1,2],[avgs(2,1),avgs(1,1)])
    %     title('Model-Based Vs. Model-Free Pupil Diameter For Exploitative Trials')
        ylim([1800,2100])
        set(gca,'XTickLabel',{'Model-Free','Model-Based'})

        subplot(2,1,2)
        bar([1,2],[avgs(2,2),avgs(1,2)])
    %     [h,p,ci,stats] = ttest2(datas{1,1}(:),datas{1,2}(:))
    %     [h,p,ci,stats] = ttest2(datas{2,1}(:),datas{2,2}(:))


        ylim([1800,2100])
        set(gca,'XTickLabel',{'Model-Free','Model-Based'})
        saveas(gcf,[pars.dirs.figDir filesep 'stimOnsetBool' num2str(gatherStimOnset) '_meas' num2str(measOpts) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize_splitByExploreOrExploit.png'])

        figure
        bar([1,2],[mean(avgs_mf(:,1),'omitnan'),mean(avgs_mb(:,1),'omitnan')])

        meanMF=mean(avgs_mf(:,1),'omitnan')
        stdMF=std(avgs_mf(:,1),'omitnan')
        meanMB=mean(avgs_mb(:,1),'omitnan')
        stdMB=std(avgs_mb(:,1),'omitnan')

        ylim([1500,2200])
        % title('model-free vs model-based pupil dilation')
        % xticklabels({'Model-Free','Model-Based'})
        set(gca,'XTickLabel',{'Model-Free','Model-Based'})
        saveas(gcf,[pars.dirs.figDir filesep 'stimOnsetBool' num2str(gatherStimOnset) '_meas' num2str(measOpts) '_cond' num2str(condition) '_modelFreeVsModelBasedAvgPupilSize.png'])

        figure
        subtightplot(2,1,1)
        hist(avgs_mf(:,1))
        xlim([0,4000])
        subtightplot(2,1,2)
        hist(avgs_mb(:,1))
        xlim([0,4000])



        if condition == 0
        %     ylim([1500,2500])
        end
        % sharedRows = avgs_mf(:,~isnan(avgs_mf(:,1)) && ~isnan(avgs_mb(:,2)))

        [h,p,ci,stats] = ttest2(avgs_mf(:,1),avgs_mb(:,1))
        % [h,p,ci,stats] = ttest2(avgs_mf(:,2),avgs_mb(:,2))
        % [h,p,ci,stats] = ttest2(avgs_mf(:,3),avgs_mb(:,3))


        figure
        pupilAvgs=[];
        accs=[];
        for ppInd = 1:length(ppData.pupilDiams)
            qDiffs = ppData.exploratoryTrialsBoolean_diffToMaxMin{ppInd};
            trials=find(qDiffs(1:150,1) == 0);
            subtightplot(7,7,ppInd);
            acc=ppData.acc{ppInd};
            pupilDiams = ppData.pupilDiams{ppInd};
            avgPup=mean(pupilDiams(:,1),'omitnan');
            pupilAvgs=[pupilAvgs,avgPup];
            accs=[accs,acc];
            hold on
            scatter(pupilDiams(trials,1),qDiffs(trials,2).^.2)
        end
        %% Grab mean and std paramVals
        allParamVals=[];
        for ppInd = 1:length(ppData.pupilDiams)
            for paramInd = 1:length(ppData.params{1})
                paramVal = ppData.params{ppInd};
                paramVal = paramVal(paramInd);
                paramValsVsDiam(ppInd,paramInd,1) = paramVal;
                allParamVals(ppInd,paramInd) = paramVal;
                paramValsVsDiam(ppCt,paramInd,2) = avgDiam;
            end
        end
        (nanmean(allParamVals,1))
        (nanstd(allParamVals,1))

        for ppType = 1:2
            if ppType == 1
                pps = modelFreePps;
            else
                pps = modelBasedPps;
            end
            allParamVals=NaN(45,4);

            for ppInd = pps
                for paramInd = 1:length(ppData.params{1})
                    paramVal = ppData.params{ppInd};
                    paramVal = paramVal(paramInd);
                    allParamVals(ppInd,paramInd) = paramVal;
                end
            end
            (nanmean(allParamVals,1))
            (nanstd(allParamVals,1))
        end
    end
end

