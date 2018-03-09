function [pupilDiams,allPupDrops] = codePupilResps(two_stage_task_data,condition,pars)
w = warning ('off','all');

% Note: 1st column is timestamp, pupil diameter is 4th col of data
pp = two_stage_task_data.subject;
load(['T:\RLmodel\dataset-eyeTracker-choices\edata' num2str(pp) num2str(condition) '.mat'])

figCount= 0;
plotPupilData = false; % boolean for visually checking performance of blink removal code
pupilDiams = zeros(size(two_stage_task_data.eton,1),3); 

for iTrial = 1:size(two_stage_task_data.eton,1)
    if iTrial > 5
        plotPupilData = false;
    end
    try
        % Have to use two_stage_task_data.eton for stimOnset (but add
        % two_stage_task_data.et0) for absolute time to access proper timestamps of
        % eye-tracker data
        if pars.gatherStimOnset
            % begOrEnd = beg of stim/stage
            begOrEndPeriodOfInt = two_stage_task_data.eton(iTrial,:) + two_stage_task_data.et0;
        elseif ~pars.gatherStimOnset || pars.grabAnticResp
            % begOrEnd = end of stimuli (when keypress happens)
            begOrEndPeriodOfInt = two_stage_task_data.etoff(iTrial,:) + two_stage_task_data.et0;
        end
        for state = 1:3
            %% First Process Stim Onset data
            if isnan(two_stage_task_data.action(iTrial,state))
                pupilDiams(iTrial,state) = nan;
                allPupDrops(iTrial,state) = nan;
            elseif begOrEndPeriodOfInt(state) < two_stage_task_data.et0 
                pupilDiams(iTrial,state) = nan; % indicates there's no eye-tracker data preceding stim (choice for stimOffset data) onset
                allPupDrops(iTrial,state) = nan;
            elseif begOrEndPeriodOfInt(state) + pars.stimDurForEyeMeans  > edata(end,1)
                pupilDiams(iTrial,state) = nan; % indicates not enough forward-going data to collect period of interest (stimDurForEyeMeans)
                allPupDrops(iTrial,state) = nan;
            else
                if pars.stimDurForEyeMeans < 0
                    if pars.gatherStimOnset
                        endPupDataRow = find(edata(:,1) == begOrEndPeriodOfInt(state));
                        begPupDataRow = endPupDataRow+pars.stimDurForEyeMeans;
                    else
                        endPupDataRow = find(edata(:,1) == begOrEndPeriodOfInt(state));
                        begPupDataRow = endPupDataRow+pars.stimDurForEyeMeans;
                        if ~pars.grabAnticResp
                            endPupDataRow = endPupDataRow - pars.antRespLength;
                        end
                    end
                else
                    begPupDataRow = find(edata(:,1) == begOrEndPeriodOfInt(state));
                    endPupDataRow = begPupDataRow+pars.stimDurForEyeMeans;
                end
                if isempty(begPupDataRow) || isempty(endPupDataRow)
                    keyboard
                else
                    try
                        pupilResp = edata(begPupDataRow:endPupDataRow,4);
                    catch
                        pupilDiams(iTrial,state) = nan;
                        allPupDrops(iTrial,state) = nan;
                        continue
                    end
                    if iTrial == 20
    %                     keyboard
                    end
                    if iTrial > 0
                        plotData = false;
                        smoothedResp = smooth(pupilResp,50);
                        veloc = smooth(diff(smoothedResp),10);
                        try
                            [~,blinkOffsets] = findpeaks(veloc,'MinPeakHeight',pars.velThresh,'MinPeakProminence',10,'MinPeakDistance',20);
                            blinkOffsets = blinkOffsets + pars.pkAdjust;
                            [~,blinkOnsets] = findpeaks(veloc*-1,'MinPeakHeight',pars.velThresh,'MinPeakProminence',10,'MinPeakDistance',20);
                            blinkOnsets = blinkOnsets - (pars.pkAdjust/2);
                            [pupVals,pupDrops]=findpeaks(veloc*-1,'MinPeakHeight',3,'MinPeakDistance',100);
                            pupDrops = pupDrops(pupVals<pars.velThresh);
                        catch why
                            keyboard
                        end
                        if ~isempty(blinkOffsets) || ~isempty(blinkOnsets)
                            plotData = true; % only plot data w/ blinks to check reliability of algorithm
                            if plotPupilData
                                blinkRemovalFig=figure;%('units','normalized','outerposition',[0 0 1 1]);

                                figCount = figCount + 1;
                                if figCount == 35 || figCount == 80
        %                             keyboard
                                end
                                numRows=4;
                                subtightplot(numRows,1,1,.05)
%                                 plot(pupilResp)
%                                 hold on
                                plot(smoothedResp,'LineWidth',5)
                                xlim([0,abs(pars.stimDurForEyeMeans)-500])
                                ylim([0,5000])
                                set(gca,'xtick',[])

                                subtightplot(numRows,1,2,.05)
                                plot([1:length(veloc)],veloc,'r','LineWidth',5)
                                ylim([-200,200])
                                set(gca,'xtick',[])
                            else
                                plotData = false;
                            end
                        end

                        if sum(blinkOnsets < 1) > 0 && ~isempty(blinkOffsets)
                            % hits this case if beginning of blink is too
                            % close to start of period of interest
                            pupilResp = pupilResp(blinkOffsets(1):end);
                            if isempty(pupilResp)
                                pupilDiams(iTrial,state) = nan;
                                allPupDrops(iTrial,state) = nan;
                                continue
                            end
                            if length(blinkOnsets) == 1
                                blinkOnsets = [];
                                blinkOffsets = [];
                            else
                                blinkOnsets = blinkOnsets(2:end)-blinkOffsets(1);
                                blinkOnsets = blinkOnsets(blinkOnsets>0); % in case no endBlink pk is found
                                blinkOffsets = blinkOffsets(2:end)-blinkOffsets(1);
                            end
                        elseif sum(blinkOnsets < 1) > 0 && isempty(blinkOffsets)
                            % eyes open very slowly, no offset detected;
                            % just leave data as is and get rid of onset
%                             if length(blinkOnsets) == 1
                            blinkOnsets = [];
                            blinkOffsets = [];
%                             else
%                                 blinkOnsets = blinkOnsets(2:end)-blinkOffsets(1);
%                                 blinkOnsets = blinkOnsets(blinkOnsets>0); % in case no endBlink pk is found
%                                 blinkOffsets = blinkOffsets(2:end);
%                             end
                        end
                        endTrialOffsets = find(blinkOnsets > (abs(pars.stimDurForEyeMeans)-100)); % doesn't catch blinks that last longer
                        if length(blinkOnsets) ~= length(blinkOffsets)
                            if length(endTrialOffsets) == 1 % indicates that blink overlaps with end of period of interest
                                pupilResp = pupilResp(1:blinkOnsets(end));
                                if length(blinkOnsets) == 1
                                    blinkOnsets = [];
                                    blinkOffsets = [];
                                else
                                    blinkOnsets = blinkOnsets(1:end-1);
                                end
                            else
                                if length(blinkOnsets) - length(blinkOffsets) == 1
                                    % hits this case if no endTrialOffset
                                    % detected due to long-lasting blink at
                                    % end of trial
                                    pupilResp = pupilResp(1:blinkOnsets(end));
                                    blinkOnsets = blinkOnsets(1:end-1);

                                elseif isempty(blinkOnsets)
                                    % hits this case if no beginning blink
                                    % detected (so no < 1 index), but still
                                    % have an offset index
                                    
                                    % We can add decent guesses for Onsets
                                    for newOnset = 1:length(blinkOffsets)
                                        offset = blinkOffsets(newOnset);
                                        onset = offset-100;
                                        if onset < 1
                                            onset = 1;
                                        end
                                        blinkOnsets(newOnset) = onset;
                                    end
                                end
                            end
                        end
                        %% After dealing with onset/offset cases, deal with rapid blinks
                        if length(blinkOnsets) ~= length(blinkOffsets)
                            % hits this case sometimes when there's
                            % rapid succession of blinks; solution:
                            % pair each onset with nearest offset and
                            % ignore remainders
                            if length(blinkOnsets) > length(blinkOffsets)
                                if length(blinkOffsets) > 1
                                    x = 1; % do nothing but check out that code works
                                end
                                diffs = [];
                                for blinkOn = 1:length(blinkOnsets)
                                    diffs(blinkOn,:) = blinkOnsets(blinkOn) - blinkOffsets;
                                end
                                pairedOnsets = [];
                                for pairedOffset = 1:length(blinkOffsets)
                                    onsetDists = diffs(:,pairedOffset);
                                    onsetDists(onsetDists > 0) = -inf; %
                                    onsetDists(onsetDists > -100) = -inf; % suggests too rapid of blink to be considered discrete from first blink
                                    [~,closestOnset] = max(onsetDists);
                                    pairedOnsets(pairedOffset) = closestOnset;
                                end
                                blinkOnsets = blinkOnsets(pairedOnsets);
                                if length(blinkOnsets) ~= length(blinkOffsets)
                                    keyboard
                                end
                            elseif length(blinkOffsets) > length(blinkOnsets)
                                diffs = [];
                                for blinkOff = 1:length(blinkOffsets)
                                    diffs(blinkOff,:) = blinkOffsets(blinkOff) - blinkOnsets;
                                end
                                pairedOffsets = [];
                                for pairedOffset = 1:length(blinkOnsets)
                                    offsetDists = diffs(:,pairedOffset);
                                    offsetDists(offsetDists < 0) = inf; %
                                    offsetDists(offsetDists < 100) = inf; % suggests too rapid of blink to be considered discrete from first blink
                                    [~,closestOffset] = min(offsetDists);
                                    pairedOffsets(pairedOffset) = closestOffset;
                                end
                                blinkOffsets = blinkOffsets(pairedOffsets);
                                if length(blinkOffsets) ~= length(blinkOffsets)
                                    keyboard
                                end
                            else
                                keyboard
                            end
                        end
                        x = 1:length(pupilResp);
                        y = pupilResp;
                        tooFarOffsets = find(blinkOffsets > length(pupilResp));
                        blinkOffsets(tooFarOffsets) = length(pupilResp);
                        tooNearOnsets = find(blinkOnsets < 1);
                        blinkOnsets(tooNearOnsets) = 1;
                        for blinkInd = 1:length(blinkOnsets)
                            y(blinkOnsets(blinkInd)+1:blinkOffsets(blinkInd)-1) = nan;
                        end
                        x = x(~isnan(y));
                        y = y(~isnan(y));
                        if plotData
                            subtightplot(numRows,1,3,.05)
                            scatter(x,y)
                            xlim([0,length(pupilResp)])
                            ylim([0,5000])
                            set(gca,'xtick',[])
                        end
                        try
                            [F] = interp1(x,y,[1:length(pupilResp)],'pchip');
                        catch why
                            keyboard
                        end
                        if plotData
                            subtightplot(numRows,1,4,.05)
                            plot([1:length(F)],smooth(F,50),'LineWidth',5)
                            xlim([0,length(F)])
                            ylim([0,5000])
                            saveas(gcf,[pars.figDir filesep 'blinkRemoval_pp' num2str(pp) 'trial_' num2str(iTrial) '.png'])
                            close all
                        end
                        pupilDiams(iTrial,state) = mean(F);
                        allPupDrops(iTrial,state) = length(pupDrops);
                        if pupilDiams(iTrial,state) == 0
                            pupilDiams(iTrial,state)=NaN;
                        end
                    end
                end
            end
        end
    catch why
        keyboard
    end
end
% keyboard
%% normalize values based on average dilation from stimulus
% Skipping this step because there are absolute differences betwen pps I
% want to caputre
zeroInds = find(pupilDiams==0);
if ~isempty(zeroInds)
    keyboard
end
% % for state = 1:3
% % %     meanDil = mean(pupilDiams(1:20,state),'omitnan');
% % %     stdDil = std(pupilDiams(1:20,state),'omitnan');
% % %     z = bsxfun(@minus, pupilDiams(:,state), meanDil)/stdDil;
% % %     pupilDiams(:,state) = z;
% %     pupilDiams(:,state) = (pupilDiams(:,state) - min(pupilDiams(:,state))) / (max(pupilDiams(:,state))-min(pupilDiams(:,state)));
% % end
if length(pupilDiams) ~= 150
    keyboard
end



