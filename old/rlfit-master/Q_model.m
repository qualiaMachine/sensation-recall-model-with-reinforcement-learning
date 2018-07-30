function out=Q_model(modelParams,choice,outcome,attentionWeights)
% basic reinforcement learning model
% calculate the action values, given model parameters
% alpha is the learning rate
N=length(outcome); %number of trials
k=length(unique(choice)); %number of options
Q=nan(N,k); %values of each choice each trial
alpha = modelParams(1);
weightsScaledToDistribution = zeros(size(attentionWeights));
attenWeight = 1;%modelParams(5); % weight parameter between TD and MB values
attenWeightBinCutoff = inf;%
try
    bin1Inds = find(attentionWeights<attenWeightBinCutoff);
    weightsScaledToDistribution(bin1Inds) = attenWeight;
    weightsScaledToDistribution(weightsScaledToDistribution==0) = 1;
%     bin2Inds = find(attentionWeights<modelParams(4) & weightsScaledToDistribution==0);
%     weightsScaledToDistribution(bin2Inds) = modelParams(5);
%     bin3Inds = find(weightsScaledToDistribution==0);
%     weightsScaledToDistribution(bin3Inds) = modelParams(6);
%     bin3Inds = find(attentionWeights<distCuts(3) & weightsScaledToDistribution==0);
%     weightsScaledToDistribution(bin3Inds) = .4;
%     bin4Inds = find(attentionWeights>distCuts(3) & weightsScaledToDistribution==0);
%     weightsScaledToDistribution(bin4Inds) = 1;
    attentionWeights = weightsScaledToDistribution;
catch why
    keyboard
end


Q(1,:) = 0; %initialize guesses

for ind = 1:(N - 1) 
    % copy forward action values to next trial
    Q(ind + 1, :) = Q(ind,:);

    % update option chosen on this trial for next trial's choice
    Q(ind + 1,choice(ind)) = Q(ind,choice(ind)) + attentionWeights(ind)*alpha*(outcome(ind)-Q(ind,choice(ind)));     
end

%return vector of action values for each trial
out=Q;

% % % maxI = 10;
% % % for ind = 1:(N - 1) 
% % % %     if ind == 1
% % % %         keyboard
% % % %     end
% % %     
% % %     for i = 1:maxI
% % %         % copy forward action values to next trial
% % %         Q(ind + 1, :) = Q(ind, :);
% % %         try
% % %         % update option chosen on this trial for next trial's choice
% % %             Q(ind + 1,choice(ind)) = Q(ind,choice(ind)) + alpha*(outcome(ind)-Q(ind,choice(ind)));     
% % %         catch why
% % %             keyboard
% % %         end
% % %         if i~= maxI
% % %             Q(ind, :) = Q(ind+1, :);
% % %         end
% % %     end
% % % end