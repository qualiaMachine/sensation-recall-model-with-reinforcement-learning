function prA=Q_model10(modelParams,action1,action2,payoff,states,pupilSizes)

% pupilSizes = ((pupilSizes-0)/(5000));
pupilSizes = (pupilSizes-min(pupilSizes))/(max(pupilSizes)-min(pupilSizes));
pupilSizes=pupilSizes.^4;
if sum(isnan(pupilSizes)) > 10
    keyboard
end
pupilSizes(pupilSizes==0)=nan;
% pupilSizes(isnan(pupilSizes)) = mean(pupilSizes,'omitnan');

% basic reinforcement learning model
% calculate the action values, given model parameters
% alpha is the learning rate
nTrials = length(payoff);
if sum(isnan(modelParams)) > 0
    keyboard
end
%% parameters
% lb = [beta,0,0,-Inf,0];%,1000,0,0]; %lower bounds
% ub = [beta,1,1,+Inf,1];%,3500,1,1]; %upper bounds
beta = modelParams(1); %inverse temperature (softmax choice parameter)
alpha = modelParams(2);
alpha2 = alpha; %learning rate
beta2 = beta;
lambda = modelParams(3); %eligibility trace decay parameter
% p = modelParams(4); %stickiness of choice parameter
w = modelParams(4); % weight parameter between TD and MB values
% pupilRecallEst=modelParams(5);
pupilExploreCoef = 0;%modelParams(5); % weight parameter between TD and MB values
pupilExploitCoef = 0;%modelParams(6); % weight parameter between TD and MB values
memoryBooleanThresh = modelParams(5);
recentMemThresh=modelParams(6);
% wInteraction = modelParams(7);
% attenWeight2 = modelParams(6); % weight parameter between TD and MB values
% attenWeightBinCutoff = modelParams(6); % weight parameter between TD and MB values
p=0;%modelParams(5);

weightsScaledToDistribution = ones(size(pupilSizes));
% recallWeights = pupilSizes;
% bin1Inds = find(recallWeights<attenWeightBinCutoff);
% weightsScaledToDistribution(bin1Inds) = recallWeights(bin1Inds)*attenWeight1;
recallWeights=weightsScaledToDistribution;

% recallWeights=recallWeights*attenWeight1;
% recallWeights=weightsScaledToDistribution;

try
    weightsScaledToDistribution = ones(size(pupilSizes));
    attentionWeights=pupilSizes;
    attentionWeights=weightsScaledToDistribution;
%     bin1Inds = find(attentionWeights<attenWeightBinCutoff);
%     weightsScaledToDistribution(bin1Inds) = attentionWeights(bin1Inds)*attenWeight1;

%     if sum(isnan(pupilSizes)) > 0
%         keyboard
%     end
    if length(pupilSizes) ~= 150 || sum(pupilSizes >1) >0
        keyboard
    end
catch why
    keyboard
end
% attentionWeights=weightsScaledToDistribution;
% attentionWeights=attenWeight1*attentionWeights;
alphaOrig=alpha;
betaOrig=beta;
lambdaOrig=lambda;
%%
% action1 = data$action1
% action2 = data$action2
% trials = length(data$action1) 
% payoff = data$payoff
% state = data$state


%% define transition prob.
N = zeros(2,2);% observed number of transitions from L (row 1) to 2 or 3 (columns)
% or R (row 2) to 2 or 3 (columns)
P2A = ones(1,length(payoff))/2; %observed probability of transition to 2 after A
P2B = ones(1,length(payoff))/2; %observed probability of transition to 2 after B

prA1 = zeros(1,nTrials);
prA2 = zeros(1,nTrials);

%% Init Q-vals for model-free portion
QA_TD = zeros(nTrials+1,3); %matrix for Q-values of A choice in 3 states
QB_TD = zeros(nTrials+1,3); %matrix for Q-values of B choice in 3 states
  
%% Init Q-vals for model-based portion
QA_MB = zeros(1,nTrials+1);%vector for Q-values for the left choice in state 1
QB_MB = zeros(1,nTrials+1);%vector for Q-values for the right choice in state 1

%% Init net Q-vals
QA_NET = zeros(1,nTrials+1);%vector for Q-values for the left choice in state 1
QB_NET = zeros(1,nTrials+1);%vector for Q-values for the right choice in state 1

% %% init values for Q-learning
% N=length(payoff); %number of trials
% a=length(unique(choice)); %number of options
% Q=nan(N,a,s); %values of each choice&state of each trial

memoryBuffer=[];
for t = 1:nTrials
    if t > 1
        origQA_TD = QA_TD; %matrix for Q-values of A choice in 3 states
        origQB_TD = QB_TD; %matrix for Q-values of B choice in 3 states

        %% Init Q-vals for model-based portion
        origQA_MB = QA_MB;%vector for Q-values for the left choice in state 1
        origQB_MB = QB_MB;%vector for Q-values for the right choice in state 1

        %% Init net Q-vals
        origQA_NET = QA_NET;%vector for Q-values for the left choice in state 1
        origQB_NET = QB_NET;%
    
        qDiff=QA_TD-QB_TD;
        for stateInd = 1:3
            if qDiff(t-1,stateInd) > 0 
                aBest(stateInd)=1;
                bBest(stateInd)=0;
            else
                aBest(stateInd)=0; 
                bBest(stateInd)=1;
            end
        end
%         if t < 50
        if (action1(t) == 0) %A was chosen in the 1st stage
            if aBest(1) >0
                exploit=true;
            else
                exploit=false;
            end
        end
%         else
%             if (action2(t) == 0) %A was chosen in the second stage
%                 if aBest() >0
%                     exploit=true;
%                 else
%                     exploit=false;
%                 end
%             end
%         end
    end
        
    origT=t;
%     alpha = alphaOrig*attentionWeights(t);
%     alpha2 = alphaOrig*attentionWeights(t); %learning rate
%     beta=betaOrig*log(attentionWeights(t));
%     beta2=betaOrig*(attentionWeights(t));

%     beta2 = beta;
    % lambda = modelParams(3); %eligibility trace decay parameter
    % p = modelParams(4); %stickiness of choice parameter
    % w = modelParams(4); 
%     lambda=lambdaOrig*attentionWeights(t);
    maxRecall = 30;
    if ~isnan(pupilSizes(t)) && pupilSizes(t) < memoryBooleanThresh && ~isempty(memoryBuffer) && exploit
        tInds = [];
        if length(memoryBuffer) >=maxRecall
            maxInds = maxRecall;
        else
            maxInds = length(memoryBuffer);
        end
        try
            tInds=memoryBuffer(end-maxInds+1:end);
        catch why
            keyboard
        end
%         for randInd = 1:maxInds
        retaintedTInds=[];
        for tInd = tInds
            if pupilSizes(tInd) < recentMemThresh
%                 tInds = [tInds, randi([1 length(memoryBuffer)],1)];
                retaintedTInds = [retaintedTInds,tInd];
            end
        end
        tInds=[retaintedTInds,t];
        if length(memoryBuffer) > maxRecall
            memoryBuffer=memoryBuffer(end-maxRecall:end); % keep most recent 15 timesteps
        end
    else
        tInds = t;
    end
    currT=tInds(end);
    for t = tInds
        if t == 1
            P2A(t) = 0.5;
            P2B(t) = 0.5;
            QA_TD(t,:) = 0;
            QB_TD(t,:) = 0;
        else
            P2A(t) = P2A(currT-1); %taking empirical probabilities from the previous trial
            P2B(t) = P2B(currT-1);
            QA_TD(t,:) = QA_TD(currT-1,:); %taking Q-values for each action from the previous trial
            QB_TD(t,:) = QB_TD(currT-1,:);
        end
    %     diffToPt5 = .5-P2A(t);
    %     P2A(t) = P2A(t)+(diffToPt5*(pupilSizes(t)/2)*pupilRecallEst);
    %     diffToPt5 = .5-P2B(t);
    %     P2B(t) = P2B(t)+(diffToPt5*(pupilSizes(t)/2)*pupilRecallEst);

        % model-based Q-values (MB)
    %     QA_MB(t) = recallWeights(t)*P2A(t)*max(QA_TD(t,2),QB_TD(t,2)) + (1-(recallWeights(t)*P2A(t)))*max(QA_TD(t,3),QB_TD(t,3));
    %     QB_MB(t) = recallWeights(t)*P2B(t)*max(QA_TD(t,2),QB_TD(t,2)) + (1-(recallWeights(t)*P2B(t)))*max(QA_TD(t,3),QB_TD(t,3));

        QA_MB(t) = P2A(t)*max(QA_TD(t,2),QB_TD(t,2)) + (1-(P2A(t)))*max(QA_TD(t,3),QB_TD(t,3));
        QB_MB(t) = P2B(t)*max(QA_TD(t,2),QB_TD(t,2)) + (1-(P2B(t)))*max(QA_TD(t,3),QB_TD(t,3));

        %weighted values
        if ~isnan(w*QA_MB(t) + (1-w)*QA_TD(t,1))
            QA_NET(t) = w*QA_MB(t) + (1-w)*QA_TD(t,1);
        else
            keyboard
        end
        if ~isnan(w*QB_MB(t) + (1-w)*QB_TD(t,1))
            QB_NET(t) = w*QB_MB(t) + (1-w)*QB_TD(t,1);
        else
            keyboard
        end

        %checking if sticky choice
        if t == 1
          sA = 0;
          sB = 0;
        else
            if (action1(t-1) == 0)
                sA = 1;
                sB = 0;
            else
                sB = 1;
                sA = 0;
            end
        end

        %first stage choice (state 1)
        prA1(t) = exp(beta*(QA_NET(t) + p*sA))/(exp(beta*(QA_NET(t)+ p*sA)) + exp(beta*(QB_NET(t)+ p*sB)));
        if prA1(t) > 1
            keyboard
        end
        %first-stage  left or right choice and updating observed number of transitions (MB)
        if (action1(t) == 0) 
            if (states(t) == 2)
                i = 2;
                N(1,1) = N(1,1) + 1;
            else 
                i = 3;
                N(1,2) = N(1,2) + 1;
            end
        else 
            if (states(t) == 3)
                i = 3; 
                N(2,2) = N(2,2) + 1;
            else
                i = 2;
                N(2,1) = N(2,1) + 1;
            end
        end


        %updating transition probabilities (MB)
        P2A(t) = (1 + N(1,1))/(2 + N(1,1) + N(1,2));
        P2B(t) = (1 + N(2,1))/(2 + N(2,1) + N(2,2));


        %second stage choice 
        prA2(t) = exp(beta2*QA_TD(t,i))/((exp(beta2*QA_TD(t,i)) + exp(beta2*QB_TD(t,i))));
        if prA1(t) > 1
            keyboard
        end
        qDiff=QA_TD-QB_TD;
        for stateInd = 1:3
            if qDiff(t,stateInd) > 0 
                aBest(stateInd)=1;
                bBest(stateInd)=0;
            else
                aBest(stateInd)=0; 
                bBest(stateInd)=1;
            end
        end
    %     paramToModulate=lambda;
        origParam = alphaOrig;
        if (action2(t) == 0) %A was chosen in the second stage
            if aBest(1) >0
                exploit=true;
                alpha = (origParam) + (attentionWeights(t)*pupilExploreCoef);
    %             alpha = (origParam/2)+(w*wInteraction*exploit*attentionWeights(t)*pupilExploreCoef*.5);
            else
                exploit=false;
    %             lambda = (origParam/2) + (attentionWeights(t)*pupilExploitCoef);
    %             alpha = (origParam)*(attentionWeights(t)*pupilExploitCoef);
    %             alpha = (origParam/2)+(w*exploit*attentionWeights(t)*pupilExploitCoef*.5);
    %             alpha = (origParam/2)+(w*wInteraction*exploit*attentionWeights(t)*pupilExploitCoef*.5);
                alpha = (origParam) + (attentionWeights(t)*pupilExploitCoef);

            end
            %updating Q-value of the first stage through eligibility trace
            %(TD)
            if (action1(t) == 0) 
                QA_TD(t,1) = QA_TD(t,1)*(1-(alpha)) + QA_TD(t,i)*alpha + (payoff(t)-QA_TD(t,i))*alpha*lambda;
            else 
                QB_TD(t,1) = QB_TD(t,1)*(1-(alpha)) + QA_TD(t,i)*alpha + (payoff(t)-QA_TD(t,i))*alpha*lambda;
            end

            if aBest(i) >0
                exploit=true;
    %             alpha2 = (origParam/2) + (exploit*attentionWeights(t)*pupilExploreCoef);
    %             alpha2 = (origParam) * (exploit*attentionWeights(t)*pupilExploreCoef);
    %             alpha2 = (origParam/2) + (w*exploit*attentionWeights(t)*pupilExploreCoef*.5);
    %             alpha2 = (origParam/2)+(w*wInteraction*exploit*attentionWeights(t)*pupilExploreCoef*.5);
                alpha2 = (origParam) + (attentionWeights(t)*pupilExploreCoef);

            else
                exploit=false;
    %             alpha2=origParam;
    %             alpha2 = (origParam/2) + (attentionWeights(t)*pupilExploitCoef);
    %             alpha2 = (origParam) * (attentionWeights(t)*pupilExploitCoef);
    %             alpha2 = (origParam/2) + (w*exploit*attentionWeights(t)*pupilExploitCoef*.5);
    %             alpha2 = (origParam/2)+(w*wInteraction*exploit*attentionWeights(t)*pupilExploitCoef*.5);
                alpha2 = (origParam) + (attentionWeights(t)*pupilExploitCoef);

            end
            %updating Q-value of the second stage action
            QA_TD(t,i) = QA_TD(t,i)*(1-(alpha2)) + payoff(t)*alpha2;


        else
            if bBest(1) >0
                exploit=true;
    %             alpha=(origParam/2)+(exploit*attentionWeights(t)*pupilExploreCoef);
    %             alpha=(origParam/2)+(w*attentionWeights(t)*pupilExploreCoef*.5);
    %             alpha=(origParam/2)+(w*wInteraction*attentionWeights(t)*pupilExploreCoef*.5);
                alpha = (origParam) + (attentionWeights(t)*pupilExploreCoef);

            else
                exploit=false;
    %             alpha = (origParam/2) + (w*attentionWeights(t)*pupilExploitCoef*.5);
    %             alpha = (origParam/2) + (w*wInteraction*attentionWeights(t)*pupilExploitCoef*.5);
                alpha = (origParam) + (attentionWeights(t)*pupilExploitCoef);

            end
            %B was chosen in the second state

            %updating Q-value of the first stage through eligibility trace
            %(TD)
            if (action1(t) == 0) 
                QA_TD(t,1) = QA_TD(t,1)*(1-(alpha)) + QB_TD(t,i)*alpha + (payoff(t)-QB_TD(t,i))*alpha*lambda;
            else 
                QB_TD(t,1) = QB_TD(t,1)*(1-(alpha)) + QB_TD(t,i)*alpha + (payoff(t)-QB_TD(t,i))*alpha*lambda;
            end

            if bBest(i) >0
                exploit=true;
    %             alpha2=(origParam/2)+(w*attentionWeights(t)*pupilExploreCoef*.5);
    %             alpha2=(origParam/2)+(w*wInteraction*attentionWeights(t)*pupilExploreCoef*.5);
                alpha2 = (origParam) + (attentionWeights(t)*pupilExploreCoef);

            else
                exploit=false;
    %             alpha2=origParam;
    %             alpha2 = (origParam/2) + (w*attentionWeights(t)*pupilExploitCoef*.5);
    %             alpha2 = (origParam/2) + (w*wInteraction*attentionWeights(t)*pupilExploitCoef*.5);
                alpha2 = (origParam) + (attentionWeights(t)*pupilExploitCoef);

            end
            %updating Q-value of the second stage action
            QB_TD(t,i) = QB_TD(t,i)*(1-(alpha2)) + payoff(t)*alpha2;
        end
    end
    %% replace Q's preceding current t with original values
    t=origT;

    if t > 1
        QA_TD(1:t-1,:) = origQA_TD(1:t-1,:) ; %matrix for Q-values of A choice in 3 states
        QB_TD(1:t-1,:) = origQB_TD(1:t-1,:) ; %matrix for Q-values of B choice in 3 states

        %% Init Q-vals for model-based portion
        QA_MB(1,1:t-1) = origQA_MB(1,1:t-1) ;%vector for Q-values for the left choice in state 1
        QB_MB(1,1:t-1) = origQB_MB(1,1:t-1) ;%vector for Q-values for the right choice in state 1

        %% Init net Q-vals
        QA_NET(1,1:t-1) = origQA_NET(1,1:t-1) ;%vector for Q-values for the left choice in state 1
        QB_NET(1,1:t-1) = origQB_NET(1,1:t-1) ;%
    end
    
    if t~=1
        memoryBuffer=[memoryBuffer,t];
    end
end
prA = [prA1,prA2];
if sum(isnan(QA_NET)) > 0 || sum(isnan(QB_NET)) > 0 
    keyboard
elseif sum(isinf(QA_NET)) > 0 || sum(isinf(QB_NET)) > 0 
    keyboard
end
stageOne=[QA_NET;QB_NET]';% firstStage q-vals
stageOne = stageOne(1:150,:);
stageTwo=[];
for trial = 1:150
    try
        stageTwo(trial,1) = QA_TD(trial,states(trial));
        stageTwo(trial,2) = QB_TD(trial,states(trial));
    catch why
        keyboard
    end
end
        
prA=[stageOne;stageTwo];
% prA = [prA;QA_TD,QB_TD];
% prA=[QA_TD(1:150,:);QB_TD(1:150,:)];

% % % % % Q(1,:,:) = 0; %initialize guesses
% % % % % 
% % % % % for ind = 1:(N - 1) 
% % % % %     % copy forward action values to next trial
% % % % %     Q(ind + 1, :,:) = Q(ind,:,:);
% % % % % 
% % % % %     % update option chosen on this trial for next trial's choice
% % % % %     Q(ind + 1,choice(ind)) = Q(ind,choice(ind)) + attentionWeights(ind)*alpha*(outcome(ind)-Q(ind,choice(ind)));     
% % % % % end
% % % % % 
% % % % % %return vector of action values for each trial
% % % % % out=Q;
% % % % % 
% % % % % % % % maxI = 10;
% % % % % % % % for ind = 1:(N - 1) 
% % % % % % % % %     if ind == 1
% % % % % % % % %         keyboard
% % % % % % % % %     end
% % % % % % % %     
% % % % % % % %     for i = 1:maxI
% % % % % % % %         % copy forward action values to next trial
% % % % % % % %         Q(ind + 1, :) = Q(ind, :);
% % % % % % % %         try
% % % % % % % %         % update option chosen on this trial for next trial's choice
% % % % % % % %             Q(ind + 1,choice(ind)) = Q(ind,choice(ind)) + alpha*(outcome(ind)-Q(ind,choice(ind)));     
% % % % % % % %         catch why
% % % % % % % %             keyboard
% % % % % % % %         end
% % % % % % % %         if i~= maxI
% % % % % % % %             Q(ind, :) = Q(ind+1, :);
% % % % % % % %         end
% % % % % % % %     end
% % % % % % % % end