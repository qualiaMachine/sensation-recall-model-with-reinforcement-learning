function prA=Q_model2(modelParams,action1,action2,payoff,states,pupilSizes)

pupilSizes(isnan(pupilSizes)) = mean(pupilSizes,'omitnan');
pupilSizes = (pupilSizes-min(pupilSizes))/(max(pupilSizes)-min(pupilSizes));

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
attenWeight1 = modelParams(5); % weight parameter between TD and MB values
% attenWeight2 = modelParams(6); % weight parameter between TD and MB values
attenWeightBinCutoff = modelParams(6); % weight parameter between TD and MB values
p=0;%modelParams(5);

weightsScaledToDistribution = ones(size(pupilSizes));

recallWeights = pupilSizes;
recallWeights=recallWeights*attenWeight1;
recallWeights=weightsScaledToDistribution;

try
    attentionWeights=pupilSizes;
    bin1Inds = find(attentionWeights<attenWeightBinCutoff);
    
    weightsScaledToDistribution(bin1Inds) = attentionWeights(bin1Inds)*attenWeight1;

        if sum(isnan(pupilSizes)) > 0
        keyboard
    end
    if length(pupilSizes) ~= 150 || sum(pupilSizes >1) >0
        keyboard
    end
catch why
    keyboard
end
attentionWeights=weightsScaledToDistribution;

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


for t = 1:nTrials
    if t == 1
        P2A(t) = 0.5;
        P2B(t) = 0.5;
        QA_TD(t,:) = 0;
        QB_TD(t,:) = 0;
    else
        P2A(t) = P2A(t-1); %taking empirical probabilities from the previous trial
        P2B(t) = P2B(t-1);
        QA_TD(t,:) = QA_TD(t-1,:); %taking Q-values for each action from the previous trial
        QB_TD(t,:) = QB_TD(t-1,:);
    end
    
    % model-based Q-values (MB)
    QA_MB(t) = recallWeights(t)*P2A(t)*max(QA_TD(t,2),QB_TD(t,2)) + (1-recallWeights(t)*P2A(t))*max(QA_TD(t,3),QB_TD(t,3));
    QB_MB(t) = recallWeights(t)*P2B(t)*max(QA_TD(t,2),QB_TD(t,2)) + (1-recallWeights(t)*P2B(t))*max(QA_TD(t,3),QB_TD(t,3));
    
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
    
    if (action2(t) == 0) %A was chosen in the second stage
        %updating Q-value of the first stage through eligibility trace
        %(TD)
        if (action1(t) == 0) 
            QA_TD(t,1) = QA_TD(t,1)*(1-(alpha*attentionWeights(t))) + QA_TD(t,i)*alpha*attentionWeights(t) + (payoff(t)-QA_TD(t,i))*alpha*lambda*attentionWeights(t);
        else 
            QB_TD(t,1) = QB_TD(t,1)*(1-(alpha*attentionWeights(t))) + QA_TD(t,i)*alpha*attentionWeights(t) + (payoff(t)-QA_TD(t,i))*alpha*lambda*attentionWeights(t);
        end
      
        %updating Q-value of the second stage action
        QA_TD(t,i) = QA_TD(t,i)*(1-(alpha2*attentionWeights(t))) + payoff(t)*alpha2*attentionWeights(t);
      
      
    else
        %B was chosen in the second state
      
        %updating Q-value of the first stage through eligibility trace
        %(TD)
        if (action1(t) == 0) 
            QA_TD(t,1) = QA_TD(t,1)*(1-(alpha*attentionWeights(t))) + QB_TD(t,i)*alpha*attentionWeights(t) + (payoff(t)-QB_TD(t,i))*alpha*lambda*attentionWeights(t);
        else 
            QB_TD(t,1) = QB_TD(t,1)*(1-(alpha*attentionWeights(t))) + QB_TD(t,i)*alpha*attentionWeights(t) + (payoff(t)-QB_TD(t,i))*alpha*lambda*attentionWeights(t);
        end

        %updating Q-value of the second stage action
        QB_TD(t,i) = QB_TD(t,i)*(1-(alpha2*attentionWeights(t))) + payoff(t)*alpha2*attentionWeights(t);

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