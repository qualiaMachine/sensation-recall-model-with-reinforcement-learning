function [beta, LL, Q] = rlfit3(Qfun, action1,action2,states,payoff,attentionWeights, lb, ub, niter)
% fits a reinforcement learning model to a multi-option choice paradigm
% inputs:
%
% Qfun is a handle to a function that accepts a vector of parameters, a
% vector of choice indices, and a vector of outcomes, and returns 
% the action values, Q
% 
% choice is a vector, one entry per trial, the index of the chosen option
%
% outcome is a set of outcomes for each trial
%
% lb and ub are vectors of upper and lower bounds on parameters
%
% niter (optional) is the number of random restarts to use in fitting
% 
% outputs:
%
% beta is the vector of fitted model parameters; first entry is the softmax
% inverse temperature, followed by parameters of the model
%
% LL is the log likelihood of the data (choice, outcome) given beta
% 
% Q is a trials x options matrix of action values


if ~exist('niter', 'var')
    niter = 10;
end

if ~exist('lb', 'var')
    lb = [];
end

if ~exist('ub', 'var')
    ub = [];
end

% rescale outcomes to offer better fit convergence
% outmean = mean(payoff(:));
% outstd = std(payoff(:));
% z = bsxfun(@minus, payoff, outmean)/outstd;
z=payoff;

% first, define a log likelihood function that takes as its input a vector
% of parameters, the first of which is the inverse temperature of the
% softmax
% LLfun = @(x, action1, action2,states,z) LL_softmax2(Qfun(x(1:end),action1,action2,z,states,attentionWeights), action1,action2);
LLfun = @(x, action1, action2,states,z,attentionWeights) LL_softmax(x(1)*Qfun(x,action1,action2,z,states,attentionWeights),states, [action1;action2]);

% then define a function to be minimized (the total negative log
% likelihood)
fitfun = @(beta)(-1)*sum(LLfun(beta, action1,action2,states, z,attentionWeights));

% now combine upper and lower bounds on softmax temp with upper and lower
% bounds on other parameters
lb = [1e-1, lb]; %lower bounds
ub = [10, ub]; %upper bounds

% optmize to fit model
w = warning ('off','all');
options = optimset('Display', 'off');
try
    [beta,fval]=multmin(fitfun, lb, ub, niter, options);
catch why
    keyboard
end
warning(w);

% return log likelihood
LL=-fval;

% get action values
% keyboard
% 0.1229626 0.7079729 0.0000000 0.1685441 0.5216178
% alpha = par[1] #learning rate
% beta = par[2] #inverse temperature (softmax choice parameter)
% alpha2 = alpha #learning rate
% beta2 = beta
% lambda = par[3] #eligibility trace decay parameter
% p = par[4] #stickiness of choice parameter
% w = par[5] # weight parameter between TD and MB values
% beta(1) = .7079729;
% beta(2) = .1229626;
% beta(3) = 0;
% beta(4) = .5216178;
% beta(5) = .1685441;
% % beta = modelParams(1); %inverse temperature (softmax choice parameter)
% % alpha = modelParams(2);
% % alpha2 = alpha; %learning rate
% % beta2 = beta;
% % lambda = modelParams(3); %eligibility trace decay parameter
% % % p = modelParams(4); %stickiness of choice parameter
% % w = c(4); % weight parameter between TD and MB values
% % attenWeight = 1;%modelParams(5); % weight parameter between TD and MB values
% % attenWeightBinCutoff = inf;%modelParams(6); % weight parameter between TD and MB values
% % p=modelParams(5);

Q = Qfun(beta(1:end), action1,action2,z,states,attentionWeights);
% x,action1,action2,z,states,attentionWeights
% undo scaling
% Q = Q*outstd + outmean; % rescale appropriately
% beta(1) = beta(2)/outstd;

