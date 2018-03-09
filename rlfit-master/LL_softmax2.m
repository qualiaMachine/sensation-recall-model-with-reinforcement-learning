function LL=LL_softmax2(prA,action1,action2)
%
%   prA = probs.short(par,data)
%   action1 = data$action1
%   action2 = data$action2
  action = [action1;action2];
%   prA_noZeros = 
  prA(prA==1) = 1-1e-7;
  LL = -sum(action'.*log(1 - prA) + (1-action)'.*log(prA));
  if isinf(LL) || isnan(LL)
      keyboard
  end
%   if (is.finite(LL) == F) 
%       LL=100000;
%   end
      

% 1-prA
% % % %calculate log of softmax function at each of the values in xx
% % % %handles cases of extreme values by Taylor series
% % % 
% % % % make choice a column vector
% % % choice = choice(:);
% % % 
% % % hh = exp(xx); %exponentiate values
% % % [N, K] = size(hh); %number of trials and options, respectively
% % % 
% % % %set limits on h
% % % hmax=1e15;
% % % hmin=1e-15;
% % % 
% % % %now do cases
% % % %first, do straightforward cases
% % % norm = sum(hh,2); %normalization (sum each row; i.e., sum over options of exp(xx))
% % % to_grab = sub2ind([N, K], (1:N)', choice); %get linear index for each option chosen
% % % lik = hh(to_grab);
% % % LL = log(lik./norm);
% % % 
% % % %now test limits and fix extreme values:
% % % 
% % % % for very small values of the likelihood, subtract logs of numerator and
% % % % denominator
% % % sel = lik < hmin;
% % % LL(sel) = xx(sel) - log(norm(sel));
% % % 
% % % % for very large likelihoods, we'll taylor expand:
% % % % log(lik/norm)=-log(1+(norm-lik)/lik), using the
% % % % fact that (norm-lik)/lik is very small if lik>>(norm-lik) (i.e., lik/norm~1)
% % % sel = lik > hmax;
% % % tiny = (norm - lik)./lik;
% % % 
% % % h=tiny(sel);
% % % LL(sel)=-h+h.^2/2-h.^3/3+h.^4/4-h.^5/5;
% % % 
% % % if isnan(LL)
% % %     keyboard
% % % end
% % % 
% % % y=LL;
% % %     
end

