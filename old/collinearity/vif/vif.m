% % function [V]=vif(X)
% % %vif() computes variance inflation coefficients  
% % %VIFs are also the diagonal elements of the inverse of the correlation matrix [1], a convenient result that eliminates the need to set up the various regressions
% % %[1] Belsley, D. A., E. Kuh, and R. E. Welsch. Regression Diagnostics. Hoboken, NJ: John Wiley & Sons, 1980.
% % 
% % R0 = corrcoef(X); % correlation matrix
% % V=diag(inv(R0))';

function V=vif(X)
% V=vif(X)
%
% Variance inflaction factor in Regression Analysis
% the variance inflation factor (VIF) quantifies the severity of 
% multicollinearity in an ordinary least squares regression analysis. 
%It provides an index that measures how much the variance of an estimated 
% regression coefficient is increased because of collinearity.
%
% INPUT:
% 
% X is the matrix n onservation x p variables  
%
% OUTPUT:
%
% V is a column vector of vif indices

[n,p]=size(X);
V=zeros(p,1);
if p>2
    for i=1:p
        pred=setdiff(1:p,i);
        [~, R2 , ~]=ols(X(:,i),X(:,pred),1);
        V(i)=1/(1-R2);
    end
else
    [~, R2 , ~]=ols(X(:,1),X(:,2),1);
    V(1)=1/(1-R2);
    V(2)=V(1);
end
