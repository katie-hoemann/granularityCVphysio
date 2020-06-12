function [R2change,Fchange,pchange,df1,df2] = twostep_hierarchical_regression(y,X,numVars)
% performs a hierarchical linear regression and derives change
% scores associated with adding the last variable(s) in X as step 2
% assumes that X is padded with a column of 1s for the intercept and
% contains variables in order of model entry
% numVars = number of variables to add in step 2

X_step1 = X(:,1:end-numVars); 
X_step2 = X;
idx = find(isnan(X));
[row,~] = ind2sub(size(X),idx);
X_step1(row,:) = [];
X_step2(row,:) = [];
y(row) = [];
[~,~,~,~,stats_step1] = regress(y,X_step1);
[~,~,~,~,stats_step2] = regress(y,X_step2);

R2_step1 = stats_step1(1);
R2_step2 = stats_step2(1);
[~,IVs_step2] = size(X_step2);
IVs_step2 = IVs_step2-1;
df1 = numVars; 
df2 = length(X)-IVs_step2-1; % N-k-1
R2change = R2_step2-R2_step1;
Fchange = ((R2_step2-R2_step1)/df1)/((1-R2_step2)/df2);
pchange = 1-fcdf(Fchange,df1,df2);
end