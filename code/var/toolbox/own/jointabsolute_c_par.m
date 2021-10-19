% JOINTABSOLUTE_C.M
% 
% Bayes estimate and joint credible set for any M x nirf matrix IRF, where 
% M is the number of posterior draws and nirf is the number of impulse 
% responses. 

function [irfabs,credibleset]=jointabsolute_c_par(IRF)
M= size(IRF,1);
c=log(M)/sqrt(M);
absolute_loss = zeros(M,1);
parfor i=1:M
   for j=1:M
%        currentJs = [currentJs, sum(abs(IRF(j,:)-IRF(i,:)))];
       absolute_loss(i,1) = absolute_loss(i,1)+sum(abs(IRF(j,:)-IRF(i,:)));
   end
%    absolute_loss(i,1) = sum(parallel.pool.Constant(currentJs));
end
absolute_loss = absolute_loss/M;
[a,I] = sort(absolute_loss);
absolute_loss = 100*(absolute_loss-a(1))/a(1);
irfabs  = IRF(absolute_loss<c,:);
credibleset = IRF(I(1:floor(0.68*M)+1),:);
