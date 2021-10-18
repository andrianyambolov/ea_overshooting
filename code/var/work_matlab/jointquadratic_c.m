% JOINTQUADRATIC_C.M

function [irfquad,credibleset]=jointquadratic_c(IRF)

M=size(IRF,1);
c=log(M)/sqrt(M);
quadratic_loss = zeros(M,1);
for i=1:M
   for j=1:M
       quadratic_loss(i,1) = quadratic_loss(i,1) + (IRF(j,:)-IRF(i,:))*(IRF(j,:)-IRF(i,:))';
   end
end
quadratic_loss = quadratic_loss/M;
[a,I] = sort(quadratic_loss);
quadratic_loss = 100*(quadratic_loss-a(1))/a(1);
irfquad  = IRF(quadratic_loss<c,:);
credibleset = IRF(I(1:floor(0.68*M)+1),:);

