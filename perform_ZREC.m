
% added W to output
function [Q,Avg_CI_size,CI_ZREC,W]=perform_ZREC(O)
    
    % O is the matrix of opinion scores
    n_pvs=size(O,1);
    n_sub=size(O,2);
    m=(mean(O,2));
    s=(std(O',0))';
    I=ones(1,n_sub);
    Z=(O-m*I)./(s*I);
    
    % compute the biases
    B=mean(Z,1);
    
    
    % compute inconsistencies
    for i=1:n_sub
        C(i)=sqrt(mean((Z(:,i)-B(i)).^2));
    end
    
    
    U=O-s*B;
    
    W=(C.^(-2))/sum((C.^(-2)));
    Q=U*W'; % Q is the R in the paper

    % compute CI avg size
    S=sqrt(  (n_sub/(n_sub-1))*  (((U-Q*ones(1,size(U,2))).^2)*W'));
    
    % mean intereval range
    Avg_CI_size=mean(2*1.96*S/sqrt(n_sub) );
    
    % 95% confidence interval
    CI_ZREC=[Q-1.96*S/sqrt(n_sub), Q+1.96*S/sqrt(n_sub)];
end