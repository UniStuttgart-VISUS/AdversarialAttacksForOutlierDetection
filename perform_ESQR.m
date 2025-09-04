function [q,Avg_CI_size,CIs,weights,Reliability]= perform_ESQR(X)

         % getting the total number of annotators
         n=size(X,2);

         % computing the ESQR contribution of each subject to each stimulus quality 
         [weights,Reliability]=compute_ESQR_Weights(X); 

         % geting the ESQR recovered quality 
         q=diag(X*weights');


         % compute the ESQR CIs and their avg size
         Std=sqrt(  (n/(n-1))*  diag(((X-q*ones(1,size(X,2))).^2)*weights'));
         CIs=[q, q]+[-ones(size(q)),ones(size(q))].*[1.96*Std/sqrt(n),1.96*Std/sqrt(n)];
         Avg_CI_size=mean(2*1.96*Std/sqrt(n) );

end


%% This funtion conpute the reliability measure and thus the weights of each annotator for each stimulus 

function [w,Relaibility]=compute_ESQR_Weights(X)

        % getting the number of annotators and objects
        num_annotators=size(X,2);
        num_objects=size(X,1);
     
        % Estimate the  ground truth distribtions P_V_i
         C=corr(X,type="Spearman");
         globalCor=fisher_z_transformation(C);
         w=abs(globalCor)/sum(abs(globalCor));

         prob=zeros(num_objects,5);
         for i=1:num_objects
             for j=1:num_annotators
                 prob(i,:)=prob(i,:)+w(j)*to_dist(X(i,j));
             end
         end

         % fix the issue with degenerate distributions by replacing 1 by
         % 1-eps;
         % this simply avoid division by zero when computing the
         % reliabilities 
         for i=1:num_objects
             for j=1:5
                 prob(i,j)=min(1-eps,prob(i,j));
             end
         end

        % computing how surprizing is each rating given the estimated ground truth distribution
        P_hat=zeros(num_objects,num_annotators); 
        for ii=1:num_annotators
            for jj=1:num_objects
                 P_hat(jj,ii)=-log(prob(jj,X(jj,ii)));
            end
        end
        
        % Computing the reliability of each opinion scores
        Relaibility=P_hat.^(-1);


        % Compute the weights as function of the Relaibility
        w=zeros(num_objects,num_annotators);
        for i=1:num_objects
            w(i,:)= (Relaibility(i,:))/sum(Relaibility(i,:));
        end
        
end




%% this function converts on opinion score into a one hot vector 
function d=to_dist(vote)
         d=zeros(1,5);
         d(vote)=1;
end