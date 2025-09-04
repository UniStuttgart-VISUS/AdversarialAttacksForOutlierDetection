

% added weights to output
function [quality, bias, inc, quality_CI, bias_CI, incCI, weights] = perform_SUREAL(Votes,threshold)

% This function run the sureal software 

% Votes is the matrix of observers rates
% the authors suggest to set the thresold=10^(-8)

  stop=0;
  quality=mean(Votes,2);
  bias=(mean((Votes-quality*ones(1,size(Votes,2)))))';


  while ~stop
       quality_prev=quality;
       residue = Votes - quality - bias';
    
       inc=(std(residue))';
       inc(inc == 0) = eps;


       weights=(1./(inc.^2))/sum((1./(inc.^2)));
       quality=(Votes-ones(size(Votes,1),1)*bias')*weights;

       bias=(mean((Votes-quality*ones(1,size(Votes,2)))))';
       stop=norm(quality_prev-quality)<threshold;
  end
  
   quality_CI=quality+[-ones(length(quality),1),ones(length(quality),1)]*1.96*( 1/sqrt( sum(inc.^(-2)) ) );
   bias_CI=bias+[-ones(length(bias),1),ones(length(bias),1)].*[inc, inc]*1.96/sqrt(size(Votes,1));
   k=size(Votes,1);
   incCI=[sqrt(k/ chi2inv(0.975,k))*inc,  sqrt(k/ chi2inv(0.025,k))*inc];
   
end
