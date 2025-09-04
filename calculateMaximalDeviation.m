
% ----
% Simulates the worst-case scenario by adding maximally biased attacker 
% ratings to a clean dataset and computing the resulting RMSE.
%
% Input: clean dataset, ground truth, number of attacker
% Output: RMSE and cleanDataset with attacker subjects added
% ----
function [rmse, data] = calculateMaximalDeviation(data,groundTruth,numAttacker)
    
    attackers = ones(numAttacker,size(data,2));
    
    for k = 1:size(data,2)
        if groundTruth(k) < 3
            attackers(:,k) = 5 .* attackers(:,k);
        end
    end
    
    data = [data ; attackers];
    attackerMOS = mean(data,1);
    
    squaredError = (groundTruth - attackerMOS).^2;
    mse = mean(squaredError);
    rmse = sqrt(mse);
end