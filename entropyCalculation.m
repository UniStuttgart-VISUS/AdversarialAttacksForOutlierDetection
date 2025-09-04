

function entropyOfSet = entropyCalculation(matrix,ratings_items,flags)
    [numSubjects, numItems] = size(matrix);
    uniqueRatings = 1:5; % ACR ratings

    % get outlierMatrix from flags and collect how the outliers rated each
    % item
    outlier = matrix(flags == 1,:);
    outlierRatings = zeros(numel(uniqueRatings),numItems);
    for k = 1:numel(uniqueRatings)      
        outlierRatings(k,:) = sum(outlier == uniqueRatings(k),1);
    end

    % update ratings for entropy calculations (remove ratings from outliers)
    ratings_items = ratings_items - outlierRatings; % "updated hashtable" (He et.al. )

    % entropy calculation
    probability = ratings_items ./ (numSubjects - sum(flags));
    logProbability = log(probability);
    logProbability(isinf(logProbability)) = 0; % handling -inf ( log(0) ) by replacing with zero. Afterwards it gets either way multiplicated by zero.

    entropyItems = -1 .* sum(probability .* logProbability,1);
    entropyOfSet = sum(entropyItems);
 
end