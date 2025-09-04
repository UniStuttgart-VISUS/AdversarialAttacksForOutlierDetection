

% Implementation of the LSA algorithm, as described in "An optimization
% model for outlier detection in categorical data" by Zengyou He et al.

function [inliers, bestResult] = perform_HB(matrix,k,seed)

    rng(seed);

    [numSubjects, numItems] = size(matrix);
    uniqueRatings = 1:5; % (ACR)

    % random initialization of k outliers and inliers
    flags = zeros(numSubjects,1);
    permutation = randperm(numSubjects);
    flags(permutation(1:k),1) = 1;

    % create matrix with number of each ACR for each item (needed for entropy calculation)
    uniqueRatingAmounts = zeros(numel(uniqueRatings),numItems);
    for k = 1:numel(uniqueRatings)
        uniqueRatingAmounts(k,:) = sum(matrix == uniqueRatings(k),1);
    end

    % save best outlier selection
    current_best_flags = flags;
    current_best_entropy = entropyCalculation(matrix, uniqueRatingAmounts, flags);
    improved = true;

    % iteration phase: 
    % iteratively try out all records as outlier (change it once with each outlier)
    % save the best outlier selection. 
    % If one run didnt improve the result --> exit
    while improved
        improved = false;

        outlierIndexes = find(flags == 1); 
        valid_indices = setdiff(1:numSubjects, outlierIndexes); 
        valid_indices = valid_indices(randperm(length(valid_indices))); % Shuffle


        for k = 1:numel(outlierIndexes)
            for c = 1:numel(valid_indices)
                flags(outlierIndexes(k)) = 0;   % exclude outlier k
                flags(valid_indices(c)) = 1;

                result = entropyCalculation(matrix,uniqueRatingAmounts,flags);

                if result < current_best_entropy
                    current_best_entropy = result;
                    current_best_flags = flags;
                    temp = outlierIndexes(k);
                    outlierIndexes(k) = valid_indices(c); 
                    valid_indices(c) = temp;
                    improved = true;
                else
                    flags(outlierIndexes(k)) = 1;   % include outlier k
                    flags(valid_indices(c)) = 0;
                end
                
            end
        end
        
    end

    inliers = matrix(current_best_flags==0,:);
    bestResult = current_best_entropy;

end