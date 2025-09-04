
% ----
% Adds 'I_new' many random rows with values 1 - 5 (ACR) to the original matrix
%
% Input: AC-Rating matrix, amount of attacker, seed, flag
% Output: AC-Rating matrix with 'I_new' many random observer (matrix with or without origData)
% ----
function result = generateAttackSet(matrix, I_new, seed, flag)
    
    % set seed
    rng(seed,"philox"); % when set the algorithm to philox, the result of parfor is the same as when just using the for loop

    [~, numCols] = size(matrix);
    if strcmp(flag, "withoutOrigData")
        result = randi([1, 5], I_new, numCols);
    elseif strcmp(flag, "withOrigData")
        newRows = randi([1, 5], I_new, numCols);
        result = [matrix; newRows];
    else
        error('Invalid flag. Choose "withoutOrigData" or "withOrigData".');
    end
end