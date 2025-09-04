

% ----
% calculates Inliers according to BT.500
% Input: AC-Rating matrix
% Output: matrix with Inliers
% ----
function newMatrix = perform_CB(matrix,preRanks,numAttacker)

    mu = mean(matrix,1);
    mu_ranks = tiedrank(mu');

    attackerRanks = tiedrank(matrix(end-numAttacker+1:end,:)');
    allRanks = [preRanks, attackerRanks];


    pearsonCorr = corr(mu', matrix', 'Type', 'Pearson')';
    spearmanCorr = corr(mu_ranks,allRanks,"Type","Pearson")';


    % subjects with no variance (same rating for all items) result in NaN
    % value
    pearsonCorr(isnan(pearsonCorr)) = 0;
    spearmanCorr(isnan(spearmanCorr)) = 0;


    c = min(pearsonCorr, spearmanCorr);


    % calulate threshold
    mean_c = mean(c);
    sigma_c = std(c);
    c_min = min(mean_c - sigma_c, 0.7);

    % discard elements that are below threshold
    subjectsToKeep = c >= c_min;
    newMatrix = matrix(subjectsToKeep,:);

end
