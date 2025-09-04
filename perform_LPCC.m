
% ----
% calculates Inliers according to ITU-T P.910
% Input: AC-Rating matrix, threshold 
% Output: matrix with Inliers
% ----
function matrix = perform_LPCC(matrix,threshold)

    while true
        mu = mean(matrix,1);
        pearsonCorr = corr(mu', matrix', 'Type', 'Pearson')';

        [minValue, minValueIndex] = min(pearsonCorr);
        if minValue < threshold
            matrix(minValueIndex,:) = [];
        else
           break;
        end
    end
end