
% ----
% Input: AC-Rating matrix
% Output: matrix with Inliers
% ----
function [newMatrix] = perform_MAZ(matrix)

    mu = mean(matrix,1);
    sigma = std(matrix);

    % Prevent division by zero
    sigma(sigma == 0) = eps;

    abs_z_scores = abs(matrix - mu) ./ sigma;  % I x J matrix
    criterion = mean(abs_z_scores, 2);  % I x 1 vector (criterion per subject)

    removalIndices = find(criterion > 1);

    newMatrix = matrix;
    newMatrix(removalIndices,:) = [];

end
