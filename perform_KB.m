

% ----
% calculates Inliers according to BT.500
% Input: AC-Rating matrix
% Output: matrix with Inliers
% ----
function newMatrix = perform_KB(matrix)

    % beta_2 Test, column wise
    beta = m_x(matrix,4) ./ m_x(matrix,2).^2;
    isNormal = (beta >= 2) & (beta <= 4);
   
    items = size(matrix,2);

    mu = mean(matrix,1);
    sigma = std(matrix, 0,1);


    upperBound = zeros(size(mu));
    lowerBound = zeros(size(mu));


    % calculate upper and lower Bound for each item different,
    % depending on the kind of distribution
    upperBound(isNormal) = mu(isNormal) + 2 * sigma(isNormal);
    lowerBound(isNormal) = mu(isNormal) - 2 * sigma(isNormal);

    
    upperBound(~isNormal) = mu(~isNormal) + sqrt(20) * sigma(~isNormal);
    lowerBound(~isNormal) = mu(~isNormal) - sqrt(20) * sigma(~isNormal);
    
   
    % compare score of i to bounds
    P = matrix >= upperBound;
    Q = matrix <= lowerBound;


    % P_i stores how many elements subject i rated above the bound, Q_i how many below
    P_i = sum(P,2);
    Q_i = sum(Q,2);

    
    % NaN values in conditionTwo do not matter, since if P_i + Q_i = 0,
    % then conditionOne is 0 and therefore the subject is not going to be
    % discarded.
    conditionOne = (P_i + Q_i) / items;
    conditionTwo = abs((P_i - Q_i) ./ (P_i + Q_i));


    % get subjects that fulfill the conditions
    subjectsToDiscard = conditionOne > 0.05 & conditionTwo < 0.3;
    
    % keep all elements that didnt fulfill the condition
    newMatrix = matrix(~subjectsToDiscard,:);

end


function result = m_x(matrix,x)
    mu = mean(matrix,1);
    deviationSquared = (matrix - mu).^x;
    result = sum(deviationSquared) / size(matrix,1);
end
