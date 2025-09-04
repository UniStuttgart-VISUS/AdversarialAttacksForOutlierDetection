
% ----
% Input: AC-Rating matrix
% Output: matrix with Inliers
% ----
function inlier = perform_NLL(matrix)
    run = true;
    currentInlier = matrix;
    while run
        run = false;
        [I, J] = size(currentInlier);
        probDistribution = zeros(J,5);
        
        % Compute probability distribution for each column
        for value = 1:5
            probDistribution(:,value) = sum(currentInlier == value, 1) / I;
        end
    
        nll = zeros(I,J);
        for i = 1:I
            for j = 1:J
                nll(i,j) = -log(probDistribution(j,currentInlier(i,j)));
            end
        end
    
        % sum all nll and divide by amount items
        nllSubject = mean(nll,2);
    
        if any(nllSubject > 1.31)
            run = true;
            [~, sortingIndices] = sort(nllSubject, 'descend');
            currentInlier(sortingIndices(1),:) = [];
        end
    end
    inlier = currentInlier;
end