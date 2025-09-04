
% ----
% Simulates a subjective rating matrix with I subjects and J items 
% based on the SUREAL model, incorporating subject bias and inconsistency. 
% ----
function [matrix,groundTruth] = simulation(I,J,seed)
    rng(seed)

    koniqData = load("brcw.mat");
    bias = koniqData.delta_brcw;
    inc = koniqData.v_brcw;
    psi = koniqData.psi_brcw;

    % get a bias and inconsistency for each subejct
    random_indices_subject = randperm(size(bias, 2), I);
    subjectBias = bias(random_indices_subject);
    subjectInc = inc(random_indices_subject);

    
    % get mean value of items
    random_indices_item = randperm(size(psi,1),J);
    itemPsi = psi(random_indices_item);
    groundTruth = itemPsi';


    % average bias = 0
    subjectBias = subjectBias - mean(subjectBias);


    % Initialize the result matrix
    matrix = zeros(I, J);

    % Generate the matrix row by row
    for i = 1:I
        for j = 1:J
            % Add mu (item mean), bias, and inconsistency noise 
            matrix(i, j) = itemPsi(j) + subjectBias(i) + normrnd(0, subjectInc(i));
        end
    end

    tau = [1.5 , 2.5 , 3.5 , 4.5, inf];

    matrix = quantize_matrix(matrix,tau);
end



function result = quantize_matrix(U, tau)
    % Converts tau to a 1x1xN array
    tau = reshape(tau, 1, 1, []);  
    
    % Compare U with each threshold in tau, results in a IxJxN logical array
    comparison = U < tau;  
    
    % get first threshold exceeded by each element in U
    [~, result] = max(comparison, [], 3); % (take only the indice)
end