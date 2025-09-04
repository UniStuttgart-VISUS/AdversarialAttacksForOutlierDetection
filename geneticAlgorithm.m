
% Genetic Algortihm with the operators: Selection, Crossover, Mutation &
% Elitism. Improves attacks by maximizing RMSE.

function [finalGeneration, allErrors, finalSizes, finalInlier, finalWeights] = geneticAlgorithm(data,generationSize, outlierDetectionMethod, numAttacker,maxIteration,percentMutations,percentElitism,groundTruth,kManyOutlier)

    % generationSize has to be even, so the size of the new generation
    % is the same as the old one (2 parents create 2 children)
    if mod(generationSize, 2) == 1  
        error("generationSize needs to be even");
    end

    % number of chromosomes that are part of Elitism
    amountElite = ceil(generationSize * (percentElitism/100));
    % preRanks needed for more efficient computation of CB
    preRanks = tiedrank(data');


    [subjects, items] = size(data);
    allErrors = zeros(generationSize,1,maxIteration);

    % set up the first generation of attackers randomly
    currentGeneration = zeros(subjects+numAttacker, items , generationSize);
    for k= 1:generationSize
        currentGeneration(:,:,k) = generateAttackSet(data,numAttacker,k,"withOrigData");
    end

    
    % calculate the fitness of the initial generation:
    [fitnessValues, ~,~,~]= fitness(groundTruth,currentGeneration,outlierDetectionMethod,preRanks,numAttacker,kManyOutlier);

    for k = 1:maxIteration
        % Temporarily save the best chromosomes from last generation (elitism)
        [eliteFitnessValue, eliteIdx] = maxk(fitnessValues, amountElite);
        eliteChromosomes = currentGeneration(:,:,eliteIdx);

        % use previous generation for GA procedure (selection --> crossover --> mutation)
        selectedParents = selection(currentGeneration,fitnessValues,k);    % passing the iteration, so that the result is reproducible (seed)
        createdChildren = crossoverBoth(selectedParents,k,data);
        currentGeneration = mutation(createdChildren,k,numAttacker,percentMutations);

        % calculate the fitness for the new generation
        [fitnessValues, newErrors, ~, ~] = fitness(groundTruth,currentGeneration,outlierDetectionMethod,preRanks,numAttacker,kManyOutlier);
        
        % replace the elite ones with the worst performing chromosomes and
        % also replace the old fitnessValues with the ones from Elite
        [~, worstIdx] = mink(fitnessValues,amountElite);
        currentGeneration(:,:,worstIdx) = eliteChromosomes;
        fitnessValues(worstIdx) = eliteFitnessValue;
        newErrors(worstIdx,1) = eliteFitnessValue; % this step of replacing the error of worst performing chromosomes with the fitnessValue, is only valid, since error and fitness is the rmse!!

        allErrors(:,1,k) = newErrors;
    end
    finalGeneration = currentGeneration;
    [~,~,finalSizes, finalInlier, finalWeights] = fitness(groundTruth,finalGeneration,outlierDetectionMethod,preRanks,numAttacker,kManyOutlier);
end


% ----
% Evaluates each chromosome (attack) by applying the given outlier method 
% and computing the RMSE to serve as its fitness value.
% ----
function [fitnessResult, newErrors, newSizes, finalInlier, finalWeights] = fitness(groundTruth,generation,outlierMethod,preRanks,numAttacker,kManyOutlier)
    generationSize = size(generation,3);
    items = size(generation(:,:,1),2);
    % preAllocate for results (fitnessResult holds fitnessValues for each attack chromosome)
    fitnessResult = zeros(1,generationSize);
    % preAllocate for the results of error Values
    newErrors = zeros(generationSize,1);
    % preAllocate for size of inlier
    newSizes = zeros(generationSize,1);
    
    finalInlier = cell(1,generationSize);
    finalWeights = cell(1,generationSize);
    for k = 1:generationSize
        if (outlierMethod == "KB")
            results = perform_KB(generation(:,:,k));
        elseif (outlierMethod == "CB")
            results = perform_CB(generation(:,:,k),preRanks,numAttacker);
        elseif (outlierMethod == "MAZ")
            results = perform_MAZ(generation(:,:,k));
        elseif (outlierMethod == "CBS")
            results = perform_CBS(generation(:,:,k),preRanks,numAttacker);
        elseif (outlierMethod == "HB")
            results = perform_HB(generation(:,:,k),kManyOutlier,k);
        elseif (outlierMethod == "SUREAL")
            [quality, ~, ~, ~, ~, ~, weights] = perform_SUREAL(generation(:,:,k)',1e-8);    % input of soft methods are transposed, because the implementation from ESQR code expect a JxI matrix
            results = quality';
        elseif (outlierMethod == "ESQR")
            [q,~,~,weights,~] = perform_ESQR(generation(:,:,k)');
            weights = sum(weights',2)/items;  % since each score is weighted, not just subjects
            results = q';
        elseif (outlierMethod == "ZREC")
            [Q,~,~, weights] = perform_ZREC(generation(:,:,k)');
            weights = weights';
            results = Q';
        elseif (outlierMethod == "NLL")
            results = perform_NLL(generation(:,:,k));
        elseif (outlierMethod == "LPCC")
            results = perform_LPCC(generation(:,:,k),0.75);
        else
            error("no valid method")
        end



        if ~ismember(outlierMethod,["ZREC", "SUREAL", "ESQR"])
            finalInlier{k} = results; 
            finalWeights{k} = [];
            MOS_inliers = mean(results, 1); % Compute mean for inliers
        else
            finalInlier{k} = []; 
            finalWeights{k} = weights;
            MOS_inliers = results; % Use direct MOS
        end

        % calculate RMSE
        squaredError = (groundTruth - MOS_inliers).^2;
        mse = mean(squaredError);
        rmse = sqrt(mse);

        % save results for record
        fitnessResult(k) = rmse;
        newErrors(k,1) = rmse;
        newSizes(k,1) = size(results,1);
    end
end



% ----
% select parents according to their fitnessValue (probability proportional
% to the fitnessValue).
% Function takes current generation, fitnessValues and current iteration
% ----
function result = selection(generation,fitnessValues,iteration)
    rng(iteration);
    
    % implementation of the roulette wheel parent selection
    generationSize = size(generation,3);
    prob = fitnessValues ./ sum(fitnessValues); % assign probability 

    selectedIndex = zeros(1,generationSize);

    % each chromosome has its "interval" (proportional to fitness)
    % depending on the random number, chromosomes are chosen if random
    % number lies in this interval
    for k = 1:generationSize
        selectedIndex(k) = 1 + sum(rand > cumsum(prob));
    end
    
    % keep only the parents that were selected
    result = generation(:,:,selectedIndex);
end


% ----
% function takes offspring and mutates numberofMutations
% many ratings, by assigning a rating value, that differs from the previous
% one.
% ----
function result = mutation(children,iteration,numAttacker,percentMutations)
   
    rng(iteration);

    [amountSubjects,amountItems,amountChromosomes] = size(children);
    
    % small percentage of all ratings get a mutation (only mutations in the attacker set)
    numberOfMutations = ceil(numAttacker*amountItems*amountChromosomes * (percentMutations/100));

    % create numberOfMutations many random coordinates
    mutationIndicesRow = randi([amountSubjects-numAttacker+1, amountSubjects], 1, numberOfMutations);
    mutationIndicesItem = randi([1, amountItems], 1, numberOfMutations);     
    mutationIndicesChromosom = randi([1, amountChromosomes], 1, numberOfMutations);

    
    % for the coordinates, assign a random integer between 1 and 5
    for k = 1:numberOfMutations
        % get new rating value, that differs from old one
        currentRating = children(mutationIndicesRow(k),mutationIndicesItem(k),mutationIndicesChromosom(k));
        possible_values = setdiff(1:5, currentRating);
        randomInteger = randsample(possible_values,1);

        children(mutationIndicesRow(k),mutationIndicesItem(k),mutationIndicesChromosom(k)) = randomInteger;
    end
    result = children;
end




% ----
% Applies crossover by swapping both items (columns) and attacker subjects (rows)
% between matched parent chromosomes to generate diverse offspring.
% ----
function result = crossoverBoth(parents,iteration,origData)
    rng(iteration);
    rows = size(origData,1);
    firstAttackerObserver = rows +1;
    lastAttackerObserver = size(parents,1);
    
    % create indexes, to then build matches that get crossed
    indexes = 1:size(parents,3);
    % randomize the indexes
    permutationIndices = randperm(length(indexes));
    permutedArray = indexes(permutationIndices);
    
    % create matches as (n x 2) matrix with the indexes of parents that together will do crossover
    matches = reshape(permutedArray, 2, []).';
    % array of the columns & rows that get crossed 
    amountCrossesPerMatchColumns = randi([1, size(parents,2)-1],1,size(matches,1)); % exclude the case to exchange all subjects/items (meaningless)
    amountCrossesPerMatchRows = randi([1, lastAttackerObserver-firstAttackerObserver],1,size(matches,1)); 
    
    % get random indexes for each element
    randomIndexesColumns = arrayfun(@(num) randperm(size(parents,2), num), amountCrossesPerMatchColumns, 'UniformOutput', false);
    randomIndexesRows = arrayfun(@(num) firstAttackerObserver - 1 + randperm(size(parents,1) - firstAttackerObserver + 1, num), amountCrossesPerMatchRows, 'UniformOutput', false);

    children = parents;
    for k = 1:size(matches,1)
        % swap the columns between the match 
        children(:,randomIndexesColumns{k},matches(k,1)) = parents(:,randomIndexesColumns{k},matches(k,2));
        children(:,randomIndexesColumns{k},matches(k,2)) = parents(:,randomIndexesColumns{k},matches(k,1));
        % swap also the rows between the match
        temp = children(randomIndexesRows{k},:,matches(k,1));
        children(randomIndexesRows{k},:,matches(k,1)) = children(randomIndexesRows{k},:,matches(k,2));
        children(randomIndexesRows{k},:,matches(k,2)) = temp;
    end
    
    result = children;
end