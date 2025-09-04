

% select outlier methods:
% (available: "KB", "CB", "MAZ", "SUREAL", "ESQR", "ZREC" ,"NLL" ,"HB"
% ,"LPCC")

methods = ["KB", "CB", "MAZ", "SUREAL", "ESQR", "ZREC" ,"NLL" ,"HB" ,"LPCC"];



% Set test parameters:
iterations = 250;
numSubjects = 30;
numItems = 20;
numAttacker = 5;

% HB method parameter: number of subjects to be rejected as outliers
kManyOutlier = 5; 



% Initialize matrices to store evaluation metrics across iterations
rmseAttackedGroundTruth = zeros(iterations,numel(methods)); % RMSE (Root Mean Square Error)
rmsdAttackedGroundTruth = zeros(iterations,numel(methods)); % RMSD (Root Mean Square Difference)
rai = zeros(iterations,numel(methods)); % RAI (Remainig Attacker Influence)

rmseNotAttacked = zeros(iterations,1); % Error between MOS of clean dataset and groundTruth (baseline error)

% Number of subjects classified as inliers after attack and outlier removal
% (used for FPR & FNR & ACC)
inlierSizes = zeros(iterations,numel(methods));
% Number of attacking subjects that remain undetected (not removed)
% (used for FPR & FNR & ACC)
remainingAttackerNumber = zeros(iterations,numel(methods));

% Stores subject weights assigned by soft outlier methods
finalWeights = zeros(iterations,numel(methods),numSubjects+numAttacker);


% robustness test across methods for various datasets
for m = 1:numel(methods)
    for i = 1:iterations
        disp(methods(m) + "   " + i+"/"+iterations)

        % simulate clean dataset
        [data, groundTruth] = simulation(numSubjects,numItems,i);

        % calculate Error between simulation's MOS and ground truth
        MOS = mean(data,1);
        rmse_MOS = sqrt(mean((groundTruth - MOS).^2));
        rmseNotAttacked(i,1) = rmse_MOS;

        
        % apply outlier method on clean dataset, in order to compute the
        % RMSD
        if (methods(m) == "KB")
            inlier = perform_KB(data);
        elseif (methods(m) == "CB")
            preRanks = tiedrank(data'); % important solely in the GA to avoid redundant computation
            inlier = perform_CB(data,preRanks,0);
        elseif (methods(m) == "MAZ")
            inlier = perform_MAZ(data);
        elseif (methods(m) == "SUREAL")
            [quality, ~, ~, ~, ~, ~, weights] = perform_SUREAL(data',1e-8);    % input of soft methods are transposed, because the implementation from ESQR code expect a JxI matrix
            prevQuality = quality';
        elseif (methods(m) == "ESQR")
            [q,~,~,~,~] = perform_ESQR(data');
            prevQuality = q';
        elseif (methods(m) == "ZREC")
            [Q,~,~, ~] = perform_ZREC(data');
            prevQuality = Q';
        elseif (methods(m) == "NLL")
            inlier = perform_NLL(data);
        elseif (methods(m) == "LPCC")
            inlier = perform_LPCC(data,0.75);
        elseif (methods(m) == "HB")
            inlier = perform_HB(data,kManyOutlier,i);
        else
            error("no valid method")
        end

        if ~ismember(methods(m),["ZREC", "SUREAL", "ESQR"])
            prevQuality = mean(inlier,1);
        end
        
   
        % run GA and take best attacker set
        [finalGeneration, allErrors,finalSizes, finalInliers, subjectWeights] = geneticAlgorithm(data,150,methods(m),numAttacker,300,0.5,3,groundTruth,kManyOutlier);
        [maxRMSE, maxIndex] = max(allErrors(:,:,end));

        % save RMSE of best attacker
        rmseAttackedGroundTruth(i,m) = maxRMSE;

        % save number of inlier subjects after attack 
        inlierSizes(i,m) = finalSizes(maxIndex);

        % differentiate between soft and hard outlier method
        if ~ismember(methods(m),["ZREC", "SUREAL", "ESQR"])

            % save number of attacker subjects that remain after applied outlier method
            isOriginal = ismember(finalInliers{maxIndex},data,'rows');
            remainingAttackerNumber(i,m) = (finalSizes(maxIndex) - sum(isOriginal));
            
            % save RAI
            rai(i,m) = (finalSizes(maxIndex) - sum(isOriginal)) / finalSizes(maxIndex);

            % save RMSD
            MOS_afterAttack = mean(finalInliers{maxIndex});
            rmsdAttackedGroundTruth(i,m) = sqrt(mean((prevQuality - MOS_afterAttack).^2));
        else

            remainingAttackerNumber(i,m) = numAttacker; % due to method being soft

            % save the assigned weight of each subject
            weights = subjectWeights{maxIndex};
            finalWeights(i,m,:) = weights;

            % save RAI
            rai(i,m) = sum(weights(numSubjects+1:end));

            % save RMSD
            if (methods(m) == "SUREAL")
                [quality, ~, ~, ~, ~, ~, weights] = perform_SUREAL(finalGeneration(:,:,maxIndex)',1e-8);
                attackedQuality = quality';
            elseif (methods(m) == "ESQR")
                [q,~,~,~,~] = perform_ESQR(finalGeneration(:,:,maxIndex)');
                attackedQuality = q';
            elseif (methods(m) == "ZREC")
                [Q,~,~, ~] = perform_ZREC(finalGeneration(:,:,maxIndex)');
                attackedQuality = Q';
            end
            rmsdAttackedGroundTruth(i,m) = sqrt(mean((prevQuality - attackedQuality).^2));
        end
    end
end



% create result matrix:

RMSE_Results = zeros(numel(methods),1);
RMSD_Results = zeros(numel(methods),1);
FPR_Results = zeros(numel(methods),1);
FNR_Results = zeros(numel(methods),1);
ACC_Results = zeros(numel(methods),1);
RAI_Results = zeros(numel(methods),1);

softMethods = ["ESQR", "SUREAL", "ZREC"];


% compute FNR & FPR
cleanRemaining = inlierSizes - remainingAttackerNumber;
FP = numSubjects - cleanRemaining;
FN = remainingAttackerNumber;
TP = numAttacker - remainingAttackerNumber;
TN = inlierSizes - remainingAttackerNumber;

FPR = FP ./ numSubjects;
FNR = FN ./ numAttacker;
ACC = (TP + TN) ./ (numSubjects + numAttacker);


% create table
for m = 1:numel(methods)
    RMSE_Results(m,1) = mean(rmseAttackedGroundTruth(:,m));
    RMSD_Results(m,1) = mean(rmsdAttackedGroundTruth(:,m));
    RAI_Results(m,1) = mean(rai(:,m));

    if ismember(methods(m), softMethods)
        FPR_Results(m,1) = NaN;
        FNR_Results(m,1) = NaN;
        ACC_Results(m,1) = NaN;
    else
        FPR_Results(m,1) = mean(FPR(:,m));
        FNR_Results(m,1) = mean(FNR(:,m));
        ACC_Results(m,1) = mean(ACC(:,m));
    end
end

table = [RMSE_Results RMSD_Results FPR_Results FNR_Results ACC_Results RAI_Results];

% add noOpt case
noOpt = zeros(iterations,2);
methods = [methods "NoOpt"];
for i = 1:iterations
    [data, groundTruth] = simulation(numSubjects,numItems,i);
    [rmse, ~] = calculateMaximalDeviation(data,groundTruth,numAttacker);
    noOpt(i,1) = rmse;
    [rmse, ~] = calculateMaximalDeviation(data,mean(data,1),numAttacker);
    noOpt(i,2) = rmse;
end
noOptResults = [mean(noOpt(:,1)) mean(noOpt(:,2)) 0 1 numSubjects/(numAttacker+numSubjects) numAttacker/(numAttacker+numSubjects)];
table = [table ; noOptResults];
[~,sortIdx] = sort(table(:,1));

finalTable_sorted = table(sortIdx,:);
methods_sorted = methods(sortIdx);





%print Latex code:

disp(" ")
disp(" ")
disp("Latex Code:")
disp(" ")


headers = {'RMSE', 'RMSD', 'FPR', 'FNR', 'ACC', 'RAI'};

% Start LaTeX table
fprintf('\\begin{table}[ht]\n\\centering\n');
fprintf('\\begin{tabular}{l%s}\n', repmat('c', 1, size(finalTable_sorted,2)));
fprintf('\\hline\n');

% Header row
fprintf('Rank');
for h = 1:numel(headers)
    fprintf(' & %s', headers{h});
end
fprintf(' \\\\\n\\hline\n');

% Table content with rank in first column
for i = 1:size(finalTable_sorted,1)
    methodName = methods_sorted{i};
    fprintf('%d. %s', i, methodName);
    for j = 1:size(finalTable_sorted,2)
        val = finalTable_sorted(i,j);
        if isnan(val)
            fprintf(' & --');
        else
            fprintf(' & %.3f', round(val, 3));
        end
    end
    fprintf(' \\\\\n');
end

fprintf('\\hline\n\\end{tabular}\n');
fprintf('\\end{table}\n');




% create density plot and save as "densityPlot.pdf" (show or hide LPCC depending on the toggle)
LPCC_toggle = "hide"; % "show" or "hide"

width = 12;
height = 0.75 * width;
figure('Units', 'centimeters', 'Position', [0, 0, width, height]);

rmseAttackedGroundTruth_withNoOpt = [rmseAttackedGroundTruth noOpt(:,1)];



hold on; % Hold the figure to plot multiple lines
for col = sortIdx'
    if methods(col) == "LPCC" && LPCC_toggle == "hide"
        continue;
    end

    [f, x] = ksdensity(rmseAttackedGroundTruth_withNoOpt(:,col)', 'Function', 'pdf');

    if col == numel(methods)
        plot(x, f, '--', 'LineWidth', 2, 'DisplayName', methods(col)); % Dashed line for NoOpt case
    else
        plot(x, f, 'LineWidth', 2, 'DisplayName', methods(col)); % Solid line
    end
end

% Labels and title
xlabel('RMSE', 'FontSize', 10);
ylabel('Probability density', 'FontSize', 10);
axis tight; 
xlim([0.07 inf]);

grid on;
legend('FontSize', 7);
hold off;

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width, height]);
set(gcf, 'PaperPosition', [0, 0, width, height]);
set(gcf, 'InvertHardcopy', 'on');

% Save the figure
print(gcf, 'densityPlot.pdf', '-dpdf', '-r600'); % PDF