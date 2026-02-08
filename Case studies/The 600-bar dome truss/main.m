clc;
clear;
close all;
format short e
rng('shuffle')

%% Problem Definition
Name = 'Dome600Bar_Frequency';

%% Information of the case studies
[problem_details] = problem_Deatails(Name);
ndim = problem_details.dim;        % Number of design ariables
lb = problem_details.lb;           % Lower bound of design variables
ub = problem_details.ub;           % Upper bound of design variables
Maxt = problem_details.MaxIt;      % Maximum number of iterations
MaxNFEs = problem_details.MaxOFE;  % Maximum number of objective function evaluations
nobs = problem_details.nPop;       % Population size
nr = problem_details.NIRs;         % Number of independent runs
FitFunc = problem_details.FitFunc; % Objective function
VarSize = [1 ndim];                % Size of design variable vector

%% Preallocate arrays for final results
numRows = Maxt;                                  % Determine the size of convergence data
Costs = zeros(numRows, nr);                      % Preallocate for Costs
Weights = zeros(numRows, nr);                    % Preallocate for Weights
Best_Costs = zeros(nr, 1);                       % Preallocate for Best_Costs
Best_Weights = zeros(nr, 1);                     % Preallocate for Best_Weights
time = zeros(nr, 1);                             % Preallocate for times
Best_Positions = struct('Bestpos', cell(1, nr)); % Preallocate for Best_Positions

%% Get optimization results
bestcostall = struct();

nto = 0;
runs = 0;

while nto < nr 
    tStart = tic;
    nto = nto + 1;
    runs = runs + 1;
    %[gbestweight, gbestfitness, gbestsol, gbestX] = DE_MEDT(nobs, lb, ub, FitFunc, Maxt, MaxNFEs, VarSize, ndim, nto);
    [gbestweight, gbestfitness, gbestsol, gbestX] = IDE(nobs, lb, ub, FitFunc, Maxt, MaxNFEs, VarSize, ndim, nto);
    bestcostall(nto).weight = gbestweight;
    bestcostall(nto).best = gbestfitness;
    bestcostall(nto).option.Best_score = gbestsol.Cost;
    bestcostall(nto).option.Best_weight = gbestsol.Weight;
    bestcostall(nto).option.Best_pos = gbestX;
    bestcostall(nto).Bestsol = gbestsol;
    Costs(:, nto) = bestcostall(nto).best;
    Weights(:, nto) = bestcostall(nto).weight;
    Best_Costs(nto, 1) = bestcostall(nto).option.Best_score;
    Best_Weights(nto, 1) = bestcostall(nto).option.Best_weight;
    Best_Positions(nto).Bestpos = bestcostall(nto).option.Best_pos(end, :);
    time(nto, 1) = toc(tStart);
    
    if Best_Weights(nto, 1) ~= Best_Costs(nto, 1)
        nto = nto - 1;
    end
    
end

%% Final Results
[FinalResults.Best, m] = min(Best_Costs);                        % Best value and its run number
FinalResults.mean = mean(Best_Costs);                            % Mean of all best scores
FinalResults.std = std(Best_Costs);                              % Standard deviation of best scores
FinalResults.worst = max(Best_Costs);                            % Worst value
FinalResults.Best_RunNumber = m;                                 % Run number for the best score
FinalResults.Bestposition_BestRun = Best_Positions(m).Bestpos;   % Best position
FinalResults.Convergence_History_BestRun = bestcostall(m).best;  % Best run history
FinalResults.Convergence_History_Mean = mean(Costs, 2);          % Mean convergence history
FinalResults.BestSol_all = bestcostall(m).Bestsol;               % Best solution of all runs