% =====================================================================
% DE-MEDT: Doppler effect–mean Euclidean distance threshold metaheuristic
% GitHub Repository: https://github.com/K-BiabaniHamedani/IDE
%
% This MATLAB implementation accompanies the paper:
% "A physics-based metaheuristic algorithm based on doppler effect 
%  phenomenon and mean euclidian distance threshold"
% Published in *Periodica Polytechnica Civil Engineering*.
%
% Authors:
%   Ali Kaveh
%       Professor of Structural Engineering, School of Civil Engineering,
%       Iran University of Science and Technology
%       ORCID: 0000-0002-0697-3263
%       Email: alikaveh@iust.ac.ir
%
%   Seyed Milad Hosseini
%       Research Fellow, School of Civil Engineering,
%       Iran University of Science and Technology
%       ORCID: 0000-0002-6059-6227
%       Email: hosseini_milad@alumni.iust.ac.ir
%
%   Ataollah Zaerreza
%       Assistant Professor of civil engineering, Ayatollah Boroujerdi University
%       ORCID: 0000-0002-0648-4072
%       Email: ataollah_zaerreza@abru.ac.ir
% =====================================================================
% DE-MEDT
% Implementation of the Doppler effect–mean Euclidean distance threshold (DE-MEDT) metaheuristic.
%
% Inputs:
% nobs        - Number of observers (population size)
% VarMin      - Minimum value of design variables (vector)
% VarMax      - Maximum value of design variables (vector)
% CostFunction- Handle to the objective function
% Maxt        - Maximum number of iterations
% MaxNFEs     - Maximum number of objective function evaluations
% VarSize     - Size of design variable vector
% ndim        - Dimension of the search space
% nto         - Run index for display purposes
%
% Outputs:
% BestWeights   - Best weight values per iteration
% BestCosts     - Best objective function values per iteration
% BestSol       - Best solution structure
% BestPositions - Best positions per iteration

%% ---------- MAIN FUNCTION ----------
function [BestWeights, BestCosts, BestSol, BestPositions] = DE_MEDT(nobs, VarMin, VarMax, CostFunction, Maxt, MaxNFEs, VarSize, ndim, nto)

% initialize iteration and objective function evaluation counters
t = 0;    % iteration counter
NFEs = 0; % objective function evaluation counter

% initialize population and evaluate initial solutions
[pop, ~, NFEs] = Initialization(nobs, VarMin, VarMax, VarSize, CostFunction, NFEs);

% preallocate arrays to store best results per iteration
BestCosts = zeros(Maxt, 1);
BestWeights = zeros(Maxt, 1);
BestPositions = zeros(Maxt, ndim);

% auxiliary containers used by the algorithm
Positions = zeros(nobs, ndim); % current positions matrix
Distance = zeros(nobs, 1);     % distance to centroid
R = zeros(Maxt, 1);            % mean distance per iteration
Normal_R = zeros(Maxt, 1);     % normalized mean distance per iteration

% template for an empty individual
empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Weight = [];

%% ---------- MAIN OPTIMIZATION LOOP ----------
while NFEs < MaxNFEs
    t = t + 1;
    newpop = repmat(empty_individual, nobs, 1); % initialize new population
    
    % collect current positions into a matrix
    for o = 1:nobs
        Positions(o, :) = pop(o).Position;
    end
    
    % compute mean position vector
    Mean_pos = mean(Positions(1:nobs, :));
    
    % compute Euclidean distance of each individual from the mean position vector
    for i = 1:nobs
        Distance(i,1) = sqrt(sum((pop(i).Position - Mean_pos).^2));
    end
    
    % loop over all observers
    for o = 1:nobs
        %% ---------- SELECT A DETERMINATIVE INDIVIDUAL ----------
        if o ~= 1 && o ~= nobs
            Selected_better = randi([1, o - 1]); % determinative individual
        elseif o == 1
            Selected_better = 1;
        else % o == nOs
            Selected_better = randi([1, nobs - 1]);
        end
        
        %% ---------- STEP-SIZE COMPUTATION ----------
        V = pop(Selected_better).Position;
        Vo = pop(Selected_better).Position - pop(o).Position;
        Vs = pop(Selected_better).Position - pop(end).Position;
        Landa_OS = rand(VarSize); % random scaling vector
        
        % compute step size (avoid division by zero)
        Sum_VVs = V + Vs;
        Sum_VVs(Sum_VVs == 0) = eps;
        stepsize = Landa_OS .* ((V + Vo) ./ Sum_VVs) .* Vs;
        
        % position update 
        newpop(o).Position = pop(o).Position + stepsize;
        
        %% ---------- Mean EUCLIDEAN DISTANCE THRESHOLD (MEDT) MECHANISM ----------
        R(t) = mean(Distance);                      % compute scatter radius index
        Normal_R(t) = R(t) ./ max(VarMax - VarMin); % compute normalized scatter radius index
        
        % compute convergence index
        if Normal_R(t) < 0.1
            Covergence_index = 10 .* Normal_R(t);
        else
            Covergence_index = 1;
        end
        
        %% ---------- DIVERSITY-ADAPTIVE STOCHASTIC PERTURBATION ----------
        if rand < 0.5 * (1 - Covergence_index)
            N = randi([1, ndim], 1, 1);
            newpop(o).Position(N) = newpop(1).Position(N) .* unifrnd(-1, 1, size(N)) .* R(t);
        end
        
        %% ---------- BOUNDARY HANDLING ----------
        newpop(o).Position = max(newpop(o).Position, VarMin);
        newpop(o).Position = min(newpop(o).Position, VarMax);
        
        %% ---------- EVALUATE NEW SOLUTION ----------
        NFEs = NFEs + 1;
        [newpop(o).Weight, newpop(o).Cost, newpop(o).sol] = CostFunction(newpop(o).Position, NFEs);
        
        % stop if maximum objective function evaluations reached
        if NFEs >= MaxNFEs
            break;
        end
    end
    
    % merge new and current populations to form the next iteration’s population
    pop = [pop; newpop];
    pop = SortPopulation(pop);
    pop = pop(1:nobs);
    
    %% ---------- UPDATE BEST SOLUTION ----------
    BestSol = pop(1);
    BestWeights(t) = BestSol.Weight;
    BestCosts(t) = BestSol.Cost;
    BestPositions(t, :) = BestSol.Position;
    
    % display progress
    disp(['Run ' num2str(nto) ': Iteration ' num2str(t) ': NFEs ' num2str(NFEs) ...
        ': Best Weight = ' num2str(BestWeights(t)) ': Best Cost = ' num2str(BestCosts(t))]);
end

end

%% ------------------- INITIALIZATION FUNCTION ------------------------
function [pop, BestSol, NFEs] = Initialization(nOs, VarMin, VarMax, VarSize, CostFunction, NFEs)
% initialize population using uniform random sampling and evaluate solutions

empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Weight = [];

pop = repmat(empty_individual, nOs, 1);

% random population
for o = 1:nOs
    pop(o).Position = unifrnd(VarMin, VarMax, VarSize);
    NFEs = NFEs + 1;
    [pop(o).Weight, pop(o).Cost, pop(o).sol] = CostFunction(pop(o).Position, NFEs);
end

% record best initial solution
pop = SortPopulation(pop);
BestSol = pop(1);

end

%% ------------------- SORT POPULATION FUNCTION ----------------------
function [pop, SortOrder] = SortPopulation(pop)
% sort a population structure array by penalized objective function values in ascending order

Costs = [pop.Cost];
[~, SortOrder] = sort(Costs);
pop = pop(SortOrder);

end