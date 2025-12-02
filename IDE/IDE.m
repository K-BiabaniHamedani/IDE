% =====================================================================
% IDE: improved Doppler effect metaheuristic
% GitHub Repository: https://github.com/K-BiabaniHamedani/IDE
%
% This MATLAB implementation accompanies the manuscript:
% "IDE: improved Doppler effect metaheuristic with adaptive learning and
% opposition-based initialization for large-scale frequency-constrained
% dome-truss optimization"
% Submitted to *Computers & Structures*.
%
% Authors:
%   Ali Kaveh
%       Professor of Structural Engineering, School of Civil Engineering,
%       Iran University of Science and Technology
%       ORCID: 0000-0002-0697-3263
%       Email: alikaveh@iust.ac.ir
%
%   Kiarash Biabani Hamedani
%       Postdoctoral Researcher, School of Civil Engineering,
%       Iran University of Science and Technology
%       ORCID: 0000-0002-0697-3263
%       Email: kiarash_biabani@alumni.iust.ac.ir
%
%   Seyed Milad Hosseini
%       Research Fellow, School of Civil Engineering,
%       Iran University of Science and Technology
%       ORCID: 0000-0002-6059-6227
%       Email: hosseini_milad@alumni.iust.ac.ir
% =====================================================================
% IDE
% Implementation of the improved Doppler effect (IDE) metaheuristic.
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
function [BestWeights, BestCosts, BestSol, BestPositions] = IDE(nobs, VarMin, VarMax, CostFunction, Maxt, MaxNFEs, VarSize, ndim, nto)

% initialize iteration and objective function evaluation counters
t = 0;    % iteration counter
NFEs = 0; % objective function evaluation counter

% initialize population and evaluate initial solutions
[pop, BestSol, NFEs] = Initialization(nobs, VarMin, VarMax, VarSize, CostFunction, NFEs);

% preallocate arrays to store best results per iteration
BestCosts = zeros(Maxt, 1);
BestWeights = zeros(Maxt, 1);
BestPositions = zeros(Maxt, ndim);

% template for an empty individual
empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Weight = [];

% diversification step interval (every 10% of MaxNFEs)
step = max(1, round(0.1 * MaxNFEs));
next_trigger = step;

% historical memory parameters
H = 30;         % maximum number of memory entries
p_memory = 0.5; % probability to use memory
memory = [];    % memory stores successful r values and success counts [r, count]

%% ---------- MAIN OPTIMIZATION LOOP ----------
while NFEs < MaxNFEs
    t = t + 1;
    newpop = repmat(empty_individual, nobs, 1); % initialize new population
    
    % loop over all observers
    for o = 1:nobs
        %% ---------- SELECT INDIVIDUALS FOR STEP SIZE ----------
        if o ~= 1 && o ~= nobs
            Selected_better = randi([1, o - 1]); % determinative individual
            Selected_worse = randi([o + 1, nobs]); % worse individual
        elseif o == 1
            Selected_better = 1;
            Selected_worse = randi([2, nobs]);
        else % o == nOs
            Selected_better = randi([1, nobs - 1]);
            Selected_worse = nobs;
        end
        
        %% ---------- STEP-SIZE COMPUTATION ----------
        V = pop(Selected_better).Position;
        Vo = pop(Selected_better).Position - pop(o).Position;
        Vs = pop(Selected_better).Position - pop(Selected_worse).Position;
        Landa_OS = rand(VarSize); % random scaling vector
        
        % avoid division by zero
        Sum_VVs = V + Vs;
        Sum_VVs(Sum_VVs == 0) = eps;
        
        % randomized amplification of step size
        if rand < 0.5
            RM = 1 + rand(1, ndim); % amplification factor
            stepsize = Landa_OS .* ((V + Vo) ./ Sum_VVs) .* Vs .* RM;
        else
            stepsize = Landa_OS .* ((V + Vo) ./ Sum_VVs) .* Vs;
        end
        
        %% ---------- DETERMINE PARAMETER r (HISTORICAL MEMORY) ----------
        use_memory = ~isempty(memory) && (rand < p_memory); % decide whether to use memory
        used_memory_idx = 0; % index of memory entry used (if any)
        
        if use_memory
            m = size(memory, 1);
            if m == 1
                idx = 1;
            else
                counts = memory(:, 2);
                probs = counts / sum(counts);
                cum = cumsum(probs);
                rand_val = rand;
                idx = find(cum >= rand_val, 1, 'first'); % roulette-wheel selection
            end
            r = memory(idx, 1); % selected historical r
            used_memory_idx = idx;
        else
            r = rand; % generate new random parameter
        end
        
        % compute adaptive weighting coefficient
        alpha = r * (1 - rand * NFEs / MaxNFEs);
        
        % select good solution for the update
        if o > 1
            upper = max(1, floor(o/2));
            GoodSol = pop(randi([1, upper]));
        else
            GoodSol = pop(1);
        end
        
        % position update (adaptive elite-guided)
        if NFEs <= MaxNFEs/2
            newpop(o).Position = alpha * pop(o).Position + (1 - alpha) * GoodSol.Position + stepsize;
        else
            newpop(o).Position = alpha * pop(o).Position + (1 - alpha) * BestSol.Position + stepsize;
        end
        
        %% ---------- ADAPTIVE BOUNCE-BACK BOUNDARY HANDLING ----------
        pos = newpop(o).Position < VarMin;
        newpop(o).Position(pos) = (pop(o).Position(pos) + VarMin(pos)) / 2;
        pos = newpop(o).Position > VarMax;
        newpop(o).Position(pos) = (pop(o).Position(pos) + VarMax(pos)) / 2;
        
        %% ---------- EVALUATE NEW SOLUTION ----------
        NFEs = NFEs + 1;
        [newpop(o).Weight, newpop(o).Cost, newpop(o).sol] = CostFunction(newpop(o).Position, NFEs);
        
        %% ---------- MEMORY UPDATE ----------
        if newpop(o).Cost < pop(o).Cost
            % successful update: replace old solution
            pop(o) = newpop(o);
            
            if used_memory_idx > 0
                % memory was used: increment success count and refresh entry
                cnt = memory(used_memory_idx, 2) + 1;
                if used_memory_idx == size(memory, 1)
                    memory(used_memory_idx, 2) = cnt;
                else
                    entry = [memory(used_memory_idx, 1), cnt];
                    memory = [memory(1:used_memory_idx - 1, :); memory(used_memory_idx + 1:end, :); entry];
                end
            else
                % memory not used: insert new r or increment existing entry
                new_r = [r, 1]; % initial success_count = 1
                if isempty(memory)
                    memory = new_r;
                else
                    [isDup, loc] = ismember(new_r(1), memory(:, 1));
                    if isDup
                        memory(loc, 2) = memory(loc, 2) + 1;
                    else
                        if size(memory, 1) < H
                            memory = [memory; new_r];
                        else
                            memory = [memory(2:end, :); new_r]; % FIFO replacement
                        end
                    end
                end
            end
        else
            % failed update: possibly remove or decrement used memory entry
            if used_memory_idx > 0 && ~isempty(memory)
                if size(memory, 1) > 1
                    memory(used_memory_idx, :) = [];
                else
                    memory(1, 2) = max(1, memory(1, 2) - 1);
                end
            end
        end
        
        % stop if maximum objective function evaluations reached
        if NFEs >= MaxNFEs
            break;
        end
    end
    
    %% ---------- PERIODIC POPULATION DIVERSIFICATION ----------
    if NFEs >= next_trigger && NFEs < MaxNFEs
        num_reinit = max(1, round(0.1 * nobs)); % reinitialize 10% of population
        pop = SortPopulation(pop); % sort population to identify worst
        for i = 1:num_reinit
            idx = nobs - i + 1; % index of a worst individual
            pop(idx).Position = unifrnd(VarMin, VarMax, VarSize); % random new position
            NFEs = NFEs + 1;
            [pop(idx).Weight, pop(idx).Cost, pop(idx).sol] = CostFunction(pop(idx).Position, NFEs);
            if NFEs >= MaxNFEs
                break;
            end
        end
        next_trigger = next_trigger + step;
    end
    
    %% ---------- UPDATE BEST SOLUTION ----------
    pop = SortPopulation(pop);
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
function [pop, BestSol, NFEs] = Initialization(nobs, VarMin, VarMax, VarSize, CostFunction, NFEs)
% initialize population using uniform random sampling, opposition-based
% learning, and small perturbations

empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Weight = [];

pop = repmat(empty_individual, nobs, 1);
popop = repmat(empty_individual, nobs, 1);
popqo = repmat(empty_individual, nobs, 1);
popsqo = repmat(empty_individual, nobs, 1);
popso = repmat(empty_individual, nobs, 1);

% random population
for i = 1:nobs
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    NFEs = NFEs + 1;
    [pop(i).Weight, pop(i).Cost, pop(i).sol] = CostFunction(pop(i).Position, NFEs);
end

% opposite population
for i = 1:nobs
    popop(i).Position = VarMax + VarMin - pop(i).Position; % opposite point
    NFEs = NFEs + 1;
    [popop(i).Weight, popop(i).Cost, popop(i).sol] = CostFunction(popop(i).Position, NFEs);
end

MID = (VarMax + VarMin) / 2; % midpoint of variable bounds

% quasi-opposite population
for i = 1:nobs
    popqo(i).Position = MID + (MID - pop(i).Position) .* rand(VarSize);
    NFEs = NFEs + 1;
    [popqo(i).Weight, popqo(i).Cost, popqo(i).sol] = CostFunction(popqo(i).Position, NFEs);
end

% secondary quasi-opposition population
for i = 1:nobs
    popsqo(i).Position = MID + (MID - popop(i).Position) .* rand(VarSize);
    NFEs = NFEs + 1;
    [popsqo(i).Weight, popsqo(i).Cost, popsqo(i).sol] = CostFunction(popsqo(i).Position, NFEs);
end

% shifted opposite population
for i = 1:nobs
    for j = 1:length(VarMin)
        if popop(i).Position(j) > MID(j)
            popso(i).Position(j) = popop(i).Position(j) + (VarMax(j) - popop(i).Position(j)) * rand;
        else
            popso(i).Position(j) = VarMin(j) + (popop(i).Position(j) - VarMin(j)) * rand;
        end
    end
    NFEs = NFEs + 1;
    [popso(i).Weight, popso(i).Cost, popso(i).sol] = CostFunction(popso(i).Position, NFEs);
end

% combine all subpopulations and keep top nobs solutions
pop = [pop; popqo; popsqo; popso; popop];
pop = SortPopulation(pop);
pop = pop(1:nobs);
BestSol = pop(1);

end

%% ------------------- SORT POPULATION FUNCTION ----------------------
function [pop, SortOrder] = SortPopulation(pop)
% sort a population structure array by penalized objective function values in ascending order

Costs = [pop.Cost];
[~, SortOrder] = sort(Costs);
pop = pop(SortOrder);

end