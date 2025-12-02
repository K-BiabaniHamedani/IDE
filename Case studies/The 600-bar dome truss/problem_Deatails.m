function [problem_details] = problem_Deatails(Name)

switch Name
    case 'Dome600Bar_Frequency'
        Problem_name = 'Dome600Bar_Frequency';
        ndim = 25;                  % Number of design variables
        lb = 0.0001.*ones(1, ndim); % Lower bound of design variables
        ub = 0.01.*ones(1, ndim);   % Upperbound of design variables
        MaxNFEs = 20000;            % Maximum number of objective function evaluations
        nobs = 30;                  % Population size
        Maxt = round(MaxNFEs/nobs); % Maximum number of Iterations
        nr = 20;                    % Number of independent runs
        FitFunc = @(x, NFEs) fobj_Dome600Bar_Frequency(x, NFEs, MaxNFEs); % Objective function
        
    case 'Dome1410Bar_Frequency'
        Problem_name = 'Dome1410Bar_Frequency';
        ndim = 47;                  % Number of design variables
        lb = 0.0001.*ones(1, ndim); % Lower bound of design variables
        ub = 0.01.*ones(1, ndim);   % Upperbound of design variables
        MaxNFEs = 20000;            % Maximum number of objective function evaluations
        nobs = 30;                  % Population size
        Maxt = round(MaxNFEs/nobs); % Maximum number of Iterations
        nr = 20;                    % Number of independent runs
        FitFunc = @(x, NFEs) fobj_Dome1410Bar_Frequency(x, NFEs, MaxNFEs); % Objective function
end

problem_details.Problem_name = Problem_name;
problem_details.dim = ndim;
problem_details.lb = lb;
problem_details.ub = ub;
problem_details.MaxIt = Maxt;
problem_details.MaxOFE = MaxNFEs;
problem_details.nPop = nobs;
problem_details.NIRs = nr;
problem_details.FitFunc = FitFunc;

end