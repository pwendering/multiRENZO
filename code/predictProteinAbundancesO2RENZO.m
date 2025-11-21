function predictProteinAbundancesO2RENZO
% Prediction of protein abundances in different O2 conditions:
% * 2%
% * 21%
% * 40%
%
% The RENZO approach, which constrains ratios between predicted protein
% abundances, is extended to multiple conditions and applied to all three
% conditions jointly.
%
% Temperature: 21 °C (day)
% Relative humidity: 60%
% Ambient p(CO2): 400 µbar
% light intensity: 100 μmol photons m−2 s−1

top_dir = regexprep(mfilename('fullpath'), '(\\[^\\]+){2}$', '');
result_dir = fullfile(top_dir, 'results');

% Result workspace name
c = clock;
datestr = sprintf('%d%02.0f%02.0f', c(1:3));
result_ws = fullfile(result_dir, ['renzo_workspace_walker' '_' datestr]);

% Parameters for the Gurobi solver
solver_params = config('solver_params');

%% ecAraCore model compatible with temperature-dependent constraints
tmp = load(config('tgemFile'));
model = tmp.TGEM;
clear tmp

n_rxns = sum(~startsWith(model.rxns, {'draw_prot', 'prot_pool'}));
n_prot_model = sum(startsWith(model.rxns, 'draw_prot'));

%% Environment parameters
model.C_a = config('C_a');  % ambient p(CO2)
I = config('I');  % light intensity
T = config('T');  % temperature
O2 = 1e4*[21, 2, 40];  % ambient p(O2) (µbar)
n_o2 = numel(O2);

%% Get temperature-dependent optimization problems

% Initialize minimum relative growth rate (RGR)
min_rgr = 1;

for i=1:n_o2
    
    % set ambient p(O2)
    tmp_model = model;
    tmp_model.O_a = O2(i);
    
    % run simulation
    tmp_sol = simulateTempEffects(tmp_model, I, 'tempRange', T);
    
    % update minimum RGR across O2 conditions
    min_rgr = min(min_rgr, tmp_sol.mu_max);
    
end

%% Add constraint to add an upper bound on the sum of fluxes

qcps = cell(1, n_o2);

for i=1:n_o2
    
    tmp_model = model;
    tmp_model.O_a = O2(i);
    
    % Fix relative growth rate
    tmp_model.lb(tmp_model.c==1) = 0.5*min_rgr;
    tmp_model.ub(tmp_model.c==1) = 0.5*min_rgr;
    
    % Solve optimization problem to obtain a flux distribution for the
    % current condition
    [tmp_sol, ~, ~, tmp_qcp] = simulateTempEffects(tmp_model, I, 'tempRange', T);
    
    % Total sum of reaction fluxes (all reactions a irreversible -> no
    % negative fluxes)
    total_flux = sum(tmp_sol.x_step2(1:n_rxns));
    
    % Add constraint to add an upper bound on the sum of fluxes
    tmp_qcp.A = [tmp_qcp.A; [ones(1, n_rxns) zeros(1, n_prot_model+1) 0 0]];
    tmp_qcp.rhs = [tmp_qcp.rhs; 1.01*total_flux];
    tmp_qcp.sense = [tmp_qcp.sense; '<'];
    tmp_qcp.constrnames = [tmp_qcp.constrnames; 'UB_sum_flux'];
    
    % Reset lower and upper bounds (will be constrained later again)
    tmp_qcp.lb(ismember(tmp_qcp.varnames, 'GROWTH')) = 0;
    tmp_qcp.ub(ismember(tmp_qcp.varnames, 'GROWTH')) = 1000;
    
    qcps{i} = tmp_qcp;
    
end


%% Load experimental data

% Mean reporter intensities and standard deviations
load(fullfile(result_dir, 'measured_protein_abundances'), ...
    'mean_intensities', ...
    'sd_intensities', ...
    'protein_IDs')

% Proteins with non-zero measurements
nz_idx = all(mean_intensities>0, 2);

%% Test relationship of objective value and lower bound on RGR
% This is to plots the objective values of the multiRENZO solutions against
% a range of RGR values, which set as lower bounds for RGR in the
% respective optimization problem.

%{
objval_test = nan(11, 1);
for i = 1:11
    
    tmp_sol = multiRENZO(qcps, ...
        mean_intensities(nz_idx, :), ...
        protein_IDs(nz_idx), ...
        'solver_params', solver_params, ...
        'pct_rgr', i/10-0.1, ...
        'prot_sd', sd_intensities(nz_idx, :));
    
    if isequal(tmp_sol.status, 'OPTIMAL')
        objval_test(i) = tmp_sol.objval;
    end
end

figure
plot(0:10:100, objval_test, 'k', 'LineWidth', 1.3)
xlabel('Percentage of optimal RGR')
ylabel({'multiRENZO objective value', '(mmol gDW^{-1})'})
set(gca, 'FontSize', 14)
%}

%% Run RENZO across all three conditions
[renzo_sol, renzo_problem] = multiRENZO(qcps, ...
    mean_intensities(nz_idx, :), ...
    protein_IDs(nz_idx), ...
    'solver_params', solver_params, ...
    'pct_rgr', 0.5, ...  % 50% of the minimum RGR value across the three conditions
    'prot_sd', sd_intensities(nz_idx, :));

% Save intermediate results
save(result_ws, ...
    'renzo_sol', 'renzo_problem', 'qcps', 'model', 'O2', 'I', 'T')

%% Variability analysis of the predicted protein abundances

% Indices of the tolerance values for the deviations of protein abundance
% ratios between the conditions.
idx_theta = startsWith(renzo_problem.varnames, 'theta_');
n_theta = sum(idx_theta);

% 1% tolerance for the sum of deviations between measured and predicted protein ratios
error_sum_tol = 0.01;

% Minimum distance between measured and predicted ratios
min_dist = renzo_problem.obj'*renzo_sol.x;

% Create problem for variability analysis with an upper bound on the sum of
% theta values.
prot_var_problem = renzo_problem;
prot_var_problem.A = [renzo_problem.A; idx_theta'];
prot_var_problem.rhs = [renzo_problem.rhs; min_dist*(1+error_sum_tol)];
prot_var_problem.sense = [renzo_problem.sense; '<'];
prot_var_problem.constrnames = [renzo_problem.constrnames; {'MIN_SUM_THETA'}];

% Find protein draw reactions in the RENZO problem
idx_prot = startsWith(renzo_problem.varnames, 'draw_prot_');
n_prot = sum(idx_prot);
idx_prot_lin = find(idx_prot);

% Initialize minimum and maximum abundances for each 
min_abundances = zeros(n_prot, 1);
max_abundances = zeros(n_prot, 1);

% Create parallel pool (if this throws an error, either comment this line
% or delete the existing parallel pool)
parpool(6);

% Minimize and maximize each protein abundance
disp('Variability analysis of protein abundances')
parfor i=1:n_prot
    
    tmp_problem = prot_var_problem;
    tmp_problem.obj(:) = 0;
    tmp_problem.obj(idx_prot_lin(i)) = 1;
    
    tmp_problem.modelsense = 'min';
    sol_tmp = gurobi(tmp_problem, solver_params);
    min_abundances(i) = sol_tmp.x(idx_prot_lin(i));
    
    tmp_problem.modelsense = 'max';
    sol_tmp = gurobi(tmp_problem, solver_params);
    max_abundances(i) = sol_tmp.x(idx_prot_lin(i));
    
end

% Append results of variability analysis to results
save(result_ws, ...
    'min_abundances', 'max_abundances', 'min_dist', '-append')

%% Sampling of protein abundances at minimal error

% Create a matrix where rows are proteins and columns are O2 conditions
min_abundances = reshape(min_abundances, [], n_o2);
max_abundances = reshape(max_abundances, [], n_o2);

% Number of samples
n_samples = 1e4;
max_tries = 5*n_samples;

% Percent of ratios to select for minimization
p_select = 0.01;

% Initialize samples
samples = cell(1, n_o2);
for i = 1:n_o2
    samples{i} = sparse(size(prot_var_problem.A, 2)-n_theta, n_samples);
end

% Initialize number of failed sampling problems and coverage values
n_failed = zeros(n_o2, 1);
coverages = nan(n_o2, 1);

disp('Sampling protein abundances')
for i = 1:n_o2
    
    t1 = tic;
    
    % Proteins of condition-specific model for current O2 concentration
    tmp_prot_idx = idx_prot_lin((i-1)*n_prot/n_o2+1:i*n_prot/n_o2);
    tmp_n_prot = n_prot/n_o2;
    
    % Maximum values of differences between random input ratios and predicted
    % ratios
    max_delta = 1000*ones(tmp_n_prot, 1);
    
    % Create sampling problem
    sampling_problem = prot_var_problem;
    
    % Add matrix to QCP for minimization of absolute distances to random
    % vector of protein ratios
    delta_mat = zeros(tmp_n_prot, size(sampling_problem.A, 2) + 2*tmp_n_prot);
    delta_mat(:, tmp_prot_idx) = -speye(tmp_n_prot);
    delta_mat(:, size(sampling_problem.A, 2)+1:end) = [speye(tmp_n_prot) -speye(tmp_n_prot)];
    
    sampling_problem.A = [ ...
        sampling_problem.A sparse(size(sampling_problem.A, 1), 2*tmp_n_prot); ...
        sparse(delta_mat) ...
        ];
    sampling_problem.lb = [sampling_problem.lb; zeros(2*tmp_n_prot, 1)];
    sampling_problem.ub = [sampling_problem.ub; max_delta; max_delta];
    sampling_problem.vtype = [sampling_problem.vtype; repmat('C', 2*tmp_n_prot, 1)];
    sampling_problem.varnames = [ ...
        sampling_problem.varnames; ...
        strcat('DELTA_PLUS_', sampling_problem.varnames(tmp_prot_idx));...
        strcat('DELTA_MINUS_', sampling_problem.varnames(tmp_prot_idx)) ...
        ];
    sampling_problem.rhs = [sampling_problem.rhs; zeros(tmp_n_prot, 1)];
    sampling_problem.sense = [sampling_problem.sense; repmat('=', tmp_n_prot, 1)];
    sampling_problem.constrnames = [ ...
        sampling_problem.constrnames; ...
        strcat('delta_constr_', sampling_problem.varnames(tmp_prot_idx))];
    
    % Initialize weights for first norm minimization
    weights = zeros(tmp_n_prot, 1);
    sampling_problem.obj = [ ...
        zeros(size(sampling_problem.A, 2), 1); ...
        weights; weights];
    
    % Set model sense to minimization
    sampling_problem.modelsense = 'min';
    
    % Update fields for quadratic constraints
    n_quad = numel(sampling_problem.quadcon);
    for j = 1:n_quad
        sampling_problem.quadcon(j).Qc = [ ...
            sampling_problem.quadcon(j).Qc ...
            zeros(size(sampling_problem.quadcon(j).Qc, 1), 2*tmp_n_prot); ...
            zeros(2*tmp_n_prot, size(sampling_problem.quadcon(j).Qc, 2)+2*tmp_n_prot)];
        sampling_problem.quadcon(j).q = [ ...
            sampling_problem.quadcon(j).q; ...
            zeros(2*tmp_n_prot, 1)];
    end
    
    % Generate as many random orderings as necessary
    all_prot_idx = 1:tmp_n_prot;
    valid_prot_idx = all_prot_idx;
    n_select = max(1, ceil(p_select*numel(valid_prot_idx)));
    rand_idx = cell2mat(arrayfun(@(i)...
        valid_prot_idx(randperm(numel(valid_prot_idx), numel(valid_prot_idx))),...
        1:ceil(n_select*n_samples/numel(valid_prot_idx)), 'un', 0));
    
    % Solve minimization problems for sampling
    tic
    j = 1;
    tries = 0;
    while j <= n_samples && tries < max_tries
        
        tries = tries + 1;
        
        % Create random flux vector
        rand_abundances = (max_abundances(:, i) - min_abundances(:, i)) ...
            .* rand(tmp_n_prot, 1) + min_abundances(:, i);
        
        % Set weights in the optimization objective
        weights = zeros(tmp_n_prot, 1);
        curr_rand_idx = rand_idx((j-1)*n_select+1:j*n_select);
        weights(curr_rand_idx) = 1./max_abundances(curr_rand_idx, i);
        weights(isnan(weights)|isinf(weights)) = 0;
        sampling_problem.obj = [ ...
            zeros(size(sampling_problem.A, 2)-2*tmp_n_prot, 1); ...
            weights; weights];
        
        % Add random abundances to minimization problem (RHS)
        sampling_problem.rhs(startsWith(sampling_problem.constrnames, ...
            'delta_constr_')) = -rand_abundances;
        
        % Solve sampling problem
        min_sol = gurobi(sampling_problem, solver_params);
        
        % Check if the solver has found an optimal solution
        if isequal(min_sol.status, 'OPTIMAL')
            % Save predicted values, except variables required for the
            % sampling objective and deviations between measured and
            % predicted protein abundance ratios
            samples{i}(:, j) = min_sol.x(~startsWith(sampling_problem.varnames, ...
                {'DELTA_', 'theta_'}));
            j = j + 1;
        else
            n_failed = n_failed + 1;
        end
        
        if j > 1 && mod(j, 100) == 1
            fprintf('Done with %d samples (%d%%)', j-1, round(100*(j-1)/n_samples))
            t2 = toc;
            fprintf(' (appr. %.f min remaining)\n',...
                (t2-t1)*(n_samples-j+1)/(j-1)/60)
        end
        
    end
    
    % Calculate coverage
    coverages(i) = getSamplingCoverage(samples{i}(tmp_prot_idx, :),...
        min_abundances(:, i), max_abundances(:, i));
    
    % Append sampling results to result workspace
    save(result_ws, ...
        'samples', 'n_failed', 'min_abundances', 'max_abundances', ...
        'min_dist', 'O2', 'coverages', ...
        '-append');
end

tend=toc;
fprintf('Total sampling time: %.2f min\n', tend/60)

end