function [solution, problem] = multiRENZO(varargin)
%MULTIRENZO Extension of the RENZO approach to multiple conditions
%
% This function constraints ratios between predicted protein abundances
% using ratios of measured abundances between the conditions.
% 
% Input
%   cell models:                cell array of enzyme-constrained metabolic models,
%                               configured for the respective conditions if applicable
%   double protein_abundances:  relative proteomics data for each condition
%   cellstr protein_ids:        protein UniProt IDs for each measured
%                               protein abundance
%   double rgr:                 relative growth rates (h^-1) for each
%                               condition
%   double bio_idx:             index of the biomass reaction (assumed to
%                               be the same in all models)
%                               If not provided, the current optimization
%                               objective is used.
%   double pct_rgr:             minimum percentage of RGR to enforce in the
%                               minimization of distances to measured
%                               ratios; default: 0.5
%   double prot_sd:             standard deviations of protein measurements;
%                               The maximum standard deviation of two conditions will
%                               be used as weight in the objective for the
%                               respective theta, similar to chi2.
%                               The standard deviations will be re-scaled, such that
%                               their mean corresponds to the mean of the
%                               predicted protein abundances from the first
%                               optimization problem.
%                               By default, all weights are equal to one.
%   struct solver_params:       parameters passed to the gurobi solver
% Output
%   struct solution:            Solution to the multiRENZO optimization,
%                               obtained from the Gurobi solver
%   struct problem:             Problem structure of the multiRENZO
%                               optimization problem (specific to Gurobi)

p = parseInput(varargin);

models = p.Results.models;
protein_abundances = p.Results.protein_abundances;
protein_ids = p.Results.protein_ids;
rgr = p.Results.rgr;
bio_idx = p.Results.bio_idx;
pct_rgr = p.Results.pct_rgr;
solver_params = p.Results.solver_params;
prot_sd = p.Results.prot_sd;

% number of conditions
n_cond = numel(models);

if isnan(rgr)
    rgr = ones(n_cond, 1);
end

% model protein IDs
model_rxns = models{1}.varnames;
model_proteins = erase(model_rxns(startsWith(model_rxns, 'draw_prot_')), 'draw_prot_');

% number of reactions and metabolites
[n_model_r, n_model_c] = size(models{1}.A);

% linear constraint matrices of the models
model_mat = cellfun(@(x)full(x.A), models, 'un', 0);
model_rhs = cell2mat(cellfun(@(x)x.rhs, models, 'un', 0)');
model_lb = cell2mat(cellfun(@(x)x.lb, models, 'un', 0)');
model_ub = cell2mat(cellfun(@(x)x.ub, models, 'un', 0)');
model_sense = cell2mat(cellfun(@(x)x.sense, models, 'un', 0)');
model_vtype = cellfun(@(x)x.vtype, models, 'un', 0);
model_vtype = vertcat(model_vtype{:});
model_varnames = cellfun(@(x)x.varnames, models, 'un', 0);
model_varnames = vertcat(model_varnames{:});
model_constrnames = cellfun(@(x)x.constrnames, models, 'un', 0);
model_constrnames = vertcat(model_constrnames{:});
model_objective = cell2mat(cellfun(@(x)x.obj, models, 'un', 0)');

if sum(bio_idx)>0
    model_objective(:) = 0;
    model_objective(bio_idx:n_model_c:end) = 1;
elseif ~any(model_objective)
    error('No objectve was set or provided.')
else
    bio_idx = find(model_objective, 1);
end

% find protein UniProt IDs in the models
idx_exp2model = cell2mat(cellfun(@(x)find(endsWith(model_rxns, x)), protein_ids, 'un', 0));
idx_exp_in_model = find(ismember(protein_ids, model_proteins));
n_matched = numel(idx_exp_in_model);

% number of ratios
n_comb = factorial(n_cond)/(2*factorial(n_cond-2));
n_theta = n_matched*n_comb;

% construct protein ratio matrix
pr_matrix = zeros(2*n_theta, n_cond*n_model_c + n_theta + 1);
pr_rhs = zeros(2*n_theta, 1);
pr_lb = zeros(n_theta, 1);
pr_ub = ones(n_theta, 1);
pr_constrnames = repmat({''}, 2*n_theta, 1);
pr_sense = repelem('<', 2*n_theta, 1);
pr_vtype = repelem('C', n_theta, 1);
pr_varnames = repmat({''}, n_theta, 1);

c = 0;
for i = 1:n_cond-1
    for j = i+1:n_cond
        ratios = protein_abundances(idx_exp_in_model, i) ./ ...
            protein_abundances(idx_exp_in_model, j);
        for k = 1:n_matched
            
            theta_idx = n_cond*n_model_c+c/2+1;
            
            prot_id = protein_ids{idx_exp_in_model(k)};
            pr_varnames{c/2+1} = sprintf('theta_%s_cond_%i_%i', ...
                prot_id, i, j);
            
            % f_ij * e_j - e_i <= theta
            c = c + 1;
            pr_matrix(c, [n_model_c*[i-1 j-i]+idx_exp2model(k) theta_idx]) = ...
                [-1 ratios(k) -1];
            pr_constrnames{c} = sprintf('ratio_%s_1_cond_%i_%i', ...
                prot_id, i, j);
            
            % e_i - f_ij * e_j <= theta
            c = c + 1;
            pr_matrix(c, [n_model_c*[i-1 j-i]+idx_exp2model(k) theta_idx]) = ...
                [1 -ratios(k) -1];
            pr_constrnames{c} = sprintf('ratio_%s_2_cond_%i_%i', ...
                prot_id, i, j);
            
        end
    end
end

% construct biomass ratio matrix
br_matrix = zeros(n_cond, n_cond*n_model_c + n_theta + 1);
br_rhs = zeros(n_cond, 1);
br_constrnames = repmat({''}, n_cond, 1);
br_sense = repelem('=', n_cond, 1);
br_varnames = {'BIOMASS_OBJ'};
br_lb = 0;
br_ub = 1;
br_vtype = 'C';

% RGR ratios relative to first condition
rgr_ratios = rgr./rgr(1);

for i = 1:n_cond
    
    br_matrix(i, n_model_c*(i-1)+bio_idx) = 1;
    br_matrix(i, end) = -rgr_ratios(i);
    br_constrnames{i} = sprintf('biomass_ratio_cond_%i', i);
    
end

% construct the optimization problem
problem = struct;
problem.A = sparse([blkdiag(model_mat{:}), zeros(n_cond*n_model_r, n_theta + 1);
    pr_matrix]);
problem.rhs = [model_rhs; pr_rhs];
problem.lb = [model_lb; pr_lb];
problem.ub = [model_ub; pr_ub];
problem.obj = [zeros(size(model_objective)); zeros(n_theta, 1); 1];
problem.sense = [model_sense; pr_sense];
problem.vtype = [model_vtype; pr_vtype];
problem.constrnames = [model_constrnames; pr_constrnames];
problem.varnames = [model_varnames; pr_varnames];
problem.modelsense = 'max';

% if necessary, combine and extend quadratic constraints
if isfield(models{1}, 'quadcon')
    
    for i=1:n_cond
        
        tmp_qc = zeros(size(problem.A, 2));
        r_idx = find(any(models{i}.quadcon.Qc, 2));
        c_idx = find(any(models{i}.quadcon.Qc, 1));
        tmp_qc((i-1)*n_model_c+r_idx, (i-1)*n_model_c+c_idx) = ...
            models{i}.quadcon.Qc(r_idx, c_idx);
        problem.quadcon(i).Qc = sparse(tmp_qc);
        
        tmp_q = zeros(size(problem.A, 2), 1);
        r_idx_q = find(models{i}.quadcon.q);
        tmp_q((i-1)*n_model_c+r_idx_q) = ...
            models{i}.quadcon.q(r_idx_q);
        problem.quadcon(i).q = sparse(tmp_q);
        problem.quadcon(i).rhs = 0;
        problem.quadcon(i).sense = '<';
        
    end
    clear tmp_qc tmp_q
end

% add biomass ratio constraints
problem.A = [problem.A; br_matrix];
problem.rhs = [problem.rhs; br_rhs];
problem.sense = [problem.sense; br_sense];
problem.constrnames = [problem.constrnames; br_constrnames];
problem.varnames = [problem.varnames; br_varnames];
problem.lb = [problem.lb; br_lb];
problem.ub = [problem.ub; br_ub];
problem.vtype = [problem.vtype; br_vtype];

% optimize biomass production
sol_bio_opt = gurobi(problem, solver_params);

% set minimum percentage for biomass production
min_rgr_pred = min(sol_bio_opt.x(model_objective>0));
problem.lb(problem.obj>0) = pct_rgr*min_rgr_pred;

% minimize sum of absolute deviations from protein ratios
problem.modelsense = 'min';

if isscalar(prot_sd) && isnan(prot_sd)
    theta_weights = ones(n_theta, 1);
else
    % calculate weights for each theta, based on standard deviations
    mean_pred_abundance = mean(mean(sol_bio_opt.x(startsWith(problem.varnames, 'draw_prot_'))));
    prot_sd_scaled = mean_pred_abundance*prot_sd./mean(prot_sd);
    theta_weights = ones(n_matched, n_cond);
    c = 0;
    for i = 1:n_cond-1
        for j = i+1:n_cond
            c = c + 1;
            theta_weights(:, c) = max(prot_sd_scaled(idx_exp_in_model, [i j]), [], 2);
        end
    end
    theta_weights = reshape(theta_weights, n_theta, 1);
end

% weighted objective
problem.obj = [zeros(n_cond*n_model_c, 1); 1./theta_weights; 0];

% sove the optimization problem
solution = gurobi(problem, solver_params);

    function p = parseInput(arguments)
        
        % validation functions
        validateModels = @(x)iscell(x)&isstruct(x{1});
        validateBioIdx = @(x)isnumeric(x)&isscalar(x);
        validateRGRprct = @(x)isnumeric(x)&isscalar(x)&x>=0&x<=1;
        
        p = inputParser;
        p.FunctionName = 'multiRENZO';
        
        addRequired(p, 'models', validateModels)
        addRequired(p, 'protein_abundances', @isnumeric)
        addRequired(p, 'protein_ids', @iscellstr)
        addParameter(p, 'rgr', NaN, @isnumeric)
        addParameter(p, 'bio_idx', NaN, validateBioIdx)
        addParameter(p, 'pct_rgr', 0.5, validateRGRprct)
        addParameter(p, 'solver_params', struct, @isstruct)
        addParameter(p, 'prot_sd', NaN, @isnumeric)
        
        parse(p, arguments{:})
        
        n_conditions = numel(p.Results.models);
        
        if size(p.Results.protein_abundances, 2) ~= n_conditions
            error('The dimension od protein abundances does not match the number of models.')
        end
        
        if numel(p.Results.protein_ids) ~= size(p.Results.protein_abundances, 1)
            error('The dimension of protein abundances does not match the number of protein IDs.')
        end
        
        if ~isnan(p.Results.rgr) && numel(p.Results.rgr) ~= n_conditions
            error('The number of growth rates does not match the number of models.')
        end
        
        if ~isscalar(p.Results.prot_sd) && ~all(all(isnan(p.Results.prot_sd)))
            if ~all(size(p.Results.prot_sd) == size(p.Results.protein_abundances))
                error('The dimension of standard deviations does not correspond to the dimension of protein abundances.')
            end
        end
        
    end

end