function [problem, adjustedParams] = addCO2UptakeConstraints(varargin)
%% [problem, adjustedParams] = addCO2UptakeConstraints(varargin)
% This function creates a quadratically-constraint optimization problem by
% adding constraints on CO2 uptake based on the C3 photosynthesis model
% described by Farquhar et al. (1980) (FvCB model).
% Input:
%       struct model:           metabolic model with additional fields for
%                               temperature-adjustment (see 'createTGEM')
%       double T:               temperature in Kelvin for adjustment of
%                               temperature-dependent parameters
%       double I:               irradiance (umol/m2/s)
%       logical saBool:         (optional) specifies whether or not a
%                               sensitivity analysis is performed (defaults
%                               to false)
%       double saPercentage:    (optional) if saBool is true, specifies the
%                               percentage for increase of the parameter
%                               given in saParameter at the given
%                               temperature
%       char saParameter:       (optional) name of the parameter that
%                               should be increased by saPercentage
% Output:
%       struct problem:         problem structure suited to be solved using
%                               the Gurobi solver
%       struct adjustedParams:  temperature-adjusted paramers
%
% Philipp Wendering, University of Potsdam (philipp.wendering@gmail.com) 23/03/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

validTemp = @(T) isnumeric(T) && (T >= 273.15) && (T <= 373.15);
p = inputParser;

addRequired(p,'model',@validateTGEM);
addRequired(p,'T',validTemp)
addRequired(p,'I',@isnumeric)
addParameter(p,'saBool',false,@islogical)
addParameter(p,'saPercentage',0,@isnumeric)
addParameter(p,'saParameter','',@ischar)

parse(p,varargin{:});

model = p.Results.model;
T = p.Results.T;
I = p.Results.I;
saBool = p.Results.saBool;
saPercentage = p.Results.saPercentage;
saParameter = p.Results.saParameter;

if ~isfield(model, 'TRANS_ID')
    model.TRANS_ID = [];
end

if ~isfield(model, 'vartype')
    model.vartype = repmat('C', size(model.S, 2), 1);
end


%% define variable indices
[ROW_DIM,COL_DIM] = size(model.S);
RXNS = numel(model.rxns);
C_IDX = findRxnIDs(model,model.C_ID);
O_IDX = findRxnIDs(model,model.O_ID);
A_IDX = COL_DIM + 1;
Z_IDX = COL_DIM + 2;
RESP_IDX = findRxnIDs(model,model.RESP_ID);
ABS_IDX = findRxnIDs(model,model.ABS_ID);
BIO_IDX = findRxnIDs(model,model.BIO_ID);
TRANS_IDX = findRxnIDs(model,model.TRANS_ID);
if TRANS_IDX == 0
    TRANS_IDX = [];
end

% estimate leaf dry mass per area (LMA) at the given temperature [g(DW) m^-2]
LMA = model.lma_model(T);

%% adjust kinetic and diffusion parameters to temperature

% increase selected parameter when performing sensitivity analysis
% increase base parameters before temperature adjustment (if possible)
if saBool && ismember(saParameter, {'V_c_max', 'K_c', 'K_o', 'k_c', ...
        'k_o', 'J_max', 'g_m', 'LMA'})
    if isequal(saParameter, 'LMA')
        LMA = (1+saPercentage) * LMA;
    else
        model.(saParameter) = (1+saPercentage) * model.(saParameter);
    end
end

% conversion factors between model units and C3 photosynthesis model parameters
fvcb_scaling = config('fvcb_scaling');
fvcb2fba = fvcb_scaling * 3600 / 1000 / LMA;

% Michaelis constant for carboxylation reaction of RuBisCO [ubar]
K_c = adj_Kc('T',T);
% Michaelis constant for oxygenation reaction of RuBisCO [ubar]
K_o = adj_Ko('T',T);
% Turnover number for carboxylation reaction of RuBisCO [s^-1]
k_c = adj_k_c('k25',model.k_c,'T',T);
% Turnover number for oxygenation reaction of RuBisCO [s^-1]
k_o = adj_k_o('k25',model.k_o,'T',T);
% Maximal reaction velocity for carboxylation [umol m-2 s-1]
V_c_max = adj_v_c_max('k25', model.V_c_max, 'T', T);
% Mesophyll conductance [mol bar-1 m-2 s-1]
g_m = adj_g_m('T',T,'g_m_ref',model.g_m,'method',config('g_m_method'));
% change g_m if we are running a sensitivity analysis
if saBool && isequal(saParameter, 'g_m') ...
        && ismember(config('g_m_method'), {'caemmerer', 'empirical'})
    % only change g_m here if the temperature adjustment does not require a
    % reference value
     g_m = (1+saPercentage) * g_m;
end
% Stomatal conductance [mol m^-2 s^-1]
g_s = adj_g_s(T);
% change g_s if we are running a sensitivity analysis
if saBool && isequal(saParameter, 'g_s')
     g_s = (1+saPercentage) * g_s;
end
% Saturation vapour pressure as proxy for intercellular water vapour pressure [bar]
[~, e_i] = vpd('T', T);
e_i = pascal2bar(e_i);
% Ambient vapour pressure [bar]
e_a = e_i - 0.01;
% atmospheric pressure [bar]
P = model.P;
% boundary layer resistance to water vapour (as in Farquhar and Wong 1984) [m^2 s mol^-1]
r_b = model.r_b;
% change r_b if we are running a sensitivity analysis
if saBool && isequal(saParameter, 'r_b')
     r_b = (1+saPercentage) * r_b;
end
% stomatal resistance [m^2 s mol^-1]
r_s = 1/g_s;
% total resistance to CO2 [m^2 s mol^-1]
r = 1.37*r_b + 1.6*r_s;
% reciprocal sum of CO2 resistances [mol m^-2 s^-1]
g = 1/(r_b + r_s);
% transpiration rate [mol m^-2 s^-1]
E = g*(e_i-e_a)/(P-0.5*e_i+e_a);
% set lower bound on H2O sink reaction
model.lb(TRANS_IDX) = E * 1e-6 * fvcb2fba;
% Specificity of RuBisCO for CO2 over O2
S_co = (k_c/K_c)*(K_o/k_o);
% ratio of oxygenation rate to carboxylation rate
gamma = model.O_a / model.C_a;
phi = gamma/S_co;
% light saturated potential rate of electron transport
J_max = adj_j_max('J_max_ref', model.J_max,'T', T);
% released CO2 per oxygenation
alpha = model.alpha;
% absorptance / quantum yield of red light
absorptance = model.absorptance;
% correction factor for the spectral quality of light
sqcf = model.sqcf;
% irradiance that can be used by the leaf
I2 = I*absorptance*(1-sqcf)/2;
% convexity factor for the irradiance-J relationship
theta = config('theta');
% additional variable Y
Y = alpha * gamma / S_co;
% potential electron transport rate
J = (I2 + J_max - sqrt((I2+J_max)^2 - 4*theta*I2*J_max)) / (2*theta);
% right-hand side for the electron transport rate-limited net CO2
% assimilation rate
A_j_rhs = J*(1-Y)/(4+8*Y); % [umol m^-2 s^-1]
A_j_rhs = A_j_rhs * fvcb2fba; % [mmol gDW^-1 h^-1]

% increase selected parameter when performing sensitivity analysis
% increase parameter if it does not have a reference value
if saBool && isequal(saParameter, 'phi')
    phi = (1+saPercentage) * phi;
end

% create output structure for adjusted parameters
adjustedParams.Kc = K_c;
adjustedParams.Ko = K_o;
adjustedParams.kc = k_c;
adjustedParams.ko = k_o;
adjustedParams.V_c_max = V_c_max;
adjustedParams.phi = phi;
adjustedParams.S_co = S_co;
adjustedParams.g_m = g_m;
adjustedParams.g_s = g_s;
adjustedParams.LMA = LMA;
adjustedParams.J_max = J_max;

%% Create Gurobi problem structure that includes constraints on CO2 uptake
% constants in constraints
c1 = 1/P/r; % mol m^-2 s^-1 bar^-1
c2 = E/2/P; % mol m^-2 s^-1 bar^-1
c3 = c1 - c2; % mol m^-2 s^-1 bar^-1
c4 = c1 + c2; % mol m^-2 s^-1 bar^-1
c5 = gamma*K_c/K_o; % unitless
c6 = (1 + c5)/ g_m; % bar m^2 s mol^-1
% allowed deviation of phi from the calculated value
tau = config('tau');

% re-scale parameters
% m^2      -> g (using LMA)
% s^-1     -> h^-1
% umol^-1  -> mmol^-1
V_c_max = V_c_max * fvcb2fba;
g_m = g_m * fvcb2fba;
c3 = c3 * fvcb2fba;
c4 = c4 * fvcb2fba;
c6 = c6 / fvcb2fba;

% initialize problem structure
problem = struct;

% ~~~~~~~~~~~~~~~~~~~~~~~ linear constraint matrix ~~~~~~~~~~~~~~~~~~~~~~ %
% constaints for limits of phi
phi_constraints = zeros(2, Z_IDX);
phi_constraints(1,[C_IDX O_IDX]) = [(1-tau)*phi -1];
phi_constraints(2,[C_IDX O_IDX]) = [-(1+tau)*phi 1];

% constraints on CO2 assimilation rate
A_c_constraint = zeros(1,Z_IDX);
A_c_constraint([C_IDX O_IDX RESP_IDX A_IDX]) = [1 -alpha -1 -1];

A_j_constraint = zeros(1,Z_IDX);
A_j_constraint([C_IDX O_IDX]) = [1 -alpha];

% define additional variable Z
Z_constraint = zeros(1,Z_IDX);
Z_constraint([A_IDX Z_IDX]) = [1 c4];

% add linear constraints to the given stoichiometric matrix
problem.A = sparse([
    [model.S, zeros(ROW_DIM,2)];
    A_c_constraint;
    A_j_constraint;
    Z_constraint;
    phi_constraints
    ]);

% ~~~~~~~~~~~~~~~~~~~~~~~ linear objective vector ~~~~~~~~~~~~~~~~~~~~~~~ %
problem.obj = [model.c;0;0];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ objective sense ~~~~~~~~~~~~~~~~~~~~~~~~~~ %
problem.modelsense = model.osenseStr;

% ~~~~~~~~~~~~~~~~~~~~ sense of linear constraints ~~~~~~~~~~~~~~~~~~~~~~ %
origSense = translateSenseCobra2Gurobi(model.csense);
problem.sense = [origSense;'=';'<';'=';'<';'<'];

% ~~~~~~~~~~~~~~~~ right-hand side of linear constraints ~~~~~~~~~~~~~~~~ %
problem.rhs = [model.b; 0; A_j_rhs; model.C_a*c3; 0; 0];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~ constraint names ~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
problem.constrnames = [...
    model.mets;
    {'A_C_CONSTRAINT'};...
    {'A_J_CONSTRAINT'};...
    {'Z_CONSTRAINT'};...
    {'PHI_CONSTRAINT_1'};...
    {'PHI_CONSTRAINT_2'}];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ lower bounds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
problem.lb = [model.lb; -1e6; 0];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ upper bounds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
problem.ub = [model.ub; 1e6; 1e6];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ variable names ~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
problem.varnames = repmat({''},COL_DIM+2,1);
problem.varnames(1:RXNS) = model.rxns;
problem.varnames(RXNS+1:RXNS+2) = {'A';'Z'};
problem.varnames(BIO_IDX) = {'GROWTH'};
problem.varnames(ABS_IDX) = {'PHOTON_UPTAKE'};
problem.varnames(TRANS_IDX) = {'TRANSPIRATION'};

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~ variable types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
problem.vtype = [model.vartype;'C';'C'];

% ~~~~~~~~~~~~~~~~~~~~~~~~~ quadratic constraints ~~~~~~~~~~~~~~~~~~~~~~~ %

% quadratic part
problem.quadcon.Qc = zeros(Z_IDX);
problem.quadcon.Qc(A_IDX, C_IDX) = -c6;
problem.quadcon.Qc(Z_IDX, C_IDX) = 1+c5;
problem.quadcon.Qc = sparse(problem.quadcon.Qc);

% linear part
problem.quadcon.q = sparse(Z_IDX,1);
problem.quadcon.q([C_IDX A_IDX Z_IDX]) = [K_c V_c_max/g_m -V_c_max];

% right-hand side
problem.quadcon.rhs = 0;

% constraint sense
problem.quadcon.sense = '<';

end