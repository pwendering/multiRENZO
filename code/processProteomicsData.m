function processProteomicsData
% This function reads and filters relative proteomics data, maps them to 
% the proteins in the metabolic model, and writes the average reporter
% intensities to a Matlab workspace.

top_dir = regexprep(mfilename('fullpath'), '(\\[^\\]+){2}$', '');
data_dir = fullfile(top_dir, 'data');
result_dir = fullfile(top_dir, 'results');

%% Load temperature-dependent ecAraCore model
tmp = load(config('tgemFile'));
model = tmp.TGEM;
clear tmp

prot_idx = find(startsWith(model.rxns, 'draw_prot_'));
protein_IDs = erase(model.rxns(prot_idx), 'draw_prot_');
n_prot = numel(prot_idx);

%% Load proteomics data
% filter out proteins with 
% * Q-value > 0.05
% * OnlyIdentifiedBySite == "+"
% * (PotentialContaminant) (no entries)
% * (Reverse) (no entries)
% * Peptide counts (razor+unique) <= 1

prot_input_file = fullfile(data_dir, 'peptide_counts.tsv');

prot_measurements = readtable(prot_input_file, ...
    'Delimiter', '\t', ...
    'FileType', 'text');

filter_idx = prot_measurements.Q_value >= 0.05 | ...
    ismember(prot_measurements.OnlyIdentifiedBySite, '+') | ...
    prot_measurements.PeptideCounts_razor_unique_ <= 1;

fprintf('Excluding %i proteins from further analysis.\n', sum(filter_idx))

prot_measurements_red = prot_measurements(~filter_idx, :);

% obtain Uniprot IDs from AGI numbers
agis = cellfun(@(x){strjoin(unique(strsplit(x, ';')), ';')}, ...
    regexprep(prot_measurements_red.ProteinIDs, '\.\d+', ''));

% read translation between AGIs and UniProt IDs
id_translation_tab = readtable(fullfile(data_dir, 'idmapping_2025_05_02.tsv'), ...
    'Delimiter', '\t', 'FileType', 'text');

up_ids = cellfun(@(x){strjoin(id_translation_tab.Entry( ...
    ismember(id_translation_tab.From, strsplit(x, ';'))), ';')}, ...
    agis);

% map protein IDs from experiment to model proteins
exp2model = cellfun(@(x)find(ismember(protein_IDs, strsplit(x, ';'))), ...
    up_ids, 'un', 0);

fprintf('%i proteins could be mapped to the model.\n', ...
    sum(~cellfun(@isempty, exp2model)))
fprintf('%i out of these proteins have multiple matches in the model.\n', ...
    sum(cellfun(@numel, exp2model)>1))

% map peptide counts to model proteins
reporter_intensities_model = zeros(n_prot, 16);
for i = 1:16
    column = sprintf('ReporterIntensityCorrected%iWalkerTMT16', i);
    tmp_intensities = cellfun(@(x)...
        max(prot_measurements_red.(column)(contains(up_ids, x))), ...
        protein_IDs, 'un', 0);
    tmp_intensities(cellfun(@isempty, tmp_intensities)) = {NaN};
    reporter_intensities_model(:, i) = cell2mat(tmp_intensities);
end

%% Calculate average intensities and standard deviations

% calculate means over the replicates
mean_intensities = arrayfun(@(i) ...
    mean(reporter_intensities_model(:, 4*(i-1)+1:4*i), 2), 1:4, 'un', 0);
mean_intensities = cell2mat(mean_intensities);

% re-order mean intensities such that the order corresponds to 21%, 2%, 40%
% O2
mean_intensities = mean_intensities(:, [3, 2, 4]);

% calculate standard deviations over the replicates
sd_intensities = arrayfun(@(i) ...
    std(reporter_intensities_model(:, 4*(i-1)+1:4*i), [], 2), 1:4, 'un', 0);
sd_intensities = cell2mat(sd_intensities);

% re-order standard deviations such that the order corresponds to 
% 21%, 2%, 40% O2
sd_intensities = sd_intensities(:, [3, 2, 4]);

% calculate scaled standard deviations over the replicates
sd_intensities_scaled = arrayfun(@(i) ...
    std(reporter_intensities_model(:, 4*(i-1)+1:4*i) ...
    ./max(reporter_intensities_model(:, 4*(i-1)+1:4*i)), [], 2), 1:4, ...
    'un', 0);
sd_intensities_scaled = cell2mat(sd_intensities_scaled);

% re-order scaled standard deviations such that the order corresponds to
% 21%, 2%, 40% O2
sd_intensities_scaled = sd_intensities_scaled(:, [3, 2, 4]);

%% Save as Matlab Workspace
save(fullfile(result_dir, 'measured_protein_abundances'),...
    'mean_intensities', ...
    'sd_intensities', ...
    'protein_IDs', ...
    'sd_intensities_scaled')

end