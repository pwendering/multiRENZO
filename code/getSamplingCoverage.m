function c = getSamplingCoverage(samples, minval, maxval)
%% c = getSamplingCoverage(samples, minval, maxval)
% Function to calculate the coverage of flux sampling (taken from GAPSPLIT)
% (Keaty & Jensen (2020), doi: 10.1093/bioinformatics/btz971).
%
% INPUT
% double samples:       matrix that contains the sampled fluxes
%                       (#variables x #samples)
% double minflux:       minima of the flux ranges, respectively
% double maxflux:       maxima of the flux ranges, respectively
%
% OUTPUT
% double c:             average coverage

nz_idx = maxval>0;

[n_row, n_col] = size(samples(nz_idx, :));

X = full([...
    reshape(minval(nz_idx), 1, n_row);...
    reshape(maxval(nz_idx), 1, n_row);...
    samples(nz_idx, :)']);
n = (n_col+2);
X = sort(X);
gaps = X(2:n,:) - X(1:n-1,:);
width = max(gaps,[],1);
rel = width ./ (X(n,:) - X(1,:));
c = 1 - mean(rel, 'omitnan');

end