function report = test_sc_batches(check_parallel)
%TEST_SC_BATCHES Check full-RHS batches against independent small references.
%   No pool is created. Pass true to also check parallel New and original Du.
if nargin < 1, check_parallel = false; end
validateattributes(check_parallel, {'logical'}, {'scalar'});
saved_path = path;
path_cleanup = onCleanup(@() path(saved_path)); %#ok<NASGU>
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'), fullfile(root, 'third_party', 'du'), ...
    fullfile(root, 'third_party', 'bini'));
if check_parallel && (isempty(ver('parallel')) || isempty(gcp('nocreate')))
    error('KemenyBatchTest:NoPool', 'Parallel checks require an already running pool.');
end
saved_rng = rng;
rng_cleanup = onCleanup(@() rng(saved_rng)); %#ok<NASGU>
rng(20261009);
A = 0.1 + rand(11);
positive = A ./ sum(A, 2);
cycle = circshift(eye(9), [0, 1]);
A = diag(1:11, 1); A = A + A'; A(1,12) = 2; A(12,1) = 2;
weighted = A ./ sum(A, 2);
% Directed positive, periodic directed, singleton blocks, and weighted graph.
% k=3 gives size-four blocks, exercising a partial final width-three batch.
cases = {positive, 3; cycle, 4; cycle, 9; weighted, 5};
configurations = 0;
max_error = 0;
du_checks = 0;
for c = 1:size(cases, 1)
    P = sparse(cases{c, 1}); k = cases{c, 2}; n = size(P, 1);
    pi_vec = [full(P')-eye(n); ones(1,n)] \ [zeros(n,1); 1];
    pi_vec = pi_vec / sum(pi_vec);
    reference = trace(inv(eye(n)-full(P)+ones(n,1)*pi_vec')) - 1; %#ok<MINV>
    edges = round(linspace(0, n, k+1));
    locals = zeros(k, 1);
    for i = 1:k
        ids = edges(i)+1:edges(i+1); other = setdiff(1:n, ids);
        S = full(P(ids,ids)) + full(P(ids,other)) * ...
            ((eye(numel(other))-full(P(other,other))) \ full(P(other,ids)));
        s = pi_vec(ids)/sum(pi_vec(ids));
        locals(i) = trace(inv(eye(numel(ids))-S+ones(numel(ids),1)*s')) - 1; %#ok<MINV>
    end
    for batch = [Inf, 1, 3]
        options = struct('UseParallel', false, 'SCBatchSize', batch, 'BatchSize', 2);
        [value, info] = kemeny_feedback(P, k, pi_vec, options);
        max_error = max(max_error, check_values(value, info, reference, locals));
        if check_parallel
            options.UseParallel = true;
            [value, info] = kemeny_feedback(P, k, pi_vec, options);
            max_error = max(max_error, check_values(value, info, reference, locals));
        end
        configurations = configurations + 1;
    end
    % Check the default SCBatchSize as well as explicit widths.
    [value, info] = kemeny_feedback(P, k, pi_vec, struct('UseParallel', false));
    max_error = max(max_error, check_values(value, info, reference, locals));
    if check_parallel
        [du_value, du_elapsed] = kemeny_k_parallel(P, k, pi_vec, struct('solver', 'lu'));
        du_error = abs(du_value-reference) / max(1, abs(reference));
        assert(isfinite(du_value) && du_error < 1e-10);
        assert(isfinite(du_elapsed) && du_elapsed >= 0);
        max_error = max(max_error, du_error);
        du_checks = du_checks + 1;
    end
end
% Unsupported options must fail explicitly.
retired_option_rejected = false;
try
    kemeny_feedback(P, k, pi_vec, struct('UseParallel', false, 'SCSolveMode', 'full-rhs'));
catch exception
    retired_option_rejected = strcmp(exception.identifier, 'KemenyFeedback:RemovedOption');
end
assert(retired_option_rejected);
report = struct('passed', true, 'configurations', configurations, ...
    'parallel_checked', check_parallel, 'original_du_checks', du_checks, ...
    'max_scaled_error', max_error);
fprintf('PASS: %d New full-RHS batch configurations; %d original Du checks; maximum scaled error %.3e.\n', ...
    configurations, du_checks, max_error);
end

function error_value = check_values(value, info, reference, locals)
expected = [reference; locals; reference-sum(locals)];
actual = [value; info.kappa_local; info.gamma];
error_value = max(abs(actual-expected) ./ max(1, abs(expected)));
assert(all(isfinite(actual)) && error_value < 1e-10);
assert(info.max_stochastic_defect < 1e-10);
assert(~isfield(info, 'timing'));
assert(abs(info.gamma-info.kappa_c-info.feedback) < 1e-10);
end
