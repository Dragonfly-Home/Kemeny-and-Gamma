function report = test_kemeny_feedback(check_parallel)
%TEST_KEMENY_FEEDBACK Independent small-matrix correctness checks.
%   Serial by default; pass true to also check an already running pool.
%   Returns a test summary. Independent references use the full fundamental
%   matrix and a weighted orthonormal contrast basis.

if nargin < 1, check_parallel = false; end
validateattributes(check_parallel, {'logical'}, {'scalar'});
saved_path = path;
path_cleanup = onCleanup(@() path(saved_path)); %#ok<NASGU>
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
saved_rng = rng;
rng_cleanup = onCleanup(@() rng(saved_rng)); %#ok<NASGU>
rng(20260929, 'twister');

if check_parallel
    if isempty(ver('parallel')) || isempty(gcp('nocreate'))
        error('KemenyFeedbackTest:NoPool', ...
            ['Parallel validation requires an already running pool. ', ...
             'Start the desired pool first, or call this test without true.']);
    end
end

% Dense, nonreversible transition matrix with a known irreducible support.
A = 0.05 + rand(9);
P_directed = A ./ sum(A, 2);

% Connected, weighted reversible graph, including unequal vertex degrees.
A = zeros(10);
for i = 1:10
    j = mod(i, 10) + 1;
    A(i, j) = 0.25 + rand();
    A(j, i) = A(i, j);
end
A(1, 6) = 1.7; A(6, 1) = 1.7;
A(2, 8) = 0.4; A(8, 2) = 0.4;
P_reversible = A ./ sum(A, 2);

% A connected bipartite random walk of exact period two.
B = 0.1 + rand(3, 5);
A = [zeros(3), B; B', zeros(5)];
P_bipartite = A ./ sum(A, 2);

% Irreducible directed cycle, with no self loops and exact period nine.
% Its sparse boundary also exercises blocks with inactive coupling columns.
P_cycle = circshift(eye(9), [0, 1]);

cases = {
    'nonreversible', P_directed, [2, 4, 9];
    'reversible weighted graph', P_reversible, [2, 3, 10];
    'period-two bipartite graph', P_bipartite, [2, 3, 8];
    'directed cycle', P_cycle, [2, 4, 9]
};

tolerance = 1e-10;
max_relative_error = 0;
configurations = 0;
assertions = 0;
fprintf('Checking Kemeny feedback on independent small-matrix references...\n');

for c = 1:size(cases, 1)
    name = cases{c, 1};
    P = cases{c, 2};
    n = size(P, 1);
    pi_vec = stationary_reference(P);
    reference_k = fundamental_kemeny(P, pi_vec);
    for k = cases{c, 3}
        label = sprintf('%s, n=%d, k=%d', name, n, k);
        [value, info] = kemeny_feedback(sparse(P), k, pi_vec, ...
            struct('UseParallel', false));
        [basis_gamma, reference_locals] = manuscript_gamma(P, pi_vec, k);
        errors = [
            check_close(value, reference_k, tolerance, [label, ': global K']);
            check_close(info.gamma, basis_gamma, tolerance, ...
                        [label, ': original contrast-basis gamma']);
            check_close(info.kappa_local, reference_locals, tolerance, ...
                        [label, ': local fundamental matrices']);
            check_close(info.gamma, reference_k - sum(reference_locals), ...
                        tolerance, [label, ': global minus local correction'])
        ];
        assertions = assertions + numel(errors);

        [value_one, info_one] = kemeny_feedback(sparse(P), k, pi_vec, ...
            struct('UseParallel', false, 'BatchSize', 1));
        errors(end + 1) = check_close(value_one, value, tolerance, ...
            [label, ': BatchSize=1 K']); %#ok<AGROW>
        errors(end + 1) = check_close(info_one.gamma, info.gamma, tolerance, ...
            [label, ': BatchSize=1 gamma']); %#ok<AGROW>
        assertions = assertions + 2;

        % A size-two batch leaves a partial final batch when k-1 is odd.
        [value_two, info_two] = kemeny_feedback(sparse(P), k, pi_vec, ...
            struct('UseParallel', false, 'BatchSize', 2));
        errors(end + 1) = check_close(value_two, reference_k, tolerance, ...
            [label, ': BatchSize=2 global K']); %#ok<AGROW>
        errors(end + 1) = check_close(info_two.gamma, basis_gamma, tolerance, ...
            [label, ': BatchSize=2 contrast-basis gamma']); %#ok<AGROW>
        assertions = assertions + 2;

        expected_sizes = diff(round(linspace(0, n, k + 1)));
        assert(isequal(info.block_sizes, expected_sizes), ...
            'KemenyFeedbackTest:Partition', '%s: incorrect block sizes.', label);
        assertions = assertions + 1;
        if k == n
            errors(end + 1) = check_close(info.kappa_local, zeros(n, 1), ...
                tolerance, [label, ': singleton constants']); %#ok<AGROW>
            errors(end + 1) = check_close(info.feedback, 0, tolerance, ...
                [label, ': singleton feedback']); %#ok<AGROW>
            assertions = assertions + 2;
        end

        if check_parallel
            [parallel_value, parallel_info] = kemeny_feedback(sparse(P), ...
                k, pi_vec, struct('UseParallel', true, 'BatchSize', 2));
            errors(end + 1) = check_close(parallel_value, value, tolerance, ...
                [label, ': parallel versus serial K']); %#ok<AGROW>
            errors(end + 1) = check_close(parallel_info.gamma, info.gamma, ...
                tolerance, [label, ': parallel versus serial gamma']); %#ok<AGROW>
            assertions = assertions + 2;
        end
        max_relative_error = max(max_relative_error, max(errors));
        configurations = configurations + 1;
        fprintf('  PASS  %-45s K = %.12g\n', label, value);
    end
end

report = struct('passed', true, 'configurations', configurations, ...
    'assertions', assertions, 'parallel_checked', check_parallel, ...
    'tolerance', tolerance, 'max_scaled_error', max_relative_error);
fprintf('PASS: %d configurations, %d checks; maximum scaled error %.3e.\n', ...
    configurations, assertions, max_relative_error);
end

function pi_vec = stationary_reference(P)
% Independent normalized stationary solve, used only for tiny test inputs.
n = size(P, 1);
pi_vec = [P' - eye(n); ones(1, n)] \ [zeros(n, 1); 1];
pi_vec = pi_vec / sum(pi_vec);
assert(all(pi_vec > 0) && norm(P' * pi_vec - pi_vec, 1) < 1e-12, ...
    'KemenyFeedbackTest:Stationary', 'Reference stationary solve failed.');
end

function value = fundamental_kemeny(P, pi_vec)
% Explicit inverse is deliberately acceptable for these small references.
n = size(P, 1);
value = trace(inv(eye(n) - P + ones(n, 1) * pi_vec')) - 1; %#ok<MINV>
end

function [gamma, local_values] = manuscript_gamma(P, pi_vec, k)
% Implement the original formula literally, one contrast a_r at a time.
n = size(P, 1);
edges = round(linspace(0, n, k + 1));
blocks = cell(k, 1);
outside = cell(k, 1);
conditionals = cell(k, 1);
local_fundamentals = cell(k, 1);
alpha = zeros(k, 1);
J = zeros(n, k);
C = zeros(k);
local_values = zeros(k, 1);
for i = 1:k
    blocks{i} = (edges(i) + 1):edges(i + 1);
    outside{i} = setdiff(1:n, blocks{i});
    alpha(i) = sum(pi_vec(blocks{i}));
    conditionals{i} = pi_vec(blocks{i}) / alpha(i);
    J(blocks{i}, i) = 1;
end
for i = 1:k
    ids = blocks{i};
    other = outside{i};
    S = P(ids, ids) + P(ids, other) * ...
        ((eye(numel(other)) - P(other, other)) \ P(other, ids));
    local_fundamentals{i} = inv(eye(numel(ids)) - S + ...
        ones(numel(ids), 1) * conditionals{i}'); %#ok<MINV>
    local_values(i) = trace(local_fundamentals{i}) - 1;
    for j = 1:k
        C(i, j) = conditionals{i}' * P(ids, blocks{j}) * ...
            ones(numel(blocks{j}), 1);
    end
end
RC = inv(eye(k) - C + ones(k, 1) * alpha'); %#ok<MINV>
% Transform a Euclidean basis perpendicular to sqrt(alpha) into a basis
% orthonormal in <a,b>_alpha, without using the production projector H.
basis = null(sqrt(alpha)');
basis = basis ./ sqrt(alpha);
assert(norm(basis' * diag(alpha) * basis - eye(k - 1), 'fro') < 1e-12, ...
    'KemenyFeedbackTest:Basis', 'Reference contrast basis is not orthonormal.');
assert(norm(alpha' * basis) < 1e-12, ...
    'KemenyFeedbackTest:Basis', 'Reference contrast basis is not centered.');
gamma = trace(RC) - 1;
for r = 1:k - 1
    a = basis(:, r);
    f = J * a;
    u = zeros(n, 1);
    for j = 1:k
        ids = blocks{j};
        other = outside{j};
        q = f(ids) + P(ids, other) * ...
            ((eye(numel(other)) - P(other, other)) \ f(other));
        u(ids) = local_fundamentals{j} * q;
    end
    b = zeros(k, 1);
    for i = 1:k
        for j = 1:k
            b(i) = b(i) + conditionals{i}' * ...
                P(blocks{i}, blocks{j}) * u(blocks{j});
        end
    end
    gamma = gamma + a' * diag(alpha) * RC * b;
end
end

function err = check_close(actual, expected, tolerance, label)
actual = actual(:);
expected = expected(:);
assert(numel(actual) == numel(expected) && all(isfinite(actual)) && ...
    all(isfinite(expected)), 'KemenyFeedbackTest:InvalidValue', ...
    '%s: invalid or differently sized results.', label);
err = max(abs(actual - expected) ./ max(1, abs(expected)));
assert(err <= tolerance, 'KemenyFeedbackTest:Mismatch', ...
    '%s: scaled error %.3e exceeds %.3e.', label, err, tolerance);
end
