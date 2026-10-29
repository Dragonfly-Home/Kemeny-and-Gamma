function [value, info] = kemeny_feedback(P, k, pi_vec, options)
%KEMENY_FEEDBACK Kemeny's constant from stochastic complements and feedback.
%   P is already ordered; pi_vec is its column stationary distribution.
%   Use RUN_KEMENY for data preparation, parallel setup and benchmarking.
%   options.UseParallel defaults to true. BatchSize (default 64) limits
%   feedback right-hand sides; SCBatchSize (default Inf) limits dense
%   stochastic-complement right-hand sides per batch.

if nargin < 4, options = struct(); end
if ~isfield(options, 'UseParallel'), options.UseParallel = true; end
if ~isfield(options, 'BatchSize'), options.BatchSize = 64; end
if isfield(options, 'SCSolveMode')
    error('KemenyFeedback:RemovedOption', ...
        'Only full-rhs is supported. Remove the SCSolveMode option.');
end
if ~isfield(options, 'SCBatchSize'), options.SCBatchSize = Inf; end
sc_batch_size = options.SCBatchSize;
validateattributes(sc_batch_size, {'numeric'}, {'scalar', 'real', 'positive', 'nonnan'});
if ~isinf(sc_batch_size) && sc_batch_size ~= fix(sc_batch_size)
    error('KemenySC:BatchSize', 'SCBatchSize must be a positive integer or Inf.');
end
sc_batch_size = double(sc_batch_size);
validateattributes(P, {'double'}, {'2d', 'square', 'real', 'finite'});
n = size(P, 1);
validateattributes(k, {'numeric'}, {'scalar', 'integer', '>=', 2, '<=', n});
validateattributes(options.BatchSize, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(options.UseParallel, {'logical', 'numeric'}, {'scalar', 'real', 'finite'});
if ~any(options.UseParallel == [0, 1])
    error('KemenyFeedback:ParallelOption', 'UseParallel must be logical or 0/1.');
end
pi_vec = pi_vec(:);
validateattributes(pi_vec, {'double'}, {'numel', n, 'real', 'finite', 'positive'});
pi_vec = pi_vec / sum(pi_vec);
P = sparse(P);
if any(nonzeros(P) < 0) || norm(sum(P, 2) - 1, inf) > 1e-10
    error('KemenyFeedback:NotStochastic', 'P must be a nonnegative row-stochastic matrix.');
end
if norm(P' * pi_vec - pi_vec, 1) > 1e-10
    error('KemenyFeedback:NotStationary', 'pi_vec is not stationary for P.');
end

% Same contiguous, nearly equal partition as the original Du implementation.
edges = round(linspace(0, n, k + 1));
block_sizes = diff(edges);
labels = repelem((1:k)', block_sizes(:));
J = sparse((1:n)', labels, 1, n, k);
alpha = full(J' * pi_vec);
W = sparse(labels, (1:n)', pi_vec ./ alpha(labels), k, n);
WP = W * P;
C = full(WP * J);
kappa_local = zeros(k, 1);
F = zeros(k);
row_defects = zeros(k, 1);
batch_size = options.BatchSize;

if options.UseParallel
    parfor i = 1:k
        [kappa_local(i), Fpart, row_defects(i)] = ...
            solve_block(P, pi_vec, J, WP, alpha, edges, i, batch_size, sc_batch_size);
        F = F + Fpart;
    end
else
    for i = 1:k
        [kappa_local(i), Fpart, row_defects(i)] = ...
            solve_block(P, pi_vec, J, WP, alpha, edges, i, batch_size, sc_batch_size);
        F = F + Fpart;
    end
end

% F maps centered block observables to feedback vectors. The weighted
% orthonormal-basis sum equals trace(R_C * F), without constructing a basis.
H = eye(k) - ones(k, 1) * alpha';
F = H * F * H;  % remove floating-point components outside the centered space
RC = (eye(k) - C + ones(k, 1) * alpha') \ eye(k);
kappa_c = trace(RC) - 1;
feedback = trace(RC * F);
gamma = kappa_c + feedback;
value = sum(kappa_local) + gamma;
if ~isfinite(value) || ~isreal(value)
    error('KemenyFeedback:InvalidResult', 'The computed Kemeny constant is not finite and real.');
end
info = struct('gamma', gamma, 'kappa_c', kappa_c, 'feedback', feedback, ...
    'kappa_local', kappa_local, 'block_sizes', block_sizes, ...
    'max_stochastic_defect', max(row_defects));
end

function [kappa, Fpart, row_defect] = solve_block(P, pi_vec, J, WP, alpha, edges, i, batch_size, sc_batch_size)
n = size(P, 1);
idx = (edges(i) + 1):edges(i + 1);
other = [1:edges(i), (edges(i + 1) + 1):n];
m = numel(idx);
k = numel(alpha);
s = pi_vec(idx) / alpha(i);
Pout = P(idx, other);
Pin = P(other, idx);
S = full(P(idx, idx));
T = zeros(m, k);
T(:, i) = 1;

% Factor once; solve active coupling columns and update S in batches.
outside = speye(numel(other)) - P(other, other);
active_cols = find(sum(spones(Pin), 1) > 0);
if isempty(active_cols)
    error('KemenyFeedback:Reducible', 'Every block must be reachable from its complement.');
end
[L, U, p, q] = lu(outside, 'vector');
sc_width = min(sc_batch_size, numel(active_cols));

% LU returns outside(p,q)=L*U. Thus Pout(:,q)*X is the same product as
% Pout*temp when temp(q,active_columns)=X, without forming temp.
Pout_q = Pout(:, q);
for sc_first = 1:sc_width:numel(active_cols)
    sc_cols = active_cols(sc_first:min(sc_first + sc_width - 1, numel(active_cols)));
    sc_rhs = full(Pin(p, sc_cols));
    sc_lower = L \ sc_rhs;
    sc_rhs = [];
    sc_upper = U \ sc_lower;
    sc_lower = [];
    sc_delta = Pout_q * sc_upper;
    S(:, sc_cols) = S(:, sc_cols) + sc_delta;
    sc_delta = [];
    sc_upper = [];
end
clear Pout_q outside Pin

% Reuse the outside LU for the k-1 nonzero block-indicator right-hand sides.
% Column i of J(other,:) is zero, so T(:,i)=1 already has its final value.
other_blocks = [1:i-1, i+1:k];
Jother = J(other, :);
for first = 1:batch_size:numel(other_blocks)
    cols = other_blocks(first:min(first + batch_size - 1, numel(other_blocks)));
    forcing = zeros(numel(other), numel(cols));
    forcing(q, :) = U \ (L \ Jother(p, cols));
    T(:, cols) = Pout * forcing;
end
clear L U p q forcing Pout Jother
row_defect = norm(full(sum(S, 2)) - 1, inf);
Q = T - sum(T, 2) * alpha';
Q = Q - ones(m, 1) * (s' * Q);

% Reuse Du's local direct inverse for both the trace and feedback solves.
e = ones(m, 1);
G = inv(speye(m) - S + e * e' / m);
kappa = trace(G) - 1;
% For s'*Q=0, G*Q solves (I-S)*u=Q; centering gives R_S*Q.
potential = G * Q;
potential = potential - e * (s' * potential);
Fpart = full(WP(:, idx) * potential);
end
