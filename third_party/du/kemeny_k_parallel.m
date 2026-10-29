function [kemeny_val, comp_time] = kemeny_k_parallel(P, k, pi_vec, options)
% KEMENY_K_PARALLEL Computes Kemeny's constant using k-partitions algorithm
    
    if nargin < 4
        options = struct();
    end
    if ~isfield(options, 'solver')
        options.solver = 'lu'; 
    end
    
    t_start = tic;
    n = size(P, 1);
    
    % =========================================================
    % 1. Fast partitioning strategy
    % =========================================================
    block_sizes = diff(round(linspace(0, n, k + 1)));
    partitions = mat2cell(1:n, 1, block_sizes);
    
    kappa_S = zeros(k, 1);
    norm_pi = zeros(k, 1);
    
    % =========================================================
    % 2. Compute local Kemeny constants (using parfor to fully utilize multi-core CPUs)
    % =========================================================
    parfor i = 1:k
        idx_i = partitions{i};
        idx_not_i = setdiff(1:n, idx_i);
        
        P_ii = P(idx_i, idx_i);
        P_i_not = P(idx_i, idx_not_i);
        P_not_i = P(idx_not_i, idx_i);
        P_not_not = P(idx_not_i, idx_not_i);
        
        I_minus_P_not = speye(length(idx_not_i)) - P_not_not;
        active_cols = find(sum(spones(P_not_i), 1) > 0); 
        
        temp = zeros(length(idx_not_i), length(idx_i));
        
        if ~isempty(active_cols)
            [L, U, p, q] = lu(I_minus_P_not, 'vector');
            temp(q, active_cols) = U \ (L \ P_not_i(p, active_cols));
        end
        
        S_ii = P_ii + P_i_not * temp;
        kappa_S(i) = kemenydirect(S_ii);
        norm_pi(i) = sum(pi_vec(idx_i));
    end
    
    % =========================================================
    % 3. Fast computation of the global interaction term Gamma
    % =========================================================
    M_global = speye(n) - P;
    M_global(1, 1) = M_global(1, 1) + 1; % Minimal perturbation to break singularity and maintain extreme sparsity
    
    [L_M, U_M, p_M, q_M] = lu(M_global, 'vector');
    
    % Batch construct right-hand sides
    C_mat = -ones(n, k-1) .* norm_pi(1:k-1)';
    for j = 1:k-1
        C_mat(partitions{j}, j) = C_mat(partitions{j}, j) + 1;
    end
    
    % Fast vectorized back-substitution to instantly solve for vectors corresponding to all j
    X_mat = zeros(n, k-1);
    X_mat(q_M, :) = U_M \ (L_M \ C_mat(p_M, :));
    
    gamma_j_vals = zeros(k-1, 1);
    norm_pi_k = norm_pi(k);
    s_k = sparse(partitions{k}, 1, pi_vec(partitions{k}) / norm_pi_k, n, 1);
    
    % O(n) vector dot products for fast merging of results
    for j = 1:k-1
        norm_pi_j = norm_pi(j);
        s_j = sparse(partitions{j}, 1, pi_vec(partitions{j}) / norm_pi_j, n, 1);
        b_j = s_j - s_k;
        
        gamma_j_vals(j) = b_j' * X_mat(:, j);
    end
    
    gamma = sum(gamma_j_vals);
    kemeny_val = sum(kappa_S) + gamma;
    
    comp_time = toc(t_start);
end