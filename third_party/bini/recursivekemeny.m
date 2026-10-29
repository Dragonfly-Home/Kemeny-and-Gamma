function k = recursivekemeny(P,pi)
%RECURSIVEKEMENY Divide and conquer with recursion for the computation of
%the Kemeny constant.
%   INPUT:  P stochastic matrix
%           pi stationary vector of the matrix P
%   OUTPUT: k Kemeny constant (approximation)

n = size(P,1);
pi = pi/sum(pi);

if n <= 500
    k = kemenydirect(P);
    return
else
    m = round(n/2);
    pi1 = pi(1:m)/sum(pi(1:m));
    pi2 = pi(m+1:n)/sum(pi(m+1:n));
    % Construction of the blocks
    P11 = P(1:m,1:m);
    P12 = P(1:m,m+1:n);
    P21 = P(m+1:n,1:m);
    P22 = P(m+1:n,m+1:n);
    I2 = speye(m,m);
    I1 = speye(n-m,n-m);
    % Computation of the Schur complements
    P1 = P11 + P12*((I1-P22)\P21);
    P2 = P22 + P21*((I2-P11)\P12);
    % Recursive calls:
    % This may be problematic if the recursion becomes to deep!
    k1 = recursivekemeny(P1,pi1);
    k2 = recursivekemeny(P2,pi2);
    % Linear-system solution with preconditioned GMRES
    options.type = "nofill";
    [L,U] = ilu(I1-P22,options);
    [L2,U2] = ilu(I2 - P11 - P12*(spdiags(diag(I1-P22),0,n-m,n-m)\P21), ...
        options);
    [x,~] = gmres(I1 - P22,ones(n-m,1),[],1e-6,n-m,L,U);
    M = @(x) x - P1*x + (pi1'*x);
    [y,~] = gmres(M,P12*x,[],1e-6,m,L2,U2);
    [y,~] = gmres(I1 - P22,P21*y,[],1e-6,n-m,L,U);
    theta = pi2'*(x+y);
    gamma = sum(pi(1:m))*theta-sum(pi(m+1:n));
    % Kemeny constant by update
    k = k1+k2+gamma;
    return
end

end