function k = kemenydirect(P)
%%KEMENYDIRECT Computes the Kemeny constant by matrix inversion

n = size(P,1);
I = speye(n);
e = ones(n,1);
k = trace(inv(I - P + e*e'/n))-1;

end