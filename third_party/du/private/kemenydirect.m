% Copyright (c) 2023, Fabio Durastante.
% Reproduced from Kemeny-and-Conquer:
% https://github.com/Cirdans-Home/Kemeny-and-Conquer
% Licensed under the BSD 3-Clause License; see
% ../LICENSE-Kemeny-and-Conquer.
function k = kemenydirect(P)
%%KEMENYDIRECT Computes the Kemeny constant by matrix inversion

n = size(P,1);
I = speye(n);
e = ones(n,1);
k = trace(inv(I - P + e*e'/n))-1;

end
