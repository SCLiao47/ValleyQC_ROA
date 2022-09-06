function lmax = func_lambdaMax(Q)
% compute the largest eigenvalues of matrix Q
    lmax = max(eig(Q));

    if lmax <= 0
        warning('No positive eigenvalue!');
    end
end