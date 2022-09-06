function [A,B] = func_genRandQuadDyn(nx)
    %% Options
    eig_base = 10;       % span size of eigenvalues as power(eig_base, ~)
    digit_A = 0;          % truncate the system to digit
    digit_B = 0;
    
    %% System dimensions
    nz = nx*(nx+1)/2;
    
    %% generate A
    % this will only generate real, negative eigenvalues
    eigVals = -power(eig_base, randn(nx, 1));
    
    T = orth(randn(nx));
    
    A = T\diag(eigVals)*T;
    
    % truncate the digit to k

    A = floor(A*power(10,digit_A)) / power(10,digit_A);
    
    % check again A is Hurwitz
    e_max = eigs(A, 1, 'largestabs');
    if e_max > 0
        A = A - 2*e_max*eye(nx);
    end    
    
    %% generate B
    % sample full B 
%     B = sqrt(eig_base) * randn(nx,nz);
% %     B = power(eig_base, randn(nx,nz));
%     B = floor(B*power(10,digit_B)) / power(10,digit_B);
    

    % sample B such that b'w is more easier to be rank-3
    B = zeros(nx,nz);
    for i = 1:nx
        bi = sqrt(eig_base)*randn(1,2);
            
        idx = randsample(nz,2);
        B(i,idx) = bi;
    end

    B = floor(B*power(10,digit_B)) / power(10,digit_B);
end