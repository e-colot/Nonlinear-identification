function H = memoryPolynomial(u, Norder, Nlag)
% H = memoryPolynomial(u, Norder, Nlag)
% Builds the regression matrix for a memory polynomial
%
% Written by E. COlot on Dec 28 2025

    H = zeros(length(u), Nlag*Norder); % Initialize regression matrix
    index = 1;
    
    for q = 0:Nlag
        % Create delayed version of u (pad with zeros at the start)
        u_delayed = [zeros(q, 1); u(1:end-q)];
        
        for k = 1:Norder
            % Use odd-only orders if desired (change to k = 1:2:Norder)
            % Standard MPM definition uses: u(n-q) * |u(n-q)|^(k-1)
            basis_vec = u_delayed .* (abs(u_delayed).^(k-1));
            
            H(:, index) = basis_vec;
            index = index + 1;
        end
    end

end
