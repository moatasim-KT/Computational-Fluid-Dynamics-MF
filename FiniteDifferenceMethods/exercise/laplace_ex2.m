% Domain setup
L = 1;
H = 1;

% Points in the domain
nx = 22;
ny = 22;

% Grid resolution
dx = L/nx;
dy = H/ny;

% Initializing solution
U_init = zeros(ny,nx);

% Boundary conditions
U_init(:,1) = 10;
U_init(:,nx) = 5;
U_init(1,:) = 4;
U_init(ny,:) = 1;

U_sol = U_init;

% Build the coefficient matrix
% Use MATLAB function spdiag() and feed in matrix coefficients
n = nx-2;
e = ones(n,1);
T = ;

% Kronecker Delta product
% Use MATLAb function kron, and take the prudct of I(x)T + T(x)I
A = ;

% Build the right hand side
% This matrix will contain your boundary conditions and is equal to
% the negative value of your stencil.
for j=2:n+1
    for i=2:n+1
        rhs(j-1,i-1) = ;
    end
end

% Reshape the rhs to a 1D vector
% Use MATLAB reshape() to [n*n, 1]


% Solve the equation x = A \ rhs


% Reshape the matrix back to a n by n grid


% Add back in the boundary conditions
for j=2:n+1
    for i=2:n+1
        
    end
end

%% Plot the solution
figure
contourf(U_sol)
title("Laplace equation solution")
xlabel("x")
ylabel("y")
colorbar