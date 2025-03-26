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
n = nx-2;
e = ones(n,1);
T = spdiags([e -2*e e],-1:1,n,n);

% Kronecker Delta product
A = kron(speye(n),T) + kron(T,speye(n));

% Build the right hand side
for j=2:n+1
    for i=2:n+1
        rhs(j-1,i-1) = -(U_init(j+1,i) + U_init(j-1,i) + U_init(j,i+1) + U_init(j,i-1));
    end
end

% Reshape the rhs to a 1D vector 
rhs = reshape(rhs,[n*n,1]);

% Solve the equation x = A \ rhs
x = A\rhs;

% Reshape the matrix back to a n by n grid
x = reshape(x,[n,n]);

% Add back in the boundary conditions
for j=2:n+1
    for i=2:n+1
        U_sol(j,i) = x(j-1,i-1); 
    end
end

%% Plot the solution
figure
contourf(U_sol)
title("Laplace equation solution")
xlabel("x")
ylabel("y")
colorbar