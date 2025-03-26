% Domain setup
L = 2;
H = 1;

% Points in the domain
nx = 44;
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
n_x = nx-2;
n_y = ny-2;

ex = ones(n_x,1);
ey = ones(n_y,1);

Tx = spdiags([ex -2*ex ex],-1:1,n_x,n_x);
Ty = spdiags([ey -2*ey ey],-1:1,n_y,n_y);

% Kronecker Delta product
A = kron(speye(n_x),Ty) + kron(Tx,speye(n_y));

% Build the right hand side
for j=2:n_y+1
    for i=2:n_x+1
        rhs(j-1,i-1) = -(U_init(j+1,i) + U_init(j-1,i) + U_init(j,i+1) + U_init(j,i-1));
    end
end

% Reshape the rhs to a 1D vector 
rhs = reshape(rhs,[n_x*n_y,1]);

% Solve the equation x = A \ rhs
x = A\rhs;

% Reshape the matrix back to a n by n grid
x = reshape(x,[n_y,n_x]);

% Add back in the boundary conditions
for j=2:n_y+1
    for i=2:n_x+1
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