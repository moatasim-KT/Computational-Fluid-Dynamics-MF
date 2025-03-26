% Points in the domain
nx = 50;
ny = 50;

% Initializing solution
U_init = zeros(ny,nx);

% Boundary conditions
U_init(:,1) = 0;
U_init(:,nx) = 0;
U_init(1,:) = 0;
U_init(ny,:) = 0;

% Corners
U_init(ny,nx) = 0.5*(U_init(ny-1,nx) + U_init(ny,nx-1));
U_init(ny,1) = 0.5*(U_init(ny-1,1) + U_init(ny,2));
U_init(1,nx) = 0.5*(U_init(2,nx) + U_init(1,nx-1));
U_init(1,1) = 0.5*(U_init(1,2) + U_init(2,1));

U_original = U_init;
U_new = U_init;

% Optimization parameter
omega = 2/(1 + (pi/nx));

% Initialize error and iteration count
maxError = 1;
iterations = 0;

% Running 3 cases,1 for each method
numCases = 0;

% Start big loop
while numCases < 3
    while maxError > 10^-3
    
        for j=2:ny-1
            for i=2:nx-1
                % Jacobi method
                if numCases == 0
                    U_new(j,i) = 0.25*(U_init(j+1,i) + U_init(j-1,i) + U_init(j,i+1) + U_init(j,i-1));
                end

                % Gauss-Seidel
                if numCases == 1
                    U_new(j,i) = 0.25*(U_new(j+1,i) + U_new(j-1,i) + U_new(j,i+1) + U_new(j,i-1));
                end

                % SOR
                if numCases == 2
                    U_new(j,i) = (1-omega)*U_init(j,i) + (0.25*omega)*(U_new(j+1,i) + U_new(j-1,i) + U_new(j,i+1) + U_new(j,i-1));
                end
            end
        end
        % Calculate the error
        
        % Stopping criteria 1
        error = U_init - U_new;
        maxError = max(max(abs(error)));

        % Stopping criteria 2
        % error = 0;
        % for j=2:ny-1
        %     for i=2:nx-1
        %         error = error + 0.25*(U_new(j+1,i) + U_new(j-1,i) + U_new(j,i+1) + U_new(j,i-1)) - U_new(j,i);
        %     end
        % end
        % maxError = abs(error);
        
        % Set U for next iteration
        U_init = U_new;
    
        % Keep track of each iteration
        iterations = iterations + 1;
    end
    caseIterations(numCases+1) = iterations;
    numCases = numCases + 1;

    U_solution = U_new;

    % Reset conditions
    U_init = U_original;
    U_new = U_init;
    maxError = 1;
    iterations = 0;
end

% Plot the solution
figure
contourf(U_solution)
title("Laplace equation solution")
xlabel("x")
ylabel("y")
colorbar