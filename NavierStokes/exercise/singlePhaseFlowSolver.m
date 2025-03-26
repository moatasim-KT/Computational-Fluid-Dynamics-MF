%% Specify domain size,# of cells in x,y and mesh resolution

xe = 1.0; % Domain length in x
ye = 1.0; % Domain length in y

nx = 8; % Number of cells in x
ny = 8; % Number of cells in y

dx = xe/nx; % Cell size x
dy = ye/ny; % Cell size y

x = repmat(dx/2:dx:xe-dx/2,ny,1); % Cell centered x location
y = repmat(dy/2:dy:ye-dy/2,nx,1)'; % Cell centered y location

%% Fluid properties

rho = 1.0; % Fluid density
mu = 0.01; % Fluid dynamic viscosity

nu = mu/rho; % Fluid kinematic viscosity

%% Initialize variables

u = zeros(ny+2,nx+2); % Face centered u velocity(to the left face)
v = zeros(ny+2,nx+2); % Face centered v velocity(to the bottom face)
p = zeros(ny+2,nx+2); % Cell centered pressure
div = zeros(ny+2,nx+2); % Cell centered divergence

u_old = zeros(ny+2,nx+2);
v_old = zeros(ny+2,nx+2);

initialVelocity = 1.0; % Top wall u velocity
Re = rho*initialVelocity*xe/mu; % Reynolds number

dt = min(0.25*dx*dx/nu,4.0*nu/initialVelocity/initialVelocity); % Time step(linear advection diffusion stability condition)

twfin = 100 * dt; % Stopping criteria for simulation

%% Real grid boundary conditions

% Parallel to wall
ut = initialVelocity; % Top wall 
ub = 0.0; % Bottom wall
vl = 0.0; % Left wall
vr = 0.0; % Right wall

% Perpendicular to wall(Not implemented - requires a different BC
% function,where outlet velocity is solved for,and slip or no slip
% conditions are used on the other walls)
%(A possible next direction for flow over a step or two-phase flow)
ul = 0.0;
ur = 0.0;
vt = 0.0;
vb = 0.0;

%% Open data files

fileID = fopen('pressureSolverConvergenceData.txt','w');
fprintf(fileID,'%8s       %8s             %8s\n','Timestep','Residual','#-of-iterations');
fclose(fileID);
%% Main loop

t = 0; % Start time for simulation
realTime = 0;
endTime = 1;
iter = endTime/dt;
while t < iter
    [u,v,u_old,v_old] = boundaryConditions(u,v,u_old,v_old,ut,ub,ul,ur,vt,vb,vl,vr,nx,ny); % Set boundary conditions
    
    [dt] = stabilityCondition(u,v,nx,ny,dx);

    [u,v] = intermediateVelocity(u,v,u_old,v_old,rho,mu,nx,ny,dx,dy,dt); % Solve for intermediate velocity condition
    [p,a_p] = pressureSolve(u,v,rho,nx,ny,dx,dy,dt,t); % Pressure iterative solver
    [u,v,u_old,v_old] = accel(p,u,v,rho,nx,ny,dx,dy,dt); % Advance time step
    
    t = t + 1;
    realTime = realTime + dt;
end
%% Post processing

U = zeros(ny,nx);
V = zeros(ny,nx);
velMag = zeros(ny,nx);

% Interpolate my cell centered U and V velocites
for i=1:nx
    for j=1:ny
        U(j,i) = 0.5 *(u(j+1,i+1) + u(j+1,i+2));
        V(j,i) = 0.5 *(v(j+1,i+1) + v(j+2,i+1));
        velMag(j,i) = sqrt(U(j,i)^2 + V(j,i)^2);
    end
end

% Calculate divergence and find max
[div] = calcDiv(div,u,v,nx,ny,dx,dy);
[m,ii,jj] = maxDiv(div,nx,ny);

% Plotting velocity vectors
figure
subplot(1,2,1)
quiver(x,y,U,V)
title('Velocity vector data')
xlabel('X-coordinate(m)')
ylabel('Y-coordinate(m)')

% Plotting stream particles
subplot(1,2,2)
[verts,averts] = streamslice(x,y,U,V);
streamline([verts averts])
title('Stream Particles')
xlabel('X-coordinate(m)')
ylabel('Y-coordinate(m)')

set(gcf,'Position',[100 100 1000 350])

%% Comparison at Re = 100
A = load('/data/Re_100.csv');
B = load('/data/Re_100y.csv');

Uq = interp2(x,y,U,0.5,y(:,1));
Vq = interp2(x,y,V,x(1,:),0.5);

figure
subplot(1,2,1)
plot(A(:,1),A(:,2),'LineWidth',1)
hold on
plot(y(:,1),Uq/initialVelocity,'o','MarkerSize',7)
title('U velocity at x=0.5,for Re=100')
xlabel('Y-coordinate(m)')
ylabel('U(m/s)')
legend('Ghia et al.',"Our simulation")

subplot(1,2,2)
plot(B(:,1),B(:,2),'LineWidth',1) 
hold on
plot(x(1,:),Vq/initialVelocity,'o','MarkerSize',7)
title('V velocity at y=0.5,for Re=100')
xlabel('X-coordinate(m)')
ylabel('V(m/s)')
legend('Ghia et al.',"Our simulation")

set(gcf,'Position',[100 100 1000 350])
%% Write your own code here!













%% Functions

function dt = stabilityCondition(u,v,nx,ny,dx)
    maxVelMag = 0.0;
    for i=1:nx+2
        for j=1:ny+2
            velMag = sqrt(u(j,i)^2 + v(j,i)^2);
            if velMag > maxVelMag
                maxVelMag = velMag;
            end
        end
    end

    dt = 0.3*dx/maxVelMag; % CFL stability condition
end

function [u,v,u_old,v_old] = boundaryConditions(u,v,u_old,v_old,ut,ub,ul,ur,vt,vb,vl,vr,nx,ny)
    u_old(:,2) = ul; % Left wall
    u_old(:,nx+2) = ur; % Right wall
    u_old(ny+2,:) = 2*ut - u_old(ny+1,:); % Top wall
    u_old(1,:) = 2*ub - u_old(2,:); % Bottom wall

    v_old(ny+2,:) = vt; % Top wall
    v_old(2,:) = vb; % Bottom wall
    v_old(:,1) = 2.0*vl - v_old(:,2); % Left wall
    v_old(:,nx+2) = 2.0*vr - v_old(:,nx+1); % Right wall

    u = u_old;
    v = v_old;
end

function [u,v] = intermediateVelocity(u,v,u_old,v_old,rho,mu,nx,ny,dx,dy,dt)
    for i=3:nx+1
        for j=2:ny+1
            % Interpolating velocities
            u_e = 0.5*(u_old(j,i) + u_old(j,i+1));
            u_w = 0.5*(u_old(j,i-1) + u_old(j,i));
            u_n = 0.5*(u_old(j,i) + u_old(j+1,i));
            u_s = 0.5*(u_old(j,i) + u_old(j-1,i));
            
            v_n = 0.5*(v_old(j+1,i-1) + v_old(j+1,i));
            v_s = 0.5*(v_old(j,i-1) + v_old(j,i));
            
            % Solving div(rho*u*u) and div(tau) 
            convection = -(rho*u_e*u_e - rho*u_w*u_w)/dx -(rho*v_n*u_n - rho*v_s*u_s)/dy;
            diffusion = mu*(u_old(j,i-1) - 2.0*u_old(j,i) + u_old(j,i+1))/dx/dx + mu*(u_old(j+1,i) - 2.0*u_old(j,i) + u_old(j-1,i))/dy/dy;
            
            % Calculate intermediate u velocity
            u(j,i) = rho*u_old(j,i) + dt*(diffusion + convection);
            u(j,i) = u(j,i)/rho;
        end
    end

    for i=2:nx+1
        for j=3:ny+1
            % Interpolating velocities
            v_e = 0.5*(v_old(j,i) + v_old(j,i+1));
            v_w = 0.5*(v_old(j,i) + v_old(j,i-1));
            v_n = 0.5*(v_old(j,i) + v_old(j+1,i));
            v_s = 0.5*(v_old(j,i) + v_old(j-1,i));

            u_e = 0.5*(u_old(j-1,i+1) + u_old(j,i+1));
            u_w = 0.5*(u_old(j-1,i) + u_old(j,i));

            % Solving div(rho*u*u) and div(tau)
            convection = -(rho*v_e*u_e - rho*v_w*u_w)/dx -(rho*v_n*v_n - rho*v_s*v_s)/dy;
            diffusion = mu*(v_old(j,i-1) - 2*v_old(j,i) + v_old(j,i+1))/dx/dx + mu*(v_old(j+1,i) - 2*v_old(j,i) + v_old(j-1,i))/dy/dy;
            
            % Calculate intermediate v velocity
            v(j,i) = rho*v_old(j,i) + dt*(diffusion + convection);
            v(j,i) = v(j,i)/rho;
        end
    end
end

function [p,a_p] = pressureSolve(u,v,rho,nx,ny,dx,dy,dt,iter)
    [p,a_p] = sor(u,v,rho,nx,ny,dx,dy,dt,iter);
end

function [p,a_p] = sor(u,v,rho,nx,ny,dx,dy,dt,iter)
    p = zeros(ny+2,nx+2);
    rhs = zeros(ny+2,nx+2);

    % Setting pressure coefficients
    a_e = ones(ny+2,nx+2)/rho/dx/dx;
    a_w = ones(ny+2,nx+2)/rho/dx/dx;
    a_n = ones(ny+2,nx+2)/rho/dy/dy;
    a_s = ones(ny+2,nx+2)/rho/dy/dy;

    a_e(:,nx+1) = 0.0;
    a_w(:,2) = 0.0;
    a_n(nx+1,:) = 0.0;
    a_s(2,:) = 0.0;

    a_p = -(a_e + a_w + a_n + a_s);

    maxError = Inf;
    TOL = 1e-2;
    counter = 0;
    beta = 2/(1 +(pi/nx));

    while abs(maxError) > TOL
        counter = counter + 1;
        maxError = 0;
        
        % Solve for pressure field
        for i=2:nx+1
            for j=2:ny+1
                rhs(j,i) =(u(j,i+1) - u(j,i))/dx +(v(j+1,i) - v(j,i))/dy;
                rhs(j,i) = rhs(j,i)/dt -(a_w(j,i)*p(j,i-1) + a_e(j,i)*p(j,i+1) + a_n(j,i)*p(j+1,i) + a_s(j,i)*p(j-1,i));
                p(j,i) = beta*rhs(j,i)/a_p(j,i) +(1-beta)*p(j,i);
            end
        end
        
        % Determine residuals
        for i=2:nx+1
            for j=2:ny+1
                rhs(j,i) =(u(j,i+1) - u(j,i))/dx +(v(j+1,i) - v(j,i))/dy;
                error = a_w(j,i)*p(j,i-1) + a_e(j,i)*p(j,i+1) + a_n(j,i)*p(j+1,i) + a_s(j,i)*p(j-1,i) + a_p(j,i)*p(j,i) - rhs(j,i)/dt;
                if abs(error) > maxError
                    maxError = abs(error);
                end
            end
        end
    end
    
    % Write convergence data
    fileID = fopen('pressureSolverConvergenceData.txt','a+');
    fprintf(fileID,'%4d            %10d          %4d\n',iter,maxError,counter);
    fclose(fileID);
end

function [u,v,u_old,v_old] = accel(p,u,v,rho,nx,ny,dx,dy,dt)
    for i=3:nx+1
        for j=2:ny+1
            % Update u velocity
            u(j,i) = rho*u(j,i) - dt*(p(j,i) - p(j,i-1))/dx;

            u(j,i) = u(j,i) / rho;
        end
    end

    for i=2:nx+1
        for j=3:ny+1
            % Update v velocity
            v(j,i) = rho*v(j,i) - dt*(p(j,i) - p(j-1,i))/dy;

            v(j,i) = v(j,i) / rho;
        end
    end
    u_old = u;
    v_old = v;
end
%% Utility functions

function div = calcDiv(div,u,v,nx,ny,dx,dy)
    for i=2:nx-1
        for j=2:ny-1
            div(j,i) =(u(j,i+1) - u(j,i))/dx +(v(j+1,i) - v(j,i))/dy;
        end
    end
end

function [max,ii,jj] = maxDiv(div,nx,ny)
    max = 0.0;
    for i=2:nx-1
        for j=2:ny-1
            if div(j,i) > max
                max = div(j,i);
                ii = i;
                jj = j;
            end
        end
    end
end
