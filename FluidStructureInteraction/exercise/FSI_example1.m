close
clear
clc
%% Specify domain size,# of cells in x,y and mesh resolution

xe = 4.0; % Domain length in x
ye = 2.0; % Domain length in y

nx = 200; % Number of cells in x
ny = 100; % Number of cells in y

nx_f = nx*2; % Number of fine grid cells in x
ny_f = ny*2; % Number of fien grid cells in y

dx = xe/nx; % Cell size x
dy = ye/ny; % Cell size y

dx_f = dx/2; % Fine grid cell size x
dy_f = dy/2; % Fine grid cell size y

x = repmat(dx/2:dx:xe-dx/2,ny,1); % Cell centered x location
y = repmat(dy/2:dy:ye-dy/2,nx,1)'; % Cell centered y location

xx = repmat(0:dx:xe,ny,1);
yy = repmat(0:dy:ye,nx,1);

%% Fluid properties

rhoL = 1.0; % Fluid density
mu = 0.01; % Fluid dynamic viscosity

nu = mu/rhoL; % Fluid kinematic viscosity
%% Initialize variables

u = zeros(ny+2,nx+2); % Face centered u velocity(to the left face)
v = zeros(ny+2,nx+2); % Face centered v velocity(to the bottom face)
p = zeros(ny+2,nx+2); % Cell centered pressure
div = zeros(ny+2,nx+2); % Cell centered divergence

u_old = zeros(ny+2,nx+2);
v_old = zeros(ny+2,nx+2);

rho = zeros(ny+2,nx+2);
rho_f = zeros(ny_f+4,nx_f+4);

U = zeros(ny,nx);
V = zeros(ny,nx);
velMag = zeros(ny,nx);

vorticity = zeros(ny+2,nx+2);
Vorticity = zeros(ny,nx);

initialVelocity = 0.2; % Inlet velocity

dt = min(0.25*dx*dx/nu,4.0*nu/initialVelocity/initialVelocity); % Time step(linear advection diffusion stability condition)

twfin = 100 * dt; % Stopping criteria for simulation

%% Initialize immersed or rigid body

points = 50;

radius = 0.25;
Re = rhoL*initialVelocity*2*radius/mu;

xCent = xe/2.0;
yCent = ye/2.0;

rhoPsi = 1.0;

FSI_x = zeros(ny+2,nx+2);
FSI_y = zeros(ny+2,nx+2);

[X,X_t,X_c] = initIBM(points,radius,xCent,yCent,dx,dy); % Initialize my immersed boundary

points = length(X);

%% Real grid boundary conditions

% Parallel to wall
ut = 0.0; % Top wall 
ub = 0.0; % Bottom wall
vl = 0.0; % Left wall
vr = 0.0; % Right wall

% Perpendicular to wall(Not implemented - requires a different BC
% function,where outlet velocity is solved for,and slip or no slip
% conditions are used on the other walls)
%(A possible next direction for flow over a step or two-phase flow)
ul = initialVelocity;
ur = 0.0;
vt = 0.0;
vb = 0.0;

% Gravity
gx = 0.0;
gy = 0.0;

%% Main loop

t = 0; % Start time for simulation
endTime = 10;
pldt = 0.05;
plt = 0.0;
fileNumber = 0;

while t < endTime
    [u,v,u_old,v_old] = boundaryConditions(u,v,u_old,v_old,ut,ub,ul,ur,vt,vb,vl,vr,nx,ny); % Set boundary conditions
    
    dt = stabilityCondition(u,v,nu,nx,ny,dx); % CFL condition

    [u,v,FSI_x,FSI_y,f_x,f_y] = ibmClassic(X_c,X,X_t,u,v,nx,ny,dx,dy);

    [u,v] = intermediateVelocity(u,v,u_old,v_old,FSI_x,FSI_y,mu,nx,ny,dx,dy,dt);

    [p,a_p] = pressureSolve(u,v,rhoL,nx,ny,dx,dy,dt); % Pressure iterative solver
    [u,v,u_old,v_old] = accel(p,u,v,rhoL,nx,ny,dx,dy,dt); % Advance time step

    X = updateFiberModel(X,u,v,dx,dy,dt); % Advance fiber model
    
    % Post
    if plt > pldt

        % figure
        % plot(X(:,1),X(:,2),'ro')
        % xlim([0 xe])
        % ylim([0 ye])
        % xlabel("x")
        % ylabel("y")
        % set(gcf,'Position',[200 200 1000,450])
        %exportgraphics(gca,'comboExample.gif','Append',true)

        for i=2:nx+1
            for j=2:ny+1
                vorticity(j,i) = sqrt(((v(j,i) - v(j,i-1))/dx -(u(j,i) - u(j-1,i))/dy)^2);
            end
        end

        for i=1:nx
            for j=1:ny
                Vorticity(j,i) = 0.5*(vorticity(j+1,i+1) + vorticity(j+1,i+2));
            end
        end

        figure(1)
        contourf(x,y,Vorticity,'LevelList',0:.2:3) 
        hold on
        colorbar
        set(gcf,'Position',[800 400 1000,380])
        xlabel("x")
        ylabel("y")
        clim([0 3])
         
        exportgraphics(gca,'Re_10_vortexField.gif','Append',true)

        fileNumber = fileNumber + 1;
        plt = 0;
    end
    dragForce = 0.0;
    liftForce = 0.0;
    ds = sqrt((X(2,1) - X(1,1))^2 +(X(2,2) - X(1,2))^2);
    
    for k=1:points
        dragForce = dragForce + -f_x(k)*ds;
        liftForce = liftForce + f_y(k)*ds;
    end
    % figure
    % plot(t,dragForce,'ro') 
    % hold on
    % xlim([0 5])
    % ylim([0 0.25])
    % xlabel("Time(s)")
    % ylabel('Drag force')
    % set(gcf,'Position',[900 500 1000,380])
    %exportgraphics(gca,'test1.gif','Append',true)
    
    % figure
    % plot(t,liftForce,'bo') 
    % hold on
    % xlim([0 5])
    % ylim([-1.5e-3 1.5e-3])
    % xlabel("Time(s)")
    % ylabel('Lift force')
    % set(gcf,'Position',[0 400 1000,380])
    %exportgraphics(gca,'test2.gif','Append',true)

    t = t + dt
    plt = plt + dt;
end
%% Post processing

dragForce = 0.0;
ds = sqrt((X(2,1) - X(1,1))^2 +(X(2,2) - X(1,2))^2);

for k=1:points
    dragForce = dragForce + f_x(k)*ds;
end

for i=1:nx
    for j=1:ny
        U(j,i) = 0.5 *(u(j+1,i+1) + u(j+1,i+2));
        V(j,i) = 0.5 *(v(j+1,i+1) + v(j+2,i+1));
        velMag(j,i) = sqrt(U(j,i)^2 + V(j,i)^2);
    end
end

% Calculate divergence and find max
div = calcDiv(div,u,v,nx,ny,dx,dy);
[m,ii,jj] = maxDiv(div,nx,ny);

%% Functions

function [dt] = stabilityCondition(u,v,nu,nx,ny,dx)
    maxVelMag = 0.0;
    for i=1:nx+2
        for j=1:ny+2
            velMag = sqrt(u(j,i)^2 + v(j,i)^2);
            if velMag > maxVelMag
                maxVelMag = velMag;
            end
        end
    end
    dt1 = min(0.25*dx*dx/nu,4.0*nu/maxVelMag/maxVelMag);
    dt2 = 0.3*dx/maxVelMag; % CFL stability condition

    dt = min(dt1,dt2);
end

function [u,v,u_old,v_old] = boundaryConditions(~,~,u_old,v_old,ut,ub,ul,ur,vt,vb,vl,vr,nx,ny)
    u_old(:,2) = ul; % Left wall
    u_old(:,nx+2) = u_old(:,nx+1); % Right wall
    u_old(ny+2,:) = 2*ut + u_old(ny+1,:); % Top wall
    u_old(1,:) = 2*ub + u_old(2,:); % Bottom wall

    v_old(ny+2,:) = vt; % Top wall
    v_old(2,:) = vb; % Bottom wall
    v_old(:,1) = 2.0*vl - v_old(:,2); % Left wall
    v_old(:,nx+2) = 2.0*vr - v_old(:,nx+1); % Right wall;

    u = u_old;
    v = v_old;
end

function [X,X_t,X_c] = initIBM(points,radius,xCent,yCent,dx,dy)
    X = zeros(points,2); % x,y position of set of lagangrian points
    X_t = zeros(points,2);

    X_c = zeros(points,points); % specified connectivity between springs

    dTheta = 2*pi/points;
    theta = 0.0;

    for j=1:points
        theta = theta + dTheta;

        x = radius*sin(theta);
        y = radius*cos(theta);

        X(j,1) = xCent + x;
        X(j,2) = yCent + y;
    end

    X_t = X;

    % Define connectivity matrix
    e = ones(points,1);
    X_c = spdiags([e 0*e e],-1:1,points,points);

    % If end points are connectead
    X_c(1,end) = 1;
    X_c(end,1) = 1;
end

function [f_x,f_y] = targetPoint(X,X_t,f_x,f_y)
    points = length(X(:,1));
    k_t = 100;

    for i=1:points
        f_x(i) = -k_t*(X(i,1) - X_t(i,1));
        f_y(i) = -k_t*(X(i,2) - X_t(i,2));
    end
end

function [f_x,f_y] = springsHookean(X_c,X,f_x,f_y)
    points = length(X(:,1));

    k_s = 20; % spring stiffness(N/m)
    r_l = 0.1; % spring resting length(m)

    for j=1:points
        for i=1:points
            % We are connected to this point
            if X_c(j,i) == 1
                x_r = X(i,1);
                x_l = X(j,1);
                y_r = X(i,2);
                y_l = X(j,2);

                tempF_x = k_s*(sqrt((x_r - x_l)^2 +(y_r - y_l)^2) - r_l)*(x_r - x_l);
                tempF_y = k_s*(sqrt((x_r - x_l)^2 +(y_r - y_l)^2) - r_l)*(y_r - y_l);

                f_x(j) = f_x(j) + tempF_x;
                f_y(j) = f_y(j) + tempF_y;
            end
        end
    end
end

function [f_x,f_y] = springsTorsional(X,f_x,f_y)
    points = length(X(:,1));
    
    k_b = 20; % bending stiffness(N/m)
    theta = 20.0*pi/180; % desired angle

    % the master is i+1
    for i=2:points-1

        X_l = X(i-1,1);
        X_m = X(i,1);
        X_r = X(i+1,1);

        Y_l = X(i-1,2);
        Y_m = X(i,2);
        Y_r = X(i+1,2);

        d_lm = sqrt((X_l - X_m)^2 +(Y_l - Y_m)^2);
        d_rm = sqrt((X_r - X_m)^2 +(Y_r - Y_m)^2);
        C = d_lm*d_rm*sin(theta);

        f_x(i) = k_b*((X_r - X_m)*(Y_m - Y_l) -(Y_r - Y_m)*(X_m - X_l) - C)*((Y_m - Y_l) +(Y_r - Y_m));
        f_y(i) = k_b*((X_r - X_m)*(Y_m - Y_l) -(Y_r - Y_m)*(X_m - X_l) - C)*(-(X_r - X_m) -(X_m - X_l));
    end
end

function [u,v,FSI_x,FSI_y,f_x,f_y] = ibmClassic(X_c,X,X_t,u,v,nx,ny,dx,dy)
    points = length(X(:,1));

    ds = dx;

    E_x = zeros(points,1);
    E_y = zeros(points,1);

    f_x = zeros(points,1);
    f_y = zeros(points,1);

    FSI_x = zeros(ny+2,nx+2);
    FSI_y = zeros(ny+2,nx+2);

    [f_x,f_y] = targetPoint(X,X_t,f_x,f_y);
    % [f_x,f_y] = springsHookean(X_c,X,f_x,f_y);
    % [f_x,f_y] = springsTorsional(X,f_x,f_y);

    for i=2:nx+1
        for j=2:ny+1
            for k=1:points
                tx_X = dx*(i-2);
                ty_X = dy*(j-2) + 0.5*dy;

                tx_Y = dx*(i-2) + 0.5*dx;
                ty_Y = dy*(j-2);

                r_x_X = abs(X(k,1) - tx_X) / dx;
                r_y_X = abs(X(k,2) - ty_X) / dy;

                r_x_Y = abs(X(k,1) - tx_Y) / dx;
                r_y_Y = abs(X(k,2) - ty_Y) / dy;

                if r_x_X <= 2
                    phi_x = delta1(r_x_X);
                else
                    phi_x = 0.0;
                end
                
                if r_y_X <= 2
                    phi_y = delta1(r_y_X);
                else
                    phi_y = 0.0;
                end

                delta = phi_x*phi_y/dx/dy;
                FSI_x(j,i) = FSI_x(j,i) + f_x(k)*delta*ds;

                if r_x_Y <= 2
                    phi_x = delta1(r_x_Y);
                else
                    phi_x = 0.0;
                end
                
                if r_y_Y <= 2
                    phi_y = delta1(r_y_Y);
                else
                    phi_y = 0.0;
                end

                delta = phi_x*phi_y/dx/dy;
                FSI_y(j,i) = FSI_y(j,i) + f_y(k)*delta*ds;
            end
        end
    end
end

function [d] = delta1(a)
    d = 0.25*(1 + cos(pi*a/2));
end

function [X] = updateFiberModel(X,u,v,dx,dy,dt)
    points = length(X(:,1));

    for k=1:points

        % Look in the vicinity to interpolate velocity of the fibers
        x_left = round(X(k,1)/dx - 2,'TieBreaker','minusinf');
        x_right = round(X(k,1)/dx + 2,'TieBreaker','plusinf');
        
        y_bot = round(X(k,2)/dy - 2,'TieBreaker','minusinf');
        y_top = round(X(k,2)/dy + 2,'TieBreaker','plusinf');

        u_fiber = 0.0;
        v_fiber = 0.0;

        for i=x_left:x_right
            for j=y_bot:y_top

                tx_X = dx*(i-2);
                ty_X = dy*(j-2) + 0.5*dy;

                tx_Y = dx*(i-2) + 0.5*dx;
                ty_Y = dy*(j-2);

                r_x_X = abs(X(k,1) - tx_X) / dx;
                r_y_X = abs(X(k,2) - ty_X) / dy;

                r_x_Y = abs(X(k,1) - tx_Y) / dx;
                r_y_Y = abs(X(k,2) - ty_Y) / dy;

                if r_x_X <= 2
                    phi_x = 0.25*(1 + cos(pi*r_x_X/2));
                else
                    phi_x = 0.0;
                end

                if r_y_X <= 2
                    phi_y = 0.25*(1 + cos(pi*r_y_X/2));
                else
                    phi_y = 0.0;
                end

                delta = phi_x*phi_y/dx/dy;
                u_fiber = u_fiber + u(j,i)*delta*dx*dy;

                if r_x_Y <= 2
                    phi_x = 0.25*(1 + cos(pi*r_x_Y/2));
                else
                    phi_x = 0.0;
                end

                if r_y_Y <= 2
                    phi_y = 0.25*(1 + cos(pi*r_y_Y/2));
                else
                    phi_y = 0.0;
                end

                delta = phi_x*phi_y/dx/dy;
                v_fiber = v_fiber + v(j,i)*delta*dx*dy;

            end
        end
        X(k,1) = X(k,1) + u_fiber*dt;
        X(k,2) = X(k,2) + v_fiber*dt;
    end
end

function [u,v] = intermediateVelocity(u,v,u_old,v_old,FSI_x,FSI_y,mu,nx,ny,dx,dy,dt)

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
            convection = -(u_e*u_e - u_w*u_w)/dx -(v_n*u_n - v_s*u_s)/dy;
            diffusion = mu*(u_old(j,i-1) - 2.0*u_old(j,i) + u_old(j,i+1))/dx/dx + mu*(u_old(j+1,i) - 2.0*u_old(j,i) + u_old(j-1,i))/dy/dy;
            
            % Calculate intermediate u velocity
            u(j,i) = u_old(j,i) + dt*(diffusion + convection + FSI_x(j,i));
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
            convection = -(v_e*u_e - v_w*u_w)/dx -(v_n*v_n - v_s*v_s)/dy;
            diffusion = mu*(v_old(j,i-1) - 2*v_old(j,i) + v_old(j,i+1))/dx/dx + mu*(v_old(j+1,i) - 2*v_old(j,i) + v_old(j-1,i))/dy/dy;
            
            % Calculate intermediate v velocity
            v(j,i) = v_old(j,i) + dt*(diffusion + convection + FSI_y(j,i));
        end
    end    
end

function [p,a_p] = pressureSolve(u,v,rho,nx,ny,dx,dy,dt)
    [p,a_p] = pressureOptimized2(u,v,rho,nx,ny,dx,dy,dt);
end

function [p,a_p] = pressureOptimized2(u,v,rho,nx,ny,dx,dy,dt)

    p = zeros(ny+2,nx+2);

    for i=2:nx+1
        for j=2:ny+1
            rhs(j,i) =((u(j,i+1) - u(j,i))/dx +(v(j+1,i) - v(j,i))/dy)/dt;
        end
    end

    e = ones(nx,1);
    Ax = spdiags([e 0*e e],-1:1,ny,ny)/dx/dx;
    Ay = spdiags([e 0*e e],-1:1,nx,nx)/dy/dy;

    Ix=speye(ny,ny);
    Iy=speye(nx,nx);

    Iy3 = speye(nx,nx);
    Iy3(1:nx-1,1:nx-1) = 0;

    a1 = kron(Iy,Ax); % This is correct
    a2 = kron(Ay,Ix); % this is correct

    A3 = a1 + a2;

    A4 = A3;
    
    middleCoeff = diag(sum(A4));
    A5 = A4 - middleCoeff;

    % Now add in boundary conditiosn
    Ax2 = spdiags([0*e -e 0*e ],-1:1,ny,ny)/dx/dx;
    Ax2(1,1) = -e(1)/dx/dx;
    Ax2(ny,ny) = -e(ny)/dx/dx;
    Iy2 = speye(nx,nx);
    Iy2(2:nx,2:nx) = 0;

    Ax22 = spdiags([0*e -e 0*e ],-1:1,ny,ny)/dx/dx;
    Ax22(ny,ny) = -e(ny)/dx/dx;

    a3 = kron(Iy2,Ax2)/rho;
    a4 = kron(Iy3,Ax22)/rho;

    A6 = A5;
    A6 = A5 + a3 + a4;

    b = zeros(ny,nx);
    for i=2:nx+1
        for j=2:ny+1
            b(j-1,i-1) = rhs(j,i);
        end
    end
    b2 = reshape(b,[ny*nx,1]);
    
    x2 = A6 \ b2;
    x = reshape(x2,[ny,nx]);

    for i=2:nx+1
        for j=2:ny+1
            p(j,i) = x(j-1,i-1);
        end
    end
    a_p = 1;
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
