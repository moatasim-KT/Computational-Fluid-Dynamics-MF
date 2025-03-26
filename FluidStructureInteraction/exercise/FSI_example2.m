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

initialVelocity = 4; % Inlet velocity

dt = min(0.25*dx*dx/nu,4.0*nu/initialVelocity/initialVelocity); % Time step(linear advection diffusion stability condition)

twfin = 100 * dt; % Stopping criteria for simulation

%% Initialize rigid body

radius = 0.25;
Re = rhoL*initialVelocity*2*radius/mu;

xCent = xe/2.0;
yCent = ye/2.0;

rhoPsi = 1.0;

psi = zeros(ny+2,nx+2);
psi_f = zeros(ny_f+4,nx_f+4);

FSI_x = zeros(ny+2,nx+2);
FSI_y = zeros(ny+2,nx+2);

psiU = zeros(ny+2,nx+2);
psiV = zeros(ny+2,nx+2);

u_RBM = 0.0;
v_RBM = 0.0;

[psi_f,totalVolume] = initPsi_f(psi_f,radius,xCent,yCent,nx_f,ny_f,dx_f,dy_f);
[u_old,v_old] = initailizePsiVelocity(psi_f,u,v,u_RBM,v_RBM,nx,ny);

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
pldt = 0.1;
plt = 0.0;
fileNumber = 1;

pp = nsidedpoly(1000,'Center',[xCent yCent],'Radius',radius);

while t < endTime
    [u,v,u_old,v_old] = boundaryConditions(u,v,u_old,v_old,ut,ub,ul,ur,vt,vb,vl,vr,nx,ny); % Set boundary conditions
    
    dt = stabilityCondition(u,v,nu,nx,ny,dx); % CFL condition

    % [psi_f,psi_old_x_f,psi_old_y_f,psiX_f,psiY_f,xCent,yCent,u_old,v_old,FSI_x,FSI_y,u_RBM,v_RBM] = conservationMomentum(FSI_x,FSI_y,psi_f,rho_f,rhoPsi,rhoL,xCent,yCent,radius,u_RBM,v_RBM,u,v,nx,ny,nx_f,ny_f,dx,dy,dx_f,dy_f,dt,t);
    [psi_f,psi_old_x_f,psi_old_y_f,psiX_f,psiY_f,xCent,yCent,u_old,v_old,FSI_x,FSI_y,u_RBM,v_RBM] = conservationMomentumStationary(FSI_x,FSI_y,psi_f,rho_f,rhoPsi,rhoL,xCent,yCent,radius,u_RBM,v_RBM,u,v,nx,ny,nx_f,ny_f,dx,dy,dx_f,dy_f,dt,t);

    [rho_f] = setProperties(psi_f,rho_f,rhoPsi,rhoL,nx_f,ny_f);
    [u,v,rhoU,rhoV] = intermediateVelocity(u,v,u_old,v_old,psi_old_x_f,psi_old_y_f,psiX_f,psiY_f,rho_f,rhoPsi,rhoL,mu,nx,ny,dx,dy,dt,t); % Solve for intermediate velocity condition

    [p,a_p] = pressureSolve(u,v,rhoL,rhoU,rhoV,nx,ny,dx,dy,dt,t); % Pressure iterative solver
    [u,v,u_old,v_old] = accel(p,u,v,rhoL,nx,ny,dx,dy,dt); % Advance time step
    
    % Post
    if plt > pldt

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
        contourf(x,y,Vorticity,'LevelList',0:2:45)
        hold on
        plot(pp,'FaceColor','red')
        colorbar
        set(gcf,'Position',[200 200 1200,400])
        clim([0 45])

        exportgraphics(gca,'oscilating_exampleTTt.gif','Append',true)
 
        plt = 0;
    
    end
    fileNumber = fileNumber + 1;
    [d(fileNumber),dragForce(fileNumber),pressureForce(fileNumber),fsi(fileNumber)] = dragCalculation(FSI_x,p,psi_f,u,v,mu,rhoL,radius,initialVelocity,nx,ny,dx,dy);

    figure(2)
    plot(t,d(fileNumber),'ro')
    hold on
    xlim([0 5])
    xlabel("Time(s)")
    ylabel('Drag force')
    set(gcf,'Position',[0 400 1000,380])
    exportgraphics(gca,'FSI_example22.gif','Append',true)

    t = t + dt
    plt = plt + dt;
end

%% Post processing

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

function [dragCoeffNum,dragForce,pressureForce,fsi] = dragCalculation(FSI_x,p,psi_f,u,v,mu,rhoL,radius,initialVelocity,nx,ny,dx ,dy)

    pressureForce = 0.0;
    dragForce1 = 0.0;
    dragForce2 = 0.0;
    convForce = 0.0;
    totalA = 0.0;

    fsi = 0.0;

    for i=2:nx-1
        for j=2:ny-1
            psi = 0.25*(psi_f(2*j,2*i-1) + psi_f(2*j-1,2*i-1) + psi_f(2*j,2*i-2) + psi_f(2*j-1,2*i-2));
            if psi > 0
                 psi = 1.0;
            end

            dragForce1 = dragForce1 + psi*(u(j,i-1) - 2*u(j,i) + u(j,i+1));
            dragForce2 = dragForce2 + psi*(u(j-1,i) - 2*u(j,i) + u(j+1,i));
            
            pressureForce = pressureForce + psi*(p(j,i) - p(j,i-1))*dx;
    
            u_e = 0.5*(u(j,i) + u(j,i+1));
            u_w = 0.5*(u(j,i) + u(j,i-1));

            u_n = 0.5*(u(j,i) + u(j+1,i));
            u_s = 0.5*(u(j,i) + u(j-1,i));
            
            v_n = 0.5*(v(j+1,i-1) + v(j+1,i));
            v_s = 0.5*(v(j,i-1) + v(j,i));
            convForce = convForce + psi*(u_e*u_e - u_w*u_w)/dx + psi*(u_n*v_n - u_s*v_s);
            
            if psi > 0 && psi < 1
                fsi = fsi + FSI_x(j,i);
            end
        end
    end
    dragForce = dragForce1 + dragForce2;
    
    dragForce_T = dragForce*mu;
    pressureForce2 = pressureForce;
    convForce2 = 0.0;
    
    totalForce = -convForce2 - pressureForce2 + dragForce_T;
    
    A =(0.5*2*radius*rhoL*initialVelocity^2);
    dragCoeffNum =  totalForce;
    fsi = fsi*dx*dx/A;
end

function dt = stabilityCondition(u,v,nu,nx,ny,dx)
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

function [u,v,u_old,v_old] = boundaryConditions(u,v,u_old,v_old,ut,ub,ul,ur,vt,vb,vl,vr,nx,ny)
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

function [rho_f] = setProperties(psi_f,rho_f,rhoPsi,rhoL,nx_f,ny_f)
    for i=1:nx_f+4
        for j=1:ny_f+4
            rho_f(j,i) = rhoPsi*psi_f(j,i) + rhoL*(1 - psi_f(j,i));
        end
    end
end

function [psi_f,totalVolume] = initPsi_f(psi_f,r,xCent,yCent,nx_f,ny_f,dx_f,dy_f)

    totalVolume = 0.0;
    % Generate volume fraction field
    for i=3:nx_f+2
        for j=3:ny_f+2
            tx =(i-2)*dx_f; % x value(left)
            ty =(j-2)*dy_f; % y value(bottom)
            txp =(i-2)*dx_f + dx_f; % x value(right)
            typ =(j-2)*dy_f + dy_f; % y value(top)
    
            if max((tx-xCent)^2,(txp-xCent)^2) +  max((ty-yCent)^2,(typ-yCent)^2) < r^2
                ftemp = 1.0;
            elseif min((tx-xCent)^2,(txp-xCent)^2) +  min((ty-yCent)^2,(typ-yCent)^2) > r*r
                ftemp = 0.0;
            else
                ftemp = 0.0;
                for l=1:40
                    for m=1:40
                        txx =(i-2)*dx_f +(l-1)*dx_f/40;
                        tyy =(j-2)*dy_f +(m-1)*dy_f/40;
                        txxp =(i-2)*dx_f + l*dx_f/40;
                        tyyp =(j-2)*dy_f + m*dy_f/40;
    
                        if max((txx-xCent)^2,(txxp-xCent)^2) + max((tyy-yCent)^2,(tyyp-yCent)^2) < r^2
                            ftemp = ftemp +(1/40)*(1/40);
                        end
                    end
                end
            end
            psi_f(j,i) = psi_f(j,i) + ftemp;
            totalVolume = totalVolume + psi_f(j,i);
        end
    end
end

function [u,v] = initailizePsiVelocity(psi_f,u,v,u_RBM,v_RBM,nx,ny)
    for i=2:nx+1
        for j=2:ny+1
            
            if psi_f(2*j,2*i-1) + psi_f(2*j,2*i-2) + psi_f(2*j-1,2*i-1) + psi_f(2*j-1,2*i-2) > 0
                u(j,i) = u_RBM;
            end

            if psi_f(2*j-1,2*i) + psi_f(2*j-1,2*i-1) + psi_f(2*j-2,2*i) + psi_f(2*j-2,2*i-1) > 0
                v(j,i) = v_RBM;
            end
        end
    end
end

function [psi_f,rho_f,xCent,yCent,totalVolume] = updatePsi(psi_f,rho_f,rhoPsi,rhoL,u_RBM,v_RBM,xCent,yCent,r,nx_f,ny_f,dx_f,dy_f,dt,direction)

    psi_f(:,:) = 0;
    if direction == 'x'
        xCent = xCent + u_RBM*dt;
    else
        yCent = yCent + v_RBM*dt;
    end

    totalVolume = 0.0;
    % Generate volume fraction field
    for i=3:nx_f+2
        for j=3:ny_f+2
            tx =(i-2)*dx_f; % x value(left)
            ty =(j-2)*dy_f; % y value(bottom)
            txp =(i-2)*dx_f + dx_f; % x value(right)
            typ =(j-2)*dy_f + dy_f; % y value(top)
    
            if max((tx-xCent)^2,(txp-xCent)^2) +  max((ty-yCent)^2,(typ-yCent)^2) < r^2
                ftemp = 1.0;
            elseif min((tx-xCent)^2,(txp-xCent)^2) +  min((ty-yCent)^2,(typ-yCent)^2) > r*r
                ftemp = 0.0;
            else
                ftemp = 0.0;
                dx = 20;
                dy = 20;
                for l=1:dx
                    for m=1:dy
                        txx =(i-2)*dx_f +(l-1)*dx_f/dx;
                        tyy =(j-2)*dy_f +(m-1)*dy_f/dx;
                        txxp =(i-2)*dx_f + l*dx_f/dy;
                        tyyp =(j-2)*dy_f + m*dy_f/dy;
    
                        if max((txx-xCent)^2,(txxp-xCent)^2) + max((tyy-yCent)^2,(tyyp-yCent)^2) < r^2
                            ftemp = ftemp +(1/dx)*(1/dy);
                        end
                    end
                end
            end
            psi_f(j,i) = psi_f(j,i) + ftemp;
            totalVolume = totalVolume + psi_f(j,i)*dx_f*dy_f;
        end
    end

    [rho_f] = setProperties(psi_f,rho_f,rhoPsi,rhoL,nx_f,ny_f);
end

function [psi_f,psi_old_x_f,psi_old_y_f,psiX_f,psiY_f,xCent,yCent,u,v,FSI_x,FSI_y,u_RBM,v_RBM] = conservationMomentum(FSI_x,FSI_y,psi_f,rho_f,rhoPsi,rhoL,xCent,yCent,radius,u_RBM,v_RBM,u,v,nx,ny,nx_f,ny_f,dx,dy,dx_f,dy_f,dt,t)

    psiX = zeros(ny+2,nx+2);
    psiY = zeros(ny+2,nx+2);

    psiX_f = zeros(ny_f+4,nx_f+4);
    psiY_f = zeros(ny_f+4,nx_f+4);

    FSI_x(:,:) = 0.0;
    FSI_y(:,:) = 0.0;

    if mod(t,2) == 0
        direction = 'x';
    else
        direction = 'y';
    end

    counter = 0;
    while(counter < 2)
        
        switch(direction)
            case 'x'
                psi_old_x_f = psi_f;
                [psi_f,rho_f,xCent,yCent,totalVolume] = updatePsi(psi_f,rho_f,rhoPsi,rhoL,u_RBM,v_RBM,xCent,yCent,radius,nx_f,ny_f,dx_f,dy_f,dt,direction);
                psiX_f = psi_f;

                psi_mass = totalVolume*rhoPsi;

                for i=2:nx+1
                    for j=2:ny+1
                        psiX(j,i) = 0.25*(psi_f(2*j,2*i-1) + psi_f(2*j-1,2*i-1) + psi_f(2*j,2*i-2) + psi_f(2*j-1,2*i-2));
                    end
                end

                rhoU = 0.0;
                for i=2:nx+1
                    for j=2:ny+1
                        if psiX(j,i) > 0.0
                            rhoU = rhoU + psiX(j,i)*u(j,i)*rhoPsi*dx*dy;
                        end
                    end
                end

                u_RBM = rhoU/psi_mass;
                u_RBM = 0.0;

                for i=2:nx+1
                    for j=2:ny+1

                        if psiX(j,i) > 0.0
                            FSI_x(j,i) = rhoPsi*(u_RBM - u(j,i))/dt;
                            u(j,i) = u_RBM;
                        end
                    end
                end

                direction = 'y';

            case 'y'
                psi_old_y_f = psi_f;
                [psi_f,rho_f,xCent,yCent,totalVolume] = updatePsi(psi_f,rho_f,rhoPsi,rhoL,u_RBM,v_RBM,xCent,yCent,radius,nx_f,ny_f,dx_f,dy_f,dt,direction);
                psiY_f = psi_f;

                psi_mass = totalVolume*rhoPsi;

                for i=2:nx+1
                    for j=2:ny+1
                        psiY(j,i) = 0.25*(psi_f(2*j-1,2*i) + psi_f(2*j-1,2*i-1) + psi_f(2*j-2,2*i) + psi_f(2*j-2,2*i-1));
                    end
                end

                rhoV = 0.0;
                for i=2:nx+1
                    for j=2:ny+1

                        if psiY(j,i) > 0.0
                            rhoV = rhoV + psiY(j,i)*v(j,i)*rhoPsi*dx*dy;
                        end
                    end
                end

                v_RBM = rhoV/psi_mass;
                %v_RBM = 0.0;

                for i=2:nx+1
                    for j=2:ny+1

                        if psiY(j,i) > 0
                            FSI_y(j,i) = rhoPsi*(v_RBM - v(j,i))/dt;
                            v(j,i) = v_RBM;
                        end
                    end
                end
                direction = 'x';
        end
        counter = counter + 1;
    end
end

function [psi_f,psi_old_x_f,psi_old_y_f,psiX_f,psiY_f,xCent,yCent,u,v,FSI_x,FSI_y,u_RBM,v_RBM] = conservationMomentumStationary(FSI_x,FSI_y,psi_f,rho_f,rhoPsi,rhoL,xCent,yCent,radius,u_RBM,v_RBM,u,v,nx,ny,nx_f,ny_f,dx,dy,dx_f,dy_f,dt,t)
    
    psiX = zeros(ny+2,nx+2);
    psiY = zeros(ny+2,nx+2);

    psi_old_x_f = psi_f;
    psi_old_y_f = psi_f;
    psiX_f = psi_f;
    psiY_f = psi_f;

    u_RBM = 0.0;
    
    for i=2:nx+1
        for j=2:ny+1
            psiX(j,i) = 0.25*(psi_f(2*j,2*i-1) + psi_f(2*j-1,2*i-1) + psi_f(2*j,2*i-2) + psi_f(2*j-1,2*i-2));
        end
    end

    for i=2:nx+1
        for j=2:ny+1
    
            if psiX(j,i) > 0.0
                u(j,i) = u_RBM;
            end
        end
    end
    
    v_RBM = 0.0;

    for i=2:nx+1
        for j=2:ny+1
            psiY(j,i) = 0.25*(psi_f(2*j-1,2*i) + psi_f(2*j-1,2*i-1) + psi_f(2*j-2,2*i) + psi_f(2*j-2,2*i-1));
        end
    end

    for i=2:nx+1
        for j=2:ny+1
    
            if psiY(j,i) > 0
                FSI_y(j,i) = rhoPsi*(v_RBM - v(j,i))/dt;
                v(j,i) = v_RBM;
            end
        end
    end
end

function [u,v,rhoU,rhoV] = intermediateVelocity(u,v,u_old,v_old,psi_old_x_f,psi_old_y_f,psiX_f,psiY_f,rho_f,rhoPsi,rhoL,mu,nx,ny,dx,dy,dt,t)

    rhoU_x = ones(ny+2,nx+2);
    rhoU_y = ones(ny+2,nx+2);
    rhoV_x = ones(ny+2,nx+2);
    rhoV_y = ones(ny+2,nx+2);

    rhoU = ones(ny+2,nx+2);
    rhoV = ones(ny+2,nx+2);

    for i=2:nx+1
        for j=2:ny+1
            psiStag = 0.25*(psiX_f(2*j,2*i-1) + psiX_f(2*j,2*i-2) + psiX_f(2*j-1,2*i-1) + psiX_f(2*j-1,2*i-2));
            rhoU_x(j,i) = rhoPsi*psiStag + rhoL*(1-psiStag);

            psiStag = 0.25*(psiY_f(2*j-1,2*i-1) + psiY_f(2*j-1,2*i-2) + psiY_f(2*j-2,2*i-1) + psiY_f(2*j-2,2*i-2));
            rhoU_y(j,i) = rhoPsi*psiStag + rhoL*(1-psiStag);
        end
    end

    for i=3:nx+1
        for j=2:ny+1
            % Interpolating velocities
            u_e = 0.5*(u_old(j,i) + u_old(j,i+1));
            u_w = 0.5*(u_old(j,i-1) + u_old(j,i));
            u_n = 0.5*(u_old(j,i) + u_old(j+1,i));
            u_s = 0.5*(u_old(j,i) + u_old(j-1,i));
            
            v_n = 0.5*(v_old(j+1,i-1) + v_old(j+1,i));
            v_s = 0.5*(v_old(j,i-1) + v_old(j,i));

            psiFace = 0.25*(psi_old_x_f(2*j,2*i-1) + psi_old_x_f(2*j,2*i-2) + psi_old_x_f(2*j-1,2*i-1) + psi_old_x_f(2*j-1,2*i-2));
            rhoFace = psiFace*rhoPsi + rhoL*(1 - psiFace);
        
            rhoU(j,i) = 0.25*(rho_f(2*j,2*i-1) + rho_f(2*j,2*i-2) + rho_f(2*j-1,2*i-1) + rho_f(2*j-1,2*i-2));

            % Solving div(rho*u*u) and div(tau) 
            convection = -rhoFace*(u_e*u_e - u_w*u_w)/dx - rhoFace*(v_n*u_n - v_s*u_s)/dy; % Simplification by sharma and patankar
            diffusion = mu*(u_old(j,i-1) - 2.0*u_old(j,i) + u_old(j,i+1))/dx/dx + mu*(u_old(j+1,i) - 2.0*u_old(j,i) + u_old(j-1,i))/dy/dy;
            
            % Calculate intermediate u velocity
            u(j,i) = rhoFace*u_old(j,i) + dt*(diffusion + convection);
            u(j,i) = u(j,i)/rhoU(j,i);
        end
    end

    for i=2:nx+1
        for j=2:ny+1
            psiStag = 0.25*(psiX_f(2*j-1,2*i-1) + psiX_f(2*j-2,2*i-2) + psiX_f(2*j-2,2*i-1) + psiX_f(2*j-2,2*i-2));
            rhoV_x(j,i) = rhoPsi*psiStag +(1 - psiStag)*rhoL;

            psiStag = 0.25*(psiY_f(2*j-1,2*i) + psiY_f(2*j-1,2*i-1) + psiY_f(2*j-2,2*i) + psiY_f(2*j-2,2*i-2));
            rhoV_y(j,i) = rhoPsi*psiStag +(1 - psiStag)*rhoL;
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

            psiFace2 = 0.25*(psi_old_y_f(2*j-1,2*i) + psi_old_y_f(2*j-1,2*i-1) + psi_old_y_f(2*j-2,2*i) + psi_old_y_f(2*j-2,2*i-2));
            rhoFace2 = psiFace2*rhoPsi +(1 - psiFace2)*rhoL;

            rhoV(j,i) = 0.25*(rho_f(2*j-1,2*i) + rho_f(2*j-1,2*i-1) + rho_f(2*j-2,2*i) + rho_f(2*j-2,2*i-2));

            % Solving div(rho*u*u) and div(tau)
            convection = -rhoFace2*(v_e*u_e - v_w*u_w)/dx - rhoFace2*(v_n*v_n - v_s*v_s)/dy;
            diffusion = mu*(v_old(j,i-1) - 2*v_old(j,i) + v_old(j,i+1))/dx/dx + mu*(v_old(j+1,i) - 2*v_old(j,i) + v_old(j-1,i))/dy/dy;
            
            % Calculate intermediate v velocity
            v(j,i) = rhoFace2*v_old(j,i) + dt*(diffusion + convection);
            v(j,i) = v(j,i)/rhoV(j,i);
        end
    end
    % FSI can be set implicitly,by just enforcing the rigid body velocity,
    % No need to add source term,also avoids confusion of rho discrepency
    
end

function [p,a_p] = pressureSolve(u,v,rho,rhoU,rhoV,nx,ny,dx,dy,dt,iter)
    [p,a_p] = pressureOptimized2(u,v,rho,rhoU,rhoV,nx,ny,dx,dy,dt,iter);
end

function [p,a_p] = pressureOptimized2(u,v,rho,rhoU,rhoV,nx,ny,dx,dy,dt,iter)
  
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

    rhoU_new = zeros(ny,nx);
    rhoV_new = zeros(ny,nx);
    
    for i=1:nx
        for j=1:ny
            rhoU_new(j,i) = 1/rhoU(j+1,i+1);
            rhoV_new(j,i) = 1/rhoV(j+1,i+1);
        end
    end

    % We can comment this out,if density of fluid=solid
    b2 = reshape(rhoV_new,[ny*nx,1]);

    c2 = reshape(rhoU_new,[ny*nx,1]);

    b2 = sparse(b2);
    c2 = sparse(c2);
    d1 = sparse(diag(b2(1:ny*nx-1,:),1));
    d2 = sparse(diag(b2(2:ny*nx,:),-1));
    d3 = sparse(diag(c2(1:ny*nx-ny,:),ny));
    d4 = sparse(diag(c2(1+ny:ny*nx,:),-ny));

    dv = d1 + d2 + d3 + d4;
    dv = sparse(dv);

    A4 = A3.*dv;
    % A4 = A3;
    
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
