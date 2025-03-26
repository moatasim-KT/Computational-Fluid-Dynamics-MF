%% Initialize variables

xe = 1.0; % Domain length in x
ye = 1.0; % Domain length in y

r = 0.15; % Sphere radius

dx = 1/200; %2*r/CPD; % Cell size x
dy = 1/200; %2*r/CPD; % Cell size y

nx = xe/dx; % Number of cells in x
ny = ye/dx; % Number of cells in y

x = repmat(dx/2:dx:xe-dx/2,ny,1); % Cell centered x location
y = repmat(dy/2:dy:ye-dy/2,nx,1)'; % Cell centered y location

u = zeros(ny+2,nx+2); % Face centered u velocity(to the left face)
v = zeros(ny+2,nx+2); % Face centered v velocity(to the bottom face)
%% Translation

% Initialize a fluid drop
xCent = 0.25;
yCent = 0.25;

f = zeros(ny+2,nx+2); % fluid volume fraction
f = initDrop(f,r,xCent,yCent,nx,ny,dx,dy);
f_init = f;

% Initialize the fluid velocity field
maxVelocity = 1.0;
u(:,:) = maxVelocity;
v(:,:) = maxVelocity;

dt = abs(0.125*dx/sqrt(2)); % CFL condition
pldt = 100;

plt = 0;
t = 0;
t_steps = round(0.5/abs(maxVelocity)/dt);

f1 = figure;
while(t < t_steps)
    
    if mod(t,2) == 0
        direction = 'x';
    else
        direction = 'y';
    end

    [f,norm,flux] = SLIC(f,u,v,direction,nx,ny,dx,dy,dt);

    % Plotting
    if plt == pldt || t == int16(t_steps-1)
        myFig = plotFigure(f);
        pause(0.1);
        plt = 0;
    end
    t = t+1;
    plt = plt+1;
end

% Post processing

massError_Translation = sum(f_init - f, "all") / sum(f_init, "all");

xCent = 0.75;
yCent = 0.75;
f_end = zeros(ny+2,nx+2);
f_end = initDrop(f_end,r,xCent,yCent,nx,ny,dx,dy);
hold on
cl = [0.5 0.5];
contour(f_end,cl,'Color','red','LineWidth',2.0)
title('F field - Translation')
colorbar
legend('SLIC advected','Exact')

%% Rotation

% Initialize a fluid drop
xCent = 0.5;
yCent = 0.5;

f = zeros(ny+2,nx+2); % fluid volume fraction
f = initDrop2(f,r,xCent,yCent,nx,ny,dx,dy);
f_init = f;

% Initialize the fluid velocity field
u(:,:) = 0;
v(:,:) = 0;
omega = 1.0;

xCent = 0.5; % Center of angular velocity omega
yCent = 0.5;

for i=2:nx+1
    for j=2:ny+1
        %Specifying u
        tx =(i-2)*dx;
        ty =(j-2)*dy + dy/2;

        xl = tx-xCent;
        yl = ty-yCent;

        u(j,i) = omega*yl;

        % Specifying v
        tx =(i-2)*dx + dx/2;
        ty =(j-2)*dy;

        xl = tx-xCent;
        yl = ty-yCent;

        v(j,i) = -omega*xl;
    end
end

dt = abs(0.125*dx/omega); % CFL condition
pldt = 750;

plt = 0;
t = 0;
t_steps = 2*pi/dt;

f2 = figure;
while t < t_steps
    
    if mod(t,2) == 0
        direction = 'x';
    else
        direction = 'y';
    end

    [f,norm,flux] = SLIC(f,u,v,direction,nx,ny,dx,dy,dt);

    % Plotting
    if plt == pldt || t == int16(t_steps-1)
        myFig = plotFigure(f);
        pause(0.1);
        plt = 0;
    end
    t = t+1;
    plt = plt+1;
end

% Post processing

massError_Rotation = sum(f_init - f, "all") / sum(f_init, "all");

hold on
contour(f_init,cl,'Color','red','LineWidth',2.0)
title('f-field Rotation')
colorbar
legend('SLIC advected','Exact')

%% Shear flow(vortex test)

% Initialize a fluid drop
xCent = 0.5;
yCent = 0.75;

f = zeros(ny+2,nx+2); % fluid volume fraction
f = initDrop(f,r,xCent,yCent,nx,ny,dx,dy);
f_init = f;

% Initialize the fluid velocity field
u(:,:) = 0;
v(:,:) = 0;
omega = 1.0;

for i=2:nx+1
    for j=2:ny+1
        %Specifying u
        tx =(i-2)*dx;
        ty =(j-2)*dy + dy/2;

        u(j,i) =(sin(pi*tx)^2)*(sin(2*pi*ty));

        % Specifying v
        tx =(i-2)*dx + dx/2;
        ty =(j-2)*dy;

        v(j,i) =(-sin(pi*ty)^2)*(sin(2*pi*tx));
    end
end

dt = abs(0.1*dx/1.0); % CFL condition
pldt = 200;

plt = 0;
t = 0;
t_end = 120*dt;
t_steps = 2560;

f3 = figure;
while t < t_steps
    
    if t == 1280
        u(:,:) = -u(:,:);
        v(:,:) = -v(:,:);
    end

    if mod(t,2) == 0
        direction = 'x';
    else
        direction = 'y';
    end
    
    [f,norm,flux,volume,v_tilde] = SLIC(f,u,v,direction,nx,ny,dx,dy,dt);

    % Plotting
    if plt == pldt || t == int16(t_steps-1)
        myFig = plotFigure(f);
        pause(0.1);
        plt = 0;
    end
    t = t+1;
    plt = plt+1;
end

% Post processing

percentVolumeLoss_Vortex = 100*sum(f_init - f, "all") / sum(f_init, "all");

div = zeros(ny,nx);
div = calcDiv(div,u,v,nx,ny,dx,dy);

xCent = 0.5;
yCent = 0.75;
f_end = zeros(ny+2,nx+2);
f_end = initDrop(f_end,r,xCent,yCent,nx,ny,dx,dy);
hold on
contour(f_end,cl,'Color','red','LineWidth',2.0)
title('f-field vortex')
colorbar
legend('SLIC advected','Exact')
%% Initialize a fluid drop

function f = initDrop(f,r,xCent,yCent,nx,ny,dx,dy)
    % Generate volume fraction field
    for i=2:nx+1
        for j=2:ny+1
            tx =(i-2)*dx; % x value(left)
            ty =(j-2)*dy; % y value(bottom)
            txp =(i-2)*dx + dx; % x value(right)
            typ =(j-2)*dy + dy; % y value(top)
    
            if max((tx-xCent)^2,(txp-xCent)^2) +  max((ty-yCent)^2,(typ-yCent)^2) <= r^2
                f(j,i) = 1.0;
            elseif min((tx-xCent)^2,(txp-xCent)^2) +  min((ty-yCent)^2,(typ-yCent)^2) > r*r
                f(j,i) = 0.0;
            else
                for l=1:40
                    for m=1:40
                        txx =(i-2)*dx +(l-1)*dx/40;
                        tyy =(j-2)*dy +(m-1)*dy/40;
                        txxp =(i-2)*dx + l*dx/40;
                        tyyp =(j-2)*dy + m*dy/40;
    
                        if max((txx-xCent)^2,(txxp-xCent)^2) + max((tyy-yCent)^2,(tyyp-yCent)^2) < r^2
                            f(j,i) = f(j,i) +(1/40)*(1/40);
                        end
                    end
                end
            end
        end
    end
end

function f = initDrop2(f,r,xCent,yCent,nx,ny,dx,dy)
    % Generate volume fraction field
    for i=2:nx+1
        for j=2:ny+1
            tx =(i-2)*dx; % x value(left)
            ty =(j-2)*dy; % y value(bottom)
            txp =(i-2)*dx + dx; % x value(right)
            typ =(j-2)*dy + dy; % y value(top)
            
            fTemp = 0.0;
            if max((tx-xCent)^2,(txp-xCent)^2) +  max((ty-yCent)^2,(typ-yCent)^2) <= r^2
                fTemp = 1.0;
            elseif min((tx-xCent)^2,(txp-xCent)^2) +  min((ty-yCent)^2,(typ-yCent)^2) > r*r
                fTemp = 0.0;
            else
                for l=1:40
                    for m=1:40
                        txx =(i-2)*dx +(l-1)*dx/40;
                        tyy =(j-2)*dy +(m-1)*dy/40;
                        txxp =(i-2)*dx + l*dx/40;
                        tyyp =(j-2)*dy + m*dy/40;
    
                        if max((txx-xCent)^2,(txxp-xCent)^2) + max((tyy-yCent)^2,(tyyp-yCent)^2) < r^2
                            fTemp = fTemp +(1/40)*(1/40);
                        end
                    end
                end
            end
            if ty < 0.5 && tx > 0.5 - 6*dx && tx < 0.5 + 6*dx
                f(j,i) = 0.0;
            else
                f(j,i) = fTemp;
            end
        end
    end
end
%% SLIC function
% Split advection scheme,which reconstructs the interface based on a
% 3-cell stencil parallel to the direction of advection.

function [f,norm,flux,volume,v_tilde] = SLIC(f,u,v,direction,nx,ny,dx,dy,dt)
    
    volume = ones(ny+2,nx+2)*dx*dy;
    v_tilde = zeros(ny+2,nx+2);
    counter = 0;
    while(counter < 2)
        norm = zeros(ny+2,nx+2);
        flux = zeros(ny+2,nx+2);
        f_old = f;
        
        switch(direction)
            case 'x'
                % Calculate norm and v_tilde
                for i=2:nx+1
                    for j=2:nx+1
                        if f(j,i) == 1.0 || f(j,i) == 0.0
                            norm(j,i) = 0;
                        elseif f(j,i+1) < f(j,i-1)
                            norm(j,i) = 1;
                        elseif f(j,i+1) > f(j,i-1)
                            norm(j,i) = -1;
                        else
                            norm(j,i) = 0;
                        end
                        
                        v_tilde(j,i) = volume(j,i) - dt*dy*(u(j,i+1) - u(j,i));
                    end
                end
                
                % Calculate flux
                for i=2:nx
                    for j=2:ny
                        if norm(j,i) == 1 || norm(j,i) == 0
                            if(u(j,i) < 0)
                                if abs(u(j,i)*dt) >= f(j,i)*dx
                                    flux(j,i) = -f(j,i)*dx;
                                else
                                    flux(j,i) = u(j,i)*dt;
                                end
                            end
                            if u(j,i+1) > 0
                                if f(j,i)*dx + u(j,i+1)*dt > dx
                                    flux(j,i+1) = f(j,i)*dx + u(j,i+1)*dt - dx;
                                else
                                    flux(j,i+1) = 0.0;
                                end
                            end
                        elseif norm(j,i) == -1
                            if u(j,i+1) > 0
                                if f(j,i)*dx <= u(j,i+1)*dt
                                    flux(j,i+1) = f(j,i)*dx;
                                else
                                    flux(j,i+1) = u(j,i+1)*dt;
                                end
                            end
                            if u(j,i) < 0
                                if f(j,i)*dx + abs(u(j,i)*dt) > dx
                                    flux(j,i) = dx - f(j,i)*dx - abs(u(j,i)*dt);
                                else
                                    flux(j,i) = 0.0;
                                end
                            end
                        end
                    end
                end
                
                % Update f
                for i=2:nx
                    for j=2:ny
                        f(j,i) = f_old(j,i)*volume(j,i) +(flux(j,i) - flux(j,i+1))/dx*volume(j,i);
                        f(j,i) = f(j,i)/v_tilde(j,i);
                    end
                end
                volume = v_tilde;

            case 'y'
                % Calculate norm
                for i=2:nx+1
                    for j=2:nx+1
                        if f(j,i) == 1.0 || f(j,i) == 0.0
                            norm(j,i) = 0;
                        elseif f(j+1,i) < f(j-1,i)
                            norm(j,i) = 1;
                        elseif f(j+1,i) > f(j-1,i)
                            norm(j,i) = -1;
                        else
                            norm(j,i) = 0;
                        end

                        v_tilde(j,i) = volume(j,i) - dt*dx*(v(j+1,i) - v(j,i));
                    end
                end
                
                % Calculate flux
                for i=2:nx
                    for j=2:ny
                        if norm(j,i) == 1 || norm(j,i) == 0
                            if v(j,i) < 0
                                if abs(v(j,i)*dt) >= f(j,i)*dy
                                    flux(j,i) = -f(j,i)*dx;
                                else
                                    flux(j,i) = v(j,i)*dt;
                                end
                            end
                            if v(j+1,i) > 0
                                if f(j,i)*dy + v(j+1,i)*dt > dy
                                    flux(j+1,i) = f(j,i)*dy + v(j+1,i)*dt - dy;
                                else
                                    flux(j+1,i) = 0.0;
                                end
                            end
                        elseif norm(j,i) == -1
                            if v(j+1,i) > 0
                                if(v(j+1,i)*dt >= f(j,i)*dy)
                                    flux(j+1,i) = f(j,i)*dy;
                                else
                                    flux(j+1,i) = v(j+1,i)*dt;
                                end
                            end
                            if v(j,i) < 0
                                if f(j,i)*dy + abs(v(j,i)*dt) > dy
                                    flux(j,i) = dy - f(j,i)*dy - abs(v(j,i)*dt);
                                else
                                    flux(j,i) = 0.0;
                                end
                            end
                        end
                    end
                end
                
                % Update f
                for i=2:nx
                    for j=2:ny
                        f(j,i) = f_old(j,i)*volume(j,i) +(flux(j,i) - flux(j+1,i))/dy*volume(j,i);
                        f(j,i) = f(j,i)/v_tilde(j,i);
                    end
                end
                volume = v_tilde;
        end

        counter = counter + 1;
        if direction == 'x'
            direction = 'y';
        else
            direction = 'x';
        end
    end
end
%% Utility functions
function div = calcDiv(div,u,v,nx,ny,dx,dy)
    for i=2:nx-1
        for j=2:ny-1
            div(j,i) =(u(j,i+1) - u(j,i))/dx +(v(j+1,i) - v(j,i))/dy;
        end
    end
end

function myFigure = plotFigure(f)
    cl = [0.5 0.5];
    myFigure = contourf(f,cl);
end
