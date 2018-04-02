%% ME 598 MPF CP2 Task2 %%

% By Cong Lin & Sushobhan Sen

clear
close all

%% Definitions %%
Lx = 0.08;        % Domain Length in x-Direction: 8 cm
Ly = 0.08;        % Domain Length in y-Direction: 8 cm
Nx = 161;      
Ny = 161;
dx = Lx/(Nx-1); 
dy = Ly/(Ny-1);

x = 0:dx:Lx;
xc = 0.045;
y = 0:dy:Ly;
yc = 0.045;
[X,Y] = meshgrid(x,y);

r = 0.005;          % Bubble radius 0.5 cm
xi = 1/r;           % curvature
rho_g = 1;          % kg/m^3
rho_l = 1000;       
sigma = 0.72;
mu_g = 1e-5;
mu_l = 1e-3;

u = zeros(Ny+2,Nx+2);
v = zeros(Ny+2,Nx+2);
uhat = zeros(Ny+2,Nx+2);
vhat = zeros(Ny+2,Nx+2);
p = zeros(Ny+2,Nx+2);

% Define IC for Phi %
phi = sqrt((X-xc).^2 + (Y-yc).^2)-r;
[c,h] = contour(X*100,Y*100,phi); 
clabel(c,h);
axis equal; 
title('\phi'); xlabel('x (cm)'); ylabel('y (cm)');
%% Compute Properties %%
epsilon = 3*dx;
H = zeros(size(phi));
H(phi < -epsilon) = 0;
H(phi > epsilon) = 1;
H(abs(phi) <= epsilon) = 0.5*(1 + phi(abs(phi) <= epsilon)/epsilon + 1/pi * sin(pi*phi(abs(phi) <= epsilon)/epsilon));
figure
contour(X*100,Y*100,H); axis equal;
title('H_\epsilon'); xlabel('x (cm)'); ylabel('y (cm)');

rho_tmp = rho_g + (rho_l - rho_g)*H;            % 2D density field
figure
contour(X*100,Y*100,rho_tmp); axis equal; colorbar;
title('Density \rho'); xlabel('x (cm)'); ylabel('y (cm)');

mu = mu_g + (mu_l - mu_g)*H;                % 2D viscosity field

delta = zeros(size(phi));
delta(phi < -epsilon) = 0;
delta(phi > epsilon) = 0;
delta(abs(phi) <= epsilon) = 0.5*(1/epsilon + 1/epsilon * cos(pi*phi(abs(phi) <= epsilon)/epsilon));
figure
contour(X*100,Y*100,delta); axis equal; colorbar;
title('\delta_\epsilon (\phi)'); xlabel('x (cm)'); ylabel('y (cm)');

Dx = 1/(dx^2)*(diag(-2*ones(Nx,1)) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1));                  % x second derivative Matrix
Dx = sparse(Dx);
Dy = 1/(dy^2)*(diag(-2*ones(Ny,1)) + diag(ones(Ny-1,1),1) + diag(ones(Ny-1,1),-1));  % y second derivative Matrix
Dy = sparse(Dy);

Cx = -1*1/(2*dx)*(diag(-1*ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1));                                       % x first derivative Matrix
Cx = sparse(Cx);
Cx(1,2) = 0; Cx(1,1) = 1; Cx(Nx,Nx-1) = 0; Cx(Nx,Nx) = 1;       % collocated zero gradient
Cy = -1*1/(2*dy)*(diag(-1*ones(Ny-1,1),1) + diag(ones(Ny-1,1),-1));                       % y first derivative Matrix
Cy = sparse(Cy);
Cy(1,2) = 0; Cy(1,1) = 1; Cy(Ny,Ny-1) = 0; Cy(Ny,Ny) = 1;       % collocated zero gradient

gradphi_x = phi*Cx';
gradphi_y = Cy*phi;
figure
contourf(X*100,Y*100,gradphi_x); axis equal; colorbar;
title('\nabla \phi _x'); xlabel('x (cm)'); ylabel('y (cm)');
figure
contourf(X*100,Y*100,gradphi_y); axis equal; colorbar;
title('\nabla \phi _y'); xlabel('x (cm)'); ylabel('y (cm)');

CSF_x = -sigma*xi*delta.*gradphi_x;
CSF_y = -sigma*xi*delta.*gradphi_y;


%% Time step %%

% compute hat velocities % 
dt = 0.01;
uhat_tmp = dt*(1./rho_tmp).*CSF_x;
vhat_tmp = dt*(1./rho_tmp).*CSF_y;
uhat(2:Ny+1,2:Nx+1) = uhat_tmp;
vhat(2:Ny+1,2:Nx+1) = vhat_tmp;

% BC routine %
uhat(1,:) = uhat(2,:);          %west BC
uhat(Ny+2,:) = uhat(Ny+1,:);    %east BC
uhat(:,1) = uhat(:,2);
uhat(:,Nx+2) = uhat(:,Nx+1);

% BC routine %
vhat(1,:) = vhat(2,:);          %west BC
vhat(Ny+2,:) = vhat(Ny+1,:);    %east BC
vhat(:,1) = vhat(:,2);
vhat(:,Nx+2) = vhat(:,Nx+1);

rho = zeros(Ny+2,Nx+2);
rho(2:Ny+1,2:Nx+1) = rho_tmp;

% BC routine %
rho(1,:) = rho(2,:);          %west BC
rho(Ny+2,:) = rho(Ny+1,:);    %east BC
rho(:,1) = rho(:,2);
rho(:,Nx+2) = rho(:,Nx+1);

omega = 0.7;
ph = p;
MaxIt = 2000;
MaxErr = 0.0001;
    for it = 1:MaxIt        
        for ii = 2:Nx+1
            for jj = 2:Ny+1                 
                
                ae = 2/(rho(jj,ii+1) + rho(jj,ii))/dx^2;
                aw = 2/(rho(jj,ii) + rho(jj,ii-1))/dx^2;     
                an = 2/(rho(jj+1,ii) + rho(jj,ii))/dy^2;
                as = 2/(rho(jj,ii) + rho(jj-1,ii))/dy^2;
                
                ap = ae + aw + an + as;
                
                source = 1/dt * ( (uhat(jj,ii+1)-uhat(jj,ii-1))/(2*dx) + (vhat(jj+1,ii)-vhat(jj-1,ii))/(2*dy) );
                term = ae*p(jj,ii+1) + aw*p(jj,ii-1) + an*p(jj+1,ii) + as*p(jj-1,ii) - source;
                term = term/ap;
                
                ph(jj,ii) = omega * term + (1-omega) * p(jj,ii);
            end
        end
        
        % BC routine %
        ph(1,:) = ph(2,:);          %west BC
        ph(Ny+2,:) = ph(Ny+1,:);    %east BC
        ph(:,1) = ph(:,2);
        ph(:,Nx+2) = ph(:,Nx+1);
        
        err = 0;
        err = sumabs(ph-p);
        if err <= MaxErr 
            break 
        end    
        p = ph;                     %update p
    end

figure
contourf(X*100,Y*100,p(2:Ny+1,2:Nx+1)); axis equal; colorbar;
title('Pressure Field (Pa)'); xlabel('x (cm)'); ylabel('y (cm)');

u_new = -dt*(1./rho_tmp).*(p(2:Ny+1,2:Nx+1)*Cx') + uhat_tmp;
v_new = -dt*(1./rho_tmp).*(Cy*p(2:Ny+1,2:Nx+1)) + vhat_tmp;

figure
contourf(X*100,Y*100,u_new); axis equal; colorbar;
title('u (m/s)'); xlabel('x (cm)'); ylabel('y (cm)');

figure
contourf(X*100,Y*100,v_new); axis equal; colorbar;
title('v (m/s)'); xlabel('x (cm)'); ylabel('y (cm)');