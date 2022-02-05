clc
clear all 
close all 
%initializing dimensions for discretization in space
nx = 2001;              %number of nodes in x
ny = 81;               %number of nodes in y
L = 1;                  %length
H = 0.1;                %height
dx = L/(nx-1);  
dy = H/(ny-1);
Re = 1e4;               %Reynold's number
%initializing meshgrid
x = linspace(0,L,nx);   
y = linspace(0,H,ny);
[X,Y]=meshgrid(x,y); 
%% Initializing v and u vectors
v = zeros(ny,nx);
u = zeros(ny,nx);
%% Setting up boundary conditions
u(:,1) = 1; % u at inlet
v(:,1) = 0; % v at inlet
u(1,:) = 0; % No slip condition at the bottom
v(1,:) = 0; % No slip condition at the bottom
u(ny,:) = 1; % Free outer flow at the top

%% Numerical Solution using Finite Difference Method
for i = 1:nx-1
    %determining u(i+1,j) from the derived FD equation
    for j = 2:ny-1
        usol(j) = u(j,i) + dx/(Re*u(j,i)*(dy^2))*(u(j+1,i)-2*u(j,i) +u(j-1,i)) - dx*v(j,i)/(2*dy*u(j,i))*(u(j+1,i)-u(j-1,i));
    end
    u(2:ny-1,i+1) = usol(2:end);
    % determining v from the derived FD equation
    for j = 2:ny
        v(j,i+1) = v(j-1,i+1) - dy/2/dx*(u(j,i+1)-u(j,i)+u(j-1,i+1)-u(j-1,i));
    end
end

%% Numerical Solution of Blasius Equation Using Runge-Kutta
h = 0.05
f1 = @(y2) y2;
f2 = @(y3) y3;
f3 = @(y1, y3) -y1*y3/2;
eta = 0:h:10;
%initial conditions of f and f'. Best f'' initial guess determined by Newton Raphson Method
y1(1) = 0;
y2(1) = 0;
y3(1) = 0.33230;   
for i = 1:(length(eta)-1)
  a = h.*[f1(y2(i)), f2(y3(i)), f3(y1(i), y3(i))];
  b = h.*[f1(y2(i)+a(2)/2), f2(y3(i)+a(3)/2), f3( y1(i)+a(1)/2, y3(i)+a(3)/2)];
  c = h.*[f1(y2(i)+b(2)/2), f2(y3(i)+b(3)/2), f3( y1(i)+b(1)/2, y3(i)+b(3)/2)];
  d = h.*[f1(y2(i)+c(2)), f2(y3(i)+c(3)), f3(y1(i)+c(1), y3(i)+c(3))];
  y3(i+1) = y3(i)+ 1/6*(a(3)+2*b(3)+2*c(3)+d(3));
  y2(i+1) = y2(i)+ 1/6*(a(2)+2*b(2)+2*c(2)+d(2));
  y1(i+1) = y1(i)+ 1/6*(a(1)+2*b(1)+2*c(1)+d(1));
end


%% Plotting
%contour of u*
figure,
h1 = subplot(321);
set(h1,'XLim',[0 L],'YLim',[0 H]);
contourf(h1,X,Y,u, [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99]);
colorbar('peer',h1,'SouthOutside');
title(h1,'u* Velocity Contour');
xlabel(h1,'x'),ylabel(h1,'y');
%contour of v*
h2 = subplot(322);
set(h2,'XLim',[0 L],'YLim',[0 H]);
contourf(h2,X,Y,v);
colorbar('peer',h2,'SouthOutside');
title(h2,'v* Velocity Contour');
xlabel(h2,'x'),ylabel(h2,'y ');

%u velocity profile FD at x = 0.5 and x = 0.0005
h3 = subplot(323);
%delta and eta at x = 0.5
delta_middle = 0.01 / sqrt(2);
eta_middle = Y(:,1001)/delta_middle;
%delta and eta at x = 0.0005
delta_lead = 0.01 / sqrt(2000);
eta_lead = Y(:,2)/delta_lead;
plot(eta_middle, u(:,1001),'-k', eta_lead, u(:,2),'-r');
xlim([0 8]);
grid on;
title(h3,'u Velocity Profile FD Method');
xlabel(h3,'{\eta}'),ylabel(h3,'u*[-]');

% v velocity profile FD at x = 0.5 and x = 0.0005
h4 = subplot(324);
plot(eta_middle, v(:,1001)/(delta_middle * 2),'-k', eta_lead, v(:,2)/(delta_lead * 2000),'-r');
xlim([0 8]);
grid on;
title(h4,'v Velocity Profile FD Method');
xlabel(h4,'{\eta}'),ylabel(h4,'v*[-]');

% u velocity profile Blasius RK4
h5 = subplot(325);
u_blasius = y2;
plot( eta, u_blasius, '-k');
xlim([0 8]);
ylim([0 1]);
grid on;
title(h5,'u Velocity Profile Blasius RK4');
xlabel(h5,'{\eta}'),ylabel(h5,'u*[-]');
% v velocity profile Blasius RK4
h6 = subplot(326);
v_blasius = ( eta .* y2 - y1 )/2;
plot( eta, v_blasius, '-k');
xlim([0 8]);
grid on;
title(h6,'v Velocity Profile Blasius RK4');
xlabel(h6,'{\eta}'),ylabel(h6,'v*[-]');