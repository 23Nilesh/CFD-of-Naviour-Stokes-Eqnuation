%% Using Central Difference Method
clc; clearvars; 

%% Defining the problem domain
n_points =101;%no_of_points
dom_length=1;
h =dom_length/(n_points-1);
x = 0:h:dom_length; %x domain span
y = 0:h:dom_length; %y domain span

% Initializing the problem
T(1,1:n_points)=1;
T(1:n_points,1)=1;
T_new(1,1:n_points)=1;
T_new(1:n_points,1)=1;

rho=1;
u=1;
v=1;
gamma=0.002;
P = rho*u*h/gamma;

%%Differencing scheme
error=1;
iterations=0;
while error>1e-7
    for i=2:n_points-1
        for j = 2:n_points-1
            a_E = gamma - rho*u*h/2;
            a_W = gamma + rho*u*h/2;
            a_N = gamma - rho*v*h/2;
            a_S = gamma + rho*v*h/2;
            a_P = rho*u*h/2 - rho*u*h/2 + rho*v*h/2 - rho*v*h/2 + gamma + gamma + gamma + gamma;
            T_new(i,j) = (a_E*T(i+1,j) + a_W*T(i-1,j) + a_N*T(i,j-1) + a_S*T(i,j+1))/a_P;
        end
    end
    iterations = iterations +1;
    error =0;
    for i = 2:n_points-1
        for j = 2:n_points-1
            error = error + abs(T(i,j)-T_new(i,j));
        end
    end
    T = T_new;
end


%% plotting
x_dom = ((1:n_points)-1).*h;
y_dom = 1-((1:n_points)-1).*h;

[X,Y] = meshgrid(x_dom,y_dom);

contourf(X,Y,T,10)
colorbar
xlabel('x')
ylabel('y')
title('u=1 v=1 Γ=0.002')

% Centerline temperature
figure;
plot(1-y,T(:,(n_points+1)/2),'--o')

title('u=1 v=1 Γ=0.002')


%% Using upwind scheme
clc; clearvars;
% Defining the problem domain
n_points =101;%no_of_points
dom_length=1;
h =dom_length/(n_points-1);
x = 0:h:dom_length; %x domain span
y = 0:h:dom_length; %y domain span

% Initializing the problem
T(1,1:n_points)=1;
T(1:n_points,1)=1;
T_new(1,1:n_points)=1;
T_new(1:n_points,1)=1;

rho=1;
u=1;
v=1;
gamma=0.01;
P = rho*u*h/gamma;

%%UPWIND
error=1;
iterations=0;
while error>1e-7
    for i=2:n_points-1
        for j = 2:n_points-1
            a_E = gamma;
            a_W = gamma + rho*u*h/2;
            a_N = gamma;
            a_S = gamma + rho*v*h/2;
            a_P = rho*u*h/2  + rho*v*h/2  + gamma + gamma + gamma + gamma;
            T_new(i,j) = (a_E*T(i+1,j) + a_W*T(i-1,j) + a_N*T(i,j-1) + a_S*T(i,j+1))/a_P;
        end
    end
    iterations = iterations +1;
    error =0;
    for i = 2:n_points-1
        for j = 2:n_points-1
            error = error + abs(T(i,j)-T_new(i,j));
        end
    end
    T = T_new;
end

%% plotting
figure;
x_dom = ((1:n_points)-1).*h;
y_dom = 1-((1:n_points)-1).*h;
[X,Y] = meshgrid(x_dom,y_dom);
contourf(X,Y,T,50)
colorbar
xlabel('x')
ylabel('y')
title('u=1 v=1 Γ=0.01');

% Centerline temperature
figure;
plot(1-y,T(:,(n_points+1)/2),'--o')
title('u=1 v=1 Γ=0.01')