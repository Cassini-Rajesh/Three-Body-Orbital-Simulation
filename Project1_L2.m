%% Project 1

% Load data from file

L2 = load("EM_L2-304P1.mat");

MU1 = L2.MU1;
T = L2.T;
tspan = linspace(0, T, 1000);
perturbation = L2.perturbation;
x0 = L2.x0;
x0 = x0(:);

%% ode45

options = odeset('RelTol',1e-9, 'AbsTol',1e-12);
[t, sol] = ode45(@(t,y) L2SS(t, y, MU1), tspan, x0, options);

%% Plot rotating frame

x_pos = sol(:,1);
y_pos = sol(:,2);

figure;
plot(x_pos, y_pos)
xlabel('x position')
ylabel('y position')
title('Lyapunov Orbit Around L2 (rotating)')
grid on;

%% Plot inertial frame

Xn_pos = sol(:,1);
Yn_pos = sol(:,2);

omega = 1;

Xn_pos = cos(omega*t) .* x_pos - sin(omega*t) .* y_pos;
Yn_pos = sin(omega*t) .* x_pos + cos(omega*t) .* y_pos;

figure;
plot(Xn_pos, Yn_pos)
xlabel('x position')
ylabel('y position')
title('Lyapunov Orbit Around L2 (inertial)')
grid on;

%% Perturbation

x0_perturbed = x0 + perturbation;

options = odeset('RelTol',1e-9, 'AbsTol',1e-12);
[t, sol_p] = ode45(@(t,y) L2SS(t, y, MU1), tspan, x0_perturbed, options);

delta_x1 = sol_p(:,1) - sol(:,1);
delta_x2 = sol_p(:,2) - sol(:,2);
delta_x3 = sol_p(:,3) - sol(:,3);
delta_x4 = sol_p(:,4) - sol(:,4);

figure;
plot(t, delta_x1, 'b');
hold on;
plot(t, delta_x2, 'r');
xlabel('Time');
ylabel('Position Deviation');
title('Position vs Time');
legend('\deltax', '\deltay');
grid on;

figure;
plot(t, delta_x3, 'b');
hold on;
plot(t, delta_x4, 'r');
xlabel('Time');
ylabel('Velocity Deviation');
title('Velocity vs Time');
legend('\deltaVx', '\deltaVy');
grid on;

%% Linearization

function jacob = LL2(t, delta_x, A)

    jacob = A * delta_x;
end 

x1 = sol(1);
x2 = sol(2);

p1 = sqrt( (x1 + MU1)^2 + x2^2 );
p2 = sqrt( (x1 + MU1 - 1)^2 + x2^2 );

Uxx = 1 - ( (1 - MU1) * ( (p1^2 - 3*(x1 + MU1)^2) / p1^5 ) ) ...
       - ( MU1 * ( (p2^2 - 3*(x1 + MU1 - 1)^2) / p2^5 ) );

Uyy = 1 - ( (1 - MU1) * ( (p1^2 - 3*x2^2) / p1^5 ) ) ...
       - ( MU1 * ( (p2^2 - 3*x2^2) / p2^5 ) );

Uxy = 3 * (1 - MU1) * ( (x1 + MU1) * x2 / p1^5 ) ...
       + 3 * MU1 * ( (x1 + MU1 - 1) * x2 / p2^5 );

A = [0 0 1 0;
     0 0 0 1;
   Uxx Uxy 0 2;
   Uxy Uyy -2 0];

[t, sol_L] = ode45(@(t, delta_x) LL2(t, delta_x, A), tspan, x0_perturbed', options);

delta_x1L = sol_L(:,1);
delta_x2L = sol_L(:,2);
delta_x3L = sol_L(:,3);
delta_x4L = sol_L(:,4);

m1 = sqrt( delta_x1.^2 + delta_x2.^2);
m2 = sqrt( delta_x1L.^2 + delta_x2L.^2);

figure;
plot(t, m1, 'b');
hold on;
plot(t, m2, 'r');
xlabel('Time')
ylabel('Magnitude of Position');
title('Magnitude of Position vs Time');
legend('Non-Linear', 'Linearized');
grid on;

































