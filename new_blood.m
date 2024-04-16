% Time settings
dt = 1; % time step are in minutes
time = 0:dt:1440; % simulation time from 0 to 24 hours

% Parameters
k_abs = 0.05; % absorption rate constant
Km = 80; % Michaelis constant for glucose absorption
ku = 0.01; % glucose utilization rate constant
ki = 0.005; % insulin response rate
kd = 0.1; % insulin decay rate
kp = 0.03; % liver response rate
k_sec = 0.2; % maximum rate of insulin secretion
Ks = 100; % half-maximum concentration for insulin secretion
n = 2; % Hill coefficient
k_ex = 0.01; % glucose excretion rate
G_thr = 180; % renal glucose threshold
Gb = 90; % basal glucose level (mg/dL)

% Initial conditions
G = zeros(1, length(time));
I = zeros(1, length(time));
D = zeros(1, length(time));
L = zeros(1, length(time));

G(1) = 100;
I(1) = 10;
D(1) = 200; % initial digestive glucose (could be a meal)
L(1) = kp * (Gb - G(1));

% Solvingn using Runge-Kutta 4th Order Method
for i = 1:(length(time) - 1)
    % Update D
    k1 = -k_abs * D(i) / (Km + D(i));
    k2 = -k_abs * (D(i) + 0.5*dt*k1) / (Km + (D(i) + 0.5*dt*k1));
    k3 = -k_abs * (D(i) + 0.5*dt*k2) / (Km + (D(i) + 0.5*dt*k2));
    k4 = -k_abs * (D(i) + dt*k3) / (Km + (D(i) + dt*k3));
    D(i+1) = D(i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    % Update G
    g_ex = k_ex * max(0, G(i) - G_thr);
    k1 = k_abs * D(i) / (Km + D(i)) - ku * G(i) * I(i) + L(i) - g_ex;
    k2 = k_abs * (D(i) + 0.5*dt*k1) / (Km + (D(i) + 0.5*dt*k1)) - ku * (G(i)+0.5*dt*k1) * (I(i)+0.5*dt*k1) + L(i) - g_ex;
    k3 = k_abs * (D(i) + 0.5*dt*k2) / (Km + (D(i) + 0.5*dt*k2)) - ku * (G(i)+0.5*dt*k2) * (I(i)+0.5*dt*k2) + L(i) - g_ex;
    k4 = k_abs * (D(i) + dt*k3) / (Km + (D(i) + dt*k3)) - ku * (G(i)+dt*k3) * (I(i)+dt*k3) + L(i) - g_ex;
    G(i+1) = G(i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    % Update I
    sec_rate = k_sec * G(i)^n / (Ks^n + G(i)^n);
    k1 = sec_rate - kd * I(i);
    k2 = sec_rate - kd * (I(i) + 0.5*dt*k1);
    k3 = sec_rate - kd * (I(i) + 0.5*dt*k2);
    k4 = sec_rate - kd * (I(i) + dt*k3);
    I(i+1) = I(i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

% Solving using Euler's Method
G_euler = zeros(1, length(time));
I_euler = zeros(1, length(time));
D_euler = zeros(1, length(time));

G_euler(1) = 100;
I_euler(1) = 10;
D_euler(1) = 200;

for i = 1:(length(time) - 1)
    % Update D
    D_euler(i+1) = D_euler(i) - k_abs * D_euler(i) / (Km + D_euler(i)) * dt;

    % Update G
    g_ex = k_ex * max(0, G_euler(i) - G_thr);
    G_euler(i+1) = G_euler(i) + (k_abs * D_euler(i) / (Km + D_euler(i)) - ku * G_euler(i) * I_euler(i) + L(i) - g_ex) * dt;

    % Update I
    sec_rate = k_sec * G_euler(i)^n / (Ks^n + G_euler(i)^n);
    I_euler(i+1) = I_euler(i) + (sec_rate - kd * I_euler(i)) * dt;
end

% Solving using Runge-Kutta 2nd Order Method (RK2)
G_rk2 = zeros(1, length(time));
I_rk2 = zeros(1, length(time));
D_rk2 = zeros(1, length(time));

G_rk2(1) = 100;
I_rk2(1) = 10;
D_rk2(1) = 200;

for i = 1:(length(time) - 1)
    % Update D
    k1 = -k_abs * D_rk2(i) / (Km + D_rk2(i));
    k2 = -k_abs * (D_rk2(i) + dt) / (Km + (D_rk2(i) + dt));
    D_rk2(i+1) = D_rk2(i) + 0.5 * (k1 + k2) * dt;

    % Update G
    g_ex = k_ex * max(0, G_rk2(i) - G_thr);
    k1 = k_abs * D_rk2(i) / (Km + D_rk2(i)) - ku * G_rk2(i) * I_rk2(i) + L(i) - g_ex;
    k2 = k_abs * (D_rk2(i) + dt) / (Km + (D_rk2(i) + dt)) - ku * (G_rk2(i) + dt * k1) * (I_rk2(i) + dt * k1) + L(i) - g_ex;
    G_rk2(i+1) = G_rk2(i) + 0.5 * (k1 + k2) * dt;

    % Update I
    sec_rate = k_sec * G_rk2(i)^n / (Ks^n + G_rk2(i)^n);
    k1 = sec_rate - kd * I_rk2(i);
    k2 = sec_rate - kd * (I_rk2(i) + dt * k1);
    I_rk2(i+1) = I_rk2(i) + 0.5 * (k1 + k2) * dt;
end

% Plot results
figure;
subplot(3,1,1);
plot(time, G, 'r', time, G_euler, 'g', time, G_rk2, 'b');
title('Glucose Concentration Over Time');
xlabel('Time (minutes)');
ylabel('Glucose (mg/dL)');
legend('RK4', 'Euler', 'RK2');


subplot(3,1,2);
plot(time, I, 'r', time, I_euler, 'g', time, I_rk2, 'b');
title('Insulin Concentration Over Time');
xlabel('Time (minutes)');
ylabel('Insulin (mU/L)');
legend('RK4', 'Euler', 'RK2');


subplot(3,1,3);
plot(time, D, 'r', time, D_euler, 'g', time, D_rk2, 'b');
title('Digestive Tract Glucose');
xlabel('Time (minutes)');
ylabel('Glucose (mg/dL)');
legend('RK4', 'Euler', 'RK2');
