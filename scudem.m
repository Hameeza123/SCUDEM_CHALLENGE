% Parameters
lambda = 0.1;    % Decay rate of persuasiveness (P)
alpha = 0.05;    % Influence of mood and messaging on behavioral intention (B)
kappa = 0.15;    % Decay rate of recall (R)
gamma = 0.2;     % Influence of persuasiveness on behavioral intention (B)
eta = 0.3;       % Effect of mood on persuasiveness (P)
tau = 20;        % Time constant for decay of behavioral intention (B)
epsilon = 0.1;   % Efficacy contribution factor (not used in this new formula)
delta = 0.05;    % Efficacy decay rate (not used in this new formula)

sigmoid = @(x) 1 ./ (1 + exp(-x));  % Define the sigmoid function

% Initial conditions
P0 = 0.5;        % Initial persuasiveness (P)
R0 = 0.7;        % Initial recall (R)
B0 = 0.6;        % Initial behavioral intention (B)
Q0 = sqrt(B0^2 + R0^2 + P0^2);  % Initial efficacy (Q) based on B, R, P

% Time span for the simulation
tspan = [0 10];  % Simulate from t = 0 to t = 5

% Range of values for m (messaging) and mu (mood)
m_values = linspace(-1, 1, 100);  % Messaging strength (m) from -1 to 1
mu_values = linspace(0, 1, 100);  % Mood (mu) from 0 to 1

% Allocate matrices to store final values of B, P, R, or Q
final_B = zeros(length(m_values), length(mu_values));  % Final behavioral intention
final_Q = zeros(length(m_values), length(mu_values));  % Final efficacy (Q)
final_P = zeros(length(m_values), length(mu_values));  % Final persuasiveness (P)
final_R = zeros(length(m_values), length(mu_values));  % Final recall (R)

for i = 1:length(m_values)
    for j = 1:length(mu_values)
        
        % Current values of m (messaging) and mu (mood)
        m = m_values(i);
        mu = mu_values(j);  
        if m  < 0
            m = 2*sigmoid(1.2*m)-1
        else
            m = 2*sigmoid(m)-1
        end
        

        %m = sigmoid(m)-1/2
        ode_system = @(t, y) [
            -lambda * y(1) + eta * mu * abs(m);           % dP/dt (persuasiveness)
            -kappa * y(2) + alpha * mu * m^2;        % dR/dt (recall)
            -alpha * exp(-t / tau) * y(3) + gamma * mu * m^2 * y(1);  % dB/dt (behavioral intention)
            0;  % dQ/dt (no longer needed as Q is directly calculated)
        ];
        initial_conditions = [P0; R0; B0; Q0];
        [t, solution] = ode45(ode_system, tspan, initial_conditions);

        % Extract the solutions for P, R, B, and Q
        P = solution(:, 1);
        R = solution(:, 2);
        B = solution(:, 3);

        % Calculate efficacy as sqrt(B^2 + R^2 + P^2)
        Q = sqrt(B.^2 + R.^2 + P.^2);

        % Store the final values of P, R, B, and Q
        final_B(i, j) = B(end);   % Store final behavioral intention (B)
        final_Q(i, j) = Q(end);   % Store final efficacy (Q)
        final_P(i, j) = P(end);   % Store final persuasiveness (P)
        final_R(i, j) = R(end);   % Store final recall (R)
    end
end

% Plot the bifurcation diagrams for final behavioral intention (B), efficacy (Q),
% persuasiveness (P), and recall (R)
figure;

% Bifurcation Diagram: Behavioral Intention (B)
subplot(2, 2, 1);
imagesc(mu_values, m_values, final_B);
colorbar;
xlabel('Mood (\mu)');
ylabel('Messaging Strength (m)');
title('Bifurcation Diagram: Behavioral Intention (B)');
axis xy;

% Bifurcation Diagram: Efficacy (Q)
subplot(2, 2, 2);
imagesc(mu_values, m_values, final_Q);
colorbar;
xlabel('Mood (\mu)');
ylabel('Messaging Strength (m)');
title('Bifurcation Diagram: Efficacy (Q)');
axis xy;

% Bifurcation Diagram: Persuasiveness (P)
subplot(2, 2, 3);
imagesc(mu_values, m_values, final_P);
colorbar;
xlabel('Mood (\mu)');
ylabel('Messaging Strength (m)');
title('Bifurcation Diagram: Persuasiveness (P)');
axis xy;

% Bifurcation Diagram: Recall (R)
subplot(2, 2, 4);
imagesc(mu_values, m_values, final_R);
colorbar;
xlabel('Mood (\mu)');
ylabel('Messaging Strength (m)');
title('Bifurcation Diagram: Recall (R)');
axis xy;

% Display the final efficacy value at the end of the simulation
disp(['Final efficacy at t = ', num2str(t(end)), ' is ', num2str(final_Q(end))]);

% Example: Plot the time series of Behavioral Intention (B) for specific m and mu
% Select a specific pair of (m, mu) values for the time series plot
m_selected = 1;   % Example value for messaging strength (m)
mu_selected = 0;  % Example value for mood (mu)
%m = m*sigmoid(1.5*m)

% Define the system of ODEs with the selected m and mu
ode_system_selected = @(t, y) [
    -lambda * y(1) + eta * mu_selected * m_selected^2;           % dP/dt (persuasiveness)
    -kappa * y(2) + alpha * mu_selected * m_selected^2;        % dR/dt (recall)
    -alpha * exp(-t / tau) * y(3) + gamma * mu_selected * m_selected^2 * y(1);  % dB/dt (behavioral intention)
    0;  % dQ/dt (no longer needed as Q is directly calculated)
];

% Initial state vector [P0; R0; B0; Q0]
initial_conditions_selected = [P0; R0; B0; Q0];

% Solve the system of ODEs using ode45 for the selected parameters
[t_selected, solution_selected] = ode45(ode_system_selected, tspan, initial_conditions_selected);

% Extract the solutions for P, R, B, and Q
P_selected = solution_selected(:, 1);
R_selected = solution_selected(:, 2);
B_selected = solution_selected(:, 3);

% Calculate efficacy as sqrt(B^2 + R^2 + P^2)
Q_selected = sqrt(B_selected.^2 + R_selected.^2 + P_selected.^2);

% Plot the time series for Behavioral Intention (B), Persuasiveness (P),
% Recall (R), and Efficacy (Q)
figure;
subplot(2,2,1);
plot(t_selected, B_selected, 'b', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('Behavioral Intention (B)');
title('Time Series of Behavioral Intention (B) for Selected m and mu');
grid on;

subplot(2,2,2);
plot(t_selected, P_selected, 'g', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('Persuasiveness (P)');
title('Time Series of Persuasiveness (P) for Selected m and mu');
grid on;

subplot(2,2,3);
plot(t_selected, R_selected, 'm', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('Recall (R)');
title('Time Series of Recall (R) for Selected m and mu');
grid on;

subplot(2,2,4);
plot(t_selected, Q_selected, 'r', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('Efficacy (Q)');
title('Time Series of Efficacy (Q) for Selected m and mu');
grid on;
