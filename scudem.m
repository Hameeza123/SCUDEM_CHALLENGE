% Parameters (fixed)
lambda = 0.1;    % Decay rate of persuasiveness (P)
alpha = 0.05;    % Influence of mood and messaging on behavioral intention (B)
kappa = 0.15;    % Decay rate of recall (R)
gamma = 0.2;     % Influence of persuasiveness on behavioral intention (B)
eta = 0.3;       % Effect of mood on persuasiveness (P)
tau = 10;        % Time constant for decay of behavioral intention (B)
epsilon = 0.1;   % Efficacy contribution factor
delta = 0.05;    % Efficacy decay rate

% Initial conditions
P0 = 0.4;        % Initial persuasiveness (P)
R0 = 0.5;        % Initial recall (R)
B0 = 0.9;        % Initial behavioral intention (B)
Q0 = 0.3;        % Initial efficacy (Q)

% Time span for the simulation
tspan = [0 100];  % Simulate from t = 0 to t = 50

% Range of values for m (messaging) and mu (mood)
m_values = linspace(-1, 1, 100);  % Messaging strength (m) from -1 to 1
mu_values = linspace(0, 1, 100);  % Mood (mu) from 0 to 1

% Allocate matrix to store final values of B or Q
final_B = zeros(length(m_values), length(mu_values));  % Final behavioral intention
final_Q = zeros(length(m_values), length(mu_values));  % Final efficacy (Q)

% Loop over m and mu to simulate the system and store final results
for i = 1:length(m_values)
    for j = 1:length(mu_values)
        
        % Current values of m (messaging) and mu (mood)
        m = m_values(i);
        mu = mu_values(j);
        
        % Define the system of ODEs with current m and mu
        ode_system = @(t, y) [
            -lambda * y(1) + eta * mu * m;           % dP/dt (persuasiveness)
            -kappa * y(2) + alpha * mu * m^2;        % dR/dt (recall)
            -alpha * exp(-t / tau) * y(3) + gamma * mu * m^2 * y(1);  % dB/dt (behavioral intention)
            epsilon * (y(1)^2 + y(2)^2 + y(3)^2) - delta * y(4);  % dQ/dt (efficacy)
        ];

        % Initial state vector [P0; R0; B0; Q0]
        initial_conditions = [P0; R0; B0; Q0];

        % Solve the system of ODEs using ode45
        [t, solution] = ode45(ode_system, tspan, initial_conditions);

        % Extract the solutions for P, R, B, and Q
        P = solution(:, 1);
        R = solution(:, 2);
        B = solution(:, 3);
        Q = solution(:, 4);

        % Store the final value of B or Q
        final_B(i, j) = B(end);   % Store final behavioral intention (B)
        final_Q(i, j) = Q(end);   % Store final efficacy (Q)
    end
end

% Plot the bifurcation diagrams for final behavioral intention (B) and efficacy (Q)
figure;

% Bifurcation Diagram: Behavioral Intention (B)
subplot(1, 2, 1);
imagesc(mu_values, m_values, final_B');
colorbar;
xlabel('Mood (\mu)');
ylabel('Messaging Strength (m)');
title('Bifurcation Diagram: Behavioral Intention (B)');
axis xy;

% Bifurcation Diagram: Efficacy (Q)
subplot(1, 2, 2);
imagesc(mu_values, m_values, final_Q');
colorbar;
xlabel('Mood (\mu)');
ylabel('Messaging Strength (m)');
title('Bifurcation Diagram: Efficacy (Q)');
axis xy;

% Display the final efficacy value at the end of the simulation
disp(['Final efficacy at t = ', num2str(t(end)), ' is ', num2str(final_Q(end))]);

% Example: Plot the time series of Behavioral Intention (B) for specific m and mu
% Select a specific pair of (m, mu) values for the time series plot
m_selected = 1;   % Example value for messaging strength (m)
mu_selected = 0;  % Example value for mood (mu)

% Define the system of ODEs with the selected m and mu
ode_system_selected = @(t, y) [
    -lambda * y(1) + eta * mu_selected * m_selected;           % dP/dt (persuasiveness)
    -kappa * y(2) + alpha * mu_selected * m_selected^2;        % dR/dt (recall)
    -alpha * exp(-t / tau) * y(3) + gamma * mu_selected * m_selected^2 * y(1);  % dB/dt (behavioral intention)
    epsilon * (y(1)^2 + y(2)^2 + y(3)^2) - delta * y(4);  % dQ/dt (efficacy)
];

% Initial state vector [P0; R0; B0; Q0]
initial_conditions_selected = [P0; R0; B0; Q0];

% Solve the system of ODEs using ode45 for the selected parameters
[t_selected, solution_selected] = ode45(ode_system_selected, tspan, initial_conditions_selected);

% Extract the solutions for P, R, B, and Q
B_selected = solution_selected(:, 3);  % Behavioral Intention (B)

% Plot the time series for Behavioral Intention (B)
figure;
plot(t_selected, B_selected, 'b', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('Behavioral Intention (B)');
title('Time Series of Behavioral Intention (B) for Selected m and mu');
grid on;
