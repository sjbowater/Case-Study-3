%% Case Study 3: Circuits
% Part 1: modeling an RC circuit.

clear all;
close all;
%% Equations of Note

t = 5e-3;      % total time of 5e-3 seconds, as seen in figure 2 of the handout. 
R = 1e3;       % 1 k ohms
C = 1e-6;      % 1 micro F
h = 8e-5;      % sampling rate in seconds per sample
hprime = 8e-4; % Sampling rate in seconds per sample as given in Task 1 Part 2.2.

timesteps  = 0:h:t;
Vinput     = 5 * ones(1, length(timesteps));
Vcapacitor = zeros(1, length(timesteps));

%% Construction and Execution of the model with h.

for k = 1:length(timesteps)
    Vcapacitor(k+1) = (1 - (h / (R * C))) * Vcapacitor(k) + (h / (R * C)) * Vinput(k); % Equation #10
end
%% Plotting of data with h

figure(1);
hold on;
plot(timesteps, Vcapacitor(1:end-1));
plot(timesteps, Vinput);
hold off;
xlabel("Time (s)");
ylabel("Voltage (V)");
title("Approximated Charge Curve vs Time");
legend("V_c", "V_i_n", "location", "best");

%% Construction and Execution of the model with h'.
timestepsprime = 0:hprime:t;
Vinputprime = 5 * ones(1, length(timestepsprime));
Vcapacitorprime = zeros(1, length(timestepsprime));
for k = 1:length(timestepsprime)
    Vcapacitorprime(k+1) = (1 - (hprime / (R * C))) * Vcapacitorprime(k) + (hprime / (R * C)) * Vinputprime(k); % Equation #10
end

%% Plotting of data with h'
figure(2);
hold on;
plot(timestepsprime, Vcapacitorprime(1:end-1));
plot(timestepsprime, Vinputprime);
hold off;
xlabel("Time (s)");
ylabel("Voltage (V)");
title("Approximated Charge Curve vs Time");
legend("V_c", "V_in", "location", "best");

%% Task 2
figure(3);
hold on;
% this is the theoretical charge function
fplot(@(k) 5*(1-exp(-k/(R*C))), [0 t]);
% this is |V_in| = 5.
fplot(5, [0 t]);
hold off;
xlabel("Time (s)");
ylabel("Voltage (V)");
title("Theoretical Charge Curve vs Time");
legend("V_c", "V_in", "location", "best");

% A comparison of the charges:
disp((5 * (1 - exp(-t / (R*C)))) - Vcapacitor(end));