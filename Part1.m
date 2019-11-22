%% Case Study 3: Circuits
% Part 1: modeling an RC circuit.

clear all;
close all;
%% Equations of Note

t = 5e-3; % total time of 5 seconds, as seen in figure 2 of the handout. 
R = 1e3; % 1000; % 1 k ohms
C = 1e-6; % 1e-6; % 1 micro F
h = 8e-5; % sampling rate in seconds per sample
hprime = 8e-4;

timesteps = 0:h:t;
Vinput = 5 * ones(1, length(timesteps));
Vresistor = zeros(1, length(timesteps));
Vcapacitor = zeros(1, length(timesteps));

%% Construction and Execution of the model with h.
for k = 1:length(timesteps)
    Vresistor(k) = Vinput(k) - Vcapacitor(k);                                          % Equation #08
    Vcapacitor(k+1) = (1 - (h / (R * C))) * Vcapacitor(k) + (h / (R * C)) * Vinput(k); % Equation #10
end
%% Plotting of data
figure(1);
hold on;
plot(Vcapacitor(:));
plot(Vinput(:));
hold off;
xlabel("Time (s)");
ylabel("Voltage (V)");
title("Approximated Charge Curve vs Time (Frequency: " + num2str(h) + " s/sample)");
legend("V_c", "V_in", "location", "best");
set(gca, 'XTick', 0:(length(timesteps)-1)/5:length(timesteps)-1)
set(gca, 'XTickLabel', 0:1e-3:t)
xlim([0 t / h + 0.01]);

%% Construction and Execution of the model with h'.
for k = 1:length(timesteps)
    Vresistor(k) = Vinput(k) - Vcapacitor(k);                                                    % Equation #08
    Vcapacitor(k+1) = (1 - (hprime / (R * C))) * Vcapacitor(k) + (hprime / (R * C)) * Vinput(k); % Equation #10
end

figure(2);
hold on;
plot(timesteps, Vcapacitor(1:end-1));
plot(timesteps, Vinput(:));
hold off;
xlabel("Time (s)");
ylabel("Voltage (V)");
title("Approximated Charge Curve vs Time (Frequency: " + num2str(hprime) + " s/sample)");
legend("V_c", "V_in", "location", "best");
set(gca, 'XTick', 0:(length(timesteps)-1)/5:length(timesteps)-1)
set(gca, 'XTickLabel', 0:1e-3:t)
xlim([0 timesteps(5)])

% %% Task 2
% % theoretical_output = 5 * (1 - exp(-t/(R*C)));
% figure;
% hold on;
% % this is the theoretical charge function
% fplot(@(k) 5*(1-exp(-k/(R*C))), [0 t]);
% % this is |V_in| = 5.
% fplot(5, [0 t]);
% hold off;
% xlabel("Time (s)");
% ylabel("Voltage (V)");
% title("Theoretical Charge Curve vs Time");
% legend("V_c", "V_in", "location", "best");
% xlim([0 t]);
% ylim([0 5]);
% 
% % A comparison of the charges:
% disp((5 * (1 - exp(-t / (R*C)))) - Vcapacitor(end));