%% Part 2: The Flipped Circuit

clear all;
close all;

%% implementation

t = 0.04; % total time of 0.04 seconds, as seen in figure 4 of the handout. 
R = 1e3;  % 1000; % 1 k ohms
C = 1e-6; % 1e-6; % 1 micro F
h = 8e-5; % sampling rate in seconds per sample

% We're currently only considering the frequency being 50 Hz.
timesteps = 0:h:t;
% timesteps, all within the time domain. 
% there will be length(timesteps) samples taken. 
Vinput50 = 5 * sin(2 * pi * 50 * timesteps);
Vresistor50 = zeros(1,length(timesteps));
Vcapacitor50 = zeros(1,length(timesteps));

%% Construction and Execution of the model with 1000Hz
for k = 1:length(timesteps)
    Vresistor50(k) = Vinput50(k) - Vcapacitor50(k);                                          % Equation #08
    Vcapacitor50(k+1) = (1 - (h / (R * C))) * Vcapacitor50(k) + (h / (R * C)) * Vinput50(k); % Equation #10
end
%% Plotting of data
figure(1);
hold on;
plot(Vcapacitor50(:));
plot(Vresistor50(:));
plot(Vinput50(:));
hold off;
xlabel("Time (s)");
ylabel("Voltage (V)");
title("Approximated Charge Curve vs Time (Freq: 50 Hz)");
legend("V_c", "V_r", "V_i_n", "location", "best");
set(gca, 'XTick', 0:(length(timesteps)-1)/4:length(timesteps)-1)
set(gca, 'XTickLabel', 0:0.01:t)
xlim([0 t / h + 0.01]);

%% Construction and Execution of the model at 1000Hz
Vinput1000 = 5 * sin(2 * pi * 1000 * timesteps);
Vresistor1000 = zeros(1,length(timesteps));
Vcapacitor1000 = zeros(1,length(timesteps));

for k = 1:length(timesteps)
    Vresistor1000(k) = Vinput1000(k) - Vcapacitor1000(k);                                          % Equation #08
    Vcapacitor1000(k+1) = (1 - (h / (R * C))) * Vcapacitor1000(k) + (h / (R * C)) * Vinput1000(k); % Equation #10
end
%% Plotting of data
figure(2);
hold on;
plot(Vcapacitor1000(:));
plot(Vresistor1000(:));
plot(Vinput1000(:));
hold off;
xlabel("Time (s)");
ylabel("Voltage (V)");
title("Approximated Charge Curve vs Time (Freq: 1000Hz)");
legend("V_c", "V_r", "V_i_n", "location", "best");
set(gca, 'XTick', 0:(length(timesteps)-1)/4:length(timesteps)-1)
set(gca, 'XTickLabel', 0:0.01:t)
xlim([0 t / h + 0.01]);