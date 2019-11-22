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
Vinput = 5 * sin(2 * pi * 50 * timesteps);
Vresistor = zeros(1,length(timesteps));
Vcapacitor = zeros(1,length(timesteps));

%% Construction and Execution of the model 
for k = 1:t/h
    Vresistor(k) = Vinput(k) - Vcapacitor(k);
    Vcapacitor(k+1) = (1 - (h / (R * C))) * Vcapacitor(k) + (h / (R * C)) * Vinput(k);
end
%% Plotting of data

figure(1);
hold on;
plot(Vcapacitor(:));
plot(Vresistor(:));
plot(Vinput(:));
hold off;
% set(gca, 'XTick', 0:k/t:k)
% set(gca, 'XTickLabel', 0:t)
% xlabel("Time (s)");
xlabel("Sample Number")
ylabel("Voltage (V)");
title("Approximated Charge Curve vs Time");
legend("V_c", "V_r", "V_i_n", "location", "best");
xlim([0 t / h]);