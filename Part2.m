%% Part 2: The Flipped Circuit

clear all;
close all;

%% implementation

t = 5; % total time of 5 seconds, as seen in figure 2 of the handout. 

R = 1; % 1000; % 1 k ohms
C = 1; % 1e-6; % 1 micro F
h = 8e-4; % sampling rate in seconds per sample

%% Construction and Execution of the model 

% |V| is the voltage in the system across all the different samples.  
% |V|_1 is V_R
% |V|_2 is V_C
% |V|_3 is V_in
r = 1;
c = 2;
in = 3;
freq = [50 1000];

V = zeros(length(freq), 3, t/h); % t/h is the total number of indicies 'k'

for f = 1:length(freq)
    for k = 1:t/h
        V(f, in, k) = 5 * sin(2 * pi * freq(f) * k);
        V(f, r, k) = V(f, in, k) - V(c, k); % Equation #8
        V(f, c, k+1) = (1 - (h / (R * C))) * V(f, c, k) + (h / (R * C)) * V(f, in, k); % Equation #10
    end
end


%% Plotting of data

figure;
hold on;
% plot(squeeze(V(1, in, :)));
 plot(squeeze(V(1, c, :)));
% plot(squeeze(V(1, r, :)));
hold off;
xlabel("Time (milliseconds)");
ylabel("Voltage (V)");
title("Approximated Charge Curve vs Time");
legend("V_c", "V_in", "location", "best");

% figure;
% hold on;
% fplot(@(t) 5 * sin(2 * pi * 50 * t), [0 0.04]);
% fplot(@(t) 5 * sin(2 * pi * 1000 * t), [0 0.04]);
% hold off;
% xlabel("Time (s)");
% ylabel("Voltage (V)");