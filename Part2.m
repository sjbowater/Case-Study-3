%% Part 2: The Flipped Circuit

clear all;
close all;

%% implementation

t = 0.005; % total time of 0.04 seconds, as seen in figure 4 of the handout. 
R = 1e3;  % 1000; % 1 k ohms
C = 1e-6; % 1e-6; % 1 micro F
h = 8e-6; % sampling rate in seconds per sample

% We're currently only considering the frequency being 50 Hz.
timesteps = 0:h:t;
% timesteps, all within the time domain. 
% there will be length(timesteps) samples taken. 
Vinput50 = 5 * sin(2 * pi * 50 * timesteps);
Vresistor50 = zeros(1,length(timesteps));
Vcapacitor50 = zeros(1,length(timesteps));

%% Construction and Execution of the model at 50Hz
for i = 1:length(timesteps)
    Vcapacitor50(i+1) = (1 - (h / (R * C))) * Vcapacitor50(i) + (h / (R * C)) * Vinput50(i); % Equation #10
end

Vresistor50 = Vinput50 - Vcapacitor50(1:end-1); % Equation #08
%% Plotting of data
figure(1);
hold on;
plot(Vcapacitor50);
plot(Vresistor50);
plot(Vinput50);
hold off;
xlabel("Time (s)");
ylabel("Voltage (V)");
title("Approximated Charge Curve vs Time");
legend("V_c", "V_r", "V_i_n", "location", "best");
set(gca, 'XTick', 0:(length(timesteps)-1)/4:length(timesteps)-1)
set(gca, 'XTickLabel', 0:0.01:t)
xlim([0 t / h + 0.01]);

%% Construction and Execution of the model at 1000Hz
% In order for the graph to be legible, we must decrease the total
% timestep. 

t = 0.001;
timesteps = 0:h:t;
Vinput1000 = 5 * sin(2 * pi * 1000 * timesteps);
Vcapacitor1000 = zeros(1,length(timesteps));

for i = 1:length(timesteps)                         
    Vcapacitor1000(i+1) = (1 - (h / (R * C))) * Vcapacitor1000(i) + (h / (R * C)) * Vinput1000(i); % Equation #10
end

Vresistor1000 = Vinput1000 - Vcapacitor1000(1:end-1);                                              % Equation #08
%% Plotting of data
figure(2);
hold on;
plot(timesteps, Vcapacitor1000(1:end-1));
plot(timesteps, Vresistor1000);
plot(timesteps, Vinput1000);
hold off;
xlabel("Time (s)");
ylabel("Voltage (V)");
title("Approximated Charge Curve vs Time");
legend("V_c", "V_r", "V_i_n", "location", "best");

%% Task 2.3: The Transfer Functions
% notice that 10 is 10^1 and 10k is 10^4. |powers| is log(10) .. log(10k)
% with steps in logarithmic units. 
t = 0.04;
timesteps = 0:h:t;
powers = 1:0.001:4;
Vc = 0;

for p = 1:length(powers)
    % we're using 10^powers(i) to decode the frequency.
    Vin = 5 * sin(2 * pi * 10^powers(p) * [0:length(timesteps)] * h);
    for i = 1:length(timesteps)
        Vc(i+1) = (1 - (h / (R * C))) * Vc(i) + (h / (R * C)) * Vin(i);
    end
    
    Vr = Vin - Vc;
    
    Hc(p) = max(Vc) / max(Vin);
    Hr(p) = max(Vr) / max(Vin);
end

figure;
hold on;
plot(powers, Hc);
plot(powers, Hr);
hold off;
xlabel("Logarithm of the Frequency in log(Hz)");
ylabel("Logarithm of the Transfer Function in log(V)");
title("Frequency vs Voltage");
legend("H_c", "H_r");

