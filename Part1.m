%% Case Study 3: Circuits
% Part 1: modeling an RC circuit.

clear all;
close all;
%% Equations of Note

t = 5; % total time of 5 seconds, as seen in figure 2 of the handout. 

%  (8) v_(R,k) = V_(in,k) - V_(C,k)
% (10) V_(C_k+1) = (1-h/RC)v_(C,k) + h/RC * v_(in,k) 

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

V = zeros(3, t/h); % t/h is the total number of indicies 'k'
V(3,:) = 5; % V_in is a constant 5V the entire time. 

for k = 1:t/h
    V(r, k) = V(in, k) - V(c, k); % Equation #8
    V(c, k+1) = (1 - (h / (R * C))) * V(c, k) + (h / (R * C)) * V(in, k); % Equation #10
end

%% Plotting of data

figure;
hold on;
plot(V(c, :));
plot(V(in, :));
hold off;
set(gca, 'XTick', 0:6250/5:k)
set(gca, 'XTickLabel', 0:5)
xlabel("Time (s)");
ylabel("Voltage (V)");
title("Approximated Charge Curve vs Time");
legend("V_c", "V_in", "location", "best");
xlim([0 k]);
ylim([0 5]);

%% Task 2
% theoretical_output = 5 * (1 - exp(-t/(R*C)));
figure;
hold on;
% this is the theoretical charge function
fplot(@(k) 5*(1-exp(-k/(R*C))), [0 5]);
% this is |V_in| = 5.
fplot(5, [0 5]);
hold off;
xlabel("Time (s)");
ylabel("Voltage (V)");
title("Theoretical Charge Curve vs Time");
legend("V_c", "V_in", "location", "best");
xlim([0 5]);
ylim([0 5]);

% A comparison of the charges:
disp((5 * (1 - exp(-t / (R*C)))) - V(c, end));