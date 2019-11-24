%% Part 3: Cascaded circuits

clear all;
close all;

%% Circuit C
t = 0.01;    % Total time simulation runs for, as seen in Figure 7. It's worth noting that if t = 0.001, the simulation appears better, but still wholly inacurate. 
h = 8e-6;    % Sample rate in seconds per second; this is the same rate from Part 1.
C = 0.68e-6; % C_1 = C_2 = C_3 = 0.68 microF. 
R = 330;     % R_1 = R_2 = R_4 = 330 Ohms
f1 = 440;    % 440 Hz
f2 = 3000;   % 3k Hz

timesteps = 0:h:t;

AC = [ 1 -1 -1  0  0  0; % i_1 - i_2 - i_3 = 0     Equation 15
       0  0  0  1  0  0; % V_in = V_in,k           Equation 19
       0  0  0  1 -1  0; % V_in - V_1 = V_c1,k     Equation 20
       0  0  0  0  1 -1; % V_1 - V_out = V_c3,k    Equation 21
       0  R  0  0 -1  0; % R_2 * i_2 - V_1 = 0     Equation 17
       0  0  R  0  0 -1; % R_4 * i_3 - V_out = 0   Equation 18
     ];

% The voltages we're tracking. 
Vc1  = zeros(1, length(timesteps));
Vc3  = zeros(1, length(timesteps));
Vout = zeros(1, length(timesteps));
Vin  = 5 * sin(2 * pi * f1 * timesteps) + sin(2 * pi * f2 * timesteps); % Equation 26.

% These could either be scalar or vector, as the values are tossed between
% each iteration of the loop. 
i1 = zeros(1, length(timesteps));
i3 = zeros(1, length(timesteps));

for k = 1:length(timesteps)
    % Solve the circuit. 
    x = linsolve(AC, [0, Vin(k), Vc1(k), Vc3(k), 0, 0]');
    % Manually storing values from x in Ax = b. These stored values are
    % used in the update equations. 
    i1(k)   = x(1);
%   i2      = x(2); i2 is never used outside of A.
    i3(k)   = x(3);
%   Vin(k)  = x(4); Vin is already calculated for all timesteps. 
    Vc1(k)  = x(4) - x(5);
    Vout(k) = x(6);
    
    % Update Equations
    Vc1(k+1) = Vc1(k) + (h / C) * i1(k); % Equation 24
    Vc3(k+1) = Vc3(k) + (h / C) * i3(k); % Equation 25
end

figure;
hold on;
plot(timesteps, Vin);
plot(timesteps, Vout);
hold off;
legend("Vin", "Vout", "location", "best");

xlabel("Time (s)");
ylabel("Voltage (V)");
title("Circuit C: Voltage (V) vs Time (s)");

%% Circuit D
AD = [ 1 -1 -1  0  0  0;
       R  0  0 -1  0  0;
       0  0  R  0  0 -1;
       0  0  0  1  0  0;
       0  0  0  1 -1  0;
       0  0  0  1 -1 -1;
     ];

% The voltages we're tracking. 
Vc2  = zeros(1, length(timesteps));
Vc3  = zeros(1, length(timesteps));
Vout = zeros(1, length(timesteps));
Vin  = 5 * sin(2 * pi * f1 * timesteps) + sin(2 * pi * f2 * timesteps); % Equation 26.

% These could either be scalar or vector, as the values are tossed between
% each iteration of the loop. 
i1 = zeros(1, length(timesteps));
i2 = zeros(1, length(timesteps));
i3 = zeros(1, length(timesteps));

for k = 1:length(timesteps)
    % Solve the circuit. 
    x = linsolve(AD, [0, 0, 0, Vin(k), Vc2(k), Vc3(k)]');
    % Manually storing values from x in Ax = b. These stored values are
    % used in the update equations. 
    i1(k)   = x(1);
    i2(k)   = x(2); 
    i3(k)   = x(3);
    Vin(k)  = x(4);
    V1(k)   = x(5);
    Vout(k) = x(6);
    
    Vc2(k) = Vin(k) - V1(k);
    
    % Update Equations
    Vc2(k+1) = Vc2(k) + (h / C) * i2(k); % Equation 24
    Vc3(k+1) = Vc3(k) + (h / C) * i3(k); % Equation 25
end

figure;
hold on;
plot(timesteps, Vin);
plot(timesteps, Vout);
hold off;
legend("Vin", "Vout", "location", "best");

xlabel("Time (s)");
ylabel("Voltage (V)");
title("Circuit D: Voltage (V) vs Time (s)");

% %% Transfer Functions
% % notice that 10 is 10^1 and 10k is 10^4. |powers| is log(10) .. log(10k)
% % with steps in logarithmic units. 
% powers = 1:0.001:4;
% 
% Vc1C  = zeros(1, length(timesteps));
% Vc3C  = zeros(1, length(timesteps));
% VoutC = zeros(1, length(timesteps));
% 
% Vc2D  = zeros(1, length(timesteps));
% Vc3D  = zeros(1, length(timesteps));
% VoutD = zeros(1, length(timesteps));
% 
% for p = 1:length(powers)
%     % we're using 10^powers(i) to decode the frequency.
%     Vin = 5 * sin(2 * pi * 10^powers(p) * [0:length(timesteps)] * h);
%     for i = 1:length(timesteps)
%         xC = linsolve(AC, [0, Vin(k), Vc1C(k), Vc3C(k),       0,      0]');
%         xD = linsolve(AD, [0,      0,       0,  Vin(k), Vc2D(k), Vc3D(k)]');
%         
%         i1C(k)   = xC(1);
%         i3C(k)   = xC(3);
%         Vc1C(k)  = xC(4) - xC(5);
%         VoutC(k) = xC(6);
%         
%         i1D(k)   = x(1);
%         i2D(k)   = x(2); 
%         i3D(k)   = x(3);
%         Vc2D(k)  = x(4) - x(5);
%         VoutD(k) = x(6);
%     
%     %   Update Equations
%         Vc1C(k+1) = Vc1C(k) + (h / C) * i1C(k); 
%         Vc3C(k+1) = Vc3C(k) + (h / C) * i3C(k); 
%         
%         Vc2D(k+1) = Vc2D(k) + (h / C) * i2D(k); 
%         Vc3D(k+1) = Vc3D(k) + (h / C) * i3D(k); 
%     end
%     
%     HC(p) = max(VoutC(k));
%     HD(p) = max(VoutD(k));
% end
% 
% figure;
% hold on;
% plot(powers, log(HC));
% plot(powers, log(HD));
% hold off;
% xlabel("Logarithm of the Frequency in log(Hz)");
% ylabel("Logarithm of the Transfer Function in log(V)");
% title("Frequency vs Voltage");