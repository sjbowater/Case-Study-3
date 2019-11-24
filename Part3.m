%% Part 3: Cascaded circuits

clear all;
close all;

%% 
t = 0.01;    % Total time simulation runs for. 
h = 8e-6;    % Sample rate in seconds per second
C1 = 0.68e-6;
C3=0.68e-6;% C_1 = C_2 = C_3 = 0.68 microF. 
R4=330;
R2=330;    % R_1 = R_2 = R_4 = 330 Ohms
f1 = 440;    % 440 Hz
f2 = 3000;   % 3k Hz

timesteps = 0:h:t;

% A = [ 1 -1 -1  0  0  0; % i_1 - i_2 - i_3 = 0     Equation 15
%       0  0  0  1  0  0; % V_in = V_in,k           Equation 19
%       0  0  0  1 -1  0; % V_in - V_1 = V_C1,k     Equation 20
%       0  0  0  0  1 -1; % V_1 - V_out = V_c3,k    Equation 21
%       0  R  0  0 -1  0; % R_2 * i_2 - V_1 = 0     Equation 17
%       0  0  R  0  0 -1; % R_4 * i_3 - V_out = 0   Equation 18
%     ];
A = [ 1 -1 -1  0  0  0;
     -1  1  1  0  0  0;
      0  0  0  1 -1  0;
      0  0  0  0 -1  1;
      0  0 -R4 0  0  1;
      0 -R2 0  0  1  0;
      0  0  0  1  0  0;
              ];

% % same as the A above but with equations in order. 
% A = [ 1 -1 -1  0  0  0;
%       0 -R  0  0  1  0;
%       0  0 -R  0  0  1;
%       0  0  0  1  0  0;
%       0  0  0  1 -1  0;
%       0  0  0  0  1 -1;
%      ];

Vc1  = zeros(1, length(timesteps));
Vc3  = zeros(1, length(timesteps));
Vout = zeros(1, length(timesteps));
Vin  = 5 * sin(2 * pi * f1 * timesteps) + sin(2 * pi * f2 * timesteps); % Equation 26.

i1 = zeros(1, length(timesteps));
i3 = zeros(1, length(timesteps));

for k=1:length(timesteps)
    
    x = linsolve(A, [0, 0, Vc1(k), -Vc3(k), 0, 0, Vin(k),]');
%     x = linsolve(A, [0 0 0 Vin(k) Vc1(k) Vout(k)]');
    i1(k)   = x(1);
    i2(k)   = x(2);
    i3(k)   = x(3);
    Vin(k)  = x(4);
    V1(k)  = x(5);
    Vout(k) = x(6);
    
    Vc1(k+1) = Vc1(k) + (h / C1) * i1(k); % Equation 24
    Vc3(k+1) = Vc3(k) + (h / C3) * i3(k); % Equation 25
end

figure;
hold on;
plot(timesteps, Vin);
plot(timesteps, Vout);
hold off;
legend("Vin", "Vout", "location", "best");

xlabel("Time (s)");
ylabel("Voltage (V)");
title("Voltage (V) vs Time (s) Circuit C");

%% Transfer Functions
