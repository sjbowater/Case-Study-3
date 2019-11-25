% function cascadedRCfilter(Vin,h) receives a time-series voltage sequence
% sampled with interval h, and returns the output voltage sequence produced
% by a circuit
%
% inputs:
% Vin - time-series vector representing the voltage input to a circuit
% h - scalar representing the sampling interval of the time series in
% seconds
%
% outputs:
% Vout - time-series vector representing the output voltage of a circuit

function Vout = RCfilter(Vin,h)
%     f_high = 1 / (2 * pi * R1 * C1); % only frequencies higher than this
%     will remain. 
%     f_low  = 1 / (2 * pi * R * C2); % only frequencies lower than this
%     will remain.
    t=length(Vin)*h;
    timesteps = 0:h:t;
    R1 = 470;
    R4 = 16;
    C1 = 33e-6;
    C2 = 47e-6;
        
    Vc2 = 0;
    Vc3 = 0;

    A = [   1  -1  -1   0   0   0;
          -R1   0   0   1  -1   0;
            0   0  R4   0   0  -1;
            0   0   0   1   0   0;
            0   0   0   0   1   0;
            0   0 -R4   0   1   0;
        ];

    for k = 1:length(timesteps)-1
        x = linsolve(A, [0, 0, 0, Vin(k), Vc2(k), Vc3(k)]');
        Vc2(k+1) = Vc2(k) + (h / C1) * x(2);
        Vc3(k+1) = Vc3(k) + (h / C2) * x(3);
    end

    Vout = (Vc2 - Vc3).';
    Vout = Vout(1:end-1);
end