function Hd = Lab9_demo_IIR_2
%LAB9_DEMO_IIR_2 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.5 and Signal Processing Toolbox 8.1.
% Generated on: 07-Apr-2023 10:08:27

% Butterworth Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are normalized to 1.

N  = 2;      % Order
Fc = 0.327;  % Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass('N,F3dB', N, Fc);
Hd = design(h, 'butter');

% [EOF]
