function Hd = IIR_2a
%IIR_2A Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.8 and Signal Processing Toolbox 8.4.
% Generated on: 25-Jun-2021 09:53:36

% Butterworth Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are normalized to 1.

N  = 2;       % Order
Fc = 0.0763;  % Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass('N,F3dB', N, Fc);
Hd = design(h, 'butter');

% [EOF]
