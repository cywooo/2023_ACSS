function Hd = Lab15_demo_IIR
%LAB15_DEMO_IIR Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.10 and Signal Processing Toolbox 8.6.
% Generated on: 29-May-2023 19:23:04

% Butterworth Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are normalized to 1.

Fpass = 0.2;         % Passband Frequency
Fstop = 0.25;        % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 80;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop);
Hd = design(h, 'butter', 'MatchExactly', match);

% [EOF]
