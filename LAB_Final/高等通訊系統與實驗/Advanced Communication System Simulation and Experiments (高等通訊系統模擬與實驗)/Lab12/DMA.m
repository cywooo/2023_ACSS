function Hd = DMA
%DMA Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.6 and Signal Processing Toolbox 8.2.
% Generated on: 29-May-2020 22:13:13

% Butterworth Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are normalized to 1.

Fpass = 0.25;        % Passband Frequency
Fstop = 0.4375;      % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 20;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop);
Hd = design(h, 'butter', 'MatchExactly', match);

% [EOF]
