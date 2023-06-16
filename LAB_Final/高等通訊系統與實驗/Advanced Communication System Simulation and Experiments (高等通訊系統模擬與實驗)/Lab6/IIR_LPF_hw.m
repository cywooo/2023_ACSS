function Hd = IIR_LPF_hw
%IIR_LPF_HW Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.6 and DSP System Toolbox 9.8.
% Generated on: 11-Apr-2020 16:29:15

% Butterworth Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are normalized to 1.

Fpass = 1/3;         % Passband Frequency
Fstop = 0.5;         % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 80;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop);
Hd = design(h, 'butter', 'MatchExactly', match);

% [EOF]
