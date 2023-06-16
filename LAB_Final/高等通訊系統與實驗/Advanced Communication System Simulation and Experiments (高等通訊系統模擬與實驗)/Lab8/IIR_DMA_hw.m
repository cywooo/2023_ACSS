function Hd = IIR_DMA_hw
%IIR_DMA_HW Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.6 and DSP System Toolbox 9.8.
% Generated on: 24-Apr-2020 11:04:06

% Butterworth Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are normalized to 1.

Fpass = 1/64;    % Passband Frequency
Fstop = 1/32+3/4/32;   % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 16;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop);
Hd = design(h, 'butter', 'MatchExactly', match);

% [EOF]
