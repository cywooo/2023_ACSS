function Hd = Lab15_HW_p2_IIR
%LAB15_HW_P2_IIR Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.10 and Signal Processing Toolbox 8.6.
% Generated on: 29-May-2023 19:57:23

% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in MHz.
Fs = 64;  % Sampling Frequency

Fstop1 = 0.5;         % First Stopband Frequency
Fpass1 = 1;           % First Passband Frequency
Fpass2 = 3;           % Second Passband Frequency
Fstop2 = 4.6;         % Second Stopband Frequency
Astop1 = 60;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 80;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);

% [EOF]
