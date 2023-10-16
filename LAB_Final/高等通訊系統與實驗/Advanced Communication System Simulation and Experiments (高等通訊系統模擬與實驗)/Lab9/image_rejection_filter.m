function Hd = image_rejection_filter
%IMAGE_REJECTION_FILTER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.6 and Signal Processing Toolbox 8.2.
% Generated on: 10-May-2020 23:26:38

% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are normalized to 1.

Fstop1 = 0.28125;     % First Stopband Frequency
Fpass1 = 0.375;       % First Passband Frequency
Fpass2 = 0.625;       % Second Passband Frequency
Fstop2 = 0.71875;     % Second Stopband Frequency
Astop1 = 50;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 60;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2);
Hd = design(h, 'butter', 'MatchExactly', match);

% [EOF]
