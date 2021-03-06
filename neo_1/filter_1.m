clear all;
close all;
clc;

Fs = 20000;  % Sampling Frequency
F=250;
Fstop1 = F-2;                     
Fpass1 = F-1;        
Fpass2 = F+1;         
Fstop2 = F+2; 
Astop1 =80;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 80;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly
% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass(Fstop1/(Fs/2), Fpass1/(Fs/2), Fpass2/(Fs/2), Fstop2/(Fs/2), Astop1, Apass, ...
                      Astop2);
Hd = butter(h, 'MatchExactly', match);

% [EOF]
for i=2:298
    b = sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\Plane wave\\one\\Console1\\Console1\\%dx_rev.txt', i);
  
    [sig]=load(b); 
    sig=filter(Hd,sig);
    sig=fliplr(sig');
    sig=sig';
    sig=sig(1:6000);
    a=zeros(1000,1);
    sig=[sig;a];
    
    fname=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\Plane wave\\TR one PW\\Console1\\Console1\\%dB.txt', i);

    save(fname,'sig','-ascii','-double');
end