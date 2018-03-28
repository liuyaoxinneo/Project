%本程序用于先滤波再反转


Fs = 20000;  % Sampling Frequency

Fstop1 =240;         % First Stopband Frequency               !300*500在2750khz附近成十字交叉
Fpass1 =241;         % First Passband Frequency
Fpass2 =243;         % Second Passband Frequency
Fstop2 =244;         % Second Stopband Frequency
Astop1 = 100;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 100;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass(Fstop1/(Fs/2), Fpass1/(Fs/2), Fpass2/(Fs/2), Fstop2/(Fs/2), Astop1, Apass, ...
                      Astop2);
Hd = butter(h, 'MatchExactly', match);


% [EOF]
for i=2:298
    b=sprintf('C:\\Documents and Settings\\Administrator\\My Documents\\Visual Studio 2005\\Projects\\Console1\\Console1\\%dx_rev.txt',i);
    [sig]=load(b);
    
    %a=zeros(2000,1);
    %sig=[a;sig];
    %B_sig=filter(Hd,sig);
    %B_sig=sig;
   
    sig=filter(Hd,sig); 
    sig=sig(1:6000);
    sig=fliplr(sig');
    sig=sig';
    a=zeros(1000,1);
    sig=[sig;a];
    fname=sprintf('C:\\Documents and Settings\\Administrator\\My Documents\\Visual Studio 2005\\Projects\\Console2\\Console2\\%dB.txt',i);
    save(fname,'sig','-ascii','-double');
end