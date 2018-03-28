function [fy,f,ft]=myfourier(amp,deltt)
[M,NN]=size(amp);
N=max(M,NN);
if(N==M)
    win=hanning(N);
else
    win=hanning(N)';
end
amp=amp.*win;
fy=fft(amp,N)/N;              %对于周期信号才是除以N，非周期乘以deltt
fs=1/deltt;
f0=fs/N;
f=0:f0:f0*(N-1);
%NN=round(750000/f0);
%ft=fy(NN+1);
ft=max(fy);
