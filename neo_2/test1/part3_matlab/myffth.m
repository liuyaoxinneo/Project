function [fy,f,ff,ffy,fs,fsy,ft,fty,fo,foy,ffif,ffify]=myffth(amp,deltt)
[M,N]=size(amp);
N=max(M,N);
win=hanning(N);
amp=amp.*win;
fy=abs(fft(amp,N)/N);
fs=1/deltt;
f0=fs/N;
f=0:f0:f0*(N-1);
f=f';

ffn=0;
ffy=0;
for j=1:N
    if abs(fy(j))>ffy
        ffy=abs(fy(j));
        ffn=j;
    end
end
ff=(ffn-1)*f0;

fs=0;
fsy=0;
for j=2*ffn-20:2*ffn+20
     if abs(fy(j))>fsy
         fsy=abs(fy(j));
         fs=j;
     end
end
fs=(fs-1)*f0;

ft=0;
fty=0;
for j=3*ffn-20:3*ffn+20
     if abs(fy(j))>fty
         fty=abs(fy(j));
         ft=j;
     end
end
ft=(ft-1)*f0;

fo=0;
foy=0;
for j=4*ffn-20:4*ffn+20
     if abs(fy(j))>foy
         foy=abs(fy(j));
         fo=j;
     end
end
fo=(fo-1)*f0;

ffif=0;
ffify=0;
for j=5*(ffn-1):5*(ffn-1)
    if abs(fy(j))>ffify
        ffify=abs(fy(j));
        ffif=j;
    end
end
ffif=ffif*f0;
ffn

         
