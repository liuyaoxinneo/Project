function [fy,f,ff,ffy,fs,fsy,ft,fty,ffif,ffify]=mypower(amp,deltt);
[M,N]=size(amp);
N=max(M,N);
win=hanning(N);

fs=1/deltt;
[fy,f]=periodogram(amp,win,N,fs);
f0=fs/N;


ffn=0;
ffy=0;
for j=1:floor(N/2)
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

ffif=0;
ffify=0;
for j=5*ffn-20:5*ffn+20
    if abs(fy(j))>ffify
        ffify=abs(fy(j));
        ffif=j;
    end
end
ffif=(ffif-1)*f0;