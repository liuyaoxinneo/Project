clear;

   
    for i=1:10
   
    b=sprintf('C:\\Documents and Settings\\Administrator\\My Documents\\Visual Studio 2005\\Projects\\nondefect\\nondefect\\%dx_rev.txt',i);
    file=load(b);
   deltt=0.5e-7;
   
   %deltt=time(2)-time(1);
   
   [fy,f,ff,ffy,fs,fsy,ft,fty,fo,foy,ffif,ffify]=myffth(file,deltt);
   [fy1,f1,ff1,ffy1,fs1,fsy1,ft1,fty1,ffif1,ffify1]=mypower(file,deltt);
    A1(i)=ffy1;
    A2(i)=fsy1;
    A3(i)=fty1;
    fname1=sprintf('sf_%d',i);
    fname2=sprintf('sfy_%d',i);
    fname3=sprintf('sf1_%d',i);
    fname4=sprintf('sfy1_%d',i);
    save(fname1,'f','-ascii','-double');
    save(fname2,'fy','-ascii','-double');
    save(fname3,'f1','-ascii','-double');
    save(fname4,'fy1','-ascii','-double');
    end