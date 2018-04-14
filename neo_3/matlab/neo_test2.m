clc
fs=1000;%����Ƶ��
T=1/fs;%����ʱ����
L=1500;%�źų���
t=(0:L-1)*T;
f1=50,f2=120;%�ź��е�����Ƶ�ʷ���
S=0.7*sin(2*pi*f1*t)+sin(2*pi*f2*t);
X=S+2*randn(size(t));%�ź��м����������ĸ���
figure(1);
plot(t(1:50),X(1:50));

Y=fft(X);
P2=abs(Y/L);%˫��Ƶ��
P1=P2(1:L/2+1);
P1(2:end-1)=2*P1(2:end-1);%����Ƶ��

f=fs*(0:(L/2))/L;
figure(2);
plot(f,P1)
xlabel('f (Hz)')
ylabel('|P1(f)|')
