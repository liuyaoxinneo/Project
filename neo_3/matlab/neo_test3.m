clc

T=5e-8;%����ʱ����
fs=1/T;%����Ƶ��
L=2e5;%�źų���
t=(0:L-1)*T;

% f1=50,f2=120;%�ź��е�����Ƶ�ʷ���
% S=0.7*sin(2*pi*f1*t)+sin(2*pi*f2*t);
% X=S+2*randn(size(t));%�ź��м����������ĸ���

X=load('E:\\Project\\neo_3\\p1\\data6\\150x_rev.txt');
figure(1);
plot(t,X);

Y=fft(X(1:1.248e5));
P2=abs(Y/L);%˫��Ƶ��
P1=P2(1:L/2+1);
P1(2:end-1)=2*P1(2:end-1);%����Ƶ��

f=fs*(0:(L/2))/L;
figure(2);
plot(f,P1)
xlabel('f (Hz)')
ylabel('|P1(f)|')