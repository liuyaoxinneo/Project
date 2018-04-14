%此工程用于画图
%----------
clear all
close all
%-----part1-----
%2~298，取中间的一个单元150，画出接收信号的波形
f=2e7;
x=0+1/f:1/f:0.01;
rev_sig=load('E:\\Project\\neo_3\\p1\\data6\\150x_rev.txt');
output=rev_sig(1:124800);
figure(1);
%plot(x,rev_sig);
plot(x(1:124800),output);
title('时域波形');
xlabel('时间 s');
ylabel('速度 m/s');
length(x)
length(rev_sig);
figure(2);
y=abs(fft(output));
semilogy(y);%使得y方向为log坐标，x为线性坐标
figure(3);
pwelch(output);
% title('左侧为50mm长的声源');


%经过滤波之后的波形
c=load('E:\\Project\\neo_2\\test1\\part2\\part2\\150B.txt');
c_length=length(c);
% subplot(3,1,3);
% plot(c)
%plot(abs(fft(rev_sig)))
% bs=fft(rev_sig,2E5);
% B=bs.*conj(bs);
% plot(B)%画出信号的功率谱