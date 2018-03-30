%此工程用于画图
%----------
clear all
close all
%-----part1-----
%2~298，取中间的一个单元150，画出接收信号的波形
f=5e6;
x=0+1/f:1/f:0.04;
rev_sig=load('E:\\Project\\neo_2\\test1\\data8\\150x_rev.txt');
figure(1);
plot(x(1:1e5),rev_sig(1:1e5));
length(x)
length(rev_sig);
figure(2);
pwelch(rev_sig);
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