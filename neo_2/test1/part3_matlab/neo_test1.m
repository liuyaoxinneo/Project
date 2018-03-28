%此工程用于画图
%----------
clear all
close all
%-----part1-----
%2~298，取中间的一个单元150，画出接收信号的波形
x=1:2E5;
b=load('E:\\Project\\neo_2\\test1\\data1\\150x_rev.txt');
a=load('E:\\Project\\neo_2\\test1\\data4\\150x_rev.txt');


%经过滤波之后的波形
c=load('E:\\Project\\neo_2\\test1\\part2\\part2\\150B.txt');
c_length=length(c)
% subplot(3,1,3);
% plot(c)
bs=fft(c,2E5);
B=bs.*conj(bs);
plot(B)%画出信号的功率谱