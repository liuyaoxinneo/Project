%�˹������ڻ�ͼ
%----------
clear all
close all
%-----part1-----
%2~298��ȡ�м��һ����Ԫ150�����������źŵĲ���
x=1:2E5;
b=load('E:\\Project\\neo_2\\test1\\data1\\150x_rev.txt');
a=load('E:\\Project\\neo_2\\test1\\data4\\150x_rev.txt');


%�����˲�֮��Ĳ���
c=load('E:\\Project\\neo_2\\test1\\part2\\part2\\150B.txt');
c_length=length(c)
% subplot(3,1,3);
% plot(c)
bs=fft(c,2E5);
B=bs.*conj(bs);
plot(B)%�����źŵĹ�����