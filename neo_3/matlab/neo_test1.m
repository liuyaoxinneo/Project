%�˹������ڻ�ͼ
%----------
clear all
close all
%-----part1-----
%2~298��ȡ�м��һ����Ԫ150�����������źŵĲ���
f=2e7;
x=0+1/f:1/f:0.01;
rev_sig=load('E:\\Project\\neo_3\\p1\\data6\\150x_rev.txt');
output=rev_sig(1:124800);
figure(1);
%plot(x,rev_sig);
plot(x(1:124800),output);
title('ʱ����');
xlabel('ʱ�� s');
ylabel('�ٶ� m/s');
length(x)
length(rev_sig);
figure(2);
y=abs(fft(output));
semilogy(y);%ʹ��y����Ϊlog���꣬xΪ��������
figure(3);
pwelch(output);
% title('���Ϊ50mm������Դ');


%�����˲�֮��Ĳ���
c=load('E:\\Project\\neo_2\\test1\\part2\\part2\\150B.txt');
c_length=length(c);
% subplot(3,1,3);
% plot(c)
%plot(abs(fft(rev_sig)))
% bs=fft(rev_sig,2E5);
% B=bs.*conj(bs);
% plot(B)%�����źŵĹ�����