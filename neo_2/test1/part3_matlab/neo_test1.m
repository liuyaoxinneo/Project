%�˹������ڻ�ͼ
%----------
clear all
close all
%-----part1-----
%2~298��ȡ�м��һ����Ԫ150�����������źŵĲ���
f=5e6;
x=0+1/f:1/f:0.04;
rev_sig=load('E:\\Project\\neo_2\\test1\\data8\\150x_rev.txt');
figure(1);
plot(x(1:1e5),rev_sig(1:1e5));
length(x)
length(rev_sig);
figure(2);
pwelch(rev_sig);
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