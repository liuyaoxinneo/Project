PW1=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\Plane wave\\TR one PW\\Console1\\Console1\\1.txt');
PW2=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\Plane wave\\TR two PW\\Console1\\Console1\\1.txt');

PA190190=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\TR one PA190190\\Console1\\Console1\\1.txt');
PA190150=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\TR one PA190150\\Console1\\Console1\\1.txt');
PA100190=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\TR one PA100190\\Console1\\Console1\\1.txt');

PA1601902=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\TR two PA160190\\Console1\\Console1\\1.txt');
PA1751902=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\TR two PA175190\\Console1\\Console1\\1.txt');
PA1901902=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\TR two PA190190\\Console1\\Console1\\1.txt');

PA190290190=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\TR two2 PA190290190\\Console1\\Console1\\1.txt');
PA190290240=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\TR two2 PA190290240\\Console1\\Console1\\1.txt');
PA190290290=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\TR two2 PA190290290\\Console1\\Console1\\1.txt');

PA190270190=sprintf('C:\\Users\\Kyle\\Documents\\Personal\\NDT\\Nonlinear Program\\TR two3 PA190270190\\Console1\\Console1\\1.txt');

x=[1 250];y=[1 150];
%[z]=load(PW1);image(x,y,z.^2*2e4);xlabel('x/mm');ylabel('y/mm');
%[z]=load(PW2);image(x,y,z.^2*3e11);xlabel('x/mm');ylabel('y/mm');

%[z]=load(PA190190);image(x,y,z.^2*4e15);xlabel('x/mm');ylabel('y/mm');
%[z]=load(PA190150);image(x,y,z.^2*4e15);xlabel('x/mm');ylabel('y/mm');
%[z]=load(PA100190);image(x,y,z.^6*6e42);xlabel('x/mm');ylabel('y/mm');

%[z]=load(PA1601902);image(x,y,z.^2*4e15);xlabel('x/mm');ylabel('y/mm');
%[z]=load(PA1751902);image(x,y,z.^2*2e15);xlabel('x/mm');ylabel('y/mm');
%[z]=load(PA1901902);image(x,y,z.^2*4e15);xlabel('x/mm');ylabel('y/mm');

%[z1]=load(PA1601902);[z2]=load(PA1901902);[z]=[z1]+[z2];image(x,y,z.^2*4e15);xlabel('x/mm');ylabel('y/mm');

%[z]=load(PA190290190);image(x,y,z.^2*5e15);xlabel('x/mm');ylabel('y/mm');
%[z]=load(PA190290240);image(x,y,z.^2*2e15);xlabel('x/mm');ylabel('y/mm');
%[z]=load(PA190290290);image(x,y,z.^2*6e15);xlabel('x/mm');ylabel('y/mm');

[z1]=load(PA190290190);[z2]=load(PA190290290);[z]=[z1]+[z2];image(x,y,z.^2*4e15);xlabel('x/mm');ylabel('y/mm');

%[z]=load(PA190270190);image(x,y,z.^2*2e17);xlabel('x/mm');ylabel('y/mm');



