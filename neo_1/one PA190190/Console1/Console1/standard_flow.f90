
!------------------------ �����ˣ�vx(m,n),vx(m,n-1),vy(m,n),vy(m-1,n);Tyy(m,n),Txy(m,n),Txx(m,n)(����not flagʱ����) -------
vxf(m,n)=vx(m,n)+delt*(2/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)-Txx(m,n)+Txy(m,n)-Txy(m-1,n))!vx��׼������ʽ

vxf(m,n-1)=vxf(m,n)+K2/K1*(vyf(m,n)-vyf(m-1,n)) !��������Txx��ʱ�䵼��=0

vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Txy(m,n-1))!vy��׼������ʽ
vyf(m-1,n)=vyf(m,n)+K2/K1*(vxf(m,n)-vxf(m,n-1)) !��������Tyy��ʱ�䵼��=0

Txxf(m,n)=Txx(m,n)+delt*K1*(vxf(m,n)-vxf(m,n-1))+delt*K2*(vyf(m,n)-vyf(m-1,n))!Txx��׼������ʽ

Tyyf(m,n)=Tyy(m,n)+delt*K1*(vyf(m,n)-vyf(m-1,n))+delt*K2*(vxf(m,n)-vxf(m,n-1))!Tyy��׼������ʽ

Txyf(m,n)=Txy(m,n)+delt*U*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n)) !Txy��׼������ʽ

