
!------------------------ 计算了：vx(m,n),vx(m,n-1),vy(m,n),vy(m-1,n);Tyy(m,n),Txy(m,n),Txx(m,n)(仅在not flag时计算) -------
vxf(m,n)=vx(m,n)+delt*(2/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)-Txx(m,n)+Txy(m,n)-Txy(m-1,n))!vx标准迭代公式

vxf(m,n-1)=vxf(m,n)+K2/K1*(vyf(m,n)-vyf(m-1,n)) !成立条件Txx的时间导数=0

vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Txy(m,n-1))!vy标准迭代公式
vyf(m-1,n)=vyf(m,n)+K2/K1*(vxf(m,n)-vxf(m,n-1)) !成立条件Tyy的时间导数=0

Txxf(m,n)=Txx(m,n)+delt*K1*(vxf(m,n)-vxf(m,n-1))+delt*K2*(vyf(m,n)-vyf(m-1,n))!Txx标准迭代公式

Tyyf(m,n)=Tyy(m,n)+delt*K1*(vyf(m,n)-vyf(m-1,n))+delt*K2*(vxf(m,n)-vxf(m,n-1))!Tyy标准迭代公式

Txyf(m,n)=Txy(m,n)+delt*U*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n)) !Txy标准迭代公式

