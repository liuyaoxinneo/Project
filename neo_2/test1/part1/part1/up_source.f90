              elseif (m==1.and.n>=100.and.n<=200)then
                    if (flag) then
					 vxf(m,n)=vx(m,n)+delt*(2/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)-Txx(m,n)+Txy(m,n))
                     vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Txy(m,n-1))
                     
                     vxf(m,n)=vxf(m,n)*exp(-pi*f*deltt/(2*Q))                   !衰减项
                     vyf(m,n)=vyf(m,n)*exp(-pi*f*deltt/(2*Q))
                       
					 else if (.NOT.flag) then
					 !t=k*deltt/2
					 Txxf(m,n)=Txx(m,n)+delt*K11(m,n)*(vxf(m,n)-vxf(m,n-1))+delt*K2*(vyf(m,n))
                     Tyyf(m,n)=A_input*sin(2*pi*fL*t)
                     Txyf(m,n)=Txy(m,n)+delt*U*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n))
					 Txxf(m,n)=Txxf(m,n)*exp(-pi*f*deltt/(2*Q))                  !衰减项
                     Tyyf(m,n)=Tyyf(m,n)*exp(-pi*f*deltt/(2*Q))
					 Txyf(m,n)=Txyf(m,n)*exp(-pi*f*deltt/(2*Q))
					 
					! K11(m,n)=K1*(1+bb*Txxf(m,n))                                 !经典非线性
					 end if