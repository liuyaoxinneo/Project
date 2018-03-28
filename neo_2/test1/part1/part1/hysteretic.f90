                else if (((n>=180.and.n<=200).and.(m>=180.and.m<=200)).or.((n>=120.and.n<=140).and.(m>=120.and.m<=140))) then                  !ÆÆËğÇøÓò
				
					 if (flag) then
					 vxf(m,n)=vx(m,n)+delt*(2.0/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)-Txx(m,n)+Txy(m,n)-Txy(m-1,n))
                    vyf(m,n)=vy(m,n)+delt*(2.0/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Txy(m,n-1))
                    
                    vxf(m,n)=vxf(m,n)*exp(-pi*f*deltt/(2*Q))                   !Ë¥¼õÏî
                    vyf(m,n)=vyf(m,n)*exp(-pi*f*deltt/(2*Q)) 

					 else if (.NOT.flag) then
					 Tvb=(Txxb(m,n)+Tyyb(m,n))/sqrt(2.0)
                     Tsb=Txyb(m,n)
                     
					 Tv=(Txx(m,n)+Tyy(m,n))/sqrt(2.0)
					 Td=(Txx(m,n)-Tyy(m,n))/sqrt(2.0)
					 Ts=Txy(m,n)
					 Tvf=Tv+delt*(Kv(m,n)/sqrt(2.0))*(vxf(m,n)-vxf(m,n-1)+vyf(m,n)-vyf(m-1,n))                        
                     Tdf=Td+delt*(Kd(m,n)/sqrt(2.0))*(vxf(m,n)-vxf(m,n-1)-vyf(m,n)+vyf(m-1,n))
                     Tsf=Ts+delt*Ks(m,n)*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n))
                      
                       Tvf=Tvf*exp(-pi*f*deltt/(2*Q))                         !Ë¥¼õÏî
                     Tdf=Tdf*exp(-pi*f*deltt/(2*Q))
					 Tsf=Tsf*exp(-pi*f*deltt/(2*Q))

					  Txxf(m,n)=(Tvf+Tdf)/SQRT(2.0)
                      Tyyf(m,n)=(Tvf-Tdf)/SQRT(2.0)
                     Txyf(m,n)=Tsf

                    K11(m,n)=K1*(1+bb*Txxf(m,n))                                 !¾­µä·ÇÏßĞÔ
                    
                    
			       call  pm_caculate(Tv,Tvf,Hem,O2,Hnum,initial_state,m,n,Kh1,Kh2,Kvf,K11,K2,r1,r2,Tvb,riseflag,1,K_test1)
			       !call  pm_caculate(Ts,Tsf,s_Hem,s_O2,s_Hnum,s_initial_state,m,n,Kh1,Kh2,Ksf,K11,K2,r1,r2,Tsb,s_riseflag,2)
					  


					  end if