module global_variables
       implicit none
       integer row,column
       parameter(row=300,column=500)
       real vxf(row,column),vx(row,column)
	   real vyf(row,column),vy(row,column)
	   real vlbxf(row),vlbx(row),vubyf(column),vuby(column)
       real Txxf(row,column),Txx(row,column),Txxb(row,column)
	   real Tyyf(row,column),Tyy(row,column),Tyyb(row,column)
	   real Txyf(row,column),Txy(row,column),Txyb(row,column)
	   real Tubxyf(column),Tlbxyf(row)
	   real Tubxy(column),Tlbxy(row)
	   real Tc(120000),Trr(120000),FF(120000),Kvv(120000),Tm(120000),vr(120000)
	   real p(row+1,column+1)
	   real Hem(row,column,100),O2(row,column,100),s_Hem(row,column,100),s_O2(row,column,100)
	   integer Hnum(row,column),s_Hnum(row,column),Hnumt

	   real K1,K2,U,w,Q
	   real deltt,deltd,r,A,A_input,f,fL,wx,wy,pp,xs,ys,pi,initial,Khr,Kh1,Kh2,r2,r1,Ktest1,Ktest2


       parameter (deltt=0.5e-7,deltd=0.5e-3)                         !离散时间间隔deltt，离散空间间隔deltd
	   parameter  (r=1e-9,Q=10,r2=1e-6,r1=2e-6)                                       !非线性参数r,衰减相关系数Q
	   real pm,xx,yy(2),x_down,x_up,y_down,y_up,xstart,xend,ystart,yend                            !pm空间密度pm，积分变量xx，上下限yy
       parameter(Kh1=0.75e15,Kh2=1e15)                                          

	   parameter (pi=3.1415926,initial=1e20)

	   

       logical flag,cornerflag,riseflag(row,column),initial_state(row,column),s_riseflag(row,column),s_initial_state(row,column)
       real Kvf(row+1,column+1),Kv(row+1,column+1),Kdf(row+1,column+1),Kd(row+1,column+1),ksf(row+1,column+1),ks(row+1,column+1),K11(row+1,column+1)
	   real B,delt,sw,bb
	   real t,x,y
	   real Fx,Fy,FA

       integer k,kk,kt,kt1,m,n,mm,nn,i,j,h,ii,jj

	   character(3) filename
end module


program  main_timeReverse
     
     use global_variables
     !use once
	 !use twice
	 !use pm_reaction
	 !use pm_density
	 implicit none
	 real Tvf,Tv,Tdf,Td,Tsf,Ts,ktv,ht,Tvb,Tsb
	 
	 
	 
	 
	 real rev_sig(row,10000),integer_sig(row,column)
	 rev_sig=0
	 integer_sig=0
	 do h=2,row-2
	 write(filename,'(I3)') h
	 open(11,file=adjustr(filename)//"B.txt",position="rewind")
	 read(11,*)(rev_sig(h,i),i=1,6000)
	 close(11)
	 end do

     
	 
	  
	  
	 
	 do h=1,1
	  K1=(5.2202+2*4.45)*1e10                !K1,K2,U为材料的模数，K11为破损区域材料模数
	 K11=(5.2202+2*4.45)*1e10

	 K2=5.2202*1e10
	 U=4.45*1e10
     bb=1e-10                                !经典非线性参数bb
	 !bb=0

	 w=2*pi*f
	
     

	 Kvf=K11+K2
	 Kv=K11+K2
     Kdf=K11-K2
     Kd=K11-K2
     Ksf=U
	 Ks=U
	 
	 Tvf=0
	 Tv=0
	 Tdf=0
	 Td=0
	 Tsf=0
	 Ts=0
	 p=7700

	 vxf=0
     vx=0
	 vyf=0
	 vy=0
	 vlbxf=0
	 vlbx=0
	 vubyf=0
	 vuby=0

	 Txxf=0
	 Txx=0
	 Tyyf=0
	 Tyy=0
	 Txyf=0
	 Txy=0
	 Tubxyf=0
	 Tubxy=0
	 Tlbxyf=0
	 Tlbxy=0

	 Tc=0
	 Trr=0
	 FF=0
	 Kvv=0
	 Tm=0
     
	

	 

	 

	 
     
     B=deltd*deltd
     delt=deltt/deltd
     
     flag=.true.

	 
     !O2=-A                                   !O2为记录的打开空间积分路径Po坐标，初始化
     Hem=initial                                  !hem为记录的终点记忆点Pc坐标，初始化
     s_Hem=initial
        
     A=1.5e8                                         !输入信号应力幅值
	  A_input=1e8
	  s_O2=-A
	  s_Hem(:,:,1)=0
	  
	  O2=-A                                !设初始状态为半关闭
	  Hem(:,:,1)=0
	  
	  Hnum=1                                       !记忆点数，初始为1
	  s_Hnum=1
	  
	  f=250e3
	  fL=4e3
	  
      ht=2*4*Q/(deltt*f)
	   ktv=1
	   kk=1
	   kt=1
	   kt1=1
     write(filename,'(I3)') h
	 open(h,file=adjustr(filename)//".txt",position="rewind")
	 open(h+1,file=adjustr(filename)//"a.txt",position="rewind")
	 open(h+2,file=adjustr(filename)//"b.txt",position="rewind")


	 initial_state=.true.
	 s_initial_state=.true.
	 do k=1,3000
	 cornerflag=.false.
	    
        
	    do n=1,column
		    do m=1,row
               if ((n==1).and.(m==1)) then
			   
			   else if ((n==1).and.(m==row)) then

			   else if ((n==column).and.(m==1)) then

			   else if ((n==column).and.(m==row)) then
			   
     

              else if ((n==1).and.(m>1).and.(m<row)) then
                    if (flag) then
                    vxf(m,n)=vx(m,n)+delt*(2/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)+Txy(m,n)-Txy(m-1,n))
                    vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Tlbxy(m)) 
                    vyf(m-1,n)=vy(m-1,n)+delt*(2/(p(m-1,n)+p(m,n)))*(Tyy(m,n)-Tyy(m-1,n)+Txy(m-1,n)-Tlbxy(m-1))   
            
                    vlbxf(m)=vxf(m,n)+K2/K1*(vyf(m,n)-vyf(m-1,n))
            
                    else if (.NOT.flag) then
					

					Tyyf(m,n)=Tyy(m,n)+delt*K1*(vyf(m,n)-vyf(m-1,n))+delt*K2*(vxf(m,n)-vlbxf(m))
					Txyf(m,n)=Txy(m,n)+delt*U*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n))
					
					end if


           
          
           else if (n==column-1.and.m>1.and.m<row)then
                   if (flag) then
					  vxf(m,n)=vxf(m,column)+K2/K1*(vyf(m,n+1)-vyf(m-1,n+1))
                     vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Txy(m,n-1))
                     
                    
                       
					 else if (.NOT.flag) then
					 Txxf(m,n)=Txx(m,n)+delt*K11(m,n)*(vxf(m,n)-vxf(m,n-1))+delt*K2*(vyf(m,n)-vyf(m-1,n))
                     Tyyf(m,n)=Tyy(m,n)+delt*K11(m,n)*(vyf(m,n)-vyf(m-1,n))+delt*K2*(vxf(m,n)-vxf(m,n-1))
                     Txyf(m,n)=Txy(m,n)+delt*U*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n))
					 !!Txxf(m,n)=Txxf(m,n)*exp(-pi*f*deltt/(2*Q))                  !衰减项
                     !!Tyyf(m,n)=Tyyf(m,n)*exp(-pi*f*deltt/(2*Q))
					 !!Txyf(m,n)=Txyf(m,n)*exp(-pi*f*deltt/(2*Q))
					 
					 K11(m,n)=K1*(1+bb*Txxf(m,n))                                 !经典非线性
					 end if
                   
                   
                   
                   
                   
                    
           else if ((n==column).and.(m>1).and.m<row)then
                    if (flag) then
					 
                     vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Txy(m,n-1))
                     
                     vxf(m,n)=rev_sig(m,kt1)
                     
                     
                     
					 else if (.NOT.flag) then
					 Tyyf(m,n)=Tyy(m,n)+delt*K1*(vyf(m,n)-vyf(m-1,n))+delt*K2*(vxf(m,n)-vxf(m,n-1))
					 end if
					 
            

				else if ((m==1).and.((n>1.and.n<column))) then
				     if (flag) then
					 vxf(m,n)=vx(m,n)+delt*(2/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)-Txx(m,n)+Txy(m,n)-Tubxy(n))
                     vxf(m,n-1)=vx(m,n-1)+delt*(2/(p(m,n-1)+p(m,n)))*(Txx(m,n)-Txx(m,n-1)+Txy(m,n-1)-Tubxy(n-1))
                     
					 vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)+Txy(m,n)-Txy(m,n-1))
                     
                     vubyf(n)=vyf(m,n)+K2/K1*(vxf(m,n)-vxf(m,n-1)) 
					 
					 else if (.NOT.flag) then
                     
                     
                     Txxf(m,n)=Txx(m,n)+delt*K1*(vxf(m,n)-vxf(m,n-1))+delt*K2*(vyf(m,n)-vubyf(n))
					 Txyf(m,n)=Txy(m,n)+delt*U*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n))
					 end if

				else if ((m==row).and.(n>1).and.(n<column)) then
				     if (flag) then
					 vxf(m,n)=vx(m,n)+delt*(2/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)-Txx(m,n)+Txy(m,n)-Txy(m-1,n))
                     
					 vxf(m,n-1)=vx(m,n-1)+delt*(2/(p(m,n-1)+p(m,n)))*(Txx(m,n)-Txx(m,n-1)+Txy(m,n-1)-Txy(m-1,n-1))
                    
					 vyf(m,n)=vyf(m-1,n)-K2/K1*(vxf(m,n)-vxf(m,n-1))

					 else if (.NOT.flag) then
                     
					 Txxf(m,n)=Txx(m,n)+delt*K1*(vxf(m,n)-vxf(m,n-1))+delt*K2*(vyf(m,n)-vyf(m-1,n))!+n11*delt*(vxf(m,n)-vx(m,n)-vxf(m,n-1)+vx(m,n-1))/deltt+n12*(vyf(m,n)-vy(m,n)-vyf(m-1,n)+vy(m-1,n))/deltt
					 
					 end if







!==========================================================================
!                     source region
!===========================================================================
                    !中间无源

!===========================================================================================

                 else 
				     if (flag) then
					 vxf(m,n)=vx(m,n)+delt*(2/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)-Txx(m,n)+Txy(m,n)-Txy(m-1,n))
                     vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Txy(m,n-1))
                     !vxf(m,n)=vxf(m,n)*exp(-pi*f*deltt/(2*Q))                   !衰减项
                     !vyf(m,n)=vyf(m,n)*exp(-pi*f*deltt/(2*Q))
                       
					 else if (.NOT.flag) then
					 Txxf(m,n)=Txx(m,n)+delt*K11(m,n)*(vxf(m,n)-vxf(m,n-1))+delt*K2*(vyf(m,n)-vyf(m-1,n))
                     Tyyf(m,n)=Tyy(m,n)+delt*K11(m,n)*(vyf(m,n)-vyf(m-1,n))+delt*K2*(vxf(m,n)-vxf(m,n-1))
                     Txyf(m,n)=Txy(m,n)+delt*U*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n))
                     
					 !Txxf(m,n)=Txxf(m,n)*exp(-pi*f*deltt/(2*Q))                  !衰减项
                     !Tyyf(m,n)=Tyyf(m,n)*exp(-pi*f*deltt/(2*Q))
					 !Txyf(m,n)=Txyf(m,n)*exp(-pi*f*deltt/(2*Q))
					 
					 K11(m,n)=K1*(1+bb*Txxf(m,n))                                 !经典非线性
					 end if
				 end if


				 if (n==column) then
				    cornerflag=.true.
				 end if

				 
!==============================================================================================				 
				 if (flag.and.cornerflag) then
				    mm=1
                    nn=1
                    vxf(mm,nn)=vx(mm,nn)+delt*(2/(p(mm,nn)+p(mm,nn+1)))*(Txx(mm,nn+1)-Txx(mm,nn)+Txy(mm,nn)-Tubxy(nn))
                    vyf(mm,nn)=vy(mm,nn)+delt*(2/(p(mm,nn)+p(mm+1,nn)))*(Tyy(mm+1,nn)-Tyy(mm,nn)+Txy(mm,nn)-Tlbxy(mm))
                
                    vubyf(nn)=vyf(mm,nn+2)-vubyf(nn+2)+vyf(mm,nn)+(2*K2/K1)*(vxf(mm,nn+1)-vxf(mm,nn))
                    vlbxf(mm)=vxf(mm,nn)+vxf(mm+2,nn)-vlbxf(mm+2)+(2*K2/K1)*(vyf(mm+1,nn)-vyf(mm,nn))

!--------------------------------------------------------------------------------------------
                    nn=1
                    mm=row
                    vxf(mm,nn)=vx(mm,nn)+delt*(2/(p(mm,nn)+p(mm,nn+1)))*(Txx(mm,nn+1)+Txy(mm,nn)-Txy(mm-1,nn))
           
                    vyf(mm,nn)=vyf(mm-1,nn)-vyf(mm,nn+2)+vyf(mm-1,nn+2)-(2*K2/K1)*(vxf(mm,nn+1)-vxf(mm,nn))
                    vlbxf(mm)=vxf(mm-2,nn)-vlbxf(mm-2)+vxf(mm,nn)+(2*K2/K1)*(vyf(mm-1,nn)-vyf(mm-2,nn))

!---------------------------------------------------------------------------------------------
                     nn=column
                     mm=1
                     vyf(mm,nn)=vy(mm,nn)+delt*(2/(p(mm,nn)+p(mm+1,nn)))*(Tyy(mm+1,nn)-Tyy(mm,nn)+Txy(mm,nn)-Txy(mm,nn-1))
            
                     vubyf(nn)=vyf(mm,nn-2)-vubyf(nn-2)+vyf(mm,nn)+(2*K2/K1)*(vxf(mm,nn-1)-vxf(mm,nn-2))                     
                     vxf(mm,nn)=vxf(mm,nn-1)-vxf(mm+2,nn)+vxf(mm+2,nn-1)-(2*K2/K1)*(vyf(mm+1,nn)-vyf(mm,nn))
            
!---------------------------------------------------------------------------------------------
                      nn=column
                      mm=row
                      vyf(mm,nn)=vyf(mm-1,nn)-vyf(mm,nn-2)+vyf(mm-1,nn-2)-(2*K2/K1)*(vxf(mm,nn-1)-vxf(mm,nn-2))
                      vxf(mm,nn)=vxf(mm,nn-1)-vxf(mm-2,nn)+vxf(mm-2,nn-1)-(2*K2/K1)*(vyf(mm-1,nn)-vyf(mm-2,nn))



!===============================================================================================
                   else if ((.NOT.flag).and.cornerflag) then
				      mm=1
                      nn=1
                       
                      Tubxyf(nn)=Tubxy(nn)+delt*U*(vxf(mm,nn)+vubyf(nn+1)-vubyf(nn))                                   
                      Tlbxyf(mm)=Tlbxy(mm)+delt*U*(vlbxf(mm+1)-vlbxf(mm)+vyf(mm,nn))

					  Txyf(mm,nn)=Txy(mm,nn)+delt*U*(vxf(mm+1,nn)-vxf(mm,nn)+vyf(mm,nn+1)-vyf(mm,nn))

                      
!-----------------------------------------------------------------------				  
					  nn=1
                      mm=row

					  Txyf(mm,nn)=Txy(mm,nn)+delt*U*(-vxf(mm,nn)+vyf(mm,nn+1)-vyf(mm,nn))

!------------------------------------------------------------------------
                      nn=column
                      mm=1
                      Txyf(mm,nn)=Txy(mm,nn)+delt*U*(vxf(mm+1,nn)-vxf(mm,nn)-vyf(mm,nn))
				   end if
           
		    end do
		  
		end do
		

		if (.Not.flag) then
		vx=vxf
        vy=vyf
		
    
        vuby=vubyf
		

        vlbx=vlbxf
		
	
        
      
       Txxb=Txx
	   Tyyb=Tyy
       Txyb=Txy
       
	   Txx=Txxf
       Tyy=Tyyf
	   
       Txy=Txyf
	  
    
   
       Tubxy=Tubxyf
	  
       
	   Tlbxy=Tlbxyf
	   
	   

	   
    
    
       Kv=Kvf
       Ks=Ksf
	  
     
	 write(*,*) k,Tyy(floor(row/2.0),column)
      

	   !if ((.NOT.flag).and.(k>ht)) then
        !if (.NOT.flag) then
        !Trr(kt)=(Txxf(50,50)+Tyyf(50,50))/sqrt(2.0)
        Trr(kt)=vxf(floor(row/2.0),column)
        Tc(kt)=vyf(floor(row/2.0),column)
		!Tm(kt)=vxf(40,50)
        !vr(ktv)=vyf(40,125)
		!ktv=ktv+1
		kt=kt+1
       end if
       
       if (flag)then
       kt1=kt1+1
       end if
       
	   flag=.NOT.flag
      if (k==1500)then
           write(h+1,"(1x,500E15.7)") ((vxf(i,j),j=1,column),i=1,row)
      elseif(k==2000)then
	       write(h+2,"(1x,500E15.7)") ((vxf(i,j),j=1,column),i=1,row)
	       !elseif(k==802)then
	       !write(h,"(1x,125E15.7)") ((vxf(i,j),j=1,125),i=1,75)
	   endif
	   
	   if(k>=2401.and.k<=2410)then                            !2301-2310,2501-2510,2801-2810,3201-3210
	         
	         integer_sig=integer_sig+array_dot_multi(vxf)
	         
	    endif 
	    
   end do

      !write(h+2,"(1x,E15.7)") Trr
      !write(h+1,"(1x,125E15.7)") ((vxf(i,j),j=1,125),i=1,75)
	 write(h,"(1x,500E15.7)") ((integer_sig(i,j)/10.0,j=1,column),i=1,row) 
	 
	 close(h)
	 close(h+1)
	 close(h+2)
	 
 end do
	 
	
	 write(*, *) "writing complete"

     pause



 
	 close(11)
	 
	 
	 contains
	    function array_dot_multi(x)
	               integer ii,jj
	               real x(300,500),array_dot_multi(300,500)
	               do jj=1,500
	                  do ii=1,300
	                      x(ii,jj)=x(ii,jj)*x(ii,jj)
	                  enddo
	               enddo
	               array_dot_multi=x
	    end function array_dot_multi
end program main_timeReverse