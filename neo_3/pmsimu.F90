module global_variables
       implicit none
       integer,parameter::row=300,column=500
       
       real vxf(row,column),vx(row,column)
	   real vyf(row,column),vy(row,column)
	   real vlbxf(row),vlbx(row),vubyf(column),vuby(column)
       real Txxf(row,column),Txx(row,column),Txxb(row,column)
	   real Tyyf(row,column),Tyy(row,column),Tyyb(row,column)
	   real Txyf(row,column),Txy(row,column),Txyb(row,column)
	   real Tubxyf(column),Tlbxyf(row)
	   real Tubxy(column),Tlbxy(row)
	   real Tc(200000),Trr(200000),FF(200000),Kvv(200000),Tm(200000),vr(200000),v_recsig_x(row,200000),v_recsig_y(row,200000)
	   real p(row+1,column+1)
	   real Hem(row,column,100),O2(row,column,100),s_Hem(row,column,100),s_O2(row,column,100)
	   integer Hnum(row,column),s_Hnum(row,column),Hnumt
       
       !全局变量
       common vxf,vx
	   common vyf,vy
	   common vlbxf,vlbx,vubyf,vuby
       common Txxf,Txx,Txxb
	   common Tyyf,Tyy,Tyyb
	   common Txyf,Txy,Txyb
	   common Tubxyf,Tlbxyf
	   common Tubxy,Tlbxy
	   common Tc,Trr,FF,Kvv,Tm,vr,v_recsig_x,v_recsig_y
	   common p
	   common Hem,O2,s_Hem,s_O2
	   common Hnum,s_Hnum,Hnumt
       
       
	   real K1,K2,U,w,Q
	   real deltt,deltd,r,A,A_input,f,fL,wx,wy,pp,xs,ys,pi,initial,Khr,Kh1,Kh2,r2,r1,K_test1(row,column),Ktest2,fatten


       parameter (deltt=5e-8,deltd=5e-4)  !离散时间间隔deltt，离散空间间隔deltd
	   parameter  (r=1e-9,Q=2,r2=1e-3,r1=2e-3)  !非线性参数r,衰减相关系数Q
	   real pm,xx,yy(2),x_down,x_up,y_down,y_up,xstart,xend,ystart,yend   !pm空间密度pm，积分变量xx，上下限yy
       parameter(Kh1=0.75e15,Kh2=1e15)                                          

	   parameter (pi=3.1415926,initial=1e20)

	   

       logical flag,cornerflag,riseflag(row,column),initial_state(row,column),s_riseflag(row,column),s_initial_state(row,column)
       real Kvf(row+1,column+1),Kv(row+1,column+1),Kdf(row+1,column+1),Kd(row+1,column+1),ksf(row+1,column+1),ks(row+1,column+1),K11(row+1,column+1)
	   real B,delt,sw,bb
	   real t,x,y
	   real Fx,Fy,FA

       integer k,kk,kt,kt1,m,n,mm,nn,i,j,h,ii,jj,sig_row

	   character(3) filename
end module


program  main_simu
     
     use global_variables
     use once
	 use twice
	 use pm_reaction_rewrite
	 use pm_density
	 implicit none
	 real Tvf,Tv,Tdf,Td,Tsf,Ts,ktv,ht,Tvb,Tsb
     
     integer remark !交替计算应力应变的因子
     
!metal 
	 K1=(5.2202+2*4.45)*1e10 !K1,K2,U为材料的模数
	 K11=(5.2202+2*4.45)*1e10 !K11为破损区域材料模数
	 K2=5.2202*1e10
	 U=4.45*1e10
	 p = 7700 !密度，kg/m^3
 
     bb=1e-10 !经典非线性参数bb
	 !bb=0

    !f=700e3 !输入外力的频率
     f = 8.5e3 !声源频率             
    !f=4e3
    !fL=4e3
    !fL=8e3
	 w=2*pi*f
	
	 Kvf=K11+K2
	 Kv=K11+K2
     Kdf=K11-K2
     Kd=K11-K2
     Ksf=2*U !根据参考文献应该为2U？？
	 Ks=2*U
	 
	 Tvf=0
	 Tv=0
	 Tdf=0
	 Td=0
	 Tsf=0
	 Ts=0
	    
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
       
     B=deltd*deltd !B:单个网格的面积
     delt=deltt/deltd
     
     !O2=-A                                   !O2为记录的打开空间积分路径Po坐标，初始化
     Hem=initial                                  !hem为记录的终点记忆点Pc坐标，初始化
     s_Hem=initial
        
      A=1.5e8  
	  A_input=1e8 !输入外力幅值
	  s_O2=-A
	  s_Hem(:,:,1)=0
	  
	  O2=-A  !设初始状态为半关闭
	  Hem(:,:,1)=0
	  
	  Hnum=1 !记忆点数，初始为1
	  s_Hnum=1
	  	  
      ht=2*4*Q/(deltt*f)
      ktv=1
	  kk=1
	  kt=1
	  kt1=1
        
	 initial_state=.true.
	 s_initial_state=.true.
	 riseflag=.true.
     
    flag=.true. !用以交替计算速度和应力
     
	 do k=1,200000 !k代表时间步数
         
         cornerflag=.false.
         
         do remark=1,2
             do n=1,column
                 do m=1,row
                     
            !排除区域角落的4个单元
               if ((n==1).and.(m==1)) then
			   
			   else if ((n==1).and.(m==row)) then

			   else if ((n==column).and.(m==1)) then

               else if ((n==column).and.(m==row)) then
			                         
!------------------------ 在紧靠源（左端，5mm~145mm）的一列单元上，正对源的单元 -------
!------------------------ 20180329：按照论文，在紧靠源（左端，50mm~100mm）的一列单元上，正对源的单元 -------       
!------------------------ 计算了：vx(m,n),vx(m,n-1),vy(m,n),vy(m-1,n);Tyy(m,n),Txy(m,n),Txx(m,n)(仅在not flag时计算) -------
               else if ((n==1).and.(m>=100).and.(m<=row-100)) then
                   if ((flag).and.(k<=150000)) then
                         vxf(m,n)=vx(m,n)+delt*(2/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)-Txx(m,n)+Txy(m,n)-Txy(m-1,n))!vx标准迭代公式
                         vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Tlbxy(m)) ! Tlbxy(m)存储的是边界上的条件
                         vyf(m-1,n)=vy(m-1,n)+delt*(2/(p(m-1,n)+p(m,n)))*(Tyy(m,n)-Tyy(m-1,n)+Txy(m-1,n)-Tlbxy(m-1)) !和上一行一样，唯一的区别在于用m-1替代m                         
                         vlbxf(m)=vxf(m,n)+K2/K1*(vyf(m,n)-vyf(m-1,n)) !成立条件：Txx时间导数为0,成立，因为没有进入第二个循环时，Txx=Txxf=0

                    else if ((.NOT.flag).and.(k<=150000)) then 
                         t=k*deltt/2
                         Fx=A_input*sin(2*pi*f*t)                                                              !输入信号Fx分量
                         Fy=0
                         Txxf(m,n)=Fx
                         Tyyf(m,n)=Tyy(m,n)+delt*K1*(vyf(m,n)-vyf(m-1,n))+delt*K2*(vxf(m,n)-vlbxf(m))!Tyy标准迭代公式
                         Txyf(m,n)=Fy !Txy(m,1)=0
                         
                    else if ((flag).and.(k>150000)) then
                         vxf(m,n)=vx(m,n)+delt*(2/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)+Txy(m,n)-Txy(m-1,n)) !成立条件：-Txx(m,n)=0
                         vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Tlbxy(m)) 
                         vyf(m-1,n)=vy(m-1,n)+delt*(2/(p(m-1,n)+p(m,n)))*(Tyy(m,n)-Tyy(m-1,n)+Txy(m-1,n)-Tlbxy(m-1))
                         vlbxf(m)=vxf(m,n)+K2/K1*(vyf(m,n)-vyf(m-1,n)) !成立条件：Txx时间导数为0，成立，因为没有进入第二个循环时，Txx=Txxf=0
            
                    else if ((.NOT.flag).and.(k>150000)) then 
                         Tyyf(m,n)=Tyy(m,n)+delt*K1*(vyf(m,n)-vyf(m-1,n))+delt*K2*(vxf(m,n)-vlbxf(m)) !Tyy标准迭代公式
                         Txyf(m,n)=Txy(m,n)+delt*U*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n)) !Txy标准迭代公式,使用了条件：U为常数
                    end if
                     
!------------------------ 在紧靠源（左端，5mm~145mm）的一列单元上，不正对源的单元 -------
                     !----- 上半部分(除角落)：n=1,m=2~9 -----
                     !------------------------ 计算了：vy(m,n);Tyy(m,n),Txy(m,n)) -------
               else if ((n==1).and.(m>1).and.(m<100)) then
                   if (flag) then
                        vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Tlbxy(m)) 
                   else if (.NOT.flag) then
                        Tyyf(m,n)=Tyy(m,n)+delt*K1*(vyf(m,n)-vyf(m-1,n))+delt*K2*(vxf(m,n)-vlbxf(m)) !Tyy标准迭代公式
                        Txyf(m,n)=Txy(m,n)+delt*U*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n)) !Txy标准迭代公式
                   end if  
                     !----- 下半部分(除角落)：n=1,m=291~299 -----
                     !------------------------ 计算了：vy(m,n);Tyy(m,n),Txy(m,n)) -------
               else if ((n==1).and.(m>row-100).and.(m<row)) then
                    if (flag) then
                        vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Tlbxy(m)) 
                    else if (.NOT.flag) then
                        Tyyf(m,n)=Tyy(m,n)+delt*K1*(vyf(m,n)-vyf(m-1,n))+delt*K2*(vxf(m,n)-vlbxf(m))
                        Txyf(m,n)=Txy(m,n)+delt*U*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n))
                    end if

!------------------------ 在最右端的一列单元上：n=500,m=2~299 -------
!------------------------ 计算了：vx(m,n),vx(m,n-1),vy(m,n),vy(m-1,n);Tyy(m,n) -------
               else if ((n==column).and.(m>1).and.(m<row)) then
                   if (flag) then
                       vxf(m,n-1)=vx(m,n-1)+delt*(2/(p(m,n-1)+p(m,n)))*(Txx(m,n)-Txx(m,n-1)+Txy(m,n-1)-Txy(m-1,n-1)) !(m,n-1)的vxf
                       vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Txy(m,n-1)) !vy标准迭代公式
                       vyf(m-1,n)=vy(m-1,n)+delt*(2/(p(m-1,n)+p(m,n)))*(Tyy(m,n)-Tyy(m-1,n)+Txy(m-1,n)-Txy(m-1,n-1)) !(m-1,n)的vyf
                       vxf(m,n)=vxf(m,n-1)-K2/K1*(vyf(m,n)-vyf(m-1,n)) !成立条件：Txx时间导数为0
                       !以上四行完成了(m,n)单元的四条边上的速度
                   else if (.NOT.flag) then
                       Tyyf(m,n)=Tyy(m,n)+delt*K1*(vyf(m,n)-vyf(m-1,n))+delt*K2*(vxf(m,n)-vxf(m,n-1)) !y方向本构方程
                   end if

!------------------------ 在最上端的一行单元上：n=2~99，201~499,m=1 -------可能是因为在50~100mm上安放着第二个声源？？
!------先将上方的第二个源去除！，改为n=2~499,m=1
!------------------------ 计算了：vx(m,n),vx(m,n-1),vy(m,n),vy(m-1,n);Txx(m,n),Txy(m,n) -------
               else if ((m==1).and.((n>1.and.n<column))) then
                   if (flag) then
                       vxf(m,n)=vx(m,n)+delt*(2/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)-Txx(m,n)+Txy(m,n)-Tubxy(n))
                       vxf(m,n-1)=vx(m,n-1)+delt*(2/(p(m,n-1)+p(m,n)))*(Txx(m,n)-Txx(m,n-1)+Txy(m,n-1)-Tubxy(n-1))
                       vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)+Txy(m,n)-Txy(m,n-1)) !成立条件：Tyy(m,n)=0
                       vubyf(n)=vyf(m,n)+K2/K1*(vxf(m,n)-vxf(m,n-1)) !成立条件：Tyy时间导数为0
                   else if (.NOT.flag) then
                       Txxf(m,n)=Txx(m,n)+delt*K1*(vxf(m,n)-vxf(m,n-1))+delt*K2*(vyf(m,n)-vubyf(n))
					   Txyf(m,n)=Txy(m,n)+delt*U*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n)) !缺少Tyy(m,n)的计算，也可能默认为0
                   end if

!------------------------ 在最下端的一行单元上：n=2~499,m=300 -------
!------------------------ 计算了：vx(m,n),vx(m,n-1),vy(m,n);Txx(m,n) -------
               else if ((m==row).and.(n>1).and.(n<column)) then
                   if (flag) then
                       vxf(m,n)=vx(m,n)+delt*(2/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)-Txx(m,n)+Txy(m,n)-Txy(m-1,n))!vx标准迭代公式
                       vxf(m,n-1)=vx(m,n-1)+delt*(2/(p(m,n-1)+p(m,n)))*(Txx(m,n)-Txx(m,n-1)+Txy(m,n-1)-Txy(m-1,n-1))!vx标准迭代公式
                       vyf(m,n)=vyf(m-1,n)-K2/K1*(vxf(m,n)-vxf(m,n-1)) !成立条件Tyy的时间导数=0
                   else if (.NOT.flag) then
                       Txxf(m,n)=Txx(m,n)+delt*K1*(vxf(m,n)-vxf(m,n-1))+delt*K2*(vyf(m,n)-vyf(m-1,n))!Txx标准迭代公式
                   end if

!------------------------ 在最上端的一行单元上：n=100~200,m=1 -------
!------------------------ 计算了：vx(m,n),vy(m,n);Txx(m,n),Txy(m,n),Tyy(m,n);考虑了衰减：对速度、应力 -------

                   
                   !将此区域先删除！2018.3.28

                   
!===========================================================================
!                         hysteretic region
!===========================================================================
!------------------------ 介质内的两个正方形区域：n=180~200,m=180~200 & n=280~？？ -------
!------------------------ 计算了：vx(m,n),vy(m,n);Txx(m,n),Txy(m,n),Tyy(m,n);考虑了衰减：对速度、应力 -------

                   
                   !将此区域先删除！2018.3.28


!------------------------ 介质内除角落、边界外的线性区域 -------
!------------------------ 计算了：vx(m,n),vy(m,n);Txx(m,n),Txy(m,n),Tyy(m,n);考虑了衰减：对速度、应力 -------
               else 
                   if (flag) then
                       vxf(m,n)=vx(m,n)+delt*(2/(p(m,n)+p(m,n+1)))*(Txx(m,n+1)-Txx(m,n)+Txy(m,n)-Txy(m-1,n)) !vx标准迭代公式
                       vyf(m,n)=vy(m,n)+delt*(2/(p(m,n)+p(m+1,n)))*(Tyy(m+1,n)-Tyy(m,n)+Txy(m,n)-Txy(m,n-1)) !vy标准迭代公式
                       vxf(m,n)=vxf(m,n)*exp(-pi*f*deltt/(2*Q)) !衰减项
                       vyf(m,n)=vyf(m,n)*exp(-pi*f*deltt/(2*Q))
                   else if (.NOT.flag) then
                       Txxf(m,n)=Txx(m,n)+delt*K11(m,n)*(vxf(m,n)-vxf(m,n-1))+delt*K2*(vyf(m,n)-vyf(m-1,n)) !此处用K11代替了标准迭代公式中的K1
                       Tyyf(m,n)=Tyy(m,n)+delt*K11(m,n)*(vyf(m,n)-vyf(m-1,n))+delt*K2*(vxf(m,n)-vxf(m,n-1)) !Tyy标准迭代公式
                       Txyf(m,n)=Txy(m,n)+delt*U*(vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n))
                       Txxf(m,n)=Txxf(m,n)*exp(-pi*f*deltt/(2*Q))!衰减项
                       Tyyf(m,n)=Tyyf(m,n)*exp(-pi*f*deltt/(2*Q))
					   Txyf(m,n)=Txyf(m,n)*exp(-pi*f*deltt/(2*Q))

					   K11(m,n)=K1*(1+bb*Txxf(m,n)) !经典非线性
                   end if
               end if
               
               if (n==column) then
				    cornerflag=.true. !即，对右侧的两个角落上的单元来说
               end if
!------------------------ 以下计算四个角落单元 -------	                  
               if (flag.and.cornerflag) then
!------------------------ 单元(1,1) -------                   
                   mm=1
                   nn=1
                   vxf(mm,nn)=vx(mm,nn)+delt*(2/(p(mm,nn)+p(mm,nn+1)))*(Txx(mm,nn+1)-Txx(mm,nn)+Txy(mm,nn)-Tubxy(nn))
                   vyf(mm,nn)=vy(mm,nn)+delt*(2/(p(mm,nn)+p(mm+1,nn)))*(Tyy(mm+1,nn)-Tyy(mm,nn)+Txy(mm,nn)-Tlbxy(mm))
                   vubyf(nn)=vyf(mm,nn+2)-vubyf(nn+2)+vyf(mm,nn)+(2*K2/K1)*(vxf(mm,nn+1)-vxf(mm,nn)) !成立条件：Tyy的时间导数为0 && 基于x、y方向速度差线性变化的设想
                   vlbxf(mm)=vxf(mm,nn)+vxf(mm+2,nn)-vlbxf(mm+2)+(2*K2/K1)*(vyf(mm+1,nn)-vyf(mm,nn))
!------------------------ 单元(300,1) -------                    
                   nn=1
                   mm=row
                   vxf(mm,nn)=vx(mm,nn)+delt*(2/(p(mm,nn)+p(mm,nn+1)))*(Txx(mm,nn+1)+Txy(mm,nn)-Txy(mm-1,nn)) !成立条件：Txx(m,n)=0                 
                   vyf(mm,nn)=vyf(mm-1,nn)-vyf(mm,nn+2)+vyf(mm-1,nn+2)-(2*K2/K1)*(vxf(mm,nn+1)-vxf(mm,nn))
                   vlbxf(mm)=vxf(mm-2,nn)-vlbxf(mm-2)+vxf(mm,nn)+(2*K2/K1)*(vyf(mm-1,nn)-vyf(mm-2,nn)) 
!------------------------ 单元(1,500) ------- 
                   nn=column
                   mm=1
                   vyf(mm,nn)=vy(mm,nn)+delt*(2/(p(mm,nn)+p(mm+1,nn)))*(Tyy(mm+1,nn)-Tyy(mm,nn)+Txy(mm,nn)-Txy(mm,nn-1)) !vy标准迭代公式            
                   vubyf(nn)=vyf(mm,nn-2)-vubyf(nn-2)+vyf(mm,nn)+(2*K2/K1)*(vxf(mm,nn-1)-vxf(mm,nn-2))                     
                   vxf(mm,nn)=vxf(mm,nn-1)-vxf(mm+2,nn)+vxf(mm+2,nn-1)-(2*K2/K1)*(vyf(mm+1,nn)-vyf(mm,nn)) 
!------------------------ 单元(300,500) -------                  
                   nn=column
                   mm=row
                   vyf(mm,nn)=vyf(mm-1,nn)-vyf(mm,nn-2)+vyf(mm-1,nn-2)-(2*K2/K1)*(vxf(mm,nn-1)-vxf(mm,nn-2))
                   vxf(mm,nn)=vxf(mm,nn-1)-vxf(mm-2,nn)+vxf(mm-2,nn-1)-(2*K2/K1)*(vyf(mm-1,nn)-vyf(mm-2,nn))
               else if ((.NOT.flag).and.cornerflag) then
!------------------------ 单元(1,1)的应力：计算了3/4 -------
                   mm=1
                   nn=1  
                   Tubxyf(nn)=Tubxy(nn)+delt*U*(vxf(mm,nn)+vubyf(nn+1)-vubyf(nn)) !vx(1,0)=0                           
                   Tlbxyf(mm)=Tlbxy(mm)+delt*U*(vlbxf(mm+1)-vlbxf(mm)+vyf(mm,nn)) !vy(1,0)=0
                   Txyf(mm,nn)=Txy(mm,nn)+delt*U*(vxf(mm+1,nn)-vxf(mm,nn)+vyf(mm,nn+1)-vyf(mm,nn)) 
!------------------------ 单元(300,1)的应力：直接幅值1/4 -------                   
                   nn=1
                   mm=row
                   Txyf(mm,nn)=Txy(mm,nn)+delt*U*(-vxf(mm,nn)+vyf(mm,nn+1)-vyf(mm,nn)) !成立条件：vx(301,1)=0
                   !源程序用的是下一行代码给Txy赋值
                   !Txyf(mm,nn) = 0
!------------------------ 单元(1,500)的应力：计算了1/4 ------- 
                   nn=column
                   mm=1
                   Txyf(mm,nn)=Txy(mm,nn)+delt*U*(vxf(mm+1,nn)-vxf(mm,nn)-vyf(mm,nn)) !成立条件：vy(1,501)=0
               end if
               
                 end do !m循环的end do
             end do !n循环的end do
              
             !write(*,*) k,remark,flag
             flag=.NOT.flag        
         end do !remark循环的end do
         remark=1
                
!------------------------ 以下为迭代(数值传递)过程 ----------------
        
!        if (.NOT.flag) then
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
            
            !write(*,*) k,Kvf(180,180),Tvf,K_test1(180,180),riseflag(180,180),initial_state(180,180),Hnum(180,180) !在控制窗口显示(180,180)单元的状态参数
            !write(*,*) "传值结束！"
            
            !Trr(kt)=vxf(100,column) !根据初始化，kt从1开始,即，Trr(kt)=vxf(100,500)

            !把最右一列除角落之外的单元(2~299,500)速度存储在v_recsig_x/y中
            do sig_row=2,row-1
                v_recsig_y(sig_row,kt)=vyf(sig_row,column)
                v_recsig_x(sig_row,kt)=vxf(sig_row,column)
            end do
            
            write(*,*) "STEP:",k,"TIME:",k*5e-8,"VX：",v_recsig_x(150,kt)
            
            !Tm(kt)=vxf(40,50)
            !vr(ktv)=vyf(40,125)
            !ktv=ktv+1
		    kt=kt+1
            !write(*,*) k,kt
!        end if
        
        !flag=.NOT.flag
        
     end do !k的end do

!==========================    
!------输出文件部分-----
!==========================
       do sig_row=2,row-2 !sig_row r=(2,298)
           write(filename,'(I3)') sig_row !以3个字符宽度把sig_row的数字r存入filename内
           open(sig_row+100,file=adjustr(filename)//"x_rev.txt",position="rewind") !打开代码为r+100，名字为“rx_rev.txt”的文件，并且把读写位置放置在文件头
           !open(sig_row+101,file=adjustr(filename)//"y_rev.txt",position="rewind") 
           write(sig_row+100,"(1x,E15.7)") (v_recsig_x(sig_row,i),i=1,200000)
           !在代码为r+100的文件中开始写文件
           !将读写位置向右移动1位，以科学记数法输出，15个字符宽度，其中小数部分占7位
           !输出的数据为v_recsig_x内(r,1~20000)
           !write(sig_row+101,"(1x,E15.7)") (v_recsig_y(sig_row,i),i=1,200000)
           close(sig_row+100) !写完后，将代码为r+100的文件关闭
           !close(sig_row+101)
       end do	  
	  !write(h+40,"(1x,E15.7)") (O2(68,54,i),i=1,500)
	 !write(*,*) h,Hnum(50,50)
      
     pause
	 write(*, *) "writing complete"
end program main_simu