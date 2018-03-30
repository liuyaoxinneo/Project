module once
      
	  use pm_density
      
	  private ep,m,n,hh,t1,t2,k,yy,h,eps
	  real ep,once_result,hh,t1,t2,yy,eps,xx
	  integer m,n,k
	  parameter(eps=1e-6)                                 !模块里不能给变量赋值，只能参数赋值？
	   
    contains
    
    !数值积分方法：梯形递归公式？？
    real function once_jifen1(x,y_down,sign)                 !下限固定，积分变量y
	real x,y(2),y_down
	integer sign
    y(1)=y_down
    y(2)=x	
	m=1
	n=1                                               
    hh=y(2)-y(1)
	t1=0.5*hh*(density(x,y(1),sign)+density(x,y(2),sign))          !
	ep=eps+1.0
    
    do while(ep>eps.and.m<=9)                           !规定最大划分次数m，避免死循环
        t2=0.5*t1
		k=1
        do while(k<=n)                                   !至少划分一次
            yy=y(1)+(k-0.5)*hh
			t2=t2+0.5*hh*density(x,yy,sign)
			k=k+1
        enddo
        if (t1==0)then
            ep=abs(t2-t1)
        elseif(t1==t2)then
            ep=0
        else
            ep=abs(1-t2/t1) 
        endif
        hh=hh/2
		n=n*2
		t2=t1
        !t1=t2
		m=m+1
    enddo
    
    once_jifen1=t2
    end function once_jifen1
      
    real function once_jifen2(x,y_up,sign)                      !上限固定,y变量
    real x,y(2),y_up
    integer sign
    m=1
	n=1
	y(1)=x
	y(2)=y_up
    hh=y(2)-y(1)
    t1=0.5*hh*(density(x,y(1),sign)+density(x,y(2),sign))          !
    ep=eps+1.0
    
    do while(ep>eps.and.m<=9)                           !规定最大划分次数m，避免死循环
        t2=0.5*t1
		k=1
        do while(k<=n)                                   !至少划分一次
            yy=y(1)+(k-0.5)*hh
			t2=t2+0.5*hh*density(x,yy,sign)
			k=k+1
        enddo
        if (t1==0)then
            ep=abs(t2-t1)
        elseif(t1==t2)then
		    ep=0
        else 
            ep=abs(1-t2/t1)
        endif
        hh=hh/2
		n=n*2
		t2=t1
        !t1=t2
		m=m+1
    enddo
    
    once_jifen2=t2  
    end function once_jifen2
    
    real function once_jifen3(x_down,y,sign)                 !下限固定,积分变量为x
    real x(2),x_down
    integer sign
    x(1)=x_down
    x(2)=y
	m=1
	n=1
    hh=x(2)-x(1)
    t1=0.5*hh*(density(x(1),y,sign)+density(x(2),y,sign))          !
	ep=eps+1.0
    
    do while(ep>eps.and.m<=9)                           !规定最大划分次数m，避免死循环
        t2=0.5*t1
		k=1
        do while(k<=n)                                   !至少划分一次
            xx=x(1)+(k-0.5)*hh
			t2=t2+0.5*hh*density(xx,y,sign)
			k=k+1
        enddo
        if (t1==0)then
            ep=abs(t2-t1)
        elseif(t1==t2)then
		    ep=0
        else 
		    ep=abs(1-t2/t1)
        endif
        hh=hh/2
		n=n*2
		t2=t1
        !t1=t2
		m=m+1
    enddo
      
    once_jifen3=t2 
    end function once_jifen3
    
    real function once_jifen4(x_up,y,sign)                      !上限固定,积分变量为x
    real y,x_up,x(2)
	integer sign
	m=1
	n=1
	x(1)=y
	x(2)=x_up
    hh=x(2)-x(1)
    t1=0.5*hh*(density(x(1),y,sign)+density(x(2),y,sign))          !
	ep=eps+1.0
    
    do while(ep>eps.and.m<=9)                           !规定最大划分次数m，避免死循环
        t2=0.5*t1
		k=1
        do while(k<=n)                                   !至少划分一次
            xx=x(1)+(k-0.5)*hh
			t2=t2+0.5*hh*density(xx,y,sign)
			k=k+1
        enddo
        if (t1==0)then
            ep=abs(t2-t1)
        elseif(t1==t2)then
		    ep=0
        else 
            ep=abs(1-t2/t1)
        endif
        hh=hh/2
		n=n*2
		t2=t1
        !t1=t2
		m=m+1
    enddo
    
    once_jifen4=t2
    end function once_jifen4


end module
