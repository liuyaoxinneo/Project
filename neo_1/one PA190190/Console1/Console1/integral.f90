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
		m=m+1
    enddo
    
    once_jifen1=t2
    end function once_jifen1