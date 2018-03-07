module twice
    
	 use pm_density
	 use once
	 implicit none
	 private hh,t1,t2,ep,s1,s2,m,k,n,x,eps,y
	 real hh,t1,t2,ep,s1,s2,x,eps,y
	 integer m,n,k
	 parameter(eps=1e-6)
	 
	 contains
	   
	    real function twice_jifen1(xstart,xend,y_down,sign)           !下限固定，积分变量为y      
		      real xstart,xend,y_down
			  integer sign
			  
			  hh=xend-xstart
			  
			  s1=once_jifen1(xstart,y_down,sign)
			  s2=once_jifen1(xend,y_down,sign)
			  t1=0.5*hh*(s1+s2)
			  ep=eps+1
			  m=1
			  n=1
			  do while(ep>eps.and.m<9)
                 t2=0.5*t1
				 k=1
				 do while(k<=n)
                    x=xstart+(k-0.5)*hh
					t2=t2+0.5*hh*once_jifen1(x,y_down,sign)
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
			  t1=t2
			  m=m+1
			  enddo
			  twice_jifen1=t2
		end function twice_jifen1
		 
		 
		 real function twice_jifen2(xstart,xend,y_up,sign)          !上限固定，积分变量为y
		      
			  real xstart,xend,y_up
			  integer sign
			  hh=xend-xstart
			  
			  s1=once_jifen2(xstart,y_up,sign)
			  s2=once_jifen2(xend,y_up,sign)
			  t1=0.5*hh*(s1+s2)
			  ep=eps+1
			  m=1
			  n=1
			  do while(ep>eps.and.m<9)
                 t2=0.5*t1
				 k=1
				 do while(k<=n)
                    x=xstart+(k-0.5)*hh
					t2=t2+0.5*hh*once_jifen2(x,y_up,sign)
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
			  t1=t2
			  m=m+1
			  enddo
			  twice_jifen2=t2
		end function twice_jifen2

		real function twice_jifen3(ystart,yend,x_down,sign)           !下限固定，积分变量为x      
		      real ystart,yend,x_down
			  integer sign
			  hh=yend-ystart
			  
			  s1=once_jifen3(x_down,ystart,sign)
			  s2=once_jifen3(x_down,yend,sign)
			  t1=0.5*hh*(s1+s2)
			  ep=eps+1
			  m=1
			  n=1
			  do while(ep>eps.and.m<9)
                 t2=0.5*t1
				 k=1
				 do while(k<=n)
                    y=ystart+(k-0.5)*hh
					t2=t2+0.5*hh*once_jifen3(x_down,y,sign)
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
			  t1=t2
			  m=m+1
			  enddo
			  twice_jifen3=t2
		end function twice_jifen3
        

		real function twice_jifen4(x_up,ystart,yend,sign)          !上限固定，积分变量为x
		      real ystart,yend,x_up
			  integer sign
			  hh=yend-ystart
			  
			  s1=once_jifen4(x_up,ystart,sign)
			  s2=once_jifen4(x_up,yend,sign)
			  t1=0.5*hh*(s1+s2)
			  ep=eps+1
			  m=1
			  n=1
			  do while(ep>eps.and.m<9)
                 t2=0.5*t1
				 k=1
				 do while(k<=n)
                    y=ystart+(k-0.5)*hh
					t2=t2+0.5*hh*once_jifen4(x_up,y,sign)
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
			  t1=t2
			  m=m+1
			  enddo
			  twice_jifen4=t2
		end function twice_jifen4
        

end module
