module pm_reaction_rewrite
  use once
  use twice
  
   private Hnumt,ii,Khr,x_down,A,x,xstart,xend,y,y_up,y_down,K_step,Hem_transition,O2_transition
  real Hnumt,Khr,x_down,A,x,xstart,xend,y,y_up,y_down,K_step,Hem_transition,O2_transition,Khr_line
   integer ii,column_i,row_i
  
  parameter(row_i=300,column_i=500)
  
  parameter(A=1.5e8)
  contains
  subroutine pm_caculate(Tv,Tvf,Hem,O2,Hnum,initial_state,m,n,Kh1,Kh2,Kvf,K11,K2,r1,r2,Tvb,riseflag,density,K_test)
  real Tv,Tvf,Hem(row_i,column_i,100),O2(row_i,column_i,100),Kh1,Kh2,Kvf(row_i+1,column_i+1),K11(row_i+1,column_i+1),K2,r1,r2,Tvb,K_test(row_i,column_i)
  integer m,n,Hnum(row_i,column_i),density
  logical initial_state(row_i,column_i),riseflag(row_i,column_i)
 
  
  
  
  if (Tvf>Tv)then       !判断应力是否增长                                       
						Khr_line=0                                                   
							  riseflag(m,n)=.true.                                    !凹点无需特殊处理
				if(Tv>A)then         !原来的应力是否超出PM空间范围
				K_step=0
				Khr_line=1/Kh2
				elseif(Tvf<-A)then  !增长后的应力是否小于PM空间范围
				K_step=0
				Khr_line=1/Kh1
				elseif(Tvf>A.and.Tv<A.and.Tv>=-A)then  !增长后的应力大于最大应力，原应力处于范围中
					   Hnum(m,n)=1
					   Hem(m,n,1)=A
					   O2(m,n,1)=-A
					   K_step=-twice_jifen1(Tv,A,-A,density)
					   Khr_line=1/Kh2
					                                                 !!
			   elseif(Tv<-A.and.Tvf>-A.and.Tvf<=A)then  !原应力小于PM空间最小应力，增长后的应力处于范围中
			          Hnum(m,n)=1
					   Hem(m,n,1)=Tvf
					   O2(m,n,1)=-A
					   K_step=-twice_jifen1(-A,Tvf,-A,density)
					   Khr_line=-(1/Kh1-1/Kh2)*twice_jifen3(Tvf,A,-A,density)-(1/Kh2)*twice_jifen1(-A,A,-A,density)
					   
			  elseif(Tv<-A.and.Tvf>A)then    !原应力小于PM空间最小应力，增长后的应力大于最大范围
			          Hnum(m,n)=1
					   Hem(m,n,1)=A
					   O2(m,n,1)=-A
					   K_step=-twice_jifen1(-A,A,-A,density)
					   Khr_line=1/Kh2
					   
			   elseif(Tv>=-A.and.Tvf<=A)then		 !原应力和增长后的应力都处于PM空间范围中
							 
							                                   								  
							          K_step=0                          !先判断是否断点抹去,连续上升无需判断最后断点,**并且最后超越断点积分路径不变
							          Hnumt=Hnum(m,n)
							          Hem_transition=Hem(m,n,1)
							           Hem(m,n,1)=Tvf
								    
					 if(Hnumt>2)then
								      
								          do while((Tvf>Hem(m,n,2)).and.Hnumt>2)                   
                                                 K_step=K_step-twice_jifen1(Hem_transition,Hem(m,n,2),O2(m,n,1),density)
                                                 O2(m,n,1)=O2(m,n,2)
								                 Hem_transition=Hem(m,n,2)
								              do ii=2,Hnumt-1                                                 !对应积分路径改变，积分路径在数组中从小到大，断点相反,由大到小
								                 O2(m,n,ii)=O2(m,n,ii+1)
								                 Hem(m,n,ii)=Hem(m,n,ii+1)
								              enddo
								                 
									          O2(m,n,Hnumt)=-A                  !废弃值统一为－A
								              Hnumt=Hnumt-1
                                               
                                          enddo
                                      if((Tvf<Hem(m,n,2)).and.Hnumt>2)then
                                        K_step=K_step-twice_jifen1(Hem_transition,Tvf,O2(m,n,1),density)
                                      endif
                                      
                                      
                                    if(Hnumt==2)then
								         if(Tvf>Hem(m,n,2))then
								            K_step=K_step-twice_jifen1(Hem_transition,Hem(m,n,2),O2(m,n,1),density)
								           
								            Hem_transition=Hem(m,n,2)                      !新增加的部分
 								            
								            O2(m,n,1)=O2(m,n,2)
								            O2(m,n,2)=-A
								            Hnumt=Hnumt-1
								          else  
								          K_step=K_step-twice_jifen1(Hem_transition,Tvf,O2(m,n,1),density)
								          endif
								    endif      
								    
								    if(Hnumt==1)then
								         K_step=K_step-twice_jifen1(Hem_transition,Tvf,O2(m,n,1),density)                
                                         
                                         O2(m,n,1)=O2(m,n,2)
								         O2(m,n,2)=-A
								       ! Hnum(m,n)=1 
								    
                                    endif
                                    
                                    
                     elseif(Hnumt==2)then
								         if(Tvf>Hem(m,n,2))then
								            K_step=K_step-twice_jifen1(Hem_transition,Hem(m,n,2),O2(m,n,1),density)
								           
								            Hem_transition=Hem(m,n,2)                 !新增加的部分
								            
								            O2(m,n,1)=O2(m,n,2)
								            O2(m,n,2)=-A
								            Hnumt=Hnumt-1
								          else  
								          K_step=K_step-twice_jifen1(Hem_transition,Tvf,O2(m,n,1),density)
								          endif
								          
								    
								    if(Hnumt==1)then
								         K_step=K_step-twice_jifen1(Hem_transition,Tvf,O2(m,n,1),density)                
                                        
                                         O2(m,n,1)=O2(m,n,2)
								         O2(m,n,2)=-A
								       ! Hnum(m,n)=1 
								    
                                    endif
                                    
                       elseif(Hnumt==1)then
								         K_step=K_step-twice_jifen1(Hem_transition,Tvf,O2(m,n,1),density)                
                                         
                                         O2(m,n,1)=O2(m,n,2)
								         O2(m,n,2)=-A
                                    
                          
                        endif
                        Hnum(m,n)=Hnumt
                  
                                    
                                       
                                     !Hnum(m,n)=Hnumt
                                 !K_test(m,n)=Hem(m,n,2)    
							   
							                                  !断点，积分路径更新完成
                               
							   !do ii=Hnum(m,n)+1,100
							    !  Hem(m,n,ii)=initial
							   !enddo
								  Khr=0
								  do ii=1,Hnum(m,n)-1
								  
                                  x_down=O2(m,n,ii)
								  Khr=Khr-(1/Kh1-1/Kh2)*twice_jifen3(Hem(m,n,ii),Hem(m,n,ii+1),x_down,density)
								  enddo  
								  
								  x_down=O2(m,n,Hnum(m,n))
								  
								  
								  Khr=Khr-(1/Kh1-1/Kh2)*twice_jifen3(Hem(m,n,Hnum(m,n)),A,x_down,density)   !积分路径最后一段缺少对应断点，特殊
								 
								 
                                  
								  x_down=O2(m,n,1)                    !有问题
								  y=Tvf
                                  y_down=-A
								  !Khr=Khr-once_jifen3(x_down,y,density)*r2-(1/Kh2)*twice_jifen1(-A,A,y_down,density)								  
								  Khr_line=Khr-(1/Kh2)*twice_jifen1(-A,A,y_down,density)
								 
			endif
								  Khr=K_step*r2/(Tvf-Tv)-Khr_line
                                  !Khr=K_step*r2/(Tvf-Tv)
                                  K_test(m,n)=Khr
								  Kvf(m,n)=1/(1/(K11(m,n)+K2)-Khr)
							      ! Kvf(m,n)=K_step*r2/(Tvf-Tv)
						          
						  
						  
 else if (Tvf<Tv) then
		Khr_line=0
		if(Tv<-A)then
		K_step=0
		Khr_line=1/Kh1
		elseif(Tvf>A)then
		K_step=0
		Khr_line=1/Kh2				       
		elseif(Tvf<-A.and.Tv>-A.and.Tv<=A)then
						       Hnum(m,n)=2
						      Hem(m,n,1)=-A
						       Hem(m,n,2)=A
						       O2(m,n,1)=-A
						       O2(m,n,2)=-A
						       K_step=-twice_jifen4(A,Tv,-A,density)
						       Khr_line=1/Kh1
						       
						       
		 elseif(Tvf<A.and.Tv>A.and.Tvf>=-A)then
		                       	 Hnum(m,n)=2
						      Hem(m,n,1)=Tvf
						       Hem(m,n,2)=A
						       O2(m,n,1)=Tvf
						       O2(m,n,2)=-A
						       K_step=-twice_jifen4(A,A,Tvf,density)
						       Khr_line=-(1/Kh1)*twice_jifen2(-A,A,A,density)-(1/Kh2-1/Kh1)*twice_jifen2(A,Tvf,A,density)
						       
		elseif(Tvf<-A.and.Tv>A)then
		                      	 Hnum(m,n)=2
						      Hem(m,n,1)=-A
						       Hem(m,n,2)=A
						       O2(m,n,1)=-A
						       O2(m,n,2)=-A
						       K_step=-twice_jifen4(A,A,-A,density)
						       Khr_line=1/Kh1
		                      
		                      
         elseif(Tvf>=-A.and.Tv<=A)then			
						       K_step=0
						       
						       
		   if((Tv==0).and.(Tvb==0).and.(Tvf<Tv).and.initial_state(m,n))then         ! 初状态,有此种情况：永远没有初始下降的初状态，因此initial_state一直为真                     
							      
							      
							      do ii=1,Hnum(m,n)
							      O2(m,n,Hnum(m,n)-ii+2)=O2(m,n,Hnum(m,n)-ii+1)
							      Hem(m,n,Hnum(m,n)-ii+2)=Hem(m,n,Hnum(m,n)-ii+1)
							      enddo
							      O2(m,n,1)=Tvf
							      Hem(m,n,1)=Tvf
							      Hnum(m,n)=Hnum(m,n)+1
								  initial_state(m,n)=.false.
							      
							      
							      K_step=K_step-twice_jifen4(Hem(m,n,2),Tv,Tvf,density)
							   
		else
							     !if((Tvf<Tv).and.(Tvb<=Tv))then         !凸点，增加断点和积分路径,必须考虑Tvb＝Tv，不然会造成下降后只有一断点，数组溢出～～
			  if(riseflag(m,n))then
								     
                                          
								  
								     O2_transition=Tv
								     do ii=1,Hnum(m,n)
							             O2(m,n,Hnum(m,n)-ii+2)=O2(m,n,Hnum(m,n)-ii+1)
							             Hem(m,n,Hnum(m,n)-ii+2)=Hem(m,n,Hnum(m,n)-ii+1)
							         enddo
								         O2(m,n,1)=Tvf
							             Hem(m,n,1)=Tvf
							             Hnum(m,n)=Hnum(m,n)+1
							         
							        
							         
							 if(Hnum(m,n)>2)then
							            do while(Tvf<O2(m,n,2).and.Hnum(m,n)>2)                !因为统一废弃值为-A，可以不考虑初始断点位1的情况：Tvf始终>=-A
							               K_step=K_step-twice_jifen4(Hem(m,n,2),O2_transition,O2(m,n,2),density)
							               O2_transition=O2(m,n,2)
							               do ii=2,Hnum(m,n)-1
							                O2(m,n,ii)=O2(m,n,ii+1)
							                Hem(m,n,ii)=Hem(m,n,ii+1)
							                enddo
							                
							                Hnum(m,n)=Hnum(m,n)-1
							             enddo
							             
							             if(Tvf>O2(m,n,2).and.Hnum(m,n)>2)then
							                K_step=K_step-twice_jifen4(Hem(m,n,2),O2_transition,Tvf,density)
							             endif   
							              
							           if(Hnum(m,n)==2)then
							              if(Tvf<O2(m,n,2))then
							               K_step=K_step-twice_jifen4(Hem(m,n,2),O2_transition,-A,density)
							               
							               O2(m,n,1)=-A
							               Hem(m,n,1)=-A
							               Hnum(m,n)=1
							               
							              else
							              K_step=K_step-twice_jifen4(Hem(m,n,2),O2_transition,Tvf,density)
							              endif
							           endif 
							            
							     elseif(Hnum(m,n)==2)then 
							             if(Tvf<O2(m,n,2))then
							               K_step=K_step-twice_jifen4(Hem(m,n,2),O2_transition,-A,density)
							               
							               O2(m,n,1)=-A
							               Hem(m,n,1)=-A
							               Hnum(m,n)=1
							               
							              else
							              K_step=K_step-twice_jifen4(Hem(m,n,2),O2_transition,Tvf,density)
							            endif
							            
							      endif  
							          
							          
							          
							          
							         
							    ! elseif((Tvf<Tv).and.(Tv<Tvb))then         !一直下降，只改变积分路径1，第一个断点即hem最后一个数据＝积分路径1，也变
			 elseif(.NOT.riseflag(m,n))then
								        O2(m,n,1)=Tvf
							            Hem(m,n,1)=Tvf
							            O2_transition=Tv
							            
							if(Hnum(m,n)>2)then
							            do while(Tvf<O2(m,n,2).and.Hnum(m,n)>2)                !因为统一废弃值为-A，可以不考虑初始断点位1的情况：Tvf始终>=-A
							               
							               K_step=K_step-twice_jifen4(Hem(m,n,2),O2_transition,O2(m,n,2),density)
							               O2_transition=O2(m,n,2)
							               
							               do ii=2,Hnum(m,n)-1
							                O2(m,n,ii)=O2(m,n,ii+1)
							                Hem(m,n,ii)=Hem(m,n,ii+1)
							                enddo
							                Hnum(m,n)=Hnum(m,n)-1
							             enddo
							         
							              if(Tvf>O2(m,n,2).and.Hnum(m,n)>2)then
							                K_step=K_step-twice_jifen4(Hem(m,n,2),O2_transition,Tvf,density)
							              endif
							         
							           if(Hnum(m,n)==2)then
							               if(Tvf<O2(m,n,2))then
							                K_step=K_step-twice_jifen4(Hem(m,n,2),O2_transition,-A,density)
							               O2(m,n,1)=-A
							               Hem(m,n,1)=-A
							               Hnum(m,n)=1
							               else
							               K_step=K_step-twice_jifen4(Hem(m,n,2),O2_transition,Tvf,density)
							               endif
							           endif  
							        
							         
							         
							 elseif(Hnum(m,n)==2)then
							            if(Tvf<O2(m,n,2))then
							                K_step=K_step-twice_jifen4(Hem(m,n,2),O2_transition,-A,density)
							               O2(m,n,1)=-A
							               Hem(m,n,1)=-A
							               O2(m,n,2)=-A
							               Hem(m,n,2)=A
							            
							            else
							            K_step=K_step-twice_jifen4(Hem(m,n,2),O2_transition,Tvf,density)
							            endif
							 
							  endif  
							            
					 endif                                               !断点,积分路径更新完成
                  endif
              
                             
                            
                               
							 !K_test(m,n)=Hem(m,n,2)
							   
							   Khr=0
                               do ii=1,Hnum(m,n)-1
							      xstart=O2(m,n,ii)
								  xend=O2(m,n,ii+1)
                                  y_up=Hem(m,n,ii)
								  Khr=Khr-(1/Kh2-1/Kh1)*twice_jifen2(xstart,xend,y_up,density)     !最后一段不算，断点定义是最后一段永远为空
								enddo  
							      Khr=Khr-(1/Kh2-1/Kh1)*twice_jifen1(-A,O2(m,n,1),-A,density)
								  
								  x=Tvf
								  
								 !! y_up=Hem(m,n,Hnum(m,n)-1)                          !此处引起数组溢出，即hnum有无可能为1：必须考虑特殊凸点，即相等造成的
								  
                                  
                                 ! xstart=-A
								  !xend=A
								 ! y_up=A
								  !Khr=Khr-once_jifen2(x,y_up,density)*r1-(1/Kh1)*twice_jifen2(-A,A,A,density)	
                                  Khr_line=Khr-(1/Kh1)*twice_jifen2(-A,A,A,density)
                                 
        endif   
                                 Khr=K_step*r1/(Tvf-Tv)-Khr_line
                                 ! Khr=K_step*r1/(Tvf-Tv)
                                  K_test(m,n)=Khr
						         
						         
						         riseflag(m,n)=.false.
						         Kvf(m,n)=1/(1/(K11(m,n)+K2)-Khr)
						         !Kvf(m,n)=0
						  end if
	end subroutine
	end module