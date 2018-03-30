module pm_density
       Private A,P_up,P_down
	   
	   real A,P_up,P_down
	   parameter(A=1.5e8,P_down=3e8,P_up=1e7)
	   contains
	   real function density(Po,Pc,sign)
	        real Po,Pc
            integer sign
			if(sign==1)then
			!density=1e-16
            density=5e-13
            elseif(sign==2)then
			!density=1e-15
            density=5e-13
			endif  	    
		end function density
	    
	    
	    
end module