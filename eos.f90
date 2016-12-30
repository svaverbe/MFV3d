        double precision function pbar(i,rhoi)

        use wp3d_h

        implicit none
        integer::i
        double precision::rhoi,k1,k2,k3
                
        k1=c0**2/rhocrit1**(2.0d0/5.0d0)
        k2=k1*rhocrit2**(1.0d0/4.0d0)
        k3=k2*rhocrit3*(1.15d0-(5.0d0/3.0d0))
                
! barotropic EOS, includes effects related to the first and second collapse phase
                
        if ( rhoi .le. rhocrit1 ) then
        pbar=c0**2*rhoi
        gam(i)=1.0d0
        endif
        
        if (( rhoi .gt. rhocrit1 ) .and. ( rhoi .le. rhocrit2 )) then
        gam(i)=7.0d0/5.0d0              
        pbar=k1*rhoi**gam(i)
        endif
                
        if (( rhoi .gt. rhocrit2 ) .and. ( rhoi .le. rhocrit3 )) then
        gam(i)=1.15d0           
        pbar=k2*rhoi**gam(i)
        endif
                
        if ( rhoi .gt. rhocrit3 ) then
        gam(i)=5.0d0/3.0d0              
        pbar=k3*rhoi**gam(i)
        endif
        
        return

        end
  