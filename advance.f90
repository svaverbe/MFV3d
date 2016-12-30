      SUBROUTINE ADVANCE(u,wprim,psi)
          
      use wp3d_h
          
      IMPLICIT NONE
          
! **************************************************************
!    Advance conservative variables by one time step 
! **************************************************************
                                     
       integer::i,in,j,k
       double precision::dt1,dth1, &
       periodic1,periodic2
       double precision::u(nmax,nv),u1(nmax,nv), &
       u2(nmax,nv),wprim(nmax,nv),wprim1(nmax,nv), &
       wprim2(nmax,nv),psi(nmax),psi1(nmax),psi2(nmax)
       external periodic1,periodic2
         
! Drift-Kick-Drift scheme for advancing the particles

        call TSTEP(wprim,psi,dt1)

        dth1=dt1/2
         
        t=t+dt1
                 
        if ( iLagrangian .eqv. .true. ) then
                 
        do i=1,n
         
! DRIFT step for the particle positions 

        x(i,1)=x(i,1)+dth1*wprim(i,2)
        x(i,2)=x(i,2)+dth1*wprim(i,3)
        x(i,3)=x(i,3)+dth1*wprim(i,4)

        enddo
        
        if ( iperiodic .eqv. .true. ) then
        
        do i=1,n 
        
        select case ( periodictypex ) 
        case (1)
        x(i,1)=periodic1(x(i,1),boxx)
        case (2) 
        x(i,1)=periodic2(x(i,1),boxx)
        end select
                
        select case ( periodictypey ) 
        case (1)
        x(i,2)=periodic1(x(i,2),boxy)
        case (2) 
        x(i,2)=periodic2(x(i,2),boxy)
        end select
                
        select case ( periodictypez ) 
        case (1)
        x(i,3)=periodic1(x(i,3),boxz)
        case (2) 
        x(i,3)=periodic2(x(i,3),boxz)
        end select
         
        enddo

        endif 
               
        CALL build_linked_list

! reconstructs the linked list

        CALL hvol(wprim)
		
! geometric quantities

        if ( igrav .eqv. .true. ) then
        call MKTREE(wprim)
        endif

        endif

        if ( iLagrangian .eqv. .true. ) then

        CALL RMDmatrices

        else

        if ((  ilagrangian .eqv. .false. ) .and. &
        (nit .eq. 1 )) then

        CALL RMDmatrices

        endif

        endif

! renormalized kernel derivatives

        if ( useprimitive .eqv. .true. ) then
        CALL find_gradients_w(wprim)
        else
        CALL find_gradients_u(u)
        endif
        
! evaluation of function gradients

        if ( useprimitive .eqv. .true. ) then
        CALL find_limiters_w(wprim)
        else
        CALL find_limiters_u(u)
        endif

        CALL CALCDOTS(u,wprim,psi) 

! KICK step for the conservative variables

        do i=1,n
                
        if ( barotropic .eqv. .true. ) then
        do k=1,4
        u1(i,k)=u(i,k)-dt1*Vudot(i,k)
        enddo
        do k=6,nv
        u1(i,k)=u(i,k)-dt1*Vudot(i,k)
        enddo
        else
        do k=1,nv
        u1(i,k)=u(i,k)-dt1*Vudot(i,k)
        enddo
        endif
                
        if (iDedner .eqv. .true. ) then
        psi1(i)=psi(i)+dt1*psidot(i)
        endif
                  
        enddo 
        
        do i=1,n
              
        call conservative_to_primitive(i,u1(i,:), &
        wprim1(i,:))
        if ( barotropic .eqv. .true. ) then
        cs(i)=dsqrt(wprim1(i,5)/wprim1(i,1))
        else
        cs(i)=dsqrt(gam(i)*wprim1(i,5)/wprim1(i,1))
        endif
                   
        enddo
                     
! renormalized kernel derivatives

        if ( useprimitive .eqv. .true. ) then
        CALL find_gradients_w(wprim1)
        else
        CALL find_gradients_u(u1)
        endif
        
! evaluation of function gradients

        if ( useprimitive .eqv. .true. ) then
        CALL find_limiters_w(wprim1)
        else
        CALL find_limiters_u(u1)
        endif

        if ( igrav .eqv. .true. ) then
        CALL MKTREE(wprim1)
        else
        do i=1,n
        am(i)=vol(i)*wprim1(i,1) 
        enddo
        endif

        CALL CALCDOTS(u1,wprim1,psi1) 

! KICK step for the conservative variables using a Runge-Kutta 

        do i=1,n
                
        if ( barotropic .eqv. .true. ) then
        do k=1,4
        u2(i,k)=0.5d0*(u1(i,k)+u(i,k)-dt1*Vudot(i,k))
        enddo
        do k=6,nv
        u2(i,k)=0.5d0*(u1(i,k)+u(i,k)-dt1*Vudot(i,k))
        enddo
        else
        do k=1,nv
        u2(i,k)=0.5d0*(u1(i,k)+u(i,k)-dt1*Vudot(i,k))
        enddo
        endif
                
        if ( iDedner .eqv. .true. ) then
        psi2(i)=0.5*(psi(i)+psi1(i)+dt1*psidot(i))
        endif  
                   
        enddo 
                
        do i=1,n
        
        call conservative_to_primitive(i,u2(i,:),wprim2(i,:))
        
        if ( barotropic .eqv. .true. ) then
        cs(i)=dsqrt(wprim2(i,5)/wprim2(i,1))
        else
        cs(i)=dsqrt(gam(i)*wprim2(i,5)/wprim2(i,1))
        endif

        enddo
                
        u(:,:)=u2(:,:)
        wprim(:,:)=wprim2(:,:)
        if ( iDedner ) then
        psi(:)=psi2(:)
        endif

        tp=tp+dt1

        if ( ilagrangian .eqv. .true. ) then
        
        do i=1,n
        
! DRIFT step for the particles

        x(i,1)=x(i,1)+dth1*wprim(i,2)
        x(i,2)=x(i,2)+dth1*wprim(i,3)
        x(i,3)=x(i,3)+dth1*wprim(i,4)
        
        enddo
                
        if ( iperiodic .eqv. .true. ) then
        
        do i=1,n 
        
        select case ( periodictypex ) 
        case (1)
        x(i,1)=periodic1(x(i,1),boxx)
        case (2) 
        x(i,1)=periodic2(x(i,1),boxx)
        end select
                
        select case ( periodictypey ) 
        case (1)
        x(i,2)=periodic1(x(i,2),boxy)
        case (2) 
        x(i,2)=periodic2(x(i,2),boxy)
        end select
                
        select case ( periodictypez ) 
        case (1)
        x(i,3)=periodic1(x(i,3),boxz)
        case (2) 
        x(i,3)=periodic2(x(i,3),boxz)
        end select
        
        enddo

        endif 

       CALL build_linked_list

! reconstructs the link list
       CALL hvol(wprim)
       if ( igrav .eqv. .true. ) then
       call MKTREE(wprim)
       endif
          
      endif
             
      if ( myrank .eq. 0 ) then
      write(2,*) 'step= ',NIT,' t=',t,' dt=',dt1         
      write(6,*) 'step= ',NIT,' t=',t,' dt=',dt1 
      endif
            
      RETURN
      END
	  
	  
      subroutine primitive_to_conservative(i,wi,ui)
                  
! converts primitive to conservative variables: (rho,vx,vy,vz,p,bx,by,bz) to 
! (rho,rho*vx,rho*vy,rho*vz,e,bx,by,bz)
                  
      use wp3d_h
                  
      double precision:: wi(nv),ui(nv)
      integer::i
                  
      ui(1)=wi(1)
      ui(2)=wi(1)*wi(2)
      ui(3)=wi(1)*wi(3)
      ui(4)=wi(1)*wi(4)
                                             
      if ( (.not. barotropic) .eqv. .true. ) then
      ui(5)=wi(5)/(gam(i)-1.0d0)+0.5d0*wi(1)*( &
      wi(2)**2+wi(3)**2+wi(4)**2)+0.5d0*(wi(6)**2+ &
      wi(7)**2+wi(8)**2)                         
      endif          
      ui(6)=wi(6)
      ui(7)=wi(7)
      ui(8)=wi(8)          
          
      return
                  
      end

                  
      subroutine conservative_to_primitive(i,ui,wi)
                  
! converts conservative to primitive variables:(rho,rho*vx,rho*vy,e,bx,by,rho*psi) to 
! (rho,vx,vy,p,bx,by)
                  
      use wp3d_h
                                  
      integer::i
      double precision:: wi(nv),ui(nv),pbar

      wi(1)=ui(1)/vol(i)          
      wi(2)=ui(2)/ui(1)
      wi(3)=ui(3)/ui(1)
      wi(4)=ui(4)/ui(1)
                 
      if ((.not. barotropic)  .eqv. .true.) then
      wi(5)=(gam(i)-1.0d0)*(ui(5)/vol(i)-0.5d0*(ui(2)**2+ &
      ui(3)**2+ui(4)**2)/ui(1)/vol(i)-0.5d0*(ui(6)**2+ui(7)**2+ &
      ui(8)**2)/vol(i)**2) 
      else
      wi(5)=pbar(i,wi(1))
      endif
                                  
      wi(6)=ui(6)/vol(i)
      wi(7)=ui(7)/vol(i)
      wi(8)=ui(8)/vol(i)
                   
      return
                  
      end