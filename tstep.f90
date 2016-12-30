      SUBROUTINE TSTEP(wprim,psi,dt)
          
 ! calculates the system time step
          
      use wp3d_h
          
      IMPLICIT NONE
      
      include 'mpif.h'
   
      double precision::dt,dti(nmax),tiny
      PARAMETER (TINY=1.E-10)

      integer::i,j,idx,iopt,ierr
      double precision::cmax,cfl,Li,dx,dy,dz,valfveni,vmhdi, &
      bij1,bij2,eps2,wprim(nmax,nv),dti1,atot,dtkin,mych,vfpsi, &
      vmhdj,valfvenj,vsig,psi(nmax)

      eps2=1.0d-10  
      cfl=0.2d0
      dt=1000.0d0
      iopt=1
          
      if ( iDedner .eqv. .true. ) then
      mych=0.0d0
      endif

      do i=n_lower,n_upper
            
      valfveni=dsqrt((wprim(i,6)**2+wprim(i,7)**2+ &
      wprim(i,8)**2+tiny)/wprim(i,1))
    
      call LLN(i,iopt,h(i))
      
      myvsigmax(i)=0.0d0
      
      do j=1,nbn

      idx=nn(j)
      
      valfvenj=dsqrt((wprim(idx,6)**2+wprim(idx,7)**2+ &
      wprim(idx,8)**2+tiny)/wprim(idx,1))

      dx=x(idx,1)-x(i,1)
      dy=x(idx,2)-x(i,2)
      dz=x(idx,3)-x(i,3)
                  
      if ( iperiodic .eqv. .true. ) then         
      call modbound(dx,dy,dz)
      endif

      if ( i .ne. idx ) then      
      bij1=(wprim(i,6)*dx+wprim(i,7)*dy+wprim(i,8)*dz) &
      /dsqrt(dx**2+dy**2+dz**2)
      else
      bij1=0.0d0
      endif
      
      if ( i .ne. idx ) then 
      bij2=-(wprim(idx,6)*dx+wprim(idx,7)*dy+wprim(idx,8)*dz) &
      /dsqrt(dx**2+dy**2+dz**2)
      else
      bij2=0.0d0
      endif
             
      vmhdi=dsqrt((cs(i)**2+valfveni**2)+ &
      dsqrt((cs(i)**2+valfveni**2)**2- &
      4*cs(i)**2*bij1**2/wprim(i,1)))/dsqrt(2.0d0) 

      vmhdj=dsqrt((cs(idx)**2+valfvenj**2)+ &
      dsqrt((cs(idx)**2+valfvenj**2)**2- &
      4*cs(idx)**2*bij2**2/wprim(idx,1)))/dsqrt(2.0d0)     
      
      vsig=vmhdi+vmhdj-dmin1(0.0d0,(-dx*(wprim(i,2)-wprim(idx,2))-dy*( &
      wprim(i,3)-wprim(idx,3))-dz*(wprim(i,4)-wprim(idx,4)))/dsqrt(dx**2+dy**2+dz**2))
            
      myvsigmax(i)=dmax1(myvsigmax(i),vsig)
   
      enddo
      
! timescale for damping the scalar potential in the Dedner scheme

      vfpsi=dsqrt(cs(i)**2+valfveni**2+(2*psi(i)/am(i)/myvsigmax(i))**2)
      
      cmax=dmax1(myvsigmax(i)/2,vfpsi)

      if ( iDedner .eqv. .true. ) then
      mych=dmax1(mych,cmax)
      endif

      Li=(3*vol(i)/pi/4)**(1.0d0/3.0d0)
      
      mytau(i)=Li/cmax/cp
                          
      mydti(i)=2*Li/myvsigmax(i)
     
      if ( igrav .eqv. .true. ) then

      atot=dsqrt(sum(Vudot(i,2:4)**2))/vol(i)
      dtkin=dsqrt(Li/atot)
      mydti(i)=dmin1(mydti(i),dtkin)

      endif

      enddo
      
      
      call MPI_ALLGATHERV(mydti(n_lower),ilen1(myrank), & 
      MPI_DOUBLE_PRECISION,dti,ilen1,idisp1,MPI_DOUBLE_PRECISION, & 
      MPI_COMM_WORLD,ierr)
      if ( iDedner .eqv. .true. ) then
      CALL MPI_ALLREDUCE(mych,ch,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
      MPI_COMM_WORLD,ierr)  
      endif
      
    
      do i=1,n

      Li=(3*vol(i)/pi/4)**(1.0d0/3.0d0)

      if ( iDedner .eqv. .true. ) then
      dti1=dmin1(dti(i),Li/ch)
      else
      dti1=dti(i)
      endif
  
      if ( dti1 .lt. dt ) then
      dt=dti1
      endif
                
!    system time step

      enddo

      dt=cfl*dt
 
      if ( dt .lt. mdt ) then
      
      if ( myrank .eq. 0 ) then

      write(2,*) 'Integration has been stopped because of a too small timestep'
      
      select case (outputtype)  
      case (1)
!     splash format 
      call PDUMP(wprim)
      case (2)
!     silo format
      call silo_3d(wprim)
      case (3)
!     vtr format 
      call plot(wprim)
      end select 
      
      endif

      stop
      endif 
    
      RETURN 
    
      END