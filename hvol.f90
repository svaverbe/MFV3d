        SUBROUTINE hvol(wprim)
      
        use wp3d_h
                
        implicit none
        
        include 'mpif.h'
               
 ! determines smoothing length and the neighbour list of of the particles, including all
 ! *interactive* pairs, and also computes the volume and mass of the particles
            
        integer::i,j,k,iopt,itab,ierr
        double precision::hpi,voli,weight,u2,dx,dy,dz, &
	    DFG,t1,t2,h1,h2,omegai,dr,wprim(nmax,nv),ni
	
        if ( hvar .eqv. .true. ) then
                
        select case ( hchoice )
                
        case (1)
                
! determine the smoothing length so that the number of neighbours is constant

        do i=n_lower,n_upper
                
        hpi=h(i)
        
        call find_smoothing_length(i,hpi)

        myh(i)=hpi

        myvol(i)=weight(i,myh(i))
       
        myam(i)=myvol(i)*wprim(i,1)
 
        enddo
                     
        case (2)
                
! determine consistent smoothing lenghts and particle weights 
! by means of a Newton-Raphson method 

        do i=n_lower,n_upper
               
        hpi=h(i)

        call findhvol(i,voli,hpi)
           
        myh(i)=hpi
        myvol(i)=voli
        
        myam(i)=myvol(i)*wprim(i,1)
       
        enddo
                
        end select 
                
        else

        do i=n_lower,n_upper
                                   
        if (( iLagrangian .eqv. .false. ) .and. &
        ( nit .eq. 1)) then

        call findhvol(i,voli,h(i))

        myvol(i)=voli
        myh(i)=h(i)
        myam(i)=myvol(i)*wprim(i,1)
        
        endif
                                       
        enddo
                
        endif 
                
! neighbour lists
! if hvar=.false. than h(i)=const 

        call MPI_ALLGATHERV(myh(n_lower),ilen1(myrank), & 
        MPI_DOUBLE_PRECISION,h,ilen1,idisp1,MPI_DOUBLE_PRECISION, & 
        MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHERV(myvol(n_lower),ilen1(myrank), & 
        MPI_DOUBLE_PRECISION,vol,ilen1,idisp1,MPI_DOUBLE_PRECISION, & 
        MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHERV(myam(n_lower),ilen1(myrank), & 
        MPI_DOUBLE_PRECISION,am,ilen1,idisp1,MPI_DOUBLE_PRECISION, & 
        MPI_COMM_WORLD,ierr)
        
        if ( igrav .eqv. .true. ) then
		
		iopt=1
	   
        do i=n_lower,n_upper
                
        hpi=h(i)
	    voli=vol(i)
		
        call LLN(i,iopt,hpi)
	   
! compute dzeta quantities to deal with spatially varying gravitational softening lengths

        t1=0.0d0
        t2=0.0d0
		
        do k=1,nbn
		
        j=nn(k) 

        if ( i .ne. j ) then 

        dx=x(i,1)-x(j,1)
        dy=x(i,2)-x(j,2)
        dz=x(i,3)-x(j,3)
		
        if ( iperiodic .eqv. .true. ) then       
        call modbound(dx,dy,dz)
        endif
		
        dr=dsqrt(dx**2+dy**2+dz**2)
		
        h1=dr/hpi
        h2=float(ndim)/hpi
		
        u2=h1**2
		
        if ( u2 .le. 4.0d0 ) then   
        itab=int(ctab*u2)+1
        t1=t1+h1*dwtab(itab)+h2*wtab(itab)
        t2=t2+am(j)*DFG(dr,hpi)
        
        endif
        
        endif

        enddo
        
        select case (hchoice)
        case (1)
        ni=float(nopt)*3.0d0/(4*pi)/hpi**3
        case (2)
        ni=1.0d0/voli
        end select 
        
        select case (hchoice) 
        case (1)
!       omegai=1.0d0-(1.0d0/float(ndim)/ni)*dfloat(nopt)*9.0d0/(4*pi)/hpi**3
        mydzeta(i)=0.0d0
        omegai=0.0d0
        case (2)
        omegai=1.0d0-(hpi/float(ndim)/ni)*t1
        mydzeta(i)=am(i)*hpi*t2/omegai/ni/float(ndim)
        end select 
       
       enddo
       
       call MPI_ALLGATHERV(mydzeta(n_lower),ilen1(myrank), & 
       MPI_DOUBLE_PRECISION,dzeta,ilen1,idisp1,MPI_DOUBLE_PRECISION, & 
       MPI_COMM_WORLD,ierr)
     
       endif
                
       return
        
       end
        
        
        double precision function weight(i,hpi)
        
        use wp3d_h
        
        implicit none 
        
        integer::i,j,k,itab
        double precision::hpi,u2,dx,dy,dz,ni
        
        ni=0.0d0
  
        do k=1,nbn
        
        j=nn(k)  

        dx=x(i,1)-x(j,1)
        dy=x(i,2)-x(j,2)
        dz=x(i,3)-x(j,3)
               
        if ( iperiodic ) then 
        call modbound(dx,dy,dz)  
        endif
                
        u2=(dx**2+dy**2+dz**2)/hpi**2
                
        if ( u2 .le. 4.0d0 ) then  
        itab=int(ctab*u2)+1    
        ni=ni+wtab(itab)/hpi**3   
        endif

        enddo
        
        weight=1.0d0/ni
        
        return

        end
        

        SUBROUTINE FINDHVOL(i,voli,hpi)
        
        use wp3d_h
     
        IMPLICIT NONE

!   this subroutine determines the smoothing length and
!   volume of each particle using an iterative procedure
                                            
        integer::i,j,it,in,iopt,itab,maxitn,maxit
          
        double precision::d2,dwdhj,hnew,hpi,h2i,h4i
        double precision::dnhi,wi,htol,u2,voli,dx,dy,dz
        double precision::tiny,hminbisec,hmaxbisec, &
        dfdh1,hpiold,func,weight,whi,hini
        logical:: bisection,converged
        parameter(tiny=1.0d-10)
                
        maxitn=25
        maxit=500
        htol=1.0d-2
        iopt=1
        hini=0.1d0
        bisection=.false. 
        hminbisec=0.0d0
        hmaxbisec=1d+6
        
        hpi=h(i)       
        hpiold=hpi
        
        it=0
        converged=.false. 
        
        do while (( it .le. maxit ) .and. ( .not. converged))
        
        h2i=hpi**2
        h4i=h2i*h2i
              
        call LLN(i,iopt,hpi)
       
! Find neighbours of this particle 

        wi=weight(i,hpi)
        
! Compute density value 

        dnhi=0.0d0

        do in=1,nbn

        j=nn(in)
        
        dx=x(i,1)-x(j,1)
        dy=x(i,2)-x(j,2)
        dz=x(i,3)-x(j,3)
        
        if ( iperiodic ) then
        call modbound(dx,dy,dz)
        endif
        
        d2=dx**2+dy**2+dz**2
        
        U2=D2/H2I

        IF ( U2 .ge. 4.0d0 ) THEN
        dwdhj=0.0d0
        ELSE
        ITAB=INT(CTAB*U2)+1
        dwdhj=DWDHTAB(ITAB)/H4I
        ENDIF 
        
        dnhi=dnhi+dwdhj
                
        enddo
        
        whi=32*pi*hpi**3/float(nopt)/3
        
! estimate of the particle weight based on the current smoothing length

        func=whi-wi
                
! this function needs to be zero 

       IF (.not. bisection) THEN
!
!--Newton-Raphson iteration
!
       dfdh1 = 32*pi*h2i/float(nopt)+dnhi*wi**2
       hnew = hpi - func/dfdh1
!   write(2,*) 'newton raphson (',it,'): hnew = ',hnew
       ELSE
!
!--Bisection iteration
!
       IF ( func .lt. 0.0d0 ) THEN
       hmaxbisec = hpi
       ELSE
       hminbisec = hpi
       ENDIF
       hnew = 0.5*(hminbisec + hmaxbisec)
!      write(6,*) hnew
       endif
       
!--Don't allow sudden jumps to huge numbers of neighbours
!  (Newton-Raphson only)
!
        IF (.not.bisection) THEN
        IF (hnew.GT.1.2*hpi) THEN
        hnew = 1.2*hpi
        ELSEIF (hnew.LT.0.8*hpi) THEN
        hnew = 0.8*hpi
        ENDIF
        ENDIF
        
        IF ( hnew.LE.0..OR.  & 
              it.EQ.maxitn .OR. nbn.LE.0  &
              .AND. .not. bisection ) THEN
!
!--switch to bisection if not converging or running into trouble
!
               WRITE(2,*) 'WARNING: switching to bisection on'// &
                 ' particle ','(',i,')',' hpi = ',hpi, &
                 ' hnew = ',hnew
               IF (nbn.LE.0) THEN
                  WRITE(2,*) '(particle has no neighbours)'
               ENDIF
               IF (it.EQ.maxit) THEN
                  WRITE(2,*) '(more than ',maxit,' iterations)'
               ENDIF
              
               bisection= .true.
               hminbisec = 0.0d0
               hmaxbisec = 1.0d6
               !--don't have to start with ridiculous h, 
               !  just something reasonably big between above limits
               hnew = 2.*hpiold
!               hnew = 0.5*(hminbisec + hmaxbisec)
!
!--otherwise check for convergence
!
               ELSEIF (DABS(hnew-hpi)/hpi .LT. htol) THEN
               converged=.true. 
               ENDIF        
            
         hpi=hnew

        it=it+1

        enddo

        if ( .not. converged ) then
        
        WRITE (2,*) 'ERROR: iteration failed'
        WRITE (2,*) i,nbn,hnew,hpi,hpiold,wi,dnhi
        stop

        else

        voli=wi

        endif
        
      RETURN
      
      END  
	  
      SUBROUTINE find_smoothing_length(i,hpi)
         
      use wp3d_h
        
      implicit none
         
!    Find the smoothing length of a particle so that its number of neighbours is exactly equal to nopt
                                             
      integer::i,j,k,ind(nmax),iopt
          
      double precision:: R2(nmax),h2,dx,dy,dz,hpi

      iopt=1
        
      CALL LLN(i,iopt,hpi)
		
      IF ( nbn .ge. nopt ) THEN
        
!  sort neighbors in increasing order of distance to particle i 
!  and find the nopt closest neighbors
!  the number of neighbors is kept exactly at nopt 

      DO j=1,nbn
          
      k=nn(j)
          
      dx=x(i,1)-x(k,1)
      dy=x(i,2)-x(k,2)
      dz=x(i,3)-x(k,3)
          
      if ( iperiodic .eqv. .true. ) then
      call modbound(dx,dy,dz)
      endif
      R2(J)=dx**2+dy**2+dz**2
      ENDDO 

      CALL indexx(nbn,r2,ind)        
      hpi=0.5d0*dsqrt(R2(IND(NOPT)))
      ELSE

      DO WHILE ( nbn .lt. NOPT ) 

      H2=2.0d0*hpi
         
      CALL LLN(i,iopt,H2) 
          
      IF ( nbn .ge. NOPT ) THEN

!  sort neighbors in increasing order of distance to particle i 
!  and find the nopt closest ones

      DO j=1,nbn
          
      k=nn(j)

      dx=x(i,1)-x(k,1)
      dy=x(i,2)-x(k,2)
      dz=x(i,3)-x(k,3)
          
      if ( iperiodic .eqv. .true. ) then
          
      call modbound(dx,dy,dz)
          
      endif

      R2(J)=dx**2+dy**2+dz**2
          
      ENDDO 

      CALL indexx(nbn,r2,ind)

      hpi=0.5*dsqrt(R2(IND(NOPT)))

      ELSE

      hpi=H2

      ENDIF 
          
      ENDDO
                  
      ENDIF
         
      RETURN
      END 
	  

      SUBROUTINE indexx(np,arr,indx)
          
      use wp3d_h
          
      implicit none
         
      integer:: np,indx(nmax),M,NSTACK
      double precision:: arr(nmax)
      PARAMETER (M=7,NSTACK=50)
      
! Indexes an array arr(1:np), i.e., outputs the array indx(1:np) such that arr(indx(j)) 
! is in ascending order for j =1, 2, ::. ;NP . The input quantities np and arr are not changed. 

      integer::i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      double precision::ah
      
      do j=1,np
      indx(j)=j
      enddo
      jstack=0
      l=1
      ir=np
10    if (ir-l.lt.M) then 
      do j=l+1,ir
        indxt=indx(j)
        ah=arr(indxt)
      do i=j-1,l,-1
        if(arr(indx(i)).le.ah) goto 20 
        indx(i+1)=indx(i)
      enddo 
      i=l-1
20    indx(i+1)=indxt 
      enddo 
      if (jstack.eq.0) return 
        ir=istack(jstack) 
        l=istack(jstack-1) 
        jstack=jstack-2 
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
      if(arr(indx(l)).gt.arr(indx(ir)))then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
      endif 
      if(arr(indx(l+1)).gt.arr(indx(ir)))then 
        itemp=indx(l+1) 
        indx(l+1)=indx(ir) 
        indx(ir)=itemp 
      endif 
      if(arr(indx(l)).gt.arr(indx(l+1)))then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
      endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        ah=arr(indxt)
30     continue 
        i=i+1 
      if(arr(indx(i)).lt.ah) goto 30 
40    continue 
        j=j-1
      if(arr(indx(j)).gt.ah) goto 40
      if(j.lt.i) goto 50 
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 30
50      indx(l+1)=indx(j) 
        indx(j)=indxt 
        jstack=jstack+2 
        if(jstack.gt.NSTACK) then
        stop
        endif 
! 'NSTACK too small in indexx’ 
      if(ir-i+1.ge.j-l)then 
        istack(jstack)=ir 
        istack(jstack-1)=i
        ir=j-1
      else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
      endif
      endif
      goto 10
      END
      