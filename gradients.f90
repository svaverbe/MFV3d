        subroutine find_gradients_w(wprim)
          
        use wp3d_h
          
        implicit none
        
        include 'mpif.h'
                  
        double precision::wprim(nmax,nv),cij(ndim)
        double precision::dumgrad(ndim,nv,nmax),dwij,u2,d2, &
        dx,dy,dz,arg,h2i,h5i
          
! this subroutine computes the gradients of the primitive 
! variables for each particle 
          
        integer::i,j,k,m,idx,iopt,ierr,itab
		
		iopt=1
                  
        do i=n_lower,n_upper
		
        call LLN(i,iopt,h(i))
                     
        do k=1,nv
          
        mygrad(:,k,i)=0.0d0
        
        do m=1,nbn
          
        idx=nn(m)
        
        if ( myconditionnumber(i) .lt. 10.0d0*ncondcrit ) then
        
! renormalized meshless derivative
                
        call coef(i,idx,cij)
                  
        mygrad(:,k,i)=mygrad(:,k,i)+(wprim(idx,k)-wprim(i,k))*cij(:)
        
        else
        
! if the condition number if bad, we go for an SPH like gradient estimate

        h2i=h(i)**2
        h5i=h2i**2*h(i)

        dx=x(idx,1)-x(i,1)
        dy=x(idx,2)-x(i,2)
        dz=x(idx,3)-x(i,3)
          
        if ( iperiodic .eqv. .true. ) then
        call modbound(dx,dy,dz)
        endif
            
        d2=dx**2+dy**2+dz**2
                  
        U2=D2/H2I

        IF ( U2 .le. 4.0d0 ) THEN
        arg=CTAB*u2
        ITAB=dint(arg)+1
        DWIJ=DWTAB(ITAB)/H5I
        else
        DWIJ=0.0d0
        endiF
        
        mygrad(1,k,i)=mygrad(1,k,i)+(wprim(idx,k)-wprim(i,k))*vol(idx)*dx*dwij
        mygrad(2,k,i)=mygrad(2,k,i)+(wprim(idx,k)-wprim(i,k))*vol(idx)*dy*dwij
        mygrad(3,k,i)=mygrad(3,k,i)+(wprim(idx,k)-wprim(i,k))*vol(idx)*dz*dwij
        
        endif
         
        enddo
                  
        enddo
               
        enddo
        
        call MPI_ALLGATHERV(mygrad(1,1,n_lower),ilen5(myrank), &
        MPI_DOUBLE_PRECISION,dumgrad,ilen5,idisp5,MPI_DOUBLE_PRECISION, &
        MPI_COMM_WORLD,ierr)
        
        do i=1,n
        do j=1,nv
        do k=1,ndim
        gradw(i,j,k)=dumgrad(k,j,i)
        enddo
        enddo
        
!        do k=1,nv
!        write(60,*) x(i,:),k,gradw(i,k,:)
!        enddo
        
        enddo
        
        return

        end

  
        subroutine find_gradients_u(u)
          
        USE wp3d_h
          
        implicit none
        
        include 'mpif.h'

        integer::i,j,m,k,idx,iopt,ierr
        double precision::u(nmax,nv),cij(ndim)
        double precision::dumgrad(ndim,nv,nmax)
       
		iopt=1
                  
        do i=n_lower,n_upper
                
        do k=1,nv
          
        mygrad(:,k,i)=0.0d0
          
        do m=1,nbn
          
        idx=nn(m)
		
        call coef(i,idx,cij)
                       
        mygrad(:,k,i)=mygrad(:,k,i)+(u(idx,k)/vol(idx)- &
        u(i,k)/vol(i))*cij(:)
         
        enddo
          
        enddo
                   
        enddo
        
        call MPI_ALLGATHERV(mygrad(1,1,n_lower),ilen5(myrank), &
        MPI_DOUBLE_PRECISION,dumgrad,ilen5,idisp5,MPI_DOUBLE_PRECISION, &
        MPI_COMM_WORLD,ierr)
        
        do i=1,n
        do j=1,nv
        do k=1,ndim
        gradu(i,j,k)=dumgrad(k,j,i)
        enddo
        enddo
        enddo
         
        return
          
        end