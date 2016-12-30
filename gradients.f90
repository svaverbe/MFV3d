        subroutine find_gradients_w(wprim)
          
        use wp3d_h
          
        implicit none
        
        include 'mpif.h'
                  
        double precision::wprim(nmax,nv),dx,dy,dz,d,cij(ndim)
        double precision::dumgrad(ndim,nv,nmax)
          
! this subroutine computes the gradients of the primitive 
! variables for each particle 
          
        integer::i,j,k,m,idx,iopt,ierr
		
		iopt=1
                  
        do i=n_lower,n_upper
		
        call LLN(i,iopt,h(i))
                     
        do k=1,nv
          
        mygrad(:,k,i)=0.0d0
        
        do m=1,nbn
          
        idx=nn(m)
                  
        dx=x(i,1)-x(idx,1)
        dy=x(i,2)-x(idx,2)
        dz=x(i,3)-x(idx,3)
                  
        if ( iperiodic .eqv. .true. ) then
        call modbound(dx,dy,dz)
        endif
		
        call coef(i,idx,cij)
                  
        mygrad(:,k,i)=mygrad(:,k,i)+(wprim(idx,k)-wprim(i,k))*cij(:)
         
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
        enddo
        
        return

        end

  
        subroutine find_gradients_u(u)
          
        USE wp3d_h
          
        implicit none
        
        include 'mpif.h'

        integer::i,j,m,k,idx,iopt,ierr
        double precision::u(nmax,nv),dx,dy,dz,d,cij(ndim)
        double precision::dumgrad(ndim,nv,nmax)
       
		iopt=1
                  
        do i=n_lower,n_upper
                
        do k=1,nv
          
        mygrad(:,k,i)=0.0d0
          
        do m=1,nbn
          
        idx=nn(m)
		
        call coef(i,idx,cij)
                  
        dx=x(i,1)-x(idx,1)
        dy=x(i,2)-x(idx,2)
        dz=x(i,3)-x(idx,3)
                  
        if ( iperiodic .eqv. .true. ) then
        call modbound(dx,dy,dz)
        endif
                  
        mygrad(:,k,i)=mygrad(:,k,i)+(u(idx,k)/vol(idx)- &
        u(i,k)/vol(idx))*cij(:)
         
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