      SUBROUTINE RMDmatrices
      
! computes matrices for renormalized meshless derivatives
! ( see Lanson,Vila,Gaburov,Hopkins)  

      use wp3d_h
      
      implicit none
      
      include 'mpif.h'

      integer::i,j,k,l,in,itab,iopt,ierr
      double precision::u2,hpi,h2i,h5i,d2,arg, &
      dx,dy,dz
      double precision::r(ndim),t1,dumA(ndim,ndim,nmax)
      double precision::dwij
	  
	  iopt=2
      
      do i=n_lower,n_upper
      
      myA(:,:,i)=0.0d0
      
      hpi=h(i)
      h2i=hpi**2
      h5i=h2i*h2i*hpi
	  
      call LLN(i,iopt,hpi)
     
      do in=1,nbn

          j=nn(in)
          
          if ( i .ne. j ) then
           
          dx=x(i,1)-x(j,1)
          dy=x(i,2)-x(j,2)
          dz=x(i,3)-x(j,3)
          
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
          endif
             
          t1=-dx**2*dwij    
          myA(1,1,i)=myA(1,1,i)+vol(j)*t1
          t1=-dx*dy*dwij
          myA(1,2,i)=myA(1,2,i)+vol(j)*t1
          t1=-dx*dz*dwij
          myA(1,3,i)=myA(1,3,i)+vol(j)*t1
          t1=-dy**2*dwij
          myA(2,2,i)=myA(2,2,i)+vol(j)*t1
          t1=-dx*dz*dwij
          myA(1,3,i)=myA(1,3,i)+vol(j)*t1
          t1=-dy*dz*dwij
          myA(2,3,i)=myA(2,3,i)+vol(j)*t1
          t1=-dz**2*dwij
          myA(3,3,i)=myA(3,3,i)+vol(j)*t1
                
         endif
       
         enddo
           
         myA(2,1,i)=myA(1,2,i)
         myA(3,1,i)=myA(1,3,i)
         myA(3,2,i)=myA(2,3,i)
                 
         r(:)=0.0d0

         call gaussj(myA(:,:,i),ndim,ndim,r,1,ndim)
         
         enddo
         
         call MPI_ALLGATHERV(myA(:,:,n_lower),ilen3(myrank), &
         MPI_DOUBLE_PRECISION,dumA,ilen3,idisp3,MPI_DOUBLE_PRECISION, &
         MPI_COMM_WORLD,ierr)
         
         do i=1,n
         A(i,:,:)=dumA(:,:,i)
         enddo
              
       RETURN
      
       end 

           
      SUBROUTINE gaussj(a,n,np,b,m,mp) 

      INTEGER::m,mp,n,np,NMAX 
      double precision a(np,np),b(np,mp) 
      PARAMETER (NMAX=50) 
      INTEGER::i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX) 
      double precision big,dum,pivinv 

      do j=1,n 
      ipiv(j)=0 
      enddo
      
      do i=1,n
      big=0.0d0
      do j=1,n 
      if( ipiv(j).ne.1) then 
      do k=1,n 
      if (ipiv(k).eq.0) then
      if (abs(a(j,k)).ge.big)then 
      big=abs(a(j,k)) 
      irow=j 
      icol=k 
      endif 
      else if  (ipiv(k).gt.1)  then 
      stop
      endif 
      enddo
      endif 
      enddo
 
      ipiv(icol)=ipiv(icol)+1 
      if (irow.ne.icol) then 
      do l=1,n 
      dum=a(irow,l) 
      a(irow,l)=a(icol,l) 
      a(icol,l)=dum 
      enddo
      do l=1,m 
      dum=b(irow,l) 
      b(irow,l)=b(icol,l) 
      b(icol,l)=dum 
      enddo
      endif 
      indxr(i)=irow
      indxc(i)=icol 
      if (a(icol,icol) .eq. 0.0d0 ) then
      stop
      endif
      pivinv=1.0d0/a(icol,icol) 
      a(icol,icol)=1.0d0
      do l=1,n 
      a(icol,l)=a(icol,l)*pivinv 
      enddo
      do l=1,m 
      b(icol,l)=b(icol,l)*pivinv 
      enddo
      do ll=1,n 
      if(ll.ne.icol)then 
      dum=a(ll,icol) 
      a(ll,icol)=0. 
      do l=1,n 
      a(ll,l)=a(ll,l)-a(icol,l)*dum 
      enddo
      do l=1,m 
      b(ll,l)=b(ll,l)-b(icol,l)*dum 
      enddo
      endif 
      enddo
      enddo
      do l=n,1,-1 
      if(indxr(l).ne.indxc(l))then 
      do k=1,n 
      dum=a(k,indxr(l))
      a(k,indxr(l))=a(k,indxc(l)) 
      a(k,indxc(l))=dum 
      enddo
      endif 
      enddo
     
      return 
      END 
	  
      subroutine coef(i,j,c)
	  
      use wp3d_h
      
      implicit none
	  
      integer::i,j,itab
      double precision::c(ndim),r(ndim),dx,dy,dz,u2, &
      d2,h2i,h5i,dwij
	 
	  dx=x(i,1)-x(j,1)
      dy=x(i,2)-x(j,2)
      dz=x(i,3)-x(j,3)
	  
      h2i=h(i)**2
      h5i=h2i**2*h(i)
          
      if ( iperiodic .eqv. .true. ) then
      call modbound(dx,dy,dz)
      endif
            
      d2=dx**2+dy**2+dz**2
                  
      U2=D2/H2I

      IF ( u2 .le. 4.0d0 ) then
      ITAB=INT(CTAB*U2)+1
      DWIJ=DWTAB(ITAB)/H5I
      else
      DWIJ=0.0d0
      endif
         
      r(1)=dx*dwij
      r(2)=dy*dwij
      r(3)=dz*dwij
       
      c=vol(j)*matmul(A(i,:,:),r(:))
         
      return
	 
      end