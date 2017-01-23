      subroutine find_limiters_w(wprim)

      USE wp3d_h

      implicit none
      
      include 'mpif.h'

      integer::i,k,l,idx,iopt,ierr
                  
      double precision::wprim(nmax,nv)
      double precision::dx,dy,dz,tiny, &
      kappa(nv),dx2,dy2,dz2,maxwreci(nv),maxwi(nv), &
      minwreci(nv),minwi(nv),wrec(nv),dwij(nv),psiij(nv)
      double precision::dumlim(nv,1:nmax),kappamin,kappamax, &
      normgradw

      kappa(1)=1.0d0
      kappa(2)=1.0d0
      kappa(3)=1.0d0
      kappa(4)=1.0d0
      kappa(5)=1.0d0
      kappa(6)=1.0d0
      kappa(7)=1.0d0
      kappa(8)=1.0d0
      kappamin=1.0d0
      kappamax=2.0d0
      
      iopt=1
 
      tiny=1.0d-10
                                 
      select case (limtype) 
                  
      case (1) 

      do i=n_lower,n_upper
                          
      maxwreci(:)=-1.0d+25
      minwreci(:)=1.0d+25
      maxwi(:)=-1.0d+25
      minwi(:)=1.0d+25
      
      call LLN(i,iopt,h(i))
         
      do k=1,nbn
          
      idx=nn(k)
	  
      dx=x(idx,1)-x(i,1)
      dy=x(idx,2)-x(i,2)
      dz=x(idx,3)-x(i,3)

      if ( iperiodic .eqv. .true. ) then
      call modbound(dx,dy,dz)
      endif

      dx2=dx/2
      dy2=dy/2
      dz2=dz/2
          
      do l=1,nv
          
      wrec(l)=wprim(i,l)+gradw(i,l,1)*dx2+ &
      gradw(i,l,2)*dy2+gradw(i,l,3)*dz2

      maxwreci(l)=dmax1(maxwreci(l),wrec(l))
      minwreci(l)=dmin1(minwreci(l),wrec(l))
          
      maxwi(l)=dmax1(maxwi(l),wprim(idx,l))
      minwi(l)=dmin1(minwi(l),wprim(idx,l))
          
      enddo
                 
      enddo

      do k=1,nv
      
      mylim(k,i)=dmin1(1.0d0,kappa(k)*(maxwi(k)-wprim(i,k))/ &
      (maxwreci(k)-wprim(i,k)+tiny),kappa(k)*(wprim(i,k)-minwi(k))/ &
      (wprim(i,k)-minwreci(k)+tiny))

      enddo

      enddo

      case (2)
         
      do i=n_lower,n_upper
      
      mylim(:,i)=1.0d+25
	  
      call LLN(i,iopt,h(i))
                 
      maxwi(:)=-1.0d+25
      minwi(:)=1.0d+25
        
      do k=1,nbn
          
      idx=nn(k)
         
      do l=1,nv
          
      maxwi(l)=dmax1(maxwi(l),wprim(idx,l))
      minwi(l)=dmin1(minwi(l),wprim(idx,l))
          
      enddo
                  
      enddo
      
                  
      do k=1,nbn
                  
      idx=nn(k)

      if ( i .ne. idx ) then

      dx=x(idx,1)-x(i,1)
      dy=x(idx,2)-x(i,2)
      dz=x(idx,3)-x(i,3)

      if ( iperiodic .eqv. .true. ) then
      call modbound(dx,dy,dz)
      endif

      dx2=dx/2
      dy2=dy/2
      dz2=dz/2

      do l=1,nv
                  
      dwij(l)=gradw(i,l,1)*dx2+gradw(i,l,2)*dy2+ &
      gradw(i,l,3)*dz2
                  
      if ( dwij(l) .eq. 0.0d0 ) then       
      psiij(l)=1.0d0       
      else        
      if ( dwij(l) .gt. 0.0d0 ) then    
      psiij(l)=(maxwi(l)-wprim(i,l))/dwij(l)    
      else    
      psiij(l)=(minwi(l)-wprim(i,l))/dwij(l)   
      endif  
      endif
      
      enddo

      mylim(:,i)=dmin1(mylim(:,i),psiij(:),1.0d0)

      endif
      
      enddo  
      
      enddo
      
      case (3)
      
      do i=n_lower,n_upper
                                
!     maxwreci(:)=-1.0d+25
!     minwreci(:)=1.0d+25
      maxwi(:)=-1.0d+25
      minwi(:)=1.0d+25
      
      call LLN(i,iopt,h(i))
         
      do k=1,nbn
          
      idx=nn(k)
	  
      dx=x(idx,1)-x(i,1)
      dy=x(idx,2)-x(i,2)
      dz=x(idx,3)-x(i,3)

      if ( iperiodic .eqv. .true. ) then
      call modbound(dx,dy,dz)
      endif

      dx2=dx/2
      dy2=dy/2
      dz2=dz/2
          
      do l=1,nv
          
!     wrec(l)=wprim(i,l)+gradw(i,l,1)*dx2+ &
!     gradw(i,l,2)*dy2+gradw(i,l,3)*dz2

!     maxwreci(l)=dmax1(maxwreci(l),wrec(l))
!     minwreci(l)=dmin1(minwreci(l),wrec(l))
          
      maxwi(l)=dmax1(maxwi(l),wprim(idx,l))
      minwi(l)=dmin1(minwi(l),wprim(idx,l))
          
      enddo
                 
      enddo

      do k=1,nv
      
      kappa(k)=dmax1(kappamin,kappamax*dmin1(1.0d0,ncondcrit/myconditionnumber(i)))
      
      normgradw=dsqrt(sum(gradw(i,k,:)**2))
     
      mylim(k,i)=dmin1(1.0d0,kappa(k)*(maxwi(k)-wprim(i,k))/ &
      (normgradw*h(i)+tiny),kappa(k)*(wprim(i,k)-minwi(k))/ &
      (normgradw*h(i)+tiny))
      
      enddo

      enddo

      end select
      
      call MPI_ALLGATHERV(mylim(1,n_lower),ilen4(myrank), &
      MPI_DOUBLE_PRECISION,dumlim,ilen4,idisp4,MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD,ierr)
      
      do i=1,n
      lim(i,:)=dumlim(:,i)
      enddo
      
      return

      end
	  

      subroutine find_limiters_u(u)

      USE wp3d_h

      implicit none
      
      include 'mpif.h'

      integer::i,k,l,idx,iopt,ierr
                  
      double precision::u(nmax,nv)
      double precision::dx,dy,dz,tiny, & 
      maxureci(nv),minureci(nv),maxui(nv), &
      minui(nv),kappa(nv),dx2,dy2,dz2,urec(nv),duij(nv),psiij(nv)
      double precision::dumlim(nv,1:nmax)

      kappa(:)=1.0d0
      tiny=1.0d-10
      
	  iopt=1

      select case (limtype) 

      case (1)

      do i=n_lower,n_upper
                      
      maxureci(:)=-1.0d+25
      minureci(:)=1.0d+25
      maxui(:)=-1.0d+25
      minui(:)=1.0d+25
	  
      call LLN(i,iopt,h(i))
          
      do k=1,nbn
          
      idx=nn(k)
          
      dx=x(idx,1)-x(i,1)
      dy=x(idx,2)-x(i,2)
      dz=x(idx,3)-x(i,3)

      if ( iperiodic .eqv. .true. ) then
      call modbound(dx,dy,dz)
      endif

      dx2=dx/2
      dy2=dy/2
      dz2=dz/2
          
      do l=1,nv
          
      urec(l)=vol(i)*(u(i,l)/vol(i)+gradu(i,l,1)*dx2+ &
      gradu(i,l,2)*dy2+gradu(i,l,3)*dz2)

      maxureci(l)=dmax1(maxureci(l),urec(l))
      minureci(l)=dmin1(minureci(l),urec(l))
          
      maxui(l)=dmax1(maxui(l),u(idx,l)/vol(idx))
      minui(l)=dmin1(minui(l),u(idx,l)/vol(idx))

      enddo
          
      enddo

      do k=1,nv

      mylim(k,i)=dmin1(1.0d0,kappa(k)*(maxui(k)-u(i,k))/ &
      (maxureci(k)-u(i,k)+tiny),kappa(k)*(u(i,k)-minui(k))/ &
      (u(i,k)-minureci(k)+tiny))

      enddo
      enddo
                  
      case (2)
                  
      do i=n_lower,n_upper
      
      mylim(:,i)=1.0d+25
	  
      call LLN(i,iopt,h(i))
                  
      maxui(:)=-1.0d+25
      minui(:)=1.0d+25
          
      do k=1,nbn
          
      idx=nn(k)
          
      do l=1,nv
          
      maxui(l)=dmax1(maxui(l),u(idx,l)/vol(idx))
      minui(l)=dmin1(minui(l),u(idx,l)/vol(idx))
          
      enddo
                  
      enddo
          
      do k=1,nbn

      idx=nn(k)

      if ( i .ne. idx ) then
                  
      dx=x(idx,1)-x(i,1)
      dy=x(idx,2)-x(i,2)
      dz=x(idx,3)-x(i,3)

      if ( iperiodic .eqv. .true. ) then
      call modbound(dx,dy,dz)
      endif

      dx2=dx/2
      dy2=dy/2
      dz2=dz/2

      do l=1,nv
                 
      duij(l)=gradu(i,l,1)*dx2+gradu(i,l,2)*dy2+ &
      gradu(i,l,3)*dz2
                  
      if ( duij(l) .eq. 0.0d0 ) then       
      psiij(l)=1.0d0        
      else       
      if ( duij(l) .gt. 0.0d0 ) then         
      psiij(l)=(maxui(l)-u(i,l)/vol(i))/duij(l)      
      else    
      psiij(l)=(minui(l)-u(i,l)/vol(i))/duij(l)     
      endif
      endif
                  
      enddo

      mylim(:,i)=dmin1(mylim(:,i),psiij(:),1.0d0)
      
      endif
                  
      enddo
      
                     
      enddo
                  
      end select 
      
      call MPI_ALLGATHERV(mylim(1,n_lower),ilen4(myrank), &
      MPI_DOUBLE_PRECISION,dumlim,ilen4,idisp4,MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD,ierr)
      
      do i=1,n
      lim(i,:)=dumlim(:,i)
      enddo

      return

      end