      subroutine reconstruct(i,j,u,wprim,uleft,uright)
          
! this subroutine computes the left and right states at the 
! midpoint between particles i and j
          
      USE wp3d_h
          
      implicit none
                 
      integer::i,j,k
          
      double precision::u(nmax,nv),wprim(nmax,nv), &
      uleft(nv),uright(nv),gradleft(nv),gradright(nv)
      double precision::wleft(nv),wright(nv),wleftlim(nv), &
      wrightlim(nv)
      double precision::dx,dy,dz,dx2,dy2,dz2

      dx=x(j,1)-x(i,1)
      dy=x(j,2)-x(i,2)
      dz=x(j,3)-x(i,3)
                  
      if ( iperiodic .eqv. .true. ) then        
      call modbound(dx,dy,dz)      
      endif
                  
!     dx2=h(i)*dx/(h(i)+h(j))
!     dy2=h(i)*dy/(h(i)+h(j))
!     dz2=h(i)*dz/(h(i)+h(j))

      dx2=dx/2
      dy2=dy/2
      dz2=dz/2
              
      if ( useprimitive .eqv. .true.) then
                  
      do k=1,nv
      
! first apply a pairwise limiter to the gradients
          
      gradleft(k)=gradw(i,k,1)*dx2+ &
      gradw(i,k,2)*dy2+gradw(i,k,3)*dz2
          
      gradright(k)=gradw(j,k,1)*dx2+ &
      gradw(j,k,2)*dy2+gradw(j,k,3)*dz2
                  
!     if ( gradleft(k)*gradright(k) .le. 0.0d0 ) then
 !    wleft(k)=wprim(i,k)
!     wright(k)=wprim(j,k)
!     else
      wleft(k)=wprim(i,k)+lim(i,k)*gradleft(k)
      wright(k)=wprim(j,k)-lim(j,k)*gradright(k)
!     endif

      call pairwiselimiter(i,j,k,wprim,wleft(k),wleftlim(k))
      call pairwiselimiter(i,j,k,wprim,wright(k),wrightlim(k))
             
      enddo
                              
!     call primitive_to_conservative(i,wleft,uleft)
!     call primitive_to_conservative(j,wright,uright)

      call primitive_to_conservative(i,wleftlim,uleft)
      call primitive_to_conservative(j,wrightlim,uright)
                  
      else
                  
      do k=1,nv
                  
      gradleft(k)=gradu(i,k,1)*dx2+ &
      gradu(i,k,2)*dy2+gradu(i,k,3)*dz2
          
      gradright(k)=gradu(j,k,1)*dx2+ &
      gradu(j,k,2)*dy2+gradu(j,k,3)*dz2
          
      if ( gradleft(k)*gradright(k) .le. 0.0d0 ) then
      uleft(k)=u(i,k)
      uright(k)=u(j,k)
      else
      uleft(k)=u(i,k)+lim(i,k)*gradleft(k)
      uright(k)=u(j,k)-lim(j,k)*gradright(k)
      endif
                  
      enddo
                                
      endif
          
      return
         
      end
      
      
      subroutine pairwiselimiter(i,j,k,wprim,w,wlim)
      
      use wp3d_h
      
      integer::i,j,k
      double precision::wprim(nmax,nv),w,wlim,phimin,phimax
      double precision::psi1,psi2,dx,dy,dz,phileft,phiright
      double precision::phiplus,phiminus,delta1,delta2,phiij
      double precision::signarg
      
      psi1=0.5d0
      psi2=0.25d0
    
      phileft=wprim(i,k)
      phiright=wprim(j,k)
      
      phimin=dmin1(phileft,phiright)
      phimax=dmax1(phileft,phiright)
      
      delta1=psi1*dabs(phileft-phiright)
      delta2=psi2*dabs(phileft-phiright)
      
      phiij=phileft+(phiright-phileft)/2
      
!     phiij=phileft+h(i)*(phiright-phileft)/(h(i)+h(j))

      if ( signarg(phimin-delta1) .eq. signarg(phimin) ) then
      phiminus=phimin-delta1
      else
      phiminus=phimin/(1.0d0+delta1/dabs(phimin))
      endif
      
      if ( signarg(phimax+delta1) .eq. signarg(phimax) ) then
      phiplus=phimax+delta1
      else
      phiplus=phimax/(1.0d0+delta1/dabs(phimax))
      endif
      
      if ( phileft .eq. phiright ) then
      wlim=phileft
      else
      if ( phileft .lt. phiright ) then
      wlim=dmax1(phiminus,dmin1(phiij+delta2,w))
      else
      wlim=dmin1(phiplus,dmax1(phiij-delta2,w))
      endif
      endif
      
      return
      end
      
      
      double precision function signarg(arg)
      
      double precision::arg
      
      if ( arg .eq. 0.0d0 ) then
      signarg=1.0d0
      else
      signarg=arg/dabs(arg)
      endif
      
      return
      
      end
	  
	  
      subroutine rotate_state(nij,ui,uir)
          
! computes the vector of conservative variables (rho,rho*vx,rho*vy,e) 
! for particle i in the rotated reference frame with rij=x-axis
          
      USE wp3d_h
          
      implicit none
          
      double precision:: ui(nv),uir(nv), & 
      nij(ndim),vec(ndim),vecr(ndim)
        
      uir(1)=ui(1)
          
      vec(1)=ui(2)
      vec(2)=ui(3)
      vec(3)=ui(4)
          
      call rotate_vector(vec,vecr,nij,1)
          
      uir(2)=vecr(1)
      uir(3)=vecr(2)
      uir(4)=vecr(3)
      uir(5)=ui(5)

      vec(1)=ui(6)
      vec(2)=ui(7)
      vec(3)=ui(8)

      call rotate_vector(vec,vecr,nij,1)
          
      uir(6)=vecr(1)
      uir(7)=vecr(2)
      uir(8)=vecr(3)
      
      return
          
      end
	  
          
      subroutine rotate_vector(vec,vecr,nij,dir)
          
      use wp3d_h
          
      implicit none
          
      double precision::vec(ndim),vecr(ndim),nij(ndim),num,tiny
      integer::dir

      tiny=1.0d-6
          
      select case (dir) 
          
      case (1)
                  
      num=dsqrt(nij(1)**2+nij(2)**2+tiny)

      vecr(1)=vec(1)*nij(1)+vec(2)*nij(2)+vec(3)*nij(3)
      vecr(2)=-nij(2)*vec(1)+nij(1)*vec(2)
      vecr(3)=-nij(1)*nij(3)*vec(1)/num-nij(2)*nij(3)* &
      vec(2)/num+num*vec(3)

      case (-1)

      num=dsqrt(nij(1)**2+nij(2)**2+tiny)

      vecr(1)=nij(1)*vec(1)-nij(2)*vec(2)/num**2- &
      nij(1)*nij(3)*vec(3)/num
      vecr(2)=nij(2)*vec(1)+nij(1)*vec(2)/num**2- &
      nij(2)*nij(3)*vec(3)/num
      vecr(3)=nij(3)*vec(1)+num*vec(3)
          
      end select 
         
      return
          
      end