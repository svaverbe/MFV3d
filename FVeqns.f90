      SUBROUTINE CALCDOTS(u,wprim,psi)
          
      use wp3d_h
          
      IMPLICIT NONE
      
      include 'mpif.h'  
      
! *****************************************************
!     This routine calculates d(Vu)/dt in the 
!     discrete conservation law(see Vila,Lanson)
! *****************************************************
          
      double precision::uleft(nv),uright(nv), &
      uleftr(nv),urightr(nv),normnij,F(nv),dx,dy,dz,d, &
      nij(ndim),ivel(ndim),ax,u(nmax,nv),S(nv),Sgrav(nv), &
      gradpsi(n_lower:n_upper,ndim),bxc,psic,wprim(nmax,nv), &
      psi(nmax),t1,u2,psiflux,eps,dgradw(ndim),gaccxi, &
      gaccyi,gacczi,dzetai, &
      h5i,h5,dwij1,dwij2,phii,help(ndim), &
	  cij(ndim),cji(ndim),dumVudot(1:nv,1:nmax), &
      dxinterface,dyinterface,dzinterface,Fpsi

      integer::i,j,k,jn,in,itab,iopt,ierr
       
! for the update of the accelerations, we use the interacting neighbour list

         eps=1.0d-3
		 iopt=2
          
         do i=n_lower,n_upper
         
! for particle i

         myVudot(:,i)=0.0d0
         mydivb(i)=0.0d0
         gradpsi(i,:)=0.0d0
         Fpsi=0.0d0
		 
         if ( igrav .eqv. .true. ) then
	     help(:)=0.0d0
         endif
		 
         call LLN(i,iopt,h(i))
		
         do in=1,nbn

         j=nn(in)
		 
          if ( i .ne. j ) then
		 
          call coef(i,j,cij)
          
          jn=0
           
          do k=1,nbn
           
          if ( nn(k) .eq. i ) then
           
          jn=k
          
          endif
          
          enddo
          
          if ( jn .ne. 0 ) then
          
          call coef(j,i,cji)
           
          dgradw(1)=vol(i)*cij(1)-vol(j)*cji(1)
          dgradw(2)=vol(i)*cij(2)-vol(j)*cji(2)
          dgradw(3)=vol(i)*cij(3)-vol(j)*cji(3)

          normnij=dsqrt(dgradw(1)**2+dgradw(2)**2+dgradw(3)**2)

          dx=x(j,1)-x(i,1)
          dy=x(j,2)-x(i,2)
          dz=x(j,3)-x(i,3)
                  
          if ( iperiodic .eqv. .true. ) then
          call modbound(dx,dy,dz)
          endif
                  
          d=dsqrt(dx**2+dy**2+dz**2)
          nij(1)=dx/d
          nij(2)=dy/d
          nij(3)=dz/d
		 
          call reconstruct(i,j,u,wprim,uleft,uright)
               
! reconstruct left and right states in the Riemann problem defined 
! in between particles i and j

          call rotate_state(nij,uleft,uleftr)
          call rotate_state(nij,uright,urightr)
                  
! rotate the state vector in the coordinate frame with rij=x-axis

          if ( iLagrangian .eqv. .true. ) then
          
          ivel(1)=(wprim(i,2)+wprim(j,2))/2
          ivel(2)=(wprim(i,3)+wprim(j,3))/2
          ivel(3)=(wprim(i,4)+wprim(j,4))/2

! velocity of the interface between particles
 
!         dxinterface=h(i)*dx/(h(i)+h(j))
!         dyinterface=h(i)*dy/(h(i)+h(j))
!         dzinterface=h(i)*dz/(h(i)+h(j))
          
!         ivel(1)=wprim(i,2)+(wprim(j,2)-wprim(i,2))*dxinterface*dx/d**2
!         ivel(2)=wprim(i,3)+(wprim(j,3)-wprim(i,3))*dyinterface*dy/d**2
!         ivel(3)=wprim(i,4)+(wprim(j,4)-wprim(i,4))*dzinterface*dz/d**2
!          
          ax=nij(1)*ivel(1)+nij(2)*ivel(2)+nij(3)*ivel(3)               
          else                  
          ax=0.0d0
          endif
                           
          select case ( solverchoice )

          case (1)

          ! HLL MHD Riemann solver
                  
           call HLLMHD(i,j,uleftr,urightr,nij,F,ax,bxc,psi,psic,psiflux)
                                   
           case (2)
		   
	   ! HLLD MHD Riemann solver

           call HLLDMHD(i,j,uleftr,urightr,nij,F,ax,bxc,psi,psic,psiflux)

           end select 
           
           select case (divbtype) 
                  
           case(1)
                  
           mydivb(i)=mydivb(i)+(wprim(j,6)-wprim(i,6))*cij(1)+ &
           (wprim(j,7)-wprim(i,7))*cij(2)+ &
           (wprim(j,8)-wprim(i,8))*cij(3)
                  
           case (2)
                  
           mydivb(i)=mydivb(i)+bxc*normnij
                  
           end select 
                       
! solve the 1D Riemann problem for the particle pair i,j in 
! the moving reference frame and rotate back to the lab frame

          if ( iDedner ) then
          Fpsi=Fpsi+normnij*psiflux
          endif
          
          ! update the vector of conservative variables 
          
          do k=1,nv
               
          t1=normnij*F(k)
          myVudot(k,i)=myVudot(k,i)+t1
		  
          if ( igrav .and. ( .not. barotropic )) then
          if ( k .eq. 1 ) then
          help(1)=help(1)-dx*t1
          help(2)=help(2)-dy*t1
          help(3)=help(3)-dz*t1
          endif
          endif
		  
          enddo
		  
          if ( igrav .eqv. .true. ) then
		  
	      h5i=h(i)**5
		  
          u2=sum(nij(:)**2)
          
          if ( U2 .ge. 4.0d0 ) THEN
          DWIJ1=0.0d0
          else
          ITAB=INT(CTAB*U2)+1
          DWIJ1=DWTAB(ITAB)/H5I
          ENDIF 
		  
	      H5=h(j)**5
		  
	      u2=(dx**2+dy**2+dz**2)/h(j)**2
          
          if ( U2 .ge. 4.0d0 ) THEN
          DWIJ2=0.0d0
          else
          ITAB=INT(CTAB*U2)+1
          DWIJ2=DWTAB(ITAB)/H5
          endif
		  
          myVudot(2,i)=myVudot(2,i)-G*(dzeta(i)*dwij1+dzeta(j)*dwij2)*dx/2
	      myVudot(3,i)=myVudot(3,i)-G*(dzeta(i)*dwij1+dzeta(j)*dwij2)*dy/2
	      myVudot(4,i)=myVudot(4,i)-G*(dzeta(i)*dwij1+dzeta(j)*dwij2)*dz/2
		  
! terms related to variable smoothing lengths and energy conservation 
! in self-gravitating flows
          
          endif
		  
		  
          if ( iDedner .eqv. .true. ) then 
                  
          select case (gradtype) 
          case (1)
          gradpsi(i,1)=gradpsi(i,1)+(psi(j)/am(j)-psi(i)/am(i))*cij(1)
          gradpsi(i,2)=gradpsi(i,2)+(psi(j)/am(j)-psi(i)/am(i))*cij(2)
          gradpsi(i,3)=gradpsi(i,3)+(psi(j)/am(j)-psi(i)/am(i))*cij(3)
          case (2)                
          gradpsi(i,1)=gradpsi(i,1)+psic*dgradw(1)
          gradpsi(i,2)=gradpsi(i,2)+psic*dgradw(2)
          gradpsi(i,3)=gradpsi(i,3)+psic*dgradw(3)
          end select
                  
          endif

         endif
         
         endif

         enddo
         
         if ( divbtype .eq. 1 ) then
         mydivb(i)=vol(i)*mydivb(i)
         endif
         
         if ( gradtype .eq. 1 ) then
         gradpsi(i,:)=vol(i)*gradpsi(i,:)
         endif
         
         enddo
         
                 
         do i=n_lower,n_upper
                
 ! include Powell's source terms  
 ! and Dedner's correction terms
 
         S(:)=0.0d0

         if ( iDedner .eqv. .true. ) then
         
 !       mypsidot(i)=-Fpsi-wprim(i,1)*myvsigmax(i)**2*mydivb(i)/4
         mypsidot(i)=-Fpsi-wprim(i,1)*ch**2*mydivb(i)
 ! Hyperbolic term
 !       mypsidot(i)=mypsidot(i)-(ch/cp)**2*psi(i)
         mypsidot(i)=mypsidot(i)-psi(i)/mytau(i)
 !  Parabolic term
 
         S(5)=-(wprim(i,6)*gradpsi(i,1)+ &
         wprim(i,7)*gradpsi(i,2)+wprim(i,8)*gradpsi(i,3))
         S(6)=-gradpsi(i,1)
         S(7)=-gradpsi(i,2)
         S(8)=-gradpsi(i,3)

         endif

         if ( iPowell .eqv. .true. ) then
 
         S(2)=S(2)-mydivb(i)*wprim(i,6)
         S(3)=S(3)-mydivb(i)*wprim(i,7)
         S(4)=S(4)-mydivb(i)*wprim(i,8)
         if ( (.not. barotropic) .eqv. .true. ) then
         S(5)=S(5)-mydivb(i)*(wprim(i,2)*wprim(i,6)+ &
         wprim(i,3)*wprim(i,7)+wprim(i,4)*wprim(i,8))
         endif
         S(6)=S(6)-mydivb(i)*wprim(i,2)
         S(7)=S(7)-mydivb(i)*wprim(i,3)
         S(8)=S(8)-mydivb(i)*wprim(i,4)
         endif
		 
         myVudot(1:4,i)=myVudot(1:4,i)-S(1:4)
         if ((.not. barotropic)  .eqv. .true.) then
         myVudot(5,i)=myVudot(5,i)-S(5)
         endif
         myVudot(6:8,i)=myVudot(6:8,i)-S(6:8)

        if ( igrav .eqv. .true. ) then
        
        Sgrav(:)=0.0d0
        
        if ( exterior(i) .eqv. .true. ) then
        gaccxi=0.0d0
        gaccyi=0.0d0
        gacczi=0.0d0
        else
        call TRG(i,gaccxi,gaccyi,gacczi,phii)
        endif
		
	    myphi(i)=phii

        Sgrav(2)=G*am(i)*gaccxi
        Sgrav(3)=G*am(i)*gaccyi
        Sgrav(4)=G*am(i)*gacczi
        if ( (.not. barotropic) .eqv. .true. ) then
        Sgrav(5)=G*(am(i)*(wprim(i,2)*gaccxi+ &
	    wprim(i,3)*gaccyi+wprim(i,4)*gacczi)+ &
	   (gaccxi*help(1)+gaccyi*help(2)+gacczi*help(3))/vol(i))
        endif
        
! source terms originating from self-gravity

        myVudot(2:4,i)=myVudot(2:4,i)-Sgrav(2:4)
        if ((.not. barotropic)  .eqv. .true.) then
        myVudot(5,i)=myVudot(5,i)-Sgrav(5)
        endif

        endif

        enddo
        
        call MPI_ALLGATHERV(myVudot(1,n_lower),ilen4(myrank), &
        MPI_DOUBLE_PRECISION,dumVudot,ilen4,idisp4,MPI_DOUBLE_PRECISION, &
        MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHERV(mydivb(n_lower),ilen1(myrank), & 
        MPI_DOUBLE_PRECISION,divb,ilen1,idisp1,MPI_DOUBLE_PRECISION, & 
        MPI_COMM_WORLD,ierr)
        if ( igrav ) then
        call MPI_ALLGATHERV(myphi(n_lower),ilen1(myrank), & 
        MPI_DOUBLE_PRECISION,phi,ilen1,idisp1,MPI_DOUBLE_PRECISION, & 
        MPI_COMM_WORLD,ierr)
        endif
        if ( iDedner ) then
        call MPI_ALLGATHERV(mypsidot(n_lower),ilen1(myrank), & 
        MPI_DOUBLE_PRECISION,psidot,ilen1,idisp1,MPI_DOUBLE_PRECISION, & 
        MPI_COMM_WORLD,ierr)
        endif
        
        do i=1,n
        Vudot(i,:)=dumVudot(:,i)
        enddo
                                
      RETURN
      END
	  
	  
     
	  
