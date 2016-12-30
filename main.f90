        MODULE wp3d_h
		     
! module file for the implementation of the Moving Finite Volume method
        
        integer::ntab,ndim,nmax   
        integer::ncmx,ncmy,ncmz,nopt,nv
        integer::nquad,nsubc,null,mxcell,incell, &
        mxnode,ncell,nstep,nhp,inittype,solverchoice, &
        imtype,nghost,root,limtype,divbtype,gradtype,maxnproc
        parameter(ncmx=250,ncmy=250,ncmz=250,nopt=32,nv=8,maxnproc=50)
        parameter(ntab=100000,nmax=500000,nghost=2)  
        parameter(ndim=3,nhp=3,nstep=10)
        parameter(nsubc=2**ndim,nquad=2*ndim-1,null=0, &
        mxnode=2*nmax)
        logical::usquad,periodic
        logical,allocatable::exterior(:)

        double precision,allocatable::x(:,:),v(:,:),h(:), &
        cs(:),vol(:),gradu(:,:,:),lim(:,:),gradw(:,:,:), &
        divb(:),eb(:),gradB(:),Vudot(:,:),psidot(:), &
        gam(:),phi(:),am(:),dzeta(:),A(:,:,:)
        
        double precision,allocatable::myh(:),myvol(:), &
        mydivb(:),myphi(:),mygrad(:,:,:),myVudot(:,:), &
        myam(:),mypsidot(:),myA(:,:,:),mydzeta(:), &
        mylim(:,:),mydti(:),mytau(:),myvsigmax(:)
        
        double precision,allocatable::rcrit2(:),quad(:,:)
        double precision,allocatable::mid(:,:),clsize(:)
        double precision::rsize,theta
        integer,allocatable::subp(:,:)
        
        double precision::norm
                         
        integer, allocatable::nn(:)
        integer::nbn,nd,n
         
        double precision::xmin,ymin,xmax,ymax,zmin,zmax, &
        dx1,dy1,dz1,boxx,boxy,boxz,hh,xcmin,xcmax,ycmin,ycmax, &
        zcmin,zcmax
        integer::nx1,ny1,nz1,ncx,ncy,ncz
        integer,allocatable::HOC(:,:,:),LL(:)
        logical,allocatable::out(:)
                
        double precision::ctab,wtab(ntab),dwtab(ntab), &
        dwdhtab(ntab),eta2,pi,t,tf,mdt,etamax, G, &
        tp,tprin,racc,rhocrit1,rhocrit2,rhocrit3,rhocrit4,ch,cp, &
        hr,c0,s0,MJR,epsa,epsr,Kpoly,trelax,ufloor,dfloor, &
        pfloor,pbet,del,rho1,rtaper,Hz,beta,T0,voltot,totmass

        integer*4::nsink,nacc,nit,ndump,periodictypex, &
        periodictypey,periodictypez,hchoice,setupchoice, &
        outputtype,typeini
        
        integer*4::myrank,n_lower,n_upper
        integer*4::ilen1(0:maxnproc-1),idisp1(0:maxnproc-1)
        integer*4::ilen2(0:maxnproc-1),idisp2(0:maxnproc-1)
        integer*4::ilen3(0:maxnproc-1),idisp3(0:maxnproc-1)
        integer*4::ilen4(0:maxnproc-1),idisp4(0:maxnproc-1)
        integer*4::ilen5(0:maxnproc-1),idisp5(0:maxnproc-1)

        logical::useFPM,iperiodic,useprimitive,  &
        hvar,ilagrangian,hevol,iPowell,iDedner, &
        barotropic,igrav,stretch,relax,restart
        
        character*20::outfn(25),fname
        
        END MODULE wp3d_h
        
	
        program MFV3d
      
  ! 3D Godunov SPH code based on Vila's discretization scheme 
  ! for hyperbolic conservation laws
  ! Hyperbolic/parabolic divergence cleaning scheme(Dedner et al.)
  ! optional coupling to self-gravity of the plasma
  ! applications: 3D protostellar collapse-simulation of accretion disks

        USE wp3d_h
   
        implicit none
        
        include 'mpif.h'
                
        double precision::u(nmax,nv),wprim(nmax,nv),psi(nmax)
        
        integer::nsim,i,p,ierr
        character*1::ndx1
                
 ! initial setting of parameters 

        pi=4.0*datan(1.0d0)
        ndump=0
        norm=1.0d0/pi
        dfloor=1.0d-6
        pfloor=1.0d-6
        etamax=0.01d0
        cp=0.1d0
        theta=0.7d0
                
        useprimitive=.true.
        hvar=.true.
        usquad=.true.
        
        open(unit=1,file='run_pars')
        
        read(1,"(8X,I6)") typeini
! choice of initial condition type: 1: test problems 2:magnetized collapse 3:Evrard collapse
        read(1,"(9X,I6)") inittype
! number of the test within the test suite of problems
        read(1,"(5X,I6)") nsim
! index number of the simulation
        read(1,"(11X,I6)") outputtype
! output format: 1: SPLASH format 2: silo format 3: vtr format 
        read(1,"(8X,I6)") hchoice
! smoothing length calculation 1: fixed number of neighbours 2: iterative procedure
        read(1,"(13X,I6)") solverchoice
! Riemann solver choice 1: HLL solver 2: HLLD solver
        read(1,"(8X,I6)") limtype 
! Choice of the limiter 1: limiter from Gaburov and Nitadori 2: limiter proposed by V. Spingel ( see AREPO paper )
        read(1,"(9X,I6)") divbtype 
! calculation of divb: 1: meshless gradient expression 2: Eqns. (42) from Gaburov and Nitadori
        read(1,"(9X,I6)") gradtype 
! calculation of divb: 1: meshless gradient expression 2: Eqns. (43) from Gaburov and Nitadori
        read(1,"(12X,L1)") iLagrangian
! toggles between Eulerian and Lagrangian mode of the formalism
        read(1,"(8X,L1)") iPowell
! includes Powell correction terms
        read(1,"(8X,L1)") iDedner
! includes Dedner cleaning scheme
        read(1,"(6X,L1)") igrav
! includes self-gravity using the tree formalism
        read(1,"(11X,L1)") barotropic
! Barotropic EOS
        read(1,"(8X,L1)") restart
! Restart option 

        close(1)
		
        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,p,ierr)
        
        call allocate_dynamic_arrays
        
        t=0.0d0

! initialisation of the kernel tables 

        CALL TABULINIT

        tp=0.0d0
        NIT=0
        mdt=1.0d-6
                
! setup of the initial conditions 

        OPEN(1,file="icinfo.dat")       

        select case (typeini)
        case (1)
        call init(u,wprim)
        case (2)
        call init_magcollapse(wprim)
        case (3)
        call evrard_ini(u,wprim)
        end select 

        write(ndx1,'(I1)') nsim 
        
        close(1)
        
        call pararange(1,n,p)
       
        OUTFN(1)=trim(fname) // "stat" // NDX1 // ".dat"
        OUTFN(2)=trim(fname) // "diag" // NDX1 // ".dat"  
        OUTFN(3)=trim(fname) // NDX1
            
        OPEN(2,file=outfn(1))
        OPEN(3,file=outfn(2))
      
        if ( myrank .eq. 0) then
        write(6,*) 'initialising...'
        endif

! First compute neighbour lists and 
! geometric quantities

        CALL build_linked_list 

        call hvol(wprim)
        
        if ( igrav .eqv. .true. ) then
        call MKTREE(wprim)
        endif

        do i=1,n

        call primitive_to_conservative(i,wprim(i,:),u(i,:))

        divb(i)=0.0d0
        if ( iDedner .eqv. .true. ) then
        psi(i)=0.0d0
        endif
        u(i,:)=vol(i)*u(i,:)
        enddo
        
        if ( myrank .eq. 0 ) then
        select case (outputtype)  
        case (1)
!       splash format 
        call PDUMP(wprim)
        case (2)
!       silo format
        call silo_3d(wprim)
        case (3)
!       vtr format 
        call plot(wprim)
        end select 
        endif
        
        call RMDmatrices
        
        CALL CALCDOTS(u,wprim,psi) 
        
        if ( myrank .eq. 0 ) then
        write(6,*) 'starting main loop...'
        endif

        do while (.true.) 

        CALL MAINIT(u,wprim,psi)
        
! Main integration loop 
       
        if ( T .gt. TF ) THEN  
        if ( myrank .eq. 0 ) then    
        write(2,*) 'End of integration has been reached:t=',t
        select case (outputtype)  
        case (1)
!       splash format 
        call PDUMP(wprim)
        case (2)
! silo format
        call silo_3d(wprim)
        case (3)
! vtr format 
        call plot(wprim)
        end select 
        endif
        exit
        endif
                        
        enddo 
 
        if ( myrank .eq. 0 ) then
        select case (outputtype)  
        case (1)
!       splash format 
        call PDUMP(wprim)
        case (2)
! silo format
        call silo_3d(wprim)
        case (3)
! vtr format 
        call plot(wprim)
        end select 
        endif
        
        call MPI_Finalize(ierr)

        close(2)
        close(3)
      
        end program MFV3d
        
		
       subroutine allocate_dynamic_arrays
          
       use wp3d_h
        
       allocate(cs(nmax),h(nmax),vol(nmax),eb(nmax), &
       lim(nmax,nv),divb(nmax),gradB(nmax),&
       gam(nmax),exterior(nmax))
         
       if ( igrav .eqv. .true. ) then
       allocate(x(mxnode,ndim))
       else
       allocate(x(nmax,ndim))
       endif
           
       if ( igrav .eqv. .true. ) then       
! tree code memory allocation
       allocate(rcrit2(1:mxnode),quad(1:mxnode,nquad))
       allocate(mid(1:mxnode,ndim),clsize(1:mxnode))
       allocate(subp(1:mxnode,nsubc),am(1:mxnode), &
       phi(mxnode),dzeta(nmax))
       allocate(LL(1:nmax),HOC(0:ncmx,0:ncmy,0:ncmz))
       else
       allocate(LL(1:nmax),HOC(0:ncmx,0:ncmy,0:ncmz),am(nmax))
       endif
           
       if ( iDedner .eqv. .true. ) then
       allocate(psidot(nmax))
       endif

       allocate(nn(nmax))
       allocate(Vudot(nmax,nv),A(nmax,ndim,ndim))
       if ( useprimitive .eqv. .true. ) then
       allocate(gradw(nmax,nv,ndim))
       else
       allocate(gradu(nmax,nv,ndim))
       endif
       allocate(out(nmax))
                       
       return
                
        end
		
  
        SUBROUTINE MAINIT(u,wprim,psi)

        use wp3d_h

        IMPLICIT NONE
                
        double precision::u(nmax,nv),wprim(nmax,nv),psi(nmax)
                
! ***********************************************************
!     One full iteration of the hydro code
! ***********************************************************   

! Advance particle positions and velocities:
        
        NIT=NIT+1

        CALL ADVANCE(u,wprim,psi)
        
        if ( myrank .eq. 0) then
        CALL ENOUT(u,wprim)
        endif
        
        if ( tp .ge. tprin ) then 
        
        tp=0.0d0
        
! Write results to file 

        if ( myrank .eq. 0 ) then      
        select case (outputtype)  
        case (1)
!       splash format 
        call PDUMP(wprim)
        case (2)
! silo format
        call silo_3d(wprim)
        case (3)
! vtr format 
        call plot(wprim)
        end select 
        endif
        
        ENDIF 

       RETURN
       END
           
      SUBROUTINE TABULINIT
! *****************************************************************
!     Calculate tabulated values of smoothing kernel for SPH summation
! **********************************************************************  

      USE wp3d_h

      integer::i
      double precision::uu2,w,dw,dwdh
                                   
! Compute tabulations of W(u) and dW(u)/du:
      DO I=1,NTAB
      UU2=4.*FLOAT(I-1)/FLOAT(NTAB-1)+1.E-15
      WTAB(I)=W(DSQRT(UU2))
      DWTAB(I)=DW(DSQRT(UU2))   
      DWDHTAB(i)=DWDH(DSQRT(UU2))
      ENDDO

      CTAB=FLOAT(NTAB-1)/4
        
      RETURN
      END

! ***********************************************************************
      double precision FUNCTION W(arg)
!     Calculate smoothing kernel function W 
!     See Rasio and Shapiro, ApJ 401, 226 (1992), Eq. 4

      use wp3d_h

      double precision::arg

      IF (arg.LT.1.) THEN
       W=1.-1.5*arg**2+0.75*arg**3
      ELSE IF (arg.LT.2.) THEN
       W=0.25*(2.-arg)**3
      ELSE
       W=0.
      ENDIF
      W=W*norm

      RETURN
      END

      double precision FUNCTION DW(arg)
      
!     Calculate the derivative of the smoothing kernel function W 

      use wp3d_h

      double precision::arg

      IF (arg.LT.1.) THEN
      DW=-3.*arg+2.25*arg**2
      ELSE IF (arg.LT.2.) THEN
      DW=-0.75*(2.-arg)**2
      ELSE
      DW=0.
      ENDIF
      DW=DW*norm/arg

      RETURN
      END
          
      double precision FUNCTION DWDH(arg)
      
!     Calculate derivative of the smoothing kernel function W with respect to the 
!     smoothing length

      use wp3d_h
          
      double precision::arg

      IF (arg.LT.1.) THEN
       DWDH=-2.0+6*arg**2-3*arg**3
      ELSE IF (arg.LT.2.) THEN
       DWDH=3*(2-arg)**2/4.-(2-arg)**3/2.
      ELSE
       DWDH=0.0d0
      ENDIF
      DWDH=DWDH*norm

      RETURN
      END