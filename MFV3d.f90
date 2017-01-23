         MODULE wp3d_h
		     
! module file for the implementation of the Moving Finite Volume method
        
        integer::ntab,ndim,nmax    
        integer::ncmx,ncmy,ncmz,maxnne,nopt,nv
        integer::nquad,nsubc,null,mxcell,incell, &
        mxnode,ncell,nstep,nhp,inittype,solverchoice, &
        imtype,nghost,root,limtype,divbtype,gradtype
        parameter(ncmx=100,ncmy=100,ncmz=100,nopt=32,nv=8)
        parameter(ntab=100000,nmax=100000,nghost=2)  
        parameter(ndim=3,maxnne=100,nhp=3,nstep=10)
        parameter(nsubc=2**ndim,nquad=2*ndim-1,null=0, &
        mxnode=2*nmax)
        logical::usquad,periodic
        logical,allocatable::exterior(:)

        double precision,allocatable:: x(:,:),v(:,:),h(:), &
        cs(:),vol(:),gradu(:,:,:),lim(:,:),gradw(:,:,:), &
        divb(:),eb(:),gradB(:),Vudot(:,:),psidot(:), &
        gam(:),phi(:),am(:),dzeta(:)
		
        double precision,allocatable:: rcrit2(:),quad(:,:)
        double precision,allocatable:: mid(:,:),clsize(:)
        double precision::rsize,theta
        integer,allocatable::subp(:,:)
        
        double precision::norm
                          
        double precision,allocatable::C(:,:,:)
                
        integer, allocatable::nni(:,:),nne(:),nn(:)
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
        pfloor,pbet,del,rho1,rtaper,Hz,beta,T0,voltot,totmass, &
        tdump

        integer::nsink,nacc,nit,ndump,periodictypex, &
        periodictypey,periodictypez,hchoice,setupchoice, &
        outputtype,typeini
        
        logical::useFPM,iperiodic,useprimitive,  &
        imhd,hvar,ilagrangian,hevol,iPowell,iDedner, &
        barotropic,igrav,stretch,relax,restart,fexists
        
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
                
        double precision::u(nmax,nv),wprim(nmax,nv),psi(nmax),dt
        
        integer::nsim,i
        character*20::restartfile
        character*1::ndx1
                
 ! initial setting of parameters 

        pi=4.0*datan(1.0d0)
        ndump=0
        norm=1.0d0/pi
        nsim=1
        inittype=1
        divbtype=2
        gradtype=2
        dfloor=1.0d-6
        pfloor=1.0d-6
        etamax=0.01d0
        outputtype=3
        typeini=1
                
        useprimitive=.true.
        imhd=.true.
        hvar=.true.
        hchoice=1
        solverchoice=2
        limtype=1
        iLagrangian=.true.
        iPowell=.true.
        iDedner=.true.     
        igrav=.false. 
        barotropic=.false. 
        usquad=.true.
		stretch=.false.
		relax=.false. 
        restart=.false. 
		
        theta=0.7d0

        call allocate_dynamic_arrays

! initialisation of the kernel tables 

        CALL TABULINIT

        tp=0.0d0
        NIT=0
        mdt=1.0d-6
                
! setup of the initial conditions 

        open(unit=1,file="icinfo.dat")

        select case (typeini)
        case (1)
        call init(u,wprim)
        case (2)
        call init_magcollapse(wprim)
        case (3)
        call evrard_ini(u,wprim)
        end select 

        write(ndx1,'(I1)') nsim 
        
        if ( restart ) then
        
        restartfile=trim(fname) // trim(".dump")
        
        open(unit=20,file=restartfile)
        
        read(20,"(ES25.15)") tdump
       
        do i=1,n
       
        read(20,"(12ES25.15)") x(i,1),x(i,2),x(i,3),wprim(i,1), &
        wprim(i,2),wprim(i,3),wprim(i,4),wprim(i,5),wprim(i,6), &
        wprim(i,7),wprim(i,8),psi(i)
        
        enddo
        
        t=tdump
        
        close(20)
        
        else
        
        t=0.0d0
        
        endif
        
        close(1)
            
        OUTFN(1)=trim(fname) // "stat" // NDX1 // ".dat"
        OUTFN(2)=trim(fname) // "ei" // NDX1 // ".dat"  
        OUTFN(3)=trim(fname) // "ek" // NDX1 // ".dat"
        if ( imhd .eqv. .true. ) then
        OUTFN(4)=trim(fname) // "eb" // NDX1 // ".dat"
        endif
        OUTFN(5)=trim(fname) // "et" // NDX1 // ".dat"
        OUTFN(6)=trim(fname) // "jt" // NDX1 // ".dat"
        OUTFN(7)=trim(fname) // "pt" // NDX1 // ".dat"
        OUTFN(8)=trim(fname) // "plots" // NDX1 // ".dat"
        OUTFN(9)=trim(fname) // "rm" // NDX1 // ".dat"
        OUTFN(10)=trim(fname) // NDX1  
        OUTFN(11)=trim(fname) // NDX1
        OUTFN(12)=trim(fname) // "restart" // NDX1
        OUTFN(14)=trim(fname) // "mt" // NDX1 // ".dat"  
        if ( imhd ) then
        OUTFN(13)=trim(fname) // "divb" // NDX1 // ".dat"
        OUTFN(15)=trim(fname) // "ep" // NDX1 // ".dat"
        endif
        
        inquire(file=outfn(1),EXIST=fexists)

        if ( fexists .eqv. .true. ) then
        
        if ( restart .eqv. .true. ) then
        open(unit=2,file=outfn(1),position='append')
        else
        open(unit=2,file=outfn(1),status='replace')
        endif
 
        else

        open(unit=2,file=outfn(1),status='new')

        endif
        
        
        inquire(file=outfn(2),EXIST=fexists)
        
        if ( fexists .eqv. .true. ) then
        
        if ( restart .eqv. .true. ) then
        open(unit=3,file=outfn(2),position='append')
        else
        open(unit=3,file=outfn(2),status='replace')
        endif
 
        else

        open(unit=3,file=outfn(2),status='new')

        endif
        
        
        inquire(file=outfn(3),EXIST=fexists)
        
        if ( fexists .eqv. .true. ) then
        
        if ( restart .eqv. .true. ) then
        open(unit=4,file=outfn(3),position='append')
        else
        open(unit=4,file=outfn(3),status='replace')
        endif
 
        else

        open(unit=4,file=outfn(3),status='new')

        endif
        
      
        if ( imhd .eqv. .true. ) then
        
        inquire(file=outfn(4),EXIST=fexists)
        
        if ( fexists .eqv. .true. ) then
        
        if ( restart .eqv. .true. ) then
        open(unit=7,file=outfn(4),position='append')
        else
        open(unit=7,file=outfn(4),status='replace')
        endif
 
        else

        open(unit=7,file=outfn(4),status='new')

        endif
        endif
        
        
        inquire(file=outfn(5),EXIST=fexists)
        
        if ( fexists .eqv. .true. ) then
        
        if ( restart .eqv. .true. ) then
        open(unit=8,file=outfn(5),position='append')
        else
        open(unit=8,file=outfn(5),status='replace')
        endif
 
        else

        open(unit=8,file=outfn(5),status='new')

        endif
        
        
        inquire(file=outfn(6),EXIST=fexists)
        
        if ( fexists .eqv. .true. ) then
        
        if ( restart .eqv. .true. ) then
        open(unit=9,file=outfn(6),position='append')
        else
        open(unit=9,file=outfn(6),status='replace')
        endif
 
        else

        open(unit=9,file=outfn(6),status='new')

        endif
        
        
        inquire(file=outfn(7),EXIST=fexists)
        
        if ( fexists .eqv. .true. ) then
        
        if ( restart .eqv. .true. ) then
        open(unit=10,file=outfn(7),position='append')
        else
        open(unit=10,file=outfn(7),status='replace')
        endif
 
        else

        open(unit=10,file=outfn(7),status='new')

        endif
        
        
        inquire(file=outfn(8),EXIST=fexists)
        
        if ( fexists .eqv. .true. ) then
        
        if ( restart .eqv. .true. ) then
        open(unit=11,file=outfn(8),position='append')
        else
        open(unit=11,file=outfn(8),status='replace')
        endif
 
        else

        open(unit=11,file=outfn(8),status='new')

        endif
        
        
        inquire(file=outfn(9),EXIST=fexists)
        
        if ( fexists .eqv. .true. ) then
        
        if ( restart .eqv. .true. ) then
        open(unit=13,file=outfn(9),position='append')
        else
        open(unit=13,file=outfn(9),status='replace')
        endif
 
        else

        open(unit=13,file=outfn(9),status='new')

        endif
        
       
        if ( imhd .eqv. .true. ) then
        
        inquire(file=outfn(15),EXIST=fexists)
        
        if ( fexists .eqv. .true. ) then
        
        if ( restart .eqv. .true. ) then
        open(unit=14,file=outfn(15),position='append')
        else
        open(unit=14,file=outfn(15),status='replace')
        endif
 
        else

        open(unit=14,file=outfn(15),status='new')

        endif
        endif
        
        
        inquire(file=outfn(14),EXIST=fexists)
        
        if ( fexists .eqv. .true. ) then
        
        if ( restart .eqv. .true. ) then
        open(unit=18,file=outfn(14),position='append')
        else
        open(unit=18,file=outfn(14),status='replace')
        endif
 
        else

        open(unit=18,file=outfn(14),status='new')

        endif
        
        
        if ( imhd .eqv. .true. ) then
        
        inquire(file=outfn(13),EXIST=fexists)
        
        if ( fexists .eqv. .true. ) then
        
        if ( restart .eqv. .true. ) then
        open(unit=15,file=outfn(13),position='append')
        else
        open(unit=15,file=outfn(13),status='replace')
        endif
 
        else

        open(unit=15,file=outfn(13),status='new')

        endif
        endif
        
       
        write(6,*) 'initialising...'

! First compute neighbour lists and 
! geometric quantities

        CALL build_linked_list  

        call neighbourlists(wprim)
        
        if ( igrav .eqv. .true. ) then
        call MKTREE(wprim)
        endif
        
        do i=1,n

        call primitive_to_conservative(i,wprim(i,:),u(i,:))

        divb(i)=0.0d0
        if ( .not. restart ) then
        if ( iDedner .eqv. .true. ) then
        psi(i)=0.0d0
        endif
        endif
        u(i,:)=vol(i)*u(i,:)
        enddo

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

        if ( igrav .eqv. .true. ) then
        call RMD_coefficients
        CALL CALCDOTS(u,wprim,psi) 
        endif

        write(6,*) 'starting main loop...'
        
        
 !       do while (.true.) 

        CALL MAINIT(u,wprim,psi,dt)
        
! Main integration loop 

        if ( dt .lt. mdt ) then
        
        write(2,fmt="(A,ES25.15)") 'Warning:time step too small,dt=',dt
        
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
        
        CALL restartdump(wprim,psi)
        
!       exit
        
        endif

        if ( t .gt. tf ) THEN   

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
        
        call restartdump(wprim,psi)
        
 !       exit
        
        endif
                        
!        enddo 

        close(2)
        close(3)
        close(4)
        if ( imhd ) then
        close(7)
        endif
        close(8)
        close(9)
        close(10)
        close(11)
        close(12)
        if ( imhd ) then
        close(15)
        endif
        if ( igrav .eqv. .true. ) then
        close(14)
        endif
        close(16)
        close(17)
        close(18)
        close(20)

        end program MFV3d
        
        
       subroutine restartdump(wprim,psi)
       
       use wp3d_h
       
       integer::i
       double precision::wprim(1:nmax,1:nv),psi(1:nmax)
       character*20::restartfile
       
       restartfile=trim(fname) // trim(".dump")
        
       open(unit=20,file=restartfile)
        
       write(20,fmt="(ES25.15)") t
       
       do i=1,n
       
       write(20,fmt="(12ES25.15)") x(i,1),x(i,2),x(i,3),wprim(i,1), &
       wprim(i,2),wprim(i,3),wprim(i,4),wprim(i,5),wprim(i,6), &
       wprim(i,7),wprim(i,8),psi(i)
        
       enddo
      
       return
       
       end
       
        
       subroutine allocate_dynamic_arrays
          
       use wp3d_h
        
       allocate(cs(nmax),h(nmax),vol(nmax), eb(nmax), &
       lim(nmax,nv),divb(nmax),gradB(nmax), &
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
       allocate(LL(1:nmax),HOC(0:ncmx,0:ncmy,0:ncmz),am(1:nmax))
       endif
           
       if ( iDedner ) then
       allocate(psidot(nmax))
       endif

       allocate(nne(nmax),nni(nmax,maxnne),nn(nmax))
       allocate(C(nmax,maxnne,ndim),Vudot(nmax,nv))
       if ( useprimitive ) then
       allocate(gradw(nmax,nv,ndim))
       else
       allocate(gradu(nmax,nv,ndim))
       endif
       allocate(out(nmax))
                       
       return
                
        end
  
        SUBROUTINE MAINIT(u,wprim,psi,dt)

        use wp3d_h

        IMPLICIT NONE
                
        double precision::u(nmax,nv),wprim(nmax,nv),psi(nmax),dt
                
! ***********************************************************
!     One full iteration of the hydro code
! ***********************************************************   

! Advance particle positions and velocities:
        
        NIT=NIT+1

        CALL ADVANCE(u,wprim,psi,dt)

        CALL ENOUT(u,wprim)
        
        if ( tp .ge. tprin ) then 
        
        tp=0.0d0
        
! Write results to file 
       
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
      double precision FUNCTION W(a)
!     Calculate smoothing kernel function W 
!     See Rasio and Shapiro, ApJ 401, 226 (1992), Eq. 4

      use wp3d_h

      double precision::a

      IF (a.LT.1.) THEN
       W=1.-1.5*a**2+0.75*a**3
      ELSE IF (a.LT.2.) THEN
       W=0.25*(2.-a)**3
      ELSE
       W=0.
      ENDIF
      W=W*norm

      RETURN
      END

      double precision FUNCTION DW(a)
      
!     Calculate the derivative of the smoothing kernel function W 

      use wp3d_h

      double precision::a

      IF (a.LT.1.) THEN
      DW=-3.*a+2.25*a**2
      ELSE IF (a.LT.2.) THEN
      DW=-0.75*(2.-a)**2
      ELSE
      DW=0.
      ENDIF
      DW=DW*norm/a

      RETURN
      END
          
      double precision FUNCTION DWDH(a)
      
!     Calculate derivative of the smoothing kernel function W with respect to the 
!     smoothing length

      use wp3d_h
          
      double precision::a 

      IF (a.LT.1.) THEN
       DWDH=-2.0+6*a**2-3*a**3
      ELSE IF (a.LT.2.) THEN
       DWDH=3*(2-a)**2/4.-(2-a)**3/2.
      ELSE
       DWDH=0.0d0
      ENDIF
      DWDH=DWDH*norm

      RETURN
      END

          
      SUBROUTINE init(u,wprim)
     
         use wp3d_h

         implicit none
         
         integer:: idx,nt,idum1,idum2, &
         idum3
         double precision::d0,b0,d,r0,p0,camb, &
         amp,wlength,wnumber,width,ran1,hini
         double precision::u(nmax,nv),wprim(nmax,nv)
                 
         select case ( inittype )
        
         case (1)
        
! strong blast wave test

        fname="blast3d"

        iperiodic=.true. 
        periodictypex=2
        periodictypey=2
        periodictypez=2
        setupchoice=1

        xcmin=0.0d0
        xcmax=1.0d0
        ycmin=0.0d0
        ycmax=1.0d0
        zcmin=0.00
        zcmax=1.0d0
                
        nx1=40
        ny1=40
        nz1=40
        d0=1.0d0
        nt=nx1*ny1*nz1
        n=nt
        tf=0.1d0
        tprin=0.01d0   
           
        r0=0.1d0
        
        gam(:)=5.0/3.0d0

        hini=0.025d0
        
        boxx=xcmax-xcmin
        boxy=ycmax-ycmin
        boxz=zcmax-zcmin
                
        call create_particle_distribution
        
        do idx=1,n
                
        wprim(idx,1)=d0
        wprim(idx,2)=0.0d0
        wprim(idx,3)=0.0d0
        wprim(idx,4)=0.0d0
         
        h(idx)=hini

        d=dsqrt((x(idx,1)-boxx/2.)**2+(x(idx,2)-boxy/2.)**2+ &
        (x(idx,3)-boxz/2)**2)
        if ( d .le. r0 ) then
        p0=10.0d0
        else
        p0=0.1d0
        endif
        wprim(idx,5)=p0
        if ( imhd ) then
        wprim(idx,6)=1.0d0/dsqrt(2.0d0)
        wprim(idx,7)=1.0d0/dsqrt(2.0d0)
        wprim(idx,8)=1.0d0/dsqrt(2.0d0)
        endif
        cs(idx)=dsqrt(gam(idx)*p0/d0)
                
        call primitive_to_conservative(idx,wprim(idx,:),u(idx,:))
        
        enddo
             
        case (2)
        
! Kelvin Helmholtz instability test

        tf=2.0d0
        tprin=0.05d0
        fname="kh3d1"

        iperiodic=.true. 
        periodictypex=1
        periodictypey=1
        periodictypez=1
        setupchoice=1

        xmin=-0.5d0
        xmax=0.5d0
        ymin=-0.5d0
        ymax=0.5d0
        zmin=-0.1d0
        zmax=0.1d0
        idum1=100
        idum2=-500
        idum3=1000
        nx1=50
        ny1=50
        nz1=10
                
        d0=1.0d0
        p0=2.5d0
        nt=nx1*ny1*nz1
        n=nt
        
        hini=0.01d0
                
! total number of particles
        
        gam(:)=7.0/5.0d0
        
        boxx=xmax-xmin
        boxy=ymax-ymin
        boxz=zmax-zmin
        
        call create_particle_distribution
                        
        do idx=1,n  
                
        if ( dabs(x(idx,2)) .gt. 0.25d0 ) then
        wprim(idx,1)=d0
        wprim(idx,2)=-0.5d0
        wprim(idx,3)=0.0d0
        wprim(idx,4)=0.0d0
        else
        wprim(idx,1)=2.0*d0
        wprim(idx,2)=0.5d0
        wprim(idx,3)=0.0d0
        wprim(idx,4)=0.0d0
        endif
        wprim(idx,2)=wprim(idx,2)+0.01*ran1(idum1)
        wprim(idx,3)=wprim(idx,3)+0.01*ran1(idum2)
        wprim(idx,4)=wprim(idx,4)+0.01*ran1(idum3)

        if ( imhd ) then
        wprim(idx,6)=0.0d0
        wprim(idx,7)=0.0d0
        wprim(idx,8)=0.0d0
        endif
        
        h(idx)=hini
        wprim(idx,5)=p0
        cs(idx)=dsqrt(gam(idx)*p0/wprim(idx,1))

        call primitive_to_conservative(idx,wprim(idx,:),u(idx,:))
        
        enddo
        
        case(3)

        iperiodic=.true. 
        periodictypex=1
        periodictypey=1
        periodictypez=1
        setupchoice=1

        tf=2.0d0
        tprin=0.05d0
        fname="kh3d2"
        
        xmin=-0.5d0
        xmax=0.5d0
        ymin=-0.5d0
        ymax=0.5d0
        zmin=-0.1d0
        zmax=0.1d0
        amp=0.1d0
        wlength=1.0d0/2.0d0
        wnumber=2*pi/wlength
        width=0.05d0/dsqrt(2.0d0)
        nx1=50
        ny1=50
        nz1=10
        d0=1.0d0
        p0=2.5d0
        nt=nx1*ny1*nz1
        n=nt
        
        hini=0.01d0
                 
! total number of particles
        
        gam(:)=5.0/3.0d0
        
        boxx=xmax-xmin
        boxy=ymax-ymin
        boxz=zmax-zmin
        
        call create_particle_distribution
                
        do idx=1,n
                        
        if ( dabs(x(idx,2)) .gt. 0.25d0 ) then
        wprim(idx,1)=d0
        wprim(idx,2)=-0.5d0
        wprim(idx,3)=0.0d0
        wprim(idx,4)=0.0d0
        else
        wprim(idx,1)=2.0d0*d0
        wprim(idx,2)=0.5d0
        wprim(idx,3)=0.0d0
        wprim(idx,4)=0.0d0
        endif

        wprim(idx,3)=wprim(idx,3)+amp*dsin(wnumber*x(idx,1))* &
       ( dexp(-(x(idx,2)+0.25d0)**2/(2*width**2))+ &
         dexp(-(x(idx,2)-0.25d0)**2/(2*width**2)))
     
        h(idx)=hini
        wprim(idx,5)=p0
         if ( imhd ) then
         wprim(idx,6)=0.0d0
         wprim(idx,7)=0.0d0
         wprim(idx,8)=0.0d0
         endif
        cs(idx)=dsqrt(gam(idx)*p0/wprim(idx,1))
        
        call primitive_to_conservative(idx,wprim(idx,:),u(idx,:))
        
        enddo
        
        case (4)
        
! blob test

        tf=2.0d0
        tprin=0.05d0
        fname="blob"
                
        iperiodic=.true. 
        periodictypex=2
        periodictypey=2
        periodictypez=2
        setupchoice=1

        xmin=0.0d0
        xmax=3.0d0
        ymin=0.0d0
        ymax=1.0d0
        zmin=0.0d0
        zmax=1.0d0
        nx1=150
        ny1=50
        nz1=50
        d0=1.0d0
        p0=1.0d0
        r0=0.1d0
        gam(:)=7.0d0/5.0d0

        camb=dsqrt(gam(1)*p0/d0)
                
        nt=nx1*ny1*nz1
        n=nt
                
        hini=0.05d0
                
! total number of particles
        
        boxx=xmax-xmin
        boxy=ymax-ymin
        boxz=zmax-zmin
                
        call create_particle_distribution
       
        do idx=1,n
                
        d=dsqrt((x(idx,1)-boxx/2.)**2+(x(idx,2)-boxy/2.)**2+ &
        (x(idx,3)-boxz/2)**2)
        if ( d .le. r0 ) then
        d0=10.0d0
        wprim(idx,2)=0.0d0
        else
        d0=1.0d0
        wprim(idx,2)=2.7d0*camb
        endif
                
        wprim(idx,1)=d0
        wprim(idx,3)=0.0d0
        wprim(idx,4)=0.0d0

         if ( imhd ) then
        wprim(idx,6)=0.0d0
        wprim(idx,7)=0.0d0
        wprim(idx,8)=0.0d0
         endif
        
        h(idx)=hini
        wprim(idx,5)=p0
        cs(idx)=dsqrt(gam(idx)*p0/wprim(idx,1))

        call primitive_to_conservative(idx,wprim(idx,:),u(idx,:))
        
        enddo
        
        case (5)
        
! advection of a current loop 

        fname="loop"

        iperiodic=.true. 
        periodictypex=2
        periodictypey=2
        periodictypez=2
        setupchoice=1
        
        xmin=0.0d0
        xmax=2.0d0
        ymin=0.0d0
        ymax=1.0d0
        zmin=0.0d0
        zmax=0.5d0
        nx1=80
        ny1=40
        nz1=20
        nt=nx1*ny1*nz1
        n=nt
        p0=1.0d0
        b0=0.001d0
        tf=10.0d0
        tprin=0.1d0   
           
        r0=0.3d0
        
        gam(:)=5.0/3.0d0

        hini=0.05d0
        
        boxx=xmax-xmin
        boxy=ymax-ymin
        boxz=zmax-zmin
                
        call create_particle_distribution
        
        do idx=1,n
           
        wprim(idx,2)=2.0d0
        wprim(idx,3)=1.0d0
        wprim(idx,4)=0.5d0
         
        h(idx)=hini

        d=dsqrt((x(idx,1)-boxx/2.)**2+(x(idx,2)-boxy/2.)**2)
          if ( d .le. r0 ) then
        d0=2.0d0
        wprim(idx,6)=(b0/d)*x(idx,2)
        wprim(idx,7)=-(b0/d)*x(idx,1)
        wprim(idx,8)=0.0d0
          else
        d0=1.0d0
        wprim(idx,6)=0.0d0
        wprim(idx,7)=0.0d0
        wprim(idx,8)=0.0d0
          endif
        wprim(idx,1)=d0
        wprim(idx,5)=p0
        
        cs(idx)=dsqrt(gam(idx)*p0/d0)
                
        call primitive_to_conservative(idx,wprim(idx,:),u(idx,:))
        
        enddo

        end select
        
        return

        end 
                
                
                
       subroutine init_magcollapse(wprim)
          
       use wp3d_h

   !------------------------------------------------------------------
   
   ! collapse of uniform density cloud embedded in 
   ! an ambient medium    
   ! uniform grid with equal mass particles

      implicit none
      integer::i,j,k,nint,ni,ncloud,nmedium,nt
      parameter(nint=100,ni=200000)
      double precision:: rho0,rcloud,r,xtest,ytest,ztest
      double precision:: gravity,boltzmann,protonmass,AU
      double precision:: msun,unitmass,unitvelocity,unittime
      double precision:: unitsurfacedensity,pbar
      double precision:: unitlength,unitdensity,unitenergy  
      double precision:: unitpressure,unitspecificenergy
      double precision:: tff,temp,radius,mass,om,mcloud
      double precision:: meanweight,Prot,mcrit,mu,c1, &
      unitBfield,mu0,dens,hini,parsec,rotratio,thetaoblique, &
      bx,bz,bmag
      double precision::wprim(nmax,nv)

       gravity=6.672d-8
       boltzmann=1.3806d-16
       protonmass=1.6726d-27
       meanweight=2.0d0*protonmass
       mu0=1.0d0
       gam(:)=1.0d0
       rhocrit1=1.0d-10
       rhocrit2=5.7d-5
       rhocrit3=1.0d0

       fname="magcol"
           
       xmin=-2.0d0
       xmax=2.0d0
       ymin=-2.0d0
       ymax=2.0d0
       zmin=-2.0d0
       zmax=2.0d0

       iperiodic=.true. 
       periodictypex=1
       periodictypey=1
       periodictypez=1
           
       boxx=xmax-xmin
       boxy=ymax-ymin
       boxz=zmax-zmin
       nmedium=50
       nt=nmedium**3
       voltot=boxx*boxy*boxz
       hini=0.1d0
       parsec=3.08567758d+18
       tf=0.01d0
       tprin=0.01d0 
	   
! temperature of the cloud in K

! physical constants in SI units

! the astronomical unit in cm

       msun=1.98892d+33

!      mass of the sun in g
!      mass of the cloud in solar mass

       mcloud=1.0d0
! mass of the cloud in units of the mass of the sun 
       mass=mcloud*msun
!  mass of the cloud in cgs units 
       temp=11.0d0
!  temperature 

       radius=0.016d0*parsec 

! Free fall time of the initially isothermal cloud

       rho0=3*mass/(4*pi*radius**3)
                  
       tff=dsqrt(radius**3/(gravity*mass)) 

       unittime=tff
       unitlength=radius
       unitmass=mcloud*msun 
       unitvelocity=unitlength/unittime
       unitdensity=unitmass/unitlength**3
       unitenergy=unitmass*unitlength**2/unittime**2
       unitpressure=unitmass/unitlength/unittime**2
       unitspecificenergy=unitenergy/unitmass
       unitsurfacedensity=unitmass/unitlength**2
       G=gravity*unitmass*unittime**2/unitlength**3
       unitBfield=dsqrt(mu0*unitdensity*unitvelocity**2)
                          
       rhocrit1=rhocrit1/unitdensity
       rhocrit2=rhocrit2/unitdensity
       rhocrit3=rhocrit3/unitdensity
	   
       rotratio=0.045d0 
       om=dsqrt(rotratio*3.0d0*(gravity*mass)/radius**3)
       Prot=2*pi/om/unittime
! rotational period of the cloud
	   
       mu=5
       c1=0.53d0
       bmag=3*mass/(2*pi*mu*c1*radius**2)* &
       (pi*gravity*mu0/5)**(1.0d0/2.0d0)
       bmag=bmag/unitbfield
       thetaoblique=0.0d0
! angle between the rotational axis and the direction of the magnetic field vector
 
       rho0=rho0/unitdensity
       c0=dsqrt(boltzmann*temp/meanweight)/unitvelocity
       rcloud=radius/unitlength
                 
! rotational period of the cloud 

      write(1,*) 'unit of time(yrs)=',unittime/(24*365*3600)
      write(1,*) 'unit of velocity=',unitvelocity
      write(1,*) 'unit of density=',unitdensity
      write(1,*) 'unit of pressure=',unitpressure
      write(1,*) 'unit of specific energy',unitspecificenergy
      write(1,*) 'gravity constant in code units=',G    
      write(1,*) 'rho0:',rho0*unitdensity
      write(1,*) 'initial temperature(Kelvin):',temp
      write(1,*) 'unitsurfacedensity:',unitsurfacedensity
      write(1,*) 'unit of magnetic field(Gauss)(µG):',unitBfield*1.0d+6
      write(1,*) 'initial magnetic field strength(µG):',bz*unitbfield*1.0d+6

! initial homogeneous density of the cloud 

       call setupem1(ni,ncloud)
      
       rhocrit4=(pi**3*(c0*unitvelocity)**6/(4*gravity**3*nopt**2))* &
       (float(ncloud)/mass)**2/unitdensity
     
       mcrit=pi**(3./2.)*(c0*unitvelocity)**3/(2*nopt*gravity**(3./2.)* &
       dsqrt(rhocrit1*unitdensity))/unitmass
        
       write(2,*) 'critical density( g cm-3)',rhocrit4*unitdensity
       write(2,*) 'critical particle mass',mcrit

       dx1=boxx/float(nmedium)
       dy1=dx1
       dz1=dx1
           
       exterior(:)=.false. 

        n=ncloud

        do i=1,ncloud
                
!       wprim(i,1)=dens(rho0,x(i,1),x(i,2))
        wprim(i,1)=rho0
        wprim(i,2)=-2*pi*x(i,2)/Prot
        wprim(i,3)=2*pi*x(i,1)/Prot
        wprim(i,4)=0.0d0
                   
        wprim(i,5)=pbar(i,wprim(i,1))
		
        bx=bmag*dsin(thetaoblique)
        bz=bmag*dcos(thetaoblique)
          
        wprim(i,6)=bx
        wprim(i,7)=0.0d0
        wprim(i,8)=bz
                   
        h(i)=hini            

       enddo
         
       do i=1,nmedium
       do j=1,nmedium
       do k=1,nmedium

       xtest=(i-1)*dx1+xmin+dx1/2
       ytest=(j-1)*dy1+ymin+dy1/2
       ztest=(k-1)*dz1+zmin+dz1/2
        
       r=dsqrt(xtest**2+ytest**2+ &
       ztest**2)
           
       if ( r .ge. rcloud ) then

       n=n+1
           
       exterior(n)=.true.

       x(n,1)=xtest
       x(n,2)=ytest
       x(n,3)=ztest
           
       wprim(n,1)=rho0/30.0d0
          
       wprim(n,2:4)=0.0d0 

       wprim(n,5)=pbar(n,wprim(n,1))        
          
       wprim(n,6)=0.0d0
       wprim(n,7)=0.0d0
       wprim(n,8)=bz
           
       h(n)=hini

      endif
          
      enddo
      enddo
      enddo     

      do i=1,n

      cs(i)=dsqrt(wprim(i,5)/wprim(i,1))
               
      enddo
               
      write(2,*) 'number of particles inside the cloud:',ncloud
  
      return
          
      end

       double precision function pbar(i,rhoi)

        use wp3d_h

        implicit none
        integer::i
        double precision::rhoi,k1,k2,k3
                
        k1=c0**2/rhocrit1**(2.0d0/5.0d0)
        k2=k1*rhocrit2**(1.0d0/4.0d0)
        k3=k2*rhocrit3*(1.15d0-(5.0d0/3.0d0))
                
! barotropic EOS, includes effects related to the first and second collapse phase
                
        if ( rhoi .le. rhocrit1 ) then
        pbar=c0**2*rhoi
        gam(i)=1.0d0
        endif
        
        if (( rhoi .gt. rhocrit1 ) .and. ( rhoi .le. rhocrit2 )) then
        gam(i)=7.0d0/5.0d0              
        pbar=k1*rhoi**gam(i)
        endif
                
        if (( rhoi .gt. rhocrit2 ) .and. ( rhoi .le. rhocrit3 )) then
        gam(i)=1.15d0           
        pbar=k2*rhoi**gam(i)
        endif
                
        if ( rhoi .gt. rhocrit3 ) then
        gam(i)=5.0d0/3.0d0              
        pbar=k3*rhoi**gam(i)
        endif
        
        return

        end
  
  
       double precision function dens(rho0,xp,yp)
        
       double precision::rho0,xp,yp,ph

! Boss and Bodenheimer initial density distribution

       if ( xp .eq. 0.0d0) then

       ph=0.0d0

       else

       ph=datan(yp/xp)

       endif

       dens=rho0*(1.0d0+0.1d0*dcos(2*ph))

       return 

       end
        

      SUBROUTINE setupem1(ni,ncloud)  

      use wp3d_h  
          
      IMPLICIT NONE

! sets up a particle distribution in 3D 
! particles are placed on a stretched grid 
                        
      double precision:: avec(3),bvec(3),cvec(3)
      double precision:: space,xp,yp,zp
      double precision:: r,rns
      INTEGER::maxn,i,idum,k,l,m,ni,ncloud
              
! setup particles on a hexagonal lattice

      rns=10.0d0
 !    mcloud=1.0d0
 !    amns=mcloud
      avec(1)=1.0d0
      avec(2)=0.0d0
      avec(3)=0.0d0
      bvec(1)=0.5d0
      bvec(2)=0.5d0*dsqrt(3.0d0)
      bvec(3)=0.0d0
      cvec(1)=0.5d0
      cvec(2)=1.0d0/dsqrt(12.0d0)
      cvec(3)=dsqrt(2.0d0/3.0d0)
      space=rns/(2.0d0*ni/2.0d0/6.0d0)**0.333d0

! The factor 1.23 makes sure to catch all points

      maxn=int(1.23*rns/space)
      write (6,*)' maxn: ',maxn,' spacing: ',space
      i=0
      idum=-2391
     
       do k=-maxn,maxn
         do l=-maxn,maxn
            do m=-maxn,maxn
               xp=space*(k*avec(1)+l*bvec(1)+m*cvec(1))
               yp=space*(k*avec(2)+l*bvec(2)+m*cvec(2))
               zp=space*(k*avec(3)+l*bvec(3)+m*cvec(3))
               r=sqrt(xp**2+yp**2+zp**2)
                  if (r.lt.rns) then
                  i=i+1
                  if(i.le.ni) then
                  
! Apply an azimuthal density perturbation(cf. Boss and Bodenheimer)
                                                    
                    x(i,1)=xp
                    x(i,2)=yp
                    x(i,3)=zp

                    else
                     write (6,*)'1es: Reached ni!'
                     stop
                  endif
               endif
            enddo
         enddo
      enddo
      ncloud=i
      write (6,*)'number of particles within rns: ',i
       
      return
                    
      END       
          
          
       
      subroutine evrard_ini(u,wprim)

   !------------------------------------------------------------------
   
   ! initial conditions 
   ! Evrard test problem 
  
        use wp3d_h

        implicit none

        integer:: i,j,k,nint,ni,nmedium,ncloud
        parameter(nint=100,ni=50000,nmedium=40)
        double precision gravity,boltzmann,protonmass,AU
        double precision msun,unitmass,unitvelocity,unittime
        double precision unitsurfacedensity
        double precision unitlength,unitdensity,unitenergy  
        double precision unitpressure,unitspecificenergy
        double precision tff,radius,mass,csound,bz
        double precision utherm,r,dens_evrard,rho0, &
        xtest,ytest,ztest,rcloud,hini
        double precision::u(nmax,nv),wprim(nmax,nv)

        gravity=6.672d-8
        boltzmann=1.3806d-16
        protonmass=1.6726d-24
        gam(:)=5.0d0/3.0d0
        hini=0.1d0
        iperiodic=.true. 
        periodictypex=1
        periodictypey=1
        periodictypez=1
          
        xmin=-2.0d0
        xmax=2.0d0
        ymin=-2.0d0
        ymax=2.0d0
        zmin=-2.0d0
        zmax=2.0d0
        boxx=xmax-xmin
        boxy=ymax-ymin
        boxz=zmax-zmin
		
        dx1=boxx/float(nmedium)
        dy1=boxy/float(nmedium)
        dz1=boxz/float(nmedium)

        fname="evrard"
        
! physical constants in cgs units

        AU=1.49598d13

! the astronomical unit in cm 

        msun=1.98892d+33

! mass of the sun in g 
! mass of the cloud in solar mass
! internal unit system of the code 

        radius=100.0*AU
        mass=msun
        tff=dsqrt(radius**3/(gravity*mass))

! Free fall time of the initially isothermal cloud 

        unittime=tff
        unitlength=radius
        unitmass=mass
        unitvelocity=unitlength/unittime
        unitdensity=unitmass/unitlength**3
        unitenergy=unitmass*unitlength**2/unittime**2
        unitpressure=unitmass/unitlength/unittime**2
        unitspecificenergy=unitenergy/unitmass
        unitsurfacedensity=unitmass/unitlength**2
        G=gravity*unitmass*unittime**2/unitlength**3
		
! system of units adopted by Evrard
! except for a constant factor in the unit of time 

         rho0=mass/(2*pi*(radius**3))/unitdensity
         utherm=0.05d0*gravity*mass/radius/unitspecificenergy 
         rcloud=radius/unitlength
         tf=2.0d0
         tprin=0.05d0

! initial temperature

         csound=dsqrt(gam(1)*(gam(1)-1.0d0)*utherm)

         write(1,*) 'unit of time(yrs)=',unittime/(24*3600*365)
         write(1,*) 'unit of velocity=',unitvelocity
         write(1,*) 'unit of density=',unitdensity
         write(1,*) 'unit of pressure=',unitpressure
         write(1,*) 'unit of specific energy',unitspecificenergy
         write(1,*) 'gravity constant in code units=',G
         write(1,*) 'specific energy in code units=',utherm
         write(1,*) 'rho0(code units):',rho0
         write(1,*) 'initial sound velocity(code units):',csound

! initial homogeneous density of the cloud 

         call setupem2(ni,ncloud)

         do i=1,ncloud

         r=dsqrt(x(i,1)**2+x(i,2)**2+x(i,3)**2)
         h(i)=hini
         wprim(i,1)=dens_evrard(r)
         wprim(i,5)=(gam(i)-1.0d0)*wprim(i,1)*utherm
         wprim(i,2:4)=0.0d0
         wprim(i,6:8)=0.0d0

         enddo

        n=ncloud

        do i=1,nmedium
        do j=1,nmedium
        do k=1,nmedium

        xtest=(i-1)*dx1+xmin+dx1/2
        ytest=(j-1)*dy1+ymin+dy1/2
        ztest=(k-1)*dz1+zmin+dz1/2
        
        r=dsqrt(xtest**2+ytest**2+ztest**2)
           
        if ( r .ge. rcloud ) then

        n=n+1
           
        exterior(n)=.true.
		
        x(n,1)=xtest
        x(n,2)=ytest
        x(n,3)=ztest
           
        wprim(n,1)=rho0/1000
 
        wprim(n,2:4)=0.0d0 

        wprim(n,5)=(gam(n)-1.0d0)*wprim(n,1)*utherm*1000.0d0
		
	    csound=dsqrt(gam(1)*(gam(1)-1.0d0)*utherm*1000.0d0)
          
        wprim(n,6:8)=0.0d0

        h(n)=hini
       
        cs(n)=csound

        endif

        enddo
        enddo
        enddo
      
        do i=1,n
        call primitive_to_conservative(i,wprim(i,:),u(i,:))
        enddo  

        return

        end


        double precision function dens_evrard(r)

        use wp3d_h

        double precision::r

! Density profile for the Evrard gas sphere 

        dens_evrard=1.0d0/r/2/pi

        return 

        end


        SUBROUTINE setupem2(ni,ncloud)

        use wp3d_h
                                          
        IMPLICIT NONE

! sets up a particle distribution in 3D 
! particles are placed on a stretched grid 

      integer ni,ncloud
      double precision r                              
      double precision avec(3),bvec(3),cvec(3)
      double precision space,xrand,xp,yp,zp
      double precision nuniform,runi,rdis,rns,amns,mpar
      integer maxn,i,idum,k,l,m
            
! Lay down particles on a hexagonal lattice

      rns=1.0d0
      amns=1.0d0
      avec(1)=1.0d0
      avec(2)=0.0d0
      avec(3)=0.0d0
      bvec(1)=0.5d0
      bvec(2)=0.5d0*dsqrt(3.0d0)
      bvec(3)=0.0d0
      cvec(1)=0.5d0
      cvec(2)=1.0d0/dsqrt(12.0d0)
      cvec(3)=dsqrt(2.0d0/3.0d0)
      space=rns/(1.0d0*ni/2.0d0/6.0d0)**0.333d0
     
! The factor 1.23 makes sure to catch all points

      maxn=int(1.23d0*rns/space)
      write (6,*)' maxn: ',maxn,' spacing: ',space
      i=0
      idum=-2391
      xrand=0.001*space

      do k=-maxn,maxn
         do l=-maxn,maxn
            do m=-maxn,maxn
               xp=space*(k*avec(1)+l*bvec(1)+m*cvec(1))
               yp=space*(k*avec(2)+l*bvec(2)+m*cvec(2))
               zp=space*(k*avec(3)+l*bvec(3)+m*cvec(3))
               r=sqrt(xp**2+yp**2+zp**2)
               if ((r.lt.rns) .and. ( r .gt. 0.0d0 )) then                              
                  i=i+1           
                  if(i.le.ni) then
                     x(i,1)=xp
                     x(i,2)=yp
                     x(i,3)=zp
                     else
                     write (6,*)'1es: Reached ni!'
                     stop
                  endif
               endif
            enddo
         enddo
       enddo
       ncloud=i
	   
       write (6,*) 'number of particles within rns: ',i
	   
       if ( stretch .eqv. .true. ) then

       mpar=amns/float(ncloud)
       
! stretching of the grid to match the required density profile   

       nuniform=3*float(ncloud)/(4*pi*rns**3)

       do i=1,ncloud
          
       runi=dsqrt(x(i,1)**2+x(i,2)**2+x(i,3)**2)
       rdis=dsqrt(2*runi**3*nuniform*mpar/3)
 
! adjusting particle positions 

       x(i,1)=x(i,1)*rdis/runi
       x(i,2)=x(i,2)*rdis/runi
       x(i,3)=x(i,3)*rdis/runi

       enddo
	   
       endif

       return
                                                                   
      END              


        subroutine create_particle_distribution
                
        use wp3d_h
                
        implicit none
                
        integer::i,j,k,idx,niter,idum
        double precision:: ran1,hini

        select case (setupchoice) 
                
        case (1)
                
         dx1=dabs(boxx)/float(nx1)
         dy1=dx1
         dz1=dx1
                
         idx=0
         do i=0,nx1-1
         do j=0,ny1-1
         do k=0,nz1-1
                
         idx=idx+1

         x(idx,1)=float(i)*dx1+xcmin+dx1/2
         x(idx,2)=float(j)*dy1+ycmin+dy1/2
         x(idx,3)=float(k)*dz1+zcmin+dz1/2

         enddo
         enddo
         enddo
                 
         n=idx

         case (2)
        
! create a random particle distribution

         idum=-100
         niter=10
         hini=0.01d0

         do i=1,n
                
         x(i,1)=xmin+(ran1(idum)+0.5d0)*boxx
         x(i,2)=ymin+(ran1(idum)+0.5d0)*boxy
         x(i,3)=zmin+(ran1(idum)+0.5d0)*boxz
         h(i)=hini
                
         enddo
 
!         if ( relax .eqv. .true. ) then
!         call relaxparticles(niter,wprim)
!         endif
                
! relax the particle distribution
                
         end select 
                
         return
                
         end
                
                
        subroutine relaxparticles(niter,wprim)
                
! subroutine to relax a random particle distribution 
                
           USE wp3d_h
                
           implicit none
                
           integer::it,i,j,in,niter,iopt,itab
           double precision:: wij,dri(nmax,ndim),deltar, &
           h3i,h2i,d2,u2,av,dx,dy,dz,periodic1,periodic2, &
           wprim(nmax,nv)
           
           iopt=1
           av=0.0001d0
                
           do it=1,niter
                
           call build_linked_list
           call neighbourlists(wprim)
            
           deltar=0.0d0
                
           do i=1,n

           dri(i,:)=0.0d0
                   
           call LLN(i,iopt,h(i))
             
           H2I=H(I)**2     
           H3I=H2I*h(i)
        
!  Compute dR for each particle

           DO IN=1,NBN
         
            J=NN(IN)

            if ( iperiodic ) then
                   
            dx=x(j,1)-x(i,1)
            dy=x(j,2)-x(i,2)
            dz=x(j,3)-x(i,3)
                 
            call modbound(dx,dy,dz)

            endif

            d2=dx**2+dy**2+dz**2
        
           U2=D2/H2I
           IF ( U2 .ge. 4.0d0 ) THEN
           WIJ=0.0d0
           ELSE
           ITAB=INT(CTAB*U2)+1
           WIJ=WTAB(ITAB)/H3I
           ENDIF 
           dri(i,1)=dri(i,1)+WIJ*dx
           dri(i,2)=dri(i,2)+WIJ*dy
           dri(i,3)=dri(i,3)+WIJ*dz

         enddo
                 
         deltar=deltar+dsqrt(dri(i,1)**2+dri(i,2)**2+dri(i,3)**2)
                 
         enddo

        do i=1,n
        x(i,:)=x(i,:)-av*dri(i,:)
        enddo

        if ( iperiodic ) then
        
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

        enddo
                
        return
                
        end
                
                
       DOUBLE PRECISION FUNCTION ran1(idum)
       INTEGER:: idum,IA,IM,IQ,IR,NTAB,NDIV
       DOUBLE PRECISION AM,EPS,RNMX
       PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,  &
       NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
       INTEGER:: j,k,iv(NTAB),iy
       SAVE iv,iy
       DATA iv /NTAB*0/, iy /0/
       if (idum.le.0.or.iy.eq.0) then 
       idum=max(-idum,1)
       do j=NTAB+8,1,-1
       k=idum/IQ
       idum=IA*(idum-k*IQ)-IR*k
       if (idum.lt.0) then 
       idum=idum+IM
       endif
       if (j.le.NTAB) then
       iv(j)=idum
       endif
       enddo
       iy=iv(1)
       endif
       k=idum/IQ
       idum=IA*(idum-k*IQ)-IR*k
       if (idum.lt.0) then 
       idum=idum+IM
       endif
       j=1+iy/NDIV
       iy=iv(j)
       iv(j)=idum
       ran1=min(AM*iy,RNMX)-0.5
       return
       END       
                
        
        SUBROUTINE PDUMP(wprim) 

        use wp3d_h 

        integer::i
        character*40::filename
        double precision::wprim(nmax,nv),rmax

        ndump=ndump+1

!       write(6,*) 'writing dump file at t=',t
        write(2,*) 'writing dump file at t=',t
        write(11,*) 'particle dump',ndump,'at t=',t
        
        write(unit=filename, fmt="(A,I0.3,A)") trim(outfn(11)), &
        ndump,".dat"

        nacc=0
        nsink=0

        open(12,FILE=filename)
        
        write(12,fmt="(3I8)") n,nacc,nsink
        write(12,fmt="(2ES25.15)") t,gam(1)
     
        do i=1,n
        
! interpolate particle positions and velocities to the current system time 
        
        write(12,fmt="(14ES25.15)") x(i,1),x(i,2),x(i,3),wprim(i,2), &
        wprim(i,3),wprim(i,4),wprim(i,6),wprim(i,7),wprim(i,8),am(i),h(i), &
        wprim(i,5),wprim(i,1),eb(i)
       
! SPH particle dump 

        enddo
        
        close(12)
        
        return 
        END  
        
                        
        SUBROUTINE plot(wprim)
        
        use wp3d_h
        
        integer:: i,j,k,dim,iopt,idx,nscal,nvec
        double precision,allocatable::xgrid(:),ygrid(:),zgrid(:)
        
        double precision:: tiny
        double precision:: wprim(nmax,nv)
        double precision, allocatable:: rhom(:,:,:), &
        velm(:,:,:,:),ptherm(:,:,:),bm(:,:,:,:)
        double precision,allocatable::divbm(:,:,:)
        double precision,allocatable::scaldata(:,:,:,:), &
        vecdata(:,:,:,:,:)
        
        character*20:: filename
        character*10,allocatable::names(:)
        
        iopt=1
        tiny=1.0d-10
        ndump=ndump+1

        allocate(rhom(nx1,ny1,nz1))
        allocate(velm(nx1,ny1,nz1,ndim))
        allocate(ptherm(nx1,ny1,nz1))
        if ( imhd ) then
        allocate(bm(nx1,ny1,nz1,ndim))
        allocate(divbm(nx1,ny1,nz1))
        endif
                
        allocate(xgrid(nx1),ygrid(ny1),zgrid(nz1))

        if ( setupchoice .eq. 2 ) then
                
        dx1=boxx/float(nx1)
        dy1=boxy/float(ny1)
        dz1=boxz/float(nz1)
                
        endif
                
        write(unit=filename, fmt="(A,I0.3,A)") trim(outfn(10)), &
       ndump,".vtr"
     
        print *, filename
                 
          idx=0
                           
          do i=0,nx1-1
          
          xgrid(i+1)=dx1*float(i)+dx1/2+xcmin
                 
          do j=0,ny1-1
          
          ygrid(j+1)=dy1*float(j)+dy1/2+ycmin
                  
          do k=0,nz1-1

          idx=idx+1
          
          zgrid(k+1)=dz1*float(k)+dz1/2+zcmin
                  
         rhom(i+1,j+1,k+1)=wprim(idx,1)
         ptherm(i+1,j+1,k+1)=wprim(idx,5)

         velm(i+1,j+1,k+1,1)=wprim(idx,2)
         velm(i+1,j+1,k+1,2)=wprim(idx,3)
         velm(i+1,j+1,k+1,3)=wprim(idx,4)

         if ( imhd ) then

         bm(i+1,j+1,k+1,1)=wprim(idx,6)
         bm(i+1,j+1,k+1,2)=wprim(idx,7)
         bm(i+1,j+1,k+1,3)=wprim(idx,8)

          if ( imhd ) then
          divbm(i+1,j+1,k+1)=dabs(divb(idx))
          else
          divbm(i+1,j+1,k+1)=0.0d0
          endif

           endif
         
          enddo
          enddo
           enddo
                  
        if ( imhd ) then
                
        nscal=3
        nvec=2
                
        allocate(names(nscal+nvec))
        allocate(scaldata(1:nx1,1:ny1,1:nz1,1:nscal))
        allocate(vecdata(1:nx1,1:ny1,1:nz1,1:ndim,1:nvec))
          
        names(1)="rho"
        names(2)="pth"
        names(3)="divb"
        names(4)="v"
        names(5)="b"
    
! write rendered plot to a vtk output file
       
        do i=1,nx1
        do j=1,ny1
        do k=1,nz1
        
! density and thermal pressure
                 
        scaldata(i,j,k,1)=dmax1(rhom(i,j,k),dfloor)
        scaldata(i,j,k,2)=dmax1(ptherm(i,j,k),pfloor)
        scaldata(i,j,k,3)=divbm(i,j,k)
     
         do dim=1,ndim

! components of the velocity field  
        
        vecdata(i,j,k,dim,1)=velm(i,j,k,dim)
       
! components of the magnetic field
        
        vecdata(i,j,k,dim,2)=bm(i,j,k,dim)
                
        enddo
                
                enddo
                enddo
                enddo
                
                else
                
                nscal=2
                nvec=1
                
                allocate(names(nscal+nvec))
                allocate(scaldata(1:nx1,1:ny1,1:nz1,1:nscal))
                allocate(vecdata(1:nx1,1:ny1,1:nz1,1:ndim,1:nvec))
         
                names(1)="rho"
        names(2)="pth"
        names(3)="v"
      
! write rendered plot to a vtk output file
       
        do i=1,nx1
        do j=1,ny1
        do k=1,nz1
        
! density and thermal pressure
                 
 !     scaldata(i,j,k,1)=dmax1(rhom(i,j,k),dfloor)
  !    scaldata(i,j,k,2)=dmax1(ptherm(i,j,k),pfloor)
        
! components of the velocity field  

        do dim=1,ndim    
        
        vecdata(i,j,k,dim,1)=velm(i,j,k,dim)
                
        enddo
                         
        enddo
        enddo
        enddo
                
        endif
        
        call vtk_3d(nx1,ny1,nz1,xgrid,ygrid,zgrid,nscal,nvec, &
        scaldata,vecdata,names,filename)
        
        deallocate(rhom,velm,ptherm,xgrid,ygrid,zgrid)
        deallocate(vecdata,scaldata)

        if ( imhd ) then

         deallocate(bm,divbm)

         endif
                                           
      return 
      END   


      SUBROUTINE render3D(wprim)
      
 ! produces a grid interpolation of the SPH particle data for use with paraview 
 ! and visit 
      
       USE wp3d_h
      
      integer:: i,j,k,in,dim,iopt,itab
      double precision::dx,dy,dz,xtry,ytry,ztry,h2,h3, &
      hsmooth,r2,u2,hsearch
      character*1 ndx1
      character*2 ndx2
      character*3 ndx3
      character(len=20):: fname1
      integer,parameter:: ngridx=60,ngridy=60,ngridz=30,nscal=2, &
      nvec=2  
      double precision,allocatable::rhom(:,:,:), &
      pm(:,:,:),velm(:,:,:,:),bm(:,:,:,:),xgrid(:),ygrid(:), &
      zgrid(:),scaldata(:,:,:,:),vecdata(:,:,:,:,:)
      double precision::wprim(nmax,nv)
      character(len=10):: names(1:nscal+nvec)
      
      allocate(xgrid(1:ngridx),ygrid(1:ngridy),zgrid(1:ngridz))
      allocate(rhom(1:ngridx,1:ngridy,1:ngridz))
      allocate(pm(1:ngridx,1:ngridy,1:ngridz))
      allocate(velm(1:ngridx,1:ngridy,1:ngridz,1:ndim))
      allocate(bm(1:ngridx,1:ngridy,1:ngridz,1:ndim))
      allocate(scaldata(1:ngridx,1:ngridy,1:ngridz,1:nscal))
      allocate(vecdata(1:ngridx,1:ngridy,1:ngridz,1:ndim,1:nvec))
      
      iopt=1
      hsearch=0.25d0
      
      xcmin=-2.0d0
      xcmax=2.0d0
      ycmin=-2.0d0
      ycmax=2.0d0
      zcmin=-1.0d0
      zcmax=1.0d0
      ndump=ndump+1
      
      if ( ndump .lt. 10 ) then

        write(ndx1,'(I1)') ndump 

        ndx3="00" // ndx1

        else

        if ( ndump .lt. 100 ) then

        write(ndx2,'(I2)') ndump 

        ndx3="0" // ndx2

        else 
        
        write(ndx3,'(I3)') ndump 

        endif 

        endif 

        fname1=trim(outfn(12)) // ndx3 // ".vtr"
        
        dx=(xcmax-xcmin)/float(ngridx)
        dy=(ycmax-ycmin)/float(ngridy)
        dz=(zcmax-zcmin)/float(ngridz)
        
        do i=1,ngridx
        
        xgrid(i)=float(i-1)*dx+xcmin+dx/2.
        
        enddo
        
        do i=1,ngridy
        
        ygrid(i)=float(i-1)*dy+ycmin+dy/2.
        
        enddo
        
        do i=1,ngridz
        
        zgrid(i)=float(i-1)*dz+zcmin+dz/2.
        
        enddo
        
    ! interpolate particle properties to a 3D grid
        
        do i=1,ngridx
        
        xtry=xgrid(i) 
        
        do j=1,ngridy
          
        ytry=ygrid(j)
                 
        do k=1,ngridz
          
        ztry=zgrid(k) 
         
        call LLN2(xtry,ytry,ztry,iopt,hsearch)
        
        rhom(i,j,k)=0.0d0  
        pm(i,j,k)=0.0d0 
        do dim=1,ndim 
        velm(i,j,k,dim)=0.0d0
        bm(i,j,k,dim)=0.0d0
        enddo
           
        hsmooth=0.0d0

        do in=1,nbn

        hsmooth=dmax1(hsmooth,h(nni(i,in)),dx/2.0d0)
        
        enddo

        h2=hsmooth**2
        h3=hsmooth*h2
        
        do in=1,nbn
          
        r2=(x(nni(i,in),1)-xtry)**2+(x(nni(i,in),2)-ytry)**2+(x(nni(i,in),3)-ztry)**2 
        
        if ( r2 .le. 4*h2 ) then 
                                         
!     We calculate the density contribution of all particles to all
!     the grid points with which their kernel function overlaps
        
        u2=r2/h2

        ITAB=INT(CTAB*U2)+1

        rhom(i,j,k)=rhom(i,j,k)+am(nni(i,in))*wtab(itab)/h3
     
        pm(i,j,k)=pm(i,j,k)+am(nni(i,in))*(wprim(nni(i,in),5)/wprim(nni(i,in),1))* &
        wtab(itab)/h3

        do dim=1,ndim

        velm(i,j,k,dim)=velm(i,j,k,dim)+am(nni(i,in))*(wprim(nni(i,in),dim+1)/wprim(nni(i,in),1))* &
        wtab(itab)/h3

        bm(i,j,k,dim)=bm(i,j,k,dim)+am(nni(i,in))*(wprim(nni(i,in),dim+5)/wprim(nni(i,in),1))* &
        wtab(itab)/h3

        enddo
      
        endif
           
        enddo
       
                              
        enddo
          
        enddo
      
        enddo
                
        names(1)="rho"
        names(2)="pth"
        names(3)="v"
        names(4)="b"
    
! write rendered plot to a vtk output file
       
        do i=1,ngridx
        do j=1,ngridy
        do k=1,ngridz
        
! density and thermal pressure
                 
 !       scaldata(i,j,k,1)=dmax1(rhom(i,j,k),dfloor)
 !       scaldata(i,j,k,2)=dmax1(pm(i,j,k),pfloor)
        
! components of the velocity field  

        do dim=1,ndim    
        
        vecdata(i,j,k,dim,1)=velm(i,j,k,dim)
       
! components of the magnetic field
        
        vecdata(i,j,k,dim,2)=bm(i,j,k,dim)
                
        enddo
                      
        enddo
        enddo
        enddo
        
        call vtk_3d(ngridx,ngridy,ngridz,xgrid,ygrid,zgrid,nscal,nvec, &
        scaldata,vecdata,names,fname1)
        
        deallocate(xgrid,ygrid,zgrid,rhom,pm,velm,bm,scaldata,vecdata)
        
             
      RETURN
      
      end
      
         
       subroutine vtk_3d(n1m,n2m,n3m,y1,y2,y3,nscal,nvec,scaldata, &
       vecdata,names,filename)
       
!      Export 3D data upon a uniform rectilinear Grid      
!      in binary vtk format for use with paraview and visit visualisation software

!      n1m,n2m,n3m: size of the data and the grid
!      y1,y2,y3: coordinates
!      nscal: number of scalar variables
!      nvec: number of vector variables
!      scal_data: scalar data to write
!      vec_data:: vector data to write
!      names: name of the data  

       use wp3d_h   

       INTEGER,INTENT(IN)::n1m,n2m,n3m
       double precision,DIMENSION(n1m),INTENT(IN)::y1
       double precision,DIMENSION(n2m),INTENT(IN)::y2
       double precision,DIMENSION(n3m),INTENT(IN)::y3
       INTEGER,INTENT(IN):: nscal,nvec
       double precision,DIMENSION(n1m,n2m,n3m,nscal),INTENT(IN)::scaldata
       double precision,DIMENSION(n1m,n2m,n3m,ndim,nvec),INTENT(IN)::vecdata
       CHARACTER(LEN=10),INTENT(IN)::names(nscal+nvec)
       CHARACTER(LEN=20),INTENT(IN)::filename
       
       integer::i,j,k,s,nfil

        nfil=41
        open(nfil,file=filename)
        write(6,*) 'writing vtk output file ',filename
        write(nfil,*)'<VTKFile type="RectilinearGrid" version="0.1"', &
             ' byte_order="LittleEndian">'
          write(nfil,*)'  <RectilinearGrid WholeExtent=', &
                    '"1 ',n1m,' 1 ',n2m,' 1 ',n3m,'">'
          write(nfil,*)'    <Piece Extent=', &
                      '"1 ',n1m,' 1 ',n2m,' 1 ',n3m,'">'
          write(nfil,*)'      <Coordinates>'
          write(nfil,*)'        <DataArray type="Float32"', &
                                ' Name="X_COORDINATES"', &
                                ' NumberOfComponents="1">' 
          write(nfil,*) (y1(i),i=1,n1m)
          write(nfil,*)'        </DataArray>'
          write(nfil,*)'        <DataArray type="Float32"', &
                                 ' Name="Y_COORDINATES"', &
                                 ' NumberOfComponents="1">'
          write(nfil,*) (y2(j),j=1,n2m)
          write(nfil,*)'        </DataArray>'
          write(nfil,*)'        <DataArray type="Float32"', &
                                 ' Name="Z_COORDINATES"', &
                                 ' NumberOfComponents="1">'
          write(nfil,*) (y3(k),k=1,n3m)
          write(nfil,*)'        </DataArray>'
          write(nfil,*)'      </Coordinates>'  
             
          write(nfil,*)'      <PointData Scalars="rho" Vectors="v">'

! output of scalar data           
          
          do s=1,nscal 
          
          write(nfil,*)'       <DataArray Name="'//trim(names(s))//'"', &
                               ' type="Float32"', &
                               ' NumberOfComponents="1"', &
                               ' format="ascii">' 
          write(nfil,*) (((scaldata(i,j,k,s),i=1,n1m),j=1,n2m),k=1,n3m)
          write(nfil,*)'        </DataArray>'
         
          enddo
          
! output of vector data   
          
          do s=1,nvec
                   
          write(nfil,*)'       <DataArray Name="'//trim(names(s+nscal))//'"', &
                               ' type="Float32"', &
                               ' NumberOfComponents="3"', &  
                               ' format="ascii">' 
          write(nfil,*) (((vecdata(i,j,k,1,s),vecdata(i,j,k,2,s), &
                  vecdata(i,j,k,3,s),i=1,n1m),j=1,n2m),k=1,n3m)
          write(nfil,*)'        </DataArray>'
          
          enddo
        
          write(nfil,*)'      </PointData>'
                  
          write(nfil,*)'    </Piece>'
          write(nfil,*)'  </RectilinearGrid>'
          write(nfil,*)'</VTKFile>'
          close(nfil)
  
       return
       end   
           
           
       subroutine silo_3d(wprim)
       
!      Export particle data on a point mesh in binary silo 
!      format for use with visit visualisation software

       use wp3d_h   
 
!       include "silo.inc"

       INTEGER:: i,len1,len2
       CHARACTER(LEN=20)::filename
       CHARACTER(LEN=30)::comment
       
       integer::ierr,dbfile,err,ndims
      
       double precision,allocatable::rho(:),pth(:),vx(:),vy(:), &
       vz(:),bx(:),by(:),bz(:),xp(:),yp(:),zp(:)
       double precision::wprim(nmax,nv)
           
       allocate(rho(n),pth(n),vx(n),vy(n),vz(n), &
       bx(n),by(n),bz(n),xp(n),yp(n),zp(n)) 
           
!       ierr=DBShowErrors(DB_ALL,0)
          
       ndump=ndump+1

        comment="point mesh with wpmhd data"
        write(unit=filename, fmt="(A,I0.3,A)") trim(outfn(10)), &
        ndump,".silo"
     
      print *, filename

      len2=len(comment)
      len1=len(filename)

!     create silo file

!     ierr = dbcreate(filename, len1,DB_CLOBBER, DB_LOCAL,  &
!     comment, len2, DB_PDB, dbfile)

!  Write output-particle distribution with cylindrical symmetry

       ndims = 3
           
       do i=1,n

           xp(i)=x(i,1)   
           yp(i)=x(i,2)    
           zp(i)=x(i,3)    
           rho(i)=wprim(i,1)
           vx(i)=wprim(i,2)
           vy(i)=wprim(i,3)
           vz(i)=wprim(i,4)
           pth(i)=wprim(i,5)
           bx(i)=wprim(i,6)
           by(i)=wprim(i,7)
           bz(i)=wprim(i,8)
      
       enddo

 !     Write a point mesh
 
 !      err = dbputpm (dbfile, "pointmesh", 9, ndims, xp,yp,zp, &
 !          n, DB_DOUBLE, DB_F77NULL, ierr)
          
 !      err = dbputpv1(dbfile, "rho", 3, "pointmesh", 9, &
 !     rho, n, DB_DOUBLE, DB_F77NULL, ierr)
 !     err = dbputpv1(dbfile, "pth", 3, "pointmesh", 9, &
 !      pth, n, DB_DOUBLE, DB_F77NULL, ierr)
 
!      if ( imhd ) then
 !      err = dbputpv1(dbfile, "eb", 2, "pointmesh", 9, &
 !     eb, n, DB_DOUBLE, DB_F77NULL, ierr)
 !     endif

!      err = dbputpv1(dbfile, "vx", 2, "pointmesh", 9, &
!       vx, n, DB_DOUBLE, DB_F77NULL, ierr)
!      err = dbputpv1(dbfile, "vy", 2, "pointmesh", 9, &
!      vy, n, DB_DOUBLE, DB_F77NULL, ierr)
!      err = dbputpv1(dbfile, "vz", 2, "pointmesh", 9, &
!      vz, n, DB_DOUBLE, DB_F77NULL, ierr)
           
!       if ( imhd ) then
!      err = dbputpv1(dbfile, "bx", 2, "pointmesh", 9, &
!       bx, n, DB_DOUBLE, DB_F77NULL, ierr)
!       err = dbputpv1(dbfile, "by", 2, "pointmesh", 9, &
!       by, n, DB_DOUBLE, DB_F77NULL, ierr)
!       err = dbputpv1(dbfile, "bz", 2, "pointmesh", 9, &
!       bz, n, DB_DOUBLE, DB_F77NULL, ierr)
!       endif
   
       deallocate(rho,pth,vx,vy,vz,bx,by,bz,xp,yp,zp)
           
 !      ierr=dbclose(dbfile)
     
       return
       end 
           
                                       
       SUBROUTINE ENOUT(u,wprim)
           
       use wp3d_h
           
       IMPLICIT NONE
           
! *****************************************************
!     Calculates energy-related quantities, and writes summary to screen
!     and file "status.sph"
! *****************************************************
                                     
       double precision::xlow,ylow,xhigh,yhigh,zlow,zhigh
       double precision::ekin,eint,etot,emag,tiny,epot
       double precision::u(nmax,nv),wprim(nmax,nv)
       double precision::amvec(ndim),pvec(ndim)
       double precision::rhomin,rhomax,hpmax,hpmin,hpi,v2i
       double precision::ptot,jtot,mtot,ebmin,ebmax,ebave
       integer::i,irhomax,irhomin,iebmin,iebmax
       integer::ihpmax,ihpmin,nnmax,nnmin,innmin,innmax
      
!     Find system box 
      XLOW=1.E30
      YLOW=1.E30
      ZLOW=1.E30
      XHIGH=-1.E30
      YHIGH=-1.E30
      ZHIGH=-1.E30
      tiny=1.0d-10

      DO I=1,n
         XLOW=MIN(X(I,1),XLOW)
         YLOW=MIN(X(I,2),YLOW)
         ZLOW=MIN(X(I,3),ZLOW)
         XHIGH=MAX(X(I,1),XHIGH)
         YHIGH=MAX(X(I,2),YHIGH) 
         ZHIGH=MAX(X(I,3),ZHIGH) 
      ENDDO

!     Calculate contributions to the total energy of the system: 
       EKIN=0.0d0
       ETOT=0.0d0
       MTOT=0.0d0
       EMAG=0.0d0
       if ( igrav .eqv. .true. ) then
       EPOT=0.0d0
       endif
 
       do i=1,n      
       V2I=wprim(i,2)**2+wprim(i,3)**2+wprim(i,4)**2
       EKIN=EKIN+0.5d0*am(i)*V2I
       ETOT=ETOT+u(i,5)
       if ( igrav .eqv. .true. ) then
       ETOT=ETOT+am(i)*phi(i)/2.0d0
       endif
       mtot=mtot+u(i,1)
       EMAG=EMAG+vol(i)*(wprim(i,6)**2+wprim(i,7)**2+wprim(i,8)**2)/2.0d0
       if ( igrav .eqv. .true. ) then
       EPOT=EPOT+am(i)*phi(i)
       endif
       ENDDO
       
       if ( igrav .eqv. .true. ) then
       EPOT=0.5d0*EPOT
       endif

! kinetic energy  

         EINT=0.0d0
           DO i=1,n
           if ( imhd ) then
           EINT=EINT+u(i,5)-0.5d0*vol(i)*(wprim(i,1)*(wprim(i,2)**2+ &
           wprim(i,3)**2+wprim(i,4)**2)+wprim(i,6)**2+wprim(i,7)**2+ &
          wprim(i,8)**2)
          else
          EINT=EINT+u(i,5)-0.5d0*vol(i)*wprim(i,1)*(wprim(i,2)**2+ &
           wprim(i,3)**2+wprim(i,4)**2)
          endif
       ENDDO
      
          amvec(:)=0.0d0
          pvec(:)=0.0d0
       
!   calculate total angular momemtum of the system

      DO I=1,n
          AMVEC(1) = AMVEC(1) + VOL(I) *  &
         (X(I,2)*WPRIM(I,4) - X(I,3)*WPRIM(I,3))
          AMVEC(2) = AMVEC(2) + VOL(I) *  &
         (X(I,3)*WPRIM(I,2) - X(I,1)*WPRIM(I,4))
          AMVEC(3) = AMVEC(3) + VOL(I) *  &
         (X(I,1)*WPRIM(I,3) - X(I,2)*WPRIM(I,2))
          PVEC(1)=PVEC(1)+VOL(I)*WPRIM(I,2)
          PVEC(2)=PVEC(2)+VOL(I)*WPRIM(I,3)
          PVEC(3)=PVEC(3)+VOL(I)*WPRIM(I,4)
      ENDDO

!     Angular momentum contributions from SPH and sink particles 

      JTOT=DSQRT(AMVEC(1)**2+AMVEC(2)**2+AMVEC(3)**2)
      PTOT=DSQRT(PVEC(1)**2+PVEC(2)**2+PVEC(3)**2)
     
!     Get min/max values of various quantities:
       RHOMIN=1.E30
       RHOMAX=0.0d0
       irhomin=1
       irhomax=1
       HPMAX=0.0d0
       HPMIN=1.E30
       ihpmin=1
       ihpmax=1
       nnmax=0
       nnmin=1e+7
       innmin=1
       innmax=1 
           
          if ( imhd ) then

          ebmax=0.0d0
          ebmin=1.0d25
          ebave=0.0d0
          iebmax=1
          iebmin=1

          endif    
      
         do i=1,n 
  
            if ( imhd ) then
            eb(i)=dabs(divb(i))*h(i)/dsqrt(wprim(i,6)**2+wprim(i,7)**2+ &
            wprim(i,8)**2+tiny)
            endif

            IF (wprim(i,1).LT.RHOMIN) THEN
              IRHOMIN=I
              RHOMIN=wprim(i,1)
            ENDIF
            IF (wprim(i,1).GT.RHOMAX) THEN
              IRHOMAX=I
               RHOMAX=wprim(i,1)
            ENDIF
            HPI=H(I)
            IF (HPI.LT.HPMIN) THEN
               IHPMIN=I
               HPMIN=HPI
            ENDIF
            IF (HPI.GT.HPMAX) THEN
               IHPMAX=I
               HPMAX=HPI
            ENDIF

            if ( imhd ) then

            IF (eb(i).LT.EBMIN) THEN
               IEBMIN=I
               EBMIN=eb(i)
            ENDIF
            IF (eb(i).GT.EBMAX) THEN
               IEBMAX=I
               EBMAX=eb(i)
            ENDIF

            ebave=ebave+eb(i)

            endif

            IF (nne(i).GT.nnmax) THEN
               innmax=i
               nnmax=nne(i)
            ENDIF
            IF (nne(i).lt.nnmin) THEN
               innmin=i
               nnmin=nne(i)
            ENDIF
      ENDDO

      ebave=ebave/float(n)
      
      WRITE (2,*) 'OUTPUT: end of iteration ',NIT,' time= ',T
      WRITE (2,*) 'System box:'
      WRITE (2,*) XLOW,XHIGH,YLOW,YHIGH,ZLOW,ZHIGH
      WRITE (2,*) 'rhomin:',RHOMIN,X(IRHOMIN,1),X(IRHOMIN,2),X(IRHOMIN,3)
      WRITE (2,*) 'rhomax:',RHOMAX,X(IRHOMAX,1),X(IRHOMAX,2),X(IRHOMAX,3)
      WRITE (2,*) 'hmin:',HPMIN,X(IHPMIN,1),X(IHPMIN,2),X(IHPMIN,3)
      WRITE (2,*) 'hmax:',HPMAX,X(IHPMAX,1),X(IHPMAX,2),X(IHPMAX,3)
      WRITE (2,*) 'nnmin:',nnmin,X(INNMIN,1),X(INNMIN,2),X(INNMIN,3)
      WRITE (2,*) 'nnmax:',nnmax,X(INNMAX,1),X(INNMAX,2),X(INNMAX,3)
      if ( imhd ) then
      WRITE (2,*) 'ebmin:',EBMIN,X(IEBMIN,1),X(IEBMIN,2),X(IEBMIN,3)
      WRITE (2,*) 'ebmax:',EBMAX,X(IEBMAX,1),X(IEBMAX,2),X(IEBMAX,3)
 
      endif
         
        write(3,*) t,eint
        write(4,*) t,ekin
        write(8,*) t,etot
        write(9,*) t,jtot
        write(10,*) t,ptot 
        write(13,*) t,rhomax
        write(18,*) t,mtot
        if ( imhd .eqv. .true. ) then
        write(15,*) t,ebmin,ebmax,ebave
        write(7,*) t,emag
        endif
        if ( igrav .eqv. .true. ) then
        write(14,*) t,epot
        endif
                
      RETURN
      END
           

      SUBROUTINE TSTEP(wprim,dt)
          
 ! calculates the system time step
          
      use wp3d_h
          
      IMPLICIT NONE
   
      double precision::dt,dti(nmax),tiny
      PARAMETER (TINY=1.E-10)

      integer::i,j,idx
      double precision::cmax,cfl,Li,dx,dy,dz,valfveni,vmaxmhdi, &
      bij,eps2,wprim(nmax,nv),dti1,atot,dtkin

      eps2=1.0d-10
          
      cfl=0.2d0
    
      dt=1000.0d0
          
       if ( iDedner .eqv. .true. ) then
       ch=0.0d0
       endif

         do i=1,n
         
              if ( imhd .eqv. .true. ) then
                  
              valfveni=dsqrt((wprim(i,6)**2+wprim(i,7)**2+ &
              wprim(i,8)**2)/wprim(i,1))

              endif

              cmax=0.0d0

              do j=1,nne(i)

              idx=nni(i,j)

              if ( imhd .eqv. .true. ) then

              dx=x(idx,1)-x(i,1)
              dy=x(idx,2)-x(i,2)
              dz=x(idx,3)-x(i,3)
                  
              if ( iperiodic .eqv. .true. ) then
                  
              call modbound(dx,dy,dz)
                  
              endif

              if ( i .ne. idx ) then
                  
              bij=(wprim(i,6)*dx+wprim(i,7)*dy+wprim(i,8)*dz) &
              /dsqrt(dx**2+dy**2+dz**2)

              else

              bij=0.0d0

              endif
                  
              vmaxmhdi=dsqrt((cs(i)**2+valfveni**2)+ &
              dsqrt((cs(i)**2+valfveni**2)**2- &
              4*cs(i)**2*bij**2/wprim(i,1)))/dsqrt(2.0d0)
           
              cmax=dmax1(cmax,vmaxmhdi)

              else
          
              cmax=dmax1(cmax,cs(idx))

              endif

          enddo

          if ( iDedner .eqv. .true. ) then
          ch=dmax1(ch,cmax)
          endif

          Li=(3*vol(i)/pi/4)**(1.0d0/3.0d0)
                          
          dti(i)=Li/cmax

          if ( igrav .eqv. .true. ) then

          atot=dsqrt(sum(Vudot(i,2:4)**2))/vol(i)

          dtkin=dsqrt(Li/atot)

          dti(i)=dmin1(dti(i),dtkin)

          endif

          enddo

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
          cp=dsqrt(0.18d0*ch)

          if ( dt .lt. mdt ) then
          write(2,*) 'Integration has been stopped because of a too small timestep'

          select case (outputtype)  
          case (1)
!         splash format 
          call PDUMP(wprim)
          case (2)
!         silo format
          call silo_3d(wprim)
          case (3)
!         vtr format 
          call plot(wprim)
          end select 

          stop
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
          
          if ( iperiodic ) then
          
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
          
          if ( iperiodic ) then
          
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
      
          
      SUBROUTINE ADVANCE(u,wprim,psi,dt1)
          
      use wp3d_h
          
      IMPLICIT NONE
          
! **************************************************************
!    Advance conservative variables by one time step 
! **************************************************************
                                     
       integer::i,in,j,k
       double precision::dt1,dth1, &
       periodic1,periodic2,resx,resy,resz
       double precision::u(nmax,nv),u1(nmax,nv), &
       u2(nmax,nv),wprim(nmax,nv),wprim1(nmax,nv), &
       wprim2(nmax,nv),psi(nmax),psi1(nmax),psi2(nmax)
       external periodic1,periodic2
         
! Drift-Kick-Drift scheme for advancing the particles

        call TSTEP(wprim,dt1)

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

        CALL neighbourlists(wprim)

        if ( igrav .eqv. .true. ) then
        call MKTREE(wprim)
        endif

        endif

! neighbour lists and geometric quantities

        if ( iLagrangian .eqv. .true. ) then

        CALL RMD_coefficients

        else

        if ((  ilagrangian .eqv. .false. ) .and. &
        (nit .eq. 1 )) then

        CALL RMD_coefficients

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
        if ( iDedner) then
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
       CALL neighbourlists(wprim)
 
       if ( igrav .eqv. .true. ) then
       call MKTREE(wprim)
       endif
      
      endif
          
      write(2,*) 'step= ',NIT,' t=',t,' dt=',dt1         
      write(6,*) 'step= ',NIT,' t=',t,' dt=',dt1 
            
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
      
      
      SUBROUTINE RMD_coefficients
      
! computes coefficients for renormalized meshless derivatives
! ( see Lanson,Vila,Gaburov)  

      use wp3d_h

      integer::i,j,k,l,in,itab
      double precision::u2,hpi,h2i,h5i,d2,A(nmax,ndim,ndim)
      double precision::B(ndim,ndim),r(ndim),t1
      double precision::dwij(nmax),rij(nmax,ndim)
          
      do i=1,n
      
      A(i,:,:)=0.0d0
                
      hpi=h(i)
      h2i=hpi**2
      h5i=h2i*h2i*hpi
      
!     Calculate gradWij's to all neighbours:

      do in=1,nne(i) 

          j=nni(i,in)

          if ( i .ne. j ) then
           
          rij(in,1)=x(i,1)-x(j,1)
          rij(in,2)=x(i,2)-x(j,2)
          rij(in,3)=x(i,3)-x(j,3)
          
          if ( iperiodic ) then
          call modbound(rij(in,1),rij(in,2),rij(in,3))
          endif
            
          d2=rij(in,1)**2+rij(in,2)**2+rij(in,3)**2
                  
          U2=D2/H2I

          IF ( U2 .ge. 4.0d0 ) THEN
          DWIJ(in)=0.0d0
          ELSE
          ITAB=INT(CTAB*U2)+1
          DWIJ(IN)=DWTAB(ITAB)/H5I
          endif
                  
         t1=-rij(in,1)**2*dwij(in)    
         A(i,1,1)=A(i,1,1)+vol(j)*t1
         t1=-rij(in,1)*rij(in,2)*dwij(in)
         A(i,1,2)=A(i,1,2)+vol(j)*t1
         t1=-rij(in,1)*rij(in,3)*dwij(in)
         A(i,1,3)=A(i,1,3)+vol(j)*t1
         t1=-rij(in,2)**2*dwij(in)
         A(i,2,2)=A(i,2,2)+vol(j)*t1
         t1=-rij(in,1)*rij(in,3)*dwij(in)
         A(i,1,3)=A(i,1,3)+vol(j)*t1
         t1=-rij(in,2)*rij(in,3)*dwij(in)
         A(i,2,3)=A(i,2,3)+vol(j)*t1
         t1=-rij(in,3)**2*dwij(in)
         A(i,3,3)=A(i,3,3)+vol(j)*t1
                 
         endif
       
         enddo
                 
         enddo
                 
         A(:,2,1)=A(:,1,2)
         A(:,3,1)=A(:,1,3)
         A(:,3,2)=A(:,2,3)
                 
         do i=1,n
                 
         r(:)=0.0d0

         call gaussj(A(i,:,:),ndim,ndim,r,1,ndim)
                 
         B=A(i,:,:)
                 
         do in=1,nne(i)
         
         r(1)=rij(in,1)*dwij(in)
         r(2)=rij(in,2)*dwij(in)
         r(3)=rij(in,3)*dwij(in)
          
         do k=1,ndim
         C(i,in,k)=0.0d0
         enddo
         
         do l=1,ndim
         do k=1,ndim
                 
         C(i,in,l)=C(i,in,l)+B(l,k)*r(k)
                 
         enddo
         enddo
         
         C(i,in,:)=C(i,in,:)*vol(j)

         enddo
                
         enddo
              
       RETURN
      
       end 

           
      SUBROUTINE gaussj(a,n,np,b,m,mp) 

      INTEGER::m,mp,n,np,NMAX 
      double precision a(np,np),b(np,mp) 
      PARAMETER (NMAX=50) 
      INTEGER::i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX), &
     ipiv(NMAX) 
      double precision big,dum,pivinv 

      do 11 j=1,n 
      ipiv(j)=0 
 11   continue 
      do 22 i=1,n
      big=0.0d0
      do 13 j=1,n 
      if( ipiv(j).ne.1) then 
      do 12 k=1,n 
      if (ipiv(k).eq.0) then
      if (abs(a(j,k)).ge.big)then 
      big=abs(a(j,k)) 
      irow=j 
      icol=k 
      endif 
      else if  (ipiv(k).gt.1)  then 
      stop
      endif 
 12   continue 
      endif 
 13   continue 
      ipiv(icol)=ipiv(icol)+1 
      if (irow.ne.icol) then 
      do 14 l=1,n 
      dum=a(irow,l) 
      a(irow,l)=a(icol,l) 
      a(icol,l)=dum 
 14   continue  
      do 15 l=1,m 
      dum=b(irow,l) 
      b(irow,l)=b(icol,l) 
      b(icol,l)=dum 
 15   continue  
      endif 
      indxr(i)=irow
      indxc(i)=icol 
      if (a(icol,icol) .eq. 0.0d0 ) then
      stop
      endif
      pivinv=1.0d0/a(icol,icol) 
      a(icol,icol)=1.0d0
      do 16 l=1,n 
      a(icol,l)=a(icol,l)*pivinv 
 16   continue 
      do  17 l=1,m 
      b(icol,l)=b(icol,l)*pivinv 
 17   continue  
      do  21 ll=1,n 
      if(ll.ne.icol)then 
      dum=a(ll,icol) 
      a(ll,icol)=0. 
      do  18 l=1,n 
      a(ll,l)=a(ll,l)-a(icol,l)*dum 
 18   continue 
      do 19 l=1,m 
      b(ll,l)=b(ll,l)-b(icol,l)*dum 
 19   continue  
      endif 
 21   continue 
 22   continue
      do 24 l=n,1,-1 
      if(indxr(l).ne.indxc(l))then 
      do 23 k=1,n 
      dum=a(k,indxr(l))
      a(k,indxr(l))=a(k,indxc(l)) 
      a(k,indxc(l))=dum 
 23   continue 
      endif 
 24   continue 
     
      return 
      END 

      
      SUBROUTINE CALCDOTS(u,wprim,psi)
          
      use wp3d_h
          
      IMPLICIT NONE
          
! *****************************************************
!     This routine calculates d(Vu)/dt in the 
!     discrete conservation law(see Vila,Lanson)
! *****************************************************
          
      double precision::uleft(nv),uright(nv), &
      uleftr(nv),urightr(nv),normnij,F(nv),dx,dy,dz,d, &
      nij(ndim),a(ndim),ax,u(nmax,nv),S(nv),gradpsi(n,ndim), &
      bxc,psic,wprim(nmax,nv),psi(nmax),t1,u2,&
      eps,dgradw(ndim),gaccxi,gaccyi,gacczi,dzetai, &
      h5i,h5,dwij1,dwij2,phii,psil,psir,help(ndim)

      integer::i,j,k,in,jn,itab
       
! for the update of the accelerations, we use the interacting neighbour list

! for particle i

         eps=1.0d-3
          
         do i=1,n
         
         Vudot(i,:)=0.0d0
         divb(i)=0.0d0
         gradpsi(i,:)=0.0d0
         if ( igrav .eqv. .true. ) then
	     help(:)=0.0d0
         endif
           
         do in=1,nne(i)

         j=nni(i,in)
           
          if ( i .ne. j ) then
                  
          jn=0
           
          do k=1,nne(j)
           
          if ( nni(j,k) .eq. i ) then
           
          jn=k
          
          endif
          
          enddo

          if ( jn .ne. 0 ) then
                  
          dgradw(1)=vol(i)*C(i,in,1)-vol(j)*C(j,jn,1)
          dgradw(2)=vol(i)*C(i,in,2)-vol(j)*C(j,jn,2)
          dgradw(3)=vol(i)*C(i,in,3)-vol(j)*C(j,jn,3)

          if ( iDedner .eqv. .true. ) then
          psil=psi(i)
          psir=psi(j)
          endif
           
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

          a(1)=(wprim(i,2)+wprim(j,2))/2
          a(2)=(wprim(i,3)+wprim(j,3))/2
          a(3)=(wprim(i,4)+wprim(j,4))/2
          ax=nij(1)*a(1)+nij(2)*a(2)+nij(3)*a(3)               
          else                  
          ax=0.0d0
          endif
                    
          if ( imhd .eqv. .true. ) then
                  
          select case ( solverchoice )

          case (1)

          ! HLL MHD Riemann solver
                  
           call HLLMHD(i,j,uleftr,urightr,nij,F,ax,bxc,psi,psic)
                                   
           case (2)
		   
	   ! HLLD MHD Riemann solver

           call HLLDMHD(i,j,uleftr,urightr,nij,F,ax,bxc,psi,psic)

           end select 
           
           select case (divbtype) 
                  
          case(1)
                  
          divb(i)=divb(i)+(wprim(j,6)-wprim(i,6))*C(i,in,1)+ &
          (wprim(j,7)-wprim(i,7))*C(i,in,2)+ &
          (wprim(j,8)-wprim(i,8))*C(i,in,3)
                  
          case (2)
                  
          divb(i)=divb(i)+bxc*normnij
                  
          end select 

                       
           else
                        
           select case (solverchoice ) 
                        
           case (2)
                        
           call HLLCHD(i,j,uleftr,urightr,nij,F,ax)
                                
           end select

          endif

! solve the 1D Riemann problem for the particle pair i,j in 
! the moving reference frame and rotate back to the lab frame

! update the vector of conservative variables 

          do k=1,nv
               
          t1=normnij*F(k)
          Vudot(i,k)=Vudot(i,k)+t1
		  
          if ( igrav .and. ( .not. barotropic )) then
          if ( .not. exterior(i) ) then
          if ( k .eq. 1 ) then
          help(1)=help(1)-dx*t1
          help(2)=help(2)-dy*t1
          help(3)=help(3)-dz*t1
          endif
          endif
          endif
		  
          enddo
		  
          if ( igrav .eqv. .true. ) then
          if ( .not. exterior(i) ) then
		  
	      h5i=h(i)**5
		  
          u2=sum(nij(:)**2)
          if ( U2 .ge. 4.0d0 ) THEN
          DWIJ1=0.0d0
          ELSE
          ITAB=INT(CTAB*U2)+1
          DWIJ1=DWTAB(ITAB)/H5I
          ENDIF 
		  
	      H5=h(j)**5
		  
	      u2=(dx**2+dy**2+dz**2)/h(j)**2
          if ( U2 .ge. 4.0d0 ) THEN
          DWIJ2=0.0d0
          ELSE
          ITAB=INT(CTAB*U2)+1
          DWIJ2=DWTAB(ITAB)/H5
          ENDIF 
		  
          Vudot(i,2)=Vudot(i,2)-G*(dzeta(i)*dwij1+dzeta(j)*dwij2)*dx/2
	      Vudot(i,3)=Vudot(i,3)-G*(dzeta(i)*dwij1+dzeta(j)*dwij2)*dy/2
	      Vudot(i,4)=Vudot(i,4)-G*(dzeta(i)*dwij1+dzeta(j)*dwij2)*dz/2
		  
! terms related to variable smoothing lengths and energy conservation 
! in self-gravitating flows

          endif          
          endif
		  
		  
          if ( iDedner ) then 
                  
          select case (gradtype) 
          case (1)
          gradpsi(i,1)=gradpsi(i,1)+(psi(j)-psi(i))*C(i,in,1)
          gradpsi(i,2)=gradpsi(i,2)+(psi(j)-psi(i))*C(i,in,2)
          gradpsi(i,3)=gradpsi(i,3)+(psi(j)-psi(i))*C(i,in,3)
          case (2)                
          gradpsi(i,1)=gradpsi(i,1)+psic*dgradw(1)
          gradpsi(i,2)=gradpsi(i,2)+psic*dgradw(2)
          gradpsi(i,3)=gradpsi(i,3)+psic*dgradw(3)
          end select
                  
          endif

         endif

         endif

         enddo
         
         

         enddo
                 
         do i=1,n
                
         S(:)=0.0d0
                 
         if ( imhd ) then

 ! include Powell's source terms  
 ! and Dedner's correction terms

         if ( iDedner ) then
                 
         if ( divbtype .eq. 2) then        
         divb(i)=divb(i)/vol(i)
         endif
         if ( gradtype  .eq. 2) then
         gradpsi(i,:)=gradpsi(i,:)/vol(i)
         endif
         psidot(i)=-ch**2*divb(i)
  ! Hyperbolic term
         psidot(i)=psidot(i)-(ch/cp)**2*psi(i)
 !  Parabolic term
 
         S(5)=S(5)-(wprim(i,6)*gradpsi(i,1)+ &
         wprim(i,7)*gradpsi(i,2)+wprim(i,8)*gradpsi(i,3))
         S(6)=S(6)-gradpsi(i,1)
         S(7)=S(7)-gradpsi(i,2)
         S(8)=S(8)-gradpsi(i,3)

         endif

         if ( iPowell ) then
 
         S(2)=S(2)-divb(i)*wprim(i,6)
         S(3)=S(3)-divb(i)*wprim(i,7)
         S(4)=S(4)-divb(i)*wprim(i,8)
         if ( .not. barotropic ) then
         S(5)=S(5)-divb(i)*(wprim(i,2)*wprim(i,6)+ &
         wprim(i,3)*wprim(i,7)+wprim(i,4)*wprim(i,8))
         endif
         S(6)=S(6)-divb(i)*wprim(i,2)
         S(7)=S(7)-divb(i)*wprim(i,3)
         S(8)=S(8)-divb(i)*wprim(i,4)

         endif

        if ( igrav .eqv. .true. ) then

        if ( exterior(i) .eqv. .true. ) then
        gaccxi=0.0d0
        gaccyi=0.0d0
        gacczi=0.0d0
        else
        call TRG(i,gaccxi,gaccyi,gacczi,phii)
        endif
		
	    phi(i)=phii

        S(2)=S(2)+G*wprim(i,1)*gaccxi
        S(3)=S(3)+G*wprim(i,1)*gaccyi
        S(4)=S(4)+G*wprim(i,1)*gacczi
        if ( .not. barotropic ) then
        S(5)=S(5)+G*(wprim(i,1)*(wprim(i,2)*gaccxi+ &
	    wprim(i,3)*gaccyi+wprim(i,4)*gacczi)+ &
	   (gaccxi*help(1)+gaccyi*help(2)+gacczi*help(3))/vol(i))
        endif

! source terms originating from self-gravity

        endif

        Vudot(i,1:4)=Vudot(i,1:4)-vol(i)*S(1:4)
        if ( .not. barotropic ) then
        Vudot(i,5)=Vudot(i,5)-vol(i)*S(5)
        endif
        if ( imhd ) then
        Vudot(i,6:8)=Vudot(i,6:8)-vol(i)*S(6:8)
        endif     
        endif
                
        enddo
                                
      RETURN
      END
          
          subroutine find_gradients_w(wprim)
          
          USE wp3d_h
          
          implicit none
                  
          double precision:: wprim(nmax,nv), &
         dx,dy,dz,d
          
! this subroutine computes the gradients of the primitive 
! variables for each particle 
          
          integer:: i,k,l,m,idx
                  
          do i=1,n
                     
          do k=1,nv
          
          gradw(i,k,:)=0.0d0
         
          do l=1,ndim
          
          do m=1,nne(i)
          
          idx=nni(i,m)
                  
          dx=x(i,1)-x(idx,1)
          dy=x(i,2)-x(idx,2)
          dz=x(i,3)-x(idx,3)
                  
          if ( iperiodic ) then
          call modbound(dx,dy,dz)
          endif
                  
          d=dsqrt(dx**2+dy**2+dz**2)
                  
          if (( i .ne. idx) .and. ( d .le. 2*h(i))) then
          
          gradw(i,k,l)=gradw(i,k,l)+(wprim(idx,k)-wprim(i,k))*C(i,m,l)

          endif
          
          enddo
         
          enddo
                  
          enddo
               
          enddo

          return

         end

  
        subroutine find_gradients_u(u)
          
          USE wp3d_h
          
          implicit none

          integer::i,m,k,l,idx 
          double precision::u(nmax,nv), &
          dx,dy,dz,d
                  
          do i=1,n
                
          do k=1,nv
          
          gradu(i,k,:)=0.0d0
         
          do l=1,ndim
          
          do m=1,nne(i)
          
          idx=nni(i,m)
                  
          dx=x(i,1)-x(idx,1)
          dy=x(i,2)-x(idx,2)
          dz=x(i,3)-x(idx,3)
                  
          if ( iperiodic ) then
          call modbound(dx,dy,dz)
          endif
                  
          d=dsqrt(dx**2+dy**2+dz**2)

          if (( i .ne. idx ) .and. ( d .le. 2*h(i))) then
          
          gradu(i,k,l)=gradu(i,k,l)+(u(idx,k)/vol(idx)- &
                  u(i,k)/vol(idx))*C(i,m,l)

          endif
          
          enddo
         
          enddo
          
          enddo
                   
          enddo
         
          return
          
          end


          subroutine find_limiters_w(wprim)

          USE wp3d_h

          implicit none

          integer::i,k,l,idx
                  
          double precision::wprim(nmax,nv)
          double precision:: dx,dy,dz,tiny, &
          kappa(nv),dx2,dy2,dz2,maxwreci(nv),maxwi(nv), &
          minwreci(nv),minwi(nv),wrec(nv),dwij(nv),psiij(nv)

          kappa(1)=1.0d0
          kappa(2)=1.0d0
          kappa(3)=1.0d0
          kappa(4)=1.0d0
          kappa(5)=1.0d0
          kappa(6)=1.0d0
          kappa(7)=1.0d0
          kappa(8)=1.0d0
 
          tiny=1.0d-10
                                 
          select case (limtype) 
                  
          case (1) 

          do i=1,n
              
          lim(i,:)=1.0d25
                                          
          maxwreci(:)=0.0d0
          minwreci(:)=1.0d+25
          maxwi(:)=0.0d0
          minwi(:)=1d+25
          
          do k=1,nne(i)
          
          idx=nni(i,k)

          dx=x(idx,1)-x(i,1)
          dy=x(idx,2)-x(i,2)
          dz=x(idx,3)-x(i,3)

          if ( iperiodic ) then
          
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
     
          lim(i,k)=dmin1(lim(i,k),1.0d0,kappa(k)*(maxwi(k)-wprim(i,k))/ &
          (maxwreci(k)-wprim(i,k)+tiny),kappa(k)*(wprim(i,k)-minwi(k))/ &
          (wprim(i,k)-minwreci(k)+tiny))

          enddo

          enddo

          case (2)
         
          do i=1,n
                 
          lim(i,:)=1.0d25

          maxwi(:)=0.0d0
          minwi(:)=1.0d+25
        
          do k=1,nne(i)
          
          idx=nni(i,k)
         
          do l=1,nv
          
          maxwi(l)=dmax1(maxwi(l),wprim(idx,l))
          minwi(l)=dmin1(minwi(l),wprim(idx,l))
          
          enddo
                  
           enddo
                  
          do k=1,nne(i)
                  
          idx=nni(i,k)

          if ( i .ne. idx ) then

          dx=x(idx,1)-x(i,1)
          dy=x(idx,2)-x(i,2)
          dz=x(idx,3)-x(i,3)

          if ( iperiodic ) then
          
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

             lim(i,:)=dmin1(lim(i,:),psiij(:),1.0d0)

               endif

                enddo
                                 
                 enddo

                end select

                return

               end

              subroutine find_limiters_u(u)

              USE wp3d_h

              implicit none

              integer::i,k,l,idx
                  
            double precision::u(nmax,nv)
            double precision::  dx,dy,dz,tiny, & 
          maxureci(nv),minureci(nv),maxui(nv), &
          minui(nv),kappa(nv),dx2,dy2,dz2,urec(nv),duij(nv), &
          psiij(nv)

            kappa(:)=1.0d0
            tiny=1.0d-10

            select case (limtype ) 

              case (1)

              do i=1,n
                      
           lim(i,:)=1.0d25

          maxureci(:)=0.0d0
          minureci(:)=1.0d+25
          maxui(:)=0.0d0
          minui(:)=1d+25
          
          do k=1,nne(i)
          
          idx=nni(i,k)
          
          dx=x(idx,1)-x(i,1)
          dy=x(idx,2)-x(i,2)
          dz=x(idx,3)-x(i,3)

          if ( iperiodic ) then
          
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

          lim(i,k)=dmin1(lim(i,k),1.0d0,kappa(k)*(maxui(k)-u(i,k))/ &
          (maxureci(k)-u(i,k)+tiny),kappa(k)*(u(i,k)-minui(k))/ &
          (u(i,k)-minureci(k)+tiny))

          enddo

          enddo
                  
          case (2)
                  
            do i=1,n
               
              lim(i,:)=1.0d25
                             
             maxui(:)=0.0d0
            minui(:)=1d+25
          
            do k=1,nne(i)
          
            idx=nni(i,k)
          
            do l=1,nv
          
            maxui(l)=dmax1(maxui(l),u(idx,l)/vol(idx))
            minui(l)=dmin1(minui(l),u(idx,l)/vol(idx))
          
            enddo
                  
            enddo
          
            do k=1,nne(i)

            idx=nni(i,k)

            if ( i .ne. idx ) then
                  
            dx=x(idx,1)-x(i,1)
            dy=x(idx,2)-x(i,2)
            dz=x(idx,3)-x(i,3)

            if ( iperiodic ) then
          
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

            endif

            lim(i,:)=dmin1(lim(i,:),psiij(:),1.0d0)
                  
            enddo
                     
            enddo
                  
            end select 

          return

          end

          
          subroutine reconstruct(i,j,u,wprim,uleft,uright)
          
! this subroutine computes the left and right states at the 
! midpoint between particles i and j
          
          USE wp3d_h
          
          implicit none
                 
          integer:: i,j,k
          
          double precision:: u(nmax,nv),wprim(nmax,nv), &
          uleft(nv),uright(nv),gradleft(nv),gradright(nv)
          double precision:: wleft(nv),wright(nv)
          double precision:: dx,dy,dz,dx2,dy2,dz2

          dx=x(j,1)-x(i,1)
          dy=x(j,2)-x(i,2)
          dz=x(j,3)-x(i,3)
                  
          if ( iperiodic ) then
                  
          call modbound(dx,dy,dz)
                  
          endif
                  
          dx2=dx/2
          dy2=dy/2
          dz2=dz/2
                                  
          if ( useprimitive ) then
                  
          do k=1,nv
                  
          gradleft(k)=gradw(i,k,1)*dx2+ &
          gradw(i,k,2)*dy2+gradw(i,k,3)*dz2
          
          gradright(k)=gradw(j,k,1)*dx2+ &
          gradw(j,k,2)*dy2+gradw(j,k,3)*dz2
                  
          if ( gradleft(k)*gradright(k) .le. 0.0d0 ) then
          wleft(k)=wprim(i,k)
          wright(k)=wprim(j,k)
          else
          wleft(k)=wprim(i,k)+lim(i,k)*gradleft(k)
          wright(k)=wprim(j,k)-lim(j,k)*gradright(k)
          endif
                  
          enddo
                              
          call primitive_to_conservative(i,wleft,uleft)
          call primitive_to_conservative(j,wright,uright)
                  
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
                 
                  if ( imhd ) then
                                  
                  if ( .not. barotropic ) then

                  ui(5)=wi(5)/(gam(i)-1.0d0)+0.5d0*wi(1)*( &
                  wi(2)**2+wi(3)**2+wi(4)**2)+0.5d0*(wi(6)**2+ &
                  wi(7)**2+wi(8)**2)
                                  
                  endif
                  
                  ui(6)=wi(6)
                  ui(7)=wi(7)
                  ui(8)=wi(8)
                  
                  else

                  if ( .not. barotropic) then

                  ui(5)=wi(5)/(gam(i)-1.0d0)+0.5d0*wi(1)*( &
                  wi(2)**2+wi(3)**2+wi(4)**2)

                  endif
                       
                  endif
                  
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
                 
                  if ( imhd ) then

                  if ( .not. barotropic ) then
                  wi(5)=(gam(i)-1.0d0)*(ui(5)/vol(i)-0.5d0*(ui(2)**2+ &
                  ui(3)**2+ui(4)**2)/ui(1)/vol(i)-0.5d0*(ui(6)**2+ui(7)**2+ &
                  ui(8)**2)/vol(i)**2) 
                  else
                  wi(5)=pbar(i,wi(1))
                  endif
                                  
                  wi(6)=ui(6)/vol(i)
                  wi(7)=ui(7)/vol(i)
                  wi(8)=ui(8)/vol(i)
                 
                  else

                  if ( .not. barotropic ) then
                  wi(5)=(gam(i)-1.0d0)*(ui(5)/vol(i)-0.5d0*(ui(2)**2+ &
                  ui(3)**2+ui(4)**2)/ui(1)/vol(i))
                  else
                  wi(5)=pbar(i,wi(1))
                  endif
                  
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

          if ( imhd ) then

          vec(1)=ui(6)
          vec(2)=ui(7)
          vec(3)=ui(8)

          call rotate_vector(vec,vecr,nij,1)
          
          uir(6)=vecr(1)
          uir(7)=vecr(2)
          uir(8)=vecr(3)

          endif
         
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
          
              
         subroutine HLLCHD(i,j,uleft,uright,nij,F,ax)
        
         USE wp3d_h
         
         implicit none
         
! HLLC Riemann solver for the system of 1D Euler equations
! see Toro chapter 10
         
         double precision:: uleft(nv),uright(nv), &
         ax,F(nv),Fv(ndim),Fvr(ndim),nij(ndim)
         double precision::vxl,vyl,vzl,el,rhol,vxr,vyr,vzr,er,rhor, &
         pr,pl,cl,cr,cmax,Sl,Sr,Sm,vxls,vxrs,vyls,vyrs, &
         vzls,vzrs,rhols,rhors,els,ers,vmax,pbar
         integer::i,j

         rhol=uleft(1)
         vxl=uleft(2)/rhol
         vyl=uleft(3)/rhol
         vzl=uleft(4)/rhol
         el=uleft(5)
                 
         if ( barotropic ) then
         pl=pbar(i,rhol)
         cl=dsqrt(pl/rhol)  
         else
         pl=(gam(i)-1.0d0)*(el-0.5d0*rhol*(vxl**2+vyl**2+vzl**2))
         cl=dsqrt(gam(i)*pl/rhol)  
         endif
        
         rhor=uright(1)
         vxr=uright(2)/rhor
         vyr=uright(3)/rhor
         vzr=uright(4)/rhor
         er=uright(5)
                 
         if ( barotropic ) then
         pr=pbar(j,rhor)
         cr=dsqrt(pr/rhor)  
         else
         pr=(gam(j)-1.0d0)*(er-0.5d0*rhor*(vxr**2+vyr**2+vzr**2))
         cr=dsqrt(gam(j)*pr/rhor)
         endif
                 
         cmax=dmax1(cl,cr)
         vmax=dmax1(vxl,vxr)
         
         Sl=vmax-cmax
         Sr=vmax+cmax
                 
         Sm=(pr-pl+rhol*vxl*(Sl-vxl)-rhor*vxr*(Sr-vxr))/ &
         (rhol*(Sl-vxl)-rhor*(Sr-vxr))
                 
         rhols=rhol*(Sl-vxl)/(Sl-Sm)
         rhors=rhor*(Sr-vxr)/(Sr-Sm)
                 
         vxls=(Sl-vxl)*Sm*rhol/(Sl-Sm)/rhols
         vxrs=(Sr-vxr)*Sm*rhor/(Sr-Sm)/rhors
                 
         vyls=(Sl-vxl)*vyl*rhol/(Sl-Sm)/rhols
         vyrs=(Sr-vxr)*vyr*rhor/(Sr-Sm)/rhors
                 
         vzls=(Sl-vxl)*vzl*rhol/(Sl-Sm)/rhols
         vzrs=(Sr-vxr)*vzr*rhor/(Sr-Sm)/rhors
                 
         els=((Sl-vxl)*rhol/(Sl-Sm))*(el/rhol+(Sm-vxl)* &
                  (Sm+pl/(rhol*(Sl-vxl))))
         ers=((Sr-vxr)*rhor/(Sr-Sm))*(er/rhor+(Sm-vxr)* &
                  (Sm+pr/(rhor*(Sr-vxr))))
                 
         if ( Sl .ge. ax ) then
         
         F(1)=rhol*vxl-ax*rhol
         Fv(1)=rhol*vxl**2+pl-ax*rhol*vxl
         Fv(2)=rhol*vxl*vyl-ax*rhol*vyl
         Fv(3)=rhol*vxl*vzl-ax*rhol*vzl
         
         call rotate_vector(Fv,Fvr,nij,-1)
         
! the flux vector related to the velocity is rotated back to 
! the lab frame to update the velocity vectors in the lab frame
         
         F(2)=Fvr(1)
         F(3)=Fvr(2)
         F(4)=Fvr(3)
         
         F(5)=(el+pl)*vxl-ax*el
         
         endif
         
         if ( Sr .le. ax ) then
         
         F(1)=rhor*vxr-ax*rhor
         Fv(1)=rhor*vxr**2+pr-ax*rhor*vxr
         Fv(2)=rhor*vxr*vyr-ax*rhor*vyr
         Fv(3)=rhor*vxr*vzr-ax*rhor*vzr
         
         call rotate_vector(Fv,Fvr,nij,-1)
         
! the flux vector related to the velocity is rotated back to 
! the lab frame to update the velocity vectors in the lab frame
         
         F(2)=Fvr(1)
         F(3)=Fvr(2)
        F(4)=Fvr(3)
         
         F(5)=(er+pr)*vxr-ax*er
         
         endif
                        
         if ( (Sl .le. ax ) .and. ( Sm .ge. ax) ) then
         
         F(1)=rhol*vxl-ax*rhols+Sl*(rhols-rhol)
                 
         Fv(1)=rhol*vxl**2+pl-ax*rhols*vxls+Sl*(rhols*vxls-rhol*vxl)
         Fv(2)=rhol*vxl*vyl-ax*rhols*vyls+Sl*(rhols*vyls-rhol*vyl)
         Fv(3)=rhol*vxl*vzl-ax*rhols*vzls+Sl*(rhols*vzls-rhol*vzl)
                       
         call rotate_vector(Fv,Fvr,nij,-1)
         
         F(2)=Fvr(1)
         F(3)=Fvr(2)
         F(4)=Fvr(3)
                 
         F(5)=(el+pl)*vxl-ax*els+Sl*(els-el)
         
         endif
                      
           if ( (Sr .ge. ax ) .and. ( Sm .le. ax) ) then
         
         F(1)=rhor*vxr-ax*rhors+Sr*(rhors-rhor)
                 
         Fv(1)=rhor*vxr**2+pr-ax*rhors*vxrs+Sr*(rhors*vxrs-rhor*vxr)
         Fv(2)=rhor*vxr*vyr-ax*rhors*vyrs+Sr*(rhors*vyrs-rhor*vyr)
         Fv(3)=rhor*vxr*vzr-ax*rhors*vzrs+Sr*(rhors*vzrs-rhor*vzr)
                 
         call rotate_vector(Fv,Fvr,nij,-1)
         
         F(2)=Fvr(1)
         F(3)=Fvr(2)
         F(4)=Fvr(3)
                 
         F(5)=(er+pr)*vxr-ax*ers+Sr*(ers-er)
         
         endif
                 
         return
         
         end
                 
                 
         subroutine HLLMHD(i,j,uleft,uright,nij,F,ax,bxc,psi,psic)
                
         use wp3d_h

         implicit none

! HLL Riemann solver for the system of 1D MHD equations

         double precision:: uleft(nv),uright(nv), &
         uLL(nv),ax,F(nv),Fv(ndim),Fvr(ndim),nij(ndim), &
         Fl(nv),Fr(nv),al,ar,bl,br,valfvenxl,valfvenxr
         double precision::vxl,vyl,vzl,el,rhol,vxr,vyr,vzr,er,rhor, &
         pr,pl,cl,cr,Sl,Sr,bxl,bxr,byl,bzl,byr,bzr,plt,prt,cmax
         double precision:: bxc,psil,psir,psic,pbar,psi(nmax)
         integer::i,j

         rhol=uleft(1)
         vxl=uleft(2)/rhol
         vyl=uleft(3)/rhol
         vzl=uleft(4)/rhol
         el=uleft(5)
         psil=psi(i)
         bxl=uleft(6)
         byl=uleft(7)
         bzl=uleft(8)
                
         if ( barotropic ) then
         pl=pbar(i,rhol)
         al=dsqrt(pl/rhol)     
         else
         pl=(gam(i)-1.0d0)*(el-0.5d0*(rhol*(vxl**2+vyl**2+vzl**2)+& 
         (bxl**2+byl**2+bzl**2)))
         al=dsqrt(gam(i)*pl/rhol)     
         endif
        
         plt=pl+(bxl**2+byl**2+bzl**2)/2

         rhor=uright(1)
         vxr=uright(2)/rhor
         vyr=uright(3)/rhor
         vzr=uright(4)/rhor
		 
         psir=psi(j)
         er=uright(5)
         bxr=uright(6)
         byr=uright(7)
         bzr=uright(8)

         if ( barotropic ) then  
         pr=pbar(j,rhor)   
         else      
         pr=(gam(j)-1.0d0)*(er-0.5d0*(rhor*(vxr**2+vyr**2+vzr**2)+ &
        (bxr**2+byr**2+bzr**2)))
         ar=dsqrt(gam(j)*pr/rhor)
         endif

         prt=pr+(bxr**2+byr**2+bzr**2)/2

         bl=dsqrt((byl**2+bxl**2+bzl**2)/rhol)
         br=dsqrt((byr**2+bxr**2+bzr**2)/rhor)
         valfvenxl=dsqrt(bxl**2/rhol)
         valfvenxr=dsqrt(bxr**2/rhor)

         cl=dsqrt(al**2+bl**2+dsqrt((al**2+bl**2)**2-  &
           4*al**2*valfvenxl**2))/dsqrt(2.0d0)

         cr=dsqrt(ar**2+br**2+dsqrt((ar**2+br**2)**2-  &
         4*ar**2*valfvenxr**2))/dsqrt(2.0d0)

         cmax=dmax1(cl,cr)
         Sl=dmin1(vxl,vxr)-cmax
         Sr=dmax1(vxl,vxr)+cmax

 ! fast magnetosonic wave speeds for the left and right state
 
         if ( iDedner ) then
         bxc=(bxl+bxr)/2-(psir-psil)/(2*cmax)
         psic=(psil+psir)/2-cmax*(bxr-bxl)/2
         else
         bxc=(bxl+bxr)/2
         endif
                 
! solution of a separate Riemann problem involving
! the evolution equation of the scalar potential

         if ( Sl .ge. ax ) then

         F(1)=rhol*vxl-ax*rhol
         Fv(1)=rhol*vxl**2+plt-bxc**2-ax*rhol*vxl
         Fv(2)=rhol*vxl*vyl-bxc*byl-ax*rhol*vyl
         Fv(3)=rhol*vxl*vzl-bxc*bzl-ax*rhol*vzl

         call rotate_vector(Fv,Fvr,nij,-1)

! the flux vector related to the velocity is rotated back to
! the lab frame to update the velocity vectors in the lab frame

         F(2)=Fvr(1)
         F(3)=Fvr(2)
         F(4)=Fvr(3)

         F(5)=(el+plt)*vxl-(vxl*bxc+vyl*byl+bzl*vzl)*bxc-ax*el

         Fv(1)=-ax*bxc
         Fv(2)=byl*vxl-bxc*vyl-ax*byl
         Fv(3)=bzl*vxl-bxc*vzl-ax*bzl

         call rotate_vector(Fv,Fvr,nij,-1)

         F(6)=Fvr(1)
         F(7)=Fvr(2)
         F(8)=Fvr(3)

         endif

         if ( Sr .le. ax ) then

         F(1)=rhor*vxr-ax*rhor
         Fv(1)=rhor*vxr**2+prt-bxc**2-ax*rhor*vxr
         Fv(2)=rhor*vxr*vyr-bxc*byr-ax*rhor*vyr
         Fv(3)=rhor*vxr*vzr-bxc*bzr-ax*rhor*vzr

         call rotate_vector(Fv,Fvr,nij,-1)

! the flux vector related to the velocity is rotated back to
! the lab frame to update the velocity vectors in the lab frame

         F(2)=Fvr(1)
         F(3)=Fvr(2)
         F(4)=Fvr(3)

         F(5)=(er+prt)*vxr-(vxr*bxc+vyr*byr+vzr*bzr)*bxc-ax*er
                 
         Fv(1)=-ax*bxc
         Fv(2)=byr*vxr-bxc*vyr-ax*byr
         Fv(3)=bzr*vxr-bxc*vzr-ax*bzr
                 
         call rotate_vector(Fv,Fvr,nij,-1)

         F(6)=Fvr(1)
         F(7)=Fvr(2)
         F(8)=Fvr(3)
                 
         endif

         if ( (Sl .le. ax ) .and. ( Sr .ge. ax) ) then

         Fl(1)=rhol*vxl-ax*rhol
         Fl(2)=rhol*vxl**2+plt-bxc**2-ax*rhol*vxl
         Fl(3)=rhol*vxl*vyl-bxc*byl-ax*rhol*vyl
         Fl(4)=rhol*vxl*vzl-bxc*bzl-ax*rhol*vzl
         Fl(5)=(el+plt)*vxl-(vxl*bxc+vyl*byl+vzl*bzl)*bxc-ax*el
         Fl(5)=-ax*bxc
         Fl(7)=byl*vxl-bxc*vyl-ax*byl
         Fl(8)=bzl*vxl-bxc*vzl-ax*bzl

         Fr(1)=rhor*vxr-ax*rhor
         Fr(2)=rhor*vxr**2+prt-bxc**2-ax*rhor*vxr
         Fr(2)=rhor*vxr*vyr-bxc*byr-ax*rhor*vyr
         Fr(3)=rhor*vxr*vzr-bxc*bzr-ax*rhor*vzr
         Fr(5)=(er+prt)*vxr-(vxr*bxc+vyr*byr+bzr*vzr)*bxc-ax*er
         Fr(6)=-ax*bxc
         Fr(7)=byr*vxr-bxc*vyr-ax*byr
         Fr(8)=bzr*vxr-bxc*vzr-ax*bzr
                 
         uLL(:)=(Sr*uright(:)-Sl*uleft(:)+Fl(:)-Fr(:))/(Sr-Sl)

         F(1)=(Sr*Fl(1)-Sl*Fr(1)+Sl*Sr*(uright(1)-uleft(1))) &
         /(Sr-Sl)-ax*uLL(1)

         Fv(1)=(Sr*Fl(2)-Sl*Fr(2)+Sl*Sr*(uright(2)-uleft(2))) &
         /(Sr-Sl)-ax*uLL(2)

         Fv(2)=(Sr*Fl(3)-Sl*Fr(3)+Sl*Sr*(uright(3)-uleft(3))) &
         /(Sr-Sl)-ax*uLL(3)
                 
         Fv(3)=(Sr*Fl(3)-Sl*Fr(3)+Sl*Sr*(uright(3)-uleft(3))) &
         /(Sr-Sl)-ax*uLL(4)
                 
         call rotate_vector(Fv,Fvr,nij,-1)

         F(2)=Fvr(1)
         F(3)=Fvr(2)
                 F(4)=Fvr(3)

         F(5)=(Sr*Fl(5)-Sl*Fr(5)+Sl*Sr*(uright(5)-uleft(5))) &
         /(Sr-Sl)-ax*uLL(5)

         Fv(1)=(Sr*Fl(6)-Sl*Fr(6)+Sl*Sr*(uright(6)-uleft(6))) &
         /(Sr-Sl)-ax*uLL(6)

         Fv(2)=(Sr*Fl(7)-Sl*Fr(7)+Sl*Sr*(uright(7)-uleft(7))) &
         /(Sr-Sl)-ax*uLL(7)
                 
                 Fv(2)=(Sr*Fl(8)-Sl*Fr(8)+Sl*Sr*(uright(8)-uleft(8))) &
         /(Sr-Sl)-ax*uLL(8)

         call rotate_vector(Fv,Fvr,nij,-1)

         F(6)=Fvr(1)
         F(7)=Fvr(2)
         F(8)=Fvr(3)

         endif

         return

         end
                 
         
         subroutine HLLDMHD(i,j,uleft,uright,nij,F,ax,bxc,psi,psic)
        
         USE wp3d_h
         
         implicit none
         
! HLLD Riemann solver for the system of 1D MHD equations
! Miyoshi et al. 
         
         double precision::uleft(nv),uright(nv), &
         ax,F(nv),Fv(ndim),Fvr(ndim),nij(ndim), &
         al,ar,bl,br,valfvenxl,valfvenxr,signbxc,test,tol
         double precision::vxl,vyl,vzl,el,rhol,vxr,vyr,vzr,er,rhor, &
         pr,pl,cl,cr,cmax,Sl,Sr,bxl,bxr,byl,byr,bzl,bzr,plt,prt
         double precision:: vxls1,vxls2,vxrs1,vxrs2,vyls1,vyls2, &
         vyrs1,vyrs2,vzls1,vzrs1,vzls2,vzrs2,rhols1,rhols2,rhors1, &
         rhors2,byls1,byls2,byrs1,byrs2,bzls1,bzls2,bzrs1,bzrs2, &
         Sm,pts,Sml,Smr,els1,els2,ers1,ers2,tiny,psi(nmax)
         double precision:: bxc,psil,psir,psic,pbar
         integer::i,j

         tol=1.0d-10
         tiny=1.0d-10
         
         rhol=uleft(1)
         vxl=uleft(2)/rhol
         vyl=uleft(3)/rhol
         vzl=uleft(4)/rhol
         psil=psi(i)
         el=uleft(5)
         bxl=uleft(6)
         byl=uleft(7)
         bzl=uleft(8)
                 
         if ( barotropic .eqv. .true. ) then
         pl=pbar(i,rhol)
         al=dsqrt(pl/rhol)     
         else
         pl=(gam(i)-1.0d0)*(el-0.5d0*(rhol*(vxl**2+vyl**2+vzl**2)+& 
        (bxl**2+byl**2+bzl**2)))
         al=dsqrt(gam(i)*pl/rhol)    
         endif
        
         plt=pl+(bxl**2+byl**2+bzl**2)/2
       
         rhor=uright(1)
         vxr=uright(2)/rhor
         vyr=uright(3)/rhor
         vzr=uright(4)/rhor
         psir=psi(j)
         er=uright(5)
         bxr=uright(6)
         byr=uright(7)
         bzr=uright(8)
                 
         if ( barotropic .eqv. .true. ) then
         pr=pbar(j,rhor)
         ar=dsqrt(pr/rhor)
         else      
         pr=(gam(j)-1.0d0)*(er-0.5d0*(rhor*(vxr**2+vyr**2+vzr**2)+ &
        (bxr**2+byr**2+bzr**2)))
         ar=dsqrt(gam(j)*pr/rhor)
         endif
                
         prt=pr+(bxr**2+byr**2+bzr**2)/2

         bl=dsqrt((byl**2+bxl**2+bzl**2)/rhol)
         br=dsqrt((byr**2+bxr**2+bzr**2)/rhor)
         valfvenxl=dsqrt(bxl**2/rhol)
         valfvenxr=dsqrt(bxr**2/rhor)

         cl=dsqrt(al**2+bl**2+dsqrt((al**2+bl**2)**2-  &
           4*al**2*valfvenxl**2))/dsqrt(2.0d0)

         cr=dsqrt(ar**2+br**2+dsqrt((ar**2+br**2)**2-  &
           4*ar**2*valfvenxr**2))/dsqrt(2.0d0)

 ! fast magnetosonic wave speeds for the left and right state

         cmax=dmax1(cl,cr)
         
         Sl=dmin1(vxl,vxr)-cmax
         Sr=dmax1(vxl,vxr)+cmax
                 
         if ( iDedner .eqv. .true. ) then
         bxc=(bxl+bxr)/2-(psir-psil)/(2*cmax)
         psic=(psil+psir)/2-cmax*(bxr-bxl)/2
         else
         bxc=(bxl+bxr)/2
         endif
      
         Sm=(prt-plt+rhol*vxl*(Sl-vxl)-rhor*vxr*(Sr-vxr))/ &
         (rhol*(Sl-vxl)-rhor*(Sr-vxr))
                 
         pts=((Sr-vxr)*rhor*plt-(Sl-vxl)*rhol*prt+ &
         rhol*rhor*(Sr-vxr)*(Sl-vxl)*(vxr-vxl))/ &
         ((Sr-vxr)*rhor-(Sl-vxl)*rhol)
                 
         rhols1=rhol*(Sl-vxl)/(Sl-Sm)
         rhors1=rhor*(Sr-vxr)/(Sr-Sm)
         rhols2=rhols1
         rhors2=rhors1

         if (( rhols1 .le. 0.0d0) .or. (rhors1 .le. 0.0d0)) then
         call HLLMHD(i,j,uleft,uright,nij,F,ax,bxc,psi,psic)
         else
                                          
         Sml=Sm-dabs(bxc)/dsqrt(rhols1)
         Smr=Sm+dabs(bxc)/dsqrt(rhors1)
                 
         vxls1=Sm
         vxrs1=Sm
         vxls2=Sm
         vxrs2=Sm
                                 
         test=rhol*(Sl-vxl)*(Sl-Sm)-bxc**2
                                 
         if ( test .le. tol ) then
         vyls1=vyl
         byls1=0.0d0
         vzls1=vzl
         bzls1=0.0d0
         else
         vyls1=vyl-bxc*byl*(Sm-vxl)/test
         vzls1=vzl-bxc*bzl*(Sm-vxl)/test
         byls1=byl*(rhol*(Sl-vxl)**2-bxc**2)/test
         bzls1=bzl*(rhol*(Sl-vxl)**2-bxc**2)/test
         endif
                 
         if ( (.not. barotropic) .eqv. .true. ) then                        
         els1=((Sl-vxl)*el-plt*vxl+pts*Sm+bxc*(vxl*bxc+vyl*byl &
         +vzl*bzl-vxls1*bxc-vyls1*byls1-vzls1*bzls1))/(Sl-Sm)
            endif
                                 
         test=rhor*(Sr-vxr)*(Sr-Sm)-bxc**2
                                 
         if ( test .le. tol ) then
         vyrs1=vyr
         byrs1=0.0d0
         vzrs1=vzr
         bzrs1=0.0d0
         else
         vyrs1=vyr-bxc*byr*(Sm-vxr)/test 
         vzrs1=vzr-bxc*bzr*(Sm-vxr)/test 
         byrs1=byr*(rhor*(Sr-vxr)**2-bxc**2)/test
         bzrs1=bzr*(rhor*(Sr-vxr)**2-bxc**2)/test
                                 
         endif
                 
         if ( (.not. barotropic) .eqv. .true. ) then                        
         ers1=((Sr-vxr)*er-prt*vxr+pts*Sm+bxc*(vxr*bxc+vyr*byr+ &
         vzr*bzr-vxrs1*bxc-vyrs1*byrs1-vzrs1*bzrs1))/(Sr-Sm)
         endif
                 
         signbxc=bxc/dabs(bxc+tiny)
                 
         vyls2=(dsqrt(rhols1)*vyls1+dsqrt(rhors1)*vyrs1+ &
         (byrs1-byls1)*signbxc)/(dsqrt(rhols1)+dsqrt(rhors1))
                                 
         vzls2=(dsqrt(rhols1)*vzls1+dsqrt(rhors1)*vzrs1+ &
         (bzrs1-bzls1)*signbxc)/(dsqrt(rhols1)+dsqrt(rhors1))
                 
          byls2=(dsqrt(rhols1)*byrs1+dsqrt(rhors1)*byls1+ &
          dsqrt(rhols1*rhors1)*(vyrs1-vyls1)*signbxc)/ &
         (dsqrt(rhols1)+dsqrt(rhors1))
                                 
          bzls2=(dsqrt(rhols1)*bzrs1+dsqrt(rhors1)*bzls1+ &
          dsqrt(rhols1*rhors1)*(vzrs1-vzls1)*signbxc)/ &
          (dsqrt(rhols1)+dsqrt(rhors1))
                                
          vyrs2=vyls2
          vzrs2=vzls2
          byrs2=byls2
          bzrs2=bzls2
                  
          if ( (.not. barotropic) .eqv. .true. ) then
                  
          els2=els1-dsqrt(rhols1)*signbxc*(vxls1*bxc+vyls1*byls1 &
          +vzls1*bzls1-vxls2*bxc-vyls2*byls2-vzls2*bzls2)  
          ers2=ers1+dsqrt(rhors1)*signbxc*(vxrs1*bxc+vyrs1*byrs1 &
          +vzrs1*bzrs1-vxrs2*bxc-vyrs2*byrs2-vzrs2*bzrs2)
                  
          endif
                  
         if ( Sl .ge. ax ) then
         
         F(1)=rhol*vxl-ax*rhol
         Fv(1)=rhol*vxl**2+plt-bxc**2-ax*rhol*vxl
         Fv(2)=rhol*vxl*vyl-bxc*byl-ax*rhol*vyl
         Fv(3)=rhol*vxl*vzl-bxc*bzl-ax*rhol*vzl
         
         call rotate_vector(Fv,Fvr,nij,-1)
         
! the flux vector related to the velocity is rotated back to 
! the lab frame to update the velocity vectors in the lab frame
         
         F(2)=Fvr(1)
         F(3)=Fvr(2)
         F(4)=Fvr(3)
                 
         if ( (.not. barotropic) .eqv. .true. ) then
         F(5)=(el+plt)*vxl-(bxc*vxl+byl*vyl+bzl*vzl)*bxc-ax*el
         endif
                 
         Fv(1)=-ax*bxc
         Fv(2)=byl*vxl-bxc*vyl-ax*byl
         Fv(3)=bzl*vxl-bxc*vzl-ax*bzl
                 
         call rotate_vector(Fv,Fvr,nij,-1)
                 
         F(6)=Fvr(1)
         F(7)=Fvr(2)
         F(8)=Fvr(3)

         endif
         
         if ( Sr .le. ax ) then
         
         F(1)=rhor*vxr-ax*rhor
         Fv(1)=rhor*vxr**2+prt-bxc**2-ax*rhor*vxr
         Fv(2)=rhor*vxr*vyr-bxc*byr-ax*rhor*vyr
         Fv(3)=rhor*vxr*vzr-bxc*bzr-ax*rhor*vzr
         
         call rotate_vector(Fv,Fvr,nij,-1)
         
! the flux vector related to the velocity is rotated back to 
! the lab frame to update the velocity vectors in the lab frame
         
         F(2)=Fvr(1)
         F(3)=Fvr(2)
         F(4)=Fvr(3)
                 
         if ( (.not. barotropic) .eqv. .true. ) then
         F(5)=(er+prt)*vxr-(vxr*bxc+vyr*byr+vzr*bzr)*bxc-ax*er
         endif
                 
         Fv(1)=-ax*bxc
         Fv(2)=byr*vxr-bxc*vyr-ax*byr
         Fv(3)=bzr*vxr-bxc*vzr-ax*bzr
                 
         call rotate_vector(Fv,Fvr,nij,-1)
                 
         F(6)=Fvr(1)
         F(7)=Fvr(2)
         F(8)=Fvr(3)
                 
         endif
         
         if (( Sl .le. ax ) .and. ( ax .le. Sml )) then
         
         F(1)=rhol*vxl+(Sl-ax)*rhols1-Sl*rhol
         Fv(1)=rhol*vxl**2+plt-bxc**2+(Sl-ax)*rhols1*vxls1-Sl*rhol*vxl
         Fv(2)=rhol*vxl*vyl-bxc*byl+(Sl-ax)*rhols1*vyls1-Sl*rhol*vyl
         Fv(3)=rhol*vxl*vzl-bxc*bzl+(Sl-ax)*rhols1*vzls1-Sl*rhol*vzl
                 
         call rotate_vector(Fv,Fvr,nij,-1)
                 
         F(2)=Fvr(1)
         F(3)=Fvr(2)
         F(4)=Fvr(3)
                 
         if ( (.not. barotropic) .eqv. .true. ) then
         F(5)=(el+plt)*vxl-(vxl*bxc+vyl*byl+vzl*bzl)*bxc+(Sl-ax)*els1 &
         -Sl*el
         endif
                 
         Fv(1)=-ax*bxc
         Fv(2)=byl*vxl-bxc*vyl+(Sl-ax)*byls1-Sl*byl
         Fv(3)=bzl*vxl-bxc*vzl+(Sl-ax)*bzls1-Sl*bzl
                       
         call rotate_vector(Fv,Fvr,nij,-1)
                 
         F(6)=Fvr(1)
         F(7)=Fvr(2)
         F(8)=Fvr(3)

         endif  
                 
         if (( Sml .le. ax ) .and. ( ax .le. Sm )) then
         
         F(1)=rhol*vxl+(Sml-ax)*rhols2-(Sml-Sl)*rhols1-Sl*rhol
         Fv(1)=rhol*vxl**2+plt-bxc**2+(Sml-ax)*rhols2*vxls2- &
         (Sml-Sl)*rhols1*vxls1-Sl*rhol*vxl
         Fv(2)=rhol*vxl*vyl-bxc*byl+(Sml-ax)*rhols2*vyls2- &
         (Sml-Sl)*rhols1*vyls1-Sl*rhol*vyl
         Fv(3)=rhol*vxl*vzl-bxc*bzl+(Sml-ax)*rhols2*vzls2- &
         (Sml-Sl)*rhols1*vzls1-Sl*rhol*vzl
                                
         call rotate_vector(Fv,Fvr,nij,-1)
                 
         F(2)=Fvr(1)
         F(3)=Fvr(2)
         F(4)=Fvr(3)
                 
         if ( (.not. barotropic) .eqv. .true. ) then
         F(5)=(el+plt)*vxl-(vxl*bxc+vyl*byl+vzl*bzl)*bxc+(Sml-ax)*els2 &
         -(Sml-Sl)*els1-Sl*el
         endif
                 
         Fv(1)=-ax*bxc
         Fv(2)=byl*vxl-bxc*vyl+(Sml-ax)*byls2-(Sml-Sl)*byls1- &
         Sl*byl
         Fv(3)=bzl*vxl-bxc*vzl+(Sml-ax)*bzls2-(Sml-Sl)*bzls1- &
         Sl*bzl
                 
         call rotate_vector(Fv,Fvr,nij,-1)
                 
         F(6)=Fvr(1)
         F(7)=Fvr(2)
         F(8)=Fvr(3)

         endif  
                 
         if (( Sm .le. ax ) .and. ( ax .le. Smr )) then
         
         F(1)=rhor*vxr+(Smr-ax)*rhors2-(Smr-Sr)*rhors1-Sr*rhor
         Fv(1)=rhor*vxr**2+prt-bxc**2+(Smr-ax)*rhors2*vxrs2- &
         (Smr-Sr)*rhors1*vxrs1-Sr*rhor*vxr
         Fv(2)=rhor*vxr*vyr-bxc*byr+(Smr-ax)*rhors2*vyrs2- &
         (Smr-Sr)*rhors1*vyrs1-Sr*rhor*vyr
         Fv(3)=rhor*vxr*vzr-bxc*bzr+(Smr-ax)*rhors2*vzrs2- &
         (Smr-Sr)*rhors1*vzrs1-Sr*rhor*vzr
         
         call rotate_vector(Fv,Fvr,nij,-1)
                 
         F(2)=Fvr(1)
         F(3)=Fvr(2)
         F(4)=Fvr(3)
                 
         if ( (.not. barotropic) .eqv. .true. ) then
         F(5)=(er+prt)*vxr-(vxr*bxc+vyr*byr+vzr*bzr)*bxc+(Smr-ax)*ers2 &
         -(Smr-Sr)*ers1-Sr*er
         endif
                 
         Fv(1)=-ax*bxc
         Fv(2)=byr*vxr-bxc*vyr+(Smr-ax)*byrs2-(Smr-Sr)*byrs1- &
         Sr*byr
         Fv(3)=bzr*vxr-bxc*vzr+(Smr-ax)*bzrs2-(Smr-Sr)*bzrs1- &
         Sr*bzr
                    
         call rotate_vector(Fv,Fvr,nij,-1)
                 
         F(6)=Fvr(1)
         F(7)=Fvr(2)
         F(8)=Fvr(3)

         endif  

         if (( Smr .le. ax ) .and. ( ax .le. Sr )) then
         
         F(1)=rhor*vxr+(Sr-ax)*rhors1-Sr*rhor
         Fv(1)=rhor*vxr**2+prt-bxc**2+(Sr-ax)*rhors1*vxrs1-Sr*rhor*vxr
         Fv(2)=rhor*vxr*vyr-bxc*byr+(Sr-ax)*rhors1*vyrs1-Sr*rhor*vyr
         Fv(3)=rhor*vxr*vzr-bxc*bzr+(Sr-ax)*rhors1*vzrs1-Sr*rhor*vzr
         
         call rotate_vector(Fv,Fvr,nij,-1)
                 
         F(2)=Fvr(1)
         F(3)=Fvr(2)
         F(4)=Fvr(3)
                 
         if ( (.not. barotropic) .eqv. .true.) then
         F(5)=(er+prt)*vxr-(vxr*bxc+vyr*byr+vzr*bzr)*bxc+(Sr-ax)*ers1 &
        -Sr*er
         endif
                 
         Fv(1)=-ax*bxc
         Fv(2)=byr*vxr-bxc*vyr+(Sr-ax)*byrs1-Sr*byr
         Fv(3)=bzr*vxr-bxc*vzr+(Sr-ax)*bzrs1-Sr*bzr
                 
         call rotate_vector(Fv,Fvr,nij,-1)
                 
         F(6)=Fvr(1)
         F(7)=Fvr(2)
         F(8)=Fvr(3)

         endif 

         endif

         return
         
         end

                 
      SUBROUTINE MKTREE(wprim)
      
! --------------------------------------------------------------
! MKTREE: initialize the tree structure for the force calculation.
! ----------------------------------------------------------------
      
      use wp3d_h

      implicit none

      integer:: I, MKCELL, K, P, IND(1:mxnode),  &
      Q, J, L, M1,M2 
      DOUBLE PRECISION:: XYZMAX, POS0(NDIM), DIST2,DR, &
      wprim(nmax,nv)
      
      mxcell=n
      incell=n+1

!     Expand root volume to enclose all particles.
!     --------------------------------------           
!     ------------------------------
!     Find maximum coordinate value.
!     ------------------------------
      
      XYZMAX = 0.0d0
      RSIZE=1.0d0
	  
      do i = 1, n
      dr=dsqrt(x(i,1)**2+x(i,2)**2+x(i,3)**2)
      xyzmax=dmax1(xyzmax,dr)
      am(i)=vol(i)*wprim(i,1)
      ENDDO
     
      DO WHILE (XYZMAX .GE. RSIZE/2.0)
      RSIZE = 2.0 * RSIZE
      ENDDO
      
!     Load bodies into the tree.
!       --------------------------
!     ---------------------------------------
!     Deallocate current tree, begin new one.
!       ---------------------------------------
      NCELL = 0
      ROOT = MKCELL()

!       ------------------------------------------
!       Initialize midpoint and size of root cell.
!       ------------------------------------------
        DO K = 1, NDIM
        MID(ROOT,K) = 0.0d0
      ENDDO
        CLSIZE(ROOT) = RSIZE
!       ---------------------------------------------
!     Load bodies into the new tree, one at a time.
!       ---------------------------------------------
      DO P = 1, n
      if ( .not. exterior(P)) then
      CALL LDBODY(P)
      endif
      ENDDO
      
!       ------------------------------------------------
!     Compute masses, center of mass coordinates, etc.
!       ------------------------------------------------
!       ---------------------------------------
!     List cells in order of decreasing size.
!       ---------------------------------------
      CALL BFLIST(IND)
!     --------------------------------------------
!    Loop processing cells from smallest to root.
!     --------------------------------------------
      DO I = NCELL, 1, -1
      P = IND(I)
!     --------------------------------------------------------------
!     Zero accumulators for this cell.  A temporary variable is used
!       for the c. of m. so as to preserve the stored midpoints.
!      --------------------------------------------------------------

!     Compute the mass and the center of mass position of the particles in cell P 

      AM(P) = 0.0d0
      DO K = 1, NDIM
      POS0(K) = 0.0d0
      ENDDO
!     -------------------------------------------------------------
!     Compute cell properties as sum of properties of its subcells.
!      -------------------------------------------------------------
      DO J = 1, NSUBC
      Q = SUBP(P,J)
!      ------------------------------
!     Only access cells which exist.
!      ------------------------------
      IF (Q .NE. NULL) THEN
!     -------------------------------------------------------
!     Sum properties of subcells to obtain values for cell P.
!     -------------------------------------------------------
      AM(P) = AM(P) + AM(Q)
      DO K = 1, NDIM
      POS0(K) = POS0(K) + AM(Q) * X(Q,K)
      ENDDO
      ENDIF
      ENDDO
!     --------------------------------------------------------
!     Normalize center of mass coordinates by total cell mass.
!      --------------------------------------------------------
      DO K = 1, NDIM
      POS0(K) = POS0(K) / AM(P)
      ENDDO 
!     -----------------------------------------------------------------
!     Check tree, compute cm-to-mid distance, and assign cell position.
!     -----------------------------------------------------------------
      DIST2 = 0.0
      DO K = 1, NDIM
      IF (POS0(K) .LT. MID(P,K) - CLSIZE(P)/2.0 .OR.  &
          POS0(K) .GE. MID(P,K) + CLSIZE(P)/2.0) THEN
      WRITE(6, '(/,1X,''TREE ERROR'',2I6,3E14.6)')  &    
       P, K, POS0(K), MID(K,P), CLSIZE(P)
      CALL OUTERR(' HACKCM: TREE STRUCTURE ERROR')
      ENDIF
      DIST2 = DIST2 + (POS0(K) - MID(P,K))**2
!       --------------------------------------------------------
!       Copy cm position to cell.  This overwrites the midpoint.
!     --------------------------------------------------------
        X(P,K) = POS0(K)
      ENDDO

!     ------------------------------------------------------------
!       Assign critical radius for cell, adding offset from midpoint
!       for more accurate forces.  This overwrites the cell size.
!       ------------------------------------------------------------
      RCRIT2(P) = (CLSIZE(P) / THETA + SQRT(DIST2))**2

      ENDDO
!       ----------------------------------------
!       Compute quadrupole moments, if required.
!       ----------------------------------------
         IF (USQUAD) THEN
!     -------------------------------
!     Loop processing cells as above.
!     -------------------------------
      DO I = NCELL, 1, -1
      P = IND(I)
!      --------------------------------------------
!     Zero accumulator for quad moments of cell P.
!     --------------------------------------------
      DO K = 1, NQUAD
      QUAD(P,K) = 0.0
      ENDDO

!     --------------------------------
!     Loop over descendents of cell P.
!      --------------------------------
        DO J = 1, NSUBC
        Q = SUBP(P,J)
        IF (Q .NE. NULL) THEN
!     --------------------------------------------------------
!     Sum properties of subcell Q to obtain values for cell P.
!     --------------------------------------------------------
      DO M1 = 1, MIN(2,NDIM)
      DO M2 = M1, NDIM
      L = (M1-1) * (NDIM-1) + M2
      QUAD(P,L) = QUAD(P,L) + 3.0 * AM(Q) * &
       (X(Q,M1) - X(P,M1)) * (X(Q,M2) - X(P,M2))
      IF (M1 .EQ. M2) THEN
      DO K = 1, NDIM
      QUAD(P,L) = QUAD(P,L) - AM(Q) *  &
       (X(Q,K) - X(P,K))**2
      ENDDO
      ENDIF
!      -------------------------------------------
!     If Q itself is a cell, add its moments too.
!     -------------------------------------------
      IF (Q .GE. INCELL) &
        QUAD(P,L) = QUAD(P,L) + QUAD(Q,L)
      ENDDO
      ENDDO
        ENDIF
      ENDDO
      ENDDO
        ENDIF
      
      RETURN
      
      END
      
      
! -------------------------------------------------------
! SBINDX: compute subcell index for node P within cell Q.
! -------------------------------------------------------
 
      integer FUNCTION SBINDX(P, Q)
      
      use wp3d_h
      
      integer::P, Q, K
 
!       ---------------------------------------------------
!       Initialize subindex to point to lower left subcell.
!       ---------------------------------------------------
      SBINDX = 1
!       ---------------------------------
!       Loop over all spatial dimensions.
!       ---------------------------------
      DO K = 1, NDIM
      IF (X(P,K) .GE. MID(Q,K))  & 
        SBINDX = SBINDX + 2 ** (NDIM - K)
      ENDDO
      END
 
! ---------------------------------------------------------
! MKCELL: function to allocate a cell, returning its index.
! ---------------------------------------------------------
 
      integer FUNCTION MKCELL()
      
      use wp3d_h
 
      integer::I
 
!     ----------------------------------------------------------
!     Terminate simulation if no remaining space for a new cell.
!     ----------------------------------------------------------

      IF (NCELL .GE. MXCELL) &  
       CALL OUTERR(' MKCELL: NO MORE MEMORY')
!       ----------------------------------------------------
!       Increment cell counter, initialize new cell pointer.
!       ----------------------------------------------------
      NCELL = NCELL + 1
      MKCELL = NCELL + n
!       --------------------------------------
!       Zero pointers to subcells of new cell.
!       --------------------------------------
      DO I = 1, NSUBC
      SUBP(MKCELL,I) = NULL
      END DO

        RETURN 
      END
      
      
      SUBROUTINE LDBODY(P)
      
 ! ----------------------------------------
! LDBODY: load particle P into tree structure.
! ----------------------------------------
        
      use wp3d_h

      integer:: P,Q, QIND, SBINDX, MKCELL, C1, K, P0

!     ---------------------------------------------
!     Start Q,QIND pair in correct subcell of root.
!     ---------------------------------------------
      Q = ROOT
      QIND = SBINDX(P, Q)
!     -----------------------------------------------------
!     Loop descending tree until an empty subcell is found.
!     -----------------------------------------------------
      DO WHILE (SUBP(Q, QIND) .NE. NULL)
!     --------------------------------------
!     On reaching another body, extend tree.
!     --------------------------------------
      IF (SUBP(Q, QIND) .LT. INCELL) THEN
!      -------------------------------------------
!     Allocate an empty cell to hold both bodies.
!      -------------------------------------------
      C1 = MKCELL()
!      ------------------------------------------------------
!     Locate midpoint of new cell wrt. parent, and set size.
!      ------------------------------------------------------
      DO K = 1, NDIM
      IF (X(P,K) .GE. MID(Q,K)) THEN
      MID(C1,K) = MID(Q,K) + CLSIZE(Q)/4.0
      ELSE
      MID(C1,K) = MID(Q,K) - CLSIZE(Q)/4.0
      ENDIF
      ENDDO
      CLSIZE(C1) = CLSIZE(Q) / 2.0
!     ------------------------------------------------------
!     Store old body in appropriate subcell within new cell.
!      ------------------------------------------------------
      P0 = SUBP(Q, QIND)
      SUBP(C1, SBINDX(P0, C1)) = P0
!     ---------------------------------------------
!     Link new cell into tree in place of old body.
!      ---------------------------------------------
      SUBP(Q, QIND) = C1
      ENDIF
!     --------------------------------------------------------
!     At this point, the node indexed by Q,QIND is known to be
!       a cell, so advance to the next level of tree, and loop.
!     --------------------------------------------------------
      Q = SUBP(Q, QIND)
      QIND = SBINDX(P, Q)
      ENDDO
!     ---------------------------------------------
!     Found place in tree for P, so store it there.
!     ---------------------------------------------
      SUBP(Q, QIND) = P

      RETURN 
      END
 
      
      SUBROUTINE OUTERR(MSG)
      CHARACTER*(*) MSG

!     Write error message to status file 

      WRITE (6, '(/,A)') MSG
     
      close(2)
      
      STOP 
      
! terminate of the execution of the program when a tree error occurs

      END
           
! -----------------------------------------------------------------
! BFLIST: list cells in breadth-first order, from largest (root) to
! smallest.  Thanks to Jun Makino for this elegant routine.
! -----------------------------------------------------------------
 
      SUBROUTINE BFLIST(IND)
      
      use wp3d_h
 
      integer:: IND(*)
 
      integer:: FACELL, LACELL, NACELL, K, I
 
!       -----------------------------------------
!       Start scan with root as only active cell.
!       -----------------------------------------
        IND(1) = ROOT
        FACELL = 1
        LACELL = 1
!       -----------------------------------
!       Loop while active cells to process.
!       -----------------------------------
      DO WHILE (FACELL .LE. LACELL)
!       ----------------------------------------------
!       Start counting active cells in next iteration.
!       ---------------------------------------------- 
        NACELL = LACELL
!       ---------------------------------------
!       Loop over subcells of each active cell.
!       ---------------------------------------
        DO K = 1, NSUBC
        DO I = FACELL, LACELL
!        -------------------------------------------
!     Add all cells on next level to active list.
!        -------------------------------------------
        IF (SUBP(IND(I),K) .GE. INCELL) THEN
        NACELL = NACELL + 1
        IND(NACELL) = SUBP(IND(I),K)
        ENDIF
      ENDDO
      ENDDO
!       ------------------------------------------------------
!       Advance first and last active cell indicies, and loop.
!       ------------------------------------------------------
        FACELL = LACELL + 1
        LACELL = NACELL
      ENDDO
!       --------------------------------------------------
!       Above loop should list all cells; check the count.
!       --------------------------------------------------
        IF (NACELL .NE. NCELL)  &
     CALL OUTERR('  BFLIST: INCONSISTENT CELL COUNT')

      RETURN 
      END
      
      
      
      SUBROUTINE TRN(i,iopt,hp)
      
      use wp3d_h

! TRN: recursive routine to walk the tree and find the neighbours of a particle 
! adapted from Joshua Barnes

        integer MXSPTR
        PARAMETER(MXSPTR = 256)
        double precision ::  DX, DY, DZ,DR2,DR,HP
        double precision ::  XP, YP, ZP, RSUPPORT
        integer I,Q, SPTR, STACK(MXSPTR), K,iopt
        
        XP=X(i,1)
        YP=X(i,2)
        ZP=X(i,3)
         
        nbn=0
                                   
!       ----------------------------------
!       Push the root cell onto the stack.
!       ----------------------------------
      SPTR = 1
      STACK(SPTR) = ROOT
!       -------------------------------------
!       Loop while nodes on stack to process.
!       -------------------------------------
      DO WHILE (SPTR .GT. 0)
!         --------------------------
!       Pop node off top of stack.
!         --------------------------
      Q = STACK(SPTR)
      SPTR = SPTR - 1

!     Compute distance to center-of-mass of node Q.

      DX = XP - X(Q,1)
      DY = YP - X(Q,2)
      DZ = ZP - X(Q,3)

  !    call modbound(dx,dy,dz)

      DR2 = DX*DX + DY*DY + DZ*DZ
      DR=DSQRT(DR2)  
            
      select case ( iopt )
             
        case (1)
             
        RSUPPORT=2*HP
           
        case (2)
          
        if ( Q .LT. INCELL ) THEN

        RSUPPORT=2*DMAX1(HP,H(Q))
          
        ELSE
          
        RSUPPORT=2*HP
          
        ENDIF
             
        end select 
          
! two options for neighbour search: 1) all neighbours within 2 smoothing lengths of the current particle and 
! 2) all *interacting* neighbours for which rij < 2*max(hi,hj), even if rij > 2*hi   

!      -------------------------------
!         Classify Q as a body or a cell.
!      -------------------------------

        IF  (Q .LT. INCELL ) THEN
!           -----------------------------------
!       A body: check if it is a neighbor 
!           -----------------------------------      

        IF ( DR .lt. RSUPPORT ) THEN
                
!     Body Q is indeed a neighbor of P

        NBN=NBN+1
        NNI(I,NBN)=Q
                             
        ENDIF 
          
      ELSE
!         
!           --------------------------------------------
!     A cell: test if the cell can contain neighbors
!           --------------------------------------------
      IF (DR2 .le. CLSIZE(Q)+RSUPPORT) THEN
!             -----------------------------------
!       examine children of cell if accepted
!             -----------------------------------
      DO K = 1, NSUBC

!               --------------------------------------
!    Push existing children onto the stack.
!               --------------------------------------
      IF (SUBP(Q,K) .NE. NULL) THEN
      IF (SPTR .GE. MXSPTR)  &
        CALL OUTERR(' TRWALK: STACK OVERFLOW')
      SPTR = SPTR + 1
      STACK(SPTR) = SUBP(Q,K)
      ENDIF
      ENDDO
        ENDIF
      ENDIF
        ENDDO

        RETURN
      END  
          
          
      SUBROUTINE TRN2(xp,yp,zp,iopt,hp)
      
      USE wp3d_h

! TRN: recursive routine to walk the tree and find the neighbours of a particle 
! adapted from Joshua Barnes

        integer::MXSPTR
        PARAMETER(MXSPTR = 256)
        double precision ::  DX, DY, DZ,DR2,DR,HP
        double precision ::  XP, YP, ZP, RSUPPORT
        integer::Q, SPTR, STACK(MXSPTR), K,iopt
        
        nbn=0
                                   
!       ----------------------------------
!       Push the root cell onto the stack.
!       ----------------------------------
      SPTR = 1
      STACK(SPTR) = ROOT
!       -------------------------------------
!       Loop while nodes on stack to process.
!       -------------------------------------
      DO WHILE (SPTR .GT. 0)
!         --------------------------
!       Pop node off top of stack.
!         --------------------------
      Q = STACK(SPTR)
      SPTR = SPTR - 1

!     Compute distance to center-of-mass of node Q.

      DX = XP - X(Q,1)
      DY = YP - X(Q,2)
      DZ = ZP - X(Q,3)
      DR2 = DX*DX + DY*DY + DZ*DZ
      DR=DSQRT(DR2)  
            
      select case ( iopt )
             
        case (1)
             
        RSUPPORT=2*HP
           
        case (2)
          
        if ( Q .LT. INCELL ) THEN

        RSUPPORT=2*DMAX1(HP,H(Q))
          
        ELSE
          
        RSUPPORT=2*HP
          
        ENDIF
             
        end select 
          
! two options for neighbour search: 1) all neighbours within 2 smoothing lengths of the current particle and 
! 2) all *interacting* neighbours for which rij < 2*max(hi,hj), even if rij > 2*hi   

!      -------------------------------
!         Classify Q as a body or a cell.
!      -------------------------------

        IF  (Q .LT. INCELL ) THEN
!           -----------------------------------
!       A body: check if it is a neighbor 
!           -----------------------------------      

        IF ( DR .lt. RSUPPORT ) THEN
                
!     Body Q is indeed a neighbor of P

        NBN=NBN+1
        NN(NBN)=Q
                             
        ENDIF 
          
      ELSE
!         
!           --------------------------------------------
!     A cell: test if the cell can contain neighbors
!           --------------------------------------------
        IF (DR2 .le. CLSIZE(Q)+RSUPPORT ) THEN
!             -----------------------------------
!       examine children of cell if accepted
!             -----------------------------------
      DO K = 1, NSUBC

!               --------------------------------------
!    Push existing children onto the stack.
!               --------------------------------------
      IF (SUBP(Q,K) .NE. NULL) THEN
      IF (SPTR .GE. MXSPTR)  &
        CALL OUTERR(' TRWALK: STACK OVERFLOW')
      SPTR = SPTR + 1
      STACK(SPTR) = SUBP(Q,K)
      ENDIF
      ENDDO
        ENDIF
      ENDIF
        ENDDO

        RETURN
      END         
          

      SUBROUTINE TRG(i,gaccxi,gaccyi,gacczi,phii)
          
!--------------------------------------------------------------
! TRG: recursive routine to walk the tree computing forces on
! particle P
! --------------------------------------------------------------

      USE wp3d_h
 
      integer:: i
      
      integer::MXSPTR
      PARAMETER(MXSPTR = 256)
      double precision::PHI0, ACC0(NDIM), POS0(NDIM)
      double precision:: DX, DY, DZ, DR2, DR2INV, DRINV, PHIM, DR5INV
      double precision:: FG,GG,DR,PHIQ,gaccxi,gaccyi,gacczi,phii
      integer:: Q, SPTR, STACK(MXSPTR), K
      LOGICAL SKPSLF    
        
!       ----------------------------------------------------------
!       Zero potential and acceleration for subsequent summations.
!       ----------------------------------------------------------
      PHI0 = 0.0
      ACC0(1) = 0.0
      ACC0(2) = 0.0
      ACC0(3) = 0.0

!     -----------------------------------------
!     Copy position of this particle for quicker reference.
!     -----------------------------------------

      POS0(1)=X(i,1)
      POS0(2)=X(i,2)
      POS0(3)=X(i,3)

!       ----------------------------------
!    Push the root cell onto the stack.
!       ----------------------------------
      SPTR = 1
      STACK(SPTR) = ROOT
!       ------------------------------------
!       Loop while nodes on stack to process.
!       -------------------------------------
      DO WHILE (SPTR .GT. 0)
!       --------------------------
!       Pop node off top of stack.
!       --------------------------
      Q = STACK(SPTR)
      SPTR = SPTR - 1
!     ---------------------------------------------
!     Compute distance to center-of-mass of node Q.
!      ---------------------------------------------
      DX = POS0(1) - X(Q,1)
      DY = POS0(2) - X(Q,2)
      DZ = POS0(3) - X(Q,3)
      DR2 = DX*DX + DY*DY + DZ*DZ
      DR=SQRT(DR2)  

!     -------------------------------
!       Classify Q as a body or a cell.
!     -------------------------------
      IF (Q .LT. INCELL) THEN
!      -----------------------------------
!       A body: check for self-interaction.
!       -----------------------------------
        IF (Q .NE. i) THEN
!     Compute body-body interaction.
!      ------------------------------
!    Spline softened gravitational potential and acceleration 
!    according to Hernquist and Katz 
      
      PHIM = AM(Q) * (FG(DR,H(I))+FG(DR,H(Q)))/2 
      PHI0 = PHI0 - PHIM
      PHIM = AM(Q) * (GG(DR,H(I))+GG(DR,H(Q)))/2
      
      ACC0(1) = ACC0(1) - PHIM * DX
      ACC0(2) = ACC0(2) - PHIM * DY
      ACC0(3) = ACC0(3) - PHIM * DZ
      
      ELSE
!       -------------------------------------------
!     Remember that self-interaction was skipped.
!             -------------------------------------------
          SKPSLF = .TRUE.
      ENDIF
      ELSE
!      --------------------------------------------
!     A cell: test if interaction can be accepted.
!      --------------------------------------------
        IF (DR2 .GE. RCRIT2(Q)) THEN
!             ----------------------------------------
!    Accepted: compute body-cell interaction.
!             ----------------------------------------
!    monopole contribution 

      DR2INV = 1.0 / DR2
      DRINV = 1.0 / DR
      PHIM = AM(Q)*DRINV
      PHI0 = PHI0 - PHIM
      PHIM = PHIM*DR2INV
      ACC0(1) = ACC0(1) - PHIM * DX
      ACC0(2) = ACC0(2) - PHIM * DY
      ACC0(3) = ACC0(3) - PHIM * DZ

!        ------------------------------------
!       Optionally include quadrupole terms.
!        ------------------------------------
        IF (USQUAD) THEN
        DRINV=1./DR
        DR2INV=1./DR2
      DR5INV = DR2INV * DR2INV * DRINV
      PHIQ = DR5INV *   &
       (0.5d0 * ((DX*DX - DZ*DZ) * QUAD(Q,1) +   &
        (DY*DY - DZ*DZ) * QUAD(Q,4)) +  & 
        DX*DY * QUAD(Q,2) + DX*DZ * QUAD(Q,3) +  &
       DY*DZ * QUAD(Q,5))
        PHI0 = PHI0 - PHIQ
        PHIQ = 5.0 * PHIQ * DR2INV
        ACC0(1) = ACC0(1) - PHIQ*DX + DR5INV *  &
       (DX*QUAD(Q,1) + DY*QUAD(Q,2) + DZ*QUAD(Q,3))
        ACC0(2) = ACC0(2) - PHIQ*DY + DR5INV * &
       (DX*QUAD(Q,2) + DY*QUAD(Q,4) + DZ*QUAD(Q,5))
        ACC0(3) = ACC0(3) - PHIQ*DZ + DR5INV *  &
       (DX*QUAD(Q,3) + DY*QUAD(Q,5) -  & 
      DZ*(QUAD(Q,1) + QUAD(Q,4)))
        ENDIF
      ELSE
!       -----------------------------------
!       Rejected: examine children of cell.
!       -----------------------------------
      DO K = 1, NSUBC
!      --------------------------------------
!      Push existing children onto the stack.
!      --------------------------------------
      IF (SUBP(Q,K) .NE. NULL) THEN
      IF (SPTR .GE. MXSPTR)  &  
       CALL OUTERR(' TRWALK: STACK OVERFLOW')
      SPTR = SPTR + 1
      STACK(SPTR) = SUBP(Q,K)
      ENDIF
      ENDDO
          ENDIF
      ENDIF
      ENDDO

!     ---------------------------------------------------
!     Check that self-interaction was explicitly skipped.
!     ---------------------------------------------------
!      IF (.NOT. SKPSLF) & 
!      CALL OUTERR(' TRWALK: MISSED SELF-INTERACTION')

!     ----------------------------------------------
!     Copy total potential and acceleration to body.
!     ----------------------------------------------
 
      PHII=PHI0
      GACCXI = ACC0(1)
      GACCYI = ACC0(2)
      GACCZI = ACC0(3)      

        RETURN
        END   


! f and g functions from Hernquist and Katz 

        double precision FUNCTION FG(r,hg)

        double precision u,r,hg

        u=r/hg

        IF ( U .lt. 1.0d0 ) THEN

        FG=u**2/3-3*u**4/20+u**5/20
        FG=-2*FG/hg+7/(5*hg)

        ELSE IF ( U .lt. 2.0d0 ) THEN

        FG=4*u**2/3-u**3+3*u**4/10-u**5/30
        FG=-FG/hg-1./(15*r)+8./(5*hg)

        ELSE

        FG=1./r

        ENDIF 

        return 

        END 


! f and g functions from Hernquist and Katz 

        double precision FUNCTION DFG(r,hg)

        double precision u,r,hg

        u=r/hg

        IF ( U .lt. 1.0d0 ) THEN

        DFG=(2*u**2-3*u**4/2.+3*u**5/5-7./5)/hg**2

        ELSE IF ( U .lt. 2.0d0 ) THEN

        DFG=(4*u**2-4*u**3+3*u**4/2.-u**5/5-8./5.)/hg**2

        ELSE

        DFG=0.0d0

        ENDIF 

        return 

        END 


        double precision FUNCTION GG(r,hg)

        double precision u,r,hg

        u=r/hg

        IF ( U .lt. 1.0d0 ) THEN

        GG=(4.0d0/3-6*u**2/5+u**3/2)/hg**3
        
        ELSE IF ( U .lt. 2.0d0 ) THEN

        GG=(-1.0d0/15+8*u**3/3-3*u**4+6*u**5/5-u**6/6)/r**3

        ELSE

        GG=1.0d0/r**3

        ENDIF 

        return 

        END 
                
        DOUBLE PRECISION FUNCTION PERIODIC1(XP,BOX)

! calculates distances taking account of the periodic boundary conditions

        use wp3d_h

        DOUBLE PRECISION::XP,box,boxhalf
         
         boxhalf=box/2.

         DO while ( XP .gt. BOXHALF )
         XP=XP-BOX
         ENDDO
         DO while ( XP .lt. -BOXHALF )
         XP=XP+BOX
         ENDDO 

         periodic1=xp
           
         RETURN
           
         END   
                                 
         DOUBLE PRECISION FUNCTION PERIODIC2(XP,BOX)

! calculates distances taking account of periodic boundary conditions

         use wp3d_h

         DOUBLE PRECISION::XP,box

         DO while ( XP .gt. BOX)
         XP=XP-BOX
         ENDDO
         DO while ( XP .lt. 0.0d0 )
         XP=XP+BOX
         ENDDO 

         periodic2=xp
           
         RETURN
           
         END  
         
         INTEGER FUNCTION modindx(idx,nc)
         
         use wp3d_h
         
         integer::idx,nc
                 
         DO while ( idx .ge. nc)
         idx=idx-nc
         ENDDO
         DO while ( idx .lt. 0 )
         idx=idx+nc
         ENDDO 
         
         modindx=idx
                 
         RETURN
         
         END 
         
         
         SUBROUTINE modbound(dx,dy,dz)
          
         USE wp3d_h
         
         implicit none
          
         double precision::dx,dy,dz
         
         IF (abs(dx).GT.0.5*boxx) then
         dx = dx - boxx*dsign(1.0d0,dx)
         endif
         IF (abs(dy).GT.0.5*boxy) then
         dy = dy - boxy*dsign(1.0d0,dy)
         endif
         IF (abs(dz).GT.0.5*boxz) then
         dz = dz - boxz*dsign(1.0d0,dz)
         endif
        
      RETURN
          
      END SUBROUTINE modbound

          
         
      subroutine build_linked_list
      
      USE wp3d_h

! construction of a Linked list for neighbour searching on a single 3D grid
         
      LOGICAL::EXC
      double precision::hpav,h1,h2,h3
    
      INTEGER::i,k,l,m,npout

! Determine system box and mean HP:
          XMIN=1.E30
          YMIN=1.E30
          XMAX=-1.E30
          YMAX=-1.E30
          ZMIN=1.0E30
          ZMAX=-1.0E30
          HPAV=0.0d0
                  
          DO I=1,N
          XMIN=MIN(XMIN,X(I,1))
          YMIN=MIN(YMIN,X(I,2))
          ZMIN=MIN(ZMIN,X(I,3))
          XMAX=MAX(XMAX,X(I,1))
          YMAX=MAX(YMAX,X(I,2))
          ZMAX=MAX(ZMAX,X(I,3))
          HPAV=HPAV+H(I)
          ENDDO
                  
      HPAV=HPAV/FLOAT(N)
      HH=HPAV*1.3d0

!     (This should give near-optimal search times)
!     6/29/98: The factor of 1.3 seems to give search times for typical
!    cases which are about 10-20% shorter than if HH=HPAV.
 
! Calculate number of cells in each dimension:

      NCX=DINT((XMAX-XMIN)/HH)+1
      NCY=DINT((YMAX-YMIN)/HH)+1
      NCZ=DINT((ZMAX-ZMIN)/HH)+1

 !    write(6,*) ncx,ncy,ncz,hh
      
! Readjust if maximum number allowed is exceeded:
! (Note that this version assumes centering on the origin!)
      EXC=.FALSE.
      IF (NCX.GT.NCMX) THEN 
!         write(6,*) &
!         'Warning! NELIST: EXC TRUE in x-direction',NCX,NCMX
         EXC=.TRUE.
         NCX=NCMX
      ENDIF
      IF (NCY.GT.NCMY) THEN 
!         write(6,*) &
!         'Warning! NELIST: EXC TRUE in y-direction',NCY,NCMY
         EXC=.TRUE.
         NCY=NCMY
      ENDIF
        IF (NCZ.GT.NCMZ) THEN
!         write(6,*) &
!         'Warning! NELIST: EXC TRUE in z-direction',NCZ,NCMZ
         EXC=.TRUE.
         NCZ=NCMZ
      ENDIF
      
      IF (EXC) THEN
      H1=dabs((xmax-xmin)/dfloat(ncx-1))
      H2=dabs((ymax-ymin)/dfloat(ncy-1))
      H3=dabs((zmax-zmin)/dfloat(ncz-1))
      HH=DMAX1(H1,H2,H3)
      ENDIF

! Initialize Head-Of-Cell array:
      DO M=0,NCMZ
        DO L=0,NCMY
          DO K=0,NCMX
            HOC(K,L,M)=0
          ENDDO
        ENDDO
      ENDDO

! Build linked lists (cf. Hockney and Eastwood):
      NPOUT=0
      IF (EXC) THEN
!   (test for particles outside of the grid)
         DO I=1,N
            K=floor((X(I,1)-XMIN)/HH)
            L=floor((X(I,2)-YMIN)/HH)
            M=floor((X(I,3)-ZMIN)/HH)
            IF ((K.LT.0).OR.(K.GT.NCMX-1).OR. &
                (L.LT.0).OR.(L.GT.NCMY-1).OR. &
                (M.LT.0).OR.(M.GT.NCMZ-1)) THEN
               OUT(I)=.TRUE.
               NPOUT=NPOUT+1
            ELSE
               OUT(I)=.FALSE.
               LL(I)=HOC(K,L,M)
               HOC(K,L,M)=I
            ENDIF
         ENDDO
      ELSE
!   (no need to check for out of bounds indices)
         DO I=1,N
          K=floor((X(I,1)-XMIN)/HH)
          L=floor((X(I,2)-YMIN)/HH)
          M=floor((X(I,3)-ZMIN)/HH)
  !         K=INT((X(I,1)-XMIN)/HH)+1
  !         L=INT((X(I,2)-YMIN)/HH)+1
  !         M=INT((X(I,3)-ZMIN)/HH)+1
            OUT(I)=.FALSE.
            LL(I)=HOC(K,L,M)
            HOC(K,L,M)=I
         ENDDO
      ENDIF
      
      IF (NPOUT.NE.0) WRITE (6,fmt="(A,I5)") &
          'NELIST: WARNING !!! NPOUT=',NPOUT

      RETURN
      END  subroutine build_linked_list
      
      
        SUBROUTINE neighbourlists(wprim)
      
        use wp3d_h
                
        implicit none
               
 ! determines smoothing length and the neighbour list of of the particles, including all
 ! *interactive* pairs, and also computes the volume and mass of the particles
            
        integer::i,j,k,iopt,itab
        double precision::hpi,voli,weight,u2,dx,dy,dz, &
	    DFG,t1,t2,h1,h2,omegai,dr,wprim(nmax,nv)
        
        if ( hvar .eqv. .true. ) then
                
        select case ( hchoice )
                
        case (1)
                
! determine the smoothing length so that the number of neighbours is constant

        do i=1,n
                
        hpi=h(i)
        
        call find_smoothing_length(i,hpi)

        h(i)=hpi

        vol(i)=weight(i,h(i))

        if ( igrav ) then
        am(i)=vol(i)*wprim(i,1)
        endif

        enddo
                     
        case (2)
                
! determine consistent smoothing lenghts and particle weights 
! by means of a Newton-Raphson method 

        do i=1,n
               
        hpi=h(i)

        call findhvol(i,voli,hpi)
                
        h(i)=hpi
        vol(i)=voli
        if ( igrav ) then
        am(i)=vol(i)*wprim(i,1)
        endif

        enddo
                
        end select 
                
        else

        do i=1,n
                                   
        if (( iLagrangian .eqv. .false. ) .and. &
        ( nit .eq. 1)) then

        call findhvol(i,voli,h(i))

        vol(i)=voli
        if ( igrav ) then
        am(i)=vol(i)*wprim(i,1)
        endif

        endif
                                       
        enddo
                
        endif 
                
        iopt=2
                
        do i=1,n
                
        call LLN(i,iopt,h(i))
            
        nne(i)=nbn

        do k=1,nne(i)
        nni(i,k)=nn(k)
        enddo 

        enddo
	   
! neighbour lists
! if hvar=.false. than h(i)=const 

        if ( igrav .eqv. .true. ) then
	   
        do i=1,n
                
        hpi=h(i)
	    voli=vol(i)
	   
! compute dzeta quantities to deal with spatially varying gravitational softening lengths

        t1=0.0d0
        t2=0.0d0
		
        do k=1,nne(i)
		
        j=nni(i,k) 

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

       omegai=1.0d0-(voli*hpi/float(ndim))*t1
       dzeta(i)=am(i)*voli*hpi*t2/omegai/float(ndim)

       enddo

       endif
                
        
       return
        
       end
        
        
        double precision function weight(i,hpi)
        
        use wp3d_h
        
        implicit none 
        
        integer:: i,j,k,itab
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

       SUBROUTINE LLN(i,iopt,hp)

        USE wp3d_h

        INTEGER::i,iopt,k1,k2,l1,l2,m1,m2
        INTEGER::kne,lne,ine,mne,KINDX,LINDX, &
        MINDX,idx1,idx2,idx3,modindx

        double precision::hp,r,dx,dy,dz,rsupport
        double precision::xp,yp,zp

! Calculate the neighbour list of a particle

          xp=x(i,1)
          yp=x(i,2)
          zp=x(i,3)
          
          nbn=0

! Find indices of cells containing potential neighbours:

         if ( iperiodic ) then
         K1=floor((XP-nhp*HP-XMIN)/HH)          
         K2=floor((XP+nhp*HP-XMIN)/HH)
         L1=floor((YP-nhp*HP-YMIN)/HH)
         L2=floor((YP+nhp*HP-YMIN)/HH)
         M1=floor((ZP-nhp*HP-ZMIN)/HH)
         M2=floor((ZP+nhp*HP-ZMIN)/HH)
            else
          K1=MAX(1,INT((xp-np*HP-XMIN)/HH)+1)
          K2=MIN(NCMX,INT((xp+np*HP-XMIN)/HH)+1)
          L1=MAX(1,INT((yp-np*HP-YMIN)/HH)+1)
          L2=MIN(NCMY,INT((yp+np*HP-YMIN)/HH)+1)
          M1=MAX(1,INT((zp-np*HP-ZMIN)/HH)+1)
          M2=MIN(NCMZ,INT((zp+np*HP-ZMIN)/HH)+1)
            endif

!     Look for neighbours in all these cells:

            DO MNE=M1,M2
            DO LNE=L1,L2
            DO KNE=K1,K2
           
                idx1=kne
                idx2=lne
                idx3=mne
                 
                 if ( iperiodic ) then

                 KINDX=modindx(idx1,ncx)
                 LINDX=modindx(idx2,ncy)
                 MINDX=modindx(idx3,ncz)
                   
                 else
                 
                 KINDX=IDX1
                 LINDX=IDX2
                 MINDX=IDX3
                 
                    endif

                    INE=HOC(KINDX,LINDX,MINDX)

                    DO WHILE (INE.NE.0)

                    DX=xp-X(INE,1)
                    DY=yp-X(INE,2)
                    DZ=zp-X(INE,3)
                                        
                    IF ( iperiodic ) then
                    call modbound(dx,dy,dz)
                    endif
                    R=dsqrt(DX**2+DY**2+DZ**2)

                    select case ( iopt )

                      case (1)
             
                      rsupport=2*HP
           
                      case (2)

                  rsupport=2*DMAX1(HP,H(INE)) 
                       
 !                    rsupport=hp+h(ine)
             
                      end select 
                              
                      IF (R.LT.rsupport) THEN
                      
!    (a new neighbour has been found)

                       NBN=NBN+1
                       NN(nbn)=INE
                    ENDIF
                    INE=LL(INE)
                 ENDDO
!    (done with that cell)
              ENDDO
           ENDDO
        ENDDO

          return

          end
          
         SUBROUTINE LLN2(xp,yp,zp,iopt,hp)

        USE wp3d_h

        INTEGER::iopt,k1,k2,l1,l2,m1,m2,LINDX,KINDX, &
        MINDX,idx1,idx2,idx3,modindx
        INTEGER::kne,lne,ine,mne

        double precision::hp,r,dx,dy,dz,rsupport
        double precision::xp,yp,zp

! Calculate the neighbour list of a particle

          nbn=0
                  
         K1=floor((XP-nhp*HP-XMIN)/HH)          
         K2=floor((XP+nhp*HP-XMIN)/HH)
         L1=floor((YP-nhp*HP-YMIN)/HH)
         L2=floor((YP+nhp*HP-YMIN)/HH)
         M1=floor((ZP-nhp*HP-ZMIN)/HH)
         M2=floor((ZP+nhp*HP-ZMIN)/HH)
 
! Find indices of cells containing potential neighbours:
          
!     Look for neighbours in all these cells:

              DO MNE=M1,M2
              DO LNE=L1,L2
              DO KNE=K1,K2
         
                idx1=kne
                idx2=lne
                idx3=mne
                 
                 if ( iperiodic ) then

                KINDX=modindx(idx1,ncx)
                LINDX=modindx(idx2,ncy)
                MINDX=modindx(idx3,ncz)
                   
                 else
                 
                 KINDX=IDX1
                 LINDX=IDX2
                 MINDX=IDX3
                 
                    endif

                 INE=HOC(KINDX,LINDX,MINDX)

                    DO WHILE (INE.NE.0)
                 
                    DX=xp-X(INE,1)
                    DY=yp-X(INE,2)
                    DZ=zp-X(INE,3)

                    if ( iperiodic ) then
                    call modbound(dx,dy,dz)
                    endif

                    R=dsqrt(DX**2+DY**2+DZ**2)

                    select case ( iopt )

                      case (1)
             
                      rsupport=2*HP
           
                      case (2)

                      rsupport=2*DMAX1(HP,H(INE)) 
             
                      end select 
                              
                      IF (R.LT.rsupport) THEN
                      
!    (a new neighbour has been found)

                       NBN=NBN+1
                       NN(nbn)=INE
                    ENDIF
                    INE=LL(INE)
                 ENDDO
!    (done with that cell)
              ENDDO
           ENDDO
        ENDDO

        return

        end   
    
            
      


 
 
        



!
