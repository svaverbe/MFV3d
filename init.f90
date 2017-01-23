        SUBROUTINE init(u,wprim)
     
        use wp3d_h

        implicit none
         
        integer:: idx,nt,idum1,idum2,idum3
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
                
        nx1=50
        ny1=50
        nz1=50
        d0=1.0d0
        nt=nx1*ny1*nz1
        n=nt
        tf=0.1d0
        tprin=0.001d0   
           
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
        wprim(idx,6)=1.0d0/dsqrt(2.0d0)
        wprim(idx,7)=1.0d0/dsqrt(2.0d0)
        wprim(idx,8)=1.0d0/dsqrt(2.0d0)
!       wprim(idx,6:8)=0.0d0
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

        wprim(idx,6)=0.0d0
        wprim(idx,7)=0.0d0
        wprim(idx,8)=0.0d0
      
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
        wprim(idx,6)=0.0d0
        wprim(idx,7)=0.0d0
        wprim(idx,8)=0.0d0
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

        wprim(idx,6)=0.0d0
        wprim(idx,7)=0.0d0
        wprim(idx,8)=0.0d0
      
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
      parameter(nint=100,ni=125000)
      double precision::rho0,rcloud,r,xtest,ytest,ztest
      double precision::gravity,boltzmann,protonmass,AU
      double precision::msun,unitmass,unitvelocity,unittime
      double precision::unitsurfacedensity,pbar
      double precision::unitlength,unitdensity,unitenergy  
      double precision::unitpressure,unitspecificenergy
      double precision::tff,temp,radius,mass,om,mcloud
      double precision::meanweight,Prot,mcrit,mu,c1, &
      unitBfield,mu0,dens,hini,parsec,rotratio,thetaoblique, &
      bx,bz,bmag
      double precision::wprim(nmax,nv)

       gravity=6.672d-8
       boltzmann=1.3806d-16
       protonmass=1.6726d-24
       meanweight=2.33d0*protonmass
       mu0=1.0d0
       gam(:)=1.0d0
       eostype=2
       
       rhocrit1=1.0d-13
       rhocrit2=5.7d-5
       rhocrit3=1.0d0

       fname="magcol"
           
       xmin=-3.0d0
       xmax=3.0d0
       ymin=-3.0d0
       ymax=3.0d0
       zmin=-3.0d0
       zmax=3.0d0

       iperiodic=.true. 
       periodictypex=1
       periodictypey=1
       periodictypez=1
           
       boxx=xmax-xmin
       boxy=ymax-ymin
       boxz=zmax-zmin
       
!      nmedium=50
       nmedium=40
       nt=nmedium**3
       voltot=boxx*boxy*boxz
       hini=0.05d0
       parsec=3.08567758d+18
       tf=1.5d0
       tprin=0.01d0
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
	   
       mu=20
       c1=0.53d0
       bmag=3*mass/(2*pi*mu*c1*radius**2)* &
       (pi*gravity*mu0/5)**(1.0d0/2.0d0)
       bmag=bmag/unitbfield
!      bmag=0.0d0
!      thetaoblique=pi/4
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
       write(1,*) 'initial magnetic field strength(µG):',bmag*unitBfield*1.0d+6

! initial homogeneous density of the cloud 

       call setupem1(ni,ncloud)
      
       rhocrit4=(pi**3*(c0*unitvelocity)**6/(4*gravity**3*nopt**2))* &
       (float(ncloud)/mass)**2/unitdensity
     
       mcrit=pi**(3./2.)*(c0*unitvelocity)**3/(2*nopt*gravity**(3./2.)* &
       dsqrt(rhocrit1*unitdensity))/unitmass
        
       write(1,*) 'critical density( g cm-3)',rhocrit4*unitdensity
       write(1,*) 'critical particle mass',mcrit

       dx1=boxx/float(nmedium)
       dy1=dx1
       dz1=dx1
           
       exterior(:)=.false. 

        n=ncloud
        
        bx=bmag*dsin(thetaoblique)/dsqrt(4*pi)
        bz=bmag*dcos(thetaoblique)/dsqrt(4*pi)

        do i=1,ncloud
                
!       wprim(i,1)=dens(rho0,x(i,1),x(i,2))
        wprim(i,1)=rho0
        wprim(i,2)=-2*pi*x(i,2)/Prot
        wprim(i,3)=2*pi*x(i,1)/Prot
        wprim(i,4)=0.0d0

        wprim(i,2:4)=0.0d0
                   
        wprim(i,5)=pbar(i,wprim(i,1))
        
        gam(i)=5.0d0/3.0d0
	
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
       
       gam(n)=5.0d0/3.0d0
          
       wprim(n,2:4)=0.0d0 

       wprim(n,5)=pbar(n,wprim(n,1))        
          
       wprim(n,6)=bx
       wprim(n,7)=0.0d0
       wprim(n,8)=bz
           
       h(n)=hini

      endif
          
      enddo
      enddo
      enddo     

      do i=1,n
      
      if ( barotropic ) then
      cs(i)=dsqrt(wprim(i,5)/wprim(i,1))
      else
      cs(i)=dsqrt(gam(i)*wprim(i,5)/wprim(i,1))
      endif
               
      enddo
               
      write(1,*) 'number of particles inside the cloud:',ncloud
  
      return
          
      end
      
    
      SUBROUTINE setupem1(ni,ncloud)  

      use wp3d_h  
          
      IMPLICIT NONE

! sets up a particle distribution in 3D 
! particles are placed on a stretched grid 
                        
      double precision::avec(3),bvec(3),cvec(3)
      double precision::space,xp,yp,zp
      double precision::r,rns
      INTEGER::maxn,i,idum,k,l,m,ni,ncloud
              
! setup particles on a hexagonal lattice

      rns=1.0d0
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
!      write (6,*)' maxn: ',maxn,' spacing: ',space
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
!      write (6,*)'number of particles within rns: ',i
       
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
        hini=0.5d0
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
         tf=0.01d0
         tprin=0.01d0

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
	   
!       write (6,*)'number of particles within rns: ',i
	   
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
 
!       if ( relax .eqv. .true. ) then
!       call relaxparticles(niter,wprim)
!       endif
                
! relax the particle distribution
                
       end select 
                
       return
                
       end
                
                
       subroutine relaxparticles(niter,wprim)
                
! subroutine to relax a random particle distribution 
                
       USE wp3d_h
                
       implicit none
                
       integer::it,i,j,in,niter,iopt,itab
       double precision::wij,dri(nmax,ndim),deltar, &
       h3i,h2i,d2,u2,av,dx,dy,dz,periodic1,periodic2, &
       wprim(nmax,nv)
           
       iopt=1
       av=0.0001d0
                
           do it=1,niter
                
           call build_linked_list
           call hvol(wprim)
            
           deltar=0.0d0
                
           do i=1,n

           dri(i,:)=0.0d0
                   
           call LLN(i,iopt,h(i))
             
           H2I=H(I)**2     
           H3I=H2I*h(i)
        
!  Compute dR for each particle

           DO IN=1,NBN
         
            J=NN(IN)

            if ( iperiodic .eqv. .true. ) then
                   
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
	  
	  

	  
	  
