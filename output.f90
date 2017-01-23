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
      integer::ihpmax,ihpmin,nnmax,nnmin
      
!     Find system box 
      XLOW=1.E30
      YLOW=1.E30
      ZLOW=1.E30
      XHIGH=-1.E30
      YHIGH=-1.E30
      ZHIGH=-1.E30
      tiny=1.0d-10

      do i=1,n
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
          EINT=EINT+u(i,5)-0.5d0*vol(i)*(wprim(i,1)*(wprim(i,2)**2+ &
          wprim(i,3)**2+wprim(i,4)**2)+wprim(i,6)**2+wprim(i,7)**2+ &
          wprim(i,8)**2)
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
        RHOMAX=0.
        irhomin=1
        irhomax=1
        HPMAX=0.
        HPMIN=1.E30
        ihpmin=1
        ihpmax=1
    
        ebmax=0.0d0
        ebmin=1.0d25
        ebave=0.0d0
        iebmax=1
        iebmin=1 
      
        do i=1,n 
  
        eb(i)=dabs(divb(i))*h(i)/dsqrt(wprim(i,6)**2+wprim(i,7)**2+ &
        wprim(i,8)**2+tiny)/vol(i)

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

        IF (eb(i).LT.EBMIN) THEN
        IEBMIN=I
        EBMIN=eb(i)
        ENDIF
        IF (eb(i).GT.EBMAX) THEN
        IEBMAX=I
        EBMAX=eb(i)
        ENDIF

        ebave=ebave+eb(i)

      ENDDO

      ebave=ebave/float(n)
      
      WRITE (2,fmt="(A,I8,A,ES12.5)") 'OUTPUT: end of iteration ',NIT,' time= ',t
      WRITE (2,*) 'System box:'
      WRITE (2,fmt="(6ES12.5)") XLOW,XHIGH,YLOW,YHIGH,ZLOW,ZHIGH
      WRITE (2,fmt="(A,4ES12.5)") 'rhomin:',RHOMIN,X(IRHOMIN,1),X(IRHOMIN,2),X(IRHOMIN,3)
      WRITE (2,fmt="(A,4ES12.5)") 'rhomax:',RHOMAX,X(IRHOMAX,1),X(IRHOMAX,2),X(IRHOMAX,3)
      WRITE (2,fmt="(A,4ES12.5)") 'hmin:',HPMIN,X(IHPMIN,1),X(IHPMIN,2),X(IHPMIN,3)
      WRITE (2,fmt="(A,4ES12.5)") 'hmax:',HPMAX,X(IHPMAX,1),X(IHPMAX,2),X(IHPMAX,3)
      WRITE (2,fmt="(A,4ES12.5)") 'ebmin:',EBMIN,X(IEBMIN,1),X(IEBMIN,2),X(IEBMIN,3)
      WRITE (2,fmt="(A,4ES12.5)") 'ebmax:',EBMAX,X(IEBMAX,1),X(IEBMAX,2),X(IEBMAX,3)
 
      if ( igrav ) then
      write(3,fmt='(12(ES12.5))') t,eint,ekin,emag,epot,jtot,ptot,rhomax,mtot,ebmin,ebmax,ebave
      else
      write(3,fmt='(11(ES12.5))') t,eint,ekin,emag,jtot,ptot,rhomax,mtot,ebmin,ebmax,ebave
      endif
                
      RETURN
      END
      
		
        SUBROUTINE PDUMP(wprim) 

        use wp3d_h 

        integer::i
        character(len=40)::filename
        double precision::wprim(nmax,nv),rmax

        ndump=ndump+1
        
        write(6,fmt="(A,ES12.5)") 'writing dump file at t=',t
        write(2,fmt="(A,ES12.5)") 'writing dump file at t=',t
        
        write(unit=filename, fmt="(A,I0.3,A)") trim(outfn(3)), &
        ndump,".dat"

        nacc=0
        nsink=0

        open(12,FILE=filename)
        
        write(12,fmt="(3I8)") n,nacc,nsink
        write(12,fmt="(2ES25.15)") t,gam(1)
     
        do i=1,n 
        
        write(12,fmt="(14ES25.15)") x(i,1),x(i,2),x(i,3),wprim(i,2), &
        wprim(i,3),wprim(i,4),wprim(i,6),wprim(i,7),wprim(i,8),am(i),h(i), &
        wprim(i,5),wprim(i,1),eb(i)
       
! SPH particle dump 

        enddo
        
        close(12)
        
        return 
        END 


        SUBROUTINE restartdump(wprim,psi) 

        use wp3d_h 

        integer::i
        character(len=40)::filename
        double precision::wprim(nmax,nv),psi(nmax),rmax

        write(6,fmt="(A,ES12.5)") 'writing restart dump file at t=',t
        write(2,fmt="(A,ES12.5)") 'writing restart dump file at t=',t
        
        filename="restartdump.dat"
        
        open(12,FILE=filename)
     
        do i=1,n

        if ( iDedner ) then
        write(12,fmt="(14ES25.15)") x(i,1),x(i,2),x(i,3),wprim(i,2), &
        wprim(i,3),wprim(i,4),wprim(i,6),wprim(i,7),wprim(i,8),am(i),h(i), &
        wprim(i,5),wprim(i,1),psi(i)
        else
        write(12,fmt="(13ES25.15)") x(i,1),x(i,2),x(i,3),wprim(i,2), &
        wprim(i,3),wprim(i,4),wprim(i,6),wprim(i,7),wprim(i,8),am(i),h(i), &
        wprim(i,5),wprim(i,1)
        endif
       
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
        
        character(len=20):: filename
        character(len=10),allocatable::names(:)
        
        iopt=1
        tiny=1.0d-10
        ndump=ndump+1

        allocate(rhom(nx1,ny1,nz1))
        allocate(velm(nx1,ny1,nz1,ndim))
        allocate(ptherm(nx1,ny1,nz1))
        allocate(bm(nx1,ny1,nz1,ndim))
        allocate(divbm(nx1,ny1,nz1))
        allocate(xgrid(nx1),ygrid(ny1),zgrid(nz1))

        if ( setupchoice .eq. 2 ) then
                
        dx1=boxx/float(nx1)
        dy1=boxy/float(ny1)
        dz1=boxz/float(nz1)
                
        endif
                
        write(unit=filename, fmt="(A,I0.3,A)") trim(outfn(3)), &
        ndump,".vtr"
     
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

         bm(i+1,j+1,k+1,1)=wprim(idx,6)
         bm(i+1,j+1,k+1,2)=wprim(idx,7)
         bm(i+1,j+1,k+1,3)=wprim(idx,8)

         divbm(i+1,j+1,k+1)=dabs(divb(idx))

         enddo
         enddo
         enddo
                  
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
                
        call vtk_3d(nx1,ny1,nz1,xgrid,ygrid,zgrid,nscal,nvec, &
        scaldata,vecdata,names,filename)
        
        deallocate(rhom,velm,ptherm,xgrid,ygrid,zgrid)
        deallocate(vecdata,scaldata,bm,divbm)
                          
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

        fname1=trim(outfn(3)) // ndx3 // ".vtr"
        
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

        hsmooth=dmax1(hsmooth,h(nn(in)),dx/2.0d0)
        
        enddo

        h2=hsmooth**2
        h3=hsmooth*h2
        
        do in=1,nbn
          
        r2=(x(nn(in),1)-xtry)**2+(x(nn(in),2)-ytry)**2+ &
		(x(nn(in),3)-ztry)**2 
        
        if ( r2 .le. 4*h2 ) then 
                                         
!     We calculate the density contribution of all particles to all
!     the grid points with which their kernel function overlaps
        
        u2=r2/h2

        ITAB=INT(CTAB*U2)+1

        rhom(i,j,k)=rhom(i,j,k)+am(nn(in))*wtab(itab)/h3
     
        pm(i,j,k)=pm(i,j,k)+am(nn(in))*(wprim(nn(in),5)/ &
		wprim(nn(in),1))*wtab(itab)/h3

        do dim=1,ndim

        velm(i,j,k,dim)=velm(i,j,k,dim)+am(nn(in))*(wprim(nn(in),dim+1)/ &
		wprim(nn(in),1))*wtab(itab)/h3

        bm(i,j,k,dim)=bm(i,j,k,dim)+am(nn(in))*(wprim(nn(in),dim+5)/ &
		wprim(nn(in),1))*wtab(itab)/h3

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
        write(unit=filename, fmt="(A,I0.3,A)") trim(outfn(3)), &
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
 
!       err = dbputpm (dbfile, "pointmesh", 9, ndims, xp,yp,zp, &
!           n, DB_DOUBLE, DB_F77NULL, ierr)
          
!       err = dbputpv1(dbfile, "rho", 3, "pointmesh", 9, &
!      rho, n, DB_DOUBLE, DB_F77NULL, ierr)
!      err = dbputpv1(dbfile, "pth", 3, "pointmesh", 9, &
!       pth, n, DB_DOUBLE, DB_F77NULL, ierr)
 
!      if ( imhd ) then
!       err = dbputpv1(dbfile, "eb", 2, "pointmesh", 9, &
!      eb, n, DB_DOUBLE, DB_F77NULL, ierr)
!      endif

!      err = dbputpv1(dbfile, "vx", 2, "pointmesh", 9, &
!       vx, n, DB_DOUBLE, DB_F77NULL, ierr)
!      err = dbputpv1(dbfile, "vy", 2, "pointmesh", 9, &
!      vy, n, DB_DOUBLE, DB_F77NULL, ierr)
!      err = dbputpv1(dbfile, "vz", 2, "pointmesh", 9, &
 !     vz, n, DB_DOUBLE, DB_F77NULL, ierr)
           
!       if ( imhd ) then
 !      err = dbputpv1(dbfile, "bx", 2, "pointmesh", 9, &
!       bx, n, DB_DOUBLE, DB_F77NULL, ierr)
!       err = dbputpv1(dbfile, "by", 2, "pointmesh", 9, &
!       by, n, DB_DOUBLE, DB_F77NULL, ierr)
!       err = dbputpv1(dbfile, "bz", 2, "pointmesh", 9, &
 !      bz, n, DB_DOUBLE, DB_F77NULL, ierr)
!       endif
   
!       deallocate(rho,pth,vx,vy,vz,bx,by,bz,xp,yp,zp)
           
 !      ierr=dbclose(dbfile)
     
       return
       end 