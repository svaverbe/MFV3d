      subroutine build_linked_list
      
      USE wp3d_h

! construction of a Linked list for neighbour searching on a single 3D grid
         
      LOGICAL::EXC
      double precision::hpav,h1,h2,h3
    
      INTEGER::i,k,l,m

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

! Readjust if maximum number allowed is exceeded:
! (Note that this version assumes centering on the origin!)
      EXC=.FALSE.
      IF (NCX.GT.NCMX) THEN 
      EXC=.TRUE.
      NCX=NCMX
      ENDIF
      IF (NCY.GT.NCMY) THEN 
      EXC=.TRUE.
      NCY=NCMY
      ENDIF
      IF (NCZ.GT.NCMZ) THEN
      EXC=.TRUE.
      NCZ=NCMZ
      ENDIF

! Initialize Head-Of-Cell array:
      DO M=0,NCMZ
        DO L=0,NCMY
          DO K=0,NCMX
            HOC(K,L,M)=0
          ENDDO
        ENDDO
      ENDDO
      
      IF (EXC) THEN
      H1=dabs((xmax-xmin)/dfloat(ncx-1))
      H2=dabs((ymax-ymin)/dfloat(ncy-1))
      H3=dabs((zmax-zmin)/dfloat(ncz-1))
      HH=DMAX1(H1,H2,H3)
      ENDIF

! Build linked lists (cf. Hockney and Eastwood):
      
!   (no need to check for out of bounds indices)
         DO I=1,N
         if ( iperiodic ) then
         K=floor((X(I,1)-XMIN)/HH)
         L=floor((X(I,2)-YMIN)/HH)
         M=floor((X(I,3)-ZMIN)/HH)
         else
         K=INT((X(I,1)-XMIN)/HH)+1
         L=INT((X(I,2)-YMIN)/HH)+1
         M=INT((X(I,3)-ZMIN)/HH)+1
         endif
         LL(I)=HOC(K,L,M)
         HOC(K,L,M)=I
         ENDDO
     
      RETURN
      END  subroutine build_linked_list
	  
	  
      SUBROUTINE LLN(i,iopt,hp)

      USE wp3d_h
      
      implicit none

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

      if ( iperiodic .eqv. .true. ) then
      K1=floor((XP-nhp*HP-XMIN)/HH)          
      K2=floor((XP+nhp*HP-XMIN)/HH)
      L1=floor((YP-nhp*HP-YMIN)/HH)
      L2=floor((YP+nhp*HP-YMIN)/HH)
      M1=floor((ZP-nhp*HP-ZMIN)/HH)
      M2=floor((ZP+nhp*HP-ZMIN)/HH)
      else
      K1=MAX(1,INT((xp-nhp*HP-XMIN)/HH)+1)
      K2=MIN(NCMX,INT((xp+nhp*HP-XMIN)/HH)+1)
      L1=MAX(1,INT((yp-nhp*HP-YMIN)/HH)+1)
      L2=MIN(NCMY,INT((yp+nhp*HP-YMIN)/HH)+1)
      M1=MAX(1,INT((zp-nhp*HP-ZMIN)/HH)+1)
      M2=MIN(NCMZ,INT((zp+nhp*HP-ZMIN)/HH)+1)
      endif

!     Look for neighbours in all these cells:

      DO MNE=M1,M2
      DO LNE=L1,L2
      DO KNE=K1,K2
           
      idx1=kne
      idx2=lne
      idx3=mne
                 
      if ( iperiodic .eqv. .true. ) then
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
                                        
      if ( iperiodic .eqv. .true. ) then
      call modbound(dx,dy,dz)
      endif
      R=dsqrt(DX**2+DY**2+DZ**2)

      select case ( iopt )
      case (1)
      rsupport=2*HP
      case (2)
      rsupport=2*DMAX1(HP,H(INE))   
      end select 
                              
      if (R.LT.rsupport) THEN
                      
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
      
      implicit none

      INTEGER::iopt,k1,k2,l1,l2,m1,m2,LINDX,KINDX, &
      MINDX,idx1,idx2,idx3,modindx
      INTEGER::kne,lne,ine,mne

      double precision::hp,r,dx,dy,dz,rsupport
      double precision::xp,yp,zp

! Calculate the neighbour list of a particle

      nbn=0
      
      if ( iperiodic ) then       
      K1=floor((XP-nhp*HP-XMIN)/HH)          
      K2=floor((XP+nhp*HP-XMIN)/HH)
      L1=floor((YP-nhp*HP-YMIN)/HH)
      L2=floor((YP+nhp*HP-YMIN)/HH)
      M1=floor((ZP-nhp*HP-ZMIN)/HH)
      M2=floor((ZP+nhp*HP-ZMIN)/HH)
      else
      K1=MAX(1,INT((xp-nhp*HP-XMIN)/HH)+1)
      K2=MIN(NCMX,INT((xp+nhp*HP-XMIN)/HH)+1)
      L1=MAX(1,INT((yp-nhp*HP-YMIN)/HH)+1)
      L2=MIN(NCMY,INT((yp+nhp*HP-YMIN)/HH)+1)
      M1=MAX(1,INT((zp-nhp*HP-ZMIN)/HH)+1)
      M2=MIN(NCMZ,INT((zp+nhp*HP-ZMIN)/HH)+1)
      endif
 
! Find indices of cells containing potential neighbours:
          
!     Look for neighbours in all these cells:

      DO MNE=M1,M2
      DO LNE=L1,L2
      DO KNE=K1,K2
         
      idx1=kne
      idx2=lne
      idx3=mne
                 
      if ( iperiodic .eqv. .true. ) then
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

      if ( iperiodic .eqv. .true. ) then
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
