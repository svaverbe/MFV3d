      SUBROUTINE MKTREE(wprim)
      
! --------------------------------------------------------------
! MKTREE: initialize the tree structure for the force calculation.
! ----------------------------------------------------------------
      
      use wp3d_h

      implicit none

      integer::I, MKCELL, K, P, IND(1:mxnode),  &
      Q, J, L, M1,M2 
      DOUBLE PRECISION::XYZMAX, POS0(NDIM), DIST2,DR, &
      wprim(nmax,nv)
      
      mxcell=n
      incell=n+1

!     Expand root volume to enclose all particles.
!     --------------------------------------           
      
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

!     ------------------------------------------
!     Initialize midpoint and size of root cell.
!     ------------------------------------------

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