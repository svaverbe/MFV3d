SUBROUTINE PARARANGE(n1,n2,nprocs)
     
USE wp3d_h
    
!  determine workload for each processor and allocate dynamic storage arrays

integer*4::irank,n1,n2,nprocs,iwork1,iwork
integer*4::ista,iend
	
iwork1=(n2-n1+1)/nprocs
iwork2=mod(n2-n1+1,nprocs)
	
 do irank=0,nprocs-1 
	
 ista=irank*iwork1+n1+min(irank,iwork2)
 iend=ista+iwork1-1
	
 if ( iwork2 .gt. irank ) then 
 iend=iend+1
 endif 
	
 ilen1(irank)=iend-ista+1
 idisp1(irank)=ista-1
 ilen2(irank)=ndim*ilen1(irank)
 idisp2(irank)=ndim*(ista-1)
 ilen3(irank)=ndim**2*ilen1(irank)
 idisp3(irank)=ndim**2*(ista-1)
 ilen4(irank)=nv*ilen1(irank)
 idisp4(irank)=nv*(ista-1)
 ilen5(irank)=ndim*nv*ilen1(irank)
 idisp5(irank)=ndim*nv*(ista-1)
 
 if ( irank .eq. myrank ) then
	
 n_lower=ista
 n_upper=iend

! range and the number of particles within this range

 allocate(myh(ista:iend),myvol(ista:iend),myam(ista:iend), &
 mydivb(ista:iend),mypsidot(ista:iend), &
 myVudot(1:nv,ista:iend),mylim(1:nv,ista:iend), &
 mygrad(1:ndim,1:nv,ista:iend), &
 myA(1:ndim,1:ndim,ista:iend),mydti(ista:iend),mytau(ista:iend), &
 myvsigmax(ista:iend),myconditionnumber(ista:iend))
 
 if ( igrav ) then
 allocate(mydzeta(ista:iend),myphi(ista:iend))
 endif

 endif 
 
 enddo
        
! allocate dynamic arrays for time stepping and local particle storage on 
! each processor
 		 		     
 return 
 end      