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
                 
                 
         subroutine HLLMHD(i,j,uleft,uright,nij,F,ax,bxc,psi,psic,psiflux)
                
         use wp3d_h

         implicit none

! HLL Riemann solver for the system of 1D MHD equations

         double precision:: uleft(nv),uright(nv), &
         uLL(nv),ax,F(nv),Fv(ndim),Fvr(ndim),nij(ndim), &
         Fl(nv),Fr(nv),al,ar,bl,br,valfvenxl,valfvenxr
         double precision::vxl,vyl,vzl,el,rhol,vxr,vyr,vzr,er,rhor, &
         pr,pl,cl,cr,Sl,Sr,bxl,bxr,byl,bzl,byr,bzr,plt,prt,cmax,psiflux
         double precision:: bxc,psil,psir,psic,pbar,psi(nmax)
         integer::i,j

         rhol=uleft(1)
         vxl=uleft(2)/rhol
         vyl=uleft(3)/rhol
         vzl=uleft(4)/rhol
         el=uleft(5)
         psil=psi(i)/am(i)
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
		 
         psir=psi(j)/am(j)
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
         
         if ( iDedner ) then
         if ( F(1) .gt. 0.0d0 ) then
         psiflux=F(1)*psil
         else
         psiflux=F(1)*psir
         endif
         endif
         
         return

         end
                 
         
         subroutine HLLDMHD(i,j,uleft,uright,nij,F,ax,bxc,psi,psic,psiflux)
        
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
         double precision:: bxc,psil,psir,psic,pbar,psiflux
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
         call HLLMHD(i,j,uleft,uright,nij,F,ax,bxc,psi,psic,psiflux)
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
         
         if ( iDedner ) then
         if ( F(1) .gt. 0.0d0 ) then
         psiflux=F(1)*psil
         else
         psiflux=F(1)*psir
         endif
         endif

         return
         
         end