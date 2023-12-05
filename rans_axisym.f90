!------------------------------------------------------------------------------------------
!-----
!----- SIMULATION OF AXISYMMETRIC JET FLOW 
!----- RANS EQUATIONS ARE DISCRITIZED USING THE FINITE-VOLUME METHOD
!----- THE SIMPLE PRESSURE CORRECTION ALGORITHM IS USED 
!----- STAGGERED GRID IS EMPLOYED FOR THE U,V VELOCITIES
!----- NO SWIRLING FLOW IS CONSIDERED (W=0)
!-----
!-------------------------------------------------------------------------------------------
  include "module_rans.f90"
  program CFD_AXIS_JET   

  use rans  
      
  double precision :: errmaxp,errmaxu,errmaxv,errmaxte
  double precision :: Qinflow,Qout,Qside 
      
  pi = 4.d0*datan(1.0d0)
  
!---READ INPUT DATA                           
      
  open(5,file='input.txt') 
  read(5,*)                                                     !   xmin,xmax: start,end of domain in x-axis                               
  read(5,*)  xmin, xmax, xuni, ymin, change_y_pos, ymax         !   ymin,ymax: start,end of domain in y-axis
  write(*,*) xmin, xmax, xuni, ymin, change_y_pos, ymax         !   xuni: length of uniform grid in x-direction
  read(5,*)                                                     !   change_y_pos: position in the vertical axis where we reverse the geometric progression to be denser at the wall                                  
  read(5,*)  ngridx2, ngrid_uni, ngridy1, ngridy2               !   ngridx1: number of nodes in uniform part in x-axis
  write(*,*) ngridx2, ngrid_uni, ngridy1, ngridy2               !   ngrid_uni: number of nodes in uniform part in y-axis  
  read(5,*)                                                     !   ngridy1: nodes in the ratio 1 part in the vertical axis                                                 
  read(5,*) ratx1,ratx2                           !   Ratio of geometrical progression in x direction: xuni-->xmax
  write(*,*) ratx1,ratx2                          !   Ratio of geometrical progression in y direction: 0-->1    
  read(5,*)                                       !                                              
  read(5,*) dvisc, Uinf                           !   Kinematic viscosity,  Free stream velocity 
  write(*,*) dvisc, Uinf
  read(5,*)                                       !   
  read(5,*) radius, ct                            !   Rotor radius, Thrust coefficient
  write(*,*) radius, ct    
  read(5,*)                                          
  read(5,*) itmax, nswp, nbackup                  !   itmax: number of maximum iterations, nswp: number of ADI sweeps 
  write(*,*) itmax, nswp, nbackup                 !   nbackup: number of iterations to store flow field in backup file
  read(5,*)                                       !                              
  read(5,*) eps, tiny                             !   eps: Convergence criterion 
  write(*,*) eps, tiny                            !   tiny: Small number
  read(5,*) 
  read(5,*) urfu,urfv,urfp,urfte,urfvis           !   urfu,urfv,urfp,urfte,urfvis: underrelaxation factors for u,v,p,k,vt
  write(*,*) urfu,urfv,urfp,urfte,urfvis    
  read(5,*)
  read(5,*) tiamb                                 !   tiamb: Ambient turbulence intensity
  write(*,*) tiamb              
  read(5,*)
  read(5,*) ibackup                               !   ibackup=0 (no backup file is written)
  write(*,*) ibackup                              !   ibackup=1 (write backup file)

!--- Grid definition
  call GRID_DEF 

  allocate (UVEL(ngridx,ngridy-1))
  allocate (VVEL(ngridx,ngridy))
  allocate (PTOT(ngridx-1,ngridy-1))
  allocate (PRE(ngridx-1,ngridy-1))
  allocate (TE(ngridx-1,ngridy-1))
  allocate (vt(ngridx-1,ngridy-1))
  allocate (UOLD(ngridx,ngridy-1))
  allocate (VOLD(ngridx,ngridy))
  allocate (TEOLD(ngridx-1,ngridy-1))
  allocate (vtold(ngridx-1,ngridy-1))
  allocate (w(ngridx),akl(ngridx))

!---- FLOW FIELD INITIALIZATION
  
  do i=1,ngridx
   do j=1,ngridy-1
     UVEL(i,j)=1.d0
   enddo
  enddo

  VVEL=0.d0
  PTOT=0.d0
  PRE=0.0d0
  TE=0.d0
  vt=0.d0
  vtold=0.d0

!----- Calculate Reynolds number                
  
  Re=Uinf*radius/dvisc
  write(*,*) 'Re=', Re

!----- Number of starting iteration
  itstart=1

!---- Read from backup files for rerun 

  if(ibackup.eq.1) then
   open(97,file='Uvelocity')
   open(98,file='Vvelocity')
   open(99,file='Pres-TKE')

    do i=1,ngridx
     do j=1,ngridy-1
      read(97,*) xx, yy, UVEL(i,j)
     enddo
     read(97,*)
    enddo

    do i=1,ngridx-1
     do j=1,ngridy
      read(98,*) xx, yy, VVEL(i,j)
     enddo
     read(98,*)
    enddo

    read(99,*) itstart
    do i=1,ngridx-1
     do j=1,ngridy-1
      read(99,*) xx,yy,PTOT(i,j),TE(i,j),vt(i,j)
     enddo
     read(99,*)
    enddo
   close(97)
   close(98)
   close(99)
  endif

!----  Open file for monitoring errors 

 open(200,file='errormax')

!----  Start basic loop of iterations
  
 do it=itstart,itmax    

!----  Calculate inflow/outflow mass

    Qinflow=0.
    do j=1,ngridy-1
     yc = y_grid_stgu(j)+0.5d0*dygridstgu(j)
     Qinflow = Qinflow + UVEL(1,j)*dygridstgu(j)*2*pi*yc
    enddo

    Qout=0.
    do j=1,ngridy-1
     yc = y_grid_stgu(j)+ 0.5*dygridstgu(j)
     Qout = Qout + UVEL(ngridx,j)*dygridstgu(j)*2*pi*yc
    enddo

    Qup=0.
    do i=1,ngridx
     yc = y_grid(ngridy)
     Qup = Qup + VVEL(i,ngridy)*dxgridstgv(i)*2*pi*yc
    enddo

!---- Calculate mass flow difference and correction factor for the outflow    

    Qcor=-Qinflow+Qout+Qup
    acor=Qcor/Qout
    write(*,*) 'Qinflow=', Qinflow
    write(*,*) 'Qout=', Qout    
    write(*,*) 'Qup=', Qup    
    write(*,*) 'Qcor=',Qcor

!---- Algebraic turbulence model for calculation of turbulent viscosity 
   if(it.eq.1) then
    call TURBVIS       
   endif
!----Momentume equation u    
   call MOMENTUM_U(errmaxu)
 
!--- Momentum equation v
   call MOMENTUM_V(errmaxv)
  
!--- Presure correction      
   call PRESCOR(errmaxp)
   
!--- Velocity correction
   call VELOCOR

   call CALCTE(errmaxte)

   call TURBVIS       

!--- Underrelaxation of pressure
!--- Define a pressure reference point
   Pref = PTOT(1,ngridy-1)
   PTOT=PTOT+urfp*PRE-Pref
   PRE=0.

!--- Write velocty and pressure field every nbackuo iterations
   
   if (mod(it,nbackup).eq.0) then
    open(97,file='Uvelocity')
    open(98,file='Vvelocity')
    open(99,file='Pres-TKE')

    do i=1,ngridx
     do j=1,ngridy-1
      write(97,'(3e18.9)') x_grid(i), 0.5*(y_grid(j)+y_grid(j+1)) , UVEL(i,j)
     enddo
     write(97,*)
    enddo

    do i=1,ngridx-1
     do j=1,ngridy
      write(98,'(3e18.9)') 0.5*(x_grid(i)+x_grid(i+1)), y_grid(j) , VVEL(i,j)
     enddo
     write(98,*)
    enddo  

    write(99,*) it
    do i=1,ngridx-1
     do j=1,ngridy-1
      write(99,'(5e18.9)') x_grid(i),y_grid(j), PTOT(i,j),TE(i,j),vt(i,j) 
     enddo
     write(99,*)
    enddo
    close(97)
    close(98)
    close(99)    
   end if 

  write(*,'(i7,3e18.9)') it,errmaxu,errmaxv,errmaxp
  write(200,'(i7,3e18.9)') it,errmaxu,errmaxv,errmaxp
   
!  write residuals
   
  if(errmaxp.lt.eps.and.errmaxu.lt.eps.and.errmaxv.lt.eps) exit
                           
  enddo
  
   
  END program CFD_AXIS_JET   
      
!--------------------------------------------------------------------------------
  Subroutine GRID_DEF
!--------------------------------------------------------------------------------
  use rans  

  double precision :: dy_grid,dx_grid1
  integer           :: i,j
 

!---    DEFINITION OF x_grid, y_grid, + staggered
      
      
  open(11,file='grid1')
  open(12,file='grid2')
  open(13,file='stggridU1')
  open(14,file='stggridU2')
  open(15,file='stggridV1')
  open(16,file='stggridV2') 

  ! 3 Areas. min -> grid_uni, uniform -> ngridy1, ngridy1 -> ngridy2
  ! In 2nd and 3rd, we give init point, last point and number of points in between and it finds the 
  ! appropriate ratio to do the job.

  ! Basic calculations
  dy_grid = 1.0d0 / real(ngrid_uni - 1)
  ngridy = ngrid_uni + ngridy1 + ngridy2 - 2

  ! Allocations
  allocate(y_grid(ngridy))

  !--- From y=ymin to y=1 (disk radius)
  y_grid(2) = y_min
  do j = 2, ngrid_uni-1
    y_grid(j + 1) = y_grid(j) + dy_grid
  enddo
  y_grid(1) = -y_grid(3)

  !--- From y=1 to y=change_y_pos
  raty1 = (change_y_pos / 1.0d0) ** (1.0d0 / real(ngridy1 - 1))
  do j = 1, ngridy1
    y_grid(ngrid_uni + j) = y_grid(ngrid_uni + j-1) * raty1
  enddo

  !--- From y=change_y_pos to y=ymax
  raty2 = (ymax / change_y_pos) ** (1.0d0 / real(ngridy2 + 1))
  do j = 0, ngridy2 - 2
    y_grid(ngrid_uni + ngridy1 + j) = change_y_pos * raty2 ** real(ngridy2 - 1 - j)
  enddo
  do j = 0, ngridy2 - 2
    y_grid(ngrid_uni + ngridy1 + j) = 2.0d0 * (ymax + change_y_pos) / 2.0d0 - y_grid(ngrid_uni + ngridy1 + j)
  enddo


!-- x-grid 
!-- 

  dx_grid2 = xuni / real(ngridx2-1, kind=8)
  dx_grid1 = dx_grid2 
  write(*,*) 'dx_grid1=', dx_grid1

  ngridx1 = int(dlog(1.d0+xmin*(1.d0-ratx1)/dx_grid1)/dlog(ratx1))

  ngridx3 = int(dlog(1.d0-(xmax-xuni)*(1.d0-ratx2)/dx_grid1)/dlog(ratx2))

  ngridx=ngridx1+ngridx2+ngridx3

  write(*,*) 'ngridx1=', ngridx1
  write(*,*) 'ngridx2=', ngridx2
  write(*,*) 'ngridx3=', ngridx3
  write(*,*) 'ngridx=', ngridx
  
  allocate(x_grid(ngridx))

  x_grid(ngridx1)=-1.d0
! x_grid(1)=xmin 
  do i=1,ngridx1-1
   x_grid(ngridx1-i)=x_grid(ngridx1-i+1)-dx_grid1*ratx1**(i-1)
  enddo

  write(*,*) 'x_grid(1)=', x_grid(1)

  do i=1,ngridx2
   x_grid(ngridx1+i)=x_grid(ngridx1+i-1)+dx_grid2
  enddo

  do i=1,ngridx3
   x_grid(ngridx1+ngridx2+i)=x_grid(ngridx1+ngridx2+i-1)+dx_grid2*ratx2**(i-1)
  enddo

  allocate(x_grid_stgu(ngridx+1),y_grid_stgu(ngridy))
  allocate(x_grid_stgv(ngridx+1),y_grid_stgv(ngridy+1))
  allocate(dxgrid(ngridx-1),dygrid(ngridy-1))
  allocate(dxgridstgu(ngridx),dygridstgu(ngridy-1))
  allocate(dxgridstgv(ngridx),dygridstgv(ngridy))
  allocate(dxgridc(ngridx-2),dygridc(ngridy-2))
  allocate(dxgridstgucc(ngridx-1),dygridstgucc(ngridy-2))
  allocate(dxgridstgvcc(ngridx-1),dygridstgvcc(ngridy-1)) 

!---Distances between grid nodes for the initial mesh (x-direction)      
  do i=1,ngridx-1
   dxgrid(i)=x_grid(i+1)- x_grid(i)
  enddo
!---Distances between grid nodes for the initial mesh (y-direction)      
  do j=1,ngridy-1
   dygrid(j)=y_grid(j+1)- y_grid(j)
  enddo

!---Distances between calculation points for the initial mesh
  do i=1,ngridx-2
   dxgridc(i)=dxgrid(i+1)/2 + dxgrid(i)/2
  enddo
        
  do j=1,ngridy-2
   dygridc(j)=dygrid(j+1)/2 + dygrid(j)/2
  enddo
          
!---  Staggered grid for u velocity 
!---  x-direction
  x_grid_stgu(1)=x_grid(1)
   
  x_grid_stgu(ngridx+1)=x_grid(ngridx)

  do i=1,ngridx-1 
   x_grid_stgu(i+1)=0.5*(x_grid(i)+x_grid(i+1))
  enddo
  
!--- y-direction     
  do j=1,ngridy
   y_grid_stgu(j)=y_grid(j)
  enddo

!--- Distances between nodes in the staggered grid of u-velocity     
  do i=1,ngridx
   dxgridstgu(i)=x_grid_stgu(i+1)- x_grid_stgu(i)
  enddo
       
  do j=1,ngridy-1
   dygridstgu(j)=y_grid_stgu(j+1) - y_grid_stgu(j)
  enddo

!--- Distances between calculation points in the staggered grid of u-velocity 
! dxgridstgucc(1)=dxgridstgu(2)/2 + dxgridstgu(1)   
! dxgridstgucc(ngridx-1)=dxgridstgu(ngridx) + dxgridstgu(ngridx-1)/2 
   
! do i=2,ngridx-2
  do i=1,ngridx-1
!  dxgridstgucc(i)=dxgridstgu(i+1)/2 + dxgridstgu(i)/2
   dxgridstgucc(i)=dxgrid(i)
  enddo
        
  do j=1,ngridy-2
!  dygridstgucc(j)=dygridstgu(j+1)/2 + dygridstgu(j)/2
   dygridstgucc(j)=dygridc(j)
  enddo

!---  Staggered grid for v velocity 

  do i=2,ngridx+1
   x_grid_stgv(i)=x_grid(i-1)
  enddo
  x_grid_stgv(1)= x_grid_stgv(2)

  y_grid_stgv(1)=y_grid(1)
  y_grid_stgv(ngridy+1)=y_grid(ngridy)
 
  do j=1,ngridy-1
   y_grid_stgv(j+1)=0.5*(y_grid(j)+y_grid(j+1))
  enddo   

!---- Distances between nodes in the staggered grid of v-velocity 
        
  do i=2,ngridx
   dxgridstgv(i)=x_grid_stgv(i+1)- x_grid_stgv(i)
  enddo
  dxgridstgv(1)=dxgridstgv(2)
       
  do j=1,ngridy
   dygridstgv(j)=y_grid_stgv(j+1) - y_grid_stgv(j)
  enddo

!---- Distances between calculation points in the staggered grid of v-velocity

! do i=2,ngridx-1
  do i=2,ngridx-2
!  dxgridstgvcc(i)=dxgridstgv(i+1)/2 + dxgridstgv(i)/2
   dxgridstgvcc(i)=dxgridc(i)
  enddo
  dxgridstgvcc(1)=0.5*dxgridstgv(2)
  dxgridstgvcc(ngridx-1)=0.5*dxgridstgv(ngridx-2)

! dygridstgvcc(1)=dygridstgv(2)/2 + dygridstgv(1)
! dygridstgvcc(ngridy-1)=dygridstgv(ngridy) + dygridstgv(ngridy-1)/2
    
! do j=2,ngridy-2
  do j=1,ngridy-1
!  dygridstgvcc(j)=dygridstgv(j+1)/2 + dygridstgv(j)/2
   dygridstgvcc(j)=dygrid(j)
  enddo
        
!---- write x,y coordinates of initial and staggered grids          

  do i=1,ngridx
   do j=1,ngridy
    write (11,*) x_grid(i), y_grid(j)
   enddo
   write (11,*)
  enddo
      
  do j=1,ngridy
   do i=1,ngridx
    write (12,*) x_grid(i), y_grid(j)
   enddo
   write (12,*)
  enddo
      
  close(12)
   
  do i=1,ngridx+1
   do j=1,ngridy
    write (13,*) x_grid_stgu(i), y_grid_stgu(j)
   enddo
   write (13,*)
  enddo
  
  close(13)
       
  do j=1,ngridy
   do i=1,ngridx+1
    write (14,*) x_grid_stgu(i), y_grid_stgu(j)
   enddo
   write (14,*)
  enddo
      
  close(14)
      
  do i=1,ngridx+1
   do j=1,ngridy+1
    write (15,*) x_grid_stgv(i), y_grid_stgv(j)
   enddo
   write (15,*)
  enddo

  do j=1,ngridy+1
   do i=1,ngridx+1
    write(16,*) x_grid_stgv(i), y_grid_stgv(j)
   enddo
   write(16,*) 
  enddo

  END Subroutine GRID_DEF
  
!--------------------------------------------------------------
  Subroutine TURBVIS
!--------------------------------------------------------------
 use rans  

 integer ::i,j
 double precision :: F,bkl,akarm,psim,aka,xi

!--- Estimation of the turbulent viscosity due to ambient turbulence   
!--- From anisotropy of atmospheric turbulence TKE=0.945*sigmax^2
!---
   Cmu=0.09
   do i=1,ngridx-1 
    do j=1,ngridy-1 
     vtold(i,j)=vt(i,j)
     if(TE(i,j).gt.1e-8) then
      vt(i,j)=Cmu*sqrt(TE(i,j))*2
     else
      vt(i,j)=1.d0/Re             
     endif
     vt(i,j)=urfvis*vt(i,j)+(1.d0-urfvis)*vtold(i,j)
    enddo
   enddo

END Subroutine TURBVIS

!--------------------------------------------------------------------------------
  Subroutine MOMENTUM_U(errmaxu)
!--------------------------------------------------------------------------------
  use rans  
  
  double precision :: err, errmaxu
  double precision :: Fe,Fw,Fn,Fs,De,Dw,Dn,Ds,Pe,Pw,Pn,Ps
  double precision :: uen,ues,uwn,uws
  double precision :: yc,alpha
      
  integer           :: i,j,isu
  
  allocate (AE(ngridx,ngridy-1), AW(ngridx,ngridy-1), AN(ngridx,ngridy-1), AS(ngridx,ngridy-1),         &
            AP(ngridx,ngridy-1), BB(ngridx,ngridy-1), DU(ngridx,ngridy-1))
          
!====INFLOW BOUNDARY CONDITIONS                        
  do j=1,ngridy-1
!--- To be completed
   DU(1,j)=1.d0
  enddo
       
!====OUTFLOW BOUNDARY CONDITIONS                        
         
  do j=2,ngridy-1
!--- To be completed
   DU(ngridx,j)=1.d0      
  enddo
      
!------ BOUNDARY CONDITIONS AT SYMMETRY AXIS   
    
  do i=2,ngridx
!--- To be completed
   DU(i,1)=1.d0     
  enddo

!------ BOUNDARY CONDITIONS AT UPPER BOUNDARY  
  do i=2,ngridx-1
!--- To be completed
   DU(i,ngridy-1)=1.d0
  end do  
    
!-------------------------------------------------------------------
!----- CALCULATION OF COEFFICIENTS FOR THE INTERNAL GRID POINTS      
!------------------------------------------------------------------- 
  do i=2,ngridx-1
   do j=2,ngridy-2

!
!---- CONVECTION TERMS
!
    yc = y_grid_stgu(j)+0.5d0*dygridstgu(j)
    Fe = 0.5*(UVEL(i,j)+UVEL(i+1,j))*dygridstgu(j)*2*pi*yc
    Fw = 0.5*(UVEL(i,j)+UVEL(i-1,j))*dygridstgu(j)*2*pi*yc
    Fn = 0.5*(VVEL(i+1,j+1)+VVEL(i,j+1))*dxgridstgu(i)*2*pi*y_grid_stgu(j+1)
    Fs = 0.5*(VVEL(i+1,j)+VVEL(i,j))*dxgridstgu(i)*2*pi*y_grid_stgu(j)
    
!
!---- DIFFUSION TERMS
!
    De=(1.d0/Re+vt(i  ,j))*dygridstgu(j)/dxgridstgucc(i)*2*pi*yc
    Dw=(1.d0/Re+vt(i-1,j))*dygridstgu(j)/dxgridstgucc(i-1)*2*pi*yc
    wce = 0.5*dxgrid(i)/dxgridstgucc(i)
    wcw = 0.5*dxgrid(i-1)/dxgridstgucc(i)
    visc =wcw*vt(i,j)  +wce*vt(i-1,j)
    visnn=wcw*vt(i,j+1)+wce*vt(i-1,j+1)
    visss=wcw*vt(i,j-1)+wce*vt(i-1,j-1)
!   visc =0.5*vt(i,j)  +0.5*vt(i-1,j)
!   visnn=0.5*vt(i,j+1)+0.5*vt(i-1,j+1)
!   visss=0.5*vt(i,j-1)+0.5*vt(i-1,j-1)
    wcn=0.5d0*dygridstgu(j)/dygridstgucc(j)
    wcs=0.5d0*dygridstgu(j)/dygridstgucc(j-1)
    wnn=0.5d0*dygridstgu(j+1)/dygridstgucc(j)
    wss=0.5d0*dygridstgu(j-1)/dygridstgucc(j-1)
    visn=wnn*visc+wcn*visnn
    viss=wss*visc+wcs*visss
    Dn=(1.d0/Re+visn)*dxgridstgu(i)/dygridstgucc(j)*2*pi*y_grid_stgu(j+1)
    Ds=(1.d0/Re+viss)*dxgridstgu(i)/dygridstgucc(j-1)*2*pi*y_grid_stgu(j)

    Pe = Fe / De
    Pw = Fw / Dw

    if(Dn.ne.0) then
     Pn = Fn / Dn
    else
     Pn = 0.d0
    endif

    if(Ds.ne.0) then
     Ps = Fs / Ds
    else
     Ps = 0.d0
    endif

!
!----- COEFFICIENTS OF DISCRETIZED EQUATION
!

    AE(i,j) = De*alpha(Pe) + dmax1(-Fe,0.d0)
    AW(i,j) = Dw*alpha(Pw) + dmax1( Fw,0.d0)
    AN(i,j) = Dn*alpha(Pn) + dmax1(-Fn,0.d0)
    AS(i,j) = Ds*alpha(Ps) + dmax1( Fs,0.d0)

    AP(i,j) = AE(i,j)+AW(i,j)+AN(i,j)+AS(i,j)
!
!---- PRESSURE GRADIENT GOES TO THE SOURCE TERM
!
    BB(i,j) = (PTOT(i-1,j)-PTOT(i,j))*dygridstgu(j)*2*pi*yc 

!
!---- INSERT WIND TURBINE THRUST
!

!--- To be completed

!
!---- FOR THE VELOCITY CORRECTION
!
    DU(i,j)=1.d0/AP(i,j)

!-----insert underrelaxation----------------------------------------------------
    AP(i,j) = AP(i,j)/urfu
    BB(i,j) = BB(i,j) + AP(i,j)*UVEL(i,j)*(1-urfu)
    DU(i,j) = DU(i,j)*urfu

   enddo
  enddo

!---- Store velocities of the previous iteration UOLD
!---- Store DU(i,j)=AREA/AE(i,j) needed for velocity correction

  do i=1,ngridx
   do j=1,ngridy-1
    UOLD(i,j) = UVEL(i,j)
    yc = y_grid(j)+0.5*dygrid(j)
    DU(i,j) = DU(i,j)*(2*pi*yc*dygrid(j))
   enddo
  enddo
  
!-----SOLVE THE SYSTEM OF EQUATIONS USING ADI
   
  do isu=1,nswp
   call LISOLV(ngridx,ngridy-1,AE,AW,AN,AS,AP,BB,UVEL)
  enddo
  
!--- Find maxiumum absolute error of u-velocity
  errmaxu=tiny
  do i=2,ngridx-1
   do j=2,ngridy-2
    err=dabs(UVEL(i,j)-UOLD(i,j))
    errmaxu=dmax1(err,errmaxu)
   enddo
  enddo

  deallocate (AE,AW,AN,AS,BB,AP)

  END Subroutine MOMENTUM_U

!--------------------------------------------------------------------------------
  Subroutine MOMENTUM_V(errmaxv)
!--------------------------------------------------------------------------------
  use rans  
  
  double precision :: err, errmaxv
  double precision :: Fe,Fw,Fn,Fs,De,Dw,Dn,Ds,Pe,Pw,Pn,Ps,Ve,Vw,Vn,Vs
  double precision :: un,up,us
  double precision :: yc,alpha

  integer          :: i,j,isv
    
  allocate  ( AE(ngridx,ngridy), AW(ngridx,ngridy), AN(ngridx,ngridy), AS(ngridx,ngridy),    & 
              AP(ngridx,ngridy), BB(ngridx,ngridy), DV(ngridx,ngridy))
              
!====INFLOW BOUNDARY CONDITIONS     
      
  do j=1,ngridy
!--- To be completed
  enddo

!====OUTFLOW BOUNDARY CONDITIONS          
  do j=2,ngridy-1
!--- To be completed
  end do
      
!-----SYMMETRY AXIS                
  do i=2,ngridx
!--- To be completed
   DV(i,2)=1.d0    

   DV(i,1)=1.d0
  enddo
         
!-----UPPER BOUNDARY CONDITIONS 
        
  do i=2,ngridx
!--- To be completed
   DV(i,ngridy)=1.d0
  end do

!--------CALCULATION OF COEFFICIENTS FOR THE INTERNAL GRID POINTS           
      
  do i=2,ngridx-1
   do j=3,ngridy-1

!
!----- CONVECTION TERMS
!
    yc=y_grid_stgv(j)+0.5d0*dygridstgv(j) 
     
    Fe = 0.5*(UVEL(i,j)+UVEL(i,j-1))*dygridstgv(j)*2*pi*yc
    Fw = 0.5*(UVEL(i-1,j)+UVEL(i-1,j-1))*dygridstgv(j)*2*pi*yc
    Fn = 0.5*(VVEL(i,j)+VVEL(i,j+1))*dxgridstgv(i)*2*pi*y_grid_stgv(j+1)
    Fs = 0.5*(VVEL(i,j)+VVEL(i,j-1))*dxgridstgv(i)*2*pi*y_grid_stgv(j) 
    
    Ve = 0.5*(VVEL(i,j)+VVEL(i+1,j))
    Vw = 0.5*(VVEL(i,j)+VVEL(i-1,j))
    Vn = 0.5*(VVEL(i,j)+VVEL(i,j+1))
    Vs = 0.5*(VVEL(i,j)+VVEL(i,j-1))

!
!----- DIFFUSION TERMS
!
    wy1 = 0.5d0*dygrid(j)/dygridstgvcc(j)
    wy2 = 0.5d0*dygrid(j-1)/dygridstgvcc(j)
    visc  = wy2*vt(i-1,j)+wy1*vt(i-1,j-1)
    visee = wy2*vt(i  ,j)+wy1*vt(i  ,j-1)
    if(i.ne.2) then
     visww = wy2*vt(i-2,j)+wy1*vt(i-2,j-1)
    else
     visww = visc
    endif
    wce  = 0.5d0*dxgridstgv(i)  /dxgridstgvcc(i)
    wee  = 0.5d0*dxgridstgv(i+1)/dxgridstgvcc(i)
    wcw  = 0.5d0*dxgridstgv(i)  /dxgridstgvcc(i-1)
    www  = 0.5d0*dxgridstgv(i-1)/dxgridstgvcc(i-1)
    vise = wce*visee + wee*visc
    visw = wcw*visww + www*visc 
    De = (1./Re + vise) * dygridstgv(j) * (1./dxgridstgvcc(i)   - Ve/yc**2) *2*pi*yc
    Dw = (1./Re + visw) * dygridstgv(j) * (1./dxgridstgvcc(i-1) - Vw/yc**2) *2*pi*yc
    Dn = (1./Re + vt(i-1,j  )) * dxgridstgv(i) * (1./dygridstgvcc(j)   - Vn/y_grid_stgv(j+1)**2)*2*pi*y_grid_stgv(j+1)
    Ds = (1./Re + vt(i-1,j-1)) * dxgridstgv(i) * (1./dygridstgvcc(j-1) - Vs/y_grid_stgv(j)**2)  *2*pi*y_grid_stgv(j)
    
    Pe=Fe/De
    Pw=Fw/Dw
    Pn=Fn/Dn
    Ps=Fs/Ds

!
!----- COEFFICIENTS OF DISCRETIZED EQUATION
!
    AE(i,j) = De*alpha(Pe) + dmax1(-Fe,0.d0)
    AW(i,j) = Dw*alpha(Pw) + dmax1( Fw,0.d0)
    AN(i,j) = Dn*alpha(Pn) + dmax1(-Fn,0.d0)
    AS(i,j) = Ds*alpha(Ps) + dmax1( Fs,0.d0)

    AP(i,j) = AE(i,j)+AW(i,j)+AN(i,j)+AS(i,j)
!
!---- PRESSURE GRADIENT GOES TO THE SOURCE TERM
!
    BB(i,j) = (PTOT(i-1,j-1)-PTOT(i-1,j))*dxgridstgv(i)*2*pi*yc
!
!---- FOR THE VELOCITY CORRECTION
!
    DV(i,j) = 1.d0/AP(i,j)

!-------insert underrelaxation-----------
    AP(i,j) = AP(i,j)/urfv    
    BB(i,j) = BB(i,j) + AP(i,j)*VVEL(i,j)*(1-urfv) 
    DV(i,j) = DV(i,j)*urfv
!------------------------------------------------------    
   enddo
  enddo

!---- Store velocities of the previous iteration VOLD
!---- Store DV(i,j)=AREA/AE(i,j) needed for velocity correction

  do i=2,ngridx
   do j=1,ngridy
    VOLD(i,j)=VVEL(i,j)
    DV(i,j) = DV(i,j)*(2*pi*y_grid(j)*dxgridstgv(i))
   enddo
  enddo
  
!--------SOLUTION OF SYSTEM OF EQUATIONS USING ADI

  do isv=1,nswp
   call LISOLV(ngridx,ngridy,AE,AW,AN,AS,AP,BB,VVEL)
  enddo 

  errmaxv=tiny
  do i=2,ngridx-1
   do j=3,ngridy-1
    err=dabs(VVEL(i,j)-VOLD(i,j))
    errmaxv=dmax1(err,errmaxv)
   enddo
  enddo

  deallocate (AE,AW,AN,AS,BB,AP)
   
  END Subroutine MOMENTUM_V      


!--------------------------------------------------------------------------------
  Subroutine PRESCOR(errmaxp)
!--------------------------------------------------------------------------------
  use rans  
 
  integer           :: i,j,isp
  
  double precision  :: err,errmaxp
  double precision  :: yc,yn,ys
  
  allocate   ( AE(ngridx-1,ngridy-1), AW(ngridx-1,ngridy-1), AN(ngridx-1,ngridy-1), AS(ngridx-1,ngridy-1),  &
               BB(ngridx-1,ngridy-1), AP(ngridx-1,ngridy-1) )
               
!-----INFLOW BOUNDARY 

   do j=1,ngridy-2 
    yc = y_grid(j)+0.5d0*dygrid(j)
    yn = y_grid(j+1)
    ys = y_grid(j)

    AN(1,j) = 0.d0                            
    AS(1,j) = 0.d0                            
    AP(1,j) = 1.d0                    
    AE(1,j) = 1.d0                    
    BB(1,j) = 0.d0                                     

!   AE(1,j) = DU(2,j)  *(2*pi*yc*dygrid(j))
!   AN(1,j) = DV(2,j+1)*(2*pi*yn*dxgrid(1))
!   AS(1,j) = DV(2,j)  *(2*pi*ys*dxgrid(1))
!   AP(1,j) = AE(1,j)+AN(1,j)+AS(1,j)
!   BB(1,j) = (UVEL(1,j)-UVEL(2,j))*dygrid(j)*2*pi*yc + (VVEL(2,j)*ys-VVEL(2,j+1)*yn)*dxgrid(1)*2*pi 

   enddo
     
!-----OUTFLOW BOUNDARY   
     
   do j=2,ngridy-2 
    yc = y_grid(j)+0.5d0*dygrid(j)
    yn = y_grid(j+1)
    ys = y_grid(j)
    AW(ngridx-1,j) = DU(ngridx-1,j) *(2*pi*yc*dygrid(j))
    AN(ngridx-1,j) = DV(ngridx,j+1)*(2*pi*yn*dxgrid(ngridx-1))
    AS(ngridx-1,j) = DV(ngridx,j)  *(2*pi*ys*dxgrid(ngridx-1))
    AP(ngridx-1,j) = AW(ngridx-1,j)+AN(ngridx-1,j)+AS(ngridx-1,j)
    BB(ngridx-1,j) =  (UVEL(ngridx-1,j)   -UVEL(ngridx,j)   )*dygrid(j)*yc    *2*pi +                    &
                      (VVEL(ngridx,j)*ys-VVEL(ngridx,j+1)*yn)*dxgrid(ngridx-1)*2*pi 

   enddo
   
   
!---- SYMMETRY AXIS 
    
   do i=2,ngridx-1
    yc = y_grid(2) + 0.5d0*dygrid(2)
    yn = y_grid(3)
    ys = y_grid(2)

    AE(i,2) = DU(i+1,2) *(2*pi*yc*dygrid(2))
    AW(i,2) = DU(i,2)   *(2*pi*yc*dygrid(2))
    AN(i,2) = DV(i+1,3)   *(2*pi*yn*dxgrid(i))
    AS(i,2) = DV(i+1,2)   *(2*pi*ys*dxgrid(i))
    AP(i,2) = AE(i,2)+AW(i,2)+AN(i,2)+AS(i,2)
    BB(i,2) = (UVEL(i,2)-UVEL(i+1,2))*dygrid(2)*2*pi*yc + VVEL(i+1,1)*dxgrid(i)*2*pi*yn

    AE(i,1) = 0.d0
    AW(i,1) = 0.d0
    AN(i,1) = 1.d0
    AP(i,1) = 1.d0   
    BB(i,1) = 0.d0
   enddo

!----UPPER BOUNDARY

   do i=1,ngridx-1
    yc = y_grid(ngridy-1) + 0.5*dygrid(ngridy-1)
    ys = y_grid(ngridy-1)
    yn = y_grid(ngridy)

    AE(i,ngridy-1) = DU(i+1,ngridy-1) *(2*pi*yc*dygrid(ngridy-1))
    AW(i,ngridy-1) = DU(i,ngridy-1)   *(2*pi*yc*dygrid(ngridy-1))
    AS(i,ngridy-1) = DV(i+1,ngridy-1) *(2*pi*ys*dxgrid(i)       )
    AP(i,ngridy-1) = AE(i,ngridy-1)+AW(i,ngridy-1)+AS(i,ngridy-1)
    BB(i,ngridy-1) = (UVEL(i,ngridy-1)     - UVEL(i+1,ngridy-1)     )*dygrid(ngridy-1)*2*pi*yc +               &
                     (VVEL(i+1,ngridy-1)*ys- VVEL(i+1,ngridy)    *yn)*dxgrid(i)       *2*pi

  enddo

!------PRESSURE  COEFFICIENTS FOR INTERNAL GRID    

  do i=2,ngridx-2
   do j=3,ngridy-2

    yc = y_grid(j)+0.5d0*dygrid(j)
    yn = y_grid(j+1)
    ys = y_grid(j)

    AE(i,j) = DU(i+1,j)*(2*pi*yc*dygrid(j))
    AW(i,j) = DU(i,j)  *(2*pi*yc*dygrid(j))
    AN(i,j) = DV(i+1,j+1)*(2*pi*yn*dxgrid(i))
    AS(i,j) = DV(i+1,j)  *(2*pi*ys*dxgrid(i))
    AP(i,j) = AE(i,j)+AW(i,j)+AN(i,j)+AS(i,j)
    BB(i,j) = (UVEL(i,j)-UVEL(i+1,j))*dygrid(j)*2*pi*yc + (VVEL(i+1,j)*ys-VVEL(i+1,j+1)*yn)*dxgrid(i)*2*pi

   enddo
  enddo
      
!--------PRESSURE SOLUTION WITH ADI
  do isp=1,nswp
   call LISOLV(ngridx-1,ngridy-1,AE,AW,AN,AS,AP,BB,PRE)
   errmaxp = tiny
   do i=2,ngridx-2
    do j=2,ngridy-2
!  do i=1,ngridx-1
!   do j=1,ngridy-1
     err=AP(i,j)*PRE(i,j)-AE(i,j)*PRE(i+1,j)-AW(i,j)*PRE(i-1,j)-AN(i,j)*PRE(i,j+1)-AS(i,j)*PRE(i,j-1)-BB(i,j)
!    errmaxp=dmax1(dabs(err),errmaxp)
     errmaxp=dmax1(dabs(BB(i,j)),errmaxp)
    enddo
   enddo
  enddo 
!------------------------------------   
   
  deallocate (AE,AW,AN,AS,BB,AP)

  END Subroutine PRESCOR
   
   
   
!--------------------------------------------------------------------------------
      Subroutine VELOCOR
!--------------------------------------------------------------------------------
      use rans  
            
!-------------------------------------------------------
!---  update solution for UVEL----------------------------    
!-------------------------------------------------------      
      do i=2,ngridx-1
        do j=2,ngridy-2
         UVEL(i,j)=UVEL(i,j)+ DU(i,j)*(PRE(i-1,j)-PRE(i,j))
        enddo
      enddo
     
!---------------------------------------------------------------        
!---  update solution for VVEL ---------------------------------
!---------------------------------------------------------------
    
      do i=2,ngridx-1
        do j=3,ngridy-1
         VVEL(i,j)=VVEL(i,j)+DV(i,j)*(PRE(i-1,j-1)-PRE(i-1,j))
        enddo
      enddo  
      
     deallocate (DU,DV)
                       
     END Subroutine VELOCOR
      
      
  Subroutine CALCTE(errmaxte)
!--------------------------------------------------------------------------------
  use rans  
  
  double precision :: err, errmaxte
  double precision :: Fe,Fw,Fn,Fs,De,Dw,Dn,Ds,Pe,Pw,Pn,Ps
  double precision :: uen,ues,uwn,uws
  double precision :: yc,alpha
      
  integer           :: i,j,isu
  
  allocate (AE(ngridx-1,ngridy-1), AW(ngridx-1,ngridy-1), AN(ngridx-1,ngridy-1), AS(ngridx-1,ngridy-1),         &
            AP(ngridx-1,ngridy-1), BB(ngridx-1,ngridy-1))
          
  sigmak = 1.d0
  cmu = 0.033d0
  cd = 1.d0

!====INFLOW BOUNDARY CONDITIONS                        
  do j=1,ngridy-2
!--- To be completed
  enddo
       
!====OUTFLOW BOUNDARY CONDITIONS                        
         
  do j=2,ngridy-2
!--- To be completed
  enddo
      
!------ BOUNDARY CONDITIONS AT SYMMETRY AXIS   
    
  do i=2,ngridx-1
!--- To be completed
  enddo

!------ BOUNDARY CONDITIONS AT UPPER BOUNDARY  
  do i=1,ngridx-1
!--- To be completed
  end do  
    
!-------------------------------------------------------------------
!----- CALCULATION OF COEFFICIENTS FOR THE INTERNAL GRID POINTS      
!------------------------------------------------------------------- 
  do i=2,ngridx-2
   do j=2,ngridy-2

!
!---- CONVECTION TERMS
!
    yc = y_grid(j)+0.5d0*dygrid(j)
    yn = y_grid(j+1)
    ys = y_grid(j)

    Fe = UVEL(i+1,j)*dygrid(j)*2*pi*yc
    Fw = UVEL(i  ,j)*dygrid(j)*2*pi*yc
    Fn = VVEL(i+1,j+1)*dxgrid(i)*2*pi*yn
    Fs = VVEL(i+1,j)*dxgrid(i)*2*pi*ys
    
!
!---- DIFFUSION TERMS
!
    vise = 0.5*(vt(i,j)+vt(i+1,j))
    visw = 0.5*(vt(i,j)+vt(i-1,j))
    visn = 0.5*(vt(i,j)+vt(i,j+1))
    viss = 0.5*(vt(i,j)+vt(i,j-1))
    De=(1.d0/Re+vise/sigmak)*dygrid(j)/dxgridc(i)*2*pi*yc
    Dw=(1.d0/Re+visw/sigmak)*dygrid(j)/dxgridc(i-1)*2*pi*yc
    Dn=(1.d0/Re+visn/sigmak)*dxgrid(i)/dygridc(j)*2*pi*yn
    Ds=(1.d0/Re+viss/sigmak)*dxgrid(i)/dygridc(j-1)*2*pi*ys

    Pe = Fe / De
    Pw = Fw / Dw

    if(Dn.ne.0) then
     Pn = Fn / Dn
    else
     Pn = 0.d0
    endif

    if(Ds.ne.0) then
     Ps = Fs / Ds
    else
     Ps = 0.d0
    endif

!
!----- COEFFICIENTS OF DISCRETIZED EQUATION
!

    AE(i,j) = De*alpha(Pe) + dmax1(-Fe,0.d0)
    AW(i,j) = Dw*alpha(Pw) + dmax1( Fw,0.d0)
    AN(i,j) = Dn*alpha(Pn) + dmax1(-Fn,0.d0)
    AS(i,j) = Ds*alpha(Ps) + dmax1( Fs,0.d0)

    AP(i,j) = AE(i,j)+AW(i,j)+AN(i,j)+AS(i,j)
!
!---- SOURCE TERMS                             
!
!----- PRODUCTION OF TKE BY REYNOLDS STRESSES
!
!----- CALCULATE VELOCITY GRADIENTS (STRAIN TENSOR)
!
    dudx = (UVEL(i+1,j  )-UVEL(i,j  ))/dxgrid(i)
    dvdy = (VVEL(i+1,j+1)-VVEL(i+1,j))/dygrid(j)
!-----
    uc   = 0.5d0*(UVEL(i,j  )+UVEL(i+1,j  ))
    unn  = 0.5d0*(UVEL(i,j+1)+UVEL(i+1,j+1))
    uss  = 0.5d0*(UVEL(i,j-1)+UVEL(i+1,j-1))
    wcn  = 0.5d0*dygrid(j)  /dygridc(j)
    wnn  = 0.5d0*dygrid(j+1)/dygridc(j)
    wcs  = 0.5d0*dygrid(j)  /dygridc(j-1)
    wss  = 0.5d0*dygrid(j-1)/dygridc(j-1)
    un   = wcn*unn + wnn*uc
    us   = wss*uc  + wcs*uss
    dudy = (un-us)/dygrid(j) 
    vc   = 0.5d0*(VVEL(i+1,j)+VVEL(i+1,j+1))
    if (i.ne.ngridx-2) then
     vee = 0.5d0*(VVEL(i+2,j)+VVEL(i+2,j+1))
    else
     vee = vc 
    endif
    vww  = 0.5d0*(VVEL(i,j)+VVEL(i,j+1))
    wce  = 0.5d0*dxgrid(i)  /dxgridc(i)
    wee  = 0.5d0*dxgrid(i+1)/dxgridc(i)
    wcw  = 0.5d0*dxgrid(i)  /dxgridc(i-1)
    www  = 0.5d0*dxgrid(i-1)/dxgridc(i-1)
    ve   = wce*vee + wee*vc
    vw   = wcw*vww + www*vc 
    dvdx = (ve-vw)/dxgrid(i) 
!------ ATTENTION!!!   INTEGRATION ON VOLUME
    vol  = 2.*pi*yc*dygrid(j)*dxgrid(i)
    prod = vt(i,j)*(2.*(dudx**2+dvdy**2)+(dudy+dvdx)**2)*vol

!----- DISSIPATION OF TKE
!------ 1st OPTION: NON-LINEARIZATION OF SOURCE TERM
!   ed = cmu*TE(i,j)**2/vt(i,j)*vol
!   ed = cd*(TE(i,j)**1.5)/2*vol
!   BB(i,j) = prod - ed 
    
!------ 2nd OPTION: LINEARIZATION OF SOURCE TERM
   
    BB(i,j) = prod  
!   AP(i,j) = AP(i,j) + cmu*TE(i,j)/vt(i,j)*vol 
    AP(i,j) = AP(i,j) + cd*TE(i,j)**0.5/2*vol 

!-----insert underrelaxation----------------------------------------------------
    AP(i,j) = AP(i,j)/urfte
    BB(i,j) = BB(i,j) + AP(i,j)*TE(i,j)*(1-urfte)

   enddo
  enddo

!---- Store TKEs of the previous iteration TKEOLD
!---- Store DU(i,j)=AREA/AE(i,j) needed for velocity correction

  do i=1,ngridx-1
   do j=1,ngridy-1
    TEOLD(i,j) = TE(i,j)
   enddo
  enddo
  
!-----SOLVE THE SYSTEM OF EQUATIONS USING ADI
   
  do isu=1,nswp
   call LISOLV(ngridx-1,ngridy-1,AE,AW,AN,AS,AP,BB,TE)
  enddo
  
!--- Find maxiumum absolute error of u-velocity
  errmaxte=tiny
  do i=2,ngridx-2
   do j=2,ngridy-2
    err=dabs(TE(i,j)-TEOLD(i,j))
    errmaxte=dmax1(err,errmaxte)
!   write(*,*) 'TE=', TE(i,j)
   enddo
  enddo

  deallocate (AE,AW,AN,AS,BB,AP)

  END Subroutine CALCTE     
            
  SUBROUTINE LISOLV(NI,NJ,AE,AW,AN,AS,AP,ASU,PHI)
    double precision :: AE(NI,NJ),AW(NI,NJ),AN(NI,NJ),AS(NI,NJ),AP(NI,NJ),ASU(NI,NJ),PHI(NI,NJ)
    double precision, allocatable :: A(:),B(:),C(:),D(:) 
    double precision TERM
    integer I,II,J,JJ,NI,NJ,NIM1,NJM1
    NIM1=NI-1 
    NJM1=NJ-1 
    allocate (A(NJ),B(NJ),C(NJ),D(NJ))
!---- COMMENCE W-E SWEEP
!---- WE INCLUDE I=1 AND I=NI IN ORDER TO CALCULATE THE CORNER POINTS (1,1),(1,NJ),(NI,1),(NI,NJ) 
!---- NO NEED TO REPEAT FOR THE N-S SWEEP 
    do I=1,NI
     A(1)=AN(I,1)
     if(I.ne.1.and.I.ne.NI) then
      C(1)=(AE(I,1)*PHI(I+1,1)+AW(I,1)*PHI(I-1,1)+ASU(I,1))/AP(I,1)
     elseif (I.eq.1) then
      C(1)=(AE(I,1)*PHI(I+1,1)+ASU(I,1))/AP(I,1)
     else
      C(1)=(AW(I,1)*PHI(I-1,1)+ASU(I,1))/AP(I,1)
     endif
!-----COMMENCE S-N TRAVERSE 
     do J=2,NJ
!-----ASSEMBLE TDMA COEFFICIENTS
      A(J)=AN(I,J)
      B(J)=AS(I,J)
      if(I.ne.1.and.I.ne.NI) then
       C(J)=AE(I,J)*PHI(I+1,J)+AW(I,J)*PHI(I-1,J)+ASU(I,J)
      elseif (I.eq.1) then
       C(J)=(AE(I,J)*PHI(I+1,J)+ASU(I,J))/AP(I,J)
      else
       C(J)=(AW(I,J)*PHI(I-1,J)+ASU(I,J))/AP(I,J)
      endif
      D(J)=AP(I,J)
!-----CALCULATE COEFFICIENTS OF RECURRENCE FORMULA
      TERM=1./(D(J)-B(J)*A(J-1))
      A(J)=A(J)*TERM
      C(J)=(C(J)+B(J)*C(J-1))*TERM
     enddo
!-----OBTAIN NEW PHI'S
     PHI(I,NJ)=C(NJ)
     do JJ=2,NJ 
      J=NJ+1-JJ 
      PHI(I,J)=A(J)*PHI(I,J+1)+C(J) 
     enddo
    enddo     
    deallocate (A,B,C,D)

    allocate (A(NI),B(NI),C(NI),D(NI))
!-----COMMENCE N-S SWEEP
    do J=2,NJM1
     A(1)=AE(1,J)
     C(1)=(AN(1,J)*PHI(1,J+1)+AS(1,J)*PHI(1,J-1)+ASU(1,J))/AP(1,J)
!-----COMMENCE W-E TRAVERSE 
!    do I=2,NIM1
     do I=2,NI
!-----ASSEMBLE TDMA COEFFICIENTS
      A(I)=AE(I,J)
      B(I)=AW(I,J)
      C(I)=AN(I,J)*PHI(I,J+1)+AS(I,J)*PHI(I,J-1)+ASU(I,J)
      D(I)=AP(I,J)
!-----CALCULATE COEFFICIENTS OF RECURRENCE FORMULA
      TERM=1./(D(I)-B(I)*A(I-1))
      A(I)=A(I)*TERM
      C(I)=(C(I)+B(I)*C(I-1))*TERM
     enddo
     PHI(NI,J)=C(NI)
!-----OBTAIN NEW PHI'S
     do II=2,NI 
      I=NI+1-II 
      PHI(I,J)=A(I)*PHI(I+1,J)+C(I) 
     enddo
    enddo     
    deallocate (A,B,C,D)
    return
    end 
 
    double precision FUNCTION alpha(x)
     double precision, intent(in) :: x
!     alpha = 1.d0
!     alpha = 1.d0-0.5*dabs(x)
      alpha = dmax1(0.d0,1.d0-0.5d0*dabs(x))

    END FUNCTION alpha
