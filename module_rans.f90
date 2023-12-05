   module rans   
!-------------------------------------------------------------------------------------
   implicit none 

   integer  :: ngridx1,ngridx2,ngridx,ngridy,ngrid_uni,ngridy1,ngridy2
   integer  :: it,itmax,nswp,nvt,istab,npos

   double precision ::  pi,eps,tiny
   double precision ::  xmin,xmax,xuni,ymin,change_y_pos,ymax,ratx1,ratx2,raty1,raty2
   double precision ::  urfu,urfv,urfp,urfte,urfvis
   double precision ::  Uinf, dvisc, radius, Re, ct 
   double precision ::  tiamb,vtamb
   double precision, allocatable :: y_grid(:),x_grid(:),x_grid_stgu(:),                           &
                                    y_grid_stgu(:),x_grid_stgv(:),y_grid_stgv(:)          

   double precision, allocatable :: dxgrid(:),dygrid(:),dxgridc(:),dygridc(:),                    &
                                    dxgridstgu(:),dygridstgu(:),dxgridstgucc(:),dygridstgucc(:),  &
                                    dxgridstgv(:),dygridstgv(:),dxgridstgvcc(:),dygridstgvcc(:)    

   double precision, allocatable :: AP(:,:),AE(:,:),AW(:,:),AN(:,:),AS(:,:),BB(:,:),DU(:,:),DV(:,:) 

   double precision, allocatable :: UVEL(:,:),VVEL(:,:),PRE(:,:),PTOT(:,:),UOLD(:,:),VOLD(:,:)
   double precision, allocatable :: TE(:,:),TEOLD(:,:),vt(:,:),vtold(:,:)

   double precision, allocatable :: w(:),akl(:)

   double precision, allocatable :: acor1,acor2

   end module rans   

