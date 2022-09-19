************************************************************************
!
!   User element for electro-chemical transport of ions 
!  and large elastic deformation of host in 2D or 3D.  This is for plane 
!  strain, axisymetric, and 3D. 
!  -PE and AX use 4-noded elements with 4 integration points.
!  -3D uses 8-noded elements with 8 integration points.
!
! Solution variables (or nodal variables) are the displacements(DOFs 1-3), 
!  the electric potential(DOFs 11) and 
!  the electro-chemical potential of moving ion(DOFs 12).
!
! In order to avoid locking for the fully-integrated element, we
!  use the F-bar method of de Souza Neto (1996).
!
!     
!              A eta (=xi_2)
!  4-node      |
!   quad       |Face 3
!        4-----------3
!        |     |     |
!        |     |     |
!  Face 4|     ------|---> xi (=xi_1)
!        |           | Face2
!        |           |
!        1-----------2
!          Face 1
!
!
!  8-node     8-----------7
!  brick     /|          /|       zeta
!           / |         / |       
!          5-----------6  |       |     eta
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi
!          1-----------2        origin at cube center
!
!     Face numbering follows:
!       Face 1 = nodes 1,2,3,4
!       Face 2 = nodes 5,8,7,6
!       Face 3 = nodes 1,5,6,2
!       Face 4 = nodes 2,6,7,3
!       Face 5 = nodes 3,7,8,4
!       Face 6 = nodes 4,8,5,1
!
! Sooraj Naryaan,   August 2020
***********************************************************************
!
! User element statement in the input file :
!
!  2D PE elements
!  *User Element,Nodes=4,Type=U1,Iproperties=2,Properties=?,Coordinates=2,Variables=4,Unsymm
!  1,2,11,12
!
!  2D AX elements
!  *User Element,Nodes=4,Type=U2,Iproperties=2,Properties=?,Coordinates=2,Variables=4,Unsymm
!  1,2,11,12
!
!  3D elements
!  *User Element,Nodes=8,Type=U3,Iproperties=2,Properties=?,Coordinates=3,Variables=4,Unsymm
!  1,2,3,11,12
!
!
!     State Variables
!     --------------------------------------------------------------
!     
!     Local SDV's (used for the solution procedure)
!       
!          svars(1+j) = cbar ---- normalized concentration at integration point
!          
!
!     In the input file, set 'User output variables'= number of global SDV's
!
!     In the input file, set 'ngSdv'= number of global SDV's
!
!     In the input file, set 'nlSdv'= number of local SDV's
!
!     In the input file, set 'varibles'=(nlSdv*nIntPt)
!
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     Nr     = props(1) ! Shear modulus parameter
!     Kbulk  = props(2) ! Bulk modulus
!     Omega  = props(3) ! Volume of a mole of fluid particles
!     chiL   = props(4) ! Chi parameter
!     chiH   = props(5) ! Chi parameter
!     thetaT = props(6) ! Transition temperature
!     Delta  = props(7) ! Width of transition in chi
!     mu0    = props(8) ! Chemical potential of pure fluid
!     D      = props(9) ! Coefficient of permeability
!     Rgas   = props(10) ! Universal gas constant
!     theta0 = props(11) ! Initial absolute temperature
!     phi0   = props(12) ! Initial polymer volume fraction
!     Cheat  = props(13) ! Specific heat
!     Ktherm = props(14) ! Thermal conductivity
!     nlSdv  = jprops(1) ! Number of local sdv's per integ pt
!     ngSdv  = jprops(2) ! Number of global sdv's per integ pt
!
!***********************************************************************

      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ElemOffset
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and 
      !    that value must be set here the same.

      integer numElem,ElemOffset,err

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
      parameter(numElem=99999)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=100000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable :: globalSdv(:,:,:)

      integer test

      end module global

***********************************************************************

      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1     LAYER,KSPT)

      use global

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION STATEV(NSTATV),COORDS(NCRDS)


      statev = 0.999
      test = max(test,noel)


      RETURN
      END

***********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.  Note that an offset of
      !  ElemOffset is used between the real mesh and the dummy mesh.
      !  If your model has more than ElemOffset UEL elements, then
      !  this will need to be modified.
     
      use global
     
      include 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.
      do i=1,nuvarm
         uvar(i) = globalSdv(noel-ElemOffset,npt,i)
      enddo
c      for example
c      uvar(2) = globalSdv(noel-ElemOffset,npt,2)
c      uvar(3) = globalSdv(noel-ElemOffset,npt,3)
c      uvar(4) = globalSdv(noel-ElemOffset,npt,4)

      return
      end subroutine uvarm

****************************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)
  

      use global
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      integer lenJobName,lenOutDir,nDim,nInt,nIntS
      character*256 jobName,outDir,fileName



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !parameter(nInt=8)  ! number of volume integration pionts
      parameter(nIntS=2) ! number of surface integration points
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




      !----------------------------------------------------------------
      ! 
      ! Perform initial checks
      !
      !
      ! Open the debug/error message file
      !
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      fileName = outDir(1:lenOutDir)//'\aaMSGS_'//
     +     jobName(1:lenJobName)//'.dat'
      open(unit=80,file=fileName,status='unknown')


      ! Check the procedure type, this should be a coupled
      !  temperature displacement or pore pressure displacement
      !  which are any of the following (64,65,72,73)
      !
      if((lflags(1).eq.64).or.(lflags(1).eq.65).or.
     +     (lflags(1).eq.72).or.(lflags(1).eq.73)) then
         !
         ! all is good
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and chekc the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         write(80,*) 'Abaqus does not have the right procedure'
         write(80,*) 'go back and chekc the procedure type'
         write(80,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Make sure Abaqus knows you are doing a large
      !  deformation problem; this only matters
      !  when it comes to output in viewer
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a small displacement analysis'
         write(80,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear purturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear purturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear purturbation step'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a linear purturbation step'
         call xit         
      endif


      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.0.0) return
      !
      ! Done with initial checks
      !
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! 
      ! Call the particular element to perform the analysis
      !
      if(jtype.eq.1) then
!         !
!         ! This is a plane strain analysis
!         !
!         nDim = 2
!         nInt = 4
!         call UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
!     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
!     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
!     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
!     +        NJPROP,PERIOD,
!     +        nDim,nInt,nIntS)
!         !
!         !
!      elseif(jtype.eq.2) then
!         !
!         ! This is an axisymmetric analysis
!         !
!         nDim = 2
!         nInt = 4
!         call UAX4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
!     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
!     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
!     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
!     +        NJPROP,PERIOD,
!     +        nDim,nInt,nIntS)
!         !
!         !
!      elseif(jtype.eq.3) then
         !
         ! This is a 3D analysis
         !
         nDim = 3
         nInt = 8
         call U3D8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,nInt,nIntS)
         !
         !
      else
         !
         ! We have a problem...
         !
         write(*,*) 'Element type not supported, jtype=',jtype
         write(80,*) 'Element type not supported, jtype=',jtype
         call xit
         !
      endif
      !
      ! Done with this element, RHS and AMATRX already returned
      !  as output from the specific element routine called
      !
      !----------------------------------------------------------------


      return
      end subroutine uel

************************************************************************

      subroutine U3D8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS)

      use global
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! All nodal DOF related variables
      real*8 u(nNode,mcrd),du(nNode,ndofel),uNew(nNode,ndofel)
      real*8 phiOld(nNode),dphi(nNode),phiNew(nNode)
      real*8 omgOld(nNode),domg(nNode),omgNew(nNode)
      real*8 uOld(nNode,ndofel),u_t(nNode,ndofel),v(nNode,3)
      real*8 coordsC(mcrd,nNode)

      ! All int. pt. DOF related variables
      real*8 omg_tau,omg_t,domgdx(3),domgdt
      real*8 F_tau(3,3),F_t(3,3),detF_tau,detF,detF_t
      real*8 phi_tau, phi_t,dphidX(3),dphidt
      
      !Integral variables mainly for loop counters
      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,lenJobName,lenOutDir,faceFlag,nIntS

      ! Misc. variables 
      real*8 Iden(3,3),Le,theta0,phi0,body(3)
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t
      real*8 cbarLmt,umeror
      real*8 qNet
      real*8 mob,dmobdu(3),dmobdphi,dmobdomg
      real*8 dcbardu(3),dcbardphi,dcbardomg
      real*8 SpTanMod(3,3,3,3),Spuomgmod(3,3),Spuphimod(3,3)
      real*8 Spomgumod(3,3,3)
      real*8 amatcu(3,9),Spphiumod(3,3,3)
      !
      real*8 delphi,delomg,delF
      real*8 F_pert(3,3),T_pert(3,3), cbar_pert
      real*8 phi_pert,omg_pert
      real*8 Csp(3,3,3,3),Ccbartan(3,3)
      real*8 Finv_pert(3,3), J_pert
      real*8 Finv_tau(3,3), J_tau
      real*8 Tk_tau(3,3),Tk_pert(3,3)
	  real*8 Tme_tau(3,3),Tes_tau(3,3)
	  !
	  real*8 T_pertp(3,3), cbar_pertp
      real*8 Finv_pertp(3,3), J_pertp
      real*8 Tk_pertp(3,3)
      !
      real*8 T_pertm(3,3), cbar_pertm
      real*8 Finv_pertm(3,3), J_pertm
      real*8 Tk_pertm(3,3)
	  
      real*8 Smat(6,1),Bmat(6,3*nNode),BodyForceRes(3*nNode,1),Qmat(9,9)
      real*8 Gmat(9,3*nNode),G0mat(9,3*nNode),Amat(9,9),wS(nIntS)
      real*8 xLocal(nIntS),yLocal(nIntS)
      real*8 Nvec(1,nNode),Fmat(6,1)
      real*8 junk1,junk2,junk3,junk4,junk5,junk6,junk7,junk8,junk9
      
      ! Properties
      real*8 Vmol,cbar0,cRmax,D,theta, Omega
      real*8 permitt,cbar_bound,z_bound,z_mob
      
      real*8 Eyoung,poisson,Gshear,Kbulk,Lambda

      ! State variables
      real*8 cbar_t,cbar_tau,cbardot
      !Outputs
      real*8 T_tau(3,3),pbar,detFs_tau,pbar_old
      

      ! Elemental shape functions and volume
      real*8 sh0(nNode),detMapJ0
      real*8 sh(nNode),detMapJ
      real*8 dsh(nNode,mcrd),detMapJC,detMapJ0C
      real*8 dshC(nNode,mcrd),xi(nInt,mcrd)
      real*8 dshxi(nNode,mcrd),dsh0(nNode,mcrd),dshC0(nNode,mcrd)
      real*8 w(nInt),ds,flux

      ! residuals and tangents
      real*8 Ru(mcrd*nNode,1),Romg(nNode,1),Rphi(nNode,1)
      real*8 Kuu(mcrd*nNode,mcrd*nNode),Komgomg(nNode,nNode)
      real*8 Kphiphi(nNode,nNode)
      real*8 Kuphi(mcrd*nNode,nNode),Kphiu(nNode,mcrd*nNode)
      real*8 Komgphi(nNode,nNode),Kphiomg(nNode,nNode)
      real*8 Komgu(nNode,mcrd*nNode),Kuomg(mcrd*nNode,nNode)
      
      
      ! BCs
      real*8 flux_FBV,J_lim, J_ox, J_red,alpha_o,alpha_r, n_el
      real*8 flux_FBV2

      real*8 zero,one,two,half,Pi,three,third,Farad,Rgas
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0,Rgas=8.314d0,Farad=96485.d0)
      
      


      ! Get element parameters
      !
      nlSdv = jprops(1) ! number of local sdv's per integ point
      ngSdv = jprops(2) ! number of global sdv's per integ point


      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
         !
         ! allocate memory for the globalSdv's
         !
         ! numElem needs to be set in the MODULE
         ! nInt needs to be set in the UEL
         !
         stat=0
c         allocate(globalSdv(numElem,nInt,ngSdv))
c         deallocate(globalSdv)
         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) 'error when allocating globalSdv'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) '   stat=',stat
            write(80,*) '  ngSdv=',ngSdv
            write(80,*) '   nInt=',nInt
            write(80,*) 'numElem=',numElem
            write(80,*) '  nNode=',nNode
            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
            write(80,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '---------- U3D8 ELEMENTS ------------------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt =',nInt
         write(*,*) '---------- nIntS=',nIntS
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv
         write(*,*) '-------------------------------------------------'
      endif


      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain initial conditions
      !
      Omega     = props(3)
      theta     = props(4)
      D         = props(6)
      z_mob     = props(7)
      cRmax     = props(8)
      cbar0     = props(9)
      permitt   = props(10)
      cbar_bound= props(11)
      z_bound   = props(12)
      
      
      ! elastic modulii 
      Eyoung         = props(1)
      poisson        = props(2)
      ! Calculate other elastic modulii
      !
      Gshear = Eyoung/(two*(one+poisson))
      Kbulk  = Eyoung/(three*(one-two*poisson))
      Lambda = Kbulk - two/third*Gshear
      

      ! Initialize the residual and tangent matrices to zero
      !
      Ru  = zero
      Romg= zero
      Rphi= zero
      
      Kphiphi = zero
      Kuu     = zero
      Komgomg = zero
      Kuomg   = zero
      Komgu   = zero
      Kphiomg = zero
      Komgphi = zero
      Kuphi   = zero
      Kphiu   = zero
      Energy  = zero


      ! Body forces
      !
      body(1:3) = zero


      ! Obtain nodal displacements and chemical potentials
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         k = k + 1
         phiNew(i) = Uall(k)
         dphi(i) = DUall(k,1)
         phiOld(i) = phiNew(i) - dphi(i)
         k = k + 1
         omgNew(i) = Uall(k)
         domg(i) = DUall(k,1)
         omgOld(i) = omgNew(i) - domg(i)
      enddo


      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


      ! Impose any time-stepping changes on the increments of
      !  chemical potential, temperature, or displacement if you wish
      !
      ! chemical potential increment
      !
      do i=1,nNode
         if(dabs(domg(i)).gt.1.d5) then
            pnewdt = 0.65
            return
         endif
      enddo
      !
      ! temperature increment
      !
      do i=1,nNode
         if(dabs(dphi(i)).gt.1.d0) then
            pnewdt = 0.6666
            return
         endif
      enddo
      !
      ! displacement increment, based on element diagonal
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,3))**two) + 
     +            ((coordsC(1,1)-coordsC(1,3))**two) + 
     +     ((coordsC(3,1)-coordsC(3,3))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.10.0*Le) then
               pnewdt = 0.55
               return
            endif
         enddo
      enddo


      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the
      !  `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures, 33, 3277-3296.
      !
      !
      ! Obtain shape functions and their local gradients at the element
      !  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.8) then
         call calcShape3DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         write(80,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape3D(nNode,dshxi,coords,dsh0,detMapJ0,stat,pnewdt)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      ! Map shape functions from local to global current coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape3D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat,pnewdt)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif



      ! Calculate the deformation gradient at the element centriod
      !  at the the begining and end of the increment for use in 
      !  the `F-bar' method. `Tau' represents the end of the increment
      !  and `t' the previous increment.
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      !
      call mdet(Fc_tau,detFc_tau)
      call mdet(Fc_t,detFc_t)
      !
      ! With the deformation gradient known at the element centriod
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over body integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.8) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.8) then
            call xint3D8pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint3D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         write(80,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! this is the first increment, of the first step
            !  give initial conditions
            !
            cbar_t  = cbar0
            pbar = zero
            !
         else
            !
            ! this is not the first increment, read old values
            !
            cbar_t  = svars(1+jj)
            pbar    = svars(2+jj)
            !
         endif


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.8) then
            call calcShape3DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.8'
            write(80,*) 'Incorrect number of nodes: nNode.ne.8'
            call xit
         endif
         

         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat,pnewdt)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! Map shape functions from local to global current coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape3D(nNode,dshxi,coordsC,dshC,detMapJC,stat,pnewdt)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! Obtain the electrochemical potential and its derivatives at 
         !  this intPt at the begining and end of the incrment
         !
         omg_tau = zero
         omg_t = zero
         domgdt = zero
         domgdx = zero
         do k=1,nNode
            omg_tau = omg_tau + omgNew(k)*sh(k)
            omg_t   = omg_t + omgOld(k)*sh(k)
            do i=1,nDim
               domgdx(i) = domgdx(i) + omgNew(k)*dshC(k,i)
            enddo
         enddo
         domgdt = (omg_tau - omg_t)/dtime


         ! Obtain the temperature and its derivative's at 
         !  this intPt at the begining and end of the incrment
         !
         phi_tau = zero
         phi_t = zero
         dphidt = zero
         dphidx = zero
         do k=1,nNode
            phi_tau = phi_tau + phiNew(k)*sh(k)
            phi_t   = phi_t + phiOld(k)*sh(k)
            do i=1,nDim
               dphidx(i) = dphidx(i) + phiNew(k)*dshC(k,i)
            enddo
         enddo
         dphidt = (phi_tau - phi_t)/dtime


         ! Obtain, and modify the deformation gradient at this integration
         !  point.  Modify the deformation gradienet for use in the `F-bar'
         !  method.  Also, take care of plane-strain or axisymetric
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 4 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.8).and.(nInt.eq.8)) then
            call mdet(F_tau,detF_tau)
            call mdet(F_t,detF_t)
            F_tau = ((detFc_tau/detF_tau)**third)*F_tau
            F_t = ((detFc_tau/detF_tau)**third)*F_t
         endif
         !
         call mdet(F_tau,detF)


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration at this integ. point 
         !
	     pbar_old = pbar		
         call integ(props,nprops,dtime,
     +        F_t, omg_t, phi_t, cbar_t,     
     +        F_tau,omg_tau,phi_tau,cbar_tau,
     +        -dphidx,
     +        T_tau,mob,pbar,detFs_tau,
     +        Tme_tau,Tes_tau,
     +        dmobdu,dmobdphi,dmobdomg,
     +        dcbardu,dcbardphi,dcbardomg
     +        )
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         svars(1+jj) = cbar_tau
         svars(2+jj) = pbar
         cbardot = (cbar_tau-cbar_t)/dtime
         jj = jj + nlSdv ! setup for the next intPt
      
         Tme_tau = T_tau - Tes_tau
		 
         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         !
         globalSdv(jelem,intPt,1) = cbar_tau   ! normalized species concentration
         globalSdv(jelem,intPt,2) = T_tau(1,1) !-\
         globalSdv(jelem,intPt,3) = T_tau(2,2) !  |-->   Stress
         globalSdv(jelem,intPt,4) = T_tau(1,2) !  |-->  Components
         globalSdv(jelem,intPt,5) = T_tau(3,3) !-/
		 globalSdv(jelem,intPt,6) = Tes_tau(1,1) !-\
         globalSdv(jelem,intPt,7) = Tes_tau(2,2) !  |-->  Electrostatic Stress
         globalSdv(jelem,intPt,8) = Tes_tau(1,2) !  |-->  Components
         globalSdv(jelem,intPt,9) = Tes_tau(3,3) !-/
		 globalSdv(jelem,intPt,10) = Tme_tau(1,1) !-\
         globalSdv(jelem,intPt,11) = Tme_tau(2,2) !  |-->  Mechanical Stress
         globalSdv(jelem,intPt,12) = Tme_tau(1,2) !  |-->  Components
         globalSdv(jelem,intPt,13) = Tme_tau(3,3) !-/
		 globalsdv(jelem,intPt,14) = -dphidx(1)
		 globalsdv(jelem,intPt,15) = -dphidx(2)
		 globalsdv(jelem,intPt,16) = -mob*domgdx(1)
		 globalsdv(jelem,intPt,17) = -mob*domgdx(2)
         globalSdv(jelem,intPt,18) = pbar
         globalSdv(jelem,intPt,19) = detFs_tau

         ! Time stepping algorithim based on the constitutive response.
         !  Here based on the change in the species concentration.
         !
         cbarLmt = 0.25d0
         umeror = dabs((cbar_tau - cbar_t)/cbarLmt)
         if(umeror.le.half) then
            pnewdt = 1.5d0
         elseif((umeror.gt.half).and.(umeror.le.0.8d0)) then
            pnewdt = 1.25d0
         elseif((umeror.gt.0.8d0).and.(umeror.le.1.25d0)) then
            pnewdt = 0.75d0
         else
            pnewdt = half
         endif


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration tangents 
         ! Finite difference based numerical evaluation of tangents

         ! First, perturbation in displacements
         Csp = zero 
         !call matInv3D(F_tau,Finv_tau,J_tau,stat)
         call m3inv(F_tau,Finv_tau)
         call mdet(F_tau,J_tau)
         
         delF= 1.d-6

         do l=1,3
          do k=1,3
		  
            F_pert = F_tau  
            F_pert(k,l) = F_pert(k,l) + delF   
            !call matInv3D(F_pert,Finv_pert,J_pert,stat)
            call m3inv(F_pert,Finv_pertp)
            call mdet(F_pert,J_pertp)
            pbar=pbar_old
            call integ(props,nprops,dtime,
     +        F_t, omg_t, phi_t, cbar_t,     
     +        F_pert,omg_tau,phi_tau,cbar_pertp,
     +        -dphidx,
     +        T_pertp,junk1,pbar,junk9,
     +        Tme_tau,Tes_tau,
     +        junk2,junk3,junk4,
     +        junk5,junk6,junk7
     +        )
	 
	        F_pert = F_tau  
            F_pert(k,l) = F_pert(k,l) - delF   
            !call matInv3D(F_pert,Finv_pert,J_pert,stat)
            call m3inv(F_pert,Finv_pertm)
            call mdet(F_pert,J_pertm)
            pbar=pbar_old
            call integ(props,nprops,dtime,
     +        F_t, omg_t, phi_t, cbar_t,     
     +        F_pert,omg_tau,phi_tau,cbar_pertm,
     +        -dphidx,
     +        T_pertm,junk1,pbar,junk9,
     +        Tme_tau,Tes_tau,
     +        junk2,junk3,junk4,
     +        junk5,junk6,junk7
     +        )


	        Tk_tau  = J_tau*matmul(T_tau,transpose(Finv_tau))
            Tk_pertp= J_pertp*matmul(T_pertp,transpose(Finv_pertp))
            Tk_pertm= J_pertm*matmul(T_pertm,transpose(Finv_pertm))
			
            Ccbartan(k,l) = (cbar_pertp - cbar_pertm)/(two*delF)
			

            do j=1,3
              do i=1,3
                Csp(i,j,k,l) =(Tk_pertp(i,j) - Tk_pertm(i,j))/(two*delF)
              enddo
            enddo
        
          enddo
        enddo


        SpTanMod = zero

        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                do m=1,3
                  do n=1,3
                      SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) + 
     +                                F_tau(j,m)*F_tau(l,n)*Csp(i,m,k,n)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      SpTanMod = SpTanMod/J_tau
      
      ! Update the tangent modulus
         !
!      SpTanMod = zero
!         do i=1,3
!            do j=1,3
!               do k=1,3
!                  do l=1,3
!                     SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) 
!     +                + Gshear*(Iden(i,k)*Iden(j,l)+Iden(i,l)*Iden(j,k)
!     +                            -(two/three)*Iden(i,j)*Iden(k,l))
!     +                + Kbulk*Iden(i,j)*Iden(k,l)
!                  enddo
!               enddo
!            enddo
!         enddo


        ! Next, perturbation in omg

          delomg = omg_tau*1.d-3+1.d-8
          omg_pert = omg_tau+delomg
	      pbar=pbar_old	   
          call integ(props,nprops,dtime,
     +        F_t, omg_t, phi_t, cbar_t,     
     +        F_tau,omg_pert,phi_tau,cbar_pertp,
     +        -dphidx,
     +        T_pertp,junk1,pbar,junk9,
     +        Tme_tau,Tes_tau,
     +        junk2,junk3,junk4,
     +        junk5,junk6,junk7
     +        )

          delomg = omg_tau*1.d-3+1.d-8
          omg_pert = omg_tau-delomg
	      pbar=pbar_old	   
          call integ(props,nprops,dtime,
     +        F_t, omg_t, phi_t, cbar_t,     
     +        F_tau,omg_pert,phi_tau,cbar_pertm,
     +        -dphidx,
     +        T_pertm,junk1,pbar,junk9,
     +        Tme_tau,Tes_tau,
     +        junk2,junk3,junk4,
     +        junk5,junk6,junk7
     +        )
	 
          Spuomgmod = (T_pertp - T_pertm)/(two*delomg)
		  dcbardomg = (cbar_pertp-cbar_pertm)/(two*delomg)
		  
        ! Next, perturbation in phi

          delphi = phi_tau*1.d-3+1.d-8
          phi_pert = phi_tau+delphi
		  pbar=pbar_old   
          call integ(props,nprops,dtime,
     +        F_t, omg_t, phi_t, cbar_t,     
     +        F_tau,omg_tau,phi_pert,cbar_pertp,
     +        -dphidx,
     +        T_pertp,junk1,pbar,junk9,
     +        Tme_tau,Tes_tau,
     +        junk2,junk3,junk4,
     +        junk5,junk6,junk7
     +        )
	 
	      delphi = phi_tau*1.d-3+1.d-8
          phi_pert = phi_tau-delphi
          pbar=pbar_old		   
          call integ(props,nprops,dtime,
     +        F_t, omg_t, phi_t, cbar_t,     
     +        F_tau,omg_tau,phi_pert,cbar_pertm,
     +        -dphidx,
     +        T_pertm,junk1,pbar,junk9,
     +        Tme_tau,Tes_tau,
     +        junk2,junk3,junk4,
     +        junk5,junk6,junk7
     +        )

          Spuphimod = (T_pertp - T_pertm)/(two*delphi)
	      dcbardphi = (cbar_pertp-cbar_pertm)/(two*delphi)
            
      dmobdomg = D*cRmax*dcbardomg*(one-two*cbar_tau)/(Rgas*theta*detF)
      dmobdphi = D*cRmax*dcbardphi*(one-two*cbar_tau)/(Rgas*theta*detF)
		 !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



         ! Compute/update the displacement residual vector
         !
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(3,3)
         Smat(4,1) = T_tau(1,2)
         Smat(5,1) = T_tau(2,3)
         Smat(6,1) = T_tau(1,3)
         do kk=1,nNode
            Bmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Bmat(2,2+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,3+nDim*(kk-1)) = dshC(kk,3)
            Bmat(4,1+nDim*(kk-1)) = dshC(kk,2)
            Bmat(4,2+nDim*(kk-1)) = dshC(kk,1)
            Bmat(5,2+nDim*(kk-1)) = dshC(kk,3)
            Bmat(5,3+nDim*(kk-1)) = dshC(kk,2)
            Bmat(6,1+nDim*(kk-1)) = dshC(kk,3)
            Bmat(6,3+nDim*(kk-1)) = dshC(kk,1)
         enddo
         !
         BodyForceRes = zero
         do kk=1,nNode
            BodyForceRes(1+nDim*(kk-1),1) = sh(kk)*body(1)
            BodyForceRes(2+nDim*(kk-1),1) = sh(kk)*body(2)
            BodyForceRes(3+nDim*(kk-1),1) = sh(kk)*body(3)
         enddo
         !
         Ru = Ru + detmapJC*w(intpt)*
     +        (
     +        -matmul(transpose(Bmat),Smat)
     +        + BodyForceRes
     +        )



         ! Compute/update the chemical potential residual vector
         !
         do i=1,nNode
            Romg(i,1) = Romg(i,1) + detmapJC*w(intpt)*
     +        (
     +        mob*(  dshC(i,1)*domgdx(1) + dshC(i,2)*domgdx(2)
     +             + dshC(i,3)*domgdx(3)  )
     +       + sh(i)*cRmax*cbardot/detF_tau
     +        )
         enddo



         ! Compute/update the temperature residual vector
         !
         qNet = (z_bound*cbar_bound+z_mob*cbar_tau)*Farad*cRmax/detF_tau

         do i=1,nNode
            Rphi(i,1) = Rphi(i,1) + detmapJC*w(intpt)*
     +        (
     +       -(   dshC(i,1)*dphidx(1) + dshC(i,2)*dphidx(2)
     +          + dshC(i,3)*dphidx(3)   )
     +       + sh(i)*qNet/permitt
     +        )
         enddo



         ! Compute/update the displacement tangent matrix
         !
         Gmat = zero
         do kk=1,nNode
            Gmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Gmat(2,2+nDim*(kk-1)) = dshC(kk,1)
            Gmat(3,3+nDim*(kk-1)) = dshC(kk,1)
            Gmat(4,1+nDim*(kk-1)) = dshC(kk,2)
            Gmat(5,2+nDim*(kk-1)) = dshC(kk,2)
            Gmat(6,3+nDim*(kk-1)) = dshC(kk,2)
            Gmat(7,1+nDim*(kk-1)) = dshC(kk,3)
            Gmat(8,2+nDim*(kk-1)) = dshC(kk,3)
            Gmat(9,3+nDim*(kk-1)) = dshC(kk,3)
         enddo

         G0mat = zero
         do kk=1,nNode
            G0mat(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(3,3+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(4,1+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(5,2+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(6,3+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(7,1+nDim*(kk-1)) = dshC0(kk,3)
            G0mat(8,2+nDim*(kk-1)) = dshC0(kk,3)
            G0mat(9,3+nDim*(kk-1)) = dshC0(kk,3)
         enddo

         Amat = zero
         Amat(1,1) = SpTanMod(1,1,1,1)
         Amat(1,2) = SpTanMod(1,1,2,1)
         Amat(1,3) = SpTanMod(1,1,3,1)
         Amat(1,4) = SpTanMod(1,1,1,2)
         Amat(1,5) = SpTanMod(1,1,2,2)
         Amat(1,6) = SpTanMod(1,1,3,2)
         Amat(1,7) = SpTanMod(1,1,1,3)
         Amat(1,8) = SpTanMod(1,1,2,3)
         Amat(1,9) = SpTanMod(1,1,3,3)
         Amat(2,1) = SpTanMod(2,1,1,1)
         Amat(2,2) = SpTanMod(2,1,2,1)
         Amat(2,3) = SpTanMod(2,1,3,1)
         Amat(2,4) = SpTanMod(2,1,1,2)
         Amat(2,5) = SpTanMod(2,1,2,2)
         Amat(2,6) = SpTanMod(2,1,3,2)
         Amat(2,7) = SpTanMod(2,1,1,3)
         Amat(2,8) = SpTanMod(2,1,2,3)
         Amat(2,9) = SpTanMod(2,1,3,3)
         Amat(3,1) = SpTanMod(3,1,1,1)
         Amat(3,2) = SpTanMod(3,1,2,1)
         Amat(3,3) = SpTanMod(3,1,3,1)
         Amat(3,4) = SpTanMod(3,1,1,2)
         Amat(3,5) = SpTanMod(3,1,2,2)
         Amat(3,6) = SpTanMod(3,1,3,2)
         Amat(3,7) = SpTanMod(3,1,1,3)
         Amat(3,8) = SpTanMod(3,1,2,3)
         Amat(3,9) = SpTanMod(3,1,3,3)
         Amat(4,1) = SpTanMod(1,2,1,1)
         Amat(4,2) = SpTanMod(1,2,2,1)
         Amat(4,3) = SpTanMod(1,2,3,1)
         Amat(4,4) = SpTanMod(1,2,1,2)
         Amat(4,5) = SpTanMod(1,2,2,2)
         Amat(4,6) = SpTanMod(1,2,3,2)
         Amat(4,7) = SpTanMod(1,2,1,3)
         Amat(4,8) = SpTanMod(1,2,2,3)
         Amat(4,9) = SpTanMod(1,2,3,3)
         Amat(5,1) = SpTanMod(2,2,1,1)
         Amat(5,2) = SpTanMod(2,2,2,1)
         Amat(5,3) = SpTanMod(2,2,3,1)
         Amat(5,4) = SpTanMod(2,2,1,2)
         Amat(5,5) = SpTanMod(2,2,2,2)
         Amat(5,6) = SpTanMod(2,2,3,2)
         Amat(5,7) = SpTanMod(2,2,1,3)
         Amat(5,8) = SpTanMod(2,2,2,3)
         Amat(5,9) = SpTanMod(2,2,3,3)
         Amat(6,1) = SpTanMod(3,2,1,1)
         Amat(6,2) = SpTanMod(3,2,2,1)
         Amat(6,3) = SpTanMod(3,2,3,1)
         Amat(6,4) = SpTanMod(3,2,1,2)
         Amat(6,5) = SpTanMod(3,2,2,2)
         Amat(6,6) = SpTanMod(3,2,3,2)
         Amat(6,7) = SpTanMod(3,2,1,3)
         Amat(6,8) = SpTanMod(3,2,2,3)
         Amat(6,9) = SpTanMod(3,2,3,3)
         Amat(7,1) = SpTanMod(1,3,1,1)
         Amat(7,2) = SpTanMod(1,3,2,1)
         Amat(7,3) = SpTanMod(1,3,3,1)
         Amat(7,4) = SpTanMod(1,3,1,2)
         Amat(7,5) = SpTanMod(1,3,2,2)
         Amat(7,6) = SpTanMod(1,3,3,2)
         Amat(7,7) = SpTanMod(1,3,1,3)
         Amat(7,8) = SpTanMod(1,3,2,3)
         Amat(7,9) = SpTanMod(1,3,3,3)
         Amat(8,1) = SpTanMod(2,3,1,1)
         Amat(8,2) = SpTanMod(2,3,2,1)
         Amat(8,3) = SpTanMod(2,3,3,1)
         Amat(8,4) = SpTanMod(2,3,1,2)
         Amat(8,5) = SpTanMod(2,3,2,2)
         Amat(8,6) = SpTanMod(2,3,3,2)
         Amat(8,7) = SpTanMod(2,3,1,3)
         Amat(8,8) = SpTanMod(2,3,2,3)
         Amat(8,9) = SpTanMod(2,3,3,3)
         Amat(9,1) = SpTanMod(3,3,1,1)
         Amat(9,2) = SpTanMod(3,3,2,1)
         Amat(9,3) = SpTanMod(3,3,3,1)
         Amat(9,4) = SpTanMod(3,3,1,2)
         Amat(9,5) = SpTanMod(3,3,2,2)
         Amat(9,6) = SpTanMod(3,3,3,2)
         Amat(9,7) = SpTanMod(3,3,1,3)
         Amat(9,8) = SpTanMod(3,3,2,3)
         Amat(9,9) = SpTanMod(3,3,3,3)


         Qmat = zero
         Qmat(1,1) = third*(Amat(1,1)+Amat(1,5)+Amat(1,9)) 
     +        - (two/three)*T_tau(1,1)
         Qmat(2,1) = third*(Amat(2,1)+Amat(2,5)+Amat(2,9))
     +        - (two/three)*T_tau(2,1)
         Qmat(3,1) = third*(Amat(3,1)+Amat(3,5)+Amat(3,9))
     +        - (two/three)*T_tau(3,1)
         Qmat(4,1) = third*(Amat(4,1)+Amat(4,5)+Amat(4,9))
     +        - (two/three)*T_tau(1,2)
         Qmat(5,1) = third*(Amat(5,1)+Amat(5,5)+Amat(5,9))
     +        - (two/three)*T_tau(2,2)
         Qmat(6,1) = third*(Amat(6,1)+Amat(6,5)+Amat(6,9))
     +        - (two/three)*T_tau(3,2)
         Qmat(7,1) = third*(Amat(7,1)+Amat(7,5)+Amat(7,9))
     +        - (two/three)*T_tau(1,3)
         Qmat(8,1) = third*(Amat(8,1)+Amat(8,5)+Amat(8,9))
     +        - (two/three)*T_tau(2,3)
         Qmat(9,1) = third*(Amat(9,1)+Amat(9,5)+Amat(9,9))
     +        - (two/three)*T_tau(3,3)
         Qmat(1,5) = Qmat(1,1)
         Qmat(2,5) = Qmat(2,1)
         Qmat(3,5) = Qmat(3,1)
         Qmat(4,5) = Qmat(4,1)
         Qmat(5,5) = Qmat(5,1)
         Qmat(6,5) = Qmat(6,1)
         Qmat(7,5) = Qmat(7,1)
         Qmat(8,5) = Qmat(8,1)
         Qmat(9,5) = Qmat(9,1)
         Qmat(1,9) = Qmat(1,1)
         Qmat(2,9) = Qmat(2,1)
         Qmat(3,9) = Qmat(3,1)
         Qmat(4,9) = Qmat(4,1)
         Qmat(5,9) = Qmat(5,1)
         Qmat(6,9) = Qmat(6,1)
         Qmat(7,9) = Qmat(7,1)
         Qmat(8,9) = Qmat(8,1)
         Qmat(9,9) = Qmat(9,1)
            

         if((nNode.eq.8).and.(nInt.eq.8)) then
            !
            ! This is the tangent using the F-bar method with the
            !  4 node fully integrated linear element
            !
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +            matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           + matmul(transpose(Gmat),matmul(Qmat,(G0mat-Gmat)))
     +           )
         else
            !
            ! This is the tangent not using the F-bar method with all
            !  other elements
            !
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           )
         endif


         ! Compute/update the electrochemical potential tangent matrix
         !
         do i=1,nNode
          do j=1,nNode

            Komgomg(i,j) = Komgomg(i,j) + detmapJC*w(intPt)*
     +        (
     +        - dmobdomg*sh(j)*( dshC(i,1)*domgdx(1)+dshC(i,2)*domgdx(2)
     +                          +dshC(i,2)*domgdx(2)  )
     +        - mob*(dshC(i,1)*dshC(j,1) + dshC(i,2)*dshC(j,2) + 
     +               dshC(i,3)*dshC(j,3))
     +        - sh(i)*sh(j)*cRmax/detF_tau*dcbardomg/dtime
     +        )
          enddo
         enddo


         ! Compute/update the electric potential tangent matrix
         !
         do i=1,nNode
          do j=1,nNode

            Kphiphi(i,j) = Kphiphi(i,j) + detmapJC*w(intPt)*
     +        (
     +    (dshC(i,1)*dshC(j,1)+dshC(i,2)*dshC(j,2)+dshC(i,3)*dshC(j,3))
     +      - sh(i)*sh(j)*(z_mob*Farad)/permitt*dcbardphi*cRmax/detF_tau
     +        )
          enddo
         enddo


         ! Compute/update the chemical potential - temperature tangent matrix.
         !
         do i=1,nNode
          do j=1,nNode

            Komgphi(i,j) = Komgphi(i,j) + detmapJC*w(intPt)*
     +        (
     +      -dmobdphi*sh(j)*( dshC(i,1)*domgdx(1) + dshC(i,2)*domgdx(2)
     +                      + dshC(i,3)*domgdx(3)  )
     +        - sh(i)*sh(j)*cRmax/detF_tau*dcbardphi/dtime
     +        )
          enddo
         enddo


         ! Compute/update the temperature - chemical potential tangent matrix.
         !
         do i=1,nNode
          do j=1,nNode

            Kphiomg(i,j) = Kphiomg(i,j) + detmapJC*w(intPt)*
     +        (
     +      -sh(i)*sh(j)*(z_mob*Farad)/permitt*dcbardomg*cRmax/detF_tau
     +        )
          enddo
         enddo
         
         
         ! do i=1,nNode
          ! do j=1,nNode         
            
            ! l=nDim*(j-1)+1
            ! Komgu(i,l) = Komgu(i,l) + detmapJC*w(intPt)*
      ! +        (
      ! +   -D/Rgas/theta*cRmax*(dshC(i,1)*domgdx(1) + dshC(i,2)*domgdx(2))
      ! +        - sh(i)*cRmax/detF_tau/dtime
      ! +        )*(Ccbartan(1,1)*dshC(j,1)+Ccbartan(1,2)*dshC(j,2))
            
            ! Komgu(i,l+1) = Komgu(i,l+1) + detmapJC*w(intPt)*
      ! +        (
      ! +   -D/Rgas/theta*cRmax*(dshC(i,1)*domgdx(1) + dshC(i,2)*domgdx(2))
      ! +        - sh(i)*cRmax/detF_tau/dtime
      ! +        )*(Ccbartan(2,1)*dshC(j,1)+Ccbartan(2,2)*dshC(j,2))
            
          ! enddo
      ! enddo

	  	  SpomgUMod = zero
         do i=1,nDim
            do k=1,nDim
               do l=1,nDim
                  SpomgUMod(i,k,l) = SpomgUMod(i,k,l)
     +                 + domgdx(k)*mob*iden(i,l)
               enddo
            enddo
         enddo
         !
         AmatCU = zero
         AmatCU(1,1) = SpomgUMod(1,1,1)
         AmatCU(1,2) = SpomgUMod(1,2,1)
         AmatCU(1,3) = SpomgUMod(1,3,1)
         AmatCU(1,4) = SpomgUMod(1,1,2)
         AmatCU(1,5) = SpomgUMod(1,2,2)
         AmatCU(1,6) = SpomgUMod(1,3,2)
         AmatCU(1,7) = SpomgUMod(1,1,3)
         AmatCU(1,8) = SpomgUMod(1,2,3)
         AmatCU(1,9) = SpomgUMod(1,3,3)
         AmatCU(2,1) = SpomgUMod(2,1,1)
         AmatCU(2,2) = SpomgUMod(2,2,1)
         AmatCU(2,3) = SpomgUMod(2,3,1)
         AmatCU(2,4) = SpomgUMod(2,1,2)
         AmatCU(2,5) = SpomgUMod(2,2,2)
         AmatCU(2,6) = SpomgUMod(2,3,2)
         AmatCU(2,7) = SpomgUMod(2,1,3)
         AmatCU(2,8) = SpomgUMod(2,2,3)
         AmatCU(2,9) = SpomgUMod(2,3,3)
         AmatCU(3,1) = SpomgUMod(3,1,1)
         AmatCU(3,2) = SpomgUMod(3,2,1)
         AmatCU(3,3) = SpomgUMod(3,3,1)
         AmatCU(3,4) = SpomgUMod(3,1,2)
         AmatCU(3,5) = SpomgUMod(3,2,2)
         AmatCU(3,6) = SpomgUMod(3,3,2)
         AmatCU(3,7) = SpomgUMod(3,1,3)
         AmatCU(3,8) = SpomgUMod(3,2,3)
         AmatCU(3,9) = SpomgUMod(3,3,3)
         
         !
         Komgu = Komgu - detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(dshC,AmatCU),Gmat)
     +        )
         ! Compute/update the displacement - electrochemical potential tangent matrix
         !  The F-bar method will have some effect, however we neglect that here.
         !
         Fmat = zero
         Fmat(1,1) = SpuomgMod(1,1)
         Fmat(2,1) = SpuomgMod(2,2)
         Fmat(3,1) = SpuomgMod(3,3)
         Fmat(4,1) = SpuomgMod(1,2)
         Fmat(5,1) = SpuomgMod(2,3)
         Fmat(6,1) = SpuomgMod(1,3)

         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo

         Kuomg = Kuomg + detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(transpose(Bmat),Fmat),Nvec)
     +        )



         ! Compute/update the displacement - electric potential tangent matrix.
         !  The F-bar method will have some effect, however we neglect that here.
         !
         Fmat = zero
         Fmat(1,1) = SpuphiMod(1,1)
         Fmat(2,1) = SpuphiMod(2,2)
         Fmat(3,1) = SpuphiMod(3,3)
         Fmat(4,1) = SpuphiMod(1,2)
         Fmat(5,1) = SpuphiMod(2,3)
         Fmat(6,1) = SpuphiMod(1,3)
         !
         Kuphi = Kuphi + detMapJC*w(intpt)*
     +        (
     +        zero+
     +        matmul(matmul(transpose(Bmat),Fmat),Nvec)
     +        )


	  	  SpphiUMod = zero
         do i=1,nDim
            do k=1,nDim
               do l=1,nDim
                  SpphiUMod(i,k,l) = SpphiUMod(i,k,l)
     +                 + dphidx(k)*iden(i,l)
               enddo
            enddo
         enddo
         !
         AmatCU = zero
         AmatCU(1,1) = SpphiUMod(1,1,1)
         AmatCU(1,2) = SpphiUMod(1,2,1)
         AmatCU(1,3) = SpphiUMod(1,3,1)
         AmatCU(1,4) = SpphiUMod(1,1,2)
         AmatCU(1,5) = SpphiUMod(1,2,2)
         AmatCU(1,6) = SpphiUMod(1,3,2)
         AmatCU(1,7) = SpphiUMod(1,1,3)
         AmatCU(1,8) = SpphiUMod(1,2,3)
         AmatCU(1,9) = SpphiUMod(1,3,3)
         AmatCU(2,1) = SpphiUMod(2,1,1)
         AmatCU(2,2) = SpphiUMod(2,2,1)
         AmatCU(2,3) = SpphiUMod(2,3,1)
         AmatCU(2,4) = SpphiUMod(2,1,2)
         AmatCU(2,5) = SpphiUMod(2,2,2)
         AmatCU(2,6) = SpphiUMod(2,3,2)
         AmatCU(2,7) = SpphiUMod(2,1,3)
         AmatCU(2,8) = SpphiUMod(2,2,3)
         AmatCU(2,9) = SpphiUMod(2,3,3)
         AmatCU(3,1) = SpphiUMod(3,1,1)
         AmatCU(3,2) = SpphiUMod(3,2,1)
         AmatCU(3,3) = SpphiUMod(3,3,1)
         AmatCU(3,4) = SpphiUMod(3,1,2)
         AmatCU(3,5) = SpphiUMod(3,2,2)
         AmatCU(3,6) = SpphiUMod(3,3,2)
         AmatCU(3,7) = SpphiUMod(3,1,3)
         AmatCU(3,8) = SpphiUMod(3,2,3)
         AmatCU(3,9) = SpphiUMod(3,3,3)
         
         !
         Kphiu = Kphiu - detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(dshC,AmatCU),Gmat)
     +        )
         ! do i=1,nNode
          ! do j=1,nNode
			
			! l=nDim*(j-1)+1
						
            ! Kphiu(i,l) = Kphiu(i,l) + detmapJC*w(intPt)*
      ! +        (
      ! +        - sh(i)*(z_mob*Farad)/permitt*cRmax/detF_tau
      ! +        )*(Ccbartan(1,1)*dsh(j,1)+Ccbartan(1,2)*dsh(j,2))

            ! Kphiu(i,l+1) = Kphiu(i,l+1) + detmapJC*w(intPt)*
      ! +        (
      ! +        - sh(i)*(z_mob*Farad)/permitt*cRmax/detF_tau
      ! +        )*(Ccbartan(2,1)*dsh(j,1)+Ccbartan(2,2)*dsh(j,2))
	 
          ! enddo
         ! enddo

            enddo
      !
      ! End the loop over body integration points
      !----------------------------------------------------------------
      
       !Rphi =zero
       !Romg =zero
       !Komgphi = zero
       !Kphiomg=zero
       !Kuomg =zero
       !Kuphi=zero
       !Kphiu=zero
       !Komgu=zero
       !Kphiphi=one
       !Komgomg=one
      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,Romg,Rphi,
     +     Kuu,Kuomg,Kuphi,Komgu,Komgomg,Komgphi,Kphiu,Kphiomg,Kphiphi,
     +     rhs,amatrx)
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine U3D8

************************************************************************
************************************************************************

      subroutine integ(props,nprops,dtime,
     +        F_t, omg_t, phi_t, cbar_t,     
     +        F_tau,omg_tau,phi_tau,cbar_tau,
     +        e_tau,
     +        T_tau,mob,pbar,detFs,
     +        Tme_tau,T_max_tau,
     +        dmobdu,dmobdphi,dmobdomg,
     +        dcbardu,dcbardphi,dcbardomg
     +        )


      ! This subroutine computes everything required for the time integration
      ! of the problem.
      !
      ! Inputs:
      !  1) material parameters, props(nprops)
      !  2) time increment, dtime
      !  3) deformation gradient, F_tau(3,3), F_t(3,3)
      !  4) electrochemical potential, omg_tau, omg_t
      !  5) electric potential, phi_tau, phi_t
      !  6) electric field,  e_tau
      !
      ! Outputs:
      !  1) Cauchy stress, T_tau(3,3)
      !  2) spatial tangent modulus, SpTanMod(3,3,3,3)
      !  3) species concentration, cbar_tau
      !  4) mobility, mob
      !  5) derivative of mob with u, omg, phi : dmobdu,dmobdphi,dmobdomg
      !  6) derivative of cbar with u, omg, phi : dcbardu,dcbardphi,dcbardomg
      !  7) derivatives of deformation residual :  SpTanMod,Spuomgmod,Spuphimod

      implicit none

      integer i,j,k,l,m,n,nprops,nargs,stat
      parameter(nargs=13)

      ! Inputs
      real*8 props(nprops),dtime
      real*8 F_tau(3,3),phi_tau,omg_tau,e_tau(3)
      real*8 F_t(3,3),phi_t,omg_t
      ! Outputs
      real*8 cbar_tau,cbar_t,T_tau(3,3),mob
      real*8 dmobdu(2),dmobdphi,dmobdomg
      real*8 dcbardu(2),dcbardphi,dcbardomg
      real*8 SpTanMod(3,3,3,3),Spuomgmod(3,3),Spuphimod(3,3)


      ! Properties
      real*8 Eyoung,anu,Gshear,Kbulk
      real*8 Omega,theta,D,mu0,z_mob,cRmax,cbar0,permitt, permitt_real


      ! In-subroutine use
      real*8 Iden(3,3)
      real*8 Finv(3,3),detF,FinvT(3,3)
      real*8 e_tau_mod
      real*8 T_max_tau(3,3),trT_max, T_ion_tau(3,3)
      real*8 Tme_tau(3,3)
      real*8 detFs,detFe, lambda_s
      real*8 Fs_tau(3,3),Fe_tau(3,3)
      real*8 Re_tau(3,3),Ue_tau(3,3),Ee_tau(3,3),trEe_tau,Ee0_tau(3,3)
      real*8 Me_tau(3,3)
      real*8 args(nargs)
      real*8 phi_per,omg_per,cbar_per
      real*8 pbar

      real*8 zero,one,two,three,third,half,Rgas,Farad
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0,Rgas = 8.314d0, Farad=96485.d0)


      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain material properties
      !
      Eyoung = props(1)
      anu    = props(2)
      Omega  = props(3)
      theta  = props(4)
      mu0    = props(5)
      D      = props(6)
      z_mob  = props(7)
      cRmax  = props(8)
      cbar0  = props(9)
      permitt= props(10)
	  !
	  permitt_real = props(13)
      !
      Gshear = Eyoung/2/(one+anu)
      Kbulk = Eyoung/three/(one-two*anu)
      
      



      ! Compute the inverse of F, its determinant, and its transpose
      !
      call matInv3D(F_tau,Finv,detF,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero'
         call xit
      endif
      FinvT = transpose(Finv)


      e_tau_mod = e_tau(1)*e_tau(1)+e_tau(2)*e_tau(2)+e_tau(3)*e_tau(3)
      T_max_tau = zero
      do i=1,3
        do j=1,3
            T_max_tau(i,j) = permitt_real*( e_tau(i)*e_tau(j)
     +                            -half*(e_tau_mod)*Iden(i,j) )
        enddo
      enddo
      !T_max_tau = zero
      trT_max = T_max_tau(1,1)+T_max_tau(2,2)+T_max_tau(3,3)
      
      ! Compute the species concentration
      !
      args(1)  = omg_tau
      args(2)  = mu0
      args(3)  = Rgas
      args(4)  = theta
      args(5)  = Omega
      args(6)  = Kbulk
      args(7)  = detF
      args(8)  = trT_max
      args(9)  = z_mob
      args(10) = Farad
      args(11) = cRmax
      args(12) = cbar0
      args(13) = phi_tau  
     
!      call solveConc(cbar_tau,args,nargs,cbar_t)
	  

      cbar_tau = dexp((omg_tau-mu0-z_mob*Farad*phi_tau-Omega*pbar)/
!      cbar_tau = dexp((omg_tau-mu0-z_mob*Farad*phi_tau)/
     +                                 (Rgas*theta))
      !
	  cbar_tau = cbar_tau/(one+cbar_tau)
	  !
	  !if (cbar_tau.gt.0.985) then
	  !cbar_tau = 0.985
	  !end if 
	 
      
      ! Compute the swelling volume ratio, detFs
      !
      detFs = 1 + Omega*cRmax*(cbar_tau-cbar0)
      ! Compute the elastic volume ratio, detFe
      !
      detFe = detF/detFs


      
      !dcbardomg = cbar_tau/Rgas/theta
      mob = D*cRmax*cbar_tau*(one-cbar_tau)/(Rgas*theta*detF)
!      dmobdomg = D*cRmax*dcbardomg*()/(Rgas*theta*detF)
!      dmobdphi = D*cRmax*dcbardphi/(Rgas*theta*detF)



      lambda_s = detFs**third
      Fs_tau = lambda_s*Iden
      Fe_tau = lambda_s**(-one)*F_tau

      call skinem(Fe_tau,Re_tau,Ue_tau,Ee_tau)
      trEe_tau = Ee_tau(1,1) + Ee_tau(2,2) + Ee_tau(3,3)
      Ee0_tau = Ee_tau - trEe_tau*third*Iden
      
      
      
      Me_tau = two*Gshear*Ee0_tau + Kbulk*trEe_tau*Iden

      ! Compute the Cauchy stress
      !
      T_tau = matmul(Re_tau,matmul(Me_tau,transpose(Re_tau)))/detFe 
     +                  + T_max_tau
      
      pbar = -third*
     +  ( T_tau(1,1)+T_tau(2,2)+T_tau(3,3) )
      
      return
      end subroutine integ

****************************************************************************

      subroutine solveConc(root,args,nargs,rootOld)

      ! This subroutine will numerically solve for the species conentration
      ! based on the current state
      ! See numerical recipies RTSAFE.

      implicit none

      integer maxit,j,nargs

      real*8 xacc,f,df,fl,fh,xl,xh,x1,x2,swap,root,dxold,one
      real*8 dx,args(nargs),zero,rootOld,temp,rootMax,rootMin

      parameter(maxit=100)
      parameter(xacc=1.d-5,zero=0.d0,one=1.d0)

      
      rootMax = 1.d2 ! 
      rootMin = 1.d-6 ! 

      x1 = rootMin
      x2 = rootMax
      call concFunc(x1,FL,DF,args,nargs)
      
      call concFunc(x2,FH,DF,args,nargs)

      if(fl*fh.ge.zero) then
         root = rootOld
         write(*,*) 'FYI, root not bracketed on conc'
         write(*,*) 'fl=',fl
         write(*,*) 'fh=',fh
         write(*,*) 'rootOld=',rootOld
         write(80,*) 'FYI, the root is not bracketed on conc'
         write(80,*) 'fl=',fl
         write(80,*) 'fh=',fh
         write(80,*) 'rootOld=',rootOld
         do j=1,nargs
         write(*,*)'args-',j,args(j)
         enddo
         call xit
         return
      endif

C
C		ORIENT THE SEARCH SO THAT F(XL) < 0.
C
      IF( FL .LT. 0.D0 ) THEN
         XL = X1
         XH = X2
      ELSE
         XH = X1
         XL = X2
         SWAP = FL
         FL = FH
         FH = SWAP
      END IF
C
C		INITIALIZE THE GUESS FOR THE ROOT, THE ''STEP SIZE
C		BEFORE LAST'', AND THE LAST STEP
C
      if(rootOld.lt.rootMin) rootOld = rootMin
      if(rootOld.gt.rootMax) rootOld = rootMax
      ROOT = rootOld !0.5D0 *( X1 + X2)
      DXOLD = DABS(X2 - X1)
      DX = DXOLD
      
      call concFunc(root,F,DF,args,nargs)

C
C			LOOP OVER ALLOWED ITERATIONS
C
      DO 10 J = 1,MAXIT
C
C			BISECT IF NEWTON OUT OF RANGE, OR NOT DECREASING
C			FAST ENOUGH.
C
         IF( ((ROOT-XH)*DF - F)*((ROOT - XL)*DF -F) .GE. 0.D0
     +        .OR. DABS(2.D0*F) .GT. DABS(DXOLD*DF) ) THEN

            DXOLD = DX
            DX = 0.5D0*(XH-XL)
            ROOT = XL + DX
            IF( XL .EQ. ROOT ) THEN
C
C			CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         ELSE
C
C			NEWTON STEP IS ACCEPTABLE. TAKE IT.
C
            DXOLD = DX
            DX = F/DF
            TEMP = ROOT
            ROOT = ROOT - DX
            IF( TEMP .EQ. ROOT) THEN
C
C			 CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         END IF
C
C		CONVERVEGENCE CRITERION
C
 
        IF( DABS(DX) .LT. XACC) RETURN

C
C			THE ONE NEW FUNCTION EVALUATION PER ITERATION
C
         call concFunc(root,F,DF,args,nargs)

C
C		MAINTAIN THE BRACKET ON THE ROOT
C
         IF( F .LT. 0.D0) THEN
            XL = ROOT
            FL = F
         ELSE
            XH = ROOT
            FH = F
         END IF

 10   CONTINUE

      WRITE(*,'(/1X,A)') 'solveConc EXCEEDING MAXIMUM ITERATIONS'
      WRITE(80,'(/1X,A)') 'solveConc EXCEEDING MAXIMUM ITERATIONS'

      return
      end subroutine solveConc

****************************************************************************

      subroutine concFunc(cbar,f,df,args,nargs)

      ! This subroutine serves as the function we would like to solve for
      !  the polymer volume fraction by finding phi such that ``f=0''

      implicit none

      integer nargs,NeoHookean,Langevin,material
      parameter(NeoHookean=1,Langevin=2)

      real*8 args(nargs),f,df,cbar
      real*8 omg_tau,mu0,Rgas,theta,Omega,Kbulk,detF,trT_max,z_mob,Farad
      real*8 cRmax,cbar0,phi_tau

      real*8 zero,one,two,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0)

      
      ! Obtain relevant quantities
      !
      omg_tau = args(1)  
      mu0     = args(2)
      Rgas    = args(3)    
      theta   = args(4)
      Omega   = args(5)
      Kbulk   = args(6)
      detF    = args(7)
      trT_max = args(8)
      z_mob   = args(9)
      Farad   = args(10)
      cRmax   = args(11)
      cbar0   = args(12)
      phi_tau = args(13)
 
	  
      ! Compute the residual
      !
      f = omg_tau - mu0 - Rgas*theta*dlog(cbar/detF) 
     +    +Omega* (
     +     Kbulk*dlog(detF/(1+Omega*cRmax*(cbar-cbar0))) 
     +   +third*trT_max*detF/(1+Omega*cRmax*(cbar-cbar0)) )
     + 	  -z_mob*Farad*phi_tau 

      ! Compute the tangent
      !
      df = - Rgas*theta/cbar 
     +   -Omega*third*(trT_max)*detF
     +                             /(1+Omega*cRmax*(cbar-cbar0))**two
     +                *(Omega*cRmax)
     +   -Omega*Kbulk/(1+Omega*cRmax*(cbar-cbar0))*Omega*cRmax    


      return
      end subroutine concFunc

************************************************************************

      subroutine AssembleElement(nDim,nNode,ndofel,
     +     Ru,Romg,Rphi,
     +     Kuu,Kuomg,Kuphi,Komgu,Komgomg,Komgphi,Kphiu,Kphiomg,Kphiphi,
     +     rhs,amatrx)
      
      !
      ! Subroutine to assemble the local elements residual and tangent
      !

      implicit none

      integer i,j,k,l,m,n,A11,A12,B11,B12,nDim,nNode,nDofEl,nDofN

      real*8 Ru(nDim*nNode,1),Romg(nNode,1),Kuu(nDim*nNode,nDim*nNode)
      real*8 Komgomg(nNode,nNode),Kuomg(nDim*nNode,nNode),rhs(ndofel,1)
      real*8 Komgu(nNode,nDim*nNode),amatrx(ndofel,ndofel),Rphi(nNode,1)
      real*8 Kphiphi(nNode,nNode),Komgphi(nNode,nNode)
      real*8 Kphiomg(nNode,nNode)
      real*8 Kuphi(nDim*nNode,nNode),Kphiu(nNode,nDim*nNode)


      ! Total number of degrees of freedom per node
      !
      nDofN = nDofEl/nNode


      ! init
      !
      rhs(:,1) = 0.d0
      amatrx = 0.d0

      if(nDim.eq.2) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            !
            ! displacement
            !
            rhs(A11,1) = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            !
            ! chemical potential
            !
            rhs(A11+2,1) = Rphi(i,1)/1.d9
            !
            ! temperature
            !
            rhs(A11+3,1) = Romg(i,1)/1.d3
         enddo
         !
         ! Assemble the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11) = Kuu(A12,B12)
               amatrx(A11,B11+1) = Kuu(A12,B12+1)
               amatrx(A11+1,B11) = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               !
               ! chemical potential
               !
               amatrx(A11+2,B11+2) = Kphiphi(i,j)/1.d9
               !
               ! temperature
               !
               amatrx(A11+3,B11+3) = Komgomg(i,j)/1.d3
               !
               ! displacement - chemical potential
               !
               amatrx(A11,B11+2) = Kuphi(A12,j)
               amatrx(A11+1,B11+2) = Kuphi(A12+1,j)
               !
               ! displacement - temperature
               !
               amatrx(A11,B11+3) = Kuomg(A12,j)
               amatrx(A11+1,B11+3) = Kuomg(A12+1,j)
               !
               ! chemical potential - displacement
               !
               amatrx(A11+2,B11) = Kphiu(i,B12)/1.d9
               amatrx(A11+2,B11+1) = Kphiu(i,B12+1)/1.d9
               !
               ! chemical potential - temperature
               !
               amatrx(A11+2,B11+3) = Kphiomg(i,j)/1.d9
               !
               ! temperature - displacement
               !
               amatrx(A11+3,B11) = Komgu(i,B12)/1.d3
               amatrx(A11+3,B11+1) = Komgu(i,B12+1)/1.d3
               !
               ! temperature - chemical potential
               !
               amatrx(A11+3,B11+2) = Komgphi(i,j)/1.d3
               !
            enddo
         enddo
         !
      elseif(nDim.eq.3) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            !
            ! displacement
            !
            rhs(A11,1)   = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            rhs(A11+2,1) = Ru(A12+2,1)
            !
            ! chemical potential
            !
            rhs(A11+3,1) = Rphi(i,1)/1.d9
            !
            ! temperature
            !
            rhs(A11+4,1) = Romg(i,1)/1.d3
            !
         enddo
         !
         ! Assembly the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11)     = Kuu(A12,B12)
               amatrx(A11,B11+1)   = Kuu(A12,B12+1)
               amatrx(A11,B11+2)   = Kuu(A12,B12+2)
               amatrx(A11+1,B11)   = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               amatrx(A11+1,B11+2) = Kuu(A12+1,B12+2)
               amatrx(A11+2,B11)   = Kuu(A12+2,B12)
               amatrx(A11+2,B11+1) = Kuu(A12+2,B12+1)
               amatrx(A11+2,B11+2) = Kuu(A12+2,B12+2)
               !
               ! chemical potential
               !
               amatrx(A11+3,B11+3) = Kphiphi(i,j)/1.d9
               !
               ! temperature
               !
               amatrx(A11+4,B11+4) = Komgomg(i,j)/1.d3
               !
               ! displacement - chemical potential
               !
               amatrx(A11,B11+3) = Kuphi(A12,j)
               amatrx(A11+1,B11+3) = Kuphi(A12+1,j)
               amatrx(A11+2,B11+3) = Kuphi(A12+2,j)
               !
               ! displacement - temperature
               !
               amatrx(A11,B11+4) = Kuomg(A12,j)
               amatrx(A11+1,B11+4) = Kuomg(A12+1,j)
               amatrx(A11+2,B11+4) = Kuomg(A12+2,j)
               !
               ! chemical potential - displacement
               !
               amatrx(A11+3,B11) = Kphiu(i,B12)/1.d9
               amatrx(A11+3,B11+1) = Kphiu(i,B12+1)/1.d9
               amatrx(A11+3,B11+2) = Kphiu(i,B12+2)/1.d9
               !
               ! chemical potential - temperature
               !
               amatrx(A11+3,B11+4) = Kphiomg(i,j)/1.d9
               !
               ! temperature - displacement
               !
               amatrx(A11+4,B11)   = Komgu(i,B12)/1.d3
               amatrx(A11+4,B11+1) = Komgu(i,B12+1)/1.d3
               amatrx(A11+4,B11+2) = Komgu(i,B12+2)/1.d3
               !
               ! temperature - chemical potential
               !
               amatrx(A11+4,B11+3) = Komgphi(i,j)/1.d3
               !
            enddo
         enddo
         !
      else
         write(*,*) 'How did you get nDim=',nDim
         call xit
      endif

      return
      end subroutine AssembleElement

!***********************************************************************
!     Dummy material for visualization purposes.
!     Returns zero stress and Jacobians.
!   

      subroutine UMAT(stress, statev, ddsdde, sse, spd, scd,
     &     rpl, ddsddt, drplde, drpldt,
     &     stran, dstran, time, dtime, temp, dtemp, predef, dpred,
     &     cmname, ndi, nshr, ntens, nstatv, props, nprops, coords,
     &     drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer,
     &     kspt, jstep, kinc)

      use global

c      implicit none

      character(len=80) cmname

!     Formal arguments to function

      real(8) :: stress(ntens), statev(nstatv), ddsdde(ntens, ntens)
      real(8) :: ddsddt(ntens), drplde(ntens), drpldt, stran(ntens)
      real(8) :: dstran(ntens), time(2), dtime, predef(1), dpred(1)
      real(8) :: props(nprops), coords(3), drot(3,3), dfgrd0(3,3)
      real(8) :: dfgrd1(3,3), sse, spd, scd, rpl, temp, dtemp, pnewdt
      real(8) :: celent

      integer :: ndi, nshr, ntens, nstatv, nprops, noel, npt, layer
      integer :: kspt, jstep(4), kinc

      integer i

      ! Parameters
      !
      real*8 zero,one,two,half,three,third,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0)


      ! Do nothing if a "dummy" step
      !
      if(dtime.eq.zero) return
      
	  
	  ! Transfer global SDV's to UMAT's statev() array for visualization.
	  do i=1,nstatv
	     statev(i) = globalSdv(noel-ElemOffset,npt,i)
	  end do

      ! No mechanics contributions
		stress = zero
		ddsdde = zero

      return
      end subroutine umat

!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint2D1pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(1,2), w(1)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w = 4.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0


      return
      end subroutine xint2D1pt
      
!************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2), w(4)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint2D4pt

************************************************************************

      subroutine xint3D1pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 2 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(1,3),w(1)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w(1) = 8.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0
      xi(1,3) = 0.d0

      return
      end subroutine xint3D1pt
     
************************************************************************

      subroutine xint3D8pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(8,3),w(8)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 8


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      w(5) = 1.d0
      w(6) = 1.d0
      w(7) = 1.d0
      w(8) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = dsqrt(1.d0/3.d0)
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = dsqrt(1.d0/3.d0)
      xi(7,1) = -dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = dsqrt(1.d0/3.d0)
      xi(8,1) = dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D8pt

************************************************************************

      subroutine xintSurf2D1pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(1),yLocal(1),w(1),zero,one,two
      parameter(zero=0.d0,one=1.d0,two=2.d0)


      ! Gauss weights
      !
      w(1) = two
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = zero
         yLocal(1) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = zero
      elseif(face.eq.3) then
         xLocal(1) = zero
         yLocal(1) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = zero
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D1pt

************************************************************************

      subroutine xintSurf2D2pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(2),yLocal(2),w(2),one,three
      parameter(one=1.d0,three=3.d0)


      ! Gauss weights
      !
      w(1) = one
      w(2) = one
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(one/three)
         xLocal(2) = one
         yLocal(2) = dsqrt(one/three)
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = dsqrt(one/three)
         xLocal(2) = -one
         yLocal(2) = -dsqrt(one/three)
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D2pt

************************************************************************

      subroutine xintSurf2D3pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(3),yLocal(3),w(3),zero,one,two,three,five,eight,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,five=5.d0,
     +     eight=8.d0,nine=9.d0)


      ! Gauss weights
      !
      w(1) = five/nine
      w(2) = eight/nine
      w(3) = five/nine
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(three/five)
         yLocal(1) = -one
         xLocal(2) = zero
         yLocal(2) = -one
         xLocal(2) = dsqrt(three/five)
         yLocal(2) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(three/five)
         xLocal(2) = one
         yLocal(2) = zero
         xLocal(3) = one
         yLocal(3) = dsqrt(three/five)
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(three/five)
         yLocal(1) = one
         xLocal(2) = zero
         yLocal(2) = one
         xLocal(3) = dsqrt(three/five)
         yLocal(3) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = dsqrt(three/five)
         xLocal(2) = -one
         yLocal(2) = zero
         xLocal(3) = -one
         yLocal(3) = -dsqrt(three/five)
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D3pt

************************************************************************

      subroutine xintSurf3D1pt(face,xLocal,yLocal,zLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 1 gauss point for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  zLocal(nIntPt): z coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(1),yLocal(1),zLocal(1),w(1),zero,one,four
      parameter(zero=0.d0,one=1.d0,four=4.d0)


      ! Gauss weights
      !
      w(1) = four
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = zero
         yLocal(1) = zero
         zLocal(1) = -one
      elseif(face.eq.2) then
         xLocal(1) = zero
         yLocal(1) = zero
         zLocal(1) = one
      elseif(face.eq.3) then
         xLocal(1) = zero
         yLocal(1) = -one
         zLocal(1) = zero
      elseif(face.eq.4) then
         xLocal(1) = one
         yLocal(1) = zero
         zLocal(1) = zero
      elseif(face.eq.5) then
         xLocal(1) = zero
         yLocal(1) = one
         zLocal(1) = zero
      elseif(face.eq.6) then
         xLocal(1) = -one
         yLocal(1) = zero
         zLocal(1) = zero
      else
         write(*,*) 'face.ne.1,2,3,4,5,6'
         write(80,*) 'face.ne.1,2,3,4,5,6'
         call xit
      endif

      end subroutine xintSurf3D1pt

************************************************************************

      subroutine xintSurf3D4pt(face,xLocal,yLocal,zLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 4 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  yLocal(nIntPt): z coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(4),yLocal(4),zLocal(4),w(4),one,three
      parameter(one=1.d0,three=3.d0)


      ! Gauss weights
      !
      w(1) = one
      w(2) = one
      w(3) = one
      w(4) = one
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -dsqrt(one/three)
         zLocal(2) = -one
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = -one
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = dsqrt(one/three)
         zLocal(4) = -one
      elseif(face.eq.2) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -dsqrt(one/three)
         zLocal(2) = one
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = one
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = dsqrt(one/three)
         zLocal(4) = one
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -one
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -one
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = -one
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = -one
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.4) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = one
         yLocal(2) = dsqrt(one/three)
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = one
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = one
         yLocal(4) = -dsqrt(one/three)
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.5) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = one
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = one
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = one
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = one
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.6) then
         xLocal(1) = -one
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = -one
         yLocal(2) = dsqrt(one/three)
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = -one
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -one
         yLocal(4) = -dsqrt(one/three)
         zLocal(4) = dsqrt(one/three)
      else
         write(*,*) 'face.ne.1,2,3,4,5,6'
         write(80,*) 'face.ne.1,2,3,4,5,6'
         call xit
      endif

      end subroutine xintSurf3D4pt
     
!************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element


      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   1-----------2        origin at center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      
      
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      

      return
      end subroutine calcShape2DLinear

************************************************************************

      subroutine calcShape3DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 8-node linear 3D element as shown
      !
      !      8-----------7
      !     /|          /|       zeta
      !    / |         / |       
      !   5-----------6  |       |     eta
      !   |  |        |  |       |   /
      !   |  |        |  |       |  /
      !   |  4--------|--3       | /
      !   | /         | /        |/
      !   |/          |/         O--------- xi
      !   1-----------2        origin at cube center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt,i,j

      real*8 xi_int(nIntPt,3),sh(8),dshxi(8,3)
      real*8 d2shxi(8,3,3),xi,eta,zeta

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      !
      ! The shape functions
      !
      sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
      sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
      sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
      sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
      sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
      sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
      sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
      sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
      !
      ! The first derivatives
      !
      dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
      dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
      dshxi(1,3) = -eighth*(one - xi)*(one - eta)
      dshxi(2,1) = eighth*(one - eta)*(one - zeta)
      dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
      dshxi(2,3) = -eighth*(one + xi)*(one - eta)
      dshxi(3,1) = eighth*(one + eta)*(one - zeta)
      dshxi(3,2) = eighth*(one + xi)*(one - zeta)
      dshxi(3,3) = -eighth*(one + xi)*(one + eta)
      dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
      dshxi(4,2) = eighth*(one - xi)*(one - zeta)
      dshxi(4,3) = -eighth*(one - xi)*(one + eta)
      dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
      dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
      dshxi(5,3) = eighth*(one - xi)*(one - eta)
      dshxi(6,1) = eighth*(one - eta)*(one + zeta)
      dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
      dshxi(6,3) = eighth*(one + xi)*(one - eta)
      dshxi(7,1) = eighth*(one + eta)*(one + zeta)
      dshxi(7,2) = eighth*(one + xi)*(one + zeta)
      dshxi(7,3) = eighth*(one + xi)*(one + eta)
      dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
      dshxi(8,2) = eighth*(one - xi)*(one + zeta)
      dshxi(8,3) = eighth*(one - xi)*(one + eta)
      !
      ! The second derivatives
      !
      d2shxi = zero
      d2shxi(1,1,2) = eighth*(one - zeta)
      d2shxi(1,2,1) = d2shxi(1,1,2)
      d2shxi(1,1,3) = eighth*(one - eta)
      d2shxi(1,3,1) = d2shxi(1,1,3)
      d2shxi(1,2,3) = eighth*(one - xi)
      d2shxi(1,3,2) = d2shxi(1,2,3)
      d2shxi(2,1,2) = -eighth*(one - zeta)
      d2shxi(2,2,1) = d2shxi(2,1,2)
      d2shxi(2,1,3) = -eighth*(one - eta)
      d2shxi(2,3,1) = d2shxi(2,1,3)
      d2shxi(2,2,3) = eighth*(one + xi)
      d2shxi(2,3,2) = d2shxi(2,2,3)
      d2shxi(3,1,2) = eighth*(one - zeta)
      d2shxi(3,2,1) = d2shxi(2,1,2)
      d2shxi(3,1,3) = -eighth*(one + eta)
      d2shxi(3,3,1) = d2shxi(2,1,3)
      d2shxi(3,2,3) = -eighth*(one + xi)
      d2shxi(3,3,2) = d2shxi(2,2,3)
      d2shxi(4,1,2) = -eighth*(one - zeta)
      d2shxi(4,2,1) = d2shxi(2,1,2)
      d2shxi(4,1,3) = eighth*(one + eta)
      d2shxi(4,3,1) = d2shxi(2,1,3)
      d2shxi(4,2,3) = -eighth*(one - xi)
      d2shxi(4,3,2) = d2shxi(2,2,3)
      d2shxi(5,1,2) = eighth*(one + zeta)
      d2shxi(5,2,1) = d2shxi(2,1,2)
      d2shxi(5,1,3) = -eighth*(one - eta)
      d2shxi(5,3,1) = d2shxi(2,1,3)
      d2shxi(5,2,3) = -eighth*(one - xi)
      d2shxi(5,3,2) = d2shxi(2,2,3)
      d2shxi(6,1,2) = eighth*(one + zeta)
      d2shxi(6,2,1) = d2shxi(2,1,2)
      d2shxi(6,1,3) = eighth*(one - eta)
      d2shxi(6,3,1) = d2shxi(2,1,3)
      d2shxi(6,2,3) = -eighth*(one + xi)
      d2shxi(6,3,2) = d2shxi(2,2,3)
      d2shxi(7,1,2) = eighth*(one + zeta)
      d2shxi(7,2,1) = d2shxi(2,1,2)
      d2shxi(7,1,3) = eighth*(one + eta)
      d2shxi(7,3,1) = d2shxi(2,1,3)
      d2shxi(7,2,3) = eighth*(one + xi)
      d2shxi(7,3,2) = d2shxi(2,2,3)
      d2shxi(8,1,2) = -eighth*(one + zeta)
      d2shxi(8,2,1) = d2shxi(2,1,2)
      d2shxi(8,1,3) = -eighth*(one + eta)
      d2shxi(8,3,1) = d2shxi(2,1,3)
      d2shxi(8,2,3) = eighth*(one - xi)
      d2shxi(8,3,2) = d2shxi(2,2,3)
      
      return
      end subroutine calcShape3DLinear

!************************************************************************


      subroutine computeSurf(xLocal,yLocal,face,coords,sh,ds)

      implicit none

      integer face

      real*8 xLocal,yLocal,ds,dshxi(4,2),sh(4),dXdXi,dXdEta,dYdXi
      real*8 dYdEta,one,coords(2,4),fourth,shape
      parameter(one=1.d0,fourth=1.d0/4.d0)

      sh(1) = fourth*(one - xLocal)*(one - yLocal)
      sh(2) = fourth*(one + xLocal)*(one - yLocal)
      sh(3) = fourth*(one + xLocal)*(one + yLocal)
      sh(4) = fourth*(one - xLocal)*(one + yLocal)
      
      dshxi(1,1) = -fourth*(one - yLocal)
      dshxi(1,2) = -fourth*(one - xLocal)
      dshxi(2,1) = fourth*(one - yLocal)
      dshxi(2,2) = -fourth*(one + xLocal)
      dshxi(3,1) = fourth*(one + yLocal)
      dshxi(3,2) = fourth*(one + xLocal)
      dshxi(4,1) = -fourth*(one + yLocal)
      dshxi(4,2) = fourth*(one - xLocal)

      dXdXi = dshxi(1,1)*coords(1,1)+dshxi(2,1)*coords(1,2)
     +     + dshxi(3,1)*coords(1,3)+dshxi(4,1)*coords(1,4)
      dXdEta = dshxi(1,2)*coords(1,1)+dshxi(2,2)*coords(1,2)
     +     + dshxi(3,2)*coords(1,3)+dshxi(4,2)*coords(1,4)
      dYdXi = dshxi(1,1)*coords(2,1)+dshxi(2,1)*coords(2,2)
     +     + dshxi(3,1)*coords(2,3)+dshxi(4,1)*coords(2,4)
      dYdEta = dshxi(1,2)*coords(2,1)+dshxi(2,2)*coords(2,2)
     +     + dshxi(3,2)*coords(2,3)+dshxi(4,2)*coords(2,4)

      if((face.eq.1).or.(face.eq.3).or.(face.eq.11).or.(face.eq.13))then
         ds = dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
      elseif((face.eq.2).or.(face.eq.4).or.
     +                                 (face.eq.12).or.(face.eq.14))then
         ds = dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
      else
         write(*,*) 'YOU HAVE A PROBLEM WITH YOUR FACE'
      endif

      end subroutine computeSurf

************************************************************************

      subroutine computeSurf3D(xLocal,yLocal,zLocal,face,coords,sh,dA)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the area dA, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 8-node brick elements

      implicit none

      integer face,stat,i,j,k

      real*8 xLocal,yLocal,zLocal,dA,dshxi(8,3),sh(8),zero,dsh(8,3),one
      real*8 coords(3,8),two,eighth,mapJ(3,3),mag,normal(3,1)

      real*8 dXdXi,dXdEta,dXdZeta,dYdXi,dYdEta,dYdZeta,dZdXi,dZdEta
      real*8 dZdZeta

      parameter(one=1.d0,two=2.d0,eighth=1.d0/8.d0,zero=0.d0)

      ! The shape functions
      !
      sh(1) = eighth*(one - xLocal)*(one - yLocal)*(one - zLocal)
      sh(2) = eighth*(one + xLocal)*(one - yLocal)*(one - zLocal)
      sh(3) = eighth*(one + xLocal)*(one + yLocal)*(one - zLocal)
      sh(4) = eighth*(one - xLocal)*(one + yLocal)*(one - zLocal)
      sh(5) = eighth*(one - xLocal)*(one - yLocal)*(one + zLocal)
      sh(6) = eighth*(one + xLocal)*(one - yLocal)*(one + zLocal)
      sh(7) = eighth*(one + xLocal)*(one + yLocal)*(one + zLocal)
      sh(8) = eighth*(one - xLocal)*(one + yLocal)*(one + zLocal)


      ! Shape function derivatives
      !
      dshxi(1,1) = -eighth*(one - yLocal)*(one - zLocal)
      dshxi(1,2) = -eighth*(one - xLocal)*(one - zLocal)
      dshxi(1,3) = -eighth*(one - xLocal)*(one - yLocal)
      dshxi(2,1) = eighth*(one - yLocal)*(one - zLocal)
      dshxi(2,2) = -eighth*(one + xLocal)*(one - zLocal)
      dshxi(2,3) = -eighth*(one + xLocal)*(one - yLocal)
      dshxi(3,1) = eighth*(one + yLocal)*(one - zLocal)
      dshxi(3,2) = eighth*(one + xLocal)*(one - zLocal)
      dshxi(3,3) = -eighth*(one + xLocal)*(one + yLocal)
      dshxi(4,1) = -eighth*(one + yLocal)*(one - zLocal)
      dshxi(4,2) = eighth*(one - xLocal)*(one - zLocal)
      dshxi(4,3) = -eighth*(one - xLocal)*(one + yLocal)
      dshxi(5,1) = -eighth*(one - yLocal)*(one + zLocal)
      dshxi(5,2) = -eighth*(one - xLocal)*(one + zLocal)
      dshxi(5,3) = eighth*(one - xLocal)*(one - yLocal)
      dshxi(6,1) = eighth*(one - yLocal)*(one + zLocal)
      dshxi(6,2) = -eighth*(one + xLocal)*(one + zLocal)
      dshxi(6,3) = eighth*(one + xLocal)*(one - yLocal)
      dshxi(7,1) = eighth*(one + yLocal)*(one + zLocal)
      dshxi(7,2) = eighth*(one + xLocal)*(one + zLocal)
      dshxi(7,3) = eighth*(one + xLocal)*(one + yLocal)
      dshxi(8,1) = -eighth*(one + yLocal)*(one + zLocal)
      dshxi(8,2) = eighth*(one - xLocal)*(one + zLocal)
      dshxi(8,3) = eighth*(one - xLocal)*(one + yLocal)


      dXdXi = zero
      dXdEta = zero
      dXdZeta = zero
      dYdXi = zero
      dYdEta = zero
      dYdZeta = zero
      dZdXi = zero
      dZdEta = zero
      dZdZeta = zero
      do k=1,8
         dXdXi = dXdXi + dshxi(k,1)*coords(1,k)
         dXdEta = dXdEta + dshxi(k,2)*coords(1,k)
         dXdZeta = dXdZeta + dshxi(k,3)*coords(1,k)
         dYdXi = dYdXi + dshxi(k,1)*coords(2,k)
         dYdEta = dYdEta + dshxi(k,2)*coords(2,k)
         dYdZeta = dYdZeta + dshxi(k,3)*coords(2,k)
         dZdXi = dZdXi + dshxi(k,1)*coords(3,k)
         dZdEta = dZdEta + dshxi(k,2)*coords(3,k)
         dZdZeta = dZdZeta + dshxi(k,3)*coords(3,k)
      enddo


      ! Jacobian of the mapping
      !
      if((face.eq.1).or.(face.eq.2)) then
         ! zeta = constant on this face
         dA = dsqrt(
     +          (dYdXi*dZdEta - dYdEta*dZdXi)**two
     +        + (dXdXi*dZdEta - dXdEta*dZdXi)**two
     +        + (dXdXi*dYdEta - dXdEta*dYdXi)**two
     +        )
      elseif((face.eq.3).or.(face.eq.5)) then
         ! eta = constant on this face
         dA = dsqrt(
     +          (dYdXi*dZdZeta - dYdZeta*dZdXi)**two
     +        + (dXdXi*dZdZeta - dXdZeta*dZdXi)**two
     +        + (dXdXi*dYdZeta - dXdZeta*dYdXi)**two
     +        )
      elseif((face.eq.4).or.(face.eq.6)) then
         ! xi = constant on this face
         dA = dsqrt(
     +          (dYdEta*dZdZeta - dYdZeta*dZdEta)**two
     +        + (dXdEta*dZdZeta - dXdZeta*dZdEta)**two
     +        + (dXdEta*dYdZeta - dXdZeta*dYdEta)**two
     +        )
         else
            write(*,*) 'never should get here'
            call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      !
      if((face.eq.1).or.(face.eq.2)) then
         ! zeta = constant on this face
         normal(1,1) = dYdXi*dZdEta - dYdEta*dZdXi
         normal(2,1) = dXdXi*dZdEta - dXdEta*dZdXi
         normal(3,1) = dXdXi*dYdEta - dXdEta*dYdXi
         if(face.eq.1) normal = -normal
      elseif((face.eq.3).or.(face.eq.5)) then
         ! eta = constant on this face
         normal(1,1) = dYdXi*dZdZeta - dYdZeta*dZdXi
         normal(2,1) = dXdXi*dZdZeta - dXdZeta*dZdXi
         normal(3,1) = dXdXi*dYdZeta - dXdZeta*dYdXi
         if(face.eq.5) normal = -normal
      elseif((face.eq.4).or.(face.eq.6)) then
         ! xi = constant on this face
         normal(1,1) = dYdEta*dZdZeta - dYdZeta*dZdEta
         normal(2,1) = dXdEta*dZdZeta - dXdZeta*dZdEta
         normal(3,1) = dXdEta*dYdZeta - dXdZeta*dYdEta
         if(face.eq.6) normal = -normal
      else
         write(*,*) 'never should get here'
         call xit
      endif
      mag = dsqrt(normal(1,1)**two+normal(2,1)**two+normal(3,1)**two)
      normal(1,1) = normal(1,1)/mag
      normal(2,1) = normal(2,1)/mag
      normal(3,1) = normal(3,1)/mag

      end subroutine computeSurf3D

************************************************************************

      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape2D'
         call xit
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2D

!*************************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a "heat transfer" and 
      !  "static" step uses MCRD=2, but for "coupled-temperature-displacement"
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape2Da'
         call xit
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2Da

************************************************************************

      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat,pnewdt)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.  This subroutine works for both 8-node
      !  linear and 20-node quadratic 3D elements.
      !
      implicit none

      integer i,j,k,nNode,ieror,stat

      real*8 dshxi(nNode,3),dsh(nNode,3),coords(3,nNode)
      real*8 mapJ(3,3),mapJ_inv(3,3),detmapJ, pnewdt

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,3
        do j=1,3
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv3D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape3D. Cutting back 66.6%.'
         !call xit
		 pnewdt = 0.333
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      ! The second derivatives may be calculated.
      !

      return
      end subroutine mapShape3D

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet
	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

****************************************************************************
      SUBROUTINE SKINEM(F,R,U,E)

C	THIS SUBROUTINE PERFORMS THE RIGHT POLAR DECOMPOSITION
C	[F] = [R][U] OF THE DEFORMATION GRADIENT [F] INTO
C	A ROTATION [R] AND THE RIGHT  STRETCH TENSOR [U].
C	THE EIGENVALUES AND EIGENVECTORS OF [U] AND
C	THE LOGARITHMIC STRAIN [E] = LN [U]
C	ARE ALSO RETURNED.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION F(3,3),FTRANS(3,3), C(3,3), OMEGA(3),
     +	          UEIGVAL(3),EIGVEC(3,3), EIGVECT(3,3), 
     +            U(3,3),E(3,3),UINV(3,3),R(3,3),TEMPM(3,3)
     
	COMMON/ERRORINFO/UMERROR     

C	F(3,3)	-- THE DEFORMATION GRADIENT MATRIX WHOSE
C		   POLAR DECOMPOSITION IS DESIRED.
C	DETF	-- THE DETRMINANT OF [F]; DETF > 0.
C	FTRANS(3,3)	-- THE TRANSPOSE OF [F].
C	R(3,3)	-- THE ROTATION MATRIX; [R]^T [R] = [I];
C		   OUTPUT.
C	U(3,3)	-- THE RIGHT STRETCH TENSOR; SYMMETRIC
C		   AND POSITIVE DEFINITE; OUTPUT.
C	UINV(3,3)	-- THE INVERSE OF [U].
C	C(3,3)	-- THE RIGHT CAUCHY-GREEN TENSOR = [U][U];
C		   SYMMETRIC AND POSITIVE DEFINITE.
C	OMEGA(3)-- THE SQUARES OF THE PRINCIPAL STRETCHES.
C 	UEIGVAL(3)	-- THE PRINCIPAL STRETCHES; OUTPUT.
C	EIGVEC(3,3)	-- MATRIX OF EIGENVECTORS OF [U];OUTPUT.
C	EIGVECT(3,3)    -- TRANSPOSE OF THE ABOVE.
C	E(3,3)	-- THE LOGARITHMIC STRAIN TENSOR, [E]=LN[U];
C		   OUTPUT.
C**********************************************************************

C	STORE THE IDENTITY MATRIX IN  [R], [U], AND [UINV]

	CALL ONEM(R)
	CALL ONEM(U)
	CALL ONEM(UINV)

C	STORE THE ZERO MATRIX IN [E]

	CALL ZEROM(E)

C      	CHECK IF THE DETERMINANT OF [F] IS GREATER THAN ZERO.
C	IF NOT, THEN PRINT DIAGNOSTIC AND STOP.

        CALL MDET(F,DETF)
        IF (DETF .LE. 0.D0) THEN
          WRITE(80,100)
          UMERROR=5.
          RETURN
        ENDIF

C      	CALCULATE THE RIGHT CAUCHY GREEN STRAIN TENSOR [C]

        CALL  MTRANS(F,FTRANS)
        CALL  MPROD(FTRANS,F,C)
 
C	CALCULATE THE EIGENVALUES AND EIGENVECTORS OF  [C]

	CALL SPECTRAL(C,OMEGA,EIGVEC)

C	CALCULATE THE PRINCIPAL VALUES OF [U] AND [E]

	UEIGVAL(1) = DSQRT(OMEGA(1))
	UEIGVAL(2) = DSQRT(OMEGA(2))
	UEIGVAL(3) = DSQRT(OMEGA(3))

	U(1,1) = UEIGVAL(1)
	U(2,2) = UEIGVAL(2)
	U(3,3) = UEIGVAL(3)

	E(1,1) = DLOG( UEIGVAL(1) )
	E(2,2) = DLOG( UEIGVAL(2) )
	E(3,3) = DLOG( UEIGVAL(3) )

C	CALCULATE THE COMPLETE TENSORS [U] AND [E]

	CALL MTRANS(EIGVEC,EIGVECT)
	CALL MPROD(EIGVEC,U,TEMPM)
	CALL MPROD(TEMPM,EIGVECT,U)
	CALL MPROD(EIGVEC,E,TEMPM)
	CALL MPROD(TEMPM,EIGVECT,E)

C	CALCULATE [UINV]

	CALL M3INV(U,UINV)

C	CALCULATE [R]

	CALL MPROD(F,UINV,R)
100     FORMAT(5X,'--ERROR IN KINEMATICS-- THE DETERMINANT OF [F]',
     +         ' IS NOT GREATER THAN 0')

	RETURN
	END

C**********************************************************************
C	THE FOLLOWING SUBROUTINES CALCULATE THE SPECTRAL
C	DECOMPOSITION OF A SYMMETRIC THREE BY THREE MATRIX
C**********************************************************************
	SUBROUTINE SPECTRAL(A,D,V)
C
C	THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
C	A SYMMETRIC 3 BY 3 MATRIX [A]. 
C
C	THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
C	EIGENVALUES IN DESCENDING ORDER, AND
C	A MATRIX [V] WHOSE COLUMNS CONTAIN THE CORRESPONDING
C	EIGENVECTORS.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER(NP=3)
	DIMENSION D(NP),V(NP,NP)
	DIMENSION A(3,3),E(NP,NP)

	DO 2 I = 1,3
          DO 1 J= 1,3
            E(I,J) = A(I,J)
1	  CONTINUE
2	CONTINUE

	CALL JACOBI(E,3,NP,D,V,NROT)
	CALL EIGSRT(D,V,3,NP)

	RETURN
	END

C**********************************************************************
	SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

C	COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C	MATRIX [A], WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL 
C	NP BY BP ARRAY. ON OUTPUT, ELEMENTS OF [A] ABOVE THE DIAGONAL 
C	ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
C	AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
C	VECTOR D RETURNS THE EIGENVALUES OF [A] IN ITS FIRST N ELEMENTS.
C	[V] IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
C	[A] WHOSE COLUMNS CONTAIN, ON OUTPUT, THE NORMALIZED
C	EIGENVECTORSOF [A]. NROT RETURNS THE NUMBER OF JACOBI ROTATIONS
C	WHICH WERE REQUIRED.

C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", PAGE 346.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (NMAX =100)
	DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

C	INITIALIZE [V] TO THE IDENTITY MATRIX

	DO 12 IP = 1,N	
	  DO 11 IQ = 1,N
	    V(IP,IQ) = 0.D0
11        CONTINUE
          V(IP,IP) = 1.D0
12	CONTINUE

C	INITIALIZE [B] AND [D] TO THE DIAGONAL OF [A], AND Z TO ZERO.
C	THE VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
C	IN EQUATION (11.1.14)

	DO 13 IP = 1,N
	  B(IP) = A(IP,IP)
	  D(IP) = B(IP)
	  Z(IP) = 0.D0
13	CONTINUE
C
	NROT = 0
	DO 24 I = 1,50

C	SUM OFF-DIAGONAL ELEMENTS

          SM = 0.D0
          DO 15 IP = 1, N-1
            DO 14 IQ = IP + 1, N
	      SM = SM + DABS ( A(IP,IQ ))
14          CONTINUE
15        CONTINUE

C	IF SUM = 0., THEN RETURN. THIS IS THE NORMAL RETURN
C	WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE 
C	UNDERFLOW.

          IF ( SM .EQ. 0.D0) RETURN
C
C	  IF ( SM .LT. 1.0D-15) RETURN

C	IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
C	|A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE, 
C	SEE EQUATION (11.1.25). THEREAFTER TRESH = 0.

          IF ( I .LT. 4) THEN
            TRESH = 0.2D0*SM/N**2
          ELSE
            TRESH = 0.D0
          ENDIF
C
          DO 22 IP = 1, N-1
            DO 21 IQ = IP+1,N
              G = 100.D0*DABS(A(IP,IQ))

C	AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
C	OFF-DIAGONAL ELEMENT IS SMALL.

	      IF ((I .GT. 4) .AND. (DABS(D(IP))+G .EQ. DABS(D(IP)))
     +            .AND. ( DABS(D(IQ))+G .EQ. DABS(D(IQ)))) THEN
                A(IP,IQ) = 0.D0
              ELSE IF ( DABS(A(IP,IQ)) .GT. TRESH) THEN
                H = D(IQ) - D(IP)
                IF (DABS(H)+G .EQ. DABS(H)) THEN

C	T = 1./(2.*THETA), EQUATION(11.1.10)

	          T =A(IP,IQ)/H
	        ELSE
	          THETA = 0.5D0*H/A(IP,IQ)
	          T =1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
	          IF (THETA .LT. 0.D0) T = -T
	        ENDIF
	        C = 1.D0/DSQRT(1.D0 + T**2)
	        S = T*C
	        TAU = S/(1.D0 + C)
	        H = T*A(IP,IQ)
	        Z(IP) = Z(IP) - H
	        Z(IQ) = Z(IQ) + H
	        D(IP) = D(IP) - H
	        D(IQ) = D(IQ) + H
	        A(IP,IQ) = 0.D0

C	CASE OF ROTATIONS 1 <= J < P
				
	        DO 16 J = 1, IP-1
	          G = A(J,IP)
	          H = A(J,IQ)
	          A(J,IP) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
16	        CONTINUE

C	CASE OF ROTATIONS P < J < Q

	        DO 17 J = IP+1, IQ-1
	          G = A(IP,J)
	          H = A(J,IQ)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
17	        CONTINUE

C	CASE OF ROTATIONS Q < J <= N

	        DO 18 J = IQ+1, N
                  G = A(IP,J)
	          H = A(IQ,J)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(IQ,J) = H + S*(G - H*TAU)
18	        CONTINUE
	        DO 19 J = 1,N
	          G = V(J,IP)
	          H = V(J,IQ)
	          V(J,IP) = G - S*(H + G*TAU)
	          V(J,IQ) = H + S*(G - H*TAU)
19	        CONTINUE
	        NROT = NROT + 1
              ENDIF
21	    CONTINUE
22	  CONTINUE

C	UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z

	  DO 23 IP = 1, N
	    B(IP) = B(IP) + Z(IP)
	    D(IP) = B(IP)
	    Z(IP) = 0.D0
23	  CONTINUE
24	CONTINUE

C	IF THE ALGORITHM HAS REACHED THIS STAGE, THEN
C	THERE ARE TOO MANY SWEEPS, PRINT A DIAGNOSTIC
C	AND STOP.

	WRITE (80,'(/1X,A/)') '50 ITERS IN JACOBI SHOULD NEVER HAPPEN'

	RETURN
	END

C**********************************************************************
	SUBROUTINE EIGSRT(D,V,N,NP)

C	GIVEN THE EIGENVALUES [D] AND EIGENVECTORS [V] AS OUTPUT FROM
C	JACOBI, THIS ROUTINE SORTS THE EIGENVALUES INTO DESCENDING ORDER, 
C	AND REARRANGES THE COLUMNS OF [V] ACCORDINGLY.

C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", P. 348.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION D(NP),V(NP,NP)

	DO 13 I = 1,N-1
	  K = I
	  P = D(I)
	  DO 11 J = I+1,N
	    IF (D(J) .GE. P) THEN
	      K = J
	      P = D(J)
	    END IF
11	  CONTINUE
	  IF (K .NE. I) THEN
	    D(K) = D(I)
	    D(I) = P
	    DO 12 J = 1,N
	      P = V(J,I)
	      V(J,I) = V(J,K)
	      V(J,K) = P
12	    CONTINUE
  	  ENDIF
13	CONTINUE

	RETURN
      END
C**********************************************************************
      	SUBROUTINE ZEROM(A)
C
C	THIS SUBROUTINE SETS ALL ENTRIES OF A 3 BY 3 MATRIX TO 0.D0.
C**********************************************************************

        REAL*8 A(3,3)

	DO 1 I=1,3
	  DO 1 J=1,3
	    A(I,J) = 0.D0
1	CONTINUE
C	
	RETURN
	END

C**********************************************************************
C**********************************************************************
	SUBROUTINE MTRANS(A,ATRANS)
 
C	THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
C	MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C**********************************************************************

	REAL*8 A(3,3),ATRANS(3,3)

	DO 1 I=1,3
 	  DO 1 J=1,3
	    ATRANS(J,I) = A(I,J)
1	CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE MPROD(A,B,C)
 
C 	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C 	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

	REAL*8 A(3,3),B(3,3),C(3,3)

	DO 2 I = 1, 3
	  DO 2 J = 1, 3
	    C(I,J) = 0.D0
	    DO 1 K = 1, 3
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1	    CONTINUE
2	CONTINUE
C
	RETURN
      END
C**********************************************************************
	SUBROUTINE M3INV(A,AINV)

C 	THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX
C	[A] AND PLACES THE RESULT IN [AINV]. 
C 	IF DET(A) IS ZERO, THE CALCULATION
C 	IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C**********************************************************************

	REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)

C	A(3,3)	        -- THE MATRIX WHOSE INVERSE IS DESIRED.
C	DET		-- THE COMPUTED DETERMINANT OF [A].
C	ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
C			   THE SIGNED MINOR (-1)**(I+J)*M_IJ
C			   IS CALLED THE COFACTOR OF A(I,J).
C	AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
C			   OBTAINED BY REPLACING EACH ELEMENT OF
C			   [A] BY ITS COFACTOR, AND THEN TAKING
C			   TRANSPOSE OF THE RESULTING MATRIX.
C	AINV(3,3)	-- RETURNED AS INVERSE OF [A].
C			   [AINV] = [AADJ]/DET.
C----------------------------------------------------------------------

	CALL MDET(A,DET)
	IF ( DET .EQ. 0.D0 ) THEN
	  WRITE(80,10)
	  write(*,10)
	  call xit
	ENDIF
	CALL MCOFAC(A,ACOFAC)
	CALL MTRANS(ACOFAC,AADJ)
	DO 1 I = 1,3
	DO 1 J = 1,3
	     AINV(I,J) = AADJ(I,J)/DET
1	CONTINUE
10	FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +         10X,'PROGRAM TERMINATED')

	RETURN
	END

C**********************************************************************
	SUBROUTINE MCOFAC(A,ACOFAC)
 
C 	THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C 	AND PLACES THE RESULT IN [ACOFAC]. 
C**********************************************************************

	REAL*8  A(3,3), ACOFAC(3,3)

	ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
	ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
	ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
	ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
	ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
	ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
	ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
	ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
	ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

	RETURN
	END

C**********************************************************************      