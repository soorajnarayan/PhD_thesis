************************************************************************
!
! User element for electro-chemical transport of ions 
!  and large elastic deformation of a host gel (polymer+solvent) in plane 
!  strain. 
!
! Solution variables (or nodal variables) are the displacements, electrostatic 
! potential, solvent chemical potential and 2 solutes electrochemical potential.
! 
! This subroutine is for the following element types
!  > two-dimensional 4 node isoparametric element as shown below
!       with 1pt (reduced) or 4pt (full) gauss integration.
!
! In order to avoid locking for the fully-integrated element, we
!  use the F-bar method of de Souza Neto (1996).
!
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
!
! Sooraj Narayan, Apr 2021.
!
***********************************************************************
!
! User element statement in the input file (set ? values as needed):
!
!  2D elements
!  *User Element,Nodes=4,Type=U1,Iproperties=2,Properties=21,Coordinates=2,Variables=12,Unsymm
!  1,2,4,11,12,13
!
!
!
!     State Variables
!     --------------------------------------------------------------
!     Global SDV's (used for visualization)
!       1) 
!       2)
!       |      
!      
!     Local SDV's (used for the solution procedure)
!       1) 
!       2)
!       | 
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
!     
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
      parameter(numElem=8000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=10000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable :: globalSdv(:,:,:)

      end module global

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
      parameter(nInt=4)  ! number of volume integration pionts
      parameter(nIntS=1) ! number of surface integration points
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
      !  deformation problem, I think this only matters
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
      ! Call the paricular element to perform the analysis
      !
      if(jtype.eq.1) then
         !
         ! This is an axisymmetric analysis
         !
         nDim = 2
c         call UAX4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
c     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
c     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
c     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
c     +        NJPROP,PERIOD,
c     +        nDim,nInt,nIntS)
         !
         !
      elseif(jtype.eq.2) then
         !
         ! This is a plane strain analysis
         !
         nDim = 2
         call UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,nInt,nIntS)
         !
         !
      elseif(jtype.eq.3) then
         !
         ! This is a 3D analysis
         !
         nDim = 3
c         call U3D8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
c     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
c     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
c     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
c     +        NJPROP,PERIOD,
c     +        nDim,nInt,nIntS)
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

      subroutine UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
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
      real*8 u(nNode,2),du(nNode,ndofel),uNew(nNode,ndofel)
      real*8 muNew(nNode), muOld(nNode),dMU(nNode)
      real*8 phiOld(nNode),dphi(nNode),phiNew(nNode)
      real*8 omgPOld(nNode),domgP(nNode),omgPNew(nNode)
      real*8 omgNOld(nNode),domgN(nNode),omgNNew(nNode)
      real*8 omgHOld(nNode),domgH(nNode),omgHNew(nNode)
      real*8 uOld(nNode,ndofel),u_t(nNode,ndofel),v(nNode,2)
      real*8 coordsC(mcrd,nNode)
      
      ! All int. pt. DOF related variables
      real*8 omgP_tau,omgP_t,domgPdx(2,1),domgPdt
      real*8 omgN_tau,omgN_t,domgNdx(2,1),domgNdt
      real*8 omgH_tau,omgH_t,domgHdx(2,1),domgHdt
      real*8 mu_tau,mu_t,dmudx(2,1),dmudt
      real*8 F_tau(3,3),F_t(3,3),detF_tau,detF,detF_t
      real*8 phi_tau,phi_t,dphidx(2,1),dphidt

      !Integral variables mainly for loop counters
      integer i,j,k,l,m,n,nn,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,lenJobName,lenOutDir,faceFlag,nIntS

      !Residuals and tangents
      real*8 Ru(2*nNode,1),Rmu(nNode,1),RomgP(nNode,1),RomgN(nNode,1)
      real*8 Rphi(nNode,1),RomgH(nNode,1)
      real*8 Kuu(2*nNode,2*nNode),Kmumu(nNode,nNode)
      real*8 Kphiphi(nNode,nNode),KomgPomgP(nNode,nNode)
      real*8 KomgNomgN(nNode,nNode),KomgHomgH(nNode,nNode)
      real*8 Kuphi(nDim*nNode,nNode),Kphiu(nNode,nDim*nNode)
      real*8 Kmuphi(nNode,nNode),Kphimu(nNode,nNode)
      real*8 Kmuu(nNode,2*nNode),Kumu(2*nNode,nNode)
      real*8 KomgPphi(nNode,nNode),KphiomgP(nNode,nNode)
      real*8 KomgPu(nNode,2*nNode),KuomgP(2*nNode,nNode)
      real*8 KomgNphi(nNode,nNode),KphiomgN(nNode,nNode)
      real*8 KomgNu(nNode,2*nNode),KuomgN(2*nNode,nNode)
      real*8 KmuomgP(nNode,nNode),KomgPmu(nNode,nNode)
      real*8 KmuomgN(nNode,nNode),KomgNmu(nNode,nNode)
      real*8 KmuomgH(nNode,nNode),KomgHmu(nNode,nNode)
      real*8 KomgNomgP(nNode,nNode),KomgPomgN(nNode,nNode)
      real*8 KomgNomgH(nNode,nNode),KomgPomgH(nNode,nNode)
      real*8 KomgHomgN(nNode,nNode),KomgHomgP(nNode,nNode)
      
      real*8 KphiomgH(nNode,nNode),KomgHphi(nNode,nNode)
      real*8 KomgHu(nNode,nDim*nNode),KuomgH(nDim*nNode,nNode)
      
      !Properties
      real*8 Omega,cbarP0,cbarN0,cRmaxP,cRmaxN,zP_mob,zN_mob
      real*8 cRf, zf
      real*8 Gshear, Kbulk
      real*8 D, DiffP, DiffN, K_rate, K1_rate, K2_rate
      
      real*8 DiffH,zH_mob,K3_rate,omgH0
      
      ! FEM quantities
      real*8 sh(nNode),detMapJ,dshC(nNode,2),w(nInt)
      real*8 sh0(nNode),detMapJ0,xi(nInt,2)
      real*8 dshxi(nNode,2),dsh0(nNode,2),dshC0(nNode,2),detMapJ0C
      real*8 wS(nIntS),xLocal(nIntS),yLocal(nIntS)
      
      real*8 flux_net
      
      !Misc. variables
      real*8 Iden(3,3),Le,theta0,zeta0,body(3)
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,ds,flux
      real*8 zeta_t,dsh(nNode,2),detMapJC,zetaLmt
      real*8 cbarP_t,cbarN_t,cbarP_tau,cbarN_tau
      real*8 cPdot,cNdot,dcPdomgP,dcNdomgN
      real*8 dcPdmu,dcNdmu,dmobPdmu,dmobNdmu
      real*8 dmobSdJ,umeror
      real*8 TR_tau(3,3),T_tau(3,3)
      real*8 SpTanMod(3,3,3,3),zeta_tau,dZdt,DzetaDmu,DzetadotDmu,mobS
      real*8 mobP,dmobPdomgP,mobN,dmobNdomgN
      real*8 Smat(3,1),Bmat(3,2*nNode),BodyForceRes(2*nNode,1),Qmat(4,4)
      real*8 dmobSdmu,Gmat(4,2*nNode),G0mat(4,2*nNode),Amat(4,4)
      real*8 Nvec(1,nNode),ResFac,AmatUC(3,1),TanFac
      real*8 dcPdphi,dcNdphi,dmobPdphi,dmobNdphi
      real*8 SpumuMod(3,3),Spmuumod(3,3,3),SpmuumodFac(3,3),AmatCU(2,4)
      ! All concentration related quantities in the UEL are in referential terms unless specified otherwise.
      real*8 cS_tau,cP_tau,cN_tau
      real*8 cS_t,cP_t,cN_t
      real*8 cP0,cN0
      !
      
      real*8 cHdot,dcHdomgH,dcHdmu,mobH,dmobHdmu,dmobHdomgH
      real*8 dcHdphi,dmobHdphi,cH_tau,cH_t,cH0
      real*8 cRf_dissoc,Ka_COOH,cH_eff
      real*8 rR_H
      real*8 dcRf_dissdcH_eff,dcH_effdomgH
      real*8 dcH_effdphi,dcH_effdzeta,dcH_effdmu
      
      
      real*8 SmatAx(4,1),BodyForceResAx(2*nNode,1)
      real*8 BmatAx(4,2*nNode),AmatAx(5,5),QmatAx(5,5)
      real*8 G0matAx(5,2*nNode),GmatAx(5,2*nNode)
      real*8 AR0,AR,ARc,AR_t

      real*8 zero,one,two,half,Pi,three,third,Rgas
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0,Rgas=8.314d0)


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
         write(*,*) '---------- UPE4 ELEMENTS ------------------------'
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
      Gshear = props(1)
      Kbulk = props(2)
      theta0 = props(8)
      zeta0   = props(9)
      cP0 = props(13)
      cN0 = props(17)
      
      cRmaxP= props(12)
      cRmaxN= props(16)
      
      zP_mob = props(14)
      zN_mob = props(18)
      
      cRf = props(19)
      zf  = props(20)
      
      D      = props(4)
      DiffP  = props(11)
      DiffN  = props(15)
      
      cH0 = props(21)
      zH_mob = props(22)
      DiffH = props(23)
      omgH0 = props(24)
      
      K_rate  =  5.d3*D/10.d-3  ! -3.d2   !5.d0!-3.d1 !*(one/zeta0-one)/Omega/Rgas/theta0
      K1_rate =  5.d3*DiffP*cP0/Rgas/theta0/10.d-3
      K2_rate =  5.d3*DiffN*cN0/Rgas/theta0/10.d-3
      K3_rate =  5.d3*DiffH*cH0/Rgas/theta0/10.d-3

      
      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rphi= zero
      Rmu = zero
      RomgP=zero
      RomgN=zero
      RomgH=zero
      
      Kuu = zero
      Kmumu=zero
      KomgPomgP=zero
      KomgNomgN=zero
      Kphiphi = zero
      KomgHomgH = zero
      
      Kuphi    = zero
      Kphiu    = zero
      Kumu     = zero
      Kmuu     = zero
      Kphimu   = zero
      Kmuphi   = zero
      
      KuomgP   = zero
      KomgPu   = zero
      KphiomgP = zero
      KomgPphi = zero
      KuomgN   = zero
      KomgNu   = zero
      KphiomgN = zero
      KomgNphi = zero
      
      KmuomgP  = zero
      KomgPmu  = zero
      KmuomgN  = zero
      KomgNmu  = zero
      KomgNomgP= zero
      KomgPomgN= zero
      KomgPomgH=zero
      KomgNomgH=zero
      KomgHomgN=zero
      KomgHomgP=zero
      
      KomgHu = zero
      KomgHphi = zero
      KphiomgH = zero
      KomgHmu = zero
      KmuomgH = zero
      
      Energy = zero


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
         omgHNew(i) = Uall(k)
         domgH(i) = DUall(k,1)
         omgHOld(i) = omgHNew(i) - domgH(i)
         k = k + 1
         muNew(i) = Uall(k)
         dMU(i) = DUall(k,1)
         muOld(i) = muNew(i) - dMU(i)
         k = k + 1
         omgPNew(i) = Uall(k)
         domgP(i) = DUall(k,1)
         omgPOld(i) = omgPNew(i) - domgP(i)
         k = k + 1
         omgNNew(i) = Uall(k)
         domgN(i) = DUall(k,1)
         omgNOld(i) = omgNNew(i) - domgN(i)
         enddo

         !omgPNew = zero
         !omgNNew = zero
         !phiNew  = zero
         !omgPOld = zero
         !omgNOld = zero
         !phiOld  = zero
         

      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


      ! Impose any time-stepping changes on the increments of
      !  chemical potential or displacement if you wish
      !
      ! chemical potential increment
      !
      do i=1,nNode
         if(dabs(dMU(i)).gt.1.d6) then
            pnewdt = 0.5
            write(*,*)'Problem in dmu'
            return
         endif
      enddo
      !
      !
      do i=1,nNode
         if(dabs(domgP(i)).gt.1.d6) then
            pnewdt = 0.5
            write(*,*)'Problem in domgP'
            return
         endif
      enddo
      !
      !
      do i=1,nNode
         if(dabs(domgN(i)).gt.1.d6) then
            pnewdt = 0.5
            write(*,*)'Problem in domgN'
            return
         endif
      enddo
      !
      !
      do i=1,nNode
         if(dabs(domgH(i)).gt.1.d6) then
            pnewdt = 0.5
            write(*,*)'Problem in domgH'
            return
         endif
      enddo
      !
      !
      do i=1,nNode
         if(dabs(dphi(i)).gt.1.d1) then
            pnewdt = 0.5
            write(*,*)'Problem in dphi'
            return
         endif
      enddo
      !
      ! displacement increment, based on element diagonal
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,3))**two) + 
     +     ((coordsC(2,1)-coordsC(2,3))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.10.0*Le) then
               pnewdt = 0.5
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
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat,pnewdt)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
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
       call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat,pnewdt)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
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

      ! For an axisymmetric problem, find the ``r'' that
      !  shows up in the integrals for axisymmetric
      !  and in the F(3,3), the factors of 2Pi are for integrals
      !  i.e., dV = 2 pi r dr dz
      !
      AR0  = zero
      ARc  = zero
      AR_t = zero
      do i=1,nNode
         ! radial coord in ref config at centroid
         AR0  = AR0 + sh0(i)*coords(1,i)
         ! radial coord in current config at centroid
         ARc  = ARc + sh0(i)*(coords(1,i) + u(i,1))
         ! radial coord in current config at centroid in previous step
         AR_t = AR_t + sh0(i)*(coords(1,i) + uOld(i,1))
         enddo
         

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
      !
      ! modify for PE
      !
      Fc_tau(3,3) = one
      Fc_t(3,3) = one
      !
      ! axisymmetric implementation detF
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
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
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
            zeta_t  = zeta0
            cP_t = cP0
            cN_t = cN0
            cH_t= cH0
            cRf_dissoc = cRf
            !
         else
            !
            ! this is not the first increment, read old values
            !
            zeta_t  = svars(1+jj)
            cP_t  = svars(2+jj)
            cN_t  = svars(3+jj)
            cH_t = svars(4+jj)
            cRf_dissoc = svars(5+jj)
            !
         endif


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            write(80,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
         

         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat,pnewdt)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
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
         call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat,pnewdt)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
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

       ! For an axisymmetric problem, find the ``r'' that
         !  shows up in the integrals for axisymmetric
         !  and in the F(3,3), the factors of 2Pi are for integrals
         !  i.e., dV = 2 pi r dr dz
         !
         !
         AR0  = zero
         AR   = zero
         AR_t = zero
         do i=1,nNode
            AR0 = AR0 + sh(i)*coords(1,i)
            AR  = AR  + sh(i)*(coords(1,i) + u(i,1))
            AR_t = AR_t + sh(i)*(coords(1,i) + uOld(i,1))
         enddo
         AR0  = two*Pi*AR0
         AR   = two*Pi*AR
         AR_t = two*Pi*AR_t
            
            
         ! Obtain the chemical potential and its derivative's at 
         !  this intPt at the begining and end of the incrment
         !
         mu_tau = zero
         mu_t = zero
         dMUdt = zero
         dMUdX = zero
         do k=1,nNode
            mu_tau = mu_tau + muNew(k)*sh(k)
            mu_t   = mu_t + muOld(k)*sh(k)
            do i=1,nDim
               dMUdX(i,1) = dMUdX(i,1) + muNew(k)*dshC(k,i)
            enddo
         enddo
         dMUdt = (mu_tau - mu_t)/dtime

         !
         omgP_tau = zero
         omgP_t = zero
         domgPdt = zero
         domgPdx = zero
         do k=1,nNode
            omgP_tau = omgP_tau + omgPNew(k)*sh(k)
            omgP_t   = omgP_t + omgPOld(k)*sh(k)
            do i=1,nDim
               domgPdx(i,1) = domgPdx(i,1) + omgPNew(k)*dshC(k,i)
            enddo
         enddo
         domgPdt = (omgP_tau - omgP_t)/dtime
         
         !
         omgN_tau = zero
         omgN_t = zero
         domgNdt = zero
         domgNdx = zero
         do k=1,nNode
            omgN_tau = omgN_tau + omgNNew(k)*sh(k)
            omgN_t   = omgN_t + omgNOld(k)*sh(k)
            do i=1,nDim
               domgNdx(i,1) = domgNdx(i,1) + omgNNew(k)*dshC(k,i)
            enddo
         enddo
         domgNdt = (omgN_tau - omgN_t)/dtime
         
         !
         omgH_tau = zero
         omgH_t = zero
         domgHdt = zero
         domgHdx = zero
         do k=1,nNode
            omgH_tau = omgH_tau + omgHNew(k)*sh(k)
            omgH_t   = omgH_t + omgHOld(k)*sh(k)
            do i=1,nDim
               domgHdx(i,1) = domgHdx(i,1) + omgHNew(k)*dshC(k,i)
            enddo
         enddo
         domgHdt = (omgH_tau - omgH_t)/dtime
         omgH_tau = omgH_tau-omgH0
         omgH_t   = omgH_t-omgH0
         
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
               dphidx(i,1) = dphidx(i,1) + phiNew(k)*dshC(k,i)
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
         ! modify F(3,3) for plane-strain 
         !
         F_tau(3,3) = one
         F_t(3,3) = one
         !
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 4 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.4).and.(nInt.eq.4)) then
            call mdet(F_tau,detF_tau)
            call mdet(F_t,detF_t)
            F_tau = ((detFc_tau/detF_tau)**third)*F_tau
            F_t = ((detFc_tau/detF_tau)**third)*F_t
         endif
         call mdet(F_tau,detF)


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration at this integ. point 
         !
         call integ(props,nprops,dtime,pnewdt,
     +        F_tau,mu_tau,omgP_tau,omgN_tau,phi_tau,omgH_tau,    
     +        zeta_t,theta0,
     +        T_tau,SpTanMod,
     +        zeta_tau,dZdt,DzetaDmu,DzetadotDmu,
     +        mobS,dmobSdmu,dmobSdJ,Omega,
     +        mobP,dmobPdomgP,mobN,dmobNdomgN,mobH,dmobHdomgH,
     +        cP_t,cN_t,cH_t,cRf,cRf_dissoc,
     +        cP_tau,cN_tau,dcPdomgP,dcNdomgN,cH_tau,dcHdomgH,
     +        dcPdmu,dcNdmu,dmobPdmu,dmobNdmu,dcHdmu,dmobHdmu,
     +        dcPdphi,dcNdphi,dmobPdphi,dmobNdphi,dcHdphi,dmobHdphi,  
     +        SpumuMod,SpmuumodFac)
         
         
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         cS_tau = (one/zeta_tau-one)/Omega
         
         Ka_COOH = 10.d0**(-4.5d0)*1.d3 
         cH_eff = cH_tau*zeta_tau/(one-zeta_tau)
         cRf_dissoc = cRf/(one+cH_eff/Ka_COOH)
         
         
         rR_H = 1/dtime*(cRf_dissoc-svars(5+jj))
         dcRf_dissdcH_eff = -cRf_dissoc/(Ka_COOH+cH_eff)
         dcH_effdomgH = dcHdomgH*zeta_tau/(one-zeta_tau)
         dcH_effdphi = dcHdphi*zeta_tau/(one-zeta_tau)
         dcH_effdmu = dcHdmu*zeta_tau/(one-zeta_tau)
         dcH_effdzeta = cH_tau*one/(one-zeta_tau)**two
         
         !cbarP_tau = cbarP_t
         !cbarN_tau = cbarN_t
         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         svars(1+jj) = zeta_tau
         svars(2+jj) = cP_tau
         svars(3+jj) = cN_tau
         svars(4+jj) = cH_tau
         svars(5+jj) = cRf_dissoc
         
         cPdot = (cP_tau-cP_t)/dtime
         cNdot = (cN_tau-cN_t)/dtime
         cHdot= (cH_tau-cH_t)/dtime
         
         
         jj = jj + nlSdv ! setup for the next intPt


         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         !
         globalSdv(jelem,intPt,1) = one/zeta_tau   ! polymer volume fraction
         globalSdv(jelem,intPt,2) = cP_tau   ! cation conc.
         globalSdv(jelem,intPt,3) = cN_tau   ! anion conc.
         globalSdv(jelem,intPt,4) = Kbulk*(dlog(detF*zeta_tau))
         globalSdv(jelem,intPt,5) = zeta_tau
         globalSdv(jelem,intPt,6) = cNdot
         globalSdv(jelem,intPt,7) = cH_tau
         
         

         ! Time stepping algorithim based on the constitutive response.
         !  Here based on the change in the polymer volume fraction change.
         !
         zetaLmt = 0.1d0
         umeror = dabs((zeta_tau - zeta_t)/zetaLmt)
         if(umeror.le.half) then
            pnewdt = 1.5d0
         elseif((umeror.gt.half).and.(umeror.le.0.8d0)) then
            pnewdt = 1.25d0
         elseif((umeror.gt.0.8d0).and.(umeror.le.1.25d0)) then
            pnewdt = 0.75d0
         else
            pnewdt = half
            write(*,*)'Problem here'
         endif

      AR = one
         ! Compute/update the displacement residual vector
         !
         ! Compute/update the displacement residual vector
         !
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(1,2)
         !
         Bmat = zero
         do kk=1,nNode
            Bmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Bmat(2,2+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,1+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,2+nDim*(kk-1)) = dshC(kk,1)
         enddo
         !
         BodyForceRes = zero
         do kk=1,nNode
            BodyForceRes(1+nDim*(kk-1),1) = sh(kk)*body(1)
            BodyForceRes(2+nDim*(kk-1),1) = sh(kk)*body(2)
         enddo
         !
         Ru = Ru + detmapJC*w(intpt)*AR*
     +        (
     +        -matmul(transpose(Bmat),Smat)
     +        + BodyForceRes
     +        )



         ! Compute/update the chemical potential residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo
         !
         ResFac = (dZdt)/(detF*Omega*zeta_tau*zeta_tau)
         !
         Rmu = Rmu + detmapJC*w(intpt)*AR*
     +        (
     +        transpose(Nvec)*ResFac 
     +        - mobS*matmul(dshC,dMUdX)
     +        - mobS/cS_tau*cP_tau*matmul(dshC,domgPdx)
     +        - mobS/cS_tau*cN_tau*matmul(dshC,domgNdx)  
     +        - mobS/cS_tau*cH_tau*matmul(dshC,domgHdx)     
     +        )
         
         
         ResFac = -cPdot/detF
         !
         RomgP  = RomgP + detmapJC*w(intpt)*AR*
     +        (
     +        transpose(Nvec)*ResFac 
     +        - mobP*matmul(dshC,domgPdX)
     +        - mobS*(cP_tau/cS_tau)**two*matmul(dshC,domgPdX)
     +        - mobS*(cP_tau/cS_tau)*
     +                matmul(dshC,dmudX)
     +        - mobS*(cP_tau*cN_tau/cS_tau**two)*
     +                  matmul(dshC,domgNdX)
     +        - mobS*(cP_tau*cH_tau/cS_tau**two)*
     +                  matmul(dshC,domgHdX)
     +        )
         !RomgP = zero
         
         ResFac = -cNdot/detF
         !
         RomgN  = RomgN + detmapJC*w(intpt)*AR*
     +        (
     +        transpose(Nvec)*ResFac 
     +        - mobN*matmul(dshC,domgNdX)
     +        - mobS*(cN_tau/cS_tau)**two*matmul(dshC,domgNdX)
     +        - mobS*(cN_tau/cS_tau)*
     +                matmul(dshC,dmudX)
     +        - mobS*(cP_tau*cN_tau/cS_tau**two)*
     +                  matmul(dshC,domgPdX)    
     +        - mobS*(cH_tau*cN_tau/cS_tau**two)*
     +                  matmul(dshC,domgHdX)    
     +        )
         !RomgN = zero
         ResFac = -cHdot/detF + rR_H/detF
         !
         RomgH  = RomgH + detmapJC*w(intpt)*AR*
     +        (
     +        transpose(Nvec)*ResFac 
     +        - mobH*matmul(dshC,domgHdX)
     +        - mobS*(cH_tau/cS_tau)**two*matmul(dshC,domgHdX)
     +        - mobS*(cH_tau/cS_tau)*
     +                matmul(dshC,dmudX)
     +        - mobS*(cP_tau*cH_tau/cS_tau**two)*
     +                  matmul(dshC,domgPdX)    
     +        - mobS*(cN_tau*cH_tau/cS_tau**two)*
     +                  matmul(dshC,domgNdX)   
     +        ) 
         
         
         
         ResFac = zP_mob*cP_tau+zN_mob*cN_tau
     +     +  zf*cRf_dissoc+ zH_mob*(cH_tau)
!     +    - zH_mob*(1.d-8/cH_tau)*(one/zeta_tau-one)**two
!     +      + zf*cRf_dissoc+ zH_mob*(cH_tau)- zH_mob*(1.d-8/cH_tau)
         !
         Rphi  = Rphi + detmapJC*w(intpt)*AR*
     +        (
     +        -transpose(Nvec)*ResFac
     +        )
!          Rphi  = Rphi + detmapJC*w(intpt)*
!     +        (
!     +        - matmul(dshC,dphidX)
!     +        )
         !Rphi = zero
         



         ! Compute/update the displacement tangent matrix
         !
         Gmat = zero
         do kk=1,nNode
            Gmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Gmat(2,2+nDim*(kk-1)) = dshC(kk,1)
            Gmat(3,1+nDim*(kk-1)) = dshC(kk,2)
            Gmat(4,2+nDim*(kk-1)) = dshC(kk,2)
         enddo

         G0mat = zero
         do kk=1,nNode
            G0mat(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(3,1+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(4,2+nDim*(kk-1)) = dshC0(kk,2)
         enddo

         Amat = zero
         Amat(1,1) = SpTanMod(1,1,1,1)
         Amat(1,2) = SpTanMod(1,1,2,1)
         Amat(1,3) = SpTanMod(1,1,1,2)
         Amat(1,4) = SpTanMod(1,1,2,2)
         Amat(2,1) = SpTanMod(2,1,1,1)
         Amat(2,2) = SpTanMod(2,1,2,1)
         Amat(2,3) = SpTanMod(2,1,1,2)
         Amat(2,4) = SpTanMod(2,1,2,2)
         Amat(3,1) = SpTanMod(1,2,1,1)
         Amat(3,2) = SpTanMod(1,2,2,1)
         Amat(3,3) = SpTanMod(1,2,1,2)
         Amat(3,4) = SpTanMod(1,2,2,2)
         Amat(4,1) = SpTanMod(2,2,1,1)
         Amat(4,2) = SpTanMod(2,2,2,1)
         Amat(4,3) = SpTanMod(2,2,1,2)
         Amat(4,4) = SpTanMod(2,2,2,2)


         Qmat = zero
         Qmat(1,1) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
         Qmat(2,1) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
         Qmat(3,1) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
         Qmat(4,1) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
         Qmat(1,4) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
         Qmat(2,4) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
         Qmat(3,4) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
         Qmat(4,4) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
            

         if((nNode.eq.4).and.(nInt.eq.4)) then
            !
            ! This is the tangent using the F-bar method with the
            !  4 node fully integrated linear element
            !
            Kuu = Kuu + detMapJC*w(intpt)*AR*
     +           (
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           + matmul(transpose(Gmat),matmul(Qmat,
     +           (G0mat-Gmat)))
     +           )
         else
            !
            ! This is the tangent not using the F-bar method with all
            !  other elements
            !
            Kuu = Kuu + detMapJC*w(intpt)*AR*
     +           (
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           )
       endif
            

         ! Compute/update the chemical potential tangent matrix
         !
         TanFac = (one/(detF*Omega*zeta_tau**two))*
     +        (two*(dZdt/zeta_tau)*DzetaDmu - DzetadotDmu)
         !
         Kmumu = Kmumu + detmapJC*w(intPt)*AR*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        + mobS*matmul(dshC,transpose(dshC))
     +        + dmobSdmu*matmul(matmul(dshC,dMUdX),Nvec)
     +        + mobS/cS_tau*dcPdmu*matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS/cS_tau*dcNdmu*matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS/cS_tau*dcHdmu*matmul(matmul(dshC,domgHdX),Nvec)
     +        )
        
         
         
         
         TanFac = 1/dtime*dcPdomgP/detF
         !
         KomgPomgP = KomgPomgP + detmapJC*w(intPt)*AR*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        + mobP*matmul(dshC,transpose(dshC))
     +        + dmobPdomgP*matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cP_tau/cS_tau)**two*matmul(dshC,transpose(dshC))
     +        + mobS*(two*cP_tau/cS_tau**two)*dcPdomgP
     +                                *matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS/cS_tau*dcPdomgP*
     +                                 matmul(matmul(dshC,dmudX),Nvec)
     +        + mobS*(cN_tau/cS_tau**two)*dcPdomgP*
     +                                 matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(cH_tau/cS_tau**two)*dcPdomgP*
     +                                matmul(matmul(dshC,domgHdX),Nvec)
     +        )
        
         
          
         
        TanFac = 1/dtime*dcNdomgN/detF
         !
         KomgNomgN = KomgNomgN + detmapJC*w(intPt)*AR*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        + mobN*matmul(dshC,transpose(dshC))
     +        + dmobNdomgN*matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(cN_tau/cS_tau)**two*matmul(dshC,transpose(dshC))
     +        + mobS*(two*cN_tau/cS_tau**two)*dcNdomgN
     +                                *matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS/cS_tau*dcNdomgN*
     +                                 matmul(matmul(dshC,dmudX),Nvec)
     +        + mobS*(cP_tau/cS_tau**two)*dcNdomgN*
     +                                 matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cH_tau/cS_tau**two)*dcNdomgN*
     +                                matmul(matmul(dshC,domgHdX),Nvec)
     +        )
         
         TanFac = 1/dtime*dcHdomgH/detF 
     +           - 1/dtime*dcRf_dissdcH_eff*dcH_effdomgH/detF
         !
         KomgHomgH = KomgHomgH + detmapJC*w(intPt)*AR*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        + mobH*matmul(dshC,transpose(dshC))
     +        + dmobHdomgH*matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS*(cH_tau/cS_tau)**two*matmul(dshC,transpose(dshC))
     +        + mobS*(two*cN_tau/cS_tau**two)*dcHdomgH
     +                              *matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS/cS_tau*dcHdomgH*
     +                                 matmul(matmul(dshC,dmudX),Nvec)
     +        + mobS*(cP_tau/cS_tau**two)*dcHdomgH*
     +                                 matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cN_tau/cS_tau**two)*dcHdomgH*
     +                                matmul(matmul(dshC,domgNdX),Nvec)
     +        )
         
         
         
         TanFac = zP_mob*dcPdphi+zN_mob*dcNdphi
     +-zf*cRf_dissoc/(Ka_COOH+cH_eff)*dcHdphi*zeta_tau/(one-zeta_tau)
     +           +zH_mob*dcHdphi 
!     +    +zH_mob/cH_tau**two*1.d-8*dcHdphi*(one/zeta_tau-one)**two
 
         !
         Kphiphi = Kphiphi + detmapJC*w(intPt)*AR*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        )
                  


        TanFac = 1/dtime*dcPdphi/detF
         !
         KomgPphi = KomgPphi + detmapJC*w(intPt)*AR*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        + dmobPdphi*matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(two*cP_tau/cS_tau**two)*dcPdphi
     +                           *matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS/cS_tau*dcPdphi*
     +                                 matmul(matmul(dshC,dmudX),Nvec)
     +        + mobS*(cN_tau/cS_tau**two)*dcPdphi*
     +                                 matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(cP_tau/cS_tau**two)*dcNdphi*
     +                                 matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(cH_tau/cS_tau**two)*dcPdphi*
     +                                matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS*(cP_tau/cS_tau**two)*dcHdphi*
     +                                matmul(matmul(dshC,domgHdX),Nvec)
     +        )
         
         
         TanFac = 1/dtime*dcNdphi/detF
         !
         KomgNphi = KomgNphi + detmapJC*w(intPt)*AR*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        + dmobNdphi*matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(two*cN_tau/cS_tau**two)*dcNdphi
     +                           *matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS/cS_tau*dcNdphi*
     +                                 matmul(matmul(dshC,dmudX),Nvec)
     +        + mobS*(cP_tau/cS_tau**two)*dcNdphi*
     +                                 matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cN_tau/cS_tau**two)*dcPdphi*
     +                                 matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cH_tau/cS_tau**two)*dcNdphi*
     +                                matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS*(cN_tau/cS_tau**two)*dcHdphi*
     +                                matmul(matmul(dshC,domgHdX),Nvec)
     +        )
         
         TanFac = 1/dtime*dcHdphi/detF
     +           - 1/dtime*dcRf_dissdcH_eff*dcH_effdphi/detF
         !
         KomgHphi = KomgHphi + detmapJC*w(intPt)*AR*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        + dmobHdphi*matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS*(two*cH_tau/cS_tau**two)*dcHdphi
     +                           *matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS/cS_tau*dcHdphi*
     +                                 matmul(matmul(dshC,dmudX),Nvec)
     +        + mobS*(cP_tau/cS_tau**two)*dcHdphi*
     +                                 matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cH_tau/cS_tau**two)*dcPdphi*
     +                                 matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cN_tau/cS_tau**two)*dcNdphi*
     +                                matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(cH_tau/cS_tau**two)*dcNdphi*
     +                                matmul(matmul(dshC,domgNdX),Nvec)
     +        )
         
         
         TanFac = zP_mob*dcPdomgP
         !
         KphiomgP = KphiomgP + detmapJC*w(intPt)*AR*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        )
         
          
          TanFac = zN_mob*dcNdomgN
         !
         KphiomgN = KphiomgN + detmapJC*w(intPt)*AR*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        )
        
         TanFac = zH_mob*dcHdomgH 
     +       -zf*cRf_dissoc/(Ka_COOH+cH_eff)*dcHdomgH
     +                             *zeta_tau/(one-zeta_tau)
!     +    +zH_mob/cH_tau**two*1.d-8*dcHdomgH*(one/zeta_tau-one)**two
         !
         KphiomgH = KphiomgH + detmapJC*w(intPt)*AR*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        )
         
         KmuomgP = KmuomgP + detmapJC*w(intPt)*AR*
     +        ( zero
     +        + mobS/cS_tau*cP_tau*matmul(dshC,transpose(dshC))
     +        + mobS/cS_tau*dcPdomgP*
     +                            matmul(matmul(dshC,domgPdX),Nvec)
     +        )
        
         
         KmuomgN = KmuomgN + detmapJC*w(intPt)*AR*
     +        ( zero
     +        + mobS/cS_tau*cN_tau*matmul(dshC,transpose(dshC))
     +        + mobS/cS_tau*dcNdomgN*
     +                            matmul(matmul(dshC,domgNdX),Nvec)
     +        )
         
         
         KmuomgH = KmuomgH + detmapJC*w(intPt)*AR*
     +        ( zero
     +        + mobS/cS_tau*cH_tau*matmul(dshC,transpose(dshC))
     +        + mobS/cS_tau*dcHdomgH*
     +                            matmul(matmul(dshC,domgHdX),Nvec)
     +        )
         
         
         Kmuphi = Kmuphi + detmapJC*w(intPt)*AR*
     +        ( zero
     +        + mobS/cS_tau*dcPdphi*
     +                            matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS/cS_tau*dcNdphi*
     +                            matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS/cS_tau*dcHdphi*
     +                            matmul(matmul(dshC,domgHdX),Nvec)
     +        )
         
         
         
         TanFac = 1/dtime*dcPdmu/detF
         KomgPmu = KomgPmu + detmapJC*w(intpt)*AR*
     +        ( zero
     +        + TanFac*matmul(transpose(Nvec),Nvec)
     +        + dmobPdmu*matmul(matmul(dshC,domgPdX),Nvec)
     +        - dmobSdmu*(cP_tau/cS_tau)**two*
     +                            matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(two*cP_tau/cS_tau**two)*dcPdmu*
     +                            matmul(matmul(dshC,domgPdX),Nvec)     
     +        + mobS/cS_tau*dcPdmu*
     +                            matmul(matmul(dshC,dmudX),Nvec)
     +        + mobS*(cP_tau/cS_tau)*
     +              matmul(dshC,transpose(dshC))
     +        - dmobSdmu*(cP_tau*cN_tau/cS_tau**two)*
     +                            matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(cN_tau/cS_tau**two)*dcPdmu*
     +                            matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(cP_tau/cS_tau**two)*dcNdmu*
     +                            matmul(matmul(dshC,domgNdX),Nvec)
     +        - dmobSdmu*(cP_tau*cH_tau/cS_tau**two)*
     +                            matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS*(cH_tau/cS_tau**two)*dcPdmu*
     +                            matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS*(cP_tau/cS_tau**two)*dcHdmu*
     +                            matmul(matmul(dshC,domgHdX),Nvec)
     +        ) 
         
         
         
         
         TanFac = 1/dtime*dcNdmu/detF
          KomgNmu = KomgNmu + detmapJC*w(intpt)*AR*
     +        ( zero
     +        + TanFac*matmul(transpose(Nvec),Nvec)
     +        + dmobNdmu*matmul(matmul(dshC,domgNdX),Nvec)
     +        - dmobSdmu*(cN_tau/cS_tau)**two*
     +                            matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(two*cN_tau/cS_tau**two)*dcNdmu*
     +                            matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS/cS_tau*dcNdmu*
     +                            matmul(matmul(dshC,dmudX),Nvec)   
     +        + mobS*(cN_tau/cS_tau)*
     +              matmul(dshC,transpose(dshC))
     +        - dmobSdmu*(cN_tau*cP_tau/cS_tau**two)*
     +                            matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cN_tau/cS_tau**two)*dcPdmu*
     +                            matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cP_tau/cS_tau**two)*dcNdmu*
     +                            matmul(matmul(dshC,domgPdX),Nvec)
     +        - dmobSdmu*(cN_tau*cH_tau/cS_tau**two)*
     +                            matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS*(cH_tau/cS_tau**two)*dcNdmu*
     +                            matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS*(cN_tau/cS_tau**two)*dcHdmu*
     +                            matmul(matmul(dshC,domgHdX),Nvec)
     +        ) 
          
          TanFac = 1/dtime*dcHdmu/detF
     +           - 1/dtime*dcRf_dissdcH_eff*dcH_effdmu/detF
     +           - 1/dtime*dcRf_dissdcH_eff*dcH_effdzeta*dzetadmu/detF
         KomgHmu = KomgHmu + detmapJC*w(intpt)*AR*
     +        ( zero
     +        + TanFac*matmul(transpose(Nvec),Nvec)
     +        + dmobHdmu*matmul(matmul(dshC,domgHdX),Nvec)
     +        - dmobSdmu*(cH_tau/cS_tau)**two*
     +                            matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS*(two*cH_tau/cS_tau**two)*dcHdmu*
     +                            matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS/cS_tau*dcHdmu*
     +                            matmul(matmul(dshC,dmudX),Nvec)   
     +        + mobS*(cH_tau/cS_tau)*
     +              matmul(dshC,transpose(dshC))
     +        - dmobSdmu*(cH_tau*cP_tau/cS_tau**two)*
     +                            matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cH_tau/cS_tau**two)*dcPdmu*
     +                            matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cP_tau/cS_tau**two)*dcHdmu*
     +                            matmul(matmul(dshC,domgPdX),Nvec)
     +        - dmobSdmu*(cH_tau*cN_tau/cS_tau**two)*
     +                            matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(cN_tau/cS_tau**two)*dcHdmu*
     +                            matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(cH_tau/cS_tau**two)*dcNdmu*
     +                            matmul(matmul(dshC,domgNdX),Nvec)
     +        ) 
         
          
          
         
         KomgPomgN = KomgPomgN + detmapJC*w(intPt)*AR*
     +        ( zero
     +        + mobS*(cP_tau/cS_tau**two)*dcNdomgN*
     +                            matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(cP_tau*cN_tau/cS_tau**two)*
     +                            matmul(dshC,transpose(dshC))    
     +        )
         
         KomgPomgH = KomgPomgH + detmapJC*w(intPt)*AR*
     +        ( zero
     +        + mobS*(cP_tau/cS_tau**two)*dcHdomgH*
     +                            matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS*(cP_tau*cH_tau/cS_tau**two)*
     +                            matmul(dshC,transpose(dshC))    
     +        )
         
         
         KomgNomgP = KomgNomgP + detmapJC*w(intPt)*AR*
     +        ( zero
     +        + mobS*(cN_tau/cS_tau**two)*dcPdomgP*
     +                            matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cN_tau*cP_tau/cS_tau**two)*
     +                            matmul(dshC,transpose(dshC))    
     +        )
         
         KomgNomgH = KomgNomgH + detmapJC*w(intPt)*AR*
     +        ( zero
     +        + mobS*(cN_tau/cS_tau**two)*dcHdomgH*
     +                            matmul(matmul(dshC,domgHdX),Nvec)
     +        + mobS*(cN_tau*cH_tau/cS_tau**two)*
     +                            matmul(dshC,transpose(dshC))    
     +        )
         
         KomgHomgP = KomgHomgP + detmapJC*w(intPt)*AR*
     +        ( zero
     +        + mobS*(cH_tau/cS_tau**two)*dcPdomgP*
     +                            matmul(matmul(dshC,domgPdX),Nvec)
     +        + mobS*(cH_tau*cP_tau/cS_tau**two)*
     +                            matmul(dshC,transpose(dshC))    
     +        )
         
         KomgHomgN = KomgHomgN + detmapJC*w(intPt)*AR*
     +        ( zero
     +        + mobS*(cH_tau/cS_tau**two)*dcNdomgN*
     +                            matmul(matmul(dshC,domgNdX),Nvec)
     +        + mobS*(cH_tau*cN_tau/cS_tau**two)*
     +                            matmul(dshC,transpose(dshC))    
     +        )
         
         
         TanFac = zP_mob*dcPdmu+zN_mob*dcNdmu
!     +           - zf*cRf_dissoc/(Kb_DMAEA+cH_tau)*dcHdmu
     +           +zH_mob*dcHdmu 
!     +    +zH_mob/cH_tau**two*1.d-8*dcHdmu*(one/zeta_tau-one)**two
!     +   -zH_mob/cH_tau*1.d-8*two*(one/zeta_tau-one)
!     +                              *(-one/zeta_tau**two)*dzetadmu 
         !
         Kphimu = Kphimu + detmapJC*w(intPt)*AR*
     +        ( zero 
     +        + TanFac*matmul(transpose(Nvec),Nvec)
     +        )
         
         
         

         ! Compute/update the chemical potential - displacement tangent matrix
         !  The F-bar method will have some effect, however we neglect that here.
         !
         Spmuumod = zero
         do i=1,nDim
            do k=1,nDim
               do l=1,nDim
                  Spmuumod(i,k,l) = Spmuumod(i,k,l)
     +                 + dMUdX(k,1)*SpmuumodFac(i,l)
               enddo
            enddo
         enddo
         
         
         AmatCU = zero
         AmatCU(1,1) = SpmuuMod(1,1,1)
         AmatCU(1,2) = SpmuuMod(1,2,1)
         AmatCU(1,3) = SpmuuMod(1,1,2)
         AmatCU(1,4) = SpmuuMod(1,2,2)
         AmatCU(2,1) = SpmuuMod(2,1,1)
         AmatCU(2,2) = SpmuuMod(2,2,1)
         AmatCU(2,3) = SpmuuMod(2,1,2)
         AmatCU(2,4) = SpmuuMod(2,2,2)
         !
         Kmuu = Kmuu - detMapJC*w(intpt)*AR*
     +        (
     +        matmul(matmul(dshC,AmatCU),Gmat)
     +        )
      
         KomgPu = KomgPu - detMapJC*w(intpt)*AR*
     +        ( zero
     +        + mobP/mobS*matmul(matmul(dshC,AmatCU),Gmat)
!     +        + matmul(matmul(dshC,AmatCU),Gmat)
     +        )
         
         KomgNu = KomgNu - detMapJC*w(intpt)*AR*
     +        ( zero
     +        + mobN/mobS*matmul(matmul(dshC,AmatCU),Gmat)
!     +        + matmul(matmul(dshC,AmatCU),Gmat)
     +        )
         
         KomgHu = KomgHu - detMapJC*w(intpt)*AR*
     +        ( zero
     +        + mobH/mobS*matmul(matmul(dshC,AmatCU),Gmat)
!     +        + matmul(matmul(dshC,AmatCU),Gmat)
     +        )
         


         ! Compute/update the displacement - chemical potential tangent matrix
         !  The F-bar method will have some effect, however we neglect that here.
         !
         AmatUC = zero
         AmatUC(1,1) = SpumuMod(1,1)
         AmatUC(2,1) = SpumuMod(2,2)
         AmatUC(3,1) = SpumuMod(1,2)
         !
         Kumu = Kumu + detMapJC*w(intpt)*AR*
     +        (
     +        matmul(matmul(transpose(Bmat),AmatUC),Nvec)
     +        )

      enddo
      !
      ! End the loop over body integration points
      !----------------------------------------------------------------
      !Rphi = zero
      !RomgP = zero
      !RomgN = zero
      !Rmu = zero
      !Ru = zero

      !----------------------------------------------------------------
      ! Start loop over surface flux terms here for flux
      !
      !
      if(ndload.gt.0) then
         !
         ! loop over faces and make proper modifications to
         !  residuals and tangents if needed
         !
         do i=1,ndload
            !
            ! 
            !
            face = jdltyp(i,1) ! label
            flux = adlmag(i,1) ! flux magnitude
            
            if(face.eq.1) then
               !
               ! flux on face 1 of the element
               !
               xLocal(1) = -dsqrt(one/three)
               yLocal(1) = -one
               xLocal(2) = dsqrt(one/three)
               yLocal(2) = -one
               w(1) = one
               w(2) = one
               !
               ! loop over integ points
               !
               do ii=1,2
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),face,coordsC,
     +                 sh,ds)


                   ! this job is axisymmetric
                     AR = zero
                     do nn=1,nNode
                        AR = AR + sh(nn)*coordsC(1,nn)
                     enddo
                     AR = two*Pi*AR
                  
                     mu_tau = zero
                     do k=1,nNode
                        mu_tau = mu_tau + muNew(k)*sh(k)
                     enddo
                     
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  
                  ! Modify the temperature residual, loop over nodes
                  !
                 
                  
                  flux_net = - K_rate * (flux - mu_tau)
                  Rmu(1,1) = Rmu(1,1) - w(ii)*ds*AR*sh(1)*flux_net
                  Rmu(2,1) = Rmu(2,1) - w(ii)*ds*AR*sh(2)*flux_net


                  ! Change on the tangent matrix
                  !
                  do k=1,nNode
                  Kmumu(1,k) = Kmumu(1,k) + w(ii)*ds*AR*sh(1)*
     +                            (K_rate)*sh(k)
                  
                  Kmumu(2,k) = Kmumu(2,k) + w(ii)*ds*AR*sh(2)*
     +                            (K_rate)*sh(k)
       
                  enddo
               enddo ! end loop over integ points 
               
               
               elseif(face.eq.2) then
               !
               ! flux on face 2 of the element
               !
               xLocal(1) = one
               yLocal(1) = -dsqrt(one/three)
               xLocal(2) = one
               yLocal(2) = dsqrt(one/three)
               w(1) = one
               w(2) = one
               !
               ! loop over integ points
               !
               do ii=1,2
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),face,coordsC,
     +                 sh,ds)


                   ! this job is axisymmetric
                     AR = zero
                     do nn=1,nNode
                        AR = AR + sh(nn)*coordsC(1,nn)
                     enddo
                     AR = two*Pi*AR
                     
                     mu_tau = zero
                     do k=1,nNode
                        mu_tau = mu_tau + muNew(k)*sh(k)
                     enddo
                  
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  
                  ! Modify the temperature residual, loop over nodes
                  !
                 
                  
                  flux_net = - K_rate * (flux - mu_tau)
                  Rmu(2,1) = Rmu(2,1) - w(ii)*ds*AR*sh(2)*flux_net
                  Rmu(3,1) = Rmu(3,1) - w(ii)*ds*AR*sh(3)*flux_net


                  ! Change on the tangent matrix
                  !
                  do k=1,nNode
                  Kmumu(2,k) = Kmumu(2,k) + w(ii)*ds*AR*sh(2)*
     +                            (K_rate)*sh(k)
                  
                  Kmumu(3,k) = Kmumu(3,k) + w(ii)*ds*AR*sh(3)*
     +                            (K_rate)*sh(k)
       
                  enddo
               enddo ! end loop over integ points 
               
               elseif(face.eq.3) then
               !
               ! flux on face 3 of the element
               !
               xLocal(1) = -dsqrt(one/three)
               yLocal(1) = one
               xLocal(2) = dsqrt(one/three)
               yLocal(2) = one
               w(1) = one
               w(2) = one
               !
               ! loop over integ points
               !
               do ii=1,2
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),face,coordsC,
     +                 sh,ds)


                   ! this job is axisymmetric
                     AR = zero
                     do nn=1,nNode
                        AR = AR + sh(nn)*coordsC(1,nn)
                     enddo
                     AR = two*Pi*AR
                     
                     mu_tau = zero
                     do k=1,nNode
                        mu_tau = mu_tau + muNew(k)*sh(k)
                     enddo
                  
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  
                  ! Modify the temperature residual, loop over nodes
                  !
                 
                  
                  flux_net = - K_rate * (flux - mu_tau)
                  Rmu(3,1) = Rmu(3,1) - w(ii)*ds*AR*sh(3)*flux_net
                  Rmu(4,1) = Rmu(4,1) - w(ii)*ds*AR*sh(4)*flux_net


                  ! Change on the tangent matrix
                  !
                  do k=1,nNode
                  Kmumu(3,k) = Kmumu(3,k) + w(ii)*ds*AR*sh(3)*
     +                            (K_rate)*sh(k)
                  
                  Kmumu(4,k) = Kmumu(4,k) + w(ii)*ds*AR*sh(4)*
     +                            (K_rate)*sh(k)
       
                  enddo
               enddo ! end loop over integ points
               
               elseif(face.eq.4) then
               !
               ! flux on face 4 of the element
               !
               xLocal(1) = -one
               yLocal(1) = dsqrt(one/three)
               xLocal(2) = -one
               yLocal(2) = -dsqrt(one/three)
               w(1) = one
               w(2) = one
               !
               ! loop over integ points
               !
               do ii=1,2
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),face,coordsC,
     +                 sh,ds)


                   ! this job is axisymmetric
                     AR = zero
                     do nn=1,nNode
                        AR = AR + sh(nn)*coordsC(1,nn)
                     enddo
                     AR = two*Pi*AR
                     
                     mu_tau = zero
                     do k=1,nNode
                        mu_tau = mu_tau + muNew(k)*sh(k)
                     enddo
                  
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  
                  ! Modify the temperature residual, loop over nodes
                  !
                 
                  
                  flux_net = -  K_rate * (flux - mu_tau)
                  Rmu(4,1) = Rmu(4,1) - w(ii)*ds*AR*sh(4)*flux_net
                  Rmu(1,1) = Rmu(1,1) - w(ii)*ds*AR*sh(1)*flux_net


                  ! Change on the tangent matrix
                  !
                  do k=1,nNode
                  Kmumu(4,k) = Kmumu(4,k) + w(ii)*ds*AR*sh(4)*
     +                            (K_rate)*sh(k)
                  
                  Kmumu(1,k) = Kmumu(1,k) + w(ii)*ds*AR*sh(1)*
     +                            (K_rate)*sh(k)
       
                  enddo
               enddo ! end loop over integ points
               
               
               elseif(face.eq.11) then
               !
               ! flux on face 1 of the element
               !
               xLocal(1) = -dsqrt(one/three)
               yLocal(1) = -one
               xLocal(2) = dsqrt(one/three)
               yLocal(2) = -one
               w(1) = one
               w(2) = one
               !
               ! loop over integ points
               !
               do ii=1,2
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),face,coordsC,
     +                 sh,ds)


                   ! this job is axisymmetric
                     AR = zero
                     do nn=1,nNode
                        AR = AR + sh(nn)*coordsC(1,nn)
                     enddo
                     AR = two*Pi*AR
                     
                     omgP_tau = zero
                     do k=1,nNode
                        omgP_tau = omgP_tau + omgPNew(k)*sh(k)
                     enddo
                  
                  
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  
                  ! Modify the temperature residual, loop over nodes
                  !
                 
                  
                  flux_net = - K1_rate * (flux - omgP_tau)
                  RomgP(1,1) = RomgP(1,1) - w(ii)*ds*AR*sh(1)*flux_net
                  RomgP(2,1) = RomgP(2,1) - w(ii)*ds*AR*sh(2)*flux_net


                  ! Change on the tangent matrix
                  !
                  do k=1,nNode
                  KomgPomgP(1,k) = KomgPomgP(1,k) + w(ii)*ds*AR*sh(1)*
     +                            (K1_rate)*sh(k)
                  
                  KomgPomgP(2,k) = KomgPomgP(2,k) + w(ii)*ds*AR*sh(2)*
     +                            (K1_rate)*sh(k)
       
                  enddo
               enddo ! end loop over integ points 
               
               
               elseif(face.eq.12) then
               !
               ! flux on face 2 of the element
               !
               xLocal(1) = one
               yLocal(1) = -dsqrt(one/three)
               xLocal(2) = one
               yLocal(2) = dsqrt(one/three)
               w(1) = one
               w(2) = one
               !
               ! loop over integ points
               !
               do ii=1,2
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),face,coordsC,
     +                 sh,ds)


                   ! this job is axisymmetric
                     AR = zero
                     do nn=1,nNode
                        AR = AR + sh(nn)*coordsC(1,nn)
                     enddo
                     AR = two*Pi*AR
                     
                     omgP_tau = zero
                     do k=1,nNode
                        omgP_tau = omgP_tau + omgPNew(k)*sh(k)
                     enddo
                  
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  
                  ! Modify the temperature residual, loop over nodes
                  !
                 
                  
                  flux_net = -K1_rate * (flux - omgP_tau)
                  RomgP(2,1) = RomgP(2,1) - w(ii)*ds*AR*sh(2)*flux_net
                  RomgP(3,1) = RomgP(3,1) - w(ii)*ds*AR*sh(3)*flux_net


                  ! Change on the tangent matrix
                  !
                  do k=1,nNode
                  KomgPomgP(2,k) = KomgPomgP(2,k) + w(ii)*ds*AR*sh(2)*
     +                            (K1_rate)*sh(k)
                  
                  KomgPomgP(3,k) = KomgPomgP(3,k) + w(ii)*ds*AR*sh(3)*
     +                            (K1_rate)*sh(k)
       
                  enddo
               enddo ! end loop over integ points 
               
               elseif(face.eq.13) then
               !
               ! flux on face 3 of the element
               !
               xLocal(1) = -dsqrt(one/three)
               yLocal(1) = one
               xLocal(2) = dsqrt(one/three)
               yLocal(2) = one
               w(1) = one
               w(2) = one
               !
               ! loop over integ points
               !
               do ii=1,2
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),face,coordsC,
     +                 sh,ds)


                   ! this job is axisymmetric
                     AR = zero
                     do nn=1,nNode
                        AR = AR + sh(nn)*coordsC(1,nn)
                     enddo
                     AR = two*Pi*AR
                     
                     omgP_tau = zero
                     do k=1,nNode
                        omgP_tau = omgP_tau + omgPNew(k)*sh(k)
                     enddo
                  
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  
                  ! Modify the temperature residual, loop over nodes
                  !
                 
                  
                  flux_net = - K1_rate * (flux - omgP_tau)
                  RomgP(3,1) = RomgP(3,1) - w(ii)*ds*AR*sh(3)*flux_net
                  RomgP(4,1) = RomgP(4,1) - w(ii)*ds*AR*sh(4)*flux_net


                  ! Change on the tangent matrix
                  !
                  do k=1,nNode
                  KomgPomgP(3,k) = KomgPomgP(3,k) + w(ii)*ds*AR*sh(3)*
     +                            (K1_rate)*sh(k)
                  
                  KomgPomgP(4,k) = KomgPomgP(4,k) + w(ii)*ds*AR*sh(4)*
     +                            (K1_rate)*sh(k)
       
                  enddo
               enddo ! end loop over integ points
               
               elseif(face.eq.14) then
               !
               ! flux on face 4 of the element
               !
               xLocal(1) = -one
               yLocal(1) = dsqrt(one/three)
               xLocal(2) = -one
               yLocal(2) = -dsqrt(one/three)
               w(1) = one
               w(2) = one
               !
               ! loop over integ points
               !
               do ii=1,2
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),face,coordsC,
     +                 sh,ds)


                   ! this job is axisymmetric
                     AR = zero
                     do nn=1,nNode
                        AR = AR + sh(nn)*coordsC(1,nn)
                     enddo
                     AR = two*Pi*AR
                     
                     omgP_tau = zero
                     do k=1,nNode
                        omgP_tau = omgP_tau + omgPNew(k)*sh(k)
                     enddo
                  
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  
                  ! Modify the temperature residual, loop over nodes
                  !
                 
                  
                  flux_net = - K1_rate * (flux - omgP_tau)
                  RomgP(4,1) = RomgP(4,1) - w(ii)*ds*AR*sh(4)*flux_net
                  RomgP(1,1) = RomgP(1,1) - w(ii)*ds*AR*sh(1)*flux_net


                  ! Change on the tangent matrix
                  !
                  do k=1,nNode
                  KomgPomgP(4,k) = KomgPomgP(4,k) + w(ii)*ds*AR*sh(4)*
     +                            (K1_rate)*sh(k)
                  
                  KomgPomgP(1,k) = KomgPomgP(1,k) + w(ii)*ds*AR*sh(1)*
     +                            (K1_rate)*sh(k)
       
                  enddo
               enddo ! end loop over integ points
               
               
               elseif(face.eq.21) then
               !
               ! flux on face 1 of the element
               !
               xLocal(1) = -dsqrt(one/three)
               yLocal(1) = -one
               xLocal(2) = dsqrt(one/three)
               yLocal(2) = -one
               w(1) = one
               w(2) = one
               !
               ! loop over integ points
               !
               do ii=1,2
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),face,coordsC,
     +                 sh,ds)


                   ! this job is axisymmetric
                     AR = zero
                     do nn=1,nNode
                        AR = AR + sh(nn)*coordsC(1,nn)
                     enddo
                     AR = two*Pi*AR
                     
                     omgN_tau = zero
                     do k=1,nNode
                        omgN_tau = omgN_tau + omgNNew(k)*sh(k)
                     enddo
                  
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  
                  ! Modify the temperature residual, loop over nodes
                  !
                 
                  
                  flux_net = -K2_rate * (flux - omgN_tau)
                  RomgN(1,1) = RomgN(1,1) - w(ii)*ds*AR*sh(1)*flux_net
                  RomgN(2,1) = RomgN(2,1) - w(ii)*ds*AR*sh(2)*flux_net


                  ! Change on the tangent matrix
                  !
                  do k=1,nNode
                  KomgNomgN(1,k) = KomgNomgN(1,k) + w(ii)*ds*AR*sh(1)*
     +                            (K2_rate)*sh(k)
                  
                  KomgNomgN(2,k) = KomgNomgN(2,k) + w(ii)*ds*AR*sh(2)*
     +                            (K2_rate)*sh(k)
       
                  enddo
               enddo ! end loop over integ points 
               
               
               elseif(face.eq.22) then
               !
               ! flux on face 2 of the element
               !
               xLocal(1) = one
               yLocal(1) = -dsqrt(one/three)
               xLocal(2) = one
               yLocal(2) = dsqrt(one/three)
               w(1) = one
               w(2) = one
               !
               ! loop over integ points
               !
               do ii=1,2
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),face,coordsC,
     +                 sh,ds)


                   ! this job is axisymmetric
                     AR = zero
                     do nn=1,nNode
                        AR = AR + sh(nn)*coordsC(1,nn)
                     enddo
                     AR = two*Pi*AR
                     
                     omgN_tau = zero
                     do k=1,nNode
                        omgN_tau = omgN_tau + omgNNew(k)*sh(k)
                     enddo
                  
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  
                  ! Modify the temperature residual, loop over nodes
                  !
                 
                  
                  flux_net = -K2_rate * (flux - omgN_tau)
                  RomgN(2,1) = RomgN(2,1) - w(ii)*ds*AR*sh(2)*flux_net
                  RomgN(3,1) = RomgN(3,1) - w(ii)*ds*AR*sh(3)*flux_net


                  ! Change on the tangent matrix
                  !
                  do k=1,nNode
                  KomgNomgN(2,k) = KomgNomgN(2,k) + w(ii)*ds*AR*sh(2)*
     +                            (K2_rate)*sh(k)
                  
                  KomgNomgN(3,k) = KomgNomgN(3,k) + w(ii)*ds*AR*sh(3)*
     +                            (K2_rate)*sh(k)
       
                  enddo
               enddo ! end loop over integ points 
               
               elseif(face.eq.23) then
               !
               ! flux on face 3 of the element
               !
               xLocal(1) = -dsqrt(one/three)
               yLocal(1) = one
               xLocal(2) = dsqrt(one/three)
               yLocal(2) = one
               w(1) = one
               w(2) = one
               !
               ! loop over integ points
               !
               do ii=1,2
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),face,coordsC,
     +                 sh,ds)


                   ! this job is axisymmetric
                     AR = zero
                     do nn=1,nNode
                        AR = AR + sh(nn)*coordsC(1,nn)
                     enddo
                     AR = two*Pi*AR
                     
                     omgN_tau = zero
                     do k=1,nNode
                        omgN_tau = omgN_tau + omgNNew(k)*sh(k)
                     enddo
                  
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  
                  ! Modify the temperature residual, loop over nodes
                  !
                 
                  
                  flux_net = - K2_rate * (flux - omgN_tau)
                  RomgN(3,1) = RomgN(3,1) - w(ii)*ds*AR*sh(3)*flux_net
                  RomgN(4,1) = RomgN(4,1) - w(ii)*ds*AR*sh(4)*flux_net


                  ! Change on the tangent matrix
                  !
                  do k=1,nNode
                  KomgNomgN(3,k) = KomgNomgN(3,k) + w(ii)*ds*AR*sh(3)*
     +                            (K2_rate)*sh(k)
                  
                  KomgNomgN(4,k) = KomgNomgN(4,k) + w(ii)*ds*AR*sh(4)*
     +                            (K2_rate)*sh(k)
       
                  enddo
               enddo ! end loop over integ points
               
               elseif(face.eq.24) then
               !
               ! flux on face 4 of the element
               !
               xLocal(1) = -one
               yLocal(1) = dsqrt(one/three)
               xLocal(2) = -one
               yLocal(2) = -dsqrt(one/three)
               w(1) = one
               w(2) = one
               !
               ! loop over integ points
               !
               do ii=1,2
                  
                  ! Compute shape functions, derivatives, and the 
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(xLocal(ii),yLocal(ii),face,coordsC,
     +                 sh,ds)


                   ! this job is axisymmetric
                     AR = zero
                     do nn=1,nNode
                        AR = AR + sh(nn)*coordsC(1,nn)
                     enddo
                     AR = two*Pi*AR
                  
                     omgN_tau = zero
                     do k=1,nNode
                        omgN_tau = omgN_tau + omgNNew(k)*sh(k)
                     enddo
                     
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  
                  ! Modify the temperature residual, loop over nodes
                  !
                 
                  
                  flux_net = -K2_rate * (flux - omgN_tau)
                  RomgN(4,1) = RomgN(4,1) - w(ii)*ds*AR*sh(4)*flux_net
                  RomgN(1,1) = RomgN(1,1) - w(ii)*ds*AR*sh(1)*flux_net


                  ! Change on the tangent matrix
                  !
                  do k=1,nNode
                  KomgNomgN(4,k) = KomgNomgN(4,k) + w(ii)*ds*AR*sh(4)*
     +                            (K2_rate)*sh(k)
                  
                  KomgNomgN(1,k) = KomgNomgN(1,k) + w(ii)*ds*AR*sh(1)*
     +                            (K2_rate)*sh(k)
       
                  enddo
               enddo ! end loop over integ points
               
         endif      
               
         enddo ! loop over ndload
               endif ! ndload.gt.0 or not
      !
      ! End loop over flux and traction terms
      !----------------------------------------------------------------
      !RomgP = zero
      !RomgN = zero
      !Rphi = zero
      !KomgNomgN = one
      !KomgPomgP = one
      !Kphiphi = one
      !KomgPmu = zero
      !KomgNmu = zero
      !KmuomgP = zero
      !KmuomgN = zero
      !Kmuphi = zero
      !Kphimu = zero
      !RomgH = zero
      !KomgHomgH = one

      if(kstep.ge.3)then
          Ru = Ru/1.d1
          Kuu = Kuu/1.d1
          Kumu = Kumu/1.d1
          Rmu = Rmu/1.d2
          Kmumu = Kmumu/1.d2
          Kmuu = Kmuu/1.d2
          KmuomgP = KmuomgP/1.d2
          KmuomgN = KmuomgN/1.d2
          KmuomgH = KmuomgH/1.d2
      endif         
               
      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,Rmu,RomgP,RomgN,Rphi,RomgH,
     +     Kuu,Kumu,Kmuu,Kmumu,KomgPomgP,KomgNomgN,Kphiphi,KomgHomgH,
     +     KomgPphi,KomgNphi,KphiomgP,KphiomgN,
     +     KmuomgN,KmuomgP,KomgPmu,KomgNmu,
     +     KomgPomgN,KomgNomgP,
     +     KomgHomgP,KomgHomgN,
     +     KomgPomgH,KomgNomgH,KmuomgH,
     +     KomgPu,KomgNu,Kphimu,Kmuphi,KomgHu,KomgHphi,KphiomgH,
     +     KomgHmu, 
     +     rhs,amatrx)
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine UPE4

************************************************************************

      subroutine integ(props,nprops,dtime,pnewdt,
     +     F_tau,mu_tau,omgP_tau,omgN_tau,phi_tau,omgH_tau,
     +     zeta_t,theta, 
     +     T_tau,SpTanMod, 
     +     zeta_tau,dZdt,DzetaDmu,DzetadotDmu,
     +     mobS,DmobsDmu,DmobsDJ,Omega,
     +     mobP,dmobPdomgP,mobN,dmobNdomgN,mobH,dmobHdomgH,
     +     cP_t,cN_t,cH_t,cRf,cRf_dissoc,
     +     cP_tau,cN_tau,dcPdomgP,dcNdomgN,cH_tau,dcHdomgH,
     +     dcPdmu,dcNdmu,dmobPdmu,dmobNdmu,dcHdmu,dmobHdmu,
     +     dcPdphi,dcNdphi,dmobPdphi,dmobNdphi,dcHdphi,dmobHdphi,
     +     SpUmuMod,SpmuUModFac)

      ! This subroutine computes everything required for the time integration
      ! of the problem.
      !
      ! Inputs:
      !  1) material parameters, props(nprops)
      !  2) time increment, dtime
      !  3) deformation gradient, F_tau(3,3)
      !  4) chemical potential, mu_tau
      !  5) old polymer volume fraction, zeta_t
      !  6) temperature, theta
      !
      ! Outputs:
      !  1) Cauchy stress, T_tau(3,3)
      !  2) spatial tangent modulus, SpTanMod(3,3,3,3)
      !  3) polymer volume fraction, zeta_tau
      !  4) time rate of polymer volume fraction, dZdt
      !  5) derivative of the zeta with mu, DzetaDmu
      !  6) derivative of the time rate of zeta with mu, DzetadotDmu
      !  7) scalar fluid permeability, mobS
      !  8) derivative of permeability with chemical potential, dmobSdmu
      !  9) volume of a mole of fluid, Omega
      ! 10) displacement - chemical potential modulus terms
      ! 11) chemical potential - displacement modulus terms

      implicit none

      integer i,j,k,l,m,n,nprops,nargs,stat
      parameter(nargs=9)

      ! Input and output variables
      real*8 props(nprops),F_tau(3,3),mu_tau,zeta_t,dtime,pnewdt
      real*8 omgP_tau,omgN_tau,phi_tau,omgH_tau
      
      real*8 zeta_tau,dZdt,DzetaDmu,DzetadotDmu,mobS
      real*8 DmobSDmu,DzetaDJ,DmobSDzeta,DmobSDJ
      real*8 mobP,dmobPdomgP,mobN,dmobNdomgN,mobH,dmobHdomgH
      real*8 cP_tau,cN_tau,dcPdomgP,dcNdomgN,cH_tau,dcHdomgH
      real*8 cP_t,cN_t,cH_t, cOH_t
      real*8 dcPdmu,dcNdmu,dmobPdmu,dmobNdmu,dcHdmu,dmobHdmu
      real*8 dcPdphi,dcNdphi,dmobPdphi,dmobNdphi,dcHdphi,dmobHdphi
      real*8 T_tau(3,3),spTanMod(3,3,3,3),SpUmuMod(3,3),SpmuUModFac(3,3)
      
      ! Properties
      real*8 theta,Gshear,Kbulk,chi,D,mu0,Omega,Rgas,Farad
      real*8 DiffP,cRmaxP,cP0,zP_mob
      real*8 DiffN,cRmaxN,cN0,zN_mob
      real*8 DiffH,cH0,zH_mob
      real*8 cRf,cRf_dissoc
      
      ! Misc. variables
      real*8 Iden(3,3)
      
      real*8 TR_tau(3,3),dTRdF(3,3,3,3)
      real*8 detF,FinvT(3,3)
      real*8 B_tau(3,3),trB_tau,C_tau(3,3),trC_tau,args(nargs),detFe
      real*8 deltaMU,dPdt_per,dPdt_m
      real*8 Finv(3,3)
      real*8 zeta_per,zeta_m
      real*8 detFs

      real*8 zero,one,two,three,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0,Farad=96485.33d0)


      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain material properties
      !
      Gshear = props(1)
      Kbulk  = props(2)
      chi    = props(3)
      D      = props(4)
      mu0    = props(5)
      Omega  = props(6)
      Rgas   = props(7)
      
      DiffP = props(11)
      cRmaxP= props(12)
      cP0   = props(13)
      zP_mob= props(14)
      
      DiffN = props(15)
      cRmaxN= props(16)
      cN0   = props(17)
      zN_mob= props(18)
      
      cH0 = props(21)
      zH_mob = props(22)
      DiffH = props(23)
      
      !chi = 5.d0*chi/4.d0-cRf_dissoc/cRf*chi/4.d0
      


      ! Compute the inverse of F, its determinant, and its transpose
      !
      call matInv3D(F_tau,Finv,detF,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero'
         pnewdt=0.37
      endif
      FinvT = transpose(Finv)


      ! Compute the left Cauchy-Green tensor and its trace
      !
      B_tau = matmul(F_tau,transpose(F_tau))
      trB_tau = B_tau(1,1) + B_tau(2,2) + B_tau(3,3)


      ! Compute the right Cauchy-Green tensor and its trace
      !
      C_tau = matmul(transpose(F_tau),F_tau)
      trC_tau = C_tau(1,1) + C_tau(2,2) + C_tau(3,3)


      cOH_t = 1.d-8/cH_t*(one/zeta_t-one)**two
      ! Compute the polymer volume fraction
      !
      args(1)  = mu_tau
      args(2)  = mu0
      args(3)  = Rgas
      args(4)  = theta
      args(5)  = chi
      args(6)  = Omega
      args(7)  = Kbulk
      args(8)  = detF
      args(9)  = cP_t + cN_t + cH_t + cOH_t
      call solveZeta(zeta_tau,args,nargs,zeta_t,pnewdt)
      !zeta_tau = zeta_t

      cP_tau = (one - zeta_tau)/(Omega*zeta_tau) !(one - zeta_tau)!one/Omega
     +           *dexp((omgP_tau-zP_mob*Farad*phi_tau)/(Rgas*theta))
      cN_tau = (one - zeta_tau)/(Omega*zeta_tau) !(one - zeta_tau)!one/Omega
     +           *dexp((omgN_tau-zN_mob*Farad*phi_tau)/(Rgas*theta))
      
      cH_tau = (one - zeta_tau)/(Omega*zeta_tau)*
     +         dexp((omgH_tau-zH_mob*Farad*phi_tau)/(Rgas*theta))
      
      dcPdomgP = cP_tau/(Rgas*theta)
      dcNdomgN = cN_tau/(Rgas*theta)
      dcPdphi  = -cP_tau/(Rgas*theta)*zP_mob*Farad
      dcNdphi  = -cN_tau/(Rgas*theta)*zN_mob*Farad
      
      dcHdomgH= cH_tau/(Rgas*theta)
      dcHdphi  = -cH_tau/(Rgas*theta)*zH_mob*Farad
      
      
      
      ! Compute the elastic volume ratio, detFe
      !
      detFe = detF*zeta_tau


      ! Compute the swelling volume ratio, detFs
      !
      detFs = one/zeta_tau


      ! Compute the time rate of the polymer volume fraction using
      !  a finite difference in time
      !
      dZdt = (zeta_tau - zeta_t)/dtime


      ! Compute the derivative of the polymer volume fraction with
      !  respect to the chemical potential.  Computed via implicit
      !  differentiation on the chemical potential equation.
      !
      DzetaDmu =  (one/(Rgas*theta))/
     +     (
     +     (one/(zeta_tau - one)) + one + two*chi*zeta_tau 
     +     - ((Omega*Kbulk)/(Rgas*theta*zeta_tau))
     +     + ((Omega*Kbulk)/(Rgas*theta*zeta_tau))*dlog(detF*zeta_tau)
     +     - Omega*(cP_t+cN_t+cH_t+cOH_t)/(one-zeta_tau)**two
     +     )
      

      dcPdmu = -1/zeta_tau*cP_tau*dzetadmu/(one-zeta_tau)
      dcNdmu = -1/zeta_tau*cN_tau*dzetadmu/(one-zeta_tau)
      dcHdmu= -1/zeta_tau*cH_tau*dzetadmu/(one-zeta_tau)
      
      

      ! Compute the derivative of the polymer volume fraction with
      !  respect to the chemical potential.  Computed via implicit
      !  differentiation on the chemical potential equation.
      !
      DzetaDJ = (
     +     (Omega*Kbulk)/(Rgas*theta*detF) 
     +     - ((Omega*Kbulk)/(Rgas*theta*detF))*dlog(detF*zeta_tau)
     +     )/
     +     (
     +     (one/(zeta_tau - one)) + one + two*chi*zeta_tau 
     +     - ((Omega*Kbulk)/(Rgas*theta*zeta_tau))
     +     + ((Omega*Kbulk)/(Rgas*theta*zeta_tau))*dlog(detF*zeta_tau)
     +     )


      ! Compute the perturbation on the chemical potential
      !
      if(dabs(mu_tau).gt.one) then
         deltaMU = dabs(mu_tau)*1.d-8
      else
         deltaMU = 1.d-8
      endif


      ! Compute a perturbed polymer volume fraction
      !
      args(1)  = mu_tau + deltaMU
      args(2)  = mu0
      args(3)  = Rgas
      args(4)  = theta
      args(5)  = chi
      args(6)  = Omega
      args(7)  = Kbulk
      args(8)  = detF
      args(9)  = cP_t + cN_t + cH_t + cOH_t
      call solveZeta(zeta_per,args,nargs,zeta_t,pnewdt)
      !zeta_tau = zeta_t

      ! Compute a perturbed polymer volume fraction
      !
      args(1)  = mu_tau - deltaMU
      args(2)  = mu0
      args(3)  = Rgas
      args(4)  = theta
      args(5)  = chi
      args(6)  = Omega
      args(7)  = Kbulk
      args(8)  = detF
      args(9)  = cP_t + cN_t + cH_t + cOH_t
      call solveZeta(zeta_m,args,nargs,zeta_t,pnewdt)
      !zeta_m = zeta_t

      ! Compute the derivative of the time rate of change of the 
      !  polymer volume fraction with respect to the chemical potential
      !
      dPdt_per = (zeta_per - zeta_t)/dtime
      dPdt_m   = (zeta_m - zeta_t)/dtime
      DzetadotDmu = 1/dtime*DzetaDmu!(dPdt_per - dPdt_m)/(two*deltaMU)


      ! Compute the fluid permeability at this integ. point
      !
      !mobS = D/(Vmol*Rgas*theta)
      !
      ! to do m = (D*cR)/(R*T), use the following line
      !mobS = (D*(one/zeta_tau - one))/(Vmol*Rgas*theta)
      !
      ! to do m = (D*c)/(R*T), use the following line
      mobS = (D*(one/zeta_tau - one))/(detF*Omega*Rgas*theta)
      mobP = (DiffP*cP_tau)/(detF*Rgas*theta)
      mobN = (DiffN*cN_tau)/(detF*Rgas*theta)
      mobH= (DiffH*cH_tau)/(detF*Rgas*theta)

      dmobPdomgP = (DiffP*dcPdomgP)/(Rgas*theta)/detF
      dmobNdomgN = (DiffN*dcNdomgN)/(Rgas*theta)/detF
      dmobHdomgH= (DiffH*dcHdomgH)/(Rgas*theta)/detF
      dmobPdphi  = (DiffP*dcPdphi)/(Rgas*theta)/detF
      dmobNdphi  = (DiffN*dcNdphi)/(Rgas*theta)/detF
      dmobHdphi = (DiffH*dcHdphi)/(Rgas*theta)/detF
      
      dmobPdmu  = (DiffP*dcPdmu)/(Rgas*theta)/detF
      dmobNdmu  = (DiffN*dcNdmu)/(Rgas*theta)/detF
      dmobHdmu = (DiffH*dcHdmu)/(Rgas*theta)/detF
      
      
      ! Compute the tangents of the fluid mobility
      !
      !DmDzeta = zero
      !
      ! to do m = (D*cR)/(R*T), use the following line
      !DmDzeta = -(D/(Vmol*zeta_tau*zeta_tau*Rgas*theta))
      !
      ! to do m = (D*c)/(R*T), use the following line
      DmobSDzeta = -(D/(detF*Omega*zeta_tau*zeta_tau*Rgas*theta))
      !
      DmobSDmu = DmobSDzeta*DzetaDmu
      DmobSDJ  = DmobSDzeta*DzetaDJ



      ! Compute the Cauchy stress
      !
      T_tau = (Gshear*(B_tau-Iden) + detFs*Kbulk*dlog(detFe)*Iden)/detF

      
      ! Compute the 1st Piola stress
      !
      TR_tau = Gshear*(F_tau - FinvT) + detFs*Kbulk*dlog(detFe)*FinvT


      ! Compute dTRdF, the so-called material tangent modulus
      !
      dTRdF = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3          
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
     +                 + Gshear*Iden(i,k)*Iden(j,l)
     +                 + Gshear*Finv(l,i)*Finv(j,k)
     +                 + detFs*Kbulk*Finv(j,i)*Finv(l,k)
     +                 - detFs*Kbulk*dlog(detFe)*Finv(l,i)*Finv(j,k)
               enddo
            enddo
         enddo
      enddo
      !
      ! Calculate the so-called spatial tangent modulus, based
      !  on the push forward of the material tangent modulus
      !
      SpTanMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) + 
     +                       (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo


      ! Compute the displacement - chemical potential modulus
      !
      SpUmuMod = (Kbulk/(detFe*zeta_tau))*Iden*DzetaDmu


      ! Compute the chemical potential - displacement modulus
      !
      SpmuUModFac = mobS*Iden
      

      return
      end subroutine integ

****************************************************************************

      subroutine solveZeta(root,args,nargs,rootOld,pnewdt)

      ! This subroutine will numerically solve for the polymer
      !  volume fraction based on the current osmotic pressure
      !  and the previous state.  See numerical recipies RTSAFE.

      implicit none

      integer maxit,j,nargs

      real*8 xacc,f,df,fl,fh,xl,xh,x1,x2,swap,root,dxold,one
      real*8 dx,args(nargs),zero,rootOld,temp,rootMax,rootMin
      real*8 pnewdt

      parameter(maxit=50)
      parameter(xacc=1.d-5,zero=0.d0,one=1.d0)

      rootMax = 0.9999d0 ! corresponds to nearly 100% dry polymer
      rootMin = 0.005d0   ! corresponds to nearly 100% fluid

      x1 = rootMin
      x2 = rootMax
      call zetaFunc(x1,FL,DF,args,nargs)
      call zetaFunc(x2,FH,DF,args,nargs)

      if(fl*fh.ge.zero) then
         root = rootOld
         write(*,*) 'FYI, root not bracketed on zeta'
         !write(*,*) 'fl=',fl
         !write(*,*) 'fh=',fh
         !write(*,*) 'rootOld=',rootOld
         !write(80,*) 'FYI, the root is not bracketed on zeta'
         !write(80,*) 'fl=',fl
         !write(80,*) 'fh=',fh
         !write(80,*) 'rootOld=',rootOld

         !write(*,*) 'mu =',args(1)
         !write(*,*) 'mu0=',args(2)
         !write(*,*) 'Rgas=',args(3)
         !write(*,*) 'theta=',args(4)
         !write(*,*) 'chi=',args(5)
         !write(*,*) 'Vmol=',args(6)
         !write(*,*) 'Kbulk=',args(7)
         !write(*,*) 'detF=',args(8)

         pnewdt = 0.15d0
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
      
      call zetaFunc(root,F,DF,args,nargs)

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
         call zetaFunc(root,F,DF,args,nargs)

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

      WRITE(*,'(/1X,A)') 'solveZeta EXCEEDING MAXIMUM ITERATIONS'
      WRITE(80,'(/1X,A)') 'solveZeta EXCEEDING MAXIMUM ITERATIONS'

      return
      end subroutine solveZeta

****************************************************************************

      subroutine zetaFunc(phi,f,df,args,nargs)

      ! This subroutine serves as the function we would like to solve for
      !  the polymer volume fraction by finding phi such that ``f=0''

      implicit none

      integer nargs,NeoHookean,Langevin,material
      parameter(NeoHookean=1,Langevin=2)

      real*8 args(nargs),f,df,mu,mu0,Rgas,theta,chi,Omega,Gshear,Kbulk
      real*8 detF,phi,RT,c_total

      real*8 zero,one,two,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0)

      
      ! Obtain relevant quantities
      !
      mu      = args(1)
      mu0     = args(2)
      Rgas    = args(3)
      theta   = args(4)
      chi     = args(5)
      Omega   = args(6)
      Kbulk   = args(7)
      detF    = args(8)
      c_total = args(9)


      ! Compute the useful quantity
      !
      RT = Rgas*theta


      ! Compute the residual
      !
      f = (mu0 - mu)/RT
     +     + dlog(one - phi) + phi + chi*phi*phi
     +     - ((Kbulk*Omega)/RT)*dlog(detF*phi)
     +     + ((Kbulk*Omega)/(two*RT))*(dlog(detF*phi)**two)
     +     - Omega*c_total/(one/phi-one)


      ! Compute the tangent
      !
      if(phi.gt.0.999d0) then
         df = zero
      else
         df = one - (one/(one - phi)) + two*chi*phi
     +        - (Kbulk*Omega)/(RT*phi)
     +        + ((Kbulk*Omega)/(RT*phi))*dlog(detF*phi)
     +        - Omega*c_total/(one-phi)**2.d0  
      endif


      return
      end subroutine zetaFunc

************************************************************************

      subroutine AssembleElement(nDim,nNode,ndofel,
     +     Ru,Rmu,RomgP,RomgN,Rphi,RomgH,
     +     Kuu,Kumu,Kmuu,Kmumu,KomgPomgP,KomgNomgN,Kphiphi,KomgHomgH,
     +     KomgPphi,KomgNphi,KphiomgP,KphiomgN,
     +     KmuomgN,KmuomgP,KomgPmu,KomgNmu,
     +     KomgPomgN,KomgNomgP,
     +     KomgHomgP,KomgHomgN,
     +     KomgPomgH,KomgNomgH,KmuomgH,
     +     KomgPu,KomgNu,Kphimu,Kmuphi,KomgHu,KomgHphi,KphiomgH,
     +     KomgHmu,
     +     rhs,amatrx)
      
      !
      ! Subroutine to assemble the local elements residual and tangent
      !

      implicit none

      integer i,j,k,l,m,n,A11,A12,B11,B12,nDim,nNode,nDofEl,nDofN

      real*8 Ru(nDim*nNode,1),Rmu(nNode,1),RomgP(nNode,1),RomgN(nNode,1)
      real*8 Rphi(nNode,1),RomgH(nNode,1),KomgHomgH(nNode,nNode)
      real*8 Kuu(nDim*nNode,nDim*nNode),Kphiphi(nNode,nNode)
      real*8 KomgPomgP(nNode,nNode),KomgNomgN(nNode,nNode)
      real*8 KomgPphi(nNode,nNode),KomgNphi(nNode,nNode)
      real*8 KphiomgP(nNode,nNode),KphiomgN(nNode,nNode)
      real*8 Kmumu(nNode,nNode),Kumu(nDim*nNode,nNode),rhs(ndofel,1)
      real*8 KmuomgN(nNode,nNode),KmuomgP(nNode,nNode)
      real*8 KomgPmu(nNode,nNode),KomgNmu(nNode,nNode)
      real*8 KomgPomgN(nNode,nNode),KomgNomgP(nNode,nNode)
      real*8 KomgPu(nNode,nDim*nNode),KomgNu(nNode,nDim*nNode) 
      real*8 Kphimu(nNode,nNode),Kmuphi(nNode,nNode)
      real*8 Kmuu(nNode,nDim*nNode),amatrx(ndofel,ndofel)
      real*8 KomgHu(nNode,nDim*nNode),KomgHphi(nNode,nNode)
      real*8 KphiomgH(nNode,nNode),KomgHmu(nNode,nNode)
      real*8 KomgPomgH(nNode,nNode),KomgNomgH(nNode,nNode)
      real*8 KmuomgH(nNode,nNode)
      real*8 KomgHomgP(nNode,nNode),KomgHomgN(nNode,nNode)


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
            ! electrical potential
            !
            rhs(A11+2,1) = Rphi(i,1)
            ! chemical potential
            !
            rhs(A11+3,1) = RomgH(i,1)
            !
            rhs(A11+4,1) = Rmu(i,1)
            !
            rhs(A11+5,1) = RomgP(i,1)
            !
            rhs(A11+6,1) = RomgN(i,1)
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
               ! electrical potential
               !
               amatrx(A11+2,B11+2) = Kphiphi(i,j)
               !
               ! chemical potential
               !
               amatrx(A11+3,B11+3) = KomgHomgH(i,j)
               !
               amatrx(A11+4,B11+4) = Kmumu(i,j)
               !
               amatrx(A11+5,B11+5) = KomgPomgP(i,j)
               !
               amatrx(A11+6,B11+6) = KomgNomgN(i,j)
               !
               ! displacement - chemical potential
               !
               amatrx(A11,B11+4) = Kumu(A12,j)
               amatrx(A11+1,B11+4) = Kumu(A12+1,j)
               !
               ! chemical potential - displacement
               !
               !
               amatrx(A11+3,B11) = KomgHu(i,B12)
               amatrx(A11+3,B11+1) = KomgHu(i,B12+1)
               !
               amatrx(A11+4,B11) = Kmuu(i,B12)
               amatrx(A11+4,B11+1) = Kmuu(i,B12+1)
               !
               !
               amatrx(A11+5,B11) = KomgPu(i,B12)
               amatrx(A11+5,B11+1) = KomgPu(i,B12+1)
               !
               !
               amatrx(A11+6,B11) = KomgNu(i,B12)
               amatrx(A11+6,B11+1) = KomgNu(i,B12+1)
               !
               ! echemical potential - electrical potential
               !
               amatrx(A11+3,B11+2) = KomgHphi(i,j)
               !
               amatrx(A11+5,B11+2) = KomgPphi(i,j)
               !
               amatrx(A11+6,B11+2) = KomgNphi(i,j)
               !
               !
               amatrx(A11+2,B11+4) = Kphimu(i,j)
               !
               amatrx(A11+4,B11+2) = Kmuphi(i,j)
               !
               amatrx(A11+2,B11+3) = KphiomgH(i,j)
               !
               amatrx(A11+2,B11+5) = KphiomgP(i,j)
               !
               amatrx(A11+2,B11+6) = KphiomgN(i,j)
               
               !
               amatrx(A11+4,B11+5) = KmuomgP(i,j)
               !
               amatrx(A11+4,B11+6) = KmuomgN(i,j)
               !
               amatrx(A11+4,B11+3) = KmuomgH(i,j)
               !
               amatrx(A11+3,B11+4) = KomgHmu(i,j)
               !
               amatrx(A11+5,B11+4) = KomgPmu(i,j)
               !
               amatrx(A11+6,B11+4) = KomgNmu(i,j)
               !
               amatrx(A11+5,B11+6) = KomgPomgN(i,j)
               !
               amatrx(A11+6,B11+5) = KomgNomgP(i,j)
               !
               amatrx(A11+5,B11+3) = KomgPomgH(i,j)
               !
               amatrx(A11+6,B11+3) = KomgNomgH(i,j)
               !
               amatrx(A11+3,B11+5) = KomgHomgP(i,j)
               !
               amatrx(A11+3,B11+6) = KomgHomgN(i,j)
            enddo
         enddo
         !
      else
         write(*,*) 'How did you get nDim=',nDim
         call xit
      endif

      return
      end subroutine AssembleElement

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

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the length ds, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 4-node quadrilateral elements

      implicit none

      integer face

      real*8 xLocal,yLocal,ds,dshxi(4,2),sh(4),dXdXi,dXdEta,dYdXi
      real*8 dYdEta,one,coords(2,4),fourth,shape,normal(2,1)
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


      ! Jacobian of the mapping
      !
      if( (face.eq.2).or.(face.eq.4).or.(face.eq.12)
     +    .or.(face.eq.14).or.(face.eq.22).or.(face.eq.24) ) then
         ds = dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
      elseif( (face.eq.1).or.(face.eq.3).or.(face.eq.11)
     +    .or.(face.eq.13).or.(face.eq.21).or.(face.eq.23) ) then
         ds = dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
      else
         write(*,*) 'never should get here'
         call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      !
      if((face.eq.2).or.(face.eq.4).or.(face.eq.12)
     +    .or.(face.eq.14).or.(face.eq.22).or.(face.eq.24)) then
         normal(1,1) = dYdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         normal(2,1) = -dXdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         if(face.eq.4) normal = -normal
         if(face.eq.14) normal = -normal
         if(face.eq.24) normal = -normal
      elseif((face.eq.1).or.(face.eq.3).or.(face.eq.11)
     +    .or.(face.eq.13).or.(face.eq.21).or.(face.eq.23)) then
         normal(1,1) = dYdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         normal(2,1) = -dXdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         if(face.eq.3) normal = -normal
         if(face.eq.13) normal = -normal
         if(face.eq.23) normal = -normal
      else
         write(*,*) 'never should get here'
         call xit
      endif

      return
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

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat,pnewdt)
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
     +  mapJ_inv(2,2),detmapJ,pnewdt
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
         pnewdt = 0.32
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2Da

************************************************************************

      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.  This subroutine works for both 8-node
      !  linear and 20-node quadratic 3D elements.
      !
      implicit none

      integer i,j,k,nNode,ieror,stat

      real*8 dshxi(nNode,3),dsh(nNode,3),coords(3,nNode)
      real*8 mapJ(3,3),mapJ_inv(3,3),detmapJ

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
         write(*,*) 'Problem: detF.lt.zero in mapShape3D'
         call xit
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