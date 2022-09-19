************************************************************************
! FeFpFg UMAT for modeling plating/stripping of lithium
! 
! Sooraj Narayan and Lallit Anand
! January 2019. Implemented in Abaqus 6.17
************************************************************************
!     State Variables
!     --------------------------------------------------------------
!      statev(1:9) = Fp     ---------- plastic distortion
!      statev(10)  = nuP    ---------- equiv. plastic shear strain rate
!      statev(11)  = S      ---------- isotropic def. resistance in shear
!      statev(12)  = Y      ---------- isotropic def. resistance in tension
!      statev(13)  = gBarP  ---------- equiv. plastic shear strain
!      statev(14)  = eBarP  ---------- equiv. plastic tensile strain
!      statev(15)  = c      ---------- concentration
!
!     Material Properties Vector
!     --------------------------------------------------------------
!      Eyoung  = props(1)  ! Elastic modulus
!      poisson = props(2)  ! Poisson ratio
!      Apre    = props(3)  ! Pre-exponential factor
!      Qact    = props(4)  ! Activation energy in kJ/mol
!      T0      = props(5)  ! Temperature in deg K
!      mRate   = props(6)  ! Rate sensitivity
!      Y0      = props(7)  ! Initial resistance in tension
!      H0      = props(8)  ! Hardening modulus  in tension
!      Ysat    = props(9)  ! Saturation resistance in tension
!      ahard   = props(10) ! Exponent in hardening relation
!      nhard   = props(11) ! Exponent in hardening relation
!      Omega   = props(12)  ! Molar volume of Li
!      alpha1  = props(13) ! Growth factor
!      alpha2  = props(14) ! Growth factor
!      alpha3  = props(15) ! Growth factor
!      cDot    = props(16) ! Rate of change of Li concentration
!**********************************************************************

      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

      include 'aba_param.inc'

      dimension stress(ntens),statev(nstatv),
     + ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     + stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     + props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      character*80 cmname,file1
      character*256 jobName,outDir,fileName
      
      integer i,j,k,l,iterError,lenJobName,lenOutDir

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),theta_t,theta_tau,T_tau(3,3)
      real*8 MatTan(3,3,3,3),S_t,S_tau,Fp_t(3,3),Fp_tau(3,3),gBarP_t
      real*8 gBarP_tau,nuP_t,nuP_tau,gBarLmt,plasticwork
      real*8 Y_t,Y_tau,eBarP_t,eBarP_tau    
      real*8 c_t, c_tau
      
      real*8 Fg_t(3,3), Fg_tau(3,3)
      
      real*8 detF_tau, F_per(3,3),dE, detF_per
      real*8 TK_tau(3,3), TK_per(3,3),T_per(3,3)
      real*8 perM(3,3)
      real*8 TK_per1(3,3),T_per1(3,3)
      
      real*8 Eyoung,poisson,Apre,Qact,T0,mRate
      real*8 Y0,S0,H0,Ysat,Ssat,ahard,nhard
      real*8 Omega,alpha1,alpha2,alpha3,cdot        
      real*8 Gshear,Kbulk,Lambda,check,eBarLmt
      
      
      ! Parameters
      !
      real*8 zero,one,two,half,three,third,nine,root_three
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0,root_three=dsqrt(3.d0))

c************************ Initialize
      !
      ! Open the debug file
      !
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      fileName = outDir(1:lenOutDir)//'\aaMSGS_'//
     +     jobName(1:lenJobName)//'.dat'
      open(unit=80,file=fileName,status='unknown')

      ! Read material properties
      !  
      Eyoung  = props(1)  ! Elastic modulus
      poisson = props(2)  ! Poisson ratio
      Apre    = props(3)  ! Preexponential factor    
      Qact    = props(4)  ! Activation energy      
      T0      = props(5)  ! Temperature
      mRate   = props(6)  ! Rate sensitivity
      Y0      = props(7)  ! Initial resistance in tension
      H0      = props(8)  ! Hardening modulus  in tension
      Ysat    = props(9)  ! Saturation resistance in tension
      ahard   = props(10)  ! Exponent in hardening relation 
      nhard   = props(11)  ! Exponent in hardening relation      
      Omega   = props(12)  ! Molar volume of Li
      alpha1  = props(13) ! Growth factor
      alpha2  = props(14) ! Growth factor
      alpha3  = props(15) ! Growth factor
      cDot    = props(16) ! Rate of change of Li concentration
      
      !write(80,*)'********** In Main'      
      !write(80,*)'Eyoung  = ',Eyoung
      !write(80,*)'poisson = ',poisson
      !write(80,*)'Apre    = ',Apre
      !write(80,*)'Qact    = ',Qact
      !write(80,*)'T0      = ',T0     
      !Write(80,*)'mrate   = ',mrate
      !write(80,*)'Y0      = ',Y0
      !Write(80,*)'Ysat    = ',Ysat
      !write(80,*)'ahard   = ',ahard
      !write(80,*)'nhard   = ',nhard    
      !Write(80,*)'Omega   = ',Omega
      !write(80,*)'alpha1  = ',alpha1
      !Write(80,*)'alpha2  = ',alpha2
      !write(80,*)'alpha3  = ',alpha3
      !Write(80,*)'cDot    = ',cDot
      
c     	When the umat is called for the first iteration in an 
c	increment, abaqus does not update  the deformation
c	gradient, and it passes in a zero strain increment.
c	in this case do not perform the integration. 
c	compute the elastic jacobian and return.
c
      !
      ! Calculate shear and bulk moduli
      !
      Gshear = Eyoung/(two*(one+poisson))
      Kbulk  = Eyoung/(three*(one-two*poisson))
      
      ddsdde  = zero   
      check=zero
      !
      do  i=1,ntens
           check = check+ dabs(dstran(i))
      enddo
      
      if(check .lt. 1.d-15) then
	    c1 = 2.d0*Gshear
            c2 = Kbulk - c1/3.d0
            
            do   i = 1, ndi
                    ddsdde(i,i)=c1
            end do
            
	    do i = ndi+1, ntens
                    ddsdde(i,i)=c1/2.d0
            enddo
            
	    do i = 1,ndi
              do  j = 1, ndi
                    ddsdde(i,j) = ddsdde(i,j) + c2
              end do
            end do      
        return
      endif

      call onem(Iden)
      
      ! Obtain old and new deformation gradients
      !
      F_t = dfgrd0
      F_tau = dfgrd1
        
      ! time(2) is  the value of the total time at 
      ! the beginning of the increment
        
      if(time(2).eq.zero) then
         Fp_t      = Iden
         nuP_t     = zero
         S_t       = Y0/root_three
         Y_t       = Y0
         gBarP_t   = zero
         eBarP_t   = zero
         c_t       = zero
         Fg_t      = Iden 
      else
         Fp_t(1,1) = statev(1)  ! plastic distortion at time t
         Fp_t(2,2) = statev(2)  ! plastic distortion at time t
         Fp_t(3,3) = statev(3)  ! plastic distortion at time t
         Fp_t(2,3) = statev(4)  ! plastic distortion at time t
         Fp_t(3,2) = statev(5)  ! plastic distortion at time t
         Fp_t(1,3) = statev(6)  ! plastic distortion at time t
         Fp_t(3,1) = statev(7)  ! plastic distortion at time t
         Fp_t(1,2) = statev(8)  ! plastic distortion at time t
         Fp_t(2,1) = statev(9)  ! plastic distortion at time t
         nuP_t     = statev(10) ! equiv. plastic shear strain rate at time t
         S_t       = statev(11) ! isotropic def. resistance in shear at time t
         Y_t       = statev(12) ! isotropic def. resistance in tension at time t
         gBarp_t   = statev(13) ! equv. plas. shear strain at time t 
         eBarp_t   = statev(14) ! equv. plas. tensile strain at time t
         c_t       = statev(15) ! Li concentration at time t
         Fg_t(1,1) = statev(16)  ! growth distortion at time t
         Fg_t(2,2) = statev(17)  ! growth distortion at time t
         Fg_t(3,3) = statev(18)  ! growth distortion at time t
         Fg_t(2,3) = statev(19)  ! growth distortion at time t
         Fg_t(3,2) = statev(20)  ! growth distortion at time t
         Fg_t(1,3) = statev(21)  ! growth distortion at time t
         Fg_t(3,1) = statev(22)  ! growth distortion at time t
         Fg_t(1,2) = statev(23)  ! growth distortion at time t
         Fg_t(2,1) = statev(24)  ! growth distortion at time t
      endif
        

            
      call integ(props,nprops,dtime,coords,cmname,time(2),
     +     F_t,F_tau,Fp_t,Fg_t,nuP_t,S_t,Y_t,gBarP_t,eBarP_t,c_t,
     +     Fp_tau,Fg_tau,nuP_tau,S_tau,Y_tau,gBarP_tau,eBarP_tau,c_tau,
     +     T_tau,plasticwork)
            
      
      ! Update state variables at the end of the increment
      !
      statev(1)  = Fp_tau(1,1)
      statev(2)  = Fp_tau(2,2)
      statev(3)  = Fp_tau(3,3)
      statev(4)  = Fp_tau(2,3)
      statev(5)  = Fp_tau(3,2)
      statev(6)  = Fp_tau(1,3)
      statev(7)  = Fp_tau(3,1)
      statev(8)  = Fp_tau(1,2)
      statev(9)  = Fp_tau(2,1)
      statev(10) = nuP_tau
      statev(11) = S_tau 
      statev(12) = Y_tau       
      statev(13) = gBarP_tau
      statev(14) = eBarP_tau 
      statev(15) = c_tau
      statev(16)  = Fg_tau(1,1)
      statev(17)  = Fg_tau(2,2)
      statev(18)  = Fg_tau(3,3)
      statev(19)  = Fg_tau(2,3)
      statev(20)  = Fg_tau(3,2)
      statev(21)  = Fg_tau(1,3)
      statev(22)  = Fg_tau(3,1)
      statev(23)  = Fg_tau(1,2)
      statev(24)  = Fg_tau(2,1)

      ! Time stepping algorithim based on the constitutive response
      !
      eBarLmt=0.02
      !
      umeror = dabs((eBarP_tau - eBarP_t)/eBarLmt)
      if(umeror.le.half) then
         pnewdt = 1.5d0
      elseif((umeror.gt.half).and.(umeror.le.0.8d0)) then
         pnewdt = 1.25d0
      elseif((umeror.gt.0.8d0).and.(umeror.le.1.25d0)) then
         pnewdt = 0.75d0
      else
         pnewdt = half
      endif

      
      umeror = dabs((c_tau - c_t)/5.d4)
      if(umeror.ge.1.d0) then
         pnewdt = 0.5d0
      endif

      ! Return the Cauchy stress
      !
      if(ntens.eq.6) then
         !
         ! 3D problem
         !
         stress(1) = T_tau(1,1)
         stress(2) = T_tau(2,2)
         stress(3) = T_tau(3,3)
         stress(4) = T_tau(1,2)
         stress(5) = T_tau(1,3)
         stress(6) = T_tau(2,3)
      elseif(ntens.eq.4) then
         !
         ! 2D problem
         !
         stress(1) = T_tau(1,1)
         stress(2) = T_tau(2,2)
         stress(3) = T_tau(3,3)
         stress(4) = T_tau(1,2)
         endif

        ! Return the plastic work dissipated
        !
        rpl = plasticwork         
 

        ! Caculate the Jacobian using
        ! a numerical perturbation scheme
        !
        MatTan = zero
        call mdet(F_tau,detF_tau)
        TK_tau = detF_tau*T_tau         
      
        ! dE is the perturbation in a strain component
        ! perM is the perturbation strain matrix
        ! Fper is the perturbed deformation gradient at
        ! the end of the step
        !
        dE = 1.d-9
        !dE = max(nuP_tau*1.d0*1.d-4,1.d-6)
        perM = zero
      
        do k=1,3
          do l=1,3              
              perM(k,l) = 1.d0
              perM = (dE/two)*(perM + transpose(perM))
              F_per = F_tau + matmul(perM,F_t)
        
          call integ(props,nprops,dtime,coords,cmname,time(2),
     +     F_t,F_per,Fp_t,Fg_t,nuP_t,S_t,Y_t,gBarP_t,eBarP_t,c_t,
     +     Fp_tau,Fg_tau,nuP_tau,S_tau,Y_tau,gBarP_tau,eBarP_tau,c_tau,
     +     T_per,plasticwork)
                 
              call mdet(F_per,detF_per)
              TK_per = detF_per*T_per
              
              
              F_per = F_tau - matmul(perM,F_t)
        
          call integ(props,nprops,dtime,coords,cmname,time(2),
     +     F_t,F_per,Fp_t,Fg_t,nuP_t,S_t,Y_t,gBarP_t,eBarP_t,c_t,
     +     Fp_tau,Fg_tau,nuP_tau,S_tau,Y_tau,gBarP_tau,eBarP_tau,c_tau,
     +     T_per1,plasticwork)
                 
              call mdet(F_per,detF_per)
              TK_per1 = detF_per*T_per1
              
              
              do i=1,3
                do j=1,3
                 MatTan(i,j,k,l) = (TK_per(i,j) -TK_per1(i,j))/(two*dE)
                enddo
              enddo
              perM = zero
              
         enddo
        enddo
        
c        
c        do k=1,3
c          do l=1,3              
c              F_per = F_tau
c              F_per(k,l) = F_tau(k,l) + dE
c        
c          call integ(props,nprops,dtime,coords,cmname,time(2),
c     +     F_t,F_per,Fp_t,Fg_t,nuP_t,S_t,Y_t,gBarP_t,eBarP_t,c_t,
c     +     Fp_tau,Fg_tau,nuP_tau,S_tau,Y_tau,gBarP_tau,eBarP_tau,c_tau,
c     +     T_per,plasticwork)
c                 
c              call mdet(F_per,detF_per)
c              TK_per = detF_per*T_per
c              
c              do i=1,3
c                do j=1,3
c                 MatTan(i,j,k,l) = (TK_per(i,j) -TK_tau(i,j))/dE
c                enddo
c              enddo
c              
c              
c         enddo
c        enddo
c        
     
      MatTan = MatTan / detF_tau
      
      !
      ! Return the stress-deformation jacobian
      !
      if(ntens.eq.6) then
          call jac3D(MatTan,ddsdde)
      elseif(ntens.eq.4) then
          call jac2D(MatTan,ddsdde)
      endif      

      


      return
      end subroutine umat

***********************************************************************

      subroutine integ(props,nprops,dtime,coords,cmname,step_time,
     +     F_t,F_tau,Fp_t,Fg_t,nuP_t,S_t,Y_t,gBarP_t,eBarP_t,c_t,
     +     Fp_tau,Fg_tau,nuP_tau,S_tau,Y_tau,gBarP_tau,eBarP_tau,c_tau,
     +     T_tau,plasticwork)
          

      implicit none
      
      integer i,j,k,l,m,n,nprops,nargs,stat
      parameter(nargs=7)
      
      character(len=*) cmname
      real*8 coords(3),step_time     
      real*8 props(nprops),dtime,F_t(3,3),F_tau(3,3)
      real*8 Fp_t(3,3),Fp_tau(3,3),nuP_t,Fg_t(3,3)
      real*8 gBarP_t,gBarP_tau,T_tau(3,3),MatTan(3,3,3,3),S_tau,S_t
      real*8 nuP_tau,T_ta9(3,3),plasticwork,Iden(3,3)
      real*8 Fp_t_inv(3,3),theta0
      real*8 Fe_tr(3,3),Re_tr(3,3),Ue_tr(3,3),Ee_tr(3,3),trEe_tr
      real*8 Ee0_tr(3,3),Me_tr(3,3),Me0_tr(3,3),pbar,tauBar_tr,Np(3,3)
      real*8 theta_tau,args(nargs),Dp_tau(3,3),Dp_vec(3,3),Dp_eig(3)
      real*8 expdtDP(3,3),Me_tau(3,3),Fe_tau(3,3),Re_tau(3,3),tmp
      real*8 detF_tau
      real*8 Ue_tau(3,3),Ee_tau(3,3),Fp_tau_inv(3,3),Ctilde(3,3,3,3)
      real*8 dSdnu,dSdgBar,factor1,factor2,ratio,tauBar_tau
      real*8 StressTempJac(3,3),PlWrkDefJac(3,3)
      real*8 Y_t,Y_tau,eBarP_t,eBarP_tau,H_t   
      
      
      real*8 c_t, c_tau,detFg_tau 
      real*8 m_g1(3),m_g2(3),m_g3(3),n_vec(3)
      real*8 lambda_g1,lambda_g2,lambda_g3       
      real*8 S_g1(3,3), S_g2(3,3),S_g3(3,3),S_g(3,3) 
      real*8 Fg_tau(3,3)      
      real*8 Fi_tr(3,3),Fi_tr_inv(3,3) 
      real*8 Fi_tau(3,3),Fi_tau_inv(3,3)      
      real*8 Finv_tau(3,3),Finv_t(3,3),Fmap_t(3,3)      
      real*8 detFe_tau,Ee_tau_dev(3,3),trEe_tau, trMe_tr
     
      real*8 Dg_tau(3,3),Dg_vec(3,3),Dg_eig(3),expdtDg(3,3)
      real*8 x_pt,R1,R2,R3,xC1,xC2,xC3,y_pt,yC1,yC2,theta_subtended
      real*8 mSlope,cInt,xSurf
      real*8 mR_g2(3)
      
      real*8 Eyoung,poisson,Apre,Qact,T0,mRate
      real*8 Y0,S0,H0,Ysat,Ssat,ahard,nhard
      real*8 Omega,alpha1,alpha2,alpha3,cdot        
      real*8 Gshear,Kbulk,Lambda,check,eBarLmt
      real*8 Stilde,fac,Hsign
      real*8 nu0    

      real*8 xacc,lowerbound,gl,dgdx,upperbound,gu      
      

      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third,six
      real*8 root_three
      real*8 Rgas
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0,
     +     six=6.d0,root_three=dsqrt(3.d0), 
     +     Rgas=8.314d-3)

      external fng ! subroutine function for calculating nuP_tau
********************************************************************* 
      ! Identity matrix
      !
      call onem(Iden)


      ! Read material properties
      !  
      Eyoung  = props(1)  ! Elastic modulus
      poisson = props(2)  ! Poisson ratio
      Apre    = props(3)  ! Preexponential factor    
      Qact    = props(4)  ! Activation energy      
      T0      = props(5)  ! Temperature
      mRate   = props(6)  ! Rate sensitivity
      Y0      = props(7)  ! Initial resistance in tension
      H0      = props(8)  ! Hardening modulus  in tension
      Ysat    = props(9)  ! Saturation resistance in tension
      ahard   = props(10)  ! Exponent in hardening relation 
      nhard   = props(11)  ! Exponent in hardening relation      
      Omega   = props(12)  ! Molar volume of Li
      alpha1  = props(13) ! Growth factor
      alpha2  = props(14) ! Growth factor
      alpha3  = props(15) ! Growth factor
      cDot    = props(16) ! Rate of change of Li concentration
      
      !write(80,*)'********** In Integ'      
      !write(80,*)'Eyoung  = ',Eyoung
      !write(80,*)'poisson = ',poisson
      !write(80,*)'Apre    = ',Apre
      !write(80,*)'Qact    = ',Qact
      !write(80,*)'T0      = ',T0     
      !Write(80,*)'mrate   = ',mrate
      !write(80,*)'Y0      = ',Y0
      !Write(80,*)'Ysat    = ',Ysat
      !write(80,*)'ahard   = ',ahard
      !write(80,*)'nhard   = ',nhard    
      !Write(80,*)'Omega   = ',Omega
      !write(80,*)'alpha1  = ',alpha1
      !Write(80,*)'alpha2  = ',alpha2
      !write(80,*)'alpha3  = ',alpha3
       !Write(80,*)'cDot    = ',cDot
      
      nu0 = Apre*dexp(-Qact/(Rgas*T0))      
      
      
      ! Calculate elastic moduli
      !
      Gshear = Eyoung/(two*(one+poisson))
      Kbulk  = Eyoung/(three*(one-two*poisson))
      Lambda = Kbulk - two_third*Gshear
      
      ! The plasticity properties are in tension, convert them to shear
      !
      nu0  = root_three*nu0
      S0   = Y0/root_three
      Ssat = Ysat/root_three
      H0   = H0/three
      
      ! The following snippet accounts for the coldspot presence
      if( ( (coords(1).ge.25.d0).and.(coords(1).le.75.d0) ) 
     &     .or.( (coords(1).ge.425.d0).and.(coords(1).le.475.d0) ) )then
          cDot = cDot/10
      endif    
      
c      if( ( (coords(1).ge.25.d0).and.(coords(1).le.75.d0) ) 
c     &     .or.( (coords(1).ge.425.d0).and.(coords(1).le.475.d0) ) )then
c          cDot = cDot/10
c      endif  
 
      ! One cycle model for 5 A current and 3mAh areal capacity
      !if(step_time.lt.2150.d0)then
      !    cDot = cDot
      !elseif((step_time.gt.2150.d0).and.(step_time.lt.2170.d0))then
      !    cDot = -cDot*(step_time-2160)/10
      !else
      !    cDot = -cDot
      !endif
      ! One cycle model for 10 A current and 3mAh areal capacity
      if(step_time.lt.1075.d0)then
          cDot = cDot
      elseif((step_time.gt.1075.d0).and.(step_time.lt.1085.d0))then
          cDot = -cDot*(step_time-1080)/5
      elseif(step_time.lt.2155.d0)then
          cDot = -cDot
      elseif((step_time.gt.2155.d0).and.(step_time.lt.2165.d0))then    
          cDot = cDot*(step_time-2160)/5
      elseif(step_time.lt.3235.d0)then
          cDot = cDot
      elseif((step_time.gt.3235.d0).and.(step_time.lt.3245.d0))then    
          cDot = -cDot*(step_time-3240)/5    
      else
          cDot = -cDot
      endif
      
       
      
      

! Volume Ratio 
      !
      call mdet(F_tau,detF_tau)
      call m3inv(F_tau,Finv_tau)
      call m3inv(F_t,Finv_t)
                

      m_g1 = zero      
      m_g2 = zero     
      m_g3 = zero
      m_g3(3) = 1.d0
      
      n_vec=zero
      n_vec(2) =1.d0 ! normal direction here is the 2-direction(in deformed space) for a flat interface
      
      Fmap_t = matmul(Fg_t,Finv_t) ! Fmap = inverse(Fe_t*Fp_t)
cc      Fmap_t = Finv_t ! Testing
          m_g2 = zero   
          do i=1,3
              do j=1,3
                  m_g2(i) = m_g2(i) + Fmap_t(i,j)*n_vec(j)
              enddo
          enddo
          
          !normalization
          m_g2 = m_g2/ dsqrt(m_g2(1)**two + m_g2(2)**two + m_g2(3)**two) 
          
          ! m_g1 orthogonal to m_g2 in the 1-2 plane
          m_g1(1) = m_g2(2)
          m_g1(2) = -m_g2(1)
      
c      !-extra    
c      m_g1 = zero
c      m_g1(1) = 1.d0      
c      m_g2 = zero
c      m_g2(2) = 1.d0      
c      m_g3 = zero
c      m_g3(3) = 1.d0    
c      !-extra
      
      ! m_g2 = zero
      ! m_g2(2) = 1
           
      do i=1,3
          do j=1,3
              S_g1(i,j) = m_g1(i)*m_g1(j)
              S_g2(i,j) = m_g2(i)*m_g2(j)
              S_g3(i,j) = m_g3(i)*m_g3(j)              
          enddo
      enddo
      
      c_tau = c_t + dtime*cDot
      if(c_tau.lt.0.d0)then
          c_tau = zero
          cDot = -c_t/dtime
      endif
      
      
      detFg_tau = one + Omega*c_tau
      
      !-extra
      lambda_g1 = detFg_tau**alpha1
      lambda_g2 = detFg_tau**alpha2      
      lambda_g3 = detFg_tau**alpha3 
      !-extra
      
      Dg_tau=(Omega*cDot/detFg_tau)*
     +               (alpha1*S_g1+alpha2*S_g2+alpha3*S_g3)
      
         call spectral(dtime*Dg_tau,Dg_eig,Dg_vec)
         expdtDg = zero
         expdtDg(1,1) = dexp(Dg_eig(1))
         expdtDg(2,2) = dexp(Dg_eig(2))
         expdtDg(3,3) = dexp(Dg_eig(3))
         expdtDg = matmul(matmul(Dg_vec,expdtDg),transpose(Dg_vec))
         Fg_tau = matmul(expdtDg,Fg_t)   
      
c      !-extra
c         do i=1,3
c        do j=1,3
c          Fg_tau(i,j) = lambda_g1*S_g1(i,j)+lambda_g2*S_g2(i,j)
c     +                  +lambda_g3*S_g3(i,j)           
c        enddo
c      enddo
c      !-extra   
         
         
      ! Compute the trial elastic deformation gradient
      !
      Fi_tr = matmul(Fp_t,Fg_tau)   
      call m3inv(Fi_tr,Fi_tr_inv)
      Fe_tr = matmul(F_tau,Fi_tr_inv)

      ! Perform kinematical calculations for the trial state
      !
      call skinem(Fe_tr,Re_tr,Ue_tr,Ee_tr)


      ! Compute the trace of the trial strain, and the deviator
      !  of the trial strain
      !
      trEe_tr = Ee_tr(1,1) + Ee_tr(2,2) + Ee_tr(3,3)
      Ee0_tr = Ee_tr - third*trEe_tr*Iden

      ! Compute the elastic volume ratio at the end of the increment
      !
      detFe_tau = detF_tau/detFg_tau
      
      ! Compute the trace of the elastic strain at the
      ! end of the increment
      !
      trEe_tau = dlog(detFe_tau)

      ! Compute the trial Mandel stress and related quantities
      !
      Me_tr = two*Gshear*Ee0_tr + Kbulk*(trEe_tau)*Iden
      !
      trMe_tr = Me_tr(1,1) + Me_tr(2,2) + Me_tr(3,3)
      !
      Me0_tr = Me_tr - (one/three)*trMe_tr*Iden


      ! Compute the trial equiv. shear stress 
      !
      tauBar_tr = dsqrt(one/two)*dsqrt(sum(Me0_tr*Me0_tr))


      ! Compute the direction of plastic flow  
      !(this is the trial, as well as the actual value at tau)
      !
      if(tauBar_tr.gt.zero) then
         Np = Me0_tr/(dsqrt(two)*tauBar_tr)
      else
         Np = Iden
      endif

      
      if (nuP_t .gt. zero) then
          Stilde = Ssat*(nuP_t/nu0)**nhard
      else 
          Stilde = Ssat
      end if
      
      fac    = one- S_t/Stilde
      Hsign = fac/dabs(fac)

      H_t = H0*(dabs(fac)**ahard)*Hsign      
    
      ! Solve the implicit relation for the equiv. plastic
      ! shearing rate at the end of the increment
      !
      if(tauBar_tr.le.zero) then
         nuP_tau = zero
      else
         args(1)   = tauBar_tr
         args(2)   = Gshear
         args(3)   = dtime
         args(4)   = S_t
         args(5)   = mRate
         args(6)   = nu0
         args(7)   = H_t     
         !
         !set the accuracy for calculating nup_tau
         !
	 xacc = 1.0d-6	       
	 !       
         !Set the lower bound
         !
         lowerbound = 1.0d-6
         !             
         call fng(lowerbound,gl,dgdx,args,nargs)
         if (gl < zero) then
             nuP_tau = zero
         else
         !	
         !Set the upper bound for rtsafe and check that 
         !a root exists
         !
         upperbound = 1.d0
         call fng(upperbound,gu,dgdx,args,nargs)
         !
         ! Make sure that the root is bracketed
         !             
         do while (gu > zero) 
            upperbound = 10.d0*upperbound
            call fng(upperbound,gu,dgdx,args,nargs)
         end do             
         !
         ! Initial guess for nup_tau
         !
         nup_tau = 1.d0 ! this is the initial guess for nup_tau
         !
         !write(80,*)'before rtsafe nuP_tau = ',nuP_tau                
         !
         ! Calculate the equiv. plastic shear strain rate                
         !
         stat = 1              
         call rtsafe(fng,lowerbound,upperbound,
     +               nup_tau,xacc,args,nargs,stat)
         !
           if (stat==0) then
             !
             ! Rtsafe was unsuccessful 
             !
             write(80,*)'Rtsafe was unsuccessful'
             write(*,*)'Rtsafe was unsuccessful, cut back the mass scaling'            
             call xit            
            endif
          endif
      endif

         
      !write(80,*)'After rtsafe nuP_tau = ',nuP_tau
       

      ! Compute the plastic stretching at the end of the increment
      !
      Dp_tau = (one/dsqrt(two))*nuP_tau*Np


      ! Compute the equiv. plastic shear and tensile strains strain at the
      ! end of the increment
      !
      gBarP_tau = gBarP_t + dtime*nuP_tau
      eBarP_tau = gBarP_tau/root_three

      ! Compute the plastic distortion at the end of
      ! the increment using the exponential map
      !
      if(nuP_tau.le.zero) then
         Fp_tau = Fp_t
      else
         call spectral(dtime*Dp_tau,Dp_eig,Dp_vec)
         expdtDp = zero
         expdtDp(1,1) = dexp(Dp_eig(1))
         expdtDp(2,2) = dexp(Dp_eig(2))
         expdtDp(3,3) = dexp(Dp_eig(3))
         expdtDp = matmul(matmul(Dp_vec,expdtDp),transpose(Dp_vec))
         Fp_tau = matmul(expdtDp,Fp_t)
      endif

      ! Check to make sure that det(Fp_tau)>0
      !
      call mdet(Fp_tau,tmp)
      if(tmp.le.zero) then
         write(*,*) 'det(Fp_tau).le.zero in INTEG'
         call xit
      endif


      ! Compute the elastic distortion at the end
      ! of the increment
      !
      Fi_tau = matmul(Fp_tau,Fg_tau)            
      call m3inv(Fi_tau,Fi_tau_inv)
      Fe_tau = matmul(F_tau,Fi_tau_inv)

      ! Perform kinematical calculations 
      !
      call skinem(Fe_tau,Re_tau,Ue_tau,Ee_tau)
      
      call devm(Ee_tau,Ee_tau_dev)
      call tracem(Ee_tau,trEe_tau)

      ! Compute the Mandel stress  the end of the increment
      !
      Me_tau = two*Gshear*Ee_tau_dev + Kbulk*(trEe_tau)*Iden      
      
      ! Compute the Cauchy stress at the end of the increment
      !
      T_tau = matmul(Re_tau,matmul(Me_tau,transpose(Re_tau)))/detFe_tau


      ! Compute the deformation resistance at the end of the increment
      !
      S_tau = S_t + H_t*dtime*nuP_tau      
      Y_tau = root_three*S_tau
     

      ! plastic power per unit volume 
      !
      plasticwork = tauBar_tau*nuP_tau


      return
      end subroutine integ

****************************************************************************

      subroutine fng(nuP,g,dg,args,nargs)

      implicit none

      integer nargs

      real*8 args(nargs),g,dg,tauBar_tr,Gshear,S_t,S_tau,dtime,mRate,H0
      real*8 ahard,H_t
      real*8 Ssat,dSdnu,dSdgBar,nu0,nuP,aux,fac1,fac2

      real*8 zero,one,two,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0)


      if(nuP.le.zero) then
         g  = one
         dg = zero
         return
      endif

      
      ! Obtain relevant quantities
      !
      tauBar_tr = args(1)
      Gshear    = args(2)
      dtime     = args(3)
      S_t       = args(4)
      mRate     = args(5)
      nu0       = args(6)
      H_t       = args(7)

      
      
      !write(80,*)'********** In fng'
      !write(80,*)'tauBar_tr = ',Taubar_tr
      !Write(80,*)'Gshear    = ',Gshear
      !write(80,*)'dtime     = ',dtime
      !write(80,*)'S_t       = ',S_t     
      !Write(80,*)'mrate     = ',mrate
      !write(80,*)'nu0       = ',nu0
      !write(80,*)'H_t        = ',H_t     
     
      

      ! Compute the resistance for this nuP 
      S_tau = S_t +dtime*H_t*nuP
      dSdnu = dtime*H_t

      ! Compute the residual
      !
      g = tauBar_tr - dtime*Gshear*nuP - S_tau*((nuP/nu0)**mRate)

      ! Compute the tangent
      !
      if((nuP/nu0).gt.zero) then
         fac1 = dSdnu*((nuP/nu0)**mRate)
         fac2 = S_tau*(mrate/nu0)*((nuP/nu0)**(mRate-one))
         dg = -Gshear*dtime - fac1 -fac2
      else
         dg = -dtime*Gshear
      endif


      return
      end subroutine fng
***************************************************************************
!****************************************************************************
!     The next subroutine solves an implicit equation for its root.
!****************************************************************************
      subroutine rtsafe(funcd,x1,x2,root,xacc,args,nargs,stat)
      
      !
      ! This is a fail-safe subroutine to solve an implicit equation
      !  f(x)=0 for its root. This routine utilizes a combination of 
      !  the bisection and the Newton-Raphson schemes. The hybrid
      !  algorithm takes a bisection step whenever Newton-Raphson
      !  would take the solution out of the bracket (x1 is the current
      !  lower bound on x, and x2 is the current upper bound on x), or
      !  whenever Newton-Raphson is not reducing the size of the
      !  brackets rapidly enough.
      !
      ! The root, returned as root, is refined until its accuracy
      !  is know within |root| < xacc.
      !
      ! funcd is a user supplied subroutine which returns both the 
      !  function valueand the first derivative of the function.
      !
      ! This subroutine is based on the one given in "Numerical
      !  Recipes."
      !
      implicit none
      !
      integer nargs,j,maxit,stat
      !
      parameter(maxit=100)
      !
      real*8 x1,x2,root,xacc,args(nargs),xl,xh,fl,fh,swap,dx,f,df,
     +  dxold,temp
      !
      external funcd


      ! Check if one of the bounds is the root
      !
      call funcd(x1,fl,df,args,nargs)
      call funcd(x2,fh,df,args,nargs)
      !
      if (fl.eq.0.d0) then
         root = x1
         return
      else if (fh.eq.0.d0) then
         root = x2
         return
      end if


      ! Verify that there is a root within the given interval
      !
      if ((fl*fh).ge.0.d0) then
         write(*,*) 'Error in bisection: Root must be bracketed.'
         call xit
      end if


      ! Orient the search so that f(xl)<0.
      !
      if (fl.lt.0.d0) then
         xl = x1
         xh = x2
      else
         xh = x1
         xl = x2
         swap = fl
         fl = fh
         fh = swap
      end if


      ! Initialize the guess for the root, the "step size before
      !  last," and the last step
      !
      root = 0.5d0*(x1 + x2)
      dxold = dabs(x2 - x1)
      dx = dxold
      call funcd(root,f,df,args,nargs)


      ! Loop over allowed iterations
      !
      do j = 1,maxit
         !
         ! Bisect if Newton is out of range or not 
         !  decreasing fast enough.
         !
         if (((((root-xh)*df - f)*((root - xl)*df -f)).ge.0.d0)
     +      .or.(dabs(2.d0*f).gt.dabs(dxold*df))) then
            !
            ! Bisection
            !
            dxold = dx
            dx = 0.5d0*(xh-xl)
            root = xl + dx
            !
            if (xl.eq.root) then
               !
               ! Change in root is negligible
               !
               return
               !
            end if
            !
         else
            !
            ! Newton step is acceptable. Take it.
            !
            dxold = dx
            dx = f/df
            temp = root
            root = root - dx
            !
            if (temp.eq.root) then
               !
               ! Change in root is negligible
               !
               return
               !
            end if
            !
         end if
         !
         ! Convergence criterion
         !
         if (dabs(dx).lt.xacc) then
            !
            ! Iteration has converged
            !
            return
         end if
         !
         ! The one new function evaluation per iteration
         !
	 call funcd(root,f,df,args,nargs)
	 !
	 ! Maintain the bracket on the root
	 !
	 if (f.lt.0.d0) then
	    xl = root
	    fl = f
	 else
	    xh = root
	    fh = f
	 end if
	 !
      end do
      
      
      ! Maximum iterations exceeded
      !
      write(*,*) 'Warning in rtsafe: Max iterations exceeded.'
      stat = 0
      
      
      return
      end subroutine rtsafe

****************************************************************************

****************************************************************************

      subroutine jac2D(SpTanMod,ddsdde)

      implicit none

      real*8 SpTanMod(3,3,3,3),ddsdde(4,4)

      ddsdde(1,1) = SpTanMod(1,1,1,1)
      ddsdde(1,2) = SpTanMod(1,1,2,2)
      ddsdde(1,3) = SpTanMod(1,1,3,3)
      ddsdde(1,4) = SpTanMod(1,1,1,2)

      ddsdde(2,1) = SpTanMod(2,2,1,1)
      ddsdde(2,2) = SpTanMod(2,2,2,2)
      ddsdde(2,3) = SpTanMod(2,2,3,3)
      ddsdde(2,4) = SpTanMod(2,2,1,2)

      ddsdde(3,1) = SpTanMod(3,3,1,1)
      ddsdde(3,2) = SpTanMod(3,3,2,2)
      ddsdde(3,3) = SpTanMod(3,3,3,3)
      ddsdde(3,4) = SpTanMod(3,3,1,2)

      ddsdde(4,1) = SpTanMod(1,2,1,1)
      ddsdde(4,2) = SpTanMod(1,2,2,2)
      ddsdde(4,3) = SpTanMod(1,2,3,3)
      ddsdde(4,4) = SpTanMod(1,2,1,2)

      end subroutine jac2D

***********************************************************************

      subroutine jac3D(SpTanMod,ddsdde)

      implicit none

      real*8 SpTanMod(3,3,3,3),ddsdde(6,6)

      ddsdde(1,1) = SpTanMod(1,1,1,1)
      ddsdde(1,2) = SpTanMod(1,1,2,2)
      ddsdde(1,3) = SpTanMod(1,1,3,3)
      ddsdde(1,4) = SpTanMod(1,1,1,2)
      ddsdde(1,5) = SpTanMod(1,1,1,3)
      ddsdde(1,6) = SpTanmod(1,1,2,3)

      ddsdde(2,1) = SpTanMod(2,2,1,1)
      ddsdde(2,2) = SpTanMod(2,2,2,2)
      ddsdde(2,3) = SpTanMod(2,2,3,3)
      ddsdde(2,4) = SpTanMod(2,2,1,2)
      ddsdde(2,5) = SpTanMod(2,2,1,3)
      ddsdde(2,6) = SpTanmod(2,2,2,3)

      ddsdde(3,1) = SpTanMod(3,3,1,1)
      ddsdde(3,2) = SpTanMod(3,3,2,2)
      ddsdde(3,3) = SpTanMod(3,3,3,3)
      ddsdde(3,4) = SpTanMod(3,3,1,2)
      ddsdde(3,5) = SpTanMod(3,3,1,3)
      ddsdde(3,6) = SpTanmod(3,3,2,3)

      ddsdde(4,1) = SpTanMod(1,2,1,1)
      ddsdde(4,2) = SpTanMod(1,2,2,2)
      ddsdde(4,3) = SpTanMod(1,2,3,3)
      ddsdde(4,4) = SpTanMod(1,2,1,2)
      ddsdde(4,5) = SpTanMod(1,2,1,3)
      ddsdde(4,6) = SpTanmod(1,2,2,3)

      ddsdde(5,1) = SpTanMod(1,3,1,1)
      ddsdde(5,2) = SpTanMod(1,3,2,2)
      ddsdde(5,3) = SpTanMod(1,3,3,3)
      ddsdde(5,4) = SpTanMod(1,3,1,2)
      ddsdde(5,5) = SpTanMod(1,3,1,3)
      ddsdde(5,6) = SpTanmod(1,3,2,3)

      ddsdde(6,1) = SpTanMod(2,3,1,1)
      ddsdde(6,2) = SpTanMod(2,3,2,2)
      ddsdde(6,3) = SpTanMod(2,3,3,3)
      ddsdde(6,4) = SpTanMod(2,3,1,2)
      ddsdde(6,5) = SpTanMod(2,3,1,3)
      ddsdde(6,6) = SpTanmod(2,3,2,3)

      end subroutine jac3D

***********************************************************************
c
c
c  The following are all utility routines used in fortran codes
c
c
c
C**********************************************************************
C	THE NEXT SUBROUTINE CALCULATES VARIOUS KINEMATICAL QUANTITIES 
C	ASSOCIATED WITH THE DEFORMATION GRADIENT
C**********************************************************************
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
          WRITE(*,100)
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
C	EIGENVALUES IN ASCENDING ORDER, AND
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

	WRITE (*,'(/1X,A/)') '50 ITERS IN JACOBI SHOULD NEVER HAPPEN'

	RETURN
	END

C**********************************************************************
	SUBROUTINE EIGSRT(D,V,N,NP)

C	GIVEN THE EIGENVALUES [D] AND EIGENVECTORS [V] AS OUTPUT FROM
C	JACOBI, THIS ROUTINE SORTS THE EIGENVALUES INTO ASCENDING ORDER, 
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
C	THE FOLLOWING SUBROUTINES ARE UTILITY ROUTINES
C**********************************************************************
	SUBROUTINE ZEROV(V,SIZE)

C	THIS SUBROUTINE STORES THE ZERO VECTOR IN A VECTOR V
C	OF SIZE SIZE.
C**********************************************************************

	INTEGER SIZE
	REAL*8 V(0:SIZE-1)

	DO 1 I=0,SIZE
	  V(I) = 0.D0
1	CONTINUE
	
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
	SUBROUTINE ONEM(A)

C	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C	3 BY 3 MATRIX [A]
C**********************************************************************

        REAL*8 A(3,3)
        DATA ZERO/0.D0/
        DATA ONE/1.D0/

	DO 1 I=1,3
	  DO 1 J=1,3
	    IF (I .EQ. J) THEN
              A(I,J) = 1.0
            ELSE
              A(I,J) = 0.0
            ENDIF
1       CONTINUE

	RETURN
	END

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
	SUBROUTINE MPROD4(A,B,C)
 
C	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C 	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

	REAL*8 A(4,4),B(4,4),C(4,4)

	DO 2 I = 1, 4
   	  DO 2 J = 1, 4
	    C(I,J) = 0.D0
	    DO 1 K = 1, 4
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1	    CONTINUE
2	CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE DOTPM(A,B,C)

C	THIS SUBROUTINE CALCULATES THE SCALAR PRODUCT OF TWO
C	3 BY 3 MATRICES [A] AND [B] AND STORES THE RESULT IN THE
C	SCALAR C.
C**********************************************************************

	REAL*8 A(3,3),B(3,3),C

	C = 0.D0
	DO 1 I = 1,3
	  DO 1 J = 1,3
            C = C + A(I,J)*B(I,J)
1	CONTINUE
C
	RETURN
	END

C**********************************************************************
	SUBROUTINE MDET(A,DET)
 
C 	THIS SUBROUTINE CALCULATES THE DETERMINANT
C 	OF A 3 BY 3 MATRIX [A].
C**********************************************************************

	REAL*8  A(3,3), DET

	DET =	  A(1,1)*A(2,2)*A(3,3) 
     +	        + A(1,2)*A(2,3)*A(3,1)
     +	        + A(1,3)*A(2,1)*A(3,2)
     +		- A(3,1)*A(2,2)*A(1,3)
     +		- A(3,2)*A(2,3)*A(1,1)
     +		- A(3,3)*A(2,1)*A(1,2)

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
	  write(*,10)
	  STOP
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
        SUBROUTINE INVAR(A,IA,IIA,IIIA)

C	THIS SUBROUTINE CALCULATES THE PRINCIPAL INVARIANTS 
C	IA, IIA, IIIA OF A TENSOR [A].
C**********************************************************************

        REAL*8 A(3,3), AD(3,3),AD2(3,3), DETA, IA,IIA,IIIA

        DO 1 I=1,3
          DO 1 J=1,3
            AD(I,J) = A(I,J)
1       CONTINUE
        IA = AD(1,1) + AD(2,2) + AD(3,3)

C	CALCULATE THE SQUARE OF [AD]

        CALL MPROD(AD,AD,AD2)
        IIA =0.5D0 * ( IA*IA - ( AD2(1,1) + AD2(2,2) + AD2(3,3) ) )

        CALL  MDET(AD,DETA)
        IIIA = DETA

        RETURN
        END

C**********************************************************************
	SUBROUTINE TRACEM(A,TRA)

C	THIS SUBROUTINE CALCULATES THE TRACE OF A 3 BY 3 MATRIX [A]
C	AND STORES THE RESULT IN THE SCALAR TRA
C**********************************************************************

	REAL*8 A(3,3),TRA

	TRA = A(1,1) + A(2,2) + A(3,3)

	RETURN 
	END

C**********************************************************************
	SUBROUTINE DEVM(A,ADEV)

C	THIS SUBROUTINE CALCULATES THE DEVIATORIC PART OF A
C	3 BY 3 MATRIX [A]
C**********************************************************************

	REAL*8 A(3,3),TRA,ADEV(3,3),IDEN(3,3)

	CALL TRACEM(A,TRA)
	CALL ONEM(IDEN)
	CALL ZEROM(ADEV)

	DO 1 I = 1,3
	  DO 1 J = 1,3
	    ADEV(I,J) = A(I,J) - (1.D0/3.D0)*TRA*IDEN(I,J)
1	CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE EQUIVS(S,SB)

C	THIS SUBROUTINE CALCULATES THE EQUIVALENT SHEAR STRESS SB
C	CORRESPONDING TO A 3 BY 3 STRESS MATRIX [S]
C**********************************************************************

	REAL*8 S(3,3),SDEV(3,3),SDOTS,SB

	SB = 0.D0
	SDOTS = 0.D0

	CALL DEVM(S,SDEV)
	CALL DOTPM(SDEV,SDEV,SDOTS)
	SB = DSQRT(1.5D0*SDOTS)

	RETURN
	END
C **********************************************************************
        SUBROUTINE PRESSURE(A,PRA)
C
C       THIS SUBROUTINE CALCULATES THE MEAN NORMAL PRESSURE
C       OF A 3 BY 3 MATRIX [A]
C       AND STORES THE RESULT IN THE SCALAR PRA
C ----------------------------------------------------------------------
C       VARIABLES
C
        REAL*8 A(3,3),PRA

        PRA = -(1.D0 / 3.D0)*( A(1,1) + A(2,2) + A(3,3) )

        RETURN 
        END
C**********************************************************************
	SUBROUTINE PRTMAT(A,M,N)
C**********************************************************************

	INTEGER M,N
	REAL*8 A(M,N)	  

	DO 10 K=1,M
	  WRITE(80,'(2X,6E12.4,2X)') (A(K,L), L=1,N)
10      CONTINUE

        RETURN
        END

C**********************************************************************
	SUBROUTINE PRTVEC(A,M)
C**********************************************************************

	INTEGER M
	REAL*8 A(M)	  

	WRITE(80,'(2X,6E12.4,2X)') (A(K), K=1,M)

        RETURN
	END
C*************************************************************************	  
c***********************************************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
C
C	Given an NxN matrix [A], with physical dimension NP, this 
C	routine replaces it by the LU decomposition of a row-wise 
C	permutation of itself. [A] and N are input. [A] is output, 
C	arranged in LU form. INDX is an output vector which records
C	the row permutation effected by the partial pivoting; 
C	D is output as +1 or -1 depending on wheter the nuber of
C	row interchanges was even or odd, respectively. This routine
C	is used in combination with LUBKSB to solve linear equations 
C	or invert a matrix.
C
	IMPLICIT REAL*8 (A-H,O-Z)
	
	PARAMETER (NMAX=100,TINY=1.0E-20)
	DIMENSION A(NP,NP),INDX(N),VV(NMAX)
	D=1.
	DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
12	CONTINUE
	DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*DABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16	CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19	CONTINUE
	IF(A(N,N).EQ.0.)A(N,N)=TINY
	RETURN
	END
C**********************************************************
       SUBROUTINE LUBKSB(A,N,NP,INDX,B)
C
C	Solves the set of N linear equations [A]{X} = {B}. 
C	Here [A] is input, not as the matrix [A], but as its LU 
C	decomposition, determined by the routine LUDCMP. INDX
C	is input as the permutation vector returned by LUDCMP. {B}
C	is input as the right-hand side vector {B}, and returns
C	with the solution vector {X}. [A], N, NP, INDX are not 
C	modified by this routine, and can be left in place
C	for succesive calls with different right-hand sides {B}.
C	This routine takes into account that {B} will begin with
C	many zero elements, so it is efficient for use in matrix 
C	inversion.
C
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION A(NP,NP),INDX(N),B(N)
       II=0
       DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12     CONTINUE
       DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14     CONTINUE
       RETURN
       END
c**********************************************************
