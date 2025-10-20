      program calc_xsect

      implicit none

c#################################################################
c  Updated by Muhammad Junaid (mjo147@uregina.ca) - 2025/May/31  #
c#################################################################

c This script computes the experimental cross section using the ratio
c DATA/MC(Yields * CS(MC)

      character*2 prv_it
      common prv_it

c      integer q2_bin
      
c     Adding variable to check physics settings
      character*20 physet
      common /phyblk/ physet

c     Get number of the previous iteration.
      
      if(iargc().ne.2) then
         print*,'*** usage: calc_xsect prv_it physet ***'
         stop
      end if
      call getarg(1,prv_it)

c     Calculate unseparated cross-sections

      if (iargc().ne.2) then
         print *, '** Running script for physics setting **', physet
         stop
      end if

      call getarg(2,physet)

      if (trim(physet) .eq. 'Q3p85_W2p62_t0p21') then
c     Set physics settings for Q2, eps
         call xsect(+1, 3.850, 0.292)
         call xsect(+1, 3.850, 0.779)
      else
         print *, 'Add physics setting details in calc_xsect script'
      end if

      stop
      end

*=======================================================================

      subroutine xsect(npol_set,q2_set,eps_set)

      implicit none

      character*2 prv_it
      common prv_it
      
      integer npol_set
      real q2_set,eps_set

      integer kset,nbin

      character*80 r_fn, kin_fn
      character*90 xunsep_fn

      character*2 pol

      integer it,ip

      integer nt,nphi
c      parameter (nt=3,nphi=8)
c      parameter (nt=2,nphi=10)

      real r,dr,w,dw,q2,dq2,tt,dtt,th_pos,th_cm
      real eps_mod,th_mod,x_mod
      real x_real,dx_real

      integer ipol
      real Eb
      real tm
      real tmn,tmx,eps,th_pq
      real x1,dx1
      real phi

      real, Dimension(9) :: t_bin_boundary

      real q2_bin

      integer t_bin, phi_bin

      character*80:: line

      integer j
      real pi
      parameter (pi=3.14159265) 

c   /*--------------------------------------------------*/
c   Constructing bname for input files
      character*20 physet
      character*50 t_bin_interval_file
      character*50 list_settings_pion
      character*50 Eb_pion_file
      common /phyblk/ physet

      t_bin_interval_file = trim(physet) // '_t_bin_interval'
      list_settings_pion = trim(physet) // '_list_settings_pion'
      Eb_pion_file = trim(physet) // '_Eb_pion'

c   Read the t and phi bins 
      open (unit=22, file= "input/"//t_bin_interval_file, action='read')
         
      read (22,*) q2_bin, t_bin, phi_bin
      nt = t_bin
      nphi = phi_bin 

      read (22, '(A)') line  
      read(line, *) (t_bin_boundary(j), j = 1, t_bin + 1)
 
      close(22)

      ipol=0
      q2=0.
      eps=0.
      tmn=0.
      tmx=0.
      kset=0

c     Open the list of settings for pion
      open(55, file= "input/"//list_settings_pion)
      
      do while(ipol.ne.npol_set.or.q2.ne.q2_set.or.eps.ne.eps_set)
         read(55,*) ipol,q2,eps,th_pq,tmn,tmx,nbin,kset
         write(6,2) ipol,q2,eps,th_pq,tmn,tmx,nbin,kset
 2       format(' list: ',i5,5f10.5,2i5)
      end do
      close(55)
      
      write(6,3)tmn,tmx,kset
 3    format(' tmn, tmx: ',2f10.5,'  for kset=',i5)
      if(tmn.eq.0..or.tmx.eq.0.) 
     *     stop '*** setting is not found in list.settings'

c     open the pion energy file
      open(57, file= 'input/' // trim(Eb_pion_file) // '.dat')     
      
      do while(.true.)
         read(57,*) Eb,q2,eps
         write(*,*) Eb,q2,eps
         if(q2.eq.q2_set.and.eps.eq.eps_set) go to 5
      end do
 5    close(57)
      Eb=Eb/1000.

      if(npol_set.lt.0) then
         pol='mn'
      else
         pol='pi'
      end if

      write(6,4)Eb,q2,eps,pol
 4    format(' xsect: Eb=',f8.5,'   at Q2=',f7.4,
     *     '  eps=',f6.4,'  pol=',a2)

c     construct ratio data file name.

      write(r_fn, 10) trim(physet), pol, nint(q2*100), nint(eps*1000)
 10   format('averages/', a, '_aver_', a2, '_', i3.3, '_', i3.3, '.dat')
      print*,'xsect: r_fn=',r_fn

      open(51,file=r_fn)

c     construct kinematics data file name.

      write(kin_fn, 20) trim(physet), nint(q2*100)
 20   format('averages/', a, '_avek_', i3.3, '.dat')

      print*,'xsect: kin_fn=',kin_fn

      open(52,file=kin_fn)

c     construct output file name.
      write(xunsep_fn,30) trim(physet), pol, nint(q2_set*100), 
     1 nint(eps_set*1000)
 30   format('xsects/', a, '_x_unsep_', a2, '_', i3.3, '_', i3.3, 
     1 '.dat')
      print *, 'xsect: xunsep_fn =', trim(xunsep_fn)

      open(61,file=xunsep_fn,status='replace')

      nbin = t_bin

      do it=1,nbin

c         tm=tmn+(it-0.5)*(tmx-tmn)/nbin
         tm = (t_bin_boundary(it) + t_bin_boundary(it+1)) / 2

c         print *," ----- "
c         print *,"t_bin: ",it,t_bin_boundary(it),t_bin_boundary(it+1),tm 
c         stop

         read(52,*) w,dw,q2,dq2,tt,th_pos
         write(6,32) w,dw,q2,dq2,tt,th_pos
 32      format('xsect: ',7f10.4)
         
         th_cm=tt*pi/180.e0
c         print*, " th_cm:", th_cm,tt

         do ip=1,nphi

            phi=(ip-0.5)*2.*pi/nphi
c            print*," "
c            print *," phi bin:",ip,phi*180/pi,"t_bin: ",it,tm 
  
            read(51,*) r,dr
c            print *, "ratio check: ", r, dr
c            stop

c            print *, "set: ", npol_set, Eb, q2_set, w, q2, tm, phi
c            print *, "bin R dr:   ",ip, it, nphi, nbin, r, dr
c            print *, "W,Q2,bin: ", w, dw, q2, dq2, th_pos

c            print*, q2_set, tm
c            stop

            call xmodel(npol_set,Eb,q2_set,w,q2,tm,phi,
     *           eps_mod,th_mod,x_mod)

cc /*--------------------------------------------------*/
cc angle check

           if (abs(th_mod-th_cm).gt.1.e-4*180./pi) then
              write(6,*)'Angle error ',th_mod*180./pi,th_cm*180./pi
              stop
           endif

c /*--------------------------------------------------*/
c cross-section
c ratio is data/simc - see GH logbook, p.55

             x_real=x_mod*r
             dx_real=x_mod*dr

             if (x_real.eq.0.0) then
                dx_real = -1000
             endif
            
c             print*, "ratio check ", r, x_mod, x_real 

c   /*--------------------------------------------------
c   This is where data output to file happens

c            print*, "xmod,eps,th ", x_mod, eps_mod, th_mod*180./pi 

            write(61,40) x_real,dx_real,x_mod,eps_mod,
     *           th_mod*180./pi,phi*180./pi,tm,w,q2
 40         format(3G15.5,f8.5,2f7.2,4f8.5)

         end do                 !phi


c         stop 
         
c         print*, x_real,dx_real,x_mod,eps_mod,phi*180./pi,tm,w,q2

c         stop

c        Write out kinematics for Henk.
c         if(npol_set.gt.0) write(99,'(5f8.3,2x,2f6.2)')
c     *   w,q2,eps_mod,th_mod*180./pi,tm,eps_set,q2_set
      end do                    !t

      close(51)
      close(52)
      close(61)
      print*,' '

      end

*=======================================================================

      subroutine xmodel(npol_set,Eb,q2_set,w,q2_gev,tm,phi,
     *     eps_mod,th_mod,x_mod)

      implicit none

      character*2 prv_it,pol
      common prv_it
      
      integer npol_set          ! SHMS polarity flag, reserved for 2H analysis
      real q2_set
      real Eb,w,q2_gev,tm,phi       ! input kinematic values
      real eps_mod,th_mod,thetacm  ! output kinematic values
      real x_mod                   ! output model cross section
      
      real sig
      real t_gev
      real phicm, pi, mp, mpig

      integer i

*-------------------------------
c var6iables needed for physics_iterate model
      character*80 p_fn
      real fitpar(16)
      real p,pe
      real wfactor
c      real f_tm,g_W,tav,f_tav
      real sigT,sigL,sigLT,sigTT
      real eps
      
*-------------------------------      
      pi = 3.1415926
      mp = 0.93827231           !mp
      mpig=139.57018/1000.      !mpi

      phicm = phi
      t_gev = abs(tm)      ! just to make sure it's positive

*     Calculate model thetacm and epsilon at first.
c      print*," enter eps_n_theta  ",npol_set,Eb,w,q2_gev,tm
      call eps_n_theta(npol_set, Eb, w, q2_gev, tm, thetacm, eps_mod)
c      print*," return eps_n_theta ",thetacm,eps_mod
   
c      stop

*-------------------------------
c     ACTIVATE THIS SECTION WHEN STARTING ITERATIONS
c     IT NEEDS TO BE EXACTLY THE SAME FUNCTIONAL FORM AS IN physics_iterate.f
      
c     Model fit parameters.
      if(npol_set.lt.0) then
         pol='mn'
      else
         pol='pl'
      end if

      write(p_fn,10) prv_it,pol,nint(100*q2_set)
 10   format('fit_params/iter',a2,'/par.',a2,'_',i3.3)

c      write(p_fn,10) pol,nint(100*q2_set),npol_set,nint(100*q2_set)
c 10   format('parameters/par.',a2,'_',i3.3,'/par.',a2,'_',i3.3)

      if (phi.lt.0.3) then
         print*, 'xmodel: p_fn=',p_fn
      endif

      open(56,file=p_fn)
      do while(.true.)
         read(56,*,end=9) p,pe,i
         fitpar(i)=p
         if (phi.lt.0.3) then
            write(6,101)fitpar(i),pe,i
 101        format(' xmodel: '2f11.4,i4)
         endif
      end do
 9    close(56)
      
c========================================================================
c     Tanja's Fpi-2 parameterization in SIMC physics_iterate.f

c      sigt = (fitpar(1)/Q2_gev)+(fitpar(2)/(Q2_gev**2))

c      sigt = (fitpar(1) / Q2_gev) * (exp(fitpar(2) * Q2_gev**2))
c     1        * (exp(fitpar(3) * abs(t_gev)))

c      sigt = (fitpar(1) / Q2_gev) * (exp(fitpar(2) * Q2_gev**2))
c     1        * (exp(fitpar(3) * abs(t_gev))/abs(t_gev))

      sigt = (fitpar(1) / Q2_gev) * (exp(fitpar(2) * Q2_gev**2))
     1        * (exp(fitpar(3) * abs(t_gev)))

c-------------------------------------------------------------------------

c      sigl = (fitpar(5)*Q2_gev
c     1	      *exp((fitpar(6)-fitpar(7)*log(Q2_gev))*abs(t_gev)))
c     2        /(1+fitpar(8)*Q2_gev+fitpar(9)*(Q2_gev**2))**2

c      sigl = (fitpar(5) + fitpar(6)/Q2_gev) * abs(t_gev) /
c     1        (abs(t_gev) + mpig**2)**2 * exp(fitpar(7) * abs(t_gev)) *
c     2        (Q2_gev / (1+fitpar(8)*Q2_gev+fitpar(9)*(Q2_gev**2))**2)

c      sigl = (fitpar(4) + fitpar(5)/Q2_gev) * (abs(t_gev) /
c     1        (abs(t_gev) + mpig**2)**2)*(exp(fitpar(6) * abs(t_gev))) *
c     2        (Q2_gev / (1+fitpar(7)*Q2_gev+fitpar(8)*(Q2_gev**2))**2)

      sigl = (fitpar(4) + fitpar(5)/Q2_gev) * (abs(t_gev) /
     1        (abs(t_gev) + mpig**2)**2)*(exp(fitpar(6) * abs(t_gev))) *
     2        (Q2_gev / (1+fitpar(7)*Q2_gev+fitpar(8)*(Q2_gev**2))**2)

c--------------------------------------------------------------------------

c      siglt=(exp(fitpar(10)+(fitpar(11)*abs(t_gev)/sqrt(Q2_gev)))
c     1       +fitpar(12)+(fitpar(13)/(Q2_gev**2)))*sin(thetacm)

c      siglt=(((fitpar(10)/(1+Q2_gev)) * exp(fitpar(11) * abs(t_gev))) + 
c     1       (fitpar(12) / (abs(t_gev)**2))) * sin(thetacm)

c      siglt=(((fitpar(10)/(1+(fitpar(11) * Q2_gev))) * 
c     1      (abs(t_gev)/(abs(t_gev)+mpig**2)**2)) + 
c     2       (fitpar(12) / (abs(t_gev)**2))) * sin(thetacm)

c      siglt=(((fitpar(10)/(Q2_gev)) + exp(fitpar(11) * abs(t_gev))) * 
c     1       (fitpar(12) / (abs(t_gev)**2))) * sin(thetacm)

c      siglt=((fitpar(9)/(Q2_gev)) + (exp(fitpar(10) * abs(t_gev))) * 
c     1       (fitpar(11)/(fitpar(12)+abs(t_gev))**2)) * sin(thetacm)

c      siglt=((fitpar(9)/(Q2_gev)) + (exp(fitpar(10) * abs(t_gev))) * 
c     1       (fitpar(11)/(abs(t_gev))**fitpar(12))) * sin(thetacm)

c      siglt=((fitpar(9)/(Q2_gev)) * (exp(fitpar(10) * abs(t_gev)))
c     1       *(abs(t_gev)/(abs(t_gev)+mpig**2)**2))*sin(thetacm)

c      siglt=((fitpar(9)/(Q2_gev)) + (fitpar(10)/abs(t_gev))) 
c     1 * sin(thetacm)

c      siglt=((fitpar(9)/(Q2_gev)) + ((exp(-fitpar(10) * abs(t_gev))) * 
c     1       (fitpar(11)/(fitpar(12)+abs(t_gev))))) * sin(thetacm)

c      siglt=((fitpar(9)/(Q2_gev)) + ((exp(-fitpar(10) * abs(t_gev))) * 
c     1     (fitpar(11)*abs(t_gev)**fitpar(12)))) * sin(thetacm)

      siglt=((fitpar(9)/(Q2_gev)) + (fitpar(10)/(abs(t_gev))) + 
     1       ((exp(fitpar(11) * abs(t_gev))) * 
     2       (fitpar(12)/(abs(t_gev))**2))) * sin(thetacm)

c--------------------------------------------------------------------------

c      sigtt=((fitpar(14)/(Q2_gev**2)) 
c     1       *(abs(t_gev)/(abs(t_gev)+mpig**2)**2))*sin(thetacm)**2

c      sigtt=((fitpar(14)/(Q2_gev**2)) 
c     1       *(abs(t_gev)/(abs(t_gev)+mpig**2)**2) * 
c     2       (exp(fitpar(15) * abs(t_gev))))*sin(thetacm)**2

c      sigtt=(((fitpar(14)/(1+Q2_gev))*exp(fitpar(15) * abs(t_gev)))+ 
c     1       (fitpar(16) / (abs(t_gev)**3))) * sin(thetacm)**2

c      sigtt=((fitpar(14)/(Q2_gev))+
c     2      (abs(t_gev)/(abs(t_gev)+mpig**2)**2) * 
c     1       (fitpar(15) / (abs(t_gev)**3))) * sin(thetacm)**2

c      sigtt=((fitpar(14)/(Q2_gev)) + 
c     1       (fitpar(15) / (abs(t_gev)**3))) * sin(thetacm)**2

c      sigtt=((fitpar(13)/(Q2_gev)) +
c     1       (fitpar(14) / (fitpar(15)+(abs(t_gev)))**3)) 
c     2      * sin(thetacm)**2

c      sigtt=((fitpar(13)/(Q2_gev)) + (fitpar(14)/(abs(t_gev)**2)) + 
c     1       ((exp(fitpar(15) * abs(t_gev))) * 
c     2       (fitpar(16)/(abs(t_gev))**3))) * sin(thetacm)**2

c      sigtt=((fitpar(13)/(Q2_gev)) + (exp(fitpar(14) * abs(t_gev))) *
c     1       (fitpar(15) / (fitpar(16)+(abs(t_gev)))**3)) 
c     2      * sin(thetacm)**2

c      sigtt=((fitpar(13)/(Q2_gev) + (fitpar(14)/(abs(t_gev)))) * 
c     1       ((exp(fitpar(15) * abs(t_gev))) * 
c     2       (fitpar(16)/(abs(t_gev))**3))) * sin(thetacm)**2

      sigtt=((fitpar(13)/(Q2_gev)) +
     1     (fitpar(14)/(abs(t_gev))) + 
     2       ((exp(fitpar(15) * abs(t_gev))) * 
     3       (fitpar(16)/(abs(t_gev))**3))) * sin(thetacm)**2

c---------------------------------------------------------------------------
      sigT=sigt
      sigL=sigl
      sigLT=siglt
      sigTT=sigtt

c We assume W scales as 1/(W^2-mp^2)^2.
      wfactor= 1.0/(w**2-mp**2)**2
c
      sig = sigT + eps_mod*sigL+(eps_mod*cos(2.*phicm)*sigTT)
     *     +sqrt(2.*eps_mod*(1.+eps_mod))*(cos(phicm))*sigLT
      sig = sig*wfactor
      sig = sig/2./pi/1.d+06    !dsig/dtdphicm in microbarns/MeV^2/rad 
c
c      write(6,*)' L ', sigL,sigL*wfactor,fitpar(5),fitpar(6),fitpar(7)
c      write(6,*)' T ', sigT,sigT*wfactor
c      write(6,*)' LT ', sigLT,sigLT*wfactor,thetacm*180./pi
c      write(6,*)' TT ', sigTT,sigTT*wfactor,thetacm*180./pi
c      write(6,*)' sig ', sig,wfactor,eps_mod,phicm*180./pi
c      
      x_mod = sig
      th_mod=thetacm
      
      if (phi.lt.0.3) then
         write(6,102) eps_mod,tm,sigL,sigT,sigTT,sigLT,x_mod
 102     format( ('xmodel: eps=',f5.3,' t=',f5.3,' sigL=',f7.2,' sigT=',
     1    f6.2,' sigTT=',f5.2,' sigLT=',f5.2,' x_mod=',e10.4) )
      endif

      return 
      end
      
      include 'eps_n_theta.f'
