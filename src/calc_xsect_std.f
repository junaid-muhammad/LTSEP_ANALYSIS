      program calc_xsect

      implicit none

c#################################################################
c  Updated by Muhammad Junaid (mjo147@uregina.ca) - 2025/May/31  #
c#################################################################

c This script computes the experimental cross section using the ratio
c DATA/MC(Yields * CS(MC)

c      character*2 prv_it
c      common prv_it

c      integer q2_bin
      
c     Get number of the previous iteration.
      
c      if(iargc().ne.1) then
c         print*,'*** usage: calc_xsect prv_it ***'
c         stop
c      end if
c      call getarg(1,prv_it)

c     Calculate unseparated cross-sections

c     Adding variable to check physics settings
      character*20 physet
      common /phyblk/ physet

      if (iargc().ne.1) then
         print *, '** Running script for physics setting **', physet
         stop
      end if

      call getarg(1,physet)

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
30    format('xsects/', a, '_x_unsep_', a2, '_', i3.3, '_', i3.3, 
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
            print*," "
            print *," phi bin:",ip,phi*180/pi,"t_bin: ",it,tm 
  
            read(51,*) r,dr
c            print *, "ratio check: ", r, dr
c            stop

c             print *, "set: ", npol_set, Eb, q2_set, w, q2, tm, phi
c             print *, "bin R dr:   ",ip, it, nphi, nbin, r, dr
c             print *, "W,Q2,bin: ", w, dw, q2, dq2, th_pos

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
            
             print*, "ratio check ", r, x_mod, x_real 

c   /*--------------------------------------------------
c   This is where data output to file happens

            print*, "xmod,eps,th ", x_mod, eps_mod, th_mod*180./pi 

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

      subroutine xmodel(npol_set,Eb,q2_set,w,q2,tm,phi,
     *     eps_mod,th_mod,x_mod)

      implicit none

      integer npol_set          ! SHMS polarity flag, reserved for 2H analysis
      real q2_set
      real Eb,w,q2,tm,phi   ! input kinematic values
      real eps_mod,th_mod,thetacm  ! output kinematic values
      real x_mod                   ! output model cross section
      
      real sig
      real tp
      real phicm, pi, mp

      integer i

*-------------------------------
c variables needed for physics_iterate model
      character*80 p_fn
      real par(14)
      real p,pe
      real wfactor
      real f_tm,g_W,tav,f_tav
      real sigT,sigL,sigLT,sigTT

*-------------------------------
c variables needed for sig_param_2021 model
c remove this section when starting iterations
      real thcm,t,wsq,eps
      real sigttv,sigltv,sigexcl3,sigexcl3p,sigexcl3m
      common/excmn/sigl,sigt,sigttv,sigltv
      real pp(17)/
     >     1.60077,
     >    -0.01523,
     >    37.08142,
     >    -4.11060,
     >    23.26192,
     >     0.00983,
     >     0.87073,
     >    -5.77115,
     >  -271.08678,
     >     0.13766,
     >    -0.00855,
     >     0.27885,
     >    -1.13212,
     >    -1.50415,
     >    -6.34766,
     >     0.55769,
     >    -0.01709/
      real pm(17)/
     >     1.75169,
     >     0.11144,
     >    47.35877,
     >    -4.69434,
     >     1.60552,
     >     0.00800,
     >     0.44194,
     >    -2.29188,
     >   -41.67194,
     >     0.69475,
     >     0.02527,
     >    -0.50178,
     >    -1.22825,
     >    -1.16878,
     >     5.75825,
     >    -1.00355,
     >     0.05055/

*-------------------------------      
      pi = 3.1415926
      mp = 0.93827231             !mp

      phicm = phi
      tp = abs(tm)      ! just to make sure it's positive

*     Calculate model thetacm and epsilon at first.
c      print*," enter eps_n_theta  ",npol_set,Eb,w,q2,tm
      call eps_n_theta(npol_set, Eb, w, q2, tm, thetacm, eps_mod)
      print*," return eps_n_theta ",thetacm,eps_mod

c      stop

*-------------------------------
c     ACTIVATE THIS SECTION WHEN STARTING ITERATIONS
c     IT NEEDS TO BE EXACTLY THE SAME FUNCTIONAL FORM AS IN physics_iterate.f
      goto 200     ! temporary skip around this section
      
*     Model fit parameters.
c      write(fn,10) prv_it,pol,nint(100*q2_set)
c 10   format('fit_params/it',a2,'/par.',a2,'_',i3.3)
*      write(fn,10) pol,nint(100*q2_set),pol,nint(100*q2_set)
* 10   format('parameters/par.',a2,'_',i3.3,'/par.',a2,'_',i3.3)
*      if (phi.lt.0.3) then
*         print*, 'xmodel: fn=',fn
*     endif

      p_fn='parameters/par.pl'

      open(56,file=p_fn)
      do while(.true.)
         read(56,*,end=9) p,pe,i
         par(i)=p
         if (phi.lt.0.3) then
            write(6,101)par(i),pe,i
 101        format(' xmodel: '2f11.4,i4)
         endif
      end do
 9    close(56)
      
c===========================
c      Vijay's pi+n SIMC parameterization in physics_iterate.f
c      
      sigL = par(5)*exp(par(6)*tp)

      tav=(0.0735+0.028*log(q2_set))*q2_set
      f_tav=(tp-tav)/tav
      sigt = par(1)+par(2)*f_tav
      
      sigLT=(par(9)/tp*exp(par(10)/tp)+par(11)/tp)*sin(thetacm) 
      
      sigTT=(par(12)/tp**3*exp(par(13)*tp)+par(14)/tp)*sin(thetacm)**2
            
c* Since I have nothing better to go on, for now I assume W scales as
c* 1/(W^2-mp^2)^2.
cc      wfactor=(2.2002**2-mp**2)**2/(w**2-mp**2)**2
      wfactor= 1.0/(w**2-mp**2)**2
c
      sig = sigT + eps_mod*sigL+(eps_mod*cos(2.*phicm)*sigTT)
     *     +sqrt(2.*eps_mod*(1.+eps_mod))*(cos(phicm))*sigLT
      sig = sig*wfactor
      sig = sig/2./pi/1.d+06    !dsig/dtdphicm in microbarns/MeV^2/rad 
c
      write(6,*)' L ', sigL*wfactor,par(5),par(6),tp
      write(6,*)' T ', sigT*wfactor
      write(6,*)' LT ', sigLT*wfactor,thetacm*180./pi
      write(6,*)' TT ', sigTT*wfactor,thetacm*180./pi
      write(6,*)' sig ', sig,wfactor,eps_mod,phicm
      phicm = phicm*180./pi
c      if (phicm .lt. 0.0) phicm = phicm + 360.0
c      
c      if (phi.lt.0.3) then
c         write(6,102) eps_mod,tm,sigL,sigT,sigTT,sigLT, x_mod
c 102     format( ('xmodel: eps=',f5.3,' t=',f5.3,' sigL=',f7.2,' sigT=',f6.2,
c     1        ' sigTT=',f5.2,' sigLT=',f5.2,' x_mod=',f10.6) )
c     endif
      goto 300

cc===========================
c remove the next section when starting iterations
c      
c SIMC MODEL FROM physics_pion.f FOR A QUICK TEST
c      
c      real*8 function sig_param_2021(thcm,phicm,t,q2,wsq,eps,
c     >      which_pion,doing_pizero)
! April 2021 fit to exclusive pi+ and pi- data from
! fpi1, fpi2, CT, pt-SIDIS, CSV-SIDIS, and KLT
! q2, t, and wsq should be in gev**2
! thcm and phicm should be in radians
! Peter Bosted, from output of ~bosted/ptc/ptb.f
 200  continue
      t=tp
      thcm=th_mod
      wsq=w*w
      eps=eps_mod
      
c      if (which_pion.eq.1.or.which_pion.eq.11.or.which_pion.eq.3) then
c         call exclfit(t,thcm,phicm,q2,
c     >           wsq,eps,sigl,sigt,sigttv,sigltv,sigexcl3,pm,1.0)
c	   else			! pi+
              call exclfit(t,thcm,phicm,q2,
     >          wsq,eps,sigl,sigt,sigttv,sigltv,sigexcl3,pp,1.0e0)
c      endif
c      sig_param_2021 = sigexcl3
       sig = sigexcl3

c      write(6,*)t,thcm,phicm,q2,wsq,eps
c      write(6,*)sigl,sigt,sigttv,sigltv
c      write(6,*)sigexcl3,pp

c===========================
 300  continue

      x_mod = sig
      th_mod=thetacm
      
      return 
      end

      subroutine exclfit(t,thetacm,phicm,q2_gev,
     >  s_gev,eps,sigl,sigt,sigtt,siglt,sig,p,fpifact)

      implicit none
      real p(17),t,thetacm,phicm,q2_gev,s_gev,eps,sigl
      real sigt,sigtt,siglt,sig,mtar_gev/0.938/
      real fpi,q2fpi2,sig219,fpifact

       fpi = fpifact / 
     >  (1.0 + p(1)*q2_gev + p(2)*q2_gev**2)

       q2fpi2 = q2_gev * fpi**2

       sigL = (p(3) + p(15)/q2_geV) * abs(t) / 
     >   (abs(t) + 0.02)**2 * q2fpi2 * 
     >   exp(p(4) * abs(t))
       sigL = sigL / (s_gev**p(11) + sqrt(s_gev)**p(17))

       sigT = p(5) / q2_gev * exp(p(6) * q2_gev**2)
       sigT = sigT / (s_gev**p(12) + sqrt(s_gev)**p(16))
       sigT = sigT * exp(p(14) * abs(t))
 
      sigLT=(p(7) / (1.0 + p(10) * q2_gev)) *
     >   exp(p(8) * abs(t)) * sin(thetacm)
       sigLT = sigLT / s_gev**p(13)

       sigTT=(p(9) / (1. + 1.0 * q2_gev)) * 
     >   exp(-7.0 * abs(t)) * sin(thetacm)**2

        sig219 = (sigt + 
     >   eps * sigl + 
     >   eps * cos(2.0*phicm) * sigtt + 
     >   sqrt(2.0 * eps * (1.0+eps)) * 
     >   cos(phicm) * siglt) / 1.0e0

* GMH: Brauel scaled all his data to W=2.19 GeV, so rescale by 1/(W**2-mp**2)**2
* HPB: factor 15.333 therefore is value of (W**2-mp**2)**2 at W=2.19
c BUT, in param_3000, everything got scaled to W=1.96

        sig = sig219 *8.539 / (s_gev - mtar_gev**2)**2
        sig = sig / 2.0 / 3.1415928 / 1.0d+06  !dsig/dtdphicm in microbarns/MeV**2/rad

      write(6,*)' L ', sigL,t
      write(6,*)' T ', sigT
      write(6,*)' LT ', sigLT,thetacm
      write(6,*)' TT ', sigTT,thetacm
      write(6,*)' sig ', sig,eps,phicm

      return
      end
      
      include 'eps_n_theta.f'
