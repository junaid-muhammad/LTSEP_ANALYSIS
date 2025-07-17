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
      real Eb,eps

      integer nt,nphi
c      parameter (nt=3,nphi=8)
c      parameter (nt=2,nphi=10)

      real r,dr,w,dw,q2,dq2,tt,dtt,th_pos,th_cm
      real tm,tmn,tmx
      real eps_mod,th_mod,x_mod
      real x_real,dx_real

      integer ipol
      real th_pq
      real*8 phi
      real x1,dx1

      real, Dimension(9) :: t_bin_boundary

      real q2_bin

      integer t_bin, phi_bin

      character*80:: line

      integer j
      real*8 pi
      parameter (pi=3.14159265) 

c   /*--------------------------------------------------*/
c   Constructing bname fo input files
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

c         read(52,*) w,dw,q2,dq2,tt,dtt,th_pos
c         read(52,*) w,dw,q2,dq2,x1,dx1,tt,th_pos
c         write(6,32) w,dw,q2,dq2,tt,dtt,th_pos
 32      format('xsect: ',7f10.4)
         
         th_cm=tt*pi/180.D0
c         print*, " th_cm:", th_cm,tt

         do ip=1,nphi

            phi=(ip-0.5)*2.*pi/nphi
c            phi=(ip-1)*2.*pi/nphi
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

cc ratio is simc/data - see GH logbook, p.55
c            if (r.ne.0) then
c               x_real=x_mod/r
c               dx_real=x_mod*dr/r**2
c            else
c               x_real=0.
c               dx_real=0.
c            endif

c           print*, "ratio check ", r, x_mod, x_real 

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

      integer npol_set
      real Eb,q2_set,w,q2,tm,eps_mod,th_mod, thetacm, x_mod

      real*8 sig
      real*8 tp
      real*8 lambda0_sq,lambdapi_sq,alphapi_t,alphapi_0
      real*8 a,b,c,d,e,f
      real*8 m_pi0sq
      real*8 phicm, phi, pi, mp
      
      real wfactor

      integer i

      real sigT,sigL,sigLT,sigTT
      real SsigT,SsigL,SsigLT,SsigTT,Ssig

      character*80 p_fn

      real par(14)
      real p,pe
      real f_tm,g_W,tav,f_tav
      
      pi = 3.1415926
      mp=0.93827231             !mp

      phicm = phi
      tp = abs(tm)      ! just to make sure it's positive

*     Calculate model thetacm and epsilon at first.
c      print*," enter eps_n_theta  ",npol_set,Eb,w,q2,tm
      call eps_n_theta(npol_set, Eb, w, q2, tm, thetacm, eps_mod)
      print*," return eps_n_theta ",thetacm,eps_mod

c      stop
      
*-------------------------------
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
c         pause
      end do
 9    close(56)
      
c===========================
c hh Vijay
c===========================
*---L
c it0
c      a =  0.25961E+02
c      b = -0.00000E+02  
c      c = -0.15838E+02
c      d =  0.00000E+00
c it1
c      a = 6.4085
c      b = -20.2917 
c      c = -15.8380
c      d = -13.0734
c it2
c      a = -3.8098
c      b = -9.8736 
c      c = -15.8380
c      d = -19.9500
c test
c      a = 6.4085
c      b = -20.2917 
c      c = -15.8380
c      d = -13.0734
c values from parfile
      a = 314.0196
      b = 0.0
      c = -33.0320
      d = 0.0
      
      sigL = ((a+b*log(q2))*exp((c+d*log(q2))*(tp)))   !GH-OK since b=d=0

* compare to SIMC
*   314.0196           14.3953    5
*   -33.0520            2.1530    6
*      sigl = fitpar(5)*exp(fitpar(6)*abs(t_gev))
      SsigL = par(5)*exp(par(6)*tp)

*---T
c it0
c      a =  0.46859E+02
c      b =  0.00000E+02 
c      c = -0.33572E+01
c      d =  0.00000E+00
c it1
c      a = 65.05584
c      b = -48.5525
c      c = 4.3149
c      d = -7.8221
c it2
c      a = 87.77741
c      b = -71.7182
c      c = 13.5436
c      d = -17.2312
c test
c      a = 65.05584
c      b = -48.5525
c      c = 4.3149
c      d = -7.8221
c values from SIMC
      a = 186.85974
      b = 0.0
      c = 3.1376
      d = 0.0
      
      sigT = a+b*log(q2)
     1     +((c+d*log(q2))*(tp-(0.0735+0.028*log(q2_set))*q2_set))
c     1     /((0.0735+0.028*log(q2_set))*q2_set))    ! SIMC doesn't have this division by tav

* compare to SIMC
*   186.85974           2.0564    1
*     3.1376            2.7112    2
*     tav=(0.0735+0.028*log(q2_set))*q2_set
*     sigt = fitpar(1)+fitpar(2)*ftav
      tav=(0.0735+0.028*log(q2_set))*q2_set
      f_tav=(tp-tav)/tav
      Ssigt = par(1)+par(2)*f_tav
      
*---LT
c it0      
c      a =  0.10000E+00
c      b = -0.28000E+02
c      c =  0.35000E+01
c it1
c      a = -100.0000
c      b = -83.8045
c      c = -4.4966
c it2
c      a = 2883.6951
c      b = -35.8699
c      c = -36.2393
c test
c      a = -100.0000
c      b = -83.8045
c      c = -4.4966
c values from SIMC
      a =  0.7270
      b =  0.0323
      c = -14.3963
      
      sigLT = ((a/tp*exp(b/(tp))+c/(tp))*sin(thetacm))

* compare to SIMC
*     0.7270            0.1883    9
*     0.0323            0.0022   10
*   -14.3963            0.5270   11
*      	   siglt=((fitpar(9)/abs(t_gev)*exp(fitpar(10)/abs(t_gev))
*	1	+fitpar(11)/abs(t_gev))*sin(thetacm))

      SsigLT=(par(9)/tp*exp(par(10)/tp)+par(11)/tp)*sin(thetacm)   ! GH-OK
      
*---TT
c     it0
c      a = -0.0008
c      b = 40.7516
c      c = 00.0000
c it1
c      a = -144.7318
c it2
c      a = -150.9774
c test
c     a = -144.7318
c values from SIMC
      a = -0.0008
      b = 40.7516
      c =  0.0000
      
      sigTT = ((a/tp**3*exp(b*tp)+c/tp)*sin(thetacm)**2)  !GH-OK
c      sigTT = sigTT*(tm/(tp+0.139570**2)**2)      !SIMC doesn't have this part

* compare to SIMC
*    -0.0008            0.0001   12
*    40.7516            0.0000   13
*     0.0000            0.0000   14
*      	   sigtt=((fitpar(12)/abs(t_gev)**3*exp(fitpar(13)*abs(t_gev))
*     1           +fitpar(14)/abs(t_gev)**1)*sin(thetacm)**2)

      SsigTT=(par(12)/tp**3*exp(par(13)*tp)+par(14)/tp)*sin(thetacm)**2
            
c===========================

* Since I have nothing better to go on, for now I assume W scales as
* 1/(W^2-mp^2)^2.
c      wfactor=(2.2002**2-mp**2)**2/(w**2-mp**2)**2
      wfactor= 1.0/(w**2-mp**2)**2             !GH-OK

      sig = sigT + eps_mod*sigL+(eps_mod*cos(2.*phicm)*sigTT)
     *     +sqrt(2.*eps_mod*(1.+eps_mod))*(cos(phicm))*sigLT
      sig = sig*wfactor
      sig = sig/2./pi/1.d+06    !dsig/dtdphicm in microbarns/MeV^2/rad  GH-OK
      
      Ssig = SsigT + eps_mod*SsigL+(eps_mod*cos(2.*phicm)*SsigTT)
     *     +sqrt(2.*eps_mod*(1.+eps_mod))*(cos(phicm))*SsigLT
      Ssig = Ssig*wfactor
      Ssig = Ssig/2./pi/1.d+06    !dsig/dtdphicm in microbarns/MeV^2/rad

      write(6,*)' L ', sigL*wfactor,SsigL*wfactor,par(5),par(6),tp
      write(6,*)' T ', sigT*wfactor,SsigT*wfactor
      write(6,*)' LT ', sigLT*wfactor,SsigLT*wfactor,thetacm*180./pi
      write(6,*)' TT ', sigTT*wfactor,SsigTT*wfactor,thetacm*180./pi
      write(6,*)' Sig ', sig,Ssig,wfactor,eps_mod,phicm*180./pi
      
c      if (phi.lt.0.3) then
c         write(6,102) eps_mod,tm,sigL,sigT,sigTT,sigLT, x_mod
c 102     format( ('xmodel: eps=',f5.3,' t=',f5.3,' sigL=',f7.2,' sigT=',f6.2,
c     1        ' sigTT=',f5.2,' sigLT=',f5.2,' x_mod=',f10.6) )
c     endif

c      x_mod = sig
      x_mod = Ssig
      th_mod=thetacm
      
      end
      include 'eps_n_theta.f'
