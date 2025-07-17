      subroutine eps_n_theta(npol,Eb,w,q2,tm,thetacm,eps)

c     To calculate model theta_pq in CM and epsilon. This subroutine is largely
c     based on theta_cm.f function, which in turn is based Jochen's script.

      implicit none

      integer npol
      real Eb,w,q2,tm,thetacm,eps

      REAL s,omega,q,tmin
      REAL p1cm,p3cm,e1cm,e3cm,p1lab

      REAL m2,m3,m4
      REAL m12,m22,m32,m42

      real mp,mp2,mpi,mpi2,mn,mn2,pi
      parameter (mp=.93827231)   !mp
      parameter (mp2=.88035493)  !mp^2
      parameter (mpi=.13956995)   !mpi
      parameter (mpi2=.01947977)  !mpi^2
      parameter (mn=.93956563)   !mn
      parameter (mn2=.88278357)  !mn^2
      parameter (pi=3.14159265)  !pi

      parameter (m3=mpi)
      parameter (m32=mpi2)

      if(npol.gt.0) then
         m2=mp
         m22=mp2
         m4=mn
         m42=mn2
      else
         m2=mn
         m22=mn2
         m4=mp
         m42=mp2
      end if

      s=w*w
c      omega=(s+q2-m22)/(2*m2)
      omega=(s+q2-m22)/(2*m2)
      q=sqrt(q2+omega**2)

*     m12=q2    !error?
      m12=-q2   !mass squared of virtual photon.

      e1cm=(s+m12-m22)/(2*w)
      e3cm=(s+m32-m42)/(2*w)
      p1lab=q
      p1cm=p1lab*m2/w
      p3cm=sqrt(e3cm*e3cm-m32)
      tmin=-((e1cm-e3cm)**2-(p1cm-p3cm)**2)

      if (tm.ge.tmin) then
         thetacm=2*asin(sqrt((tm-tmin)/(4*p1cm*p3cm)))
      else
         thetacm=-1.
         print*, 'eps_n_theta: *** tm=',tm,' <  tmin=',tmin,' ! ***'
      endif

      eps=1.+2.*(q2+omega**2)/(4.*Eb*(Eb-omega)-q2)
      eps=1./eps

c      write(*,'(a13,7(F8.5,1x))')
c     *     'eps_n_theta: ',w,q2,t,tmin,thetacm,eps,omega

      end
