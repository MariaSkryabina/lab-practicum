c*NVVMIE2.for ******************************************
c*NVVMIE2 * Light Scattering by Spherical Particles ****
c                                                      *
c  Calculations of extinction, scattering, absorption, *
c  etc. efficiency factors and scattering matrix       *
c  for homogeneous spheres (Mie theory)                *
c......................................................*
c  Input data:          filename:  nvvmie2.dat         *
c   ri = n+k*i: complex index of refraction            *
c           x1: minumum size                           *
c           dx: step over size                         *
c           x2: maximum size                           *
c       dtheta: step over scattering angle             *
c......................................................*
c  Output data:         filename:  nvvmie2.out         *
c   Qext : extinction factor                           *
c   Qsca : scattering factor                           *
c   Qabs : absorption factor                           *
c   Qbk  : backscattering factor                       *
c   Qpr  : radiation pressure factor                   *
c  albedo: albedo                                      *
c   g    : asymmetry factor                            *
c P1,... : elements of scattering matrix               *
c            for Stokes vector (I_r,I_l,U,V)           *
c    I   : intensity of scattered light                *
c            element (1,1) for Stokes vector (I,Q,U,V) *
c  I/I0  : normalized intensity                        *
c  Pol   : polarization degree                         *
c......................................................*
c NB! In order to treat very large particles,          *
c     one needs to enlarge the parameter NTERMS.       *
c     In order to consider more scattering angles      *
c     one needs to enlarge the parameter NTHETA.       *
c......................................................*
c created by N.V. Voshchinnikov                        *
c (c) 2002 Astronomical Institute, St. Petersburg Univ.*
c*******************************************************
c
      parameter(ntheta=361)
      implicit real*8(a-h,o-z)
      complex*16, allocatable:: r_ri_array(:)
      complex*16 :: ri
      real(8), allocatable::r_lambda(:), r_Qpr_array(:)
      DIMENSION
     *          theta(ntheta), P1(ntheta), P2(ntheta),
     *          P3(ntheta), P4(ntheta)
    1 format(3X,'SPHERES: homogeneous',/,
     *       3x,'THEORY:  exact')
    2 format(3x,'m = ',2F14.10,'*i')
    3 format(1X,91('.'))
    4 format(6X,'x',9X,'Qext',7X,'Qsca',7X,'Qabs',9X,
     *       'Qbk',11x,'Qpr',6X,'albedo',7x,'g')
    5 format(1X,91('-'))
    6 format(1X,F10.2,3(F11.6),f15.6,3(f11.5))
C    7 format(2d14.10)
    7 format(d13.3,2d11.4)
    8 format(d16.6)
      common /par/ pi180, dtheta, k
      common /theta/ theta
c* input
      print *,'start nvvmie2'

      open(unit=1, file='silicate.RI', status='old')
      i_row_number = 837
      allocate(r_lambda(0:i_row_number-1) )
      allocate(r_ri_array(0:i_row_number-1) )
      do i_counter=0,i_row_number-1
         read(1,7) r_lambda(i_counter), r_ri_array(i_counter)
      enddo
      close(1)   
      
      
      write(*,*) r_ri_array(0)
      
C      open (unit=05,file='nvvmie2.dat',status='old',
C     *      access='sequential')
C      open (unit=07,file='nvvmie2.out',status='unknown',
C     *      access='append')
      open (unit=07,file='nvvmie2.out',status='replace')
c     .........................................................................
c     .   set program constants
c     .........................................................................
       pi = 4d0 * datan(1d0)
       dpi = 2d0 * pi
       PI180 = DPI / 3.6D2                                                  
C     read (5,7) ri
C      read (5,8) x1, dx, x2, dtheta
C      nx = (x2 - x1) / dx + 1

      theta(1) = 0d0
      dtheta = 0.5d0
      DO I = 2, 100000
       if(i.gt.ntheta) then
         write(*,*) 'ntheta, i=', ntheta, i
         pause
         stop
       end if
      theta(I) =  theta(i-1) + dtheta
      if(dabs(theta(i)-180d0).lt.1d-6) go to 654
      end do
  654 continue
      k = i
      IF(K.GT.ntheta) STOP 181

c      dh = dtheta * pi180

C      write (7,1)
C      write (7,2) ri
C      write (7,3)
C      write (7,4)
C      write (7,5)
c*  calculations
C      do i = 1, nx
C        x = x1 + dx * (i - 1)
      r_a_max = 1.0D0
      r_a_step = 0.05D0
      r_a_min = r_a_step
      i_a_number = int((r_a_max - r_a_min) / r_a_step) + 1
      allocate(r_Qpr_array(0:i_a_number-1) )
    
      write(*,*) "Welcome to Diluc's tavern!"
      do i_lambda_index = 0, i_row_number - 1
       
       do i_a_index = 0,i_a_number - 1
        r_a = r_a_min + i_a_index*r_a_step
        x = 2.0D0*pi*r_a/r_lambda(i_lambda_index)
        ri = r_ri_array(i_lambda_index)
        call shexqi(ri, x, Qext, Qsca, Qabs, Qbk, Qpr,
     *               albedo, g, P1, P2, p3, p4, ier)
c* output
       ! if (ier.eq.0) then
           r_Qpr_array(i_a_index) = Qpr
           
        
        
C          write (7,6) x, Qext, Qsca, Qabs, Qbk, Qpr,
C     *                albedo, g
C          write (*,*) x, Qext, Qabs

   24 FORMAT(1X,'Theta',5X,' P1',9X,' P2',9x,' P3',9x,' P4',10x,
     *       'I',10x,'I/I0',8X,'Pol')
   40 FORMAT(1X,91('-'))
*  440 FORMAT(1X,78('-'))
   41 FORMAT(1X,91('.'))
C      write(7, 41)
C      write(7, 24)
C      write(7, 41)
      DO 12 J=1,K 
      fi=(P1(J)+P2(J))/2d0
      if(j.eq.1) pfi=fi
      ffi=fi/pfi
      P=(P1(J)-P2(J))/(P1(J)+P2(J))*1D2
c  250 FORMAT(1X,f5.1,1p5d12.3,0pF9.3,2X,F9.3)
  250 FORMAT(1X,f5.1,1p6d12.3,2X,0pF9.3)
C     write(7,250) Theta(j), P1(J), P2(J), P3(J), P4(J), fi, ffi, P
 12   CONTINUE                                  
*      PRINT 440
C      write(7, 40)
   11 CONTINUE                                  

        !end if
c        if(ier.eq.1) then
c          write(7,*) 'ier=1'
c        end if
c        if(ier.eq.2) then
c          write(7,*) 'ier=2'
          !pause
          !stop
c        end if
       end do 
       write(7,*) r_Qpr_array
      end do
      write(*,*) "You are drunk!"
c*
C*      write (7,5)
      close (7)
      stop
      end
c-------------------------------------------------------
c **********   shexqi -  Spheres: homogeneous
c                        Theory: exact
c                        Results: efficiency factors + matrix
c-------------------------------------------------------
      subroutine shexqi(ri, x, Qext, Qsca, Qabs, Qbk,
     *                   Qpr, albedo, g, p1, p2, p3, p4, ier)
      parameter(nterms=5000, ntheta=361)
      implicit real*8(a-h,o-q,t-z), complex*16(r-s)
      dimension ru(2*nterms)
      dimension fact(0:1)
      DIMENSION
     *          ss1(ntheta), ss2(ntheta),
     *          P1(ntheta), P2(ntheta),
     *          P3(ntheta), P4(ntheta),
     *          Pii(nterms,ntheta), tau(nterms,ntheta)
      data factor / 1d250 /
      data eps    / 1d-15 /
      data fact / 1d0, 1d250 /
      common /par/ pi180, dd, k
      common /theta/ theta(ntheta)
c* efficiency factors
      ier  = 0
      Qext = 0d0
      Qsca = 0d0
      Qabs = 0d0
      Qbk  = 0d0
      Qpr  = 0d0
      albedo = 0d0
      g    = 0d0
c* null argument
      if(x.le.1d-6) then
       ier=1
       return
      end if
c*
      pi = 4d0 * datan(1d0)
      ax = 1d0 / x
      b = 2D0 * ax**2
      ss = (0d0,0d0)
      s3 = (0d0,-1d0)
      an = 3d0
c* choice of number for subroutine aa [Loskutov (1971)]
       y = cdsqrt(RI*dconjg(ri))*x
       num = 1.25 * y + 15.5
       if(y.lt.1d0) go to 11
       if(y.gt.100d0.and.y.lt.50000d0) go to 12
       if(y.ge.50000d0) go to 13
       go to 14
   11  num = 7.5 * y + 9.0
       go to 14
   12  num = 1.0625 * y + 28.5
       go to 14
   13  num=1.005*y+50.5
   14  continue
       if(num.gt.2*nterms) then
         ier = 2
         write(*,*) '2*nterms, num=', 2*nterms, num
         return
       end if
c* logarithmic derivative to Bessel function
c* (complex argument)
      call aa(ax,ri,num,ru)
c* Bessel functions (first terms)
      ass = dsqrt(pi / 2d0 * ax)
      w1 =  2d0 / pi * ax
        Si = dsin(x)/x
        Co = dcos(x)/x
c n=0
          besJ0 = Si / ass
          besY0 = - Co / ass
          iu0 = 0
c n=1
          besJ1 = (Si * ax - Co) / ass
          besY1 = (- Co * ax - Si) / ass
          iu1 = 0
          iu2 = 0
c* Mie coefficients (first term)
      s = ru(1) / ri + ax
      s1 = s * besJ1 - besJ0
      s2 = s * besY1 - besY0
      ra0 = s1 / (s1 - s3 * s2)
      s = ru(1) * ri + ax
      s1 = s * besJ1 - besJ0
      s2 = s * besY1 - besY0
      rb0 = s1 / (s1 - s3 * s2)
c* efficiency factors (first term)
      r = -1.5d0*(ra0-rb0)
      Qext = an * (ra0 + rb0)
      Qsca = an * (ra0 * dconjg(ra0) + rb0 *
     *       dconjg(rb0))
c* amplitude functions (first term)
      sSs=(0.0D0,0.5D0)
      
      write(*,*) 'take a drink for Childe!'
      DO N = 1, K
        PIi(1,n) = 1D0
        C=dcos(pi180*theta(N))
        TAU(1,n) = C
        Ss1(n) = an / 2d0 * (ra0 * pii(1,n) + rb0 * tau(1,n))
        Ss2(n) = an / 2d0 * (ra0 * tau(1,n) + rb0 * pii(1,n))
c        write(7,*) n, ss1(n)
      end do
c*

      write(*,*) 'take a drink for Zhongli!'
      z = -1d0
      do I = 2, 5001
          an = an + 2d0
          an2 = an - 2d0
c* Bessel functions
         if(iu1.eq.iu0) then
            besY2 = an2 * ax * besY1 - besY0
         else
            besY2 = an2 * ax * besY1 - besY0 / factor
         end if
         if(dabs(besY2).gt.1d300) then
           besY2 = besY2 / factor
           iu2 = iu1 + 1
         end if
         besJ2 = (w1 + besY2 * besJ1) / besY1
c* Mie coefficients
      s = ru(I) / ri + I * AX
      s1 = s * besJ2 / fact(iu2) - besJ1 / fact(iu1)
      s2 = s * besY2 * fact(iu2) - besY1 * fact(iu1)
      ra1 = s1 / (s1 - s3 * s2)
c
      s = ru(I) * ri + i * ax
      s1 = s * besJ2 / fact(iu2) - besJ1 / fact(iu1)
      s2 = s * besY2 * fact(iu2) - besY1 * fact(iu1)
      rb1 = s1 / (s1 - s3 * s2)
c*
       if(i.gt.nterms) then
c        write(*,*) 'nterms, i=', nterms, i     
      	 write(*,*) 'No, Venti, put ur clothes back!!!'
      	 Qpr = -10
c         pause
        return
       end if
c* amplitude functions
      DO N = 1, K
       C = dcos(pi180*theta(N))
       if(i.eq.2) then
         PIi(i,n) = 3D0 * C
         TAU(i,n) = 6D0 * C**2 - 3D0
c         write(*,*) i, n, pii(i,n)
         go to 2
       end if
      PIi(I,n) = (C * (2d0 * I - 1d0) * pii(I-1,n) -
     *         I * pii(I-2,n))/(I - 1d0)
      TAU(I,n) = C * (PIi(I,n) - pii(I-2,n)) - (2d0 * I - 1d0) *
     *         (1-C**2) * PIi(I-1,n) + TAU(I-2,n)
    2  continue
        Ss1(n) = ss1(n) + an / (i * (i + 1d0)) *
     *           (ra1 * pii(i,n) + rb1 * tau(i,n))
        Ss2(n) = ss2(n) + an / (i * (i + 1d0)) *
     *           (ra1 * tau(i,n) + rb1 * pii(i,n))
      end do
      
c* efficiency factors
      z = -z
      rr = z * (i + 0.5d0) * (ra1 - rb1)
      r = r + rr
      ss = ss + (i - 1d0) * (i + 1d0) / i *
     *     (ra0 * dconjg(ra1) + rb0 * dconjg(rb1)) +
     *     an2 / i / (i - 1d0) * (ra0 * dconjg(rb0))
      qq = an * (ra1 + Rb1)
      Qext = Qext + qq
      Qsca = Qsca + an * (ra1 * dconjg(ra1) +
     *       rb1 * dconjg(rb1))
      if(dabs(qq / qext).lt.eps) go to 1
      
c* Bessel functions
          besJ0 = besJ1
          besJ1 = besJ2
          besY0 = besY1
          besY1 = besY2
            iu0 =   iu1
            iu1 =   iu2
            ra0 =   ra1
            rb0 =   rb1
      end do
      
      write(*,*) 'take a drink for hilichurl!'
c* end of cycle
   1  continue
c* efficiency factors (final calculations)
      write(*,*) 'take a drink for Klee!'
      Qext = b * Qext
      Qsca = b * Qsca
      Qbk = 2d0 * b * r * dconjg(r)
      Qpr = Qext - 2d0 * b * ss
      Qabs = Qext - Qsca
      albedo = Qsca / Qext
      g = (Qext - Qpr) / Qsca
c* scattering matrix (final calculations)
      DO N = 1, K
c        write(7,*) n, ss1(n)
c        write(7,*) n, ss2(n)
       sS3=DCONJG(Ss1(n))
       Ss4=DCONJG(Ss2(n))
       P1(N)=Ss1(n) * sS3
       P2(N)=Ss2(n) * Ss4
       P3(N)=(Ss1(n) * Ss4 + Ss2(n) * Ss3) / 2D0
       P4(N)=(Ss1(n) * Ss4 - Ss2(n) * Ss3) * sSs
      end do
      write(*,*) 'take a drink for Kaeya!'
      return
      end
c------------------------------------------------------
c aa-subroutine for calculations of the ratio of
c    derivative to the function for Bessel functions
c    of half order with complex argument: J'(n)/J(n).
c    The calculations are given by the recursive
c    expression ``from top to bottom'' beginning
c    from n=num.
c    ru-array of results.
c    a=1/x (a=2*pi*a(particle radius)/lambda -
c    size parameter).
c    ri - complex refractive index.
c August 1989, AO LGU
c------------------------------------------------------
      subroutine aa(a, ri, num, ru)
      implicit real*8 (a-h,o-q,t-z), complex*16 (r-s)
      dimension ru(num)
      s=  a / ri
      ru(num) = (num + 1D0) * s
      num1 = num - 1
      do j = 1, num1
        i = num - j
        i1 = i + 1
        s1 = i1 * s
        ru(i) = s1 - 1D0 / (ru(i1) + s1)
      end do
      return
      end
