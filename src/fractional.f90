module fractional
   use util, only: to_lower => struc2lc
   implicit none
   private
   integer, parameter :: GENERAL_EXACT_EVALUATION = 1, &
                         INTERPOLATION_EVALUATION = 2, &
                         SPECIAL_EXACT_EVALUATION = 3
   integer, parameter :: CAPUTO_FD = 1, &
                         CAPUTO_FABRIZIO_FD = 2
   ! normally, it should be ceiling(alpha)
   integer :: p = 2000
   ! level of fractional derivative
   real(8) :: alpha_
   ! used algorithm
   integer :: algorithm_ = -1
   ! used fractional derivative definition
   integer :: definition_ = -1
   ! interface to C function; see 2F2.c
   ! note, this routine returns 0 < hyp2F2 <= 1
   ! z should be <= 0
   interface
      real(8) function hyp2F2(a, b, c, d, z) bind(C, name="hyp2F2")
         real(8), intent(in), value :: a, b, c, d, z
      end function hyp2F2
   end interface
   ! internal global arrays for interpolation algorithm
   integer, parameter :: n_2F2 = 101, k_2F2 = 10
   real(8) eprr_2F2(0:n_2F2 - 1)
   real(8) f_2F2(0:n_2F2 - 1, 0:8), Bs_2F2(0:n_2F2 - 1, 0:8)
   real(8) t_2F2(n_2F2 + k_2F2)
   !
   public :: set_alpha_level, GTO_fractional_integral
   public :: p
contains
   subroutine set_alpha_level()
      character(len=80) :: tempstr
      write (*, *) "Type `get` for getting information about current fractional derivative parameters"
      write (*, *) "Type `set_p` for updating fractional derivative parameters with specific level of internal derivatives"
      write (*, *) "Type `set_auto` for updating fractional derivative parameters with automatic level of internal derivatives"
      do while (.true.)
         read (*, *) tempstr
         call to_lower(tempstr)
         select case (trim(tempstr))
         case ("get")
            call get_alpha_level()
            return
         case ("set_p")
            write (*, *) "Input alpha parameter; it should be <= 1"
            read (*, *) alpha_
            write (*, *) "Input level of internal derivative; it can be 0 or 1"
            read (*, *) p
            if (alpha_ > p) then
               write (*, *) "Alpha can not be more than p"
               write (*, *) "Try again..."
               cycle
            end if
         case ("set_auto")
            write (*, *) "Input alpha parameter; it should be <= 1"
            read (*, *) alpha_
            p = ceiling(alpha_)
            write (*, '(A,I0)') "Automatically chosen level is ", p
         case default
            write (*, *) "Try again..."
            cycle
         end select
         exit
      end do
      write (*, *)
      write (*, *) "There is two possible definition of fractional derivatives: `Caputo` and `Caputo-Fabrizio`"
      write (*, *) "Type `C` or `Caputo` for choosing Caputo fractional derivative defition"
      write (*, *) "Type `CF` or `Caputo-Fabrizio` or `CaputoFabrizio` for choosing Caputo-Fabrizio fractional derivative defition"
      do while (.true.)
         read (*, *) tempstr
         call to_lower(tempstr)
         select case (trim(tempstr))
         case ("c", "caputo")
            definition_ = CAPUTO_FD
         case ("cf", "caputofabrizio", "caputo-fabrizio")
            definition_ = CAPUTO_FABRIZIO_FD
         case default
            write (*, *) "Try again..."
            cycle
         end select
         exit
      end do
      write (*, *)
      write (*, *) "There is several ways for evaluation of fractional derivatives:"
      write (*, *) "Type `exact` for enabling evaluation of general analytic form of fractional derivative"
      if (definition_ == CAPUTO_FD) then
         write (*, *) "Type `approx` or `approximate` for interpolation of fractional derivative"
         write (*, *) "  This option is helpful for Caputo fractional derivative for avoiding of recomputing computationally-costly hypergeometric function"
         write (*, *) "Type `special` for enabling evaluation of specific analytic form of fractional derivative"
         write (*, *) "For special, only several pairs are available:"
         write (*, *) "    p      alpha"
         write (*, *) "    1       1.0"
         write (*, *) "    0       0.0"
      end if
      do while (.true.)
         read (*, *) tempstr
         call to_lower(tempstr)
         select case (trim(tempstr))
         case ("exact")
            algorithm_ = GENERAL_EXACT_EVALUATION
         case ("approx", "approximate")
            algorithm_ = INTERPOLATION_EVALUATION
            call prepare_interpolation()
            if (definition_ == CAPUTO_FABRIZIO_FD) algorithm_ = GENERAL_EXACT_EVALUATION
         case ("special")
            algorithm_ = SPECIAL_EXACT_EVALUATION
            if (definition_ == CAPUTO_FABRIZIO_FD) algorithm_ = GENERAL_EXACT_EVALUATION
         case default
            write (*, *) "Try again..."
            cycle
         end select
         exit
      end do
      call get_alpha_level()
   end subroutine set_alpha_level
   subroutine get_alpha_level()
      select case (definition_)
      case (CAPUTO_FD)
         write (*, *) "Caputo definion of fractional derivative will be used"
      case (CAPUTO_FABRIZIO_FD)
         write (*, *) "Caputo-Fabrizio definion of fractional derivative will be used"
      case default
         write (*, *) "Unknown type of fractional derivative"
      end select
      write (*, *) "Fractional derivative parameters:"
      write (*, *) "Alpha: ", alpha_
      write (*, '(A,I0,A,I0)') "Level of derivative is ", p, "; should be ", ceiling(alpha_)
      select case (algorithm_)
      case (GENERAL_EXACT_EVALUATION)
         write (*, *) "Evaluation of fractional derivative using general analytic form"
      case (INTERPOLATION_EVALUATION)
         write (*, *) "Evaluation of fractional derivative using interpolation"
      case (SPECIAL_EXACT_EVALUATION)
         write (*, *) "Evaluation of fractional derivative using specific analytic form"
      case default
         write (*, *) "Unknown algorithm"
      end select
      write (*, *)
   end subroutine get_alpha_level
   real(8) function GTO_fractional_integral(n, x)
      integer, intent(in) :: n
      real(8), intent(in) :: x
      select case (definition_)
      case (CAPUTO_FD)
         GTO_fractional_integral = GTO_fractional_integral_Caputo(n, x)
      case (CAPUTO_FABRIZIO_FD)
         GTO_fractional_integral = GTO_fractional_integral_CaputoFabrizio(n, x)
      case default
         error stop "Unknown fractional derivative definition"
      end select
   end function GTO_fractional_integral
   ! This routine evaluates the following expression:
   ! \frac{1}{\Gamma(p-\alpha)} \int\limits_0^1 \frac{t^n}{(1-t)^{p-\alpha-1}} exp(x t^2) dt
   ! using general, approximational of special algorithms
   ! n is integer
   ! x is real(8) <= 0
   real(8) function GTO_fractional_integral_Caputo(n, x)
      use bspline_sub_module, only: db1val
      real(8), parameter :: one = 1e0_8, zero = 0e0_8
      integer, intent(in) :: n
      real(8), intent(in) :: x
      real(8) :: top, down
      real(8) :: gtop, gdown
      real(8) :: val2F2
      integer :: iout, inbvx
      if (x > zero) then
         write (*, *) "x is ", x
         error stop "x greater than zero!"
      end if
      val2F2 = -one
      GTO_fractional_integral_Caputo = 0e0_8
      if (algorithm_ == GENERAL_EXACT_EVALUATION) then
         top = real(n + 1, kind=8)
         down = top + real(p, kind=8) - alpha_
         gtop = gamma(top)
         gdown = gamma(down)
         val2F2 = hyp2F2(top, top + one, down, down + one, x)
         GTO_fractional_integral_Caputo = gtop/gdown*val2F2
      else if (algorithm_ == INTERPOLATION_EVALUATION) then
         top = real(n + 1, kind=8)
         down = top + real(p, kind=8) - alpha_
         gtop = gamma(top)
         gdown = gamma(down)
         inbvx = 1
         ! normally, x shall be passed, but here some logic is broken since Bspline expects increasing values of x
         call db1val(-x, 0, t_2F2, n_2F2, k_2F2, Bs_2F2(:, n), val2F2, iout, inbvx)
         GTO_fractional_integral_Caputo = gtop/gdown*val2F2
      else if (algorithm_ == SPECIAL_EXACT_EVALUATION) then
         select case (p)
         case (1)
            if (alpha_ == one) then
               GTO_fractional_integral_Caputo = exp(x)
            end if
         case (0)
            if (alpha_ == zero) then
               GTO_fractional_integral_Caputo = exp(x)
            end if
         end select
         if (GTO_fractional_integral_Caputo /= -one) return
         write (*, *) "Alpha: ", alpha_, " p: ", p
         error stop "No specialization for given Alpha and p"
      else
         write (*, '(A,I0)') "Algorithm internal value is ", algorithm_
         error stop "Algorithm is not implemented!"
      end if
   end function GTO_fractional_integral_Caputo
   real(8) function GTO_fractional_integral_CaputoFabrizio(n, x)
      real(8), parameter :: two = 2e0_8, one = 1e0_8, zero = 0e0_8
      integer, intent(in) :: n
      real(8), intent(in) :: x
      integer, parameter :: factorials(0:8) = (/1, 1, 2, 6, 24, 120, 720, 5040, 40320/)
      integer :: Cnk(0:8)
      integer :: i, signi
      real(8) :: prefactor, a, b, u, w
      real(8) :: intsum
      if (x > zero) then
         write (*, *) "x is ", x
         error stop "x greater than zero!"
      end if
      do i = 0, n
         Cnk(i) = factorials(n)/factorials(i)/factorials(n - i)
      end do
      GTO_fractional_integral_CaputoFabrizio = -one
      intsum = zero
      select case (algorithm_)
      case (GENERAL_EXACT_EVALUATION)
         a = alpha_/(alpha_ - one)
         b = -x
         u = a/two/sqrt(b)
         w = sqrt(b)
         prefactor = one/(one - alpha_)*exp(a)*exp(u*u)/w**n/two
         if (mod(n, 2) == 1) then
            signi = -1
         else
            signi = 1
         end if
         do i = 0, n
            signi = -signi
            intsum = intsum + signi*Cnk(i)*u**(n - i)*(gammad((i + one)/two, u*u) - gammad((i + one)/two, (u + w)**2))
         end do
         GTO_fractional_integral_CaputoFabrizio = prefactor*intsum
      case default
         write (*, '(A,I0)') "Algorithm internal value is ", algorithm_
         error stop "Algorithm is not implemented!"
      end select
   end function GTO_fractional_integral_CaputoFabrizio
   subroutine prepare_interpolation()
      select case (definition_)
      case (CAPUTO_FD)
         call prepare_interpolation_Caputo()
      case (CAPUTO_FABRIZIO_FD)
         error stop "Caputo-Fabrizio fractional derivative can not be interpolated!"
      case default
         error stop "Unknown definition of fractional derivative in prepare_interpolation!"
      end select
   end subroutine prepare_interpolation
   ! this routine fills arrays for interpolation evaluation of 2F2 function
   subroutine prepare_interpolation_Caputo()
      use bspline_sub_module, only: db1ink
      real(8), parameter :: one = 1e0_8, zero = 0e0_8
      real(8) :: top, down
      integer :: i, j
      integer :: iout
      do i = 0, n_2F2 - 1
         eprr_2F2(i) = exp(0.1_8*real(i, 8)) - one
         do j = 0, 8
            if (i > 0 .and. f_2F2(max(i - 1, 0), j) < 1e-40_8) then
               f_2F2(i, j) = zero
            else
               top = real(j + 1, kind=8)
               down = top + real(p, kind=8) - alpha_
               ! we evaluate negative eprr
               f_2F2(i, j) = hyp2F2(top, top + one, down, down + one, -eprr_2F2(i))
            end if
         end do
         write (*, *) i, eprr_2F2(i), f_2F2(i, 0:4)
      end do
      do j = 0, 8
         ! here, eprr_2F2 is from 0 to +inf
         call db1ink(eprr_2F2, n_2F2, f_2F2(:, j), k_2F2, 0, t_2F2, Bs_2F2(:, j), iout)
      end do
   end subroutine prepare_interpolation_Caputo
   subroutine prepare_interpolation_CaputoFabrizio()
      error stop "prepare_interpolation_Caputo_Fabrizio is not implemented!"
   end subroutine prepare_interpolation_CaputoFabrizio
   function gammad(x, p, ifault)
!*****************************************************************************80
!
!! GAMMAD computes the Incomplete Gamma Integral
!
!  Auxiliary functions:
!
!    ALOGAM = logarithm of the gamma function,
!    ALNORM = algorithm AS66
!
!  Modified:
!
!    20 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by B Shea.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    B Shea,
!    Algorithm AS 239:
!    Chi-squared and Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 466-473.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, P, the parameters of the incomplete
!    gamma ratio.  0 <= X, and 0 < P.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error.
!    1, X < 0 or P <= 0.
!
!    Output, real ( kind = 8 ) GAMMAD, the value of the incomplete
!    Gamma integral.
!
      implicit none

      real(kind=8) a
      real(kind=8) an
      real(kind=8) arg
      real(kind=8) b
      real(kind=8) c
      real(kind=8), parameter :: elimit = -88.0D+00
      real(kind=8) gammad
! originally, ifault is not optional argument
      integer(kind=4), optional :: ifault
      real(kind=8), parameter :: oflo = 1.0D+37
      real(kind=8) p
      real(kind=8), parameter :: plimit = 1000.0D+00
      real(kind=8) pn1
      real(kind=8) pn2
      real(kind=8) pn3
      real(kind=8) pn4
      real(kind=8) pn5
      real(kind=8) pn6
      real(kind=8) rn
      real(kind=8), parameter :: tol = 1.0D-14
      logical upper
      real(kind=8) x
      real(kind=8), parameter :: xbig = 1.0D+08

      gammad = 0.0D+00
!
!  Check the input.
!
      if (x < 0.0D+00) then
         ifault = 1
         return
      end if

      if (p <= 0.0D+00) then
         ifault = 1
         return
      end if

      ifault = 0

      if (x == 0.0D+00) then
         gammad = 0.0D+00
         return
      end if
!
!  If P is large, use a normal approximation.
!
      if (plimit < p) then

         pn1 = 3.0D+00*sqrt(p)*((x/p)**(1.0D+00/3.0D+00) &
                                + 1.0D+00/(9.0D+00*p) - 1.0D+00)

         upper = .false.
         gammad = alnorm(pn1, upper)
         return

      end if
!
!  If X is large set GAMMAD = 1.
!
      if (xbig < x) then
         gammad = 1.0D+00
         return
      end if
!
!  Use Pearson's series expansion.
!  (Note that P is not large enough to force overflow in ALOGAM).
!  No need to test IFAULT on exit since P > 0.
!
      if (x <= 1.0D+00 .or. x < p) then

         arg = p*log(x) - x - alngam(p + 1.0D+00, ifault)
         c = 1.0D+00
         gammad = 1.0D+00
         a = p

         do

            a = a + 1.0D+00
            c = c*x/a
            gammad = gammad + c

            if (c <= tol) then
               exit
            end if

         end do

         arg = arg + log(gammad)

         if (elimit <= arg) then
            gammad = exp(arg)
         else
            gammad = 0.0D+00
         end if
!
!  Use a continued fraction expansion.
!
      else

         arg = p*log(x) - x - alngam(p, ifault)
         a = 1.0D+00 - p
         b = a + x + 1.0D+00
         c = 0.0D+00
         pn1 = 1.0D+00
         pn2 = x
         pn3 = x + 1.0D+00
         pn4 = x*b
         gammad = pn3/pn4

         do

            a = a + 1.0D+00
            b = b + 2.0D+00
            c = c + 1.0D+00
            an = a*c
            pn5 = b*pn3 - an*pn1
            pn6 = b*pn4 - an*pn2

            if (pn6 /= 0.0D+00) then

               rn = pn5/pn6

               if (abs(gammad - rn) <= min(tol, tol*rn)) then
                  exit
               end if

               gammad = rn

            end if

            pn1 = pn3
            pn2 = pn4
            pn3 = pn5
            pn4 = pn6
!
!  Re-scale terms in continued fraction if terms are large.
!
            if (oflo <= abs(pn5)) then
               pn1 = pn1/oflo
               pn2 = pn2/oflo
               pn3 = pn3/oflo
               pn4 = pn4/oflo
            end if

         end do

         arg = arg + log(gammad)

         if (elimit <= arg) then
            gammad = 1.0D+00 - exp(arg)
         else
            gammad = 1.0D+00
         end if

      end if

      return
   contains
      function alnorm(x, upper)
!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
         implicit none

         real(kind=8), parameter :: a1 = 5.75885480458D+00
         real(kind=8), parameter :: a2 = 2.62433121679D+00
         real(kind=8), parameter :: a3 = 5.92885724438D+00
         real(kind=8) alnorm
         real(kind=8), parameter :: b1 = -29.8213557807D+00
         real(kind=8), parameter :: b2 = 48.6959930692D+00
         real(kind=8), parameter :: c1 = -0.000000038052D+00
         real(kind=8), parameter :: c2 = 0.000398064794D+00
         real(kind=8), parameter :: c3 = -0.151679116635D+00
         real(kind=8), parameter :: c4 = 4.8385912808D+00
         real(kind=8), parameter :: c5 = 0.742380924027D+00
         real(kind=8), parameter :: c6 = 3.99019417011D+00
         real(kind=8), parameter :: con = 1.28D+00
         real(kind=8), parameter :: d1 = 1.00000615302D+00
         real(kind=8), parameter :: d2 = 1.98615381364D+00
         real(kind=8), parameter :: d3 = 5.29330324926D+00
         real(kind=8), parameter :: d4 = -15.1508972451D+00
         real(kind=8), parameter :: d5 = 30.789933034D+00
         real(kind=8), parameter :: ltone = 7.0D+00
         real(kind=8), parameter :: p = 0.398942280444D+00
         real(kind=8), parameter :: q = 0.39990348504D+00
         real(kind=8), parameter :: r = 0.398942280385D+00
         logical up
         logical upper
         real(kind=8), parameter :: utzero = 18.66D+00
         real(kind=8) x
         real(kind=8) y
         real(kind=8) z

         up = upper
         z = x

         if (z < 0.0D+00) then
            up = .not. up
            z = -z
         end if

         if (ltone < z .and. ((.not. up) .or. utzero < z)) then

            if (up) then
               alnorm = 0.0D+00
            else
               alnorm = 1.0D+00
            end if

            return

         end if

         y = 0.5D+00*z*z

         if (z <= con) then

            alnorm = 0.5D+00 - z*(p - q*y &
                                  /(y + a1 + b1 &
                                    /(y + a2 + b2 &
                                      /(y + a3))))

         else

            alnorm = r*exp(-y) &
                     /(z + c1 + d1 &
                       /(z + c2 + d2 &
                         /(z + c3 + d3 &
                           /(z + c4 + d4 &
                             /(z + c5 + d5 &
                               /(z + c6))))))

         end if

         if (.not. up) then
            alnorm = 1.0D+00 - alnorm
         end if

         return
      end function alnorm
      function alngam(xvalue, ifault)
!*****************************************************************************80
!
!! ALNGAM computes the logarithm of the gamma function.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Allan Macleod.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Allan Macleod,
!    Algorithm AS 245,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the Gamma function.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!    1, XVALUE is less than or equal to 0.
!    2, XVALUE is too big.
!
!    Output, real ( kind = 8 ) ALNGAM, the logarithm of the gamma function of X.
!
         implicit none

         real(kind=8) alngam
         real(kind=8), parameter :: alr2pi = 0.918938533204673D+00
         integer(kind=4) ifault
         real(kind=8), dimension(9) :: r1 = (/ &
                                       -2.66685511495D+00, &
                                       -24.4387534237D+00, &
                                       -21.9698958928D+00, &
                                       11.1667541262D+00, &
                                       3.13060547623D+00, &
                                       0.607771387771D+00, &
                                       11.9400905721D+00, &
                                       31.4690115749D+00, &
                                       15.2346874070D+00/)
         real(kind=8), dimension(9) :: r2 = (/ &
                                       -78.3359299449D+00, &
                                       -142.046296688D+00, &
                                       137.519416416D+00, &
                                       78.6994924154D+00, &
                                       4.16438922228D+00, &
                                       47.0668766060D+00, &
                                       313.399215894D+00, &
                                       263.505074721D+00, &
                                       43.3400022514D+00/)
         real(kind=8), dimension(9) :: r3 = (/ &
                                       -2.12159572323D+05, &
                                       2.30661510616D+05, &
                                       2.74647644705D+04, &
                                       -4.02621119975D+04, &
                                       -2.29660729780D+03, &
                                       -1.16328495004D+05, &
                                       -1.46025937511D+05, &
                                       -2.42357409629D+04, &
                                       -5.70691009324D+02/)
         real(kind=8), dimension(5) :: r4 = (/ &
                                       0.279195317918525D+00, &
                                       0.4917317610505968D+00, &
                                       0.0692910599291889D+00, &
                                       3.350343815022304D+00, &
                                       6.012459259764103D+00/)
         real(kind=8) x
         real(kind=8) x1
         real(kind=8) x2
         real(kind=8), parameter :: xlge = 5.10D+05
         real(kind=8), parameter :: xlgst = 1.0D+30
         real(kind=8) xvalue
         real(kind=8) y

         x = xvalue
         alngam = 0.0D+00
!
!  Check the input.
!
         if (xlgst <= x) then
            ifault = 2
            return
         end if

         if (x <= 0.0D+00) then
            ifault = 1
            return
         end if

         ifault = 0
!
!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
!
         if (x < 1.5D+00) then

            if (x < 0.5D+00) then

               alngam = -log(x)
               y = x + 1.0D+00
!
!  Test whether X < machine epsilon.
!
               if (y == 1.0D+00) then
                  return
               end if

            else

               alngam = 0.0D+00
               y = x
               x = (x - 0.5D+00) - 0.5D+00

            end if

            alngam = alngam + x*(((( &
                                   r1(5)*y &
                                   + r1(4))*y &
                                   + r1(3))*y &
                                  + r1(2))*y &
                                 + r1(1))/(((( &
                                             y &
                                             + r1(9))*y &
                                             + r1(8))*y &
                                            + r1(7))*y &
                                           + r1(6))

            return

         end if
!
!  Calculation for 1.5 <= X < 4.0.
!
         if (x < 4.0D+00) then

            y = (x - 1.0D+00) - 1.0D+00

            alngam = y*(((( &
                          r2(5)*x &
                          + r2(4))*x &
                          + r2(3))*x &
                         + r2(2))*x &
                        + r2(1))/(((( &
                                    x &
                                    + r2(9))*x &
                                    + r2(8))*x &
                                   + r2(7))*x &
                                  + r2(6))
!
!  Calculation for 4.0 <= X < 12.0.
!
         else if (x < 12.0D+00) then

            alngam = (((( &
                        r3(5)*x &
                        + r3(4))*x &
                        + r3(3))*x &
                       + r3(2))*x &
                      + r3(1))/(((( &
                                  x &
                                  + r3(9))*x &
                                  + r3(8))*x &
                                 + r3(7))*x &
                                + r3(6))
!
!  Calculation for 12.0 <= X.
!
         else

            y = log(x)
            alngam = x*(y - 1.0D+00) - 0.5D+00*y + alr2pi

            if (x <= xlge) then

               x1 = 1.0D+00/x
               x2 = x1*x1

               alngam = alngam + x1*(( &
                                     r4(3)* &
                                     x2 + r4(2))* &
                                     x2 + r4(1))/(( &
                                                  x2 + r4(5))* &
                                                  x2 + r4(4))

            end if

         end if

         return
      end function alngam
   end function gammad
end module fractional
