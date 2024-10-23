module moi_aux
   private
   real(kind=8), allocatable :: dxyzsel(:, :), val(:)
   integer :: ndmoi
   public :: dxyzsel, val, ndmoi, moi
contains

   subroutine moi(sqdfmoi, sqdfdistmoi, outputidxmoi, tminmoi, tmaxmoi, domsizemoi, &
                  dxmoi, dymoi, dzmoi, wtsmoi, moimatmoi, eigvecmoi, eigvalmoi, perc99valmoi, debugmoi)
      !-----------------------------------------------------------------------
      !
      ! Subroutine to find the moment of inertia(MOI) matrix. References of equations
      ! can be at check in Hadavand and Deustch (2020) and Hassanpour and Deustch (2007).
      !
      !                      MOI Subroutine
      !                      **************
      !
      !
      ! INPUT PARAMETERS:
      ! * sqdfmoi: array with pseudocorrelations of pairs.
      ! * tminmoi, and tmaxmoi: trimming limits (they are equal of the initial declaration)
      ! * outputidxmoi : number of pairs.
      ! * dxmoi, dymoi, dzmoi: arrays with distances (as vectors: including negatives)
      !
      !
      ! OUTPUT PARAMETERS:
      ! * Eigenvectors, eigenvalues and MOI matrix
      !
      !
      ! INTERNAL VARIABLES
      ! * nd, i: integers to keep counts and loops.
      !
      ! Jacobi subroutine used. Included at the end.
      !
      !-----------------------------------------------------------------------
      !************* Variables *************
      use sortem

      implicit none

      logical :: testfl

      ! Subroutine parameters:
      integer, intent(in) :: outputidxmoi
      real(kind=8), intent(in) :: tminmoi, tmaxmoi, domsizemoi
      real(kind=8), dimension(outputidxmoi), intent(in) :: sqdfmoi, &
                                                           sqdfdistmoi, dxmoi, dymoi, dzmoi, wtsmoi
      logical, intent(in) :: debugmoi
      integer, dimension(outputidxmoi) :: index_
      integer :: ncellx, ncelly, ncellz, ncellt, icellx, icelly, icellz, &
                 icell
      real(kind=8), allocatable :: npoint(:)
      real(kind=8) :: dxmin, dymin, dzmin, dxmax, dymax, dzmax, &
                      power, perc_min, perc_max
      integer :: p_min, p_max
      real(kind=8), parameter :: BIGDBLE = 1d21, SMALLDBLE = 1d-6

      ! Calculate the eigenvectors and eigenvalues
      double precision, parameter:: abserr = 1.0e-09
      integer i, j, k, nrot
      ! Output variables:
      real(kind=8), intent(out) :: eigvecmoi(3, 3), moimatmoi(3, 3), eigvalmoi(3)
      real(kind=8), intent(out) :: perc99valmoi
      real(kind=8) :: min_corr, max_corr
      ! Internal variables:

      integer :: nd, perc, nn, jj, ss, nb
      !real(kind=4), dimension(outputidxmoi) :: val
      real(kind=8), dimension(outputidxmoi) ::  temp_arr, temp_arr2
      real(kind=8), allocatable :: distx(:), disty(:), distz(:), wts(:), &
                                   sum_moi(:), moi_xyz(:, :), dist(:)
      real(kind=8) :: max_rad_s, dataval, mod99, zero

      ! End of variable declaration

      !
      ! Read through the arrays and storage in convenient way:
      !
      nd = 0

      !Get P99 to restrict the values that enters MOI calculation

      ! Order the squared differences, find the 99 th percentile, create an array with only those percentiles
      ! Find the mean and variance, then standardize every difference of pairs.
      nn = 0

      do jj = 1, outputidxmoi
         nn = nn + 1
         temp_arr(nn) = sqdfmoi(jj)
      end do

      call dblemodsortem(temp_arr, outputidxmoi, 0, temp_arr2)
      ss = size(temp_arr)
      ! Find P99:
      p_max = 99

      ! If the array is divisible by 99:
      mod99 = mod(ss, 99)
      zero = 0.0

      if (mod99 .eq. zero) then
         perc_max = int(0.99*ss)
         max_corr = temp_arr(perc_max)

         ! If it is not:
      else
         perc_max = int(0.99*ss)
         max_corr = (temp_arr(perc_max) + temp_arr(perc_max + 1))/2

      end if
      perc_max = floor(0.99*ss)
      perc99valmoi = max_corr

      p_min = 1
      ! If the array is divisible by 99:
      if (mod(ss, 1) .eq. 0) then
         perc_min = int(0.01*ss)
         min_corr = temp_arr(perc_min)
         ! If it is not:
      else
         perc_min = int(0.01*ss)
         min_corr = (temp_arr(perc_min) + temp_arr(perc_min + 1))/2
      end if

      ! Just consider values below P99
      do i = 1, outputidxmoi
         if (sqdfmoi(i) .le. max_corr) then
            nd = nd + 1
         else
            continue
         end if
      end do

      !Allocate
      allocate (wts(nd))
      allocate (val(nd))
      allocate (distx(nd))
      allocate (disty(nd))
      allocate (distz(nd))
      allocate (dxyzsel(3, nd))
      allocate (dist(nd))

      ! Allocate the selected values:
      ndmoi = nd
      nd = 0
      do i = 1, outputidxmoi
         if (sqdfmoi(i) .le. max_corr) then
            nd = nd + 1
            val(nd) = sqdfmoi(nd)
            distx(nd) = abs(dxmoi(nd))
            disty(nd) = abs(dymoi(nd))
            distz(nd) = abs(dzmoi(nd))
            dxyzsel(1, nd) = dxmoi(nd)
            dxyzsel(2, nd) = dymoi(nd)
            dxyzsel(3, nd) = dzmoi(nd)
            dist(nd) = sqdfdistmoi(nd)
            wts(nd) = 1
         else
            continue
         end if
      end do

      if (nd .lt. 2) then
         write (*, *) ' too few data for calculating MOI', nd
         stop
      end if

      i = 0
      j = 0
      k = 0

      !
      ! Main loop to calculate the elements of MOI matrix
      ! Set initial values to zero in case of 2D data.

      moimatmoi = 0
      power = 2.5

      do i = 1, nd - 1
         ! for elements off the diagonal
         moimatmoi(2, 1) = moimatmoi(2, 1) - val(i)*wts(i)*dxmoi(i)*dymoi(i)/dist(i)**power
         moimatmoi(3, 1) = moimatmoi(3, 1) - val(i)*wts(i)*dxmoi(i)*dzmoi(i)/dist(i)**power
         moimatmoi(3, 2) = moimatmoi(3, 2) - val(i)*wts(i)*dymoi(i)*dzmoi(i)/dist(i)**power
         ! for elements in the diagonal
         moimatmoi(1, 1) = moimatmoi(1, 1) + val(i)*wts(i)*(dymoi(i)**2 + dzmoi(i)**2)/dist(i)**power
         moimatmoi(2, 2) = moimatmoi(2, 2) + val(i)*wts(i)*(dxmoi(i)**2 + dzmoi(i)**2)/dist(i)**power
         moimatmoi(3, 3) = moimatmoi(3, 3) + val(i)*wts(i)*(dxmoi(i)**2 + dymoi(i)**2)/dist(i)**power
      end do

      ! assign the values to similar terms:
      moimatmoi(1, 2) = moimatmoi(2, 1)
      moimatmoi(1, 3) = moimatmoi(3, 1)
      moimatmoi(2, 3) = moimatmoi(3, 2)

      ! show results in screen
      if (debugmoi) then
         do i = 1, 3
            write (*, '(3f15.3)') (moimatmoi(i, j), j=1, 3)
         end do
      end if

      ! Get the eigenvectors and eigenvalues:
      call jacobi(moimatmoi, 3, 3, eigvalmoi, eigvecmoi, nrot)
      deallocate (wts, val, distx, disty, distz, dxyzsel, dist)
      return

   end subroutine moi

   subroutine jacobi(a, n, np, d, v, nrot)
      !-----------------------------------------------------------------------
      ! Code from Numerical Recipes in Fortran (Press, Teukolsky, Vetterling,
      !                                         Flannery) 1992
      ! Used under the original license included with book and code purchase
      !-----------------------------------------------------------------------
      ! Computes all eigenvalues and eigenvectors of a real symmetric
      ! matrix a, which is of size n by n, stored in a physical np by
      ! np array. On output, elements of a above the diagonal are dest-
      ! royed. d returns eigenvalues of a in its first elements. v is
      ! a matrix with the same logical and physical dimensions as a,
      ! whose columns contain, on output, the normalized eigenvectors
      ! of a. nrot returns the number of Jacobi rotationsthat were required.
      implicit none

      integer, intent(in) :: n, np
      integer, intent(out) :: nrot
      integer :: NMAX
      real(kind=8):: a(np, np)
      real(kind=8), intent(out) :: d(np), v(np, np)
      parameter(NMAX=500)

      integer i, ip, iq, j
      real c, g, h, s, sm, t, tau, theta, tresh, b(NMAX), z(NMAX)

      do ip = 1, n          ! Initialize to the Identity matrix
         do iq = 1, n
            v(ip, iq) = 0.
         end do
         v(ip, ip) = 1.
      end do

      do ip = 1, n
         b(ip) = a(ip, ip)  ! Initialize b and d to the diagonal a.
         d(ip) = b(ip)
         z(ip) = 0.
      end do
      nrot = 0
      do i = 1, 100
         sm = 0.
         do ip = 1, n - 1
            do iq = ip + 1, n
               sm = sm + abs(a(ip, iq))
            end do
         end do
         if (sm == 0.) return
         if (i < 4) then
            tresh = 0.2*sm/n**2
         else
            tresh = 0.
         end if
         do ip = 1, n - 1
            do iq = ip + 1, n
               g = 100.*abs(a(ip, iq))

               if ((i > 4) .and. (abs(d(ip)) + g .eq. abs(d(ip))) .and. (abs(d(iq)) + g .eq. abs(d(iq)))) then
                  a(ip, iq) = 0.
               else if (abs(a(ip, iq)) > tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h) + g == abs(h)) then
                     t = a(ip, iq)/h
                  else
                     theta = 0.5*h/a(ip, iq)
                     t = 1./(abs(theta) + sqrt(1.+theta**2))
                     if (theta < 0.) t = -t
                  end if

                  c = 1./sqrt(1 + t**2)
                  s = t*c
                  tau = s/(1.+c)
                  h = t*a(ip, iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip, iq) = 0.
                  do j = 1, ip - 1
                     g = a(j, ip)
                     h = a(j, iq)
                     a(j, ip) = g - s*(h + g*tau)
                     a(j, iq) = h + s*(g - h*tau)
                  end do

                  do j = ip + 1, iq - 1
                     g = a(ip, j)
                     h = a(j, iq)
                     a(ip, j) = g - s*(h + g*tau)
                     a(j, iq) = h + s*(g - h*tau)
                  end do

                  do j = iq + 1, n
                     g = a(ip, j)
                     h = a(iq, j)
                     a(ip, j) = g - s*(h + g*tau)
                     a(iq, j) = h + s*(g - h*tau)
                  end do

                  do j = 1, n
                     g = v(j, ip)
                     h = v(j, iq)
                     v(j, ip) = g - s*(h + g*tau)
                     v(j, iq) = h + s*(g - h*tau)
                  end do
                  nrot = nrot + 1
               end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.
         end do
      end do
      write (*, *) 'ERROR: too many iterations in jacobi subroutine'

      return
   end subroutine jacobi

end module moi_aux
