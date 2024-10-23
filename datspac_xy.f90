module spacing
   private
   public :: acorni1, datspac_horz, datspac_vert, samplespac_vert
contains

   double precision function acorni1(idum) result(acorni11)
      !-----------------------------------------------------------------------
      !
      ! Fortran implementation of ACORN random number generator of order less
      ! than or equal to 12 (higher orders can be obtained by increasing the
      ! parameter value MAXORD).
      !
      !
      ! NOTES: 1. The variable idum is a dummy variable. The common block
      !           IACO is used to transfer data into the function.
      !
      !        2. Before the first call to ACORN the common block IACO must
      !           be initialised by the user, as follows. The values of
      !           variables in the common block must not subsequently be
      !           changed by the user.
      !
      !             KORDEI - order of generator required ( must be =< MAXORD)
      !
      !             MAXINT - modulus for generator, must be chosen small
      !                      enough that 2*MAXINT does not overflow
      !
      !             ixv(1) - seed for random number generator
      !                      require 0 < ixv(1) < MAXINT
      !
      !             (ixv(I+1),I=1,KORDEI)
      !                    - KORDEI initial values for generator
      !                      require 0 =< ixv(I+1) < MAXINT
      !
      !        3. After initialisation, each call to ACORN generates a single
      !           random number between 0 and 1.
      !
      !        4. An example of suitable values for parameters is
      !
      !             KORDEI   = 10
      !             MAXINT   = 2**30
      !             ixv(1)   = an odd integer in the (approximate) range
      !                        (0.001 * MAXINT) to (0.999 * MAXINT)
      !             ixv(I+1) = 0, I=1,KORDEI
      !
      !
      !
      ! Author: R.S.Wikramaratna,                           Date: October 1990
      !-----------------------------------------------------------------------
      implicit double precision(a - h, o - z)
      parameter(KORDEI=12, MAXOP1=KORDEI + 1, MAXINT=2**30)
      common/iaco/ixv(MAXOP1)
      do i = 1, KORDEI
         ixv(i + 1) = (ixv(i + 1) + ixv(i))
         if (ixv(i + 1) .ge. MAXINT) ixv(i + 1) = ixv(i + 1) - MAXINT
      end do
      acorni11 = dble(ixv(KORDEI + 1))/MAXINT
      return
   end function acorni1

   subroutine datspac_horz(xyzsh, ndsh, idcolsh, idsh, avgminsh, avgavgsh, avgamaxsh, debugsh)
      !-----------------------------------------------------------------------
      ! Modified from datspac program as subroutine

      ! Compute the distribution of data spacings in a 2D domain
      ! The intent is as to enable selection of a cell size for cell declus
      ! Assuming (mostly) vertical boreholes, a 2D treatment is sufficient

      ! Author: D. Khan                                     Date: January 2014
      ! Modified by: Luis Davila                            Date: February 2023
      !-----------------------------------------------------------------------

      use sortem
      implicit none

      ! Input variables from main
      integer, intent(in) :: ndsh, idcolsh!number of data and idsh column number
      real(kind=8), dimension(ndsh, 3), intent(in):: xyzsh
      real(kind=8), dimension(ndsh), intent(in) :: idsh
      logical, intent(in) :: debugsh

      !Internal variables
      real, parameter :: EPSLON = 1.0e-2, MAXCOL = 50
      real, parameter::KORDEI = 12, MAXOP1 = KORDEI + 1, MAXINT = 2**30! ACORN parameters
      common/iaco/ixv(MAXOP1)

      real(kind=8) rnd, acorni11, zdist
      real(kind=8), allocatable :: distmin(:), distmax(:), avgdist(:), dist(:), tmpdist(:), dh(:), tmpdh(:)
      real(kind=8), allocatable :: xloc(:), yloc(:), zloc(:), tmpxloc(:), tmpyloc(:), tmpzloc(:)
      real(kind=8) maxspac, maxspac2, dx, dy, h2, minspac, minspac2, max_dist
      real(kind=4) :: var(MAXCOL)

      integer icx, icy, ndist, na, nb, zcount
      character*512 :: str
      logical testfl

      integer :: idh, idd, ixv, i, j, ndh, ndmax, idat, ndt, idf, j_acc
      real(kind=8) :: tmp, midx, midy
      integer :: idum
      real(kind=8) :: sumdist
      !outputs
      real(kind=8), intent(out) :: avgminsh, avgavgsh, avgamaxsh

      ! Assign input variables to internal variables
      ! The minspacing will be set to 1m assuming sample spaces are greater than this.
      minspac = 1
      minspac2 = minspac*minspac
      ! The max spacing will be calculated using the furthest coordinates
      max_dist = sqrt((maxval(xyzsh(:, 1)) - minval(xyzsh(:, 1)))**2 + &
                      (maxval(xyzsh(:, 2)) - minval(xyzsh(:, 2)))**2 + &
                      (maxval(xyzsh(:, 3)) - minval(xyzsh(:, 3)))**2)
      maxspac = max_dist
      maxspac2 = maxspac*maxspac
      !drillhole
      if (idcolsh .ne. 0) then
         idh = idcolsh
         idd = 1 !if 1 to use DHID
      else
         idd = 0
      end if

      ! random number seed for despiking
      ixv(1) = 69069
      do i = 1, 1000
         rnd = acorni1(idum)
      end do

      ! Read data file
      !

      if (idd == 1) then
         allocate (dh(ndsh))
         dh = idsh
         call dblemodsortem(dh, ndsh, 0)
         ndh = 1
         do i = 2, ndsh
            if (dh(i) /= dh(i - 1)) ndh = ndh + 1
         end do
      end if

      ! Initialize x,yloc arrays
      allocate (xloc(ndsh), yloc(ndsh), zloc(ndsh))
      allocate (distmin(ndsh), avgdist(ndsh), distmax(ndsh), tmpdist(ndsh))
      ! Setting the number of ndmax as 4 to look the closest neighbors
      ndmax = 4
      ndist = min(ndsh, ndmax)
      allocate (dist(ndist))

      ! Read in the values

      xloc = xyzsh(:, 1)
      yloc = xyzsh(:, 2)
      zloc = xyzsh(:, 3)

      if (idd == 1) then
         allocate (tmpxloc(ndh), tmpyloc(ndh), tmpdh(ndh))
         call dblemodsortem(dh, ndsh, 0)
         call dblemodsortem(dh, ndsh, 2, xloc, yloc)
         idh = 1
         do i = 2, ndsh
            if (dh(i) /= dh(i - 1)) idh = idh + 1
            tmpxloc(idh) = xloc(i)
            tmpyloc(idh) = yloc(i)
         end do
         deallocate (xloc, yloc)
         allocate (xloc(ndh), yloc(ndh))
         xloc = tmpxloc
         yloc = tmpyloc
         deallocate (tmpxloc, tmpyloc)
         ndt = ndh
      end if
      if (xloc(1) .lt. 0) xloc(1) = minval(xyzsh(:, 1))
      if (yloc(1) .lt. 0) yloc(1) = minval(xyzsh(:, 2))

      ! Loop over all the data locations in the map:
      distmin = 0.0
      avgdist = 0.0
      do idat = 1, ndt

         ! Find the nearest samples up to ndmax neighbors. First initialize the counter arrays:

         na = 0
         nb = 0
         dist = 1.0e+20

         dx = 0.0
         dy = 0.0
         h2 = 0.0

         ! Scan all the samples (up to maxspac):

         do idf = 1, ndt
            if (idf == idat) cycle
            dx = xloc(idf) - xloc(idat)
            dy = yloc(idf) - yloc(idat)
            h2 = dx*dx + dy*dy
            if (h2 .gt. maxspac2) cycle
            if (h2 .lt. minspac2) cycle
            na = na + 1
            tmpdist(na) = sqrt(h2)
         end do

         ! sort the distances and remove duplicate locations
         call dblemodsortem(tmpdist, na, 0)
         do i = 1, na
            ! Consider up to ndmax neighboring sample locations:
            if (tmpdist(i + 1) - tmpdist(i) .lt. minspac) cycle
            if (nb == ndmax) exit
            nb = nb + 1
            if (nb .eq. 0) exit
            dist(nb) = tmpdist(i)
         end do
         if (nb .eq. 0) go to 3
         ! Get the minimum spacing from this datum
         distmin(idat) = dist(1)

         ! Get the maximum spacing from this datum
         distmax(idat) = dist(nb)

         ! Find the average spacing from this datum
         sumdist = 0.0
         do i = 1, nb
            sumdist = sumdist + dist(i)
         end do
         avgdist(idat) = sumdist/nb
         ! End Main Loop
3        continue
      end do

      ! Sort min, max and avg spacings to turn them into a CDF:

      call dblemodsortem(distmin, ndt, 0)
      call dblemodsortem(avgdist, ndt, 0)
      call dblemodsortem(distmax, ndt, 0)
      ! despike the values

      do i = ndt, 2, -1
         if (distmin(i - 1) == distmin(i)) then
            rnd = acorni1(idum)
            if (rnd .gt. 0.5) then
               distmin(i) = distmin(i) + rnd*EPSLON
               tmp = real(rnd)*EPSLON
            else
               distmin(i) = distmin(i) - rnd*EPSLON
            end if
         end if
         if (avgdist(i - 1) == avgdist(i)) then
            rnd = acorni1(idum)
            if (rnd .gt. 0.5) then
               avgdist(i) = avgdist(i) + rnd*EPSLON
            else
               avgdist(i) = avgdist(i) - rnd*EPSLON
            end if
         end if
         if (distmax(i - 1) == distmax(i)) then
            rnd = acorni1(idum)
            if (rnd .gt. 0.5) then
               distmax(i) = distmax(i) + rnd*EPSLON
               tmp = real(rnd)*EPSLON
            else
               distmax(i) = distmax(i) - rnd*EPSLON
            end if
         end if
      end do

      ! average of average data spacing array

      avgavgsh = 0.0
      do i = 1, ndt
         avgavgsh = avgavgsh + avgdist(i)
      end do
      avgavgsh = avgavgsh/ndt

      ! average of minimum data spacing array

      avgminsh = 0.0
      do i = 1, ndt
         avgminsh = avgminsh + distmin(i)
      end do
      avgminsh = avgminsh/ndt

      ! average of maximum data spacing array

      avgamaxsh = 0.0
      do i = 1, ndt
         avgamaxsh = avgamaxsh + distmax(i)
      end do
      avgamaxsh = avgamaxsh/ndt
      if (debugsh) then
         write (*, *) 'Inferred spacing:'
         write (*, *) ' Avg. Min. spacing, Avg. avg. spacing, and Avg. max. data spacing:'
         write (*, '(3f18.3)') avgminsh, avgavgsh, avgamaxsh
         write (*, *)
      end if
      return

   end subroutine datspac_horz

   subroutine datspac_vert(xyzsv, ndsv, idcolsv, idsv, zavrsv, debugsv, lavrsv)
      !-----------------------------------------------------------------------
      ! Modified from datspac program as subroutine (D. Khan)

      ! Get the average sampling spacing in the Z direction and the average
      ! drillhole lenght in the Z direction

      ! Author: Luis Davila                           Date: February 2023
      !-----------------------------------------------------------------------

      use sortem
      implicit none

      ! Input variables from main
      integer, intent(in) :: ndsv, idcolsv!number of data and idsv column number
      real(kind=8), dimension(ndsv, 3), intent(in):: xyzsv
      real(kind=8), dimension(ndsv), intent(in) :: idsv
      logical, intent(in):: debugsv

      !Internal variables
      real(kind=8), dimension(ndsv, 4) :: xyzsv_idsv
      real, parameter :: EPSLON = 1.0e-2, MAXCOL = 50
      real, parameter::KORDEI = 12, MAXOP1 = KORDEI + 1, MAXINT = 2**30! ACORN parameters
      common/iaco/ixv(MAXOP1)

      real(kind=8) rnd, acorni1, zdist, ztemp, lenght_temp
      real(kind=8), allocatable :: distmin(:), distmax(:), avgdist(:), dist(:), tmpdist(:), dh(:), tmpdh(:)
      real(kind=8), allocatable :: xloc(:), yloc(:), zloc(:), tmpxloc(:), tmpyloc(:), tmpzloc(:)
      real(kind=8) maxspac, maxspac2, dx, dy, h2, minspac, minspac2, max_dist
      real(kind=4) :: var(MAXCOL)
      real(kind=8), allocatable :: zdh(:), zdistdh(:), lenght(:)

      integer icx, icy, ndist, na, nb, zcount, meannsamp, nsamples
      character*512 :: str
      logical testfl

      integer :: idh, idd, ixv, i, j, k, l, ndh, ndmax, idat, ndt, idf, n, j_acc
      real(kind=8) :: tmp
      integer :: idum
      real(kind=8) :: sumdist
      !outputs
      real(kind=8), intent(out) :: zavrsv, lavrsv

      ! Assign input variables to internal variables
      ! The minspacing will be set to 2m since it is not likely having drilholes less
      ! than 2m apart for exploration stages
      minspac = 2
      minspac2 = minspac*minspac
      ! The max spacing will be set to 0.1 the domain size
      max_dist = sqrt((maxval(xyzsv(:, 1)) - minval(xyzsv(:, 1)))**2 + &
                      (maxval(xyzsv(:, 2)) - minval(xyzsv(:, 2)))**2 + &
                      (maxval(xyzsv(:, 3)) - minval(xyzsv(:, 3)))**2)
      maxspac = max_dist
      maxspac2 = maxspac*maxspac
      !drillhole
      if (idcolsv .ne. 0) then
         idh = idcolsv
         idd = 1 !if 1 to use DHID
      else
         idd = 0
      end if
      !Find the number of drillholes and allocate the matrix

      if (idd == 1) then
         allocate (dh(ndsv))
         dh = idsv
         call dblemodsortem(dh, ndsv, 0)
         ndh = 1
         ! sort first and then count the number of changes to find the
         ! number of dh
         do i = 2, ndsv
            if (dh(i) /= dh(i - 1)) ndh = ndh + 1
         end do
      end if
      allocate (zdistdh(ndh)) ! to save the average distance spacing of each ddh
      allocate (lenght(ndh)) ! to save the lenght of each ddh

      ! Initialize xyzsv_idsv
      ! Setting the number of ndmax as 2 to look the closest neighbor
      ndmax = 2
      ndist = min(ndsv, ndmax)
      allocate (dist(ndist))

      ! Read in the values

      xyzsv_idsv(:, 1) = xyzsv(:, 1)
      xyzsv_idsv(:, 2) = xyzsv(:, 2)
      xyzsv_idsv(:, 3) = xyzsv(:, 3)
      xyzsv_idsv(:, 4) = idsv
      ! Find the average distances in Z direction
      ! Grouping by DH
      nsamples = 0 ! to start the count of number of samples
      if (idd == 1) then
         call dblemodsortem(xyzsv_idsv(:, 4), ndsv, 3, xyzsv_idsv(:, 3), xyzsv_idsv(:, 2), xyzsv_idsv(:, 1))
         idh = 1
         j_acc = 0
         ndh = 1
         do i = 2, ndsv
            if (i .lt. j_acc + 1) cycle
            j = 0
            if (dh(i) /= dh(i - 1)) ndh = ndh + 1
            do n = 1, ndsv
               if (dh(n) .gt. dh(i)) exit
               if (dh(i) .eq. dh(n)) then
                  j = j + 1
               end if
            end do
            allocate (zdh(j))
            j_acc = j_acc + j
            zdh = xyzsv_idsv(j_acc - j + 1:j_acc, 3)
            ! count the number of samples per drillhole
            ! to calculate and report as average.
            nsamples = nsamples + j
            !sort the z values of this dh
            call dblemodsortem(zdh, j, 0)
            zdist = 0
            zcount = 0
            do l = 2, j
               zdist = zdist + abs(zdh(l) - zdh(l - 1))
               zcount = zcount + 1
            end do
            zdistdh(ndh) = zdist/zcount
            ! Find the lenght of this drillhole
            lenght(ndh) = abs(maxval(zdh) - minval(zdh))
            deallocate (zdh)

         end do
         ztemp = 0
         lenght_temp = 0
         do i = 1, ndh
            ztemp = ztemp + zdistdh(i)
            lenght_temp = lenght_temp + lenght(i)
         end do
         zavrsv = ztemp/ndh
         lavrsv = lenght_temp/ndh
         meannsamp = int(nsamples/ndh)
      end if

      if (debugsv) then
         write (*, *) ' Average Z spacing '
         write (*, '(1f18.3)') zavrsv
         write (*, *) ' Average drillhole Z lenght '
         write (*, '(1f18.3)') lavrsv
         write (*, *) ' Average number of samples '
         write (*, '(i18)') meannsamp
         write (*, *) ' Number of Drillholes '
         write (*, '(i18)') ndh
         write (*, *) ''
      end if
      return
   end subroutine datspac_vert

   subroutine samplespac_vert(xyzsv, ndsv, idcolsv, idsv, xyzavr, debugsv, lavrsv)
      !-----------------------------------------------------------------------
      ! Get the average sampling spacing downhole and the average drillhole
      ! lenght, whether inclined or vertical drillhole.

      ! Author: Luis Davila                           Date: March 2024
      !-----------------------------------------------------------------------

      use sortem
      implicit none

      ! Input variables from main
      integer, intent(in) :: ndsv, idcolsv
      real(kind=8), dimension(ndsv, 3), intent(in):: xyzsv
      real(kind=8), dimension(ndsv), intent(in) :: idsv
      logical, intent(in):: debugsv

      !Internal variables
      real(kind=8), dimension(ndsv, 4) :: xyzsv_idsv
      real, parameter :: EPSLON = 1.0e-2, MAXCOL = 50
      real, parameter::KORDEI = 12, MAXOP1 = KORDEI + 1, MAXINT = 2**30! ACORN parameters
      common/iaco/ixv(MAXOP1)

      real(kind=8) rnd, acorni1, xdist, ydist, zdist, xtemp, ytemp, ztemp, lenght_temp, &
         xavrsv, yavrsv, zavrsv
      real(kind=8), allocatable :: distmin(:), distmax(:), avgdist(:), dist(:), tmpdist(:), dh(:), tmpdh(:)
      real(kind=8), allocatable :: xloc(:), yloc(:), zloc(:), tmpxloc(:), tmpyloc(:), tmpzloc(:)
      real(kind=8) maxspac, maxspac2, dx, dy, h2, minspac, minspac2, max_dist
      real(kind=4) :: var(MAXCOL)
      real(kind=8), allocatable :: xdh(:), ydh(:), zdh(:), xdistdh(:), ydistdh(:), &
                                   zdistdh(:), lenght(:)

      integer icx, icy, ndist, na, nb, xcount, ycount, zcount, meannsamp, nsamples
      character*512 :: str
      logical testfl

      integer :: idh, idd, ixv, i, j, k, l, ndh, ndmax, idat, ndt, idf, n, j_acc
      real(kind=8) :: tmp
      integer :: idum
      real(kind=8) :: sumdist
      !outputs
      real(kind=8), intent(out) :: xyzavr, lavrsv

      ! Assign input variables to internal variables
      minspac = 2
      minspac2 = minspac*minspac
      max_dist = sqrt((maxval(xyzsv(:, 1)) - minval(xyzsv(:, 1)))**2 + &
                      (maxval(xyzsv(:, 2)) - minval(xyzsv(:, 2)))**2 + &
                      (maxval(xyzsv(:, 3)) - minval(xyzsv(:, 3)))**2)
      maxspac = max_dist
      maxspac2 = maxspac*maxspac
      !drillhole
      if (idcolsv .ne. 0) then
         idh = idcolsv
         idd = 1
      else
         idd = 0
      end if
      !Find the number of drillholes and allocate the matrix

      if (idd == 1) then
         allocate (dh(ndsv))
         dh = idsv
         call dblemodsortem(dh, ndsv, 0)
         ndh = 1
         do i = 2, ndsv
            if (dh(i) /= dh(i - 1)) ndh = ndh + 1
         end do
      end if
      allocate (xdistdh(ndh))
      allocate (ydistdh(ndh))
      allocate (zdistdh(ndh))
      allocate (lenght(ndh))

      ! Setting the number of ndmax as 2 to look the closest neighbor
      ndmax = 2
      ndist = min(ndsv, ndmax)
      allocate (dist(ndist))

      ! Read in the values

      xyzsv_idsv(:, 1) = xyzsv(:, 1)
      xyzsv_idsv(:, 2) = xyzsv(:, 2)
      xyzsv_idsv(:, 3) = xyzsv(:, 3)
      xyzsv_idsv(:, 4) = idsv
      ! Find the average distances in Z direction
      ! Grouping by DH
      nsamples = 0 ! to start the count of number of samples
      if (idd == 1) then
         call dblemodsortem(xyzsv_idsv(:, 4), ndsv, 3, xyzsv_idsv(:, 3), xyzsv_idsv(:, 2), xyzsv_idsv(:, 1))
         idh = 1
         j_acc = 0
         ndh = 1
         do i = 2, ndsv
            if (i .lt. j_acc + 1) cycle
            j = 0
            if (dh(i) /= dh(i - 1)) ndh = ndh + 1
            do n = 1, ndsv
               if (dh(n) .gt. dh(i)) exit
               if (dh(i) .eq. dh(n)) then
                  j = j + 1
               end if
            end do
            allocate (zdh(j))
            j_acc = j_acc + j
            xdh = xyzsv_idsv(j_acc - j + 1:j_acc, 1)
            ydh = xyzsv_idsv(j_acc - j + 1:j_acc, 2)
            zdh = xyzsv_idsv(j_acc - j + 1:j_acc, 3)
            ! count the number of samples per drillhole
            ! to calculate and report as average.
            nsamples = nsamples + j
            !sort the x values of this dh
            call dblemodsortem(xdh, j, 0)
            xdist = 0
            xcount = 0
            do l = 2, j
               xdist = xdist + abs(xdh(l) - xdh(l - 1))
               xcount = xcount + 1
            end do
            !sort the y values of this dh
            call dblemodsortem(ydh, j, 0)
            ydist = 0
            ycount = 0
            do l = 2, j
               ydist = ydist + abs(ydh(l) - ydh(l - 1))
               ycount = ycount + 1
            end do
            !sort the z values of this dh
            call dblemodsortem(zdh, j, 0)
            zdist = 0
            zcount = 0
            do l = 2, j
               zdist = zdist + abs(zdh(l) - zdh(l - 1))
               zcount = zcount + 1
            end do
            xdistdh(ndh) = xdist/xcount
            ydistdh(ndh) = ydist/ycount
            zdistdh(ndh) = zdist/zcount
            ! Find the lenght of this drillhole
            lenght(ndh) = sqrt((maxval(xdh) - minval(xdh))**2 + &
                               (maxval(ydh) - minval(ydh))**2 + &
                               (maxval(zdh) - minval(zdh))**2)
            deallocate (xdh, ydh, zdh)

         end do
         ! Get the average of all drillholes in x,y,z
         xtemp = 0
         ytemp = 0
         ztemp = 0
         lenght_temp = 0
         do i = 1, ndh
            xtemp = xtemp + xdistdh(i)
            ytemp = ytemp + ydistdh(i)
            ztemp = ztemp + zdistdh(i)
            lenght_temp = lenght_temp + lenght(i)
         end do
         xavrsv = xtemp/ndh
         yavrsv = ytemp/ndh
         zavrsv = ztemp/ndh
         xyzavr = sqrt(xavrsv**2 + yavrsv**2 + zavrsv**2)
         lavrsv = lenght_temp/ndh
         meannsamp = int(nsamples/ndh)
      end if

      if (debugsv) then
         write (*, *) ' Average sample spacing '
         write (*, '(1f18.3)') xyzavr
         write (*, *) ' Average drillhole lenght '
         write (*, '(1f18.3)') lavrsv
         write (*, *) ''
      end if
      return
   end subroutine samplespac_vert

end module spacing
