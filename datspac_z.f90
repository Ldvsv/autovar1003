module sp3d
   private
   public :: spac3d
contains

   subroutine spac3d(nd, xyz, hspac, vspac, nlaghorz, nlagvert, deptype, &
                     majspc, minspc, terspc)
      ! Subroutine to find the lag spacing for 3D datasets, modified from varcalc
      implicit none
      real(kind=8), parameter :: PI = 4*atan(1.0d0), DEG2RAD = PI/180d0, RAD2DEG = 180d0/PI
      real(kind=8), parameter :: BIGDBLE = 1d21, SMALLDBLE = 1d-6
      integer, intent(in) :: nd
      integer, intent(in) :: nlaghorz, nlagvert, deptype
      real(kind=8), intent(in) :: hspac, vspac
      real(kind=8), dimension(nd, 3), intent(in) :: xyz
      real(kind=8), intent(out) :: majspc, minspc, terspc
      ! internal variables
      real(kind=8), parameter :: azmtol = 22.5*DEG2RAD, diptol = 22.5*DEG2RAD
      integer :: i, j, k, l, nl
      real(kind=8) :: lagtol, azmtoldist, diptoldist, xmax, xmin, ymax, ymin, zmax, zmin, &
                      xdist, h, angle
      real(kind=8), allocatable, dimension(:) :: majdist, majnum, mindist, &
                                                 minnum, terdist, ternum, tmajdist, tmajnum, tmindist, tminnum, tterdist, tternum

      !allocate arrays for referential lag distances:

      allocate (majdist(nlaghorz + 1))
      allocate (mindist, majnum, minnum, tmajdist, tmajnum, &
                tmindist, tminnum, mold=majdist)
      allocate (terdist(nlagvert + 1))
      allocate (ternum, tterdist, tternum, mold=terdist)
      majdist = 0.0
      majnum = 0.0
      mindist = 0.0
      minnum = 0.0
      terdist = 0.0
      ternum = 0.0
      ! ----------------------------------
      ! Spacing in major direction
      ! ----------------------------------
      azmtoldist = 0.5*hspac/tan(azmtol)
      diptoldist = 0.5*vspac/tan(diptol)
      do i = 1, nd
         ! precalculate filters
         lagtol = 0.5*hspac
         if (deptype .eq. 1) then
            ymax = xyz(i, 2) + hspac
            ymin = xyz(i, 2) - hspac
            zmax = xyz(i, 3) + vspac
            zmin = xyz(i, 3) - vspac
         end if

         if (deptype .eq. 2) then
            ymax = xyz(i, 2) + 0.5*hspac
            ymin = xyz(i, 2) - 0.5*hspac
            zmax = xyz(i, 3) + 0.5*vspac
            zmin = xyz(i, 3) - 0.5*vspac
         end if
         do j = 1, nd
            ! Check horizontal bandwidths
            if (xyz(j, 2) .gt. ymax) cycle
            if (xyz(j, 2) .lt. ymin) cycle
            ! Check vertical bandwidths
            if (xyz(j, 3) .gt. zmax) cycle
            if (xyz(j, 3) .lt. zmin) cycle
            ! Check if the points are at the same point on the x' axis
            xdist = abs(xyz(j, 1) - xyz(i, 1))
            if (i .eq. j) then
               cycle
            else
               ! Check azimuth tolerance if it matters at this distance
               if (xdist .lt. azmtoldist) then
                  ! Check absolute azimuth tolerance
                  angle = abs(atan((xyz(j, 2) - xyz(i, 2))/(max(xdist, SMALLDBLE))))
                  if (angle .gt. azmtol) cycle
               end if
               ! Check dip tolerance if it matters at this distance
               if (xdist .lt. diptoldist) then
                  ! Check absolute dip tolerance
                  angle = abs(atan((xyz(j, 3) - xyz(i, 3))/(max(xdist, SMALLDBLE))))
                  if (angle .gt. diptol) cycle
               end if
            end if
            ! All constrains passed at this point
            ! Set the point on the right lag
            h = sqrt((xyz(i, 1) - xyz(j, 1))**2 + &
                     (xyz(i, 2) - xyz(j, 2))**2 + &
                     (xyz(i, 3) - xyz(j, 3))**2)
            ! Calculating the sum of wts per lag
            do l = 0, nlaghorz
               if ((h .ge. (l*hspac - lagtol)) .and. &
                   (h .le. (l*hspac + lagtol))) then
                  majnum(l + 1) = majnum(l + 1) + 1
                  majdist(l + 1) = majdist(l + 1) + xdist
               end if
            end do
         end do
      end do
      ! Left only values different than zero:
      tmajnum = majnum
      tmajdist = majdist
      deallocate (majnum, majdist)
      nl = 0
      do l = 1, nlaghorz + 1
         if (tmajnum(l) .gt. 0.0 .or. tmajdist(l) .gt. 0.0) then
            nl = nl + 1
         end if
      end do
      ! Allocate only for valid values
      allocate (majnum(nl), majdist(nl))
      !empty nl to use it again
      nl = 0
      ! Fill the valid values:
      do l = 1, nlaghorz + 1
         if (tmajnum(l) .gt. 0.0 .or. tmajdist(l) .gt. 0.0) then
            nl = nl + 1
            majnum(nl) = tmajnum(l)
            majdist(nl) = tmajdist(l)
         end if
      end do

      ! Get the distances for the lags:
      do l = 1, nl
         majdist(l) = majdist(l)/majnum(l)
      end do
      ! Get an average of lag distance
      majspc = 0.0
      do k = 1, nl - 1
         majspc = majspc + (majdist(k + 1) - majdist(k))
      end do
      ! Get final value:
      majspc = majspc/(nl - 1)
      !set nl to zero to use it again
      nl = 0

      !-----------------------------------
      ! Spacing in minor direction
      !-----------------------------------
      azmtoldist = 0.5*hspac/tan(azmtol)
      diptoldist = 0.5*vspac/tan(diptol)
      do i = 1, nd
         ! precalculate filters
         lagtol = 0.5*hspac
         if (deptype .eq. 1) then
            ymax = xyz(i, 1) + hspac
            ymin = xyz(i, 1) - hspac
            zmax = xyz(i, 3) + vspac
            zmin = xyz(i, 3) - vspac
         end if

         if (deptype .eq. 2) then
            ymax = xyz(i, 1) + 0.5*hspac
            ymin = xyz(i, 1) - 0.5*hspac
            zmax = xyz(i, 3) + 0.5*vspac
            zmin = xyz(i, 3) - 0.5*vspac
         end if
         do j = 1, nd
            ! Check horizontal bandwidths
            if (xyz(j, 1) .gt. ymax) cycle
            if (xyz(j, 1) .lt. ymin) cycle
            ! Check vertical bandwidths
            if (xyz(j, 3) .gt. zmax) cycle
            if (xyz(j, 3) .lt. zmin) cycle
            ! Check if the points are at the same point on the x' axis
            xdist = abs(xyz(j, 2) - xyz(i, 2))
            if (i .eq. j) then
               cycle
            else
               ! Check azimuth tolerance if it matters at this distance
               if (xdist .lt. azmtoldist) then
                  ! Check absolute azimuth tolerance
                  angle = abs(atan((xyz(j, 1) - xyz(i, 1))/(max(xdist, SMALLDBLE))))
                  if (angle .gt. azmtol) cycle
               end if
               ! Check dip tolerance if it matters at this distance
               if (xdist .lt. diptoldist) then
                  ! Check absolute dip tolerance
                  angle = abs(atan((xyz(j, 3) - xyz(i, 3))/(max(xdist, SMALLDBLE))))
                  if (angle .gt. diptol) cycle
               end if
            end if
            ! All constrains passed at this point
            ! Set the point on the right lag
            h = sqrt((xyz(i, 1) - xyz(j, 1))**2 + &
                     (xyz(i, 2) - xyz(j, 2))**2 + &
                     (xyz(i, 3) - xyz(j, 3))**2)
            ! Calculating the sum of wts per lag
            do l = 0, nlaghorz
               if ((h .ge. (l*hspac - lagtol)) .and. &
                   (h .le. (l*hspac + lagtol))) then
                  minnum(l + 1) = minnum(l + 1) + 1
                  mindist(l + 1) = mindist(l + 1) + xdist
               end if
            end do
         end do
      end do
      ! Left only values different than zero:
      tminnum = minnum
      tmindist = mindist
      deallocate (minnum, mindist)
      nl = 0
      do l = 1, nlaghorz + 1
         if (tminnum(l) .gt. 0.0 .or. tmindist(l) .gt. 0.0) then
            nl = nl + 1
         end if
      end do
      ! Allocate only for valid values
      allocate (minnum(nl), mindist(nl))
      !empty nl to use it again
      nl = 0
      ! Fill the valid values:
      do l = 1, nlaghorz + 1
         if (tminnum(l) .gt. 0.0 .or. tmindist(l) .gt. 0.0) then
            nl = nl + 1
            minnum(nl) = tminnum(l)
            mindist(nl) = tmindist(l)
         end if
      end do

      ! Get the distances for the lags:
      do l = 1, nl
         mindist(l) = mindist(l)/minnum(l)
      end do
      ! Get an average of lag distance
      minspc = 0.0
      do k = 1, nl - 1
         minspc = minspc + (mindist(k + 1) - mindist(k))
      end do
      ! Get final value:
      minspc = minspc/(nl - 1)

      !set nl to zero to use it again
      nl = 0

      !-----------------------------------
      ! Spacing in tertiary direction
      !-----------------------------------
      azmtoldist = 0.5*vspac/tan(azmtol)
      diptoldist = 0.5*vspac/tan(diptol)
      do i = 1, nd
         ! precalculate filters
         lagtol = 0.5*vspac
         if (deptype .eq. 1) then
            ymax = xyz(i, 1) + 0.5*hspac
            ymin = xyz(i, 1) - 0.5*hspac
            zmax = xyz(i, 2) + 4*vspac
            zmin = xyz(i, 2) - 4*vspac
         end if

         if (deptype .eq. 2) then
            ymax = xyz(i, 1) + 0.5*hspac
            ymin = xyz(i, 1) - 0.5*hspac
            zmax = xyz(i, 2) + 0.5*vspac
            zmin = xyz(i, 2) - 0.5*vspac
         end if
         do j = 1, nd
            ! Check horizontal bandwidths
            if (xyz(j, 1) .gt. ymax) cycle
            if (xyz(j, 1) .lt. ymin) cycle
            ! Check vertical bandwidths
            if (xyz(j, 2) .gt. zmax) cycle
            if (xyz(j, 2) .lt. zmin) cycle
            ! Check if the points are at the same point on the x' axis
            xdist = abs(xyz(j, 3) - xyz(i, 3))
            if (i .eq. j) then
               cycle
            else
               ! Check azimuth tolerance if it matters at this distance
               if (xdist .lt. azmtoldist) then
                  ! Check absolute azimuth tolerance
                  angle = abs(atan((xyz(j, 1) - xyz(i, 1))/(max(xdist, SMALLDBLE))))
                  if (angle .gt. azmtol) cycle
               end if
               ! Check dip tolerance if it matters at this distance
               if (xdist .lt. diptoldist) then
                  ! Check absolute dip tolerance
                  angle = abs(atan((xyz(j, 2) - xyz(i, 2))/(max(xdist, SMALLDBLE))))
                  if (angle .gt. diptol) cycle
               end if
            end if
            ! All constrains passed at this point
            ! Set the point on the right lag
            h = sqrt((xyz(i, 1) - xyz(j, 1))**2 + &
                     (xyz(i, 2) - xyz(j, 2))**2 + &
                     (xyz(i, 3) - xyz(j, 3))**2)
            ! Calculating the sum of wts per lag
            do l = 0, nlagvert
               if ((h .ge. (l*vspac - lagtol)) .and. &
                   (h .le. (l*vspac + lagtol))) then
                  ternum(l + 1) = ternum(l + 1) + 1
                  terdist(l + 1) = terdist(l + 1) + xdist
               end if
            end do
         end do
      end do
      tternum = ternum
      tterdist = terdist
      deallocate (ternum, terdist)
      nl = 0
      do l = 1, nlagvert + 1
         if (tternum(l) .gt. 0.0 .or. tterdist(l) .gt. 0.0) then
            nl = nl + 1
         end if
      end do
      ! Allocate only for valid values
      allocate (ternum(nl), terdist(nl))
      !empty nl to use it again
      nl = 0
      ! Fill the valid values:
      do l = 1, nlagvert + 1
         if (tternum(l) .gt. 0.0 .or. tterdist(l) .gt. 0.0) then
            nl = nl + 1
            ternum(nl) = tternum(l)
            terdist(nl) = tterdist(l)
         end if
      end do

      ! Get the distances for the lags:
      do l = 1, nl
         terdist(l) = terdist(l)/ternum(l)
      end do
      ! Get an average of lag distance
      terspc = 0.0
      do k = 1, nl - 1
         terspc = terspc + (terdist(k + 1) - terdist(k))
      end do
      ! Get final value:
      terspc = terspc/(nl - 1)
      !set nl to zero to use it again
      nl = 0
      return
   end subroutine spac3d

end module sp3d
