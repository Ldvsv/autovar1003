module alldir

   public ::  omnipair
contains

   ! Subroutine to calculate the omnidirectional pairs inside the cube
   ! Modified from varcalc.
   subroutine omnipair(ndv, maxdistv, tminv, tmaxv, xv, yv, zv, vv1, &
                       geoloc, wts, outputidxv, sqdfdistv, dxv1, dyv1, dzv1, &
                       sqdfv, idfv, wtsomniv, debugv)
      ! Variable declaration
      use sortem
      use gen_subs
      implicit none
      ! input variables
      integer, intent(in) :: ndv, geoloc
      real(kind=8), dimension(ndv), intent(in) :: xv, yv, zv, vv1, wts
      real(kind=8), intent(in) :: maxdistv, tminv, tmaxv
      logical, intent(in) :: debugv
      ! output variables
      real(kind=8), allocatable, intent(out) ::dxv1(:), dyv1(:), dzv1(:), sqdfv(:), sqdfdistv(:), &
                                                wtsomniv(:)
      integer, allocatable, intent(out) :: idfv(:)
      integer, intent(out) :: outputidxv
      ! internal variables
      real(kind=8) :: h, srch_rnge, xdist, ydist, zdist, ymin, zmin, ymax, zmax, &
                      farthest, azmtolrad, diptolrad, dist, maxx_dist, av, ss, num_mass, den_mass
      integer, allocatable :: introtindex(:)
      real(kind=8), allocatable :: dblerotindex(:), xyz(:, :)
      integer :: i, j, ntotal, nps, startidx
      logical :: omni
      real(kind=8), parameter :: PI = 4*atan(1.0d0)
      real(kind=8), parameter :: DEG2RAD = PI/180d0, RAD2DEG = 180d0/PI
      ! End of variable declaration

      srch_rnge = 0.2*maxdistv
      nps = 0  ! Number of pairs selected

      ! Get the average and variance of the selected points
      av = 0.0
      ss = 0.0
      do i = 1, ndv
         av = av + vv1(i)
         ss = ss + vv1(i)*vv1(i)
      end do
      av = av/max(real(ndv), 1.0)
      ss = (ss/max(real(ndv), 1.0)) - av*av

      ! Allocate in convenient way and do calculations
      allocate (xyz(ndv, 3))
      xyz(:, 1) = xv(:)
      xyz(:, 2) = yv(:)
      xyz(:, 3) = zv(:)

      ! Loop over points to calculate distances
      do i = 1, ndv - 1  ! This point, i, is the "tail" variable in GSLIB notation
         do j = 1, ndv  ! This point, j, is the "head" variable in GSLIB notation
            ! Check if the points are at the same point on the xv' axis
            xdist = xyz(j, 1) - xyz(i, 1)
            ydist = xyz(j, 2) - xyz(i, 2)
            zdist = xyz(j, 3) - xyz(i, 3)
            ! Calculate the distance between the points
            dist = sqrt(xdist**2 + ydist**2 + zdist**2)
            if (dist .le. srch_rnge) then
               nps = nps + 1
            else
               continue
            end if
         end do
      end do

      ! Note: The calculation of the pairs comes from a modification of varcalc.
      ! Allocate the needed memory
      ntotal = nps

      allocate (sqdfv(ntotal), dxv1(ntotal), dyv1(ntotal), dzv1(ntotal), idfv(ntotal), &
                sqdfdistv(ntotal), wtsomniv(ntotal))

      if (debugv) write (*, *) 'Working on selection of pairs'
      if (debugv) write (*, '(A, 1I15)') 'Total of pairs ', ntotal
      if (debugv) write (*, '(A, 11f18.3)') 'maxdistv ', maxdistv
      if (debugv) write (*, '(A, 11f18.3)') 'srch_rnge ', srch_rnge

      ! Initialize index tracking for values
      allocate (dblerotindex(ndv))
      do i = 1, ndv
         dblerotindex(i) = i
      end do

      ! Set omnidirectional search parameters
      azmtolrad = 90d0*DEG2RAD
      diptolrad = 90d0*DEG2RAD
      omni = .true.

      ! Sort along X
      call dblemodsortem(xyz(:, 1), ndv, 3, xyz(:, 2), xyz(:, 3), dblerotindex)

      ! Track the sorted indices for computing values
      allocate (introtindex(ndv))
      do i = 1, ndv
         introtindex(i) = int(dblerotindex(i))
      end do

      ! Calculate the maximum distance between pairs
      ! Restrict to selected part of the anchor:
      if (geoloc .eq. 1) then
         call get_max_dist(ndv, 3, xyz, maxx_dist)
      end if

      if (geoloc .eq. 0) then
         ! Calculate the maximum distance between pairs to only calculate points inside this range
         ! Restrict to 0.2 domain size:
         call get_max_dist(ndv, 3, xyz, maxx_dist)
         maxx_dist = 0.2*maxx_dist

      end if

      ! Main variogram loop over all points indexed by i
      outputidxv = 0

      do i = 1, ndv - 1  ! This point, i, is the "tail" variable in GSLIB notation

         ! Loop over all possible pairs forward of the current point
         startidx = 1

         do j = startidx, ndv - 1  ! This point, j, is the "head" variable in GSLIB notation

            ! Skip if the points are at the same location on the xv' axis
            if (i .eq. j) cycle

            ! Calculate distances between points
            xdist = xyz(j, 1) - xyz(i, 1)
            ydist = xyz(j, 2) - xyz(i, 2)
            zdist = xyz(j, 3) - xyz(i, 3)

            ! Calculate the distance between the points
            h = sqrt((xyz(i, 1) - xyz(j, 1))**2 + &
                     (xyz(i, 2) - xyz(j, 2))**2 + &
                     (xyz(i, 3) - xyz(j, 3))**2)

            ! Check if the distance is within the search range
            if (h .le. srch_rnge) then
               ! Inside the search radius, proceed with calculations

               ! Are both values within the allowed range?
               if (vv1(introtindex(i)) .le. tminv .or. vv1(introtindex(i)) .gt. tmaxv) cycle
               if (vv1(introtindex(j)) .le. tminv .or. vv1(introtindex(j)) .gt. tmaxv) cycle

               ! Not trimming, proceed with calculations
               outputidxv = outputidxv + 1

               ! Set the first table output
               idfv(outputidxv) = outputidxv
               sqdfdistv(outputidxv) = h
               dxv1(outputidxv) = xdist
               dyv1(outputidxv) = ydist
               dzv1(outputidxv) = zdist
               wtsomniv(outputidxv) = abs(wts(i)*wts(j))

               ! Pseudo Correlogram calculation
               sqdfv(outputidxv) = (vv1(introtindex(i))*vv1(introtindex(j)) - av**2)/ss
            else
               continue
            end if

         end do

      end do

      return

   end subroutine omnipair

end module alldir
