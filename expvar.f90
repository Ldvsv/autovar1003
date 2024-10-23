module varcalcs
   private
   public :: varexp, varexpflx
contains

   ! Subroutine to calculate the experimental correlogram and then
   ! get the experimental variogram, it is modified from the klp program.
   subroutine varexp(xev, yev, zev, vev, lrotmatev, ndev, usemajev, useminev, useterev, &
                     nmajlagev, nminlagev, nterlagev, dmajlagev, dminlagev, dterlagev, zavr_anis_in, &
                     usexdirev, useydirev, usezdirev, debugev, deptypeev, wtsvar, &
                     lvmajnumev, lvmajcovev, lvmajmuev, lvmajmuhev, lvmajsuev, lvmajsuhev, lvmajdisev, lvmajwtev, &
                     lvminnumev, lvmincovev, lvminmuev, lvminmuhev, lvminsuev, lvminsuhev, lvmindisev, lvminwtev, &
                     lvternumev, lvtercovev, lvtermuev, lvtermuhev, lvtersuev, lvtersuhev, lvterdisev, lvterwtev, &
                     anchx, anchy, anchz, anx, geoloc, strattransev)

      ! Variable declaration
      use sortem
      implicit none

      real, parameter :: EPSLON = 1.0e-6
      real(kind=8), parameter :: PI = 4*atan(1.0d0), DEG2RAD = PI/180d0, RAD2DEG = 180d0/PI
      real(kind=8), parameter :: BIGDBLE = 1d21, SMALLDBLE = 1d-6
      ! inputs
      integer, intent(in) ::  ndev, deptypeev, strattransev, anx, geoloc
      integer, intent(inout) :: nmajlagev, nminlagev, nterlagev
      real(kind=8), dimension(ndev), intent(in):: xev, yev, zev, vev, wtsvar
      real(kind=8), dimension(3, 3), intent(in) :: lrotmatev
      logical, intent(in) :: usemajev, useminev, useterev, usexdirev, useydirev, &
                             usezdirev, debugev
      real(kind=8), intent(in) :: dmajlagev, dminlagev, dterlagev, zavr_anis_in
      real(kind=8), allocatable, dimension(:) :: majorwt, minorwt, tertiwt
      real(kind=8), intent(in) :: anchx, anchy, anchz ! coordinates of current anchor
      ! outputs:
      real(kind=8), allocatable, dimension(:), intent(inout) :: &
         lvmajnumev, lvmajcovev, lvmajmuev, lvmajmuhev, &
         lvmajsuev, lvmajsuhev, lvmajdisev, lvmajwtev
      real(kind=8), allocatable, dimension(:), intent(inout) :: &
         lvminnumev, lvmincovev, lvminmuev, lvminmuhev, &
         lvminsuev, lvminsuhev, lvmindisev, lvminwtev
      real(kind=8), allocatable, dimension(:), intent(inout) :: &
         lvternumev, lvtercovev, lvtermuev, lvtermuhev, &
         lvtersuev, lvtersuhev, lvterdisev, lvterwtev
      ! temporary arrays to select outputs:
      real(kind=8), dimension(nmajlagev + 1) :: &
         tlvmajnumev, tlvmajcovev, tlvmajmuev, tlvmajmuhev, &
         tlvmajsuev, tlvmajsuhev, tlvmajdisev, tlvmajwtev
      real(kind=8), dimension(nminlagev + 1) :: &
         tlvminnumev, tlvmincovev, tlvminmuev, tlvminmuhev, &
         tlvminsuev, tlvminsuhev, tlvmindisev, tlvminwtev
      real(kind=8), dimension(nterlagev + 1) :: &
         tlvternumev, tlvtercovev, tlvtermuev, tlvtermuhev, &
         tlvtersuev, tlvtersuhev, tlvterdisev, tlvterwtev
      ! internal variables
      real(kind=8) :: xyz(ndev, 3), xyzrot(ndev, 3), wts(ndev)
      real(kind=8), dimension(ndev) :: dblerotindex
      real(kind=8) :: farthest, xmax, xmin, ymax, ymin, zmax, zmin, &
                      xdist, ydist, zdist, angle, azmtoldist, diptoldist, lagtol, h
      real(kind=8), parameter :: azmtol = 22.5*DEG2RAD, diptol = 22.5*DEG2RAD

      real(kind=8) :: dmaj, dmin, dter, dxv, dyv, dzv, tolfrac, dmajtolh, dmajtolv, &
                      dmintolh, dmintolv, dtertolh, dtertolv, reftanh, reftanv, genwtsum, genmean, genvar, &
                      maxcovmaj, zavranisev
      integer :: i, j, k, l, ncountmaj, ncountmin, ncountter

      real(kind=8) :: ama, bma, ami, bmi, ate, bte, wmaj, wmin, wter
      logical :: searchdif
      ! temporary variables for checking weights
      character(len=8) :: fmt
      character(len=100) :: output_line
      integer :: i_exp, tmp_anx, nmaj, nmin, nter! number of pairs per direction
      integer :: majdir, mindir, terdir
      ! all lags considered.
      character(256) :: ii
      real(kind=8) :: midx, midy, midz, danx
      real(kind=8) :: wt, sp
      real(kind=8), allocatable, dimension(:) :: danxs, wtss, hdist, lagns, sumwt
      real(kind=8), allocatable, dimension(:) :: accwtmaj, accwtmin, accwtter, gwts
      real(kind=8) :: varmajglc, varminglc, varterglc, daxs, sdaxs

      zavranisev = zavr_anis_in
      ! code
      !
      fmt = '(I2.2)'
      write (ii, fmt) anx

      ! Get the sp
      ! the smoothing parameter is going to be the same for all directions.
      sp = max(abs(minval(xev) - maxval(xev)), abs(minval(yev) - maxval(yev)))
      sp = sp

      ! Deallocate dummy input arrays
      deallocate (lvmajnumev, lvmajcovev, lvmajmuev, lvmajmuhev, &
                  lvmajsuev, lvmajsuhev, lvmajdisev, lvmajwtev)
      deallocate (lvternumev, lvtercovev, lvtermuev, lvtermuhev, &
                  lvtersuev, lvtersuhev, lvterdisev, lvterwtev)
      deallocate (lvminnumev, lvmincovev, lvminmuev, lvminmuhev, &
                  lvminsuev, lvminsuhev, lvmindisev, lvminwtev)
      !Allocate arrays
      allocate (lvmajnumev(nmajlagev + 1))
      allocate (lvmajcovev, lvmajmuev, lvmajmuhev, lvmajsuev, lvmajsuhev, &
                lvmajdisev, lvmajwtev, majorwt, accwtmaj, mold=lvmajnumev)
      allocate (lvminnumev(nminlagev + 1))
      allocate (lvmincovev, lvminmuev, lvminmuhev, lvminsuev, lvminsuhev, &
                lvmindisev, lvminwtev, minorwt, accwtmin, mold=lvminnumev)
      allocate (lvternumev(nterlagev + 1))
      allocate (lvtercovev, lvtermuev, lvtermuhev, lvtersuev, lvtersuhev, &
                lvterdisev, lvterwtev, tertiwt, accwtter, mold=lvternumev)
      ! Set all arrays to zero to start
      lvmajnumev = 0; lvmajcovev = 0.0; lvmajmuev = 0.0; lvmajmuhev = 0.0
      lvmajsuev = 0.0; lvmajsuhev = 0.0; lvmajdisev = 0.0; lvmajwtev = 0.0
      lvminnumev = 0; lvmincovev = 0.0; lvminmuev = 0.0; lvminmuhev = 0.0
      lvminsuev = 0.0; lvminsuhev = 0.0; lvmindisev = 0.0; lvminwtev = 0.0
      lvternumev = 0; lvtercovev = 0.0; lvtermuev = 0.0; lvtermuhev = 0.0
      lvtersuev = 0.0; lvtersuhev = 0.0; lvterdisev = 0.0; lvterwtev = 0.0
      majorwt = 0.0; minorwt = 0.0; tertiwt = 0.0
      accwtmaj = 1.0; accwtmin = 1.0; accwtter = 1.0
      ! Set temporary arrays to zero
      tlvmajnumev = 0; tlvmajcovev = 0.0; tlvmajmuev = 0.0; tlvmajmuhev = 0.0
      tlvmajsuev = 0.0; tlvmajsuhev = 0.0; tlvmajdisev = 0.0; tlvmajwtev = 0.0
      tlvminnumev = 0; tlvmincovev = 0.0; tlvminmuev = 0.0; tlvminmuhev = 0.0
      tlvminsuev = 0.0; tlvminsuhev = 0.0; tlvmindisev = 0.0; tlvminwtev = 0.0
      tlvternumev = 0; tlvtercovev = 0.0; tlvtermuev = 0.0; tlvtermuhev = 0.0
      tlvtersuev = 0.0; tlvtersuhev = 0.0; tlvterdisev = 0.0; tlvterwtev = 0.0

      ! Set tolerances
      reftanh = tan(22d0*DEG2RAD)
      reftanv = tan(5d0*DEG2RAD)
      ! If strat transformation has been done, another approach of restrictions
      ! will be set in place
      searchdif = .false.
      if (strattransev .eq. 1) searchdif = .true.
      ! Get the rotated coordinates
      xyz(:, 1) = xev
      xyz(:, 2) = yev
      xyz(:, 3) = zev
      xyzrot = matmul(xyz, lrotmatev)
      ! Assign array of weights
      wts = wtsvar
      if (geoloc .eq. 1) then
         ! Get the weights for the rotated coordinates
         allocate (gwts(ndev))
         do j = 1, ndev
            daxs = sqrt((anchx - xyzrot(j, 1))**2 + &
                        (anchy - xyzrot(j, 2))**2 + &
                        (anchz - xyzrot(j, 3))**2)
            sdaxs = daxs*daxs
            gwts(j) = exp(-(sdaxs/(2*sp**2))) + EPSLON
         end do

         ! Get the general weighted variance in case the variance for anchors are
         ! smaller than the covariance values
         genwtsum = 0.0
         do j = 1, ndev
            genwtsum = genwtsum + gwts(j)
         end do
         genmean = 0.0
         genvar = 0.0

         ! Get the weighted mean
         do j = 1, ndev
            genmean = genmean + (gwts(j)/genwtsum)*vev(j)
         end do

         ! Get the weighted variance
         do j = 1, ndev
            genvar = genvar + (gwts(j)/genwtsum)*(vev(j) - genmean)**2
         end do
      end if

      ! index tracking
      do i = 1, ndev
         dblerotindex(i) = i
      end do

      majdir = 1
      mindir = 2
      terdir = 3
      !
      !First calculate the sum of weights if anchors selected
      !
      if (geoloc .eq. 1) then
         ! Major direction:
         if (usemajev) then
            call wtpairs(xyzrot, vev, gwts, dmajlagev, nmajlagev, zavranisev, ndev, deptypeev, &
                         sp, anchx, anchy, anchz, azmtol, diptol, majdir, danxs, wtss, hdist, lagns, majorwt)
         end if

         ! Minor direction:
         if (useminev) then
            call wtpairs(xyzrot, vev, gwts, dminlagev, nminlagev, zavranisev, ndev, deptypeev, &
                         sp, anchx, anchy, anchz, azmtol, diptol, mindir, danxs, wtss, hdist, lagns, minorwt)

         end if

         ! Tertiary direction:
         if (useterev) then
            call wtpairs(xyzrot, vev, gwts, dterlagev, nterlagev, zavranisev, ndev, deptypeev, &
                         sp, anchx, anchy, anchz, azmtol, diptol, terdir, danxs, wtss, hdist, lagns, tertiwt)

         end if
      end if
      !
      ! Calculate the values for the variogram
      !

      ! Major direction:
      if (usemajev) then
         call corrc(xyzrot, vev, gwts, majorwt, dmajlagev, nmajlagev, zavranisev, ndev, deptypeev, &
                    sp, anchx, anchy, anchz, azmtol, diptol, majdir, geoloc, lvmajnumev, lvmajdisev, lvmajcovev, lvmajmuev, &
                    lvmajmuhev, lvmajsuev, lvmajsuhev, varmajglc)
      end if

      ! Minor direction:
      if (useminev) then
         call corrc(xyzrot, vev, gwts, minorwt, dminlagev, nminlagev, zavranisev, ndev, deptypeev, &
                    sp, anchx, anchy, anchz, azmtol, diptol, mindir, geoloc, lvminnumev, lvmindisev, lvmincovev, lvminmuev, &
                    lvminmuhev, lvminsuev, lvminsuhev, varminglc)
      end if

      ! Tertiary direction:
      if (useterev) then
         call corrc(xyzrot, vev, gwts, tertiwt, dterlagev, nterlagev, zavranisev, ndev, deptypeev, &
                    sp, anchx, anchy, anchz, azmtol, diptol, terdir, geoloc, lvternumev, lvterdisev, lvtercovev, lvtermuev, &
                    lvtermuhev, lvtersuev, lvtersuhev, varterglc)
      end if
      !
      ! Deallocate all temporary weight arrays for the next anchor
      !
      deallocate (majorwt, minorwt, tertiwt)
      !
      ! Calculate the correlations for each lag
      !

      if (geoloc .eq. 1) then ! Calculations for weighted variograms:
         do l = 1, nmajlagev + 1
            if (usemajev .and. lvmajnumev(l) .ge. 5) then
               lvmajdisev(l) = lvmajdisev(l)/real(lvmajnumev(l))
               lvmajcovev(l) = lvmajcovev(l) - lvmajmuev(l)*lvmajmuhev(l)

               lvmajwtev(l) = 1

            else
               lvmajnumev(l) = 0
               lvmajcovev(l) = 0.0
               lvmajdisev(l) = 0.0
               lvmajwtev(l) = 0.0
            end if
         end do
         do l = 1, nminlagev + 1
            if (useminev .and. lvminnumev(l) .ge. 5) then
               lvmindisev(l) = lvmindisev(l)/real(lvminnumev(l))
               lvmincovev(l) = lvmincovev(l) - lvminmuev(l)*lvminmuhev(l)

               lvminwtev(l) = 1
            else
               lvminnumev(l) = 0
               lvmincovev(l) = 0.0
               lvmindisev(l) = 0.0
               lvminwtev(l) = 0.0
            end if
         end do
         do l = 1, nterlagev + 1
            if (useterev .and. lvternumev(l) .ge. 5) then
               lvterdisev(l) = lvterdisev(l)/real(lvternumev(l))
               lvtercovev(l) = lvtercovev(l) - lvtermuev(l)*lvtermuhev(l)

               lvterwtev(l) = 1
            else
               lvternumev(l) = 0
               lvtercovev(l) = 0.0
               lvterdisev(l) = 0.0
               lvterwtev(l) = 0.0
            end if
         end do
      else ! Calculations for general variograms:
         do l = 1, nmajlagev + 1
            if (usemajev .and. lvmajnumev(l) .ge. 5) then
               lvmajcovev(l) = lvmajcovev(l)/real(lvmajnumev(l))
               lvmajdisev(l) = lvmajdisev(l)/real(lvmajnumev(l))
               lvmajmuev(l) = lvmajmuev(l)/real(lvmajnumev(l))
               lvmajmuhev(l) = lvmajmuhev(l)/real(lvmajnumev(l))
               lvmajsuev(l) = lvmajsuev(l)/real(lvmajnumev(l))
               lvmajsuhev(l) = lvmajsuhev(l)/real(lvmajnumev(l))
               lvmajcovev(l) = lvmajcovev(l) - lvmajmuev(l)*lvmajmuhev(l)
               lvmajsuev(l) = lvmajsuev(l) - lvmajmuev(l)*lvmajmuev(l)
               lvmajsuhev(l) = lvmajsuhev(l) - lvmajmuhev(l)*lvmajmuhev(l)
               if (lvmajsuev(l) .lt. EPSLON) lvmajsuev(l) = EPSLON
               if (lvmajsuhev(l) .lt. EPSLON) lvmajsuhev(l) = EPSLON
               lvmajcovev(l) = lvmajcovev(l)/sqrt(lvmajsuev(l)*lvmajsuhev(l))

               lvmajwtev(l) = 1

            else
               lvmajnumev(l) = 0
               lvmajcovev(l) = 0.0
               lvmajdisev(l) = 0.0
               lvmajwtev(l) = 0.0
            end if
         end do
         do l = 1, nminlagev + 1
            if (useminev .and. lvminnumev(l) .ge. 5) then
               lvmincovev(l) = lvmincovev(l)/real(lvminnumev(l))
               lvmindisev(l) = lvmindisev(l)/real(lvminnumev(l))
               lvminmuev(l) = lvminmuev(l)/real(lvminnumev(l))
               lvminmuhev(l) = lvminmuhev(l)/real(lvminnumev(l))
               lvminsuev(l) = lvminsuev(l)/real(lvminnumev(l))
               lvminsuhev(l) = lvminsuhev(l)/real(lvminnumev(l))
               lvmincovev(l) = lvmincovev(l) - lvminmuev(l)*lvminmuhev(l)
               lvminsuev(l) = lvminsuev(l) - lvminmuev(l)*lvminmuev(l)
               lvminsuhev(l) = lvminsuhev(l) - lvminmuhev(l)*lvminmuhev(l)
               if (lvminsuev(l) .lt. EPSLON) lvminsuev(l) = EPSLON
               if (lvminsuhev(l) .lt. EPSLON) lvminsuhev(l) = EPSLON
               lvmincovev(l) = lvmincovev(l)/sqrt(lvminsuev(l)*lvminsuhev(l))

               lvminwtev(l) = 1
            else
               lvminnumev(l) = 0
               lvmincovev(l) = 0.0
               lvmindisev(l) = 0.0
               lvminwtev(l) = 0.0
            end if
         end do
         do l = 1, nterlagev + 1
            if (useterev .and. lvternumev(l) .ge. 5) then
               lvtercovev(l) = lvtercovev(l)/real(lvternumev(l))
               lvterdisev(l) = lvterdisev(l)/real(lvternumev(l))
               lvtermuev(l) = lvtermuev(l)/real(lvternumev(l))
               lvtermuhev(l) = lvtermuhev(l)/real(lvternumev(l))
               lvtersuev(l) = lvtersuev(l)/real(lvternumev(l))
               lvtersuhev(l) = lvtersuhev(l)/real(lvternumev(l))
               lvtercovev(l) = lvtercovev(l) - lvtermuev(l)*lvtermuhev(l)
               lvtersuev(l) = lvtersuev(l) - lvtermuev(l)*lvtermuev(l)
               lvtersuhev(l) = lvtersuhev(l) - lvtermuhev(l)*lvtermuhev(l)
               if (lvtersuev(l) .lt. EPSLON) lvtersuev(l) = EPSLON
               if (lvtersuhev(l) .lt. EPSLON) lvtersuhev(l) = EPSLON
               lvtercovev(l) = lvtercovev(l)/sqrt(lvtersuev(l)*lvtersuhev(l))

               lvterwtev(l) = 1
            else
               lvternumev(l) = 0
               lvtercovev(l) = 0.0
               lvterdisev(l) = 0.0
               lvterwtev(l) = 0.0
            end if
         end do

      end if

      if (geoloc .eq. 1) then ! Calculations for weighted Covariance:
         do l = 1, nmajlagev + 1
            maxcovmaj = maxval(lvmajcovev)
            if (maxcovmaj .gt. varmajglc) varmajglc = genvar
            lvmajcovev(l) = (varmajglc - lvmajcovev(l))/varmajglc
         end do
         do l = 1, nminlagev + 1
            lvmincovev(l) = (varminglc - lvmincovev(l))/varminglc
         end do
         do l = 1, nterlagev + 1
            lvtercovev(l) = (varterglc - lvtercovev(l))/varterglc
         end do
      else ! if not calculating geoloc. var, we are getting the correlogram:
         do l = 1, nmajlagev + 1
            lvmajcovev(l) = 1 - lvmajcovev(l)
         end do
         do l = 1, nminlagev + 1
            lvmincovev(l) = 1 - lvmincovev(l)
         end do
         do l = 1, nterlagev + 1
            lvtercovev(l) = 1 - lvtercovev(l)
         end do
      end if

      !
      ! Remove where the lag distance is zero
      !
      ncountmaj = 0; ncountmin = 0; ncountter = 0
      do l = 1, nmajlagev + 1
         if (lvmajdisev(l) .gt. 0.0) ncountmaj = ncountmaj + 1
      end do
      do l = 1, nminlagev + 1
         if (lvmindisev(l) .gt. 0.0) ncountmin = ncountmin + 1
      end do
      do l = 1, nterlagev + 1
         if (useterev) then
            if (lvterdisev(l) .gt. 0.0) ncountter = ncountter + 1
         else
            ncountter = 1
         end if

      end do
      ! Allocate temporary arrays to deallocate outputs array
      tlvmajnumev = lvmajnumev; tlvmajcovev = lvmajcovev
      tlvmajmuev = lvmajmuev; tlvmajmuhev = lvmajmuhev
      tlvmajsuev = lvmajsuev; tlvmajsuhev = lvmajsuhev
      tlvmajdisev = lvmajdisev; tlvmajwtev = lvmajwtev

      tlvminnumev = lvminnumev; tlvmincovev = lvmincovev
      tlvminmuev = lvminmuev; tlvminmuhev = lvminmuhev
      tlvminsuev = lvminsuev; tlvminsuhev = lvminsuhev
      tlvmindisev = lvmindisev; tlvminwtev = lvminwtev

      tlvternumev = lvternumev; tlvtercovev = lvtercovev
      tlvtermuev = lvtermuev; tlvtermuhev = lvtermuhev
      tlvtersuev = lvtersuev; tlvtersuhev = lvtersuhev
      tlvterdisev = lvterdisev; tlvterwtev = lvterwtev

      ! Deallocate output arrays
      deallocate (lvmajnumev, lvmajcovev, lvmajmuev, lvmajmuhev, &
                  lvmajsuev, lvmajsuhev, lvmajdisev, lvmajwtev)
      deallocate (lvternumev, lvtercovev, lvtermuev, lvtermuhev, &
                  lvtersuev, lvtersuhev, lvterdisev, lvterwtev)
      deallocate (lvminnumev, lvmincovev, lvminmuev, lvminmuhev, &
                  lvminsuev, lvminsuhev, lvmindisev, lvminwtev)

      ! Allocate for valid values only
      allocate (lvmajnumev(ncountmaj))
      allocate (lvmajcovev, lvmajmuev, lvmajmuhev, lvmajsuev, &
                lvmajsuhev, lvmajdisev, lvmajwtev, mold=lvmajnumev)
      allocate (lvminnumev(ncountmin))
      allocate (lvmincovev, lvminmuev, lvminmuhev, lvminsuev, &
                lvminsuhev, lvmindisev, lvminwtev, mold=lvminnumev)
      allocate (lvternumev(ncountter))
      allocate (lvtercovev, lvtermuev, lvtermuhev, lvtersuev, &
                lvtersuhev, lvterdisev, lvterwtev, mold=lvternumev)

      ! fill the arrays
      j = 0
      do i = 1, nmajlagev + 1
         if (tlvmajdisev(i) .gt. 0.0) then
            j = j + 1
            lvmajnumev(j) = tlvmajnumev(i)
            lvmajcovev(j) = tlvmajcovev(i)
            lvmajmuev(j) = tlvmajmuev(i)
            lvmajmuhev(j) = tlvmajmuhev(i)
            lvmajsuev(j) = tlvmajsuev(i)
            lvmajsuhev(j) = tlvmajsuhev(i)
            lvmajdisev(j) = tlvmajdisev(i)
            lvmajwtev(j) = tlvmajwtev(i)
         end if
      end do
      j = 0
      do i = 1, nminlagev + 1
         if (tlvmindisev(i) .gt. 0.0) then
            j = j + 1
            lvminnumev(j) = tlvminnumev(i)
            lvmincovev(j) = tlvmincovev(i)
            lvminmuev(j) = tlvminmuev(i)
            lvminmuhev(j) = tlvminmuhev(i)
            lvminsuev(j) = tlvminsuev(i)
            lvminsuhev(j) = tlvminsuhev(i)
            lvmindisev(j) = tlvmindisev(i)
            lvminwtev(j) = tlvminwtev(i)
         end if
      end do
      j = 0
      if (useterev) then
         do i = 1, nterlagev + 1
            if (tlvterdisev(i) .gt. 0.0) then
               j = j + 1
               lvternumev(j) = tlvternumev(i)
               lvtercovev(j) = tlvtercovev(i)
               lvtermuev(j) = tlvtermuev(i)
               lvtermuhev(j) = tlvtermuhev(i)
               lvtersuev(j) = tlvtersuev(i)
               lvtersuhev(j) = tlvtersuhev(i)
               lvterdisev(j) = tlvterdisev(i)
               lvterwtev(j) = tlvterwtev(i)
            end if
         end do
      end if

      ! update number of lags for output

      nmajlagev = ncountmaj
      nminlagev = ncountmin
      if (useterev) nterlagev = ncountter

      !
      ! Check first point and set the rigth weight: if the variogram value is
      ! greater than the second point and if the number of pairs is less
      !
      i = 0; j = 0; k = 0
      if (lvmajcovev(1) .gt. lvmajcovev(2) .and. lvmajnumev(1) .lt. lvmajnumev(2)) &
         lvmajwtev(1) = 0.0
      if (lvmincovev(1) .gt. lvmincovev(2) .and. lvminnumev(1) .lt. lvminnumev(2)) &
         lvminwtev(1) = 0.0
      if (useterev) then
         if (lvtercovev(1) .gt. lvtercovev(2) .and. lvternumev(1) .lt. lvternumev(2)) &
            lvterwtev(1) = 0.0
      end if

      ! Checking only for 2D datasets
      if (.not. useterev) then
         do i = nminlagev, 2, -1
            if (lvmincovev(i) .lt. 1. .and. lvmindisev(i) .gt. 0.5*maxval(lvmindisev)) then
               if (lvmincovev(i) .lt. lvmincovev(i - 1)) lvminwtev(i) = 0.
            end if
         end do
      end if
      ! Checking only for 3D datasets.
      if (useterev) then
         do i = nterlagev, 2, -1
            if (lvtercovev(i) .lt. 1. .and. lvterdisev(i) .gt. 0.5*maxval(lvterdisev)) then
               if (lvtercovev(i) .lt. lvtercovev(i - 1)) lvterwtev(i) = 0.
            end if
         end do
      end if
      ! Show results on screen
      if (debugev) then
         write (*, *) ' Inferred experimental variogram: '
         if (usemajev) then
            write (*, *)
            write (*, *) ' Major variogram '
            do l = 1, nmajlagev
               write (*, 120) lvmajnumev(l), lvmajdisev(l), lvmajcovev(l)
            end do
         end if
         if (useminev) then
            write (*, *)
            write (*, *) ' Minor variogram '
            do l = 1, nminlagev
               write (*, 120) lvminnumev(l), lvmindisev(l), lvmincovev(l)
            end do
         end if
         if (useterev) then
            write (*, *)
            write (*, *) ' Tertiary variogram '
            do l = 1, nterlagev
               write (*, 120) lvternumev(l), lvterdisev(l), lvtercovev(l)
            end do
         end if
      end if
120   format(4x, 4f12.4, 2x, 4f12.4)
      return
   end subroutine varexp

   subroutine varexpflx(xev, yev, zev, vev, lrotmatev, ndev, usemajev, useminev, useterev, &
                        nmajlagev, nminlagev, nterlagev, dmajlagev, dminlagev, dterlagev, zavr_anis, &
                        usexdirev, useydirev, usezdirev, debugev, deptypeev, wtsvar, &
                        lvmajnumev, lvmajcovev, lvmajmuev, lvmajmuhev, lvmajsuev, lvmajsuhev, lvmajdisev, lvmajwtev, &
                        lvminnumev, lvmincovev, lvminmuev, lvminmuhev, lvminsuev, lvminsuhev, lvmindisev, lvminwtev, &
                        lvternumev, lvtercovev, lvtermuev, lvtermuhev, lvtersuev, lvtersuhev, lvterdisev, lvterwtev, &
                        anchx, anchy, anchz, anx, geoloc, strattransev)

      ! Variable declaration
      use sortem
      implicit none

      real, parameter :: EPSLON = 1.0e-6
      real(kind=8), parameter :: PI = 4*atan(1.0d0), DEG2RAD = PI/180d0, RAD2DEG = 180d0/PI
      real(kind=8), parameter :: BIGDBLE = 1d21, SMALLDBLE = 1d-6
      ! inputs
      integer, intent(in) ::  ndev, deptypeev, strattransev, anx, geoloc
      integer, intent(inout) :: nmajlagev, nminlagev, nterlagev
      real(kind=8), dimension(ndev), intent(in):: xev, yev, zev, vev, wtsvar
      real(kind=8), dimension(3, 3), intent(in) :: lrotmatev
      logical, intent(in) :: usemajev, useminev, useterev, usexdirev, useydirev, &
                             usezdirev, debugev
      real(kind=8), intent(in) :: dmajlagev, dminlagev, dterlagev, zavr_anis
      real(kind=8) :: zavranisev
      real(kind=8), allocatable, dimension(:) :: majwt, minwt, terwt
      real(kind=8), intent(in) :: anchx, anchy, anchz ! coordinates of current anchor
      ! outputs:
      real(kind=8), allocatable, dimension(:), intent(inout) :: &
         lvmajnumev, lvmajcovev, lvmajmuev, lvmajmuhev, &
         lvmajsuev, lvmajsuhev, lvmajdisev, lvmajwtev
      real(kind=8), allocatable, dimension(:), intent(inout) :: &
         lvminnumev, lvmincovev, lvminmuev, lvminmuhev, &
         lvminsuev, lvminsuhev, lvmindisev, lvminwtev
      real(kind=8), allocatable, dimension(:), intent(inout) :: &
         lvternumev, lvtercovev, lvtermuev, lvtermuhev, &
         lvtersuev, lvtersuhev, lvterdisev, lvterwtev
      ! temporary arrays to select outputs:
      real(kind=8), dimension(nmajlagev) :: &
         tlvmajnumev, tlvmajcovev, tlvmajmuev, tlvmajmuhev, &
         tlvmajsuev, tlvmajsuhev, tlvmajdisev, tlvmajwtev
      real(kind=8), dimension(nminlagev) :: &
         tlvminnumev, tlvmincovev, tlvminmuev, tlvminmuhev, &
         tlvminsuev, tlvminsuhev, tlvmindisev, tlvminwtev
      real(kind=8), dimension(nterlagev) :: &
         tlvternumev, tlvtercovev, tlvtermuev, tlvtermuhev, &
         tlvtersuev, tlvtersuhev, tlvterdisev, tlvterwtev
      ! internal variables
      real(kind=8) :: xyz(ndev, 3), xyzrot(ndev, 3), wts(ndev)
      real(kind=8), dimension(ndev) :: dblerotindex
      real(kind=8), allocatable, dimension(:) :: gwts, wmajacc, wminacc, wteracc
      real(kind=8) :: farthest, xmax, xmin, ymax, ymin, zmax, zmin, &
                      xdist, ydist, zdist, angle, azmtoldist, diptoldist, lagtol, h
      real(kind=8), parameter :: azmtol = 22.5*DEG2RAD, diptol = 22.5*DEG2RAD

      real(kind=8) :: dmaj, dmin, dter, dxv, dyv, dzv, tolfrac, dmajtolh, dmajtolv, &
                      dmintolh, dmintolv, dtertolh, dtertolv, reftanh, reftanv, sdaxs, daxs
      integer :: i, j, k, l, ncountmaj, ncountmin, ncountter

      real(kind=8) :: ama, bma, ami, bmi, ate, bte, wmaj, wmin, wter
      logical :: searchdif
      ! temporary variables for checking weights
      character(len=8) :: fmt
      character(len=100) :: output_line
      integer :: i_exp, tmp_anx, nmaj, nmin, nter! number of pairs per direction
      integer :: majdir, mindir, terdir
      ! all lags considered.
      character(256) :: ii
      real(kind=8) :: midx, midy, midz, danx
      real(kind=8) :: wt, sp ! to calculate Gaussian weighting for each pair
      real(kind=8), allocatable, dimension(:) :: danxs, wtss, hdist, lagns, sumwt
      real(kind=8), allocatable, dimension(:) :: accwtmaj, accwtmin, accwtter
      !
      ! code
      !
      zavranisev = zavr_anis
      ! Specify the outputs file name
      fmt = '(I2.2)'
      write (ii, fmt) anx

      ! Get the sp
      sp = min(abs(minval(xev) - maxval(xev)), abs(minval(yev) - maxval(yev)))
      sp = sp

      ! Deallocate dummy input arrays
      deallocate (lvmajnumev, lvmajcovev, lvmajmuev, lvmajmuhev, &
                  lvmajsuev, lvmajsuhev, lvmajdisev, lvmajwtev)
      deallocate (lvternumev, lvtercovev, lvtermuev, lvtermuhev, &
                  lvtersuev, lvtersuhev, lvterdisev, lvterwtev)
      deallocate (lvminnumev, lvmincovev, lvminmuev, lvminmuhev, &
                  lvminsuev, lvminsuhev, lvmindisev, lvminwtev)
      !Allocate arrays
      allocate (lvmajnumev(nmajlagev + 1))
      allocate (lvmajcovev, lvmajmuev, lvmajmuhev, lvmajsuev, lvmajsuhev, &
                lvmajdisev, lvmajwtev, majwt, accwtmaj, mold=lvmajnumev)
      allocate (lvminnumev(nminlagev + 1))
      allocate (lvmincovev, lvminmuev, lvminmuhev, lvminsuev, lvminsuhev, &
                lvmindisev, lvminwtev, minwt, accwtmin, mold=lvminnumev)
      allocate (lvternumev(nterlagev + 1))
      allocate (lvtercovev, lvtermuev, lvtermuhev, lvtersuev, lvtersuhev, &
                lvterdisev, lvterwtev, terwt, accwtter, mold=lvternumev)
      ! set all to zero to start
      lvmajnumev = 0; lvmajcovev = 0.0; lvmajmuev = 0.0; lvmajmuhev = 0.0
      lvmajsuev = 0.0; lvmajsuhev = 0.0; lvmajdisev = 0.0; lvmajwtev = 0.0
      lvminnumev = 0; lvmincovev = 0.0; lvminmuev = 0.0; lvminmuhev = 0.0
      lvminsuev = 0.0; lvminsuhev = 0.0; lvmindisev = 0.0; lvminwtev = 0.0
      lvternumev = 0; lvtercovev = 0.0; lvtermuev = 0.0; lvtermuhev = 0.0
      lvtersuev = 0.0; lvtersuhev = 0.0; lvterdisev = 0.0; lvterwtev = 0.0
      majwt = 0.0; minwt = 0.0; terwt = 0.0
      accwtmaj = 1.0; accwtmin = 1.0; accwtter = 1.0
      ! Set tolerances
      reftanh = tan(22d0*DEG2RAD)
      reftanv = tan(5d0*DEG2RAD)
      !Set same weights
      wmaj = 1; wmin = 1; wter = 1

      if (geoloc .eq. 0) then
         searchdif = .false.
         if (strattransev .eq. 1) searchdif = .true.
         !Setting the pair look up and calculations.
         dxv = 0.0; dyv = 0.0; dzv = 0.0; tolfrac = 0.25
         do i = 1, ndev
            do j = 1, ndev
               if (usexdirev) dxv = xev(i) - xev(j)
               if (useydirev) dyv = yev(i) - yev(j)
               if (usezdirev) dzv = zev(i) - zev(j)
               ! Getting the components using the rotation matrix
               dmaj = abs(lrotmatev(1, 1)*dxv + lrotmatev(1, 2)*dyv + lrotmatev(1, 3)*dzv)
               dmin = abs(lrotmatev(2, 1)*dxv + lrotmatev(2, 2)*dyv + lrotmatev(2, 3)*dzv)
               dter = abs(lrotmatev(3, 1)*dxv + lrotmatev(3, 2)*dyv + lrotmatev(3, 3)*dzv)
               wmaj = 1
               wmin = 1
               wter = 1
               !
               !     Major direction variogram
               !
               if (usemajev .and. i .ne. j) then
                  l = 1 + IDNINT(dmaj/dmajlagev)
                  !set tolerances based on deposit type
                  if (deptypeev .eq. 1) dmajtolh = 2*tolfrac*dmaj
                  if (deptypeev .eq. 1) dmajtolv = 4*zavranisev
                  if (deptypeev .eq. 2) dmajtolh = dmajlagev
                  if (deptypeev .eq. 2) dmajtolv = 1.5*zavranisev

                  if (deptypeev .eq. 2 .and. searchdif) dmajtolh = tolfrac*dmaj
                  if (deptypeev .eq. 2 .and. searchdif) dmajtolv = 1.5*zavranisev
                  if (l .le. nmajlagev) then
                     if (useminev .and. dmin .gt. dmajtolh) go to 220
                     if (useterev .and. dter .gt. dmajtolv) go to 220
                     ! To consider az tol of 22 and dip tol of 5
                     if (deptypeev .eq. 2) then
                        if (reftanh .lt. (dmin/dmaj)) go to 220
                        if (reftanv .lt. (dter/dmaj)) go to 220
                     end if
                     lvmajnumev(l) = lvmajnumev(l) + 1
                     lvmajwtev(l) = lvmajwtev(l) + wmaj
                     lvmajdisev(l) = lvmajdisev(l) + dmaj
                     lvmajcovev(l) = lvmajcovev(l) + vev(i)*vev(j)*wmaj
                     lvmajmuev(l) = lvmajmuev(l) + vev(i)*wmaj
                     lvmajmuhev(l) = lvmajmuhev(l) + vev(j)*wmaj
                     lvmajsuev(l) = lvmajsuev(l) + vev(i)*vev(i)*wmaj
                     lvmajsuhev(l) = lvmajsuhev(l) + vev(j)*vev(j)*wmaj
220                  continue
                  end if
               end if
               !
               !     Minor direction variogram
               !
               if (useminev .and. i .ne. j) then
                  l = 1 + IDNINT(dmin/dminlagev)
                  if (deptypeev .eq. 1) dmintolh = 2*tolfrac*dmin
                  if (deptypeev .eq. 1) dmintolv = 4*zavranisev
                  if (deptypeev .eq. 2) dmintolh = dminlagev
                  if (deptypeev .eq. 2) dmintolv = 1.5*zavranisev

                  if (deptypeev .eq. 2 .and. searchdif) dmintolh = tolfrac*dmin
                  if (deptypeev .eq. 2 .and. searchdif) dmintolv = 1.5*zavranisev
                  if (l .le. nminlagev) then
                     if (usemajev .and. dmaj .gt. dmintolh) go to 221
                     if (useterev .and. dter .gt. dmintolv) go to 221
                     ! To consider az tol of 22 and dip tol of 5
                     if (deptypeev .eq. 2) then
                        if (reftanh .lt. (dmaj/dmin)) go to 221
                        if (reftanv .lt. (dter/dmin)) go to 221
                     end if
                     lvminnumev(l) = lvminnumev(l) + 1
                     lvminwtev(l) = lvminwtev(l) + wmin
                     lvmindisev(l) = lvmindisev(l) + dmin
                     lvmincovev(l) = lvmincovev(l) + vev(i)*vev(j)*wmin
                     lvminmuev(l) = lvminmuev(l) + vev(i)*wmin
                     lvminmuhev(l) = lvminmuhev(l) + vev(j)*wmin
                     lvminsuev(l) = lvminsuev(l) + vev(i)*vev(i)*wmin
                     lvminsuhev(l) = lvminsuhev(l) + vev(j)*vev(j)*wmin
221                  continue
                  end if
               end if
               !
               !     Tertiary direction variogram
               !
               if (useterev .and. i .ne. j) then
                  l = 1 + IDNINT(dter/dterlagev)
                  if (deptypeev .eq. 1) dtertolh = 8*tolfrac*dterlagev
                  if (deptypeev .eq. 1) dtertolv = 4*zavranisev
                  if (deptypeev .eq. 2) dtertolh = 1.5*tolfrac*dterlagev
                  if (deptypeev .eq. 2) dtertolv = 3*zavranisev

                  if (deptypeev .eq. 2 .and. searchdif) dtertolh = 0.5*tolfrac*dterlagev
                  if (deptypeev .eq. 2 .and. searchdif) dtertolv = 0.5*zavranisev
                  if (l .le. nterlagev) then
                     if (usemajev .and. dmaj .gt. dtertolh) go to 222
                     if (useminev .and. dmin .gt. dtertolv) go to 222
                     lvternumev(l) = lvternumev(l) + 1
                     lvterwtev(l) = lvterwtev(l) + wter
                     lvterdisev(l) = lvterdisev(l) + dter
                     lvtercovev(l) = lvtercovev(l) + vev(i)*vev(j)*wter
                     lvtermuev(l) = lvtermuev(l) + vev(i)*wter
                     lvtermuhev(l) = lvtermuhev(l) + vev(j)*wter
                     lvtersuev(l) = lvtersuev(l) + vev(i)*vev(i)*wter
                     lvtersuhev(l) = lvtersuhev(l) + vev(j)*vev(j)*wter
222                  continue
                  end if
               end if
            end do
         end do

         do l = 1, nmajlagev + 1
            ! Calculations for general correlogram:
            if (usemajev .and. lvmajnumev(l) .ge. 5) then
               lvmajcovev(l) = lvmajcovev(l)/real(lvmajnumev(l))
               lvmajdisev(l) = lvmajdisev(l)/real(lvmajnumev(l))
               lvmajmuev(l) = lvmajmuev(l)/real(lvmajnumev(l))
               lvmajmuhev(l) = lvmajmuhev(l)/real(lvmajnumev(l))
               lvmajsuev(l) = lvmajsuev(l)/real(lvmajnumev(l))
               lvmajsuhev(l) = lvmajsuhev(l)/real(lvmajnumev(l))
               lvmajcovev(l) = lvmajcovev(l) - lvmajmuev(l)*lvmajmuhev(l)
               lvmajsuev(l) = lvmajsuev(l) - lvmajmuev(l)*lvmajmuev(l)
               lvmajsuhev(l) = lvmajsuhev(l) - lvmajmuhev(l)*lvmajmuhev(l)
               if (lvmajsuev(l) .lt. EPSLON) lvmajsuev(l) = EPSLON
               if (lvmajsuhev(l) .lt. EPSLON) lvmajsuhev(l) = EPSLON
               lvmajcovev(l) = lvmajcovev(l)/sqrt(lvmajsuev(l)*lvmajsuhev(l))
               lvmajwtev(l) = 1

            else
               lvmajnumev(l) = 0
               lvmajcovev(l) = 0.0
               lvmajdisev(l) = 0.0
               lvmajwtev(l) = 0.0
            end if
         end do
         do l = 1, nminlagev + 1
            if (useminev .and. lvminnumev(l) .ge. 5) then
               lvmincovev(l) = lvmincovev(l)/real(lvminnumev(l))
               lvmindisev(l) = lvmindisev(l)/real(lvminnumev(l))
               lvminmuev(l) = lvminmuev(l)/real(lvminnumev(l))
               lvminmuhev(l) = lvminmuhev(l)/real(lvminnumev(l))
               lvminsuev(l) = lvminsuev(l)/real(lvminnumev(l))
               lvminsuhev(l) = lvminsuhev(l)/real(lvminnumev(l))
               lvmincovev(l) = lvmincovev(l) - lvminmuev(l)*lvminmuhev(l)
               lvminsuev(l) = lvminsuev(l) - lvminmuev(l)*lvminmuev(l)
               lvminsuhev(l) = lvminsuhev(l) - lvminmuhev(l)*lvminmuhev(l)
               if (lvminsuev(l) .lt. EPSLON) lvminsuev(l) = EPSLON
               if (lvminsuhev(l) .lt. EPSLON) lvminsuhev(l) = EPSLON
               lvmincovev(l) = lvmincovev(l)/sqrt(lvminsuev(l)*lvminsuhev(l))
               lvminwtev(l) = 1
            else
               lvminnumev(l) = 0
               lvmincovev(l) = 0.0
               lvmindisev(l) = 0.0
               lvminwtev(l) = 0.0
            end if
         end do

         do l = 1, nterlagev + 1
            ! Calculations for general correlogram:
            if (useterev .and. lvternumev(l) .ge. 5) then
               lvtercovev(l) = lvtercovev(l)/real(lvternumev(l))
               lvterdisev(l) = lvterdisev(l)/real(lvternumev(l))
               lvtermuev(l) = lvtermuev(l)/real(lvternumev(l))
               lvtermuhev(l) = lvtermuhev(l)/real(lvternumev(l))
               lvtersuev(l) = lvtersuev(l)/real(lvternumev(l))
               lvtersuhev(l) = lvtersuhev(l)/real(lvternumev(l))
               lvtercovev(l) = lvtercovev(l) - lvtermuev(l)*lvtermuhev(l)
               lvtersuev(l) = lvtersuev(l) - lvtermuev(l)*lvtermuev(l)
               lvtersuhev(l) = lvtersuhev(l) - lvtermuhev(l)*lvtermuhev(l)
               if (lvtersuev(l) .lt. EPSLON) lvtersuev(l) = EPSLON
               if (lvtersuhev(l) .lt. EPSLON) lvtersuhev(l) = EPSLON
               lvtercovev(l) = lvtercovev(l)/sqrt(lvtersuev(l)*lvtersuhev(l))
               lvterwtev(l) = 1
            else
               lvternumev(l) = 0
               lvtercovev(l) = 0.0
               lvterdisev(l) = 0.0
               lvterwtev(l) = 0.0
            end if
         end do

         ! Convert the correlogram to variogram values, since it is already
         ! normalized: var = 1 - corr
         do l = 1, nmajlagev + 1
            lvmajcovev(l) = 1 - lvmajcovev(l)
         end do
         do l = 1, nminlagev + 1
            lvmincovev(l) = 1 - lvmincovev(l)
         end do
         do l = 1, nterlagev + 1
            lvtercovev(l) = 1 - lvtercovev(l)
         end do
         !
      end if

      if (geoloc .eq. 1) then
         ! Get the weights for the rotated coordinates
         allocate (gwts(ndev))
         gwts = 0.0
         do j = 1, ndev
            daxs = sqrt((anchx - xev(j))**2 + &
                        (anchy - yev(j))**2 + &
                        (anchz - zev(j))**2)
            sdaxs = daxs*daxs
            gwts(j) = exp(-(sdaxs/(2*sp**2))) + EPSLON
         end do
         allocate (wmajacc(nmajlagev + 1))
         allocate (wminacc(nminlagev + 1))
         if (usezdirev) allocate (wteracc(nterlagev + 1))
         wmajacc = 0.0
         wminacc = 0.0
         if (usezdirev) wteracc = 0.0
         searchdif = .false.
         if (strattransev .eq. 1) searchdif = .true.
         !Setting the pair look up and calculations.
         dxv = 0.0; dyv = 0.0; dzv = 0.0; tolfrac = 0.25
         ! Get the weights for each lag.
         do i = 1, ndev
            do j = 1, ndev
               if (usexdirev) dxv = xev(i) - xev(j)
               if (useydirev) dyv = yev(i) - yev(j)
               if (usezdirev) dzv = zev(i) - zev(j)
               ! Getting the components using the rotation matrix
               dmaj = abs(lrotmatev(1, 1)*dxv + lrotmatev(1, 2)*dyv + lrotmatev(1, 3)*dzv)
               dmin = abs(lrotmatev(2, 1)*dxv + lrotmatev(2, 2)*dyv + lrotmatev(2, 3)*dzv)
               dter = abs(lrotmatev(3, 1)*dxv + lrotmatev(3, 2)*dyv + lrotmatev(3, 3)*dzv)
               !
               !     Major direction variogram
               !
               if (usemajev .and. i .ne. j) then
                  l = 1 + IDNINT(dmaj/dmajlagev)
                  if (deptypeev .eq. 1) dmajtolh = 2*tolfrac*dmaj
                  if (deptypeev .eq. 1) dmajtolv = 4*zavranisev
                  if (deptypeev .eq. 2) dmajtolh = dmajlagev
                  if (deptypeev .eq. 2) dmajtolv = 1.5*zavranisev

                  if (deptypeev .eq. 2 .and. searchdif) dmajtolh = tolfrac*dmaj
                  if (deptypeev .eq. 2 .and. searchdif) dmajtolv = 1.5*zavranisev
                  if (l .le. nmajlagev) then
                     if (useminev .and. dmin .gt. dmajtolh) go to 230
                     if (useterev .and. dter .gt. dmajtolv) go to 230
                     ! To consider az tol of 22 and dip tol of 5
                     if (deptypeev .eq. 2) then
                        if (reftanh .lt. (dmin/dmaj)) go to 230
                        if (reftanv .lt. (dter/dmaj)) go to 230
                     end if
                     wmaj = sqrt(gwts(i)*gwts(j))
                     wmajacc(l) = wmajacc(l) + wmaj
230                  continue
                  end if

               end if
               !
               !     Minor direction variogram
               !
               if (useminev .and. i .ne. j) then
                  l = 1 + IDNINT(dmin/dminlagev)
                  if (deptypeev .eq. 1) dmintolh = 2*tolfrac*dmin
                  if (deptypeev .eq. 1) dmintolv = 4*zavranisev
                  if (deptypeev .eq. 2) dmintolh = dminlagev
                  if (deptypeev .eq. 2) dmintolv = 1.5*zavranisev

                  if (deptypeev .eq. 2 .and. searchdif) dmintolh = tolfrac*dmin
                  if (deptypeev .eq. 2 .and. searchdif) dmintolv = 1.5*zavranisev

                  if (l .le. nminlagev) then
                     if (usemajev .and. dmaj .gt. dmintolh) go to 231
                     if (useterev .and. dter .gt. dmintolv) go to 231
                     ! To consider az tol of 22 and dip tol of 5
                     if (deptypeev .eq. 2) then
                        if (reftanh .lt. (dmaj/dmin)) go to 231
                        if (reftanv .lt. (dter/dmin)) go to 231
                     end if
                     wmin = sqrt(gwts(i)*gwts(j))
                     wminacc(l) = wminacc(l) + wmin
231                  continue
                  end if
               end if
               !
               !     Tertiary direction variogram
               !
               if (useterev .and. i .ne. j) then
                  l = 1 + IDNINT(dter/dterlagev)
                  if (deptypeev .eq. 1) dtertolh = 8*tolfrac*dterlagev
                  if (deptypeev .eq. 1) dtertolv = 4*zavranisev
                  if (deptypeev .eq. 2) dtertolh = 1.5*tolfrac*dterlagev
                  if (deptypeev .eq. 2) dtertolv = 3*zavranisev

                  if (deptypeev .eq. 2 .and. searchdif) dtertolh = 0.5*tolfrac*dterlagev
                  if (deptypeev .eq. 2 .and. searchdif) dtertolv = 0.5*zavranisev

                  if (l .le. nterlagev) then
                     if (usemajev .and. dmaj .gt. dtertolh) go to 232
                     if (useminev .and. dmin .gt. dtertolv) go to 232
                     wter = sqrt(gwts(i)*gwts(j))
                     wteracc(l) = wteracc(l) + wter
232                  continue
                  end if
               end if
            end do
         end do
         ! Get the pairs accumulation
         do i = 1, ndev
            do j = 1, ndev
               if (usexdirev) dxv = xev(i) - xev(j)
               if (useydirev) dyv = yev(i) - yev(j)
               if (usezdirev) dzv = zev(i) - zev(j)
               ! Getting the components using the rotation matrix
               dmaj = abs(lrotmatev(1, 1)*dxv + lrotmatev(1, 2)*dyv + lrotmatev(1, 3)*dzv)
               dmin = abs(lrotmatev(2, 1)*dxv + lrotmatev(2, 2)*dyv + lrotmatev(2, 3)*dzv)
               dter = abs(lrotmatev(3, 1)*dxv + lrotmatev(3, 2)*dyv + lrotmatev(3, 3)*dzv)
               !
               !     Major direction variogram
               !
               if (usemajev .and. i .ne. j) then
                  l = 1 + IDNINT(dmaj/dmajlagev)
                  !set tolerances based on deposit type
                  if (deptypeev .eq. 1) dmajtolh = 2*tolfrac*dmaj
                  if (deptypeev .eq. 1) dmajtolv = 4*zavranisev
                  if (deptypeev .eq. 2) dmajtolh = dmajlagev
                  if (deptypeev .eq. 2) dmajtolv = 1.5*zavranisev

                  if (deptypeev .eq. 2 .and. searchdif) dmajtolh = tolfrac*dmaj
                  if (deptypeev .eq. 2 .and. searchdif) dmajtolv = 1.5*zavranisev
                  if (l .le. nmajlagev) then
                     if (useminev .and. dmin .gt. dmajtolh) go to 240
                     if (useterev .and. dter .gt. dmajtolv) go to 240
                     ! To consider az tol of 22 and dip tol of 5
                     if (deptypeev .eq. 2) then
                        if (reftanh .lt. (dmin/dmaj)) go to 240
                        if (reftanv .lt. (dter/dmaj)) go to 240
                     end if
                     wmaj = sqrt(gwts(i)*gwts(j))
                     lvmajnumev(l) = lvmajnumev(l) + 1
                     lvmajwtev(l) = lvmajwtev(l) + wmaj
                     lvmajdisev(l) = lvmajdisev(l) + dmaj
                     lvmajcovev(l) = lvmajcovev(l) + vev(i)*vev(j)*wmaj/wmajacc(l)
                     lvmajmuev(l) = lvmajmuev(l) + vev(i)*wmaj/wmajacc(l)
                     lvmajmuhev(l) = lvmajmuhev(l) + vev(j)*wmaj/wmajacc(l)
                     lvmajsuev(l) = lvmajsuev(l) + vev(i)*vev(i)*wmaj/wmajacc(l)
                     lvmajsuhev(l) = lvmajsuhev(l) + vev(j)*vev(j)*wmaj/wmajacc(l)
240                  continue
                  end if

               end if
               !
               !     Minor direction variogram
               !
               if (useminev .and. i .ne. j) then
                  l = 1 + IDNINT(dmin/dminlagev)
                  if (deptypeev .eq. 1) dmintolh = 2*tolfrac*dmin
                  if (deptypeev .eq. 1) dmintolv = 4*zavranisev
                  if (deptypeev .eq. 2) dmintolh = dminlagev
                  if (deptypeev .eq. 2) dmintolv = 1.5*zavranisev

                  if (deptypeev .eq. 2 .and. searchdif) dmintolh = tolfrac*dmin
                  if (deptypeev .eq. 2 .and. searchdif) dmintolv = 1.5*zavranisev

                  if (l .le. nminlagev) then
                     if (usemajev .and. dmaj .gt. dmintolh) go to 241
                     if (useterev .and. dter .gt. dmintolv) go to 241
                     if (deptypeev .eq. 2) then
                        if (reftanh .lt. (dmaj/dmin)) go to 241
                        if (reftanv .lt. (dter/dmin)) go to 241
                     end if
                     wmin = sqrt(gwts(i)*gwts(j))
                     lvminnumev(l) = lvminnumev(l) + 1
                     lvminwtev(l) = lvminwtev(l) + wmin
                     lvmindisev(l) = lvmindisev(l) + dmin
                     lvmincovev(l) = lvmincovev(l) + vev(i)*vev(j)*wmin/wminacc(l)
                     lvminmuev(l) = lvminmuev(l) + vev(i)*wmin/wminacc(l)
                     lvminmuhev(l) = lvminmuhev(l) + vev(j)*wmin/wminacc(l)
                     lvminsuev(l) = lvminsuev(l) + vev(i)*vev(i)*wmin/wminacc(l)
                     lvminsuhev(l) = lvminsuhev(l) + vev(j)*vev(j)*wmin/wminacc(l)
241                  continue
                  end if
               end if
               !
               !     Tertiary direction variogram
               !
               if (useterev .and. i .ne. j) then
                  l = 1 + IDNINT(dter/dterlagev)
                  if (deptypeev .eq. 1) dtertolh = 8*tolfrac*dterlagev
                  if (deptypeev .eq. 1) dtertolv = 4*zavranisev
                  if (deptypeev .eq. 2) dtertolh = 1.5*tolfrac*dterlagev
                  if (deptypeev .eq. 2) dtertolv = 3*zavranisev

                  if (deptypeev .eq. 2 .and. searchdif) dtertolh = 0.5*tolfrac*dterlagev
                  if (deptypeev .eq. 2 .and. searchdif) dtertolv = 0.5*zavranisev

                  if (l .le. nterlagev) then
                     if (usemajev .and. dmaj .gt. dtertolh) go to 242
                     if (useminev .and. dmin .gt. dtertolv) go to 242
                     wter = sqrt(gwts(i)*gwts(j))
                     lvternumev(l) = lvternumev(l) + 1
                     lvterwtev(l) = lvterwtev(l) + wter
                     lvterdisev(l) = lvterdisev(l) + dter
                     lvtercovev(l) = lvtercovev(l) + vev(i)*vev(j)*wter/wteracc(l)
                     lvtermuev(l) = lvtermuev(l) + vev(i)*wter/wteracc(l)
                     lvtermuhev(l) = lvtermuhev(l) + vev(j)*wter/wteracc(l)
                     lvtersuev(l) = lvtersuev(l) + vev(i)*vev(i)*wter/wteracc(l)
                     lvtersuhev(l) = lvtersuhev(l) + vev(j)*vev(j)*wter/wteracc(l)
242                  continue
                  end if
               end if
            end do
         end do
         ! Get the weighted correlogram
         do l = 1, nmajlagev + 1
            if (usemajev .and. lvmajnumev(l) .ge. 5) then
               lvmajdisev(l) = lvmajdisev(l)/real(lvmajnumev(l))
               lvmajcovev(l) = lvmajcovev(l) - lvmajmuev(l)*lvmajmuhev(l)! wt cov
               lvmajsuev(l) = lvmajsuev(l) - lvmajmuev(l)*lvmajmuev(l)
               lvmajsuhev(l) = lvmajsuhev(l) - lvmajmuhev(l)*lvmajmuhev(l)
               if (lvmajsuev(l) .lt. EPSLON) lvmajsuev(l) = EPSLON
               if (lvmajsuhev(l) .lt. EPSLON) lvmajsuhev(l) = EPSLON
               lvmajcovev(l) = lvmajcovev(l)/sqrt(lvmajsuev(l)*lvmajsuhev(l))!wt corr
               lvmajwtev(l) = 1
            else
               lvmajnumev(l) = 0
               lvmajcovev(l) = 0.0
               lvmajdisev(l) = 0.0
               lvmajwtev(l) = 0.0
            end if
         end do
         do l = 1, nminlagev + 1
            if (useminev .and. lvminnumev(l) .ge. 5) then
               lvmindisev(l) = lvmindisev(l)/real(lvminnumev(l))
               lvmincovev(l) = lvmincovev(l) - lvminmuev(l)*lvminmuhev(l)
               lvminsuev(l) = lvminsuev(l) - lvminmuev(l)*lvminmuev(l)
               lvminsuhev(l) = lvminsuhev(l) - lvminmuhev(l)*lvminmuhev(l)
               if (lvminsuev(l) .lt. EPSLON) lvminsuev(l) = EPSLON
               if (lvminsuhev(l) .lt. EPSLON) lvminsuhev(l) = EPSLON
               lvmincovev(l) = lvmincovev(l)/sqrt(lvminsuev(l)*lvminsuhev(l))
               lvminwtev(l) = 1
            else
               lvminnumev(l) = 0
               lvmincovev(l) = 0.0
               lvmindisev(l) = 0.0
               lvminwtev(l) = 0.0
            end if
         end do

         do l = 1, nterlagev + 1
            if (useterev .and. lvternumev(l) .ge. 5) then
               lvterdisev(l) = lvterdisev(l)/real(lvternumev(l))
               lvtercovev(l) = lvtercovev(l) - lvtermuev(l)*lvtermuhev(l)
               lvtersuev(l) = lvtersuev(l) - lvtermuev(l)*lvtermuev(l)
               lvtersuhev(l) = lvtersuhev(l) - lvtermuhev(l)*lvtermuhev(l)
               if (lvtersuev(l) .lt. EPSLON) lvtersuev(l) = EPSLON
               if (lvtersuhev(l) .lt. EPSLON) lvtersuhev(l) = EPSLON
               lvtercovev(l) = lvtercovev(l)/sqrt(lvtersuev(l)*lvtersuhev(l))
               lvterwtev(l) = 1
            else
               lvternumev(l) = 0
               lvtercovev(l) = 0.0
               lvterdisev(l) = 0.0
               lvterwtev(l) = 0.0
            end if
         end do
         !
         ! Convert the correlogram to variogram values, since it is already
         ! normalized: var = 1 - corr
         do l = 1, nmajlagev + 1
            lvmajcovev(l) = 1 - lvmajcovev(l)
         end do
         do l = 1, nminlagev + 1
            lvmincovev(l) = 1 - lvmincovev(l)
         end do
         do l = 1, nterlagev + 1
            lvtercovev(l) = 1 - lvtercovev(l)
         end do

      end if

      !
      ! Remove where the lag distance is zero
      !
      ncountmaj = 0; ncountmin = 0; ncountter = 0
      do l = 1, nmajlagev
         if (lvmajdisev(l) .gt. 0.0) ncountmaj = ncountmaj + 1
      end do
      do l = 1, nminlagev
         if (lvmindisev(l) .gt. 0.0) ncountmin = ncountmin + 1
      end do
      do l = 1, nterlagev
         if (useterev) then
            if (lvterdisev(l) .gt. 0.0) ncountter = ncountter + 1
         else
            ncountter = 1
         end if

      end do
      ! Allocate temporaty arrays to deallocate outputs array
      tlvmajnumev = lvmajnumev; tlvmajcovev = lvmajcovev
      tlvmajmuev = lvmajmuev; tlvmajmuhev = lvmajmuhev
      tlvmajsuev = lvmajsuev; tlvmajsuhev = lvmajsuhev
      tlvmajdisev = lvmajdisev; tlvmajwtev = lvmajwtev

      tlvminnumev = lvminnumev; tlvmincovev = lvmincovev
      tlvminmuev = lvminmuev; tlvminmuhev = lvminmuhev
      tlvminsuev = lvminsuev; tlvminsuhev = lvminsuhev
      tlvmindisev = lvmindisev; tlvminwtev = lvminwtev

      tlvternumev = lvternumev; tlvtercovev = lvtercovev
      tlvtermuev = lvtermuev; tlvtermuhev = lvtermuhev
      tlvtersuev = lvtersuev; tlvtersuhev = lvtersuhev
      tlvterdisev = lvterdisev; tlvterwtev = lvterwtev

      ! Deallocate output arrays
      deallocate (lvmajnumev, lvmajcovev, lvmajmuev, lvmajmuhev, &
                  lvmajsuev, lvmajsuhev, lvmajdisev, lvmajwtev)
      deallocate (lvternumev, lvtercovev, lvtermuev, lvtermuhev, &
                  lvtersuev, lvtersuhev, lvterdisev, lvterwtev)
      deallocate (lvminnumev, lvmincovev, lvminmuev, lvminmuhev, &
                  lvminsuev, lvminsuhev, lvmindisev, lvminwtev)

      ! Allocate for valid values only
      allocate (lvmajnumev(ncountmaj))
      allocate (lvmajcovev, lvmajmuev, lvmajmuhev, lvmajsuev, &
                lvmajsuhev, lvmajdisev, lvmajwtev, mold=lvmajnumev)
      allocate (lvminnumev(ncountmin))
      allocate (lvmincovev, lvminmuev, lvminmuhev, lvminsuev, &
                lvminsuhev, lvmindisev, lvminwtev, mold=lvminnumev)
      allocate (lvternumev(ncountter))
      allocate (lvtercovev, lvtermuev, lvtermuhev, lvtersuev, &
                lvtersuhev, lvterdisev, lvterwtev, mold=lvternumev)

      ! fill the arrays
      j = 0
      do i = 1, nmajlagev
         if (tlvmajdisev(i) .gt. 0.) then
            j = j + 1
            lvmajnumev(j) = tlvmajnumev(i); lvmajcovev(j) = tlvmajcovev(i)
            lvmajmuev(j) = tlvmajmuev(i); lvmajmuhev(j) = tlvmajmuhev(i)
            lvmajsuev(j) = tlvmajsuev(i); lvmajsuhev(j) = tlvmajsuhev(i)
            lvmajdisev(j) = tlvmajdisev(i); lvmajwtev(j) = tlvmajwtev(i)
         end if
      end do
      j = 0
      do i = 1, nminlagev
         if (tlvmindisev(i) .gt. 0.) then
            j = j + 1
            lvminnumev(j) = tlvminnumev(i); lvmincovev(j) = tlvmincovev(i)
            lvminmuev(j) = tlvminmuev(i); lvminmuhev(j) = tlvminmuhev(i)
            lvminsuev(j) = tlvminsuev(i); lvminsuhev(j) = tlvminsuhev(i)
            lvmindisev(j) = tlvmindisev(i); lvminwtev(j) = tlvminwtev(i)
         end if
      end do
      j = 0
      if (useterev) then
         do i = 1, nterlagev
            if (tlvterdisev(i) .gt. 0.) then
               j = j + 1
               lvternumev(j) = tlvternumev(i); lvtercovev(j) = tlvtercovev(i)
               lvtermuev(j) = tlvtermuev(i); lvtermuhev(j) = tlvtermuhev(i)
               lvtersuev(j) = tlvtersuev(i); lvtersuhev(j) = tlvtersuhev(i)
               lvterdisev(j) = tlvterdisev(i); lvterwtev(j) = tlvterwtev(i)
            end if
         end do
      end if

      ! update number of lags for output

      nmajlagev = ncountmaj
      nminlagev = ncountmin
      if (useterev) nterlagev = ncountter

      !
      ! Check first point and set the rigth weight: if the variogram value is
      ! greater than the second point and if the number of pairs is less
      !
      i = 0; j = 0; k = 0
      if (lvmajcovev(1) .gt. lvmajcovev(2) .and. lvmajnumev(1) .lt. lvmajnumev(2)) &
         lvmajwtev(1) = 0.0
      if (lvmincovev(1) .gt. lvmincovev(2) .and. lvminnumev(1) .lt. lvminnumev(2)) &
         lvminwtev(1) = 0.0
      if (useterev) then
         if (lvtercovev(1) .gt. lvtercovev(2) .and. lvternumev(1) .lt. lvternumev(2)) &
            lvterwtev(1) = 0.0
      end if
      ! Checking only for 2D datasets
      if (.not. useterev) then
         do i = nminlagev, 2, -1
            if (lvmincovev(i) .lt. 1. .and. lvmindisev(i) .gt. 0.5*maxval(lvmindisev)) then
               if (lvmincovev(i) .lt. lvmincovev(i - 1)) lvminwtev(i) = 0.
            end if
         end do
      end if
      ! Checking only for 3D datasets.
      if (useterev) then
         do i = nterlagev, 2, -1
            if (lvtercovev(i) .lt. 1. .and. lvterdisev(i) .gt. 0.5*maxval(lvterdisev)) then
               if (lvtercovev(i) .lt. lvtercovev(i - 1)) lvterwtev(i) = 0.
            end if
         end do
      end if
      !
      ! Show results on screen
      !
      if (debugev) then
         write (*, *) ' Inferred experimental variogram: '
         if (usemajev) then
            write (*, *)
            write (*, *) ' Major variogram '
            do l = 1, nmajlagev
               write (*, 120) lvmajnumev(l), lvmajdisev(l), lvmajcovev(l)
            end do
         end if
         if (useminev) then
            write (*, *)
            write (*, *) ' Minor variogram '
            do l = 1, nminlagev
               write (*, 120) lvminnumev(l), lvmindisev(l), lvmincovev(l)
            end do
         end if
         if (useterev) then
            write (*, *)
            write (*, *) ' Tertiary variogram '
            do l = 1, nterlagev
               write (*, 120) lvternumev(l), lvterdisev(l), lvtercovev(l)
            end do
         end if
      end if
120   format(4x, 4f12.4, 2x, 4f12.4)
      return
   end subroutine varexpflx

   subroutine wtpairs(xyz, v, wts, lagsize, nlag, vspc, nd, deptype, sp, anchx, &
                      anchy, anchz, azmtol, diptol, dir, danxs, wtss, hdist, lagns, majwt)
      !
      ! Subroutine to check all possible pairs and calculate the weights.
      !
      ! Variable declaration
      use sortem
      implicit none

      real, parameter :: EPSLON = 1.0e-6
      real(kind=8), parameter :: PI = 4*atan(1.0d0), DEG2RAD = PI/180d0, RAD2DEG = 180d0/PI
      real(kind=8), parameter :: BIGDBLE = 1d21, SMALLDBLE = 1d-6
      ! inputs
      integer, intent(in) :: nd, deptype, nlag, dir
      real(kind=8), dimension(nd, 3), intent(in):: xyz
      real(kind=8), dimension(nd), intent(in):: v
      real(kind=8), dimension(:), intent(in):: wts
      real(kind=8), intent(in) :: sp, azmtol, diptol, lagsize, vspc, anchx, &
                                  anchy, anchz
      ! outputs:
      real(kind=8), allocatable, dimension(:), intent(out) :: danxs, wtss, hdist, lagns, majwt
      real(kind=8), allocatable, dimension(:) :: xxxnum

      ! internal variables
      character(len=100) :: output_line
      integer :: i_exp, tmp_anx, nmaj, i, j, l
      real(kind=8) :: azmtoldist, diptoldist, ymax, ymin, zmax, zmin, lagtol, &
                      xdist, angle, h, midx, midy, midz, danx, wt, wmaj, sp2
      real(kind=8), parameter ::g = 2.0, gg = 2.0

      ! code
      !Allocate the arrays for the weights
      allocate (majwt(nlag + 1))
      allocate (xxxnum(nlag + 1))
      xxxnum = 0.0
      majwt = 0.0
      i_exp = 9
      nmaj = 0
      !
      ! Major direction
      !

      if (dir .eq. 1) then
         azmtoldist = 0.5*lagsize/tan(azmtol)
         diptoldist = 0.5*vspc/tan(diptol)
         do i = 1, nd
            ! precalculate filters
            lagtol = 0.5*lagsize
            if (deptype .eq. 1) then
               ymax = xyz(i, 2) + lagsize
               ymin = xyz(i, 2) - lagsize
               zmax = xyz(i, 3) + vspc
               zmin = xyz(i, 3) - vspc
            end if

            if (deptype .eq. 2) then
               ymax = xyz(i, 2) + 0.5*lagsize
               ymin = xyz(i, 2) - 0.5*lagsize
               zmax = xyz(i, 3) + 0.5*vspc
               zmin = xyz(i, 3) - 0.5*vspc
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
               ! Get average of weights for this pair
               midx = (xyz(i, 1) + xyz(j, 1))/2
               midy = (xyz(i, 2) + xyz(j, 2))/2
               midz = (xyz(i, 3) + xyz(j, 3))/2
               ! additionally weighted base on the distance.
               danx = sqrt((midx - anchx)**2 + &
                           (midy - anchy)**2 + &
                           (midz - anchz)**2)
               wt = exp(-(danx**g/(2*(sp**gg)))) + EPSLON
               wmaj = sqrt(wts(i)*wts(j))
               ! Calculating the sum of wts per lag
               do l = 0, nlag
                  if ((h .ge. (l*lagsize - lagtol)) .and. &
                      (h .le. (l*lagsize + lagtol))) then
                     majwt(l + 1) = majwt(l + 1) + wmaj
                     xxxnum(l + 1) = xxxnum(l + 1) + 1
                     nmaj = nmaj + 1
                  end if
               end do
               ! For calculating the acc values for the correlogram.
            end do
         end do
         allocate (danxs(nmaj))
         allocate (wtss, hdist, lagns, mold=danxs)
         danxs = 1; wtss = 0; lagns = 1; hdist = 5
      end if
      !
      ! Minor direction
      !
      if (dir .eq. 2) then
         azmtoldist = 0.5*lagsize/tan(azmtol)
         diptoldist = 0.5*vspc/tan(diptol)
         do i = 1, nd
            ! precalculate filters
            lagtol = 0.5*lagsize
            if (deptype .eq. 1) then
               ymax = xyz(i, 1) + lagsize
               ymin = xyz(i, 1) - lagsize
               zmax = xyz(i, 3) + vspc
               zmin = xyz(i, 3) - vspc
            end if

            if (deptype .eq. 2) then
               ymax = xyz(i, 1) + 0.5*lagsize
               ymin = xyz(i, 1) - 0.5*lagsize
               zmax = xyz(i, 3) + 0.5*vspc
               zmin = xyz(i, 3) - 0.5*vspc
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
               ! Get average of weights for this pair
               midx = (xyz(i, 1) + xyz(j, 1))/2
               midy = (xyz(i, 2) + xyz(j, 2))/2
               midz = (xyz(i, 3) + xyz(j, 3))/2
               ! additionally weighted base on the distance.
               danx = sqrt((midx - anchx)**2 + &
                           (midy - anchy)**2 + &
                           (midz - anchz)**2)
               !if (danx .gt. 30.0) cycle
               wt = exp(-(danx**g/(2*(sp**gg)))) + EPSLON
               wmaj = sqrt(wts(i)*wts(j))
               ! Calculating the sum of wts per lag
               do l = 0, nlag
                  if ((h .ge. (l*lagsize - lagtol)) .and. &
                      (h .le. (l*lagsize + lagtol))) then
                     majwt(l + 1) = majwt(l + 1) + wmaj
                     nmaj = nmaj + 1
                  end if
               end do
            end do
         end do
         allocate (danxs(nmaj))
         allocate (wtss, hdist, lagns, mold=danxs)
         danxs = 1; wtss = 0; lagns = 1; hdist = 5
      end if
      !
      ! Tertiary direction
      !

      if (dir .eq. 3) then
         azmtoldist = 0.5*vspc/tan(azmtol)
         diptoldist = 0.5*vspc/tan(diptol)
         do i = 1, nd
            ! precalculate filters
            lagtol = 0.5*vspc
            if (deptype .eq. 1) then
               ymax = xyz(i, 1) + lagsize
               ymin = xyz(i, 1) - lagsize
               zmax = xyz(i, 2) + 4*vspc
               zmin = xyz(i, 2) - 4*vspc
            end if

            if (deptype .eq. 2) then
               ymax = xyz(i, 1) + 0.5*lagsize
               ymin = xyz(i, 1) - 0.5*lagsize
               zmax = xyz(i, 2) + 0.5*vspc
               zmin = xyz(i, 2) - 0.5*vspc
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
               ! Get average of weights for this pair
               midx = (xyz(i, 1) + xyz(j, 1))/2
               midy = (xyz(i, 2) + xyz(j, 2))/2
               midz = (xyz(i, 3) + xyz(j, 3))/2
               ! additionally weighted base on the distance.
               danx = sqrt((midx - anchx)**2 + &
                           (midy - anchy)**2 + &
                           (midz - anchz)**2)
               wt = exp(-(danx**g/(2*(sp**gg)))) + EPSLON
               wmaj = sqrt(wts(i)*wts(j))
               ! Calculating the sum of wts per lag
               do l = 0, nlag
                  if ((h .ge. (l*lagsize - lagtol)) .and. &
                      (h .le. (l*lagsize + lagtol))) then
                     majwt(l + 1) = majwt(l + 1) + wmaj
                     nmaj = nmaj + 1
                  end if
               end do
            end do
         end do
         allocate (danxs(nmaj))
         allocate (wtss, hdist, lagns, mold=danxs)
         danxs = 1; wtss = 0; lagns = 1; hdist = 5
      end if
      close (i_exp)
   end subroutine wtpairs

   ! Subroutine to calculate the values for the correlogram

   subroutine corrc(xyz, v, wts, accwt, lagsize, nlag, vspc, nd, deptype, sp, anchx, &
                    anchy, anchz, azmtol, diptol, dir, glc, lvxxxnum, lvxxxdis, lvxxxcov, lvxxxmu, &
                    lvxxxmuh, lvxxxsu, lvxxxsuh, varxxxglc)
      !
      ! Subroutine to check all possible pairs and calculate the weights.
      !
      ! Variable declaration
      use sortem
      implicit none

      real, parameter :: EPSLON = 1.0e-6
      real(kind=8), parameter :: PI = 4*atan(1.0d0), DEG2RAD = PI/180d0, RAD2DEG = 180d0/PI
      real(kind=8), parameter :: BIGDBLE = 1d21, SMALLDBLE = 1d-6
      ! inputs
      integer, intent(in) :: nd, deptype, nlag, dir, glc
      real(kind=8), dimension(nd, 3), intent(in):: xyz
      real(kind=8), dimension(nd), intent(in):: v
      real(kind=8), dimension(:), intent(in):: wts, accwt
      real(kind=8), intent(in) :: sp, azmtol, diptol, lagsize, vspc, anchx, &
                                  anchy, anchz
      ! outputs:
      real(kind=8), allocatable, dimension(:), intent(out) ::lvxxxnum, lvxxxdis, &
                                                              lvxxxcov, lvxxxmu, lvxxxmuh, lvxxxsu, lvxxxsuh
      real(kind=8), intent(out) :: varxxxglc
      ! internal variables
      integer :: tmp_anx, nmaj, i, j, l
      real(kind=8) :: azmtoldist, diptoldist, ymax, ymin, zmax, zmin, lagtol, &
                      xdist, angle, h, midx, midy, midz, danx, wt, wmaj, sp2, sumweight, partial
      real(kind=8), allocatable, dimension(:) :: xxxwt
      real(kind=8), parameter ::g = 2.0, gg = 2.0

      ! code
      !Allocate the arrays for the weights
      allocate (xxxwt(nlag + 1))
      ! Allocate arrays for correlogram calculation
      allocate (lvxxxnum, lvxxxdis, lvxxxcov, lvxxxmu, &
                lvxxxmuh, lvxxxsu, lvxxxsuh, mold=xxxwt)
      ! Set all to zero to start
      lvxxxnum = 0.0; lvxxxdis = 0.0; lvxxxcov = 0.0; lvxxxmu = 0.0
      lvxxxmuh = 0.0; lvxxxsu = 0.0; lvxxxsuh = 0.0; xxxwt = 0.0
      ! If the weights were already calculated:
      xxxwt = accwt
      varxxxglc = 0.0
      !
      ! Major direction
      !

      if (dir .eq. 1) then
         azmtoldist = 0.5*lagsize/tan(azmtol)
         diptoldist = 0.5*vspc/tan(diptol)
         do i = 1, nd
            ! precalculate filters
            lagtol = 0.5*lagsize
            if (deptype .eq. 1) then
               ymax = xyz(i, 2) + lagsize
               ymin = xyz(i, 2) - lagsize
               zmax = xyz(i, 3) + vspc
               zmin = xyz(i, 3) - vspc
            end if

            if (deptype .eq. 2) then
               ymax = xyz(i, 2) + 0.5*lagsize
               ymin = xyz(i, 2) - 0.5*lagsize
               zmax = xyz(i, 3) + 0.5*vspc
               zmin = xyz(i, 3) - 0.5*vspc
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
               if (glc .eq. 1) then
                  ! Get average of weights for this pair
                  midx = (xyz(i, 1) + xyz(j, 1))/2
                  midy = (xyz(i, 2) + xyz(j, 2))/2
                  midz = (xyz(i, 3) + xyz(j, 3))/2
                  ! additionally weighted base on the distance.
                  danx = sqrt((midx - anchx)**2 + &
                              (midy - anchy)**2 + &
                              (midz - anchz)**2)
                  wt = exp(-(danx**g/(2*(sp**gg)))) + EPSLON
                  wmaj = sqrt(wts(i)*wts(j))
                  ! Calculating the acc values for the correlogram.
                  do l = 0, nlag
                     if ((h .ge. (l*lagsize - lagtol)) .and. &
                         (h .le. (l*lagsize + lagtol))) then
                        lvxxxnum(l + 1) = lvxxxnum(l + 1) + 1
                        lvxxxdis(l + 1) = lvxxxdis(l + 1) + xdist
                        lvxxxcov(l + 1) = lvxxxcov(l + 1) + ((v(i)*v(j)))*wmaj/xxxwt(l + 1)
                        lvxxxmu(l + 1) = lvxxxmu(l + 1) + v(i)*wmaj/xxxwt(l + 1)
                        lvxxxmuh(l + 1) = lvxxxmuh(l + 1) + v(j)*wmaj/xxxwt(l + 1)
                        varxxxglc = varxxxglc + wmaj*((v(i) - v(j))**2)
                        lvxxxsu(l + 1) = lvxxxsu(l + 1) + v(i)*v(i)*wmaj/xxxwt(l + 1)
                        lvxxxsuh(l + 1) = lvxxxsuh(l + 1) + v(j)*v(j)*wmaj/xxxwt(l + 1)
                     end if
                  end do
               else
                  do l = 0, nlag
                     if ((h .ge. (l*lagsize - lagtol)) .and. &
                         (h .le. (l*lagsize + lagtol))) then
                        lvxxxnum(l + 1) = lvxxxnum(l + 1) + 1
                        lvxxxdis(l + 1) = lvxxxdis(l + 1) + xdist
                        lvxxxcov(l + 1) = lvxxxcov(l + 1) + v(i)*v(j)
                        lvxxxmu(l + 1) = lvxxxmu(l + 1) + v(i)
                        lvxxxmuh(l + 1) = lvxxxmuh(l + 1) + v(j)
                        lvxxxsu(l + 1) = lvxxxsu(l + 1) + v(i)*v(i)
                        lvxxxsuh(l + 1) = lvxxxsuh(l + 1) + v(j)*v(j)
                     end if
                  end do
               end if
            end do
         end do
         if (glc .eq. 1) then
            ! Get the sum of weights per lag
            sumweight = 0.0
            do l = 0, nlag
               partial = xxxwt(l + 1)
               if (partial .gt. 0) sumweight = sumweight + partial
            end do
            varxxxglc = varxxxglc/(2*sumweight) ! geoloc variance
         end if
      end if
      !
      ! Minor direction
      !
      if (dir .eq. 2) then
         azmtoldist = 0.5*lagsize/tan(azmtol)
         diptoldist = 0.5*vspc/tan(diptol)
         do i = 1, nd
            ! precalculate filters
            lagtol = 0.5*lagsize
            if (deptype .eq. 1) then
               ymax = xyz(i, 1) + lagsize
               ymin = xyz(i, 1) - lagsize
               zmax = xyz(i, 3) + vspc
               zmin = xyz(i, 3) - vspc
            end if

            if (deptype .eq. 2) then
               ymax = xyz(i, 1) + 0.5*lagsize
               ymin = xyz(i, 1) - 0.5*lagsize
               zmax = xyz(i, 3) + 0.5*vspc
               zmin = xyz(i, 3) - 0.5*vspc
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
               if (glc .eq. 1) then
                  ! Get average of weights for this pair
                  midx = (xyz(i, 1) + xyz(j, 1))/2
                  midy = (xyz(i, 2) + xyz(j, 2))/2
                  midz = (xyz(i, 3) + xyz(j, 3))/2
                  ! additionally weighted base on the distance.
                  danx = sqrt((midx - anchx)**2 + &
                              (midy - anchy)**2 + &
                              (midz - anchz)**2)
                  wt = exp(-(danx**g/(2*(sp**gg)))) + EPSLON
                  wmaj = sqrt(wts(i)*wts(j))
                  ! Calculating the acc values for the correlogram.
                  do l = 0, nlag
                     if ((h .ge. (l*lagsize - lagtol)) .and. &
                         (h .le. (l*lagsize + lagtol))) then
                        lvxxxnum(l + 1) = lvxxxnum(l + 1) + 1
                        lvxxxdis(l + 1) = lvxxxdis(l + 1) + xdist
                        lvxxxcov(l + 1) = lvxxxcov(l + 1) + ((v(i)*v(j)))*wmaj/xxxwt(l + 1)
                        lvxxxmu(l + 1) = lvxxxmu(l + 1) + v(i)*wmaj/xxxwt(l + 1)
                        lvxxxmuh(l + 1) = lvxxxmuh(l + 1) + v(j)*wmaj/xxxwt(l + 1)
                        varxxxglc = varxxxglc + wmaj*((v(i) - v(j))**2)
                        lvxxxsu(l + 1) = lvxxxsu(l + 1) + v(i)*v(i)*wmaj/xxxwt(l + 1)
                        lvxxxsuh(l + 1) = lvxxxsuh(l + 1) + v(j)*v(j)*wmaj/xxxwt(l + 1)
                     end if
                  end do
               else
                  do l = 0, nlag
                     if ((h .ge. (l*lagsize - lagtol)) .and. &
                         (h .le. (l*lagsize + lagtol))) then
                        lvxxxnum(l + 1) = lvxxxnum(l + 1) + 1
                        lvxxxdis(l + 1) = lvxxxdis(l + 1) + xdist
                        lvxxxcov(l + 1) = lvxxxcov(l + 1) + v(i)*v(j)
                        lvxxxmu(l + 1) = lvxxxmu(l + 1) + v(i)
                        lvxxxmuh(l + 1) = lvxxxmuh(l + 1) + v(j)
                        lvxxxsu(l + 1) = lvxxxsu(l + 1) + v(i)*v(i)
                        lvxxxsuh(l + 1) = lvxxxsuh(l + 1) + v(j)*v(j)
                     end if
                  end do
               end if
            end do
         end do
         if (glc .eq. 1) then
            ! Get the sum of weights per lag
            sumweight = 0.0
            do l = 0, nlag
               partial = xxxwt(l + 1)
               if (partial .gt. 0) sumweight = sumweight + partial
            end do
            varxxxglc = varxxxglc/(2*sumweight) ! geoloc variance
         end if
      end if
      !
      ! Tertiary direction
      !

      if (dir .eq. 3) then
         azmtoldist = 0.5*vspc/tan(azmtol)
         diptoldist = 0.5*vspc/tan(diptol)
         do i = 1, nd
            ! precalculate filters
            lagtol = 0.5*vspc
            if (deptype .eq. 1) then
               ymax = xyz(i, 1) + lagsize
               ymin = xyz(i, 1) - lagsize
               zmax = xyz(i, 2) + 4*vspc
               zmin = xyz(i, 2) - 4*vspc
            end if

            if (deptype .eq. 2) then
               ymax = xyz(i, 1) + 0.5*lagsize
               ymin = xyz(i, 1) - 0.5*lagsize
               zmax = xyz(i, 2) + 0.5*vspc
               zmin = xyz(i, 2) - 0.5*vspc
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
               if (glc .eq. 1) then
                  ! Get average of weights for this pair
                  midx = (xyz(i, 1) + xyz(j, 1))/2
                  midy = (xyz(i, 2) + xyz(j, 2))/2
                  midz = (xyz(i, 3) + xyz(j, 3))/2
                  ! additionally weighted base on the distance.
                  danx = sqrt((midx - anchx)**2 + &
                              (midy - anchy)**2 + &
                              (midz - anchz)**2)
                  wt = exp(-(danx**g/(2*(sp**gg)))) + EPSLON
                  wmaj = sqrt(wts(i)*wts(j))
                  ! Calculating the acc values for the correlogram.
                  do l = 0, nlag
                     if ((h .ge. (l*lagsize - lagtol)) .and. &
                         (h .le. (l*lagsize + lagtol))) then
                        lvxxxnum(l + 1) = lvxxxnum(l + 1) + 1
                        lvxxxdis(l + 1) = lvxxxdis(l + 1) + xdist
                        lvxxxcov(l + 1) = lvxxxcov(l + 1) + ((v(i)*v(j)))*wmaj/xxxwt(l + 1)
                        lvxxxmu(l + 1) = lvxxxmu(l + 1) + v(i)*wmaj/xxxwt(l + 1)
                        lvxxxmuh(l + 1) = lvxxxmuh(l + 1) + v(j)*wmaj/xxxwt(l + 1)
                        varxxxglc = varxxxglc + wmaj*((v(i) - v(j))**2)
                        lvxxxsu(l + 1) = lvxxxsu(l + 1) + v(i)*v(i)*wmaj/xxxwt(l + 1)
                        lvxxxsuh(l + 1) = lvxxxsuh(l + 1) + v(j)*v(j)*wmaj/xxxwt(l + 1)
                     end if
                  end do
               else
                  do l = 0, nlag
                     if ((h .ge. (l*lagsize - lagtol)) .and. &
                         (h .le. (l*lagsize + lagtol))) then
                        lvxxxnum(l + 1) = lvxxxnum(l + 1) + 1
                        lvxxxdis(l + 1) = lvxxxdis(l + 1) + xdist
                        lvxxxcov(l + 1) = lvxxxcov(l + 1) + v(i)*v(j)
                        lvxxxmu(l + 1) = lvxxxmu(l + 1) + v(i)
                        lvxxxmuh(l + 1) = lvxxxmuh(l + 1) + v(j)
                        lvxxxsu(l + 1) = lvxxxsu(l + 1) + v(i)*v(i)
                        lvxxxsuh(l + 1) = lvxxxsuh(l + 1) + v(j)*v(j)
                     end if
                  end do
               end if
            end do
         end do
         if (glc .eq. 1) then
            ! Get the sum of weights per lag
            sumweight = 0.0
            do l = 0, nlag
               partial = xxxwt(l + 1)
               if (partial .gt. 0) sumweight = sumweight + partial
            end do
            varxxxglc = varxxxglc/(2*sumweight) ! geoloc variance
         end if
      end if
   end subroutine corrc
end module varcalcs
