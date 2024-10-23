!  Autovar.f90
!
!****************************************************************************
!
! Notes:
! This program takes as input the data file and produce the experimental
! variogramand its model. The only input need it is the information of
! the data file; additionally, the practitioner picks: (1)if the geolocated
! option is used, (2) the anchor information and (3) the deposit type.
! When using anchors each pair will get a weight base on the distance
! to the anchor (Gaussian Wt function).
!
! This program adapted or take programs/subroutines from other authors:
! Clayton V. Deutsch, Jared Deutsch, Matthew V. Deutsch, John Manchuck,
! D. Khan, etc. Includes GSLIB subroutines like varcalc, varmodel, datspac,
! sortem, cova, stratcoord, etc.
!
! This is a research software and the final results need to be evaluated
! by the practitioner. The author and contributors disclaim any liability
! for the use of this software.
!
!  Autovar structure
! 1. Initial set up base on parameter file
! 2. Identify directions of continuity using MOI
! 3. Experimental variogram calculation on the directions from 2
! 4. Automatic variogram modeling
! 5. Evaluation of the model base on anisotropy relationships(if selected
!    by user)
! 6. Output writing: experimental variogram and model, rotation matrix
!    (directions of continuity),weights per anchor (if selected)
! If the geolocated option is selected, the program will loop from 2 to 6.
!
! Author: Luis Davila
! Initial date Oct 2022
! Finish date June 30 2023
! Version 1.000: October 6, 2023
! Version 1.001: November 14, 2023
! Version 1.002: January 30, 2024
! Version 1.003: March 23, 2024
!****************************************************************************
program Autovar
   use sortem ! Quicksort module
   use gen_subs! for general subroutines in modules
   use spacing ! To calculate horizontal and vertical space.
   use alldir ! omnidirectional calculation of pairs
   use moi_aux ! MOI algorithm
   use sp3d ! Module to calculate 3D spacing
   use varcalcs! for pair calculation and correlogram calc
   use random ! Mersenne-Twister random number generator
   !****************************************************************************
   ! Variables
   !****************************************************************************
   implicit none
   ! Parameters
   integer, parameter :: MAXROT = 1
   real, parameter :: EPSLON = 1.0e-6
   real(kind=8), parameter :: MAXCAT = 32, VERSION = 1.003, PI = 4*atan(1.0d0)
   real(kind=8), parameter :: BIGDBLE = 1d21, SMALLDBLE = 1d-6
   real(kind=8), parameter :: DEG2RAD = PI/180d0, RAD2DEG = 180d0/PI
   ! Variables declared by user
   !   file names
   character*512 :: datafl, outfl_mod, outfl_exp, str, dbgfl, varmodgs, wtsfl, anxfl
   integer :: i_exp, i_mod, lout, lun, ndir
   !   Columns and values
   integer :: idcol, xcol, ycol, zcol, vcol, wcol, geolocvar, nx, ny, nz, &
              dep_type, anis, varimp
   real(kind=8) :: rma, rmi, rte, tmin, tmax, xmn, ymn, zmn, xsiz, ysiz, zsiz, &
                   rtmajmin, rtmajter
   real(kind=8), allocatable :: x(:), y(:), z(:), v(:), id(:), w(:), &
                                xt(:), yt(:), zt(:), wt(:), iidd(:), vv(:), xxx(:), yyy(:), zzz(:)
   !   related logicals:
   logical :: not_z, not_dhid, dat2d, dat3d
   ! Variables used if geolocated option selected
   integer :: nanchor, nach, anx, nax
   real(kind=8) :: sdist, rdist, radsqd, ddax, ddaxx, ddaxy, ddaxz, sr
   real(kind=8), dimension(:), allocatable :: anchx, anchy, anchz, sumw, xx, yy, zz, &
                                              vvax, wtsax
   character(512), dimension(:), allocatable :: modelnames, expnames, dbgnames, &
                                                modgsnames, wtsflnames
   character(256) :: ii
   character(len=8) :: fmt ! format descriptor
   real(kind=8), allocatable :: wtsanch(:, :), wtscurrent(:, :)
   real(kind=8), allocatable :: wtsomni(:), wtsomni1(:)
   integer :: axcol, aycol, azcol, nanchororig
   logical :: geoloc
   ! Internal variables
   !   Pair selection:
   integer :: outputidx
   integer, allocatable :: idf(:)
   real(kind=8), allocatable :: sq_df(:), dx(:), dy(:), dz(:), sqdfdist(:)
   real(kind=8), allocatable :: xs(:), ys(:), zs(:), vs(:), xyzs(:, :)
   real(kind=8), dimension(:, :), allocatable :: xyz, xyz_anis, xyzref
   real(kind=8) :: max_dist, h
   !   MOI:
   integer :: outputidx1, z_count, nd_tmp
   real(kind=8) :: a1, a2, a3, dom_size, perc_99val, eigt1, eigt2
   real(kind=8) :: eigveco(3, 3), introtindexx(3), dblerotindexx(3), &
                   eigvec(3, 3), eigval(3), moimat(3, 3), lrotmat(3, 3), right_eiveco(3, 3), &
                   lrotmatrob(3, 3)
   real(kind=8) :: aa, bb, tt, azm_opt, dip_opt, tilt_opt, err, besterr, &
                   rotmat_t(3, 3), revrotmat_t(3, 3), azmff, dipff, tilff, azmff1, dipff1, &
                   azmff2, dipff2, tilff2, err2, besterr2!find GSLIB angles
   integer :: mm, nn, pp, o_maj, o_min, o_ter! find GSLIB angles.
   real(kind=8), allocatable :: sq_df1(:), sqdfdist1(:), dx1(:), dy1(:), dz1(:)
   !       Variables for MOI cube
   real(kind=8) :: maxx_zero, maxy_zero, maxz_zero, range_x, range_y, range_z, &
                   cube_side, midx, midy, midz, lower_x, lower_y, lower_z, upper_x, upper_y, upper_z
   real(kind=8), allocatable :: x_tmp(:), y_tmp(:), z_tmp(:), v_tmp(:), w_tmp(:)
   real, allocatable :: geolc(:)
   !Other internal variables
   integer :: i, j, k, lin, MAXDAT, nvari, test, nd, strattrans
   integer :: size_maj_arr, size_min_arr, size_ter_arr
   logical :: testfl, scaling
   real(kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax, av, ss, dist, d_x, d_y, d_z
   real:: var(512), vrt
   ! Variables for distance standardization
   real(kind=8) :: avgmin, avgavg, avgmax, zavr, ssavr, ddhlavr, lavr, zavrs, lavrs, rt_dh
   ! Angle related
   integer, allocatable :: npoints(:)
   real(kind=8) :: az_major, az_minor, az_ter, dip_major, dip_minor, dip_ter, &
                   azmf, dipf, tilf
   real(kind=8), allocatable :: pointazm(:), pointdip(:), pointsep(:)
   ! Experimental variogram calc
   integer :: il, ih, nxc, nyc, nzc, iv, l
   logical :: usexdir, useydir, usezdir, debug, usemaj, usemin, useter
   real(kind=8) :: deltax, deltay, deltaz, dims, xcs, ycs, zcs, avgsiz, &
                   debug_onscr, dxv, dyv, dzv, dmaj, dmin, dter
   real(kind=8), parameter :: MAXCOVS = 4
   real(kind=8), allocatable, dimension(:) :: &
      lvmajnum, lvmajcov, lvmajmu, lvmajmuh, lvmajsu, lvmajsuh, lvmajdis, lvmajwt, &
      lvminnum, lvmincov, lvminmu, lvminmuh, lvminsu, lvminsuh, lvmindis, lvminwt, &
      lvternum, lvtercov, lvtermu, lvtermuh, lvtersu, lvtersuh, lvterdis, lvterwt
   integer :: nlaghorz, nlagvert
   real(kind=8):: majspc, minspc, terspc
   ! Temporary arrays in case of attempting to improve the variogram
   real(kind=8), allocatable, dimension(:) :: &
      tlvmajnum, tlvmajcov, tlvmajmu, tlvmajmuh, tlvmajsu, tlvmajsuh, tlvmajdis, tlvmajwt, &
      tlvminnum, tlvmincov, tlvminmu, tlvminmuh, tlvminsu, tlvminsuh, tlvmindis, tlvminwt, &
      tlvternum, tlvtercov, tlvtermu, tlvtermuh, tlvtersu, tlvtersuh, tlvterdis, tlvterwt
   ! Anisotropy ratios to evaluate variograms
   real(kind=8):: optmajmin1, optmajter1, optmajmin2, optmajter2, difratio1, difratio2, &
                  difratio3, difratio4
   real(kind=8), allocatable :: dmajco(:), dminco(:), dterco(:)
   !   Optimizing the ranges and lag selection
   integer :: nmaj_lag, nmin_lag, nter_lag
   real(kind=8) :: dmajmin, dmajmax, dminmin, dminmax, dtermin, dtermax, &
                   dmajlag, dminlag, dterlag
   real(kind=8) :: minx_zero, miny_zero, minz_zero
   real(kind=8) :: rt_lag, avgmin_anish, avgavg_anish, avgmax_anish, zavr_anis, &
                   max_maj_dist, max_min_dist, max_ter_dist, dmaj_ref, dmin_ref, dter_ref, &
                   avgmin_vv, avgavg_vv, avgmax_vv
   real(kind=8), allocatable::bin_edges(:), bin(:)
   !   Variogram values selection to model
   integer :: nvarmaj, nvarmin, nvarter
   ! Optimization function
   integer :: it(3), opt_it(3), it_(3, 2), minpairs, nst, ist, idum
   integer :: maxiter, rseed, nvargs
   integer, allocatable :: varnumpairs(:), varindex(:)
   integer, dimension(3000) :: varind
   logical :: npairswt, invdistwt, fixhmaxvert, fixhminhmax
   real(kind=8), dimension(3, 3, 3) :: rotmat
   real(kind=8) :: silln, anismax, anismin, hmaxvert, hminhmax, opt_c0, c0
   real(kind=8), dimension(3, 2) :: cc_, azm, dip, tilt, ahmax, ahmin, avert
   real(kind=8), dimension(2) :: sillf, c0_
   real(kind=8), allocatable :: varlagdist(:), varvalue(:), &
                                varazm(:), vardip(:), varazm_ref(:), vardip_ref(:), varwts(:)
   real(kind=8), dimension(3) :: opt_cc, opt_azm, opt_dip, opt_tilt, &
                                 opt_ahmax, opt_ahmin, opt_avert, cc, aa1, aa2, aa3
   real(kind=8), dimension(3000) :: varmodelsph_, lagdist_, model_varazm, &
                                    model_vardip
   !****************************************************************************
   ! Start
   !****************************************************************************
   !********* Read parameter file ******************************
   write (*, 999) VERSION
999 format(/' Autovar Version: ', f7.3/)

   ! Get the name of the par file - try default name if no input:
   do i = 1, 512
      str(i:i) = ' '
   end do

   call getarg(1, str)
   if (str(1:1) .eq. ' ') then
      write (*, *) 'Which parameter file do you want to use?'
      read (*, '(a)') str
   end if
   if (str(1:1) .eq. ' ') str(1:20) = 'Autovar.par        '
   inquire (file=str, exist=testfl)
   if (.not. testfl) then
      write (*, *) 'ERROR - the parameter file does not exist,'
      write (*, *) '        check for the file and try again  '
      write (*, *)
      if (str(1:20) .eq. 'Autovar.par        ') then
         write (*, *) '        creating a blank parameter file'
         call makepar
         write (*, *)
      end if
      stop
   end if
   lin = 1
   open (lin, file=str, status='OLD')

   ! Find Start of parameters
1  read (lin, '(a4)', end=98) str(1:4)
   if (str(1:4) .ne. 'STAR') go to 1

   ! Read input of parameters
   read (lin, '(a512)', err=98) datafl
   call chknam(datafl, 512)
   write (*, *) ' data file = ', datafl(1:40)

   read (lin, *, err=98) xcol, ycol, zcol, vcol, idcol
   write (*, *) ' columns = ', xcol, ycol, zcol, vcol, idcol
   ! Check if Z column is provided, X and Y must be provided
   if (xcol .eq. 0) then
      write (*, *) 'X column is necessary'
      go to 99
   end if
   if (ycol .eq. 0) then
      write (*, *) 'Y column is necessary'
      go to 99
   end if
   if (zcol .eq. 0) then
      not_z = .true.
   end if
   if ((xcol .ne. 0) .and. (ycol .ne. 0) .and. (zcol .eq. 0)) then
      write (*, *) ' 2D Dataset'
      dat2d = .true.
      ndir = 2
   else if ((xcol .ne. 0) .and. (ycol .ne. 0) .and. (zcol .ne. 0)) then
      write (*, *) ' 3D Dataset'
      dat3d = .true.
      ndir = 3
   end if
   ! Check if ID column is provided
   if (idcol .eq. 0) then
      not_dhid = .true.
      if (.not. not_z) then
         write (*, *) 'DH ID is necessary to calculate vertical spacing.'
         write (*, *) 'This calculations will be skipped, high chance of getting wrong directions.'
         write (*, *)
      end if
   end if
   read (lin, *, err=98) tmin, tmax
   write (*, *) ' trimming limits = ', tmin, tmax

   read (lin, *, err=98) geolocvar, nanchor
   if (geolocvar .eq. 1) then
      write (*, *) ' considering geolocated variogram  '
      write (*, *) '    Number of anchors = ', nanchor
   else
      write (*, *) ' Not considering geolocated variogram '
      nanchor = 1
   end if

   read (lin, '(a512)', err=98) anxfl
   call chknam(anxfl, 512)
   write (*, *) ' data file = ', anxfl(1:40)

   read (lin, *, err=98) axcol, aycol, azcol
   write (*, '(A,3f18.3)') '  Columns for anchor coordinates = ', axcol, aycol, azcol

   read (lin, *, err=98) varimp
   if (varimp .eq. 1) write (*, '(A)') ' Attemp to improve variogram'

   read (lin, *, err=98) anis, rtmajmin, rtmajter
   if (anis .eq. 1) write (*, '(A,1f18.3)') ' Using anisotropy ratio  = ', rtmajmin, &
      rtmajter

   read (lin, *, err=98) nx, xmn, xsiz
   write (*, '(A,I18,2f18.3)') '  nx, xmn, xsiz = ', nx, xmn, xsiz

   read (lin, *, err=98) ny, ymn, ysiz
   write (*, '(A,I18,2f18.3)') '  ny, ymn, ysiz = ', ny, ymn, ysiz

   read (lin, *, err=98) nz, zmn, zsiz
   write (*, '(A,I18,2f18.3)') '  nz, zmn, zsiz = ', nz, zmn, zsiz

   read (lin, '(a512)', err=98) outfl_mod
   call chknam(outfl_mod, 512)
   if (geolocvar .eq. 0) then
      write (*, *) ' Variogram model output = ', outfl_mod(1:40)
   else
      allocate (modelnames(nanchor))
      outfl_mod = outfl_mod(:index(outfl_mod, '.') - 1)
      write (*, *) ' There will be as many models as anchors = ', nanchor
      do i = 1, nanchor
         if (i .le. 9) fmt = '(I1.1)'
         if (i .gt. 9) fmt = '(I2.1)'
         write (ii, fmt) i
         modelnames(i) = trim(adjustl(outfl_mod))//trim(ii)//".out"
         write (*, *) "   File", i, ":", trim(modelnames(i))
      end do
   end if

   read (lin, '(a512)', err=98) outfl_exp
   call chknam(outfl_exp, 512)

   if (geolocvar .eq. 0) then
      write (*, *) 'Experimental variogram output = ', outfl_exp(1:40)
   else
      allocate (expnames(nanchor))
      outfl_exp = outfl_exp(:index(outfl_exp, '.') - 1)
      write (*, *) ' There will be as many exp. vars. as anchors = ', nanchor
      do i = 1, nanchor
         if (i .le. 9) fmt = '(I1.1)'
         if (i .gt. 9) fmt = '(I2.1)'
         write (ii, fmt) i
         expnames(i) = trim(adjustl(outfl_exp))//trim(ii)//".out"
         write (*, *) "   File", i, ":", trim(expnames(i))
      end do
   end if

   read (lin, '(a512)', err=98) dbgfl
   call chknam(dbgfl, 512)

   if (geolocvar .eq. 0) then
      write (*, *) 'Debug file = ', dbgfl(1:40)
   else
      allocate (dbgnames(nanchor))
      dbgfl = dbgfl(:index(dbgfl, '.') - 1)
      write (*, *) ' There will be as many debug files as anchors = ', nanchor
      do i = 1, nanchor
         if (i .le. 9) fmt = '(I1.1)'
         if (i .gt. 9) fmt = '(I2.1)'
         write (ii, fmt) i
         dbgnames(i) = trim(adjustl(dbgfl))//trim(ii)//".out"
         write (*, *) "   File", i, ":", trim(dbgnames(i))
      end do
   end if

   read (lin, '(a512)', err=98) varmodgs
   call chknam(varmodgs, 512)

   if (geolocvar .eq. 0) then
      write (*, *) 'Text file variogram model= ', varmodgs(1:40)
   else
      allocate (modgsnames(nanchor))
      varmodgs = varmodgs(:index(varmodgs, '.') - 1)
      write (*, *) ' There will be as many vars. models  as anchors = ', nanchor
      do i = 1, nanchor
         if (i .le. 9) fmt = '(I1.1)'
         if (i .gt. 9) fmt = '(I2.1)'
         write (ii, fmt) i
         modgsnames(i) = trim(adjustl(varmodgs))//trim(ii)//".out"
         write (*, '(a,i5,a,a)') "    File", i, ":", trim(modgsnames(i))
      end do
   end if

   read (lin, '(a512)', err=98) wtsfl
   call chknam(wtsfl, 512)

   if (geolocvar .eq. 0) then
      write (*, *) 'File for weights= ', wtsfl(1:40)
   else
      allocate (wtsflnames(nanchor))
      wtsfl = wtsfl(:index(wtsfl, '.') - 1)
      write (*, *) ' There will be as many weight files as anchors = ', nanchor
      do i = 1, nanchor
         if (i .le. 9) fmt = '(I1.1)'
         if (i .gt. 9) fmt = '(I2.1)'
         write (ii, fmt) i
         wtsflnames(i) = trim(adjustl(wtsfl))//trim(ii)//".out"
         write (*, *) "   File", i, ":", trim(wtsflnames(i))
      end do
   end if

   read (lin, *, err=98) debug_onscr
   if (debug_onscr .eq. 1) then
      write (*, *) ' debug on screen'
      debug = .true.
   else
      write (*, *) ' not debug on screen'
      debug = .false.
   end if

   read (lin, *, err=98) dep_type
   if (dep_type .eq. 1) then
      write (*, *) ' Disseminated deposited type'
   else if (dep_type .eq. 2) then
      write (*, *) ' Stratigraphic deposited type'
   else
      write (*, *) ' Working as disseminated deposit type'
      dep_type = 1
   end if
   write (*, *) ''
   write (*, *) ''

   !
   ! Only keep data within the grid limits - find grid limits:
   !
   xmin = xmn - xsiz/2.0
   xmax = xmn + (real(nx) + 0.5)*xsiz
   ymin = ymn - ysiz/2.0
   ymax = ymn + (real(ny) + 0.5)*ysiz
   zmin = zmn - zsiz/2.0
   zmax = zmn + (real(nz) + 0.5)*zsiz
   dom_size = dmax1(xmax, ymax, zmax)

   !
   ! Count the number of input data:
   !
   inquire (file=datafl, exist=testfl)
   if (.not. testfl) then
      write (*, *) 'ERROR data file ', datafl, ' does not exist!'
      stop
   end if
   open (lin, file=datafl, status='OLD')
   read (lin, *)
   read (lin, *, err=99) nvari
   do i = 1, nvari
      read (lin, *)
   end do
   MAXDAT = 0
22 read (lin, *, end=55, err=99) (var(j), j=1, nvari)
   if (xcol .gt. 0) then
      if (var(xcol) .lt. xmin .or. var(xcol) .gt. xmax) go to 22
   end if
   if (ycol .gt. 0) then
      if (var(ycol) .lt. ymin .or. var(ycol) .gt. ymax) go to 22
   end if
   if (zcol .gt. 0) then
      if (var(zcol) .lt. zmin .or. var(zcol) .gt. zmax) go to 22
   end if
   if (var(vcol) .lt. tmin .or. var(vcol) .ge. tmax) go to 22
   MAXDAT = MAXDAT + 1
   go to 22
55 continue
   !
   !
   ! Allocate the needed memory:
   !
   allocate (x(MAXDAT), stat=test)
   if (test .ne. 0) stop 'ERROR: Allocation failed'
   allocate (y(MAXDAT), stat=test)
   if (test .ne. 0) stop 'ERROR: Allocation failed'
   allocate (z(MAXDAT), stat=test)
   if (test .ne. 0) stop 'ERROR: Allocation failed'
   allocate (v(MAXDAT), stat=test)
   if (test .ne. 0) stop 'ERROR: Allocation failed'
   allocate (id(MAXDAT), stat=test)
   if (test .ne. 0) stop 'ERROR: Allocation failed'
   !allocate temporary arrays
   allocate (xt(MAXDAT), yt(MAXDAT), zt(MAXDAT))
   !
   ! Read the input data:
   !
   rewind (lin)
   read (lin, *, err=99)
   read (lin, *, err=99) nvari
   nd = 0
   av = 0.0
   ss = 0.0
   do i = 1, nvari
      read (lin, *, err=99)
   end do
   !
   ! Some tests on column numbers:
   !
   if (xcol .gt. nvari .or. ycol .gt. nvari .or. zcol .gt. nvari .or. vcol .gt. nvari) &
      stop ' column specification error '

   !
   ! Read all the data until the end of the file:
   !
2  read (lin, *, end=3, err=99) (var(j), j=1, nvari)
   vrt = var(vcol)
   if (vrt .lt. tmin .or. vrt .ge. tmax) go to 2
   nd = nd + 1
   if (xcol .le. 0) then
      x(nd) = xmn
   else
      x(nd) = var(xcol)
   end if
   if (ycol .le. 0) then
      y(nd) = ymn
   else
      y(nd) = var(ycol)
   end if
   if (zcol .le. 0) then
      z(nd) = zmn
   else
      z(nd) = var(zcol)
   end if
   if (idcol .le. 0) then
      id(nd) = nd
   else
      id(nd) = var(idcol)
   end if
   if (vcol .le. 0) then
      v(nd) = 0
   else
      v(nd) = var(vcol)
   end if

   av = av + vrt
   ss = ss + vrt*vrt
   go to 2
3  close (lin)
   ! Save initial versions of original data, since they are sorted and modified
   !  more than once
   allocate (iidd(nd))
   allocate (xxx, yyy, zzz, vv, mold=iidd)
   xxx = x
   yyy = y
   zzz = z
   vv = v
   iidd = id
   !Not considering declustering weights
   allocate (wt(nd))
   wt = 1
   !
   ! Report the number of data as a check:
   !
   av = av/max(real(nd), 1.0)
   ss = (ss/max(real(nd), 1.0)) - av*av
   write (*, *) ' Data for Autovar: Variable number ', vcol
   write (*, '(A,1I15)') '  Number   = ', nd
   write (*, '(A,1f18.3)') '   Average  = ', av
   write (*, '(A,1f18.3)') '   Variance = ', ss
   if (nd .lt. 1) stop ' ERROR: there are no data'
   write (*, *)
   write (*, *)

   !
   ! Read anchor file
   !
   allocate (anchx(nanchor))
   allocate (anchy, anchz, mold=anchx)
   if (geolocvar .eq. 1) then
      inquire (file=anxfl, exist=testfl)
      if (.not. testfl) then
         write (*, *) 'ERROR data file ', anxfl, ' does not exist!'
         stop
      end if
      open (lin, file=anxfl, status='OLD')
      read (lin, *)
      read (lin, *, err=99) nvari
      do i = 1, nvari
         read (lin, *)
      end do

      ! read until the end of the file
      i = 0
4     read (lin, *, end=5, err=99) (var(j), j=1, nvari)
      i = i + 1
      if (axcol .le. 0) then
         anchx(i) = (minval(x) + maxval(x))/2
      else
         anchx(i) = var(axcol)
      end if
      if (ycol .le. 0) then
         anchy(i) = (minval(y) + maxval(y))/2
      else
         anchy(i) = var(aycol)
      end if
      if (zcol .le. 0) then
         anchz(i) = (minval(z) + maxval(z))/2
      else
         anchz(i) = var(azcol)
      end if
      go to 4
5     close (lin)
   else
      anchx(1) = (minval(x) + maxval(x))/2
      anchy(1) = (minval(y) + maxval(y))/2
      anchz(1) = (minval(z) + maxval(z))/2
   end if

   !*****************************************************************
   ! Calculate the weights for each anchor if using geolocated option
   !*****************************************************************
   if (geolocvar .eq. 0) then  ! not using geolocated option
      allocate (wtsanch(nd, 1))
      allocate (wtscurrent(nd, 1))! to keep format in loop
      wtsanch = 1
      wtscurrent = 1
   else ! using geolocated option
      allocate (wtsanch(nd, nanchor))
      allocate (wtscurrent(nd, 1))! for separating the weights per iteration
      allocate (sumw(nanchor))
      minx_zero = minval(x(:))
      maxx_zero = maxval(x(:))
      miny_zero = minval(y(:))
      maxy_zero = maxval(y(:))
      range_x = maxx_zero - minx_zero
      range_y = maxy_zero - miny_zero
      radsqd = max(range_x, range_y)
      do i = 1, nanchor
         do j = 1, nd
            sdist = sqrt((anchx(i) - x(j))**2 + &
                         (anchy(i) - y(j))**2 + &
                         (anchz(i) - z(j))**2)
            rdist = real(sdist*sdist)
            wtsanch(j, i) = 1*(exp(-(rdist/2*radsqd)) + EPSLON)*wt(i)
            sumw(i) = sumw(i) + wtsanch(j, i)
         end do
      end do
   end if
   !*************************************************************
   ! Calculate the rotation matrix: MOI, eigenvectors-eigenvalues
   !*************************************************************

   !
   !  Get all the pairs inside 0.20 of the max domain size
   !
   allocate (xyz(nd, 3))
   allocate (xyzref(nd, 3))
   xyz(:, 1) = x(:)
   xyz(:, 2) = y(:)
   xyz(:, 3) = z(:)
   if (ndir .eq. 2) xyzref = xyz ! to ensure that xyzref exist for 2D
   !If only 2D data being used, set Z as constant
   if (ndir .eq. 2) then
      xyz(:, 3) = 1
      z(:) = 1
   end if
   ! Set the default and assume no stratigraphic transformation
   strattrans = 0
   ! Analize data spacing. MOI sensitivity to data geometry.
   call datspac_horz(xyz, nd, idcol, id, avgmin, avgavg, avgmax, debug)
   if (.not. debug) then
      write (*, *) 'Inferred spacing:'
      write (*, *) '  Avg. Min. spacing, Avg. avg. spacing, and Avg. max. data spacing:'
      write (*, '(3f18.3)') avgmin, avgavg, avgmax
      write (*, *) ''
   end if
   minx_zero = minval(x(:))
   maxx_zero = maxval(x(:))
   miny_zero = minval(y(:))
   maxy_zero = maxval(y(:))
   range_x = maxx_zero - minx_zero
   range_y = maxy_zero - miny_zero
   ! Run only for 3d data and only if DHID available
   if (dat3d) then
      if (.not. not_dhid) then
         ! Get the averages:
         if (ndir .eq. 3) then
            call datspac_vert(xyz, nd, idcol, id, zavr, debug, lavr)
            call samplespac_vert(xyz, nd, idcol, id, ssavr, debug, ddhlavr)
            if (.not. debug) then
               write (*, *) '  Average Z spacing '
               write (*, '(1f18.3)') zavr
               write (*, *) ''
               write (*, *) ''
            end if
         end if
         scaling = .false.
         ! Standardize in Z direction for disseminated deposits:
         if (dep_type .eq. 1) then
            xyz(:, 3) = (xyz(:, 3))*((avgmin)/(zavr))
            scaling = .true.
            ! save for reuse in this order
            xyzref = xyz
         end if
      else
         continue
      end if
   end if
   ! Calculate a general variogram first, consider geolocvar=0
   if (geolocvar .eq. 1) then
      nanchororig = nanchor
      geolocvar = 0
      geoloc = .True.! This logical saves the use of the geoloc opt.
      nanchor = 1
   else
      geoloc = .False.
   end if
   ! Do loop for each anchor. If no anchor, nanchor = 1
111 do anx = 1, nanchor
      ! Save referential arrays of original data.
      x(:) = xxx
      y(:) = yyy
      z(:) = zzz
      v = vv
      id = iidd
      !If only 2D data being used, set Z as constant
      if (ndir .eq. 2) then
         xyz(:, 3) = 1
         z(:) = 1
      end if
      !Set the right weights to analyze per anchor
      wtscurrent(:, 1) = wtsanch(:, anx)
      !
      ! For 2D datasets
      !
      if (ndir .eq. 2) then
         call get_max_dist(nd, ndir, xyz, max_dist)
         if (geolocvar .eq. 1) then ! Select points around the current anchor
            ! Select weights corresponding to those points
            nax = 0
            sr = 0.3
            do i = 1, nd
               ddaxx = abs((x(i) - anchx(anx)))
               ddaxy = abs((y(i) - anchy(anx)))
               if ((ddaxx .le. (sr*max_dist)) .and. (ddaxy .le. (sr*max_dist))) &
                  nax = nax + 1
            end do
            allocate (xx(nax))
            allocate (yy, zz, vvax, wtsax, mold=xx)
            j = 0
            do i = 1, nd
               ddaxx = abs((x(i) - anchx(anx)))
               ddaxy = abs((y(i) - anchy(anx)))
               if ((ddaxx .le. (sr*max_dist)) .and. (ddaxy .le. (sr*max_dist))) then
                  j = j + 1
                  xx(j) = x(i)
                  yy(j) = y(i)
                  zz(j) = z(i)
                  vvax(j) = v(i)
                  wtsax(j) = wtscurrent(i, 1)
               end if
            end do
         else ! Take all the points if no anchors being considered
            allocate (xx(nd))
            allocate (yy, zz, vvax, wtsax, mold=xx)
            xx = x
            yy = y
            zz = z
            vvax = v
            wtsax = wtscurrent(i, 1)
            nax = nd
         end if

         call omnipair(nax, max_dist, tmin, tmax, xx, yy, zz, vvax, geolocvar, &
                       wtsax, outputidx, sqdfdist, dx, dy, dz, sq_df, idf, wtsomni, debug)
      end if
      !
      ! For 3D datasets
      !
      if (ndir .eq. 3) then
         call get_max_dist(nd, ndir, xyz, max_dist)
         !First select the points for the anchor
         if (geolocvar .eq. 1) then ! Select points and weights for each anchor
            nax = 0
            sr = 0.3
            do i = 1, nd
               ddaxx = abs((x(i) - anchx(anx)))
               ddaxy = abs((y(i) - anchy(anx)))
               ddaxz = abs((z(i) - anchz(anx)))
               if ((ddaxx .le. (sr*max_dist)) .and. (ddaxy .le. (sr*max_dist)) &
                   .and. (ddaxz .le. (sr*max_dist))) nax = nax + 1
            end do
            allocate (xx(nax))
            allocate (yy, zz, vvax, wtsax, mold=xx)
            j = 0
            do i = 1, nd
               ddaxx = abs((x(i) - anchx(anx)))
               ddaxy = abs((y(i) - anchy(anx)))
               ddaxz = abs((z(i) - anchz(anx)))
               if ((ddaxx .le. (sr*max_dist)) .and. (ddaxy .le. (sr*max_dist)) &
                   .and. (ddaxz .le. (sr*max_dist))) then
                  j = j + 1
                  xx(j) = x(i)
                  yy(j) = y(i)
                  zz(j) = z(i)
                  vvax(j) = v(i)
                  wtsax(j) = wtscurrent(i, 1)
               end if
            end do
            ! Find the range in X and Y for the selected points
            maxx_zero = maxval(xx)
            minx_zero = minval(xx)
            maxy_zero = maxval(yy)
            miny_zero = minval(yy)
            range_x = maxx_zero - minx_zero
            range_y = maxy_zero - miny_zero
         else ! Take all the points if no anchors being considered
            allocate (xx(nd))
            allocate (yy, zz, vvax, wtsax, mold=xx)
            xx = x
            yy = y
            zz = z
            vvax = v
            wtsax = wtscurrent(i, 1)
            ! Take the original scaled values as well
         end if
         !For disseminated deposits calculated direction inside cube
         if (dep_type .eq. 1) then
            if (geolocvar .eq. 1) then  !Reevaluate the Z extends since it's
               ! for data close to the anchors
               minz_zero = minval(zz)
               maxz_zero = maxval(zz)
               range_z = maxz_zero - minz_zero
               cube_side = min(range_x, range_y)
               midx = (maxx_zero + minx_zero)/2
               midy = (maxy_zero + miny_zero)/2
               midz = (maxz_zero + minz_zero)/2
               lower_x = midx - cube_side/2
               lower_y = midy - cube_side/2
               lower_z = midz - cube_side/2
               upper_x = midx + cube_side/2
               upper_y = midy + cube_side/2
               upper_z = midz + cube_side/2
                    !! Count points that satisfy this condition
               nd_tmp = 0
               ! make sure that the number of data is valid
               if (geolocvar .eq. 0) nax = nd
               do i = 1, nax
                  if ((xx(i) .gt. lower_x) .and. (xx(i) .lt. upper_x)) then
                     if ((yy(i) .gt. lower_y) .and. (yy(i) .lt. upper_y)) then
                        if ((zz(i) .gt. lower_z) .and. (zz(i) .lt. upper_z)) then
                           nd_tmp = nd_tmp + 1
                        end if
                     end if
                  end if
               end do
               ! Allocate arrays for points inside this condition
               allocate (w_tmp(nd_tmp))
               allocate (x_tmp, y_tmp, z_tmp, v_tmp, mold=w_tmp)
               ! Fill the arrays
               j = 0
               do i = 1, nax
                  if ((xx(i) .gt. lower_x) .and. (xx(i) .lt. upper_x)) then
                     if ((yy(i) .gt. lower_y) .and. (yy(i) .lt. upper_y)) then
                        if ((zz(i) .gt. lower_z) .and. (zz(i) .lt. upper_z)) then
                           j = j + 1
                           x_tmp(j) = xx(i)
                           y_tmp(j) = yy(i)
                           z_tmp(j) = zz(i)
                           v_tmp(j) = vvax(i)
                           w_tmp(j) = wtsax(i)
                        end if
                     end if
                  end if
               end do
               max_dist = sqrt((maxval(x_tmp) - minval(x_tmp))**2 + &
                               (maxval(y_tmp) - minval(y_tmp))**2 + &
                               (maxval(z_tmp) - minval(z_tmp))**2)
               call omnipair(nd_tmp, max_dist, tmin, tmax, x_tmp, y_tmp, z_tmp, v_tmp, &
                             geolocvar, w_tmp, outputidx, sqdfdist, dx, dy, dz, sq_df, idf, wtsomni, debug)
            end if

            if (geolocvar .eq. 0) then    ! if no geolocated option used
               ! Reevaluate the z extends since it's been scaled
               minz_zero = minval(xyz(:, 3))
               maxz_zero = maxval(xyz(:, 3))
               range_z = maxz_zero - minz_zero
               cube_side = min(range_x, range_y)
               midx = (maxx_zero + minx_zero)/2
               midy = (maxy_zero + miny_zero)/2
               midz = (maxz_zero + minz_zero)/2
               lower_x = midx - cube_side/2
               lower_y = midy - cube_side/2
               lower_z = midz - cube_side/2
               upper_x = midx + cube_side/2
               upper_y = midy + cube_side/2
               upper_z = midz + cube_side/2
                    !! Count points that satisfy this condition
               nd_tmp = 0
               do i = 1, nd
                  if ((xyz(i, 1) .gt. lower_x) .and. (xyz(i, 1) .lt. upper_x)) then
                     if ((xyz(i, 2) .gt. lower_y) .and. (xyz(i, 2) .lt. upper_y)) then
                        if ((xyz(i, 3) .gt. lower_z) .and. (xyz(i, 3) .lt. upper_z)) then
                           nd_tmp = nd_tmp + 1
                        end if
                     end if
                  end if
               end do
               ! Allocate arrays for points inside this condition
               allocate (w_tmp(nd_tmp))
               allocate (x_tmp, y_tmp, z_tmp, v_tmp, mold=w_tmp)
               ! Fill the arrays
               j = 0
               do i = 1, nd
                  if ((xyz(i, 1) .gt. lower_x) .and. (xyz(i, 1) .lt. upper_x)) then
                     if ((xyz(i, 2) .gt. lower_y) .and. (xyz(i, 2) .lt. upper_y)) then
                        if ((xyz(i, 3) .gt. lower_z) .and. (xyz(i, 3) .lt. upper_z)) then
                           j = j + 1
                           x_tmp(j) = xx(i)
                           y_tmp(j) = yy(i)
                           z_tmp(j) = zz(i)
                           v_tmp(j) = vv(i)
                           w_tmp(j) = wtsax(i)
                        end if
                     end if
                  end if
               end do
               max_dist = sqrt((maxval(x_tmp) - minval(x_tmp))**2 + &
                               (maxval(y_tmp) - minval(y_tmp))**2 + &
                               (maxval(z_tmp) - minval(z_tmp))**2)
               call omnipair(nd_tmp, max_dist, tmin, tmax, x_tmp, y_tmp, z_tmp, v_tmp, &
                             geolocvar, w_tmp, outputidx, sqdfdist, dx, dy, dz, sq_df, idf, wtsomni, debug)
            end if

         end if
         !For tabular deposits calculated directions using all pairs
         if (dep_type .eq. 2) then
            call get_max_dist(nd, ndir, xyz, max_dist)
            call omnipair(nd, max_dist, tmin, tmax, x, y, z, v, geolocvar, wtscurrent, &
                          outputidx, sqdfdist, dx, dy, dz, sq_df, idf, wtsomni, debug)
         end if

      end if

      !
      ! Getting all values gt 0 to calculate MOI
      !
      z_count = 0
      do i = 1, outputidx
         if (sq_df(i) .lt. 0.0) then
            sq_df(i) = 0.0
            z_count = z_count + 1
         end if
      end do
      outputidx1 = 0
      do i = 1, outputidx
         if (sq_df(i) .gt. 0.0) then
            outputidx1 = outputidx1 + 1
         end if
      end do

      ! Allocate for values above 0 to speed calculations.
      allocate (sq_df1(outputidx1))
      allocate (sqdfdist1(outputidx1))
      allocate (dx1(outputidx1))
      allocate (dy1(outputidx1))
      allocate (dz1(outputidx1))
      allocate (wtsomni1(outputidx1))
      ! Fill the arrays
      j = 0
      do i = 1, outputidx
         if (sq_df(i) .gt. 0.0) then
            j = j + 1
            sq_df1(j) = sq_df(i)
            sqdfdist1(j) = sqdfdist(i)
            dx1(j) = dx(i)
            dy1(j) = dy(i)
            dz1(j) = dz(i)
            wtsomni1(j) = wtsomni(i)
         else
            continue
         end if
      end do
      !
      ! Calculate MOI
      !
      if (debug) then
         write (*, *) ' Calculating MOI matrix... '
      end if
      call moi(sq_df1, sqdfdist1, outputidx1, tmin, tmax, dom_size, dx1, dy1, dz1, &
               wtsomni1, moimat, eigvec, eigval, perc_99val, debug)

      if (debug) then
         write (*, *) ('Eigenvectors: ')
         do i = 1, 3
            write (*, '(3f18.3)') (eigvec(i, k), k=1, 3)
         end do
         write (*, *) ('Eigenvalues: ')
         write (*, '(3f18.3)') (eigval(i), i=1, 3)
      end if
      ! Save referential order:
      do i = 1, 3
         dblerotindexx(i) = i
      end do

      ! Order the eigenvectors according to the eigenvalues in ascendent order:
      call dblemodsortem(eigval(:), 3, 1, dblerotindexx)

      do i = 1, 3
         introtindexx(i) = int(dblerotindexx(i))
      end do

        !! Definitive order: having the theory as basis
      if (ndir .eq. 3 .and. strattrans .eq. 0) then
         eigveco(1, 1) = eigvec(1, dblerotindexx(1))! xmaj
         eigveco(2, 1) = eigvec(2, dblerotindexx(1))! ymaj
         eigveco(3, 1) = eigvec(3, dblerotindexx(1))! zmaj
         eigveco(1, 2) = eigvec(1, dblerotindexx(2))! xmin
         eigveco(2, 2) = eigvec(2, dblerotindexx(2))! ymin
         eigveco(3, 2) = eigvec(3, dblerotindexx(2))! zmin
         eigveco(1, 3) = eigvec(1, dblerotindexx(3))! xter
         eigveco(2, 3) = eigvec(2, dblerotindexx(3))! yter
         eigveco(3, 3) = eigvec(3, dblerotindexx(3))! zter
         if (debug) then
            write (*, *) ('Ordered Eigenvectors: ')
            do i = 1, 3
               write (*, '(3f18.3)') (eigveco(i, k), k=1, 3)
            end do
         end if
      end if

      if (ndir .eq. 2) then
         eigveco(1, 1) = eigvec(1, dblerotindexx(1))! xmaj
         eigveco(2, 1) = eigvec(2, dblerotindexx(1))! ymaj
         eigveco(3, 1) = 0! zmaj
         eigveco(1, 2) = eigvec(1, dblerotindexx(2))! xmin
         eigveco(2, 2) = eigvec(2, dblerotindexx(2))! ymin
         eigveco(3, 2) = 0! zmin
         eigveco(1, 3) = eigvec(1, dblerotindexx(3))! xter
         eigveco(2, 3) = eigvec(2, dblerotindexx(3))! yter
         eigveco(3, 3) = 1! zter
         ! Get azimuth and dip
         a1 = PI/2 - atan(eigveco(2, 1)/eigveco(1, 1)) !Az
         a2 = 0 !Dip
         a3 = 0 !Tilt
         a1 = a1*RAD2DEG
         a2 = a2*RAD2DEG
         a3 = a3*RAD2DEG
         if (debug) then
            write (*, *) ('Ordered Eigenvectors: ')
            do i = 1, 3
               write (*, '(3f18.3)') (eigveco(i, k), k=1, 3)
            end do
            write (*, *)
            write (*, *) '2D dataset inferred angles:'
            write (*, '(3f15.3)') a1, a2, a3
            write (*, *)
         end if

      end if

      ! Referential Az and Dip of each vector of the principal directions
      az_major = atan2(eigveco(1, 1), eigveco(2, 1))*RAD2DEG !Az
      dip_major = atan2(eigveco(3, 1), sqrt(eigveco(1, 1)**2 + eigveco(2, 1)**2)) &
                  *RAD2DEG
      az_minor = atan2(eigveco(1, 2), eigveco(2, 2))*RAD2DEG !Az
      dip_minor = atan2(eigveco(3, 2), sqrt(eigveco(1, 2)**2 + eigveco(2, 2)**2)) &
                  *RAD2DEG
      az_ter = atan2(eigveco(1, 3), eigveco(2, 3))*RAD2DEG !Az
      dip_ter = atan2(eigveco(3, 3), sqrt(eigveco(1, 3)**2 + eigveco(2, 3)**2)) &
                *RAD2DEG

      !*****************************************************************
      ! Calculate the variogram(correlogram in the inferred directions)
      !*****************************************************************
      ! Assign the rotation matrix values - transpose them:

      lrotmat(1, 1) = eigveco(1, 1)
      lrotmat(1, 2) = eigveco(2, 1)
      lrotmat(1, 3) = eigveco(3, 1)
      lrotmat(2, 1) = eigveco(1, 2)
      lrotmat(2, 2) = eigveco(2, 2)
      lrotmat(2, 3) = eigveco(3, 2)
      lrotmat(3, 1) = eigveco(1, 3)
      lrotmat(3, 2) = eigveco(2, 3)
      lrotmat(3, 3) = eigveco(3, 3)

      ! Ensure the right set of eigenvectors, according to GSLIB angle convention
      call get_octant(eigveco, o_maj, o_min, o_ter, right_eiveco)

      ! Get the angles from right_eiveco:
      lrotmatrob(1, 1) = right_eiveco(1, 1)
      lrotmatrob(1, 2) = right_eiveco(2, 1)
      lrotmatrob(1, 3) = right_eiveco(3, 1)
      lrotmatrob(2, 1) = right_eiveco(1, 2)
      lrotmatrob(2, 2) = right_eiveco(2, 2)
      lrotmatrob(2, 3) = right_eiveco(3, 2)
      lrotmatrob(3, 1) = right_eiveco(1, 3)
      lrotmatrob(3, 2) = right_eiveco(2, 3)
      lrotmatrob(3, 3) = right_eiveco(3, 3)

      if (lrotmatrob(1, 3) .eq. 1) then ! if dip=90
         azmf = atan(-1*lrotmatrob(2, 1)/lrotmatrob(3, 1))*RAD2DEG !Az
         dipf = 90
         tilf = 0
      else ! if dip != 90
         bb = asin(-lrotmatrob(1, 3))
         aa = atan(lrotmatrob(1, 2)/lrotmatrob(1, 1))
         tt = atan(-1*lrotmatrob(2, 3)/lrotmatrob(3, 3))
         dipf = -bb*RAD2DEG !Dip
         azmf = 90 - aa*RAD2DEG !Az
         tilf = -tt*RAD2DEG! Tilt
      end if
      ! Do an optimization check to find the best angles possible
      ! Get the reference error with the previously calculated angles.
      err2 = 0.0
      do i = 1, 3
         do j = 1, 3
            err2 = err2 + abs( &
                   (eigveco(i, j))**2 - (right_eiveco(i, j)**2) &
                   )
         end do
      end do
      rotmat_t = 0.0
      revrotmat_t = 0.0
      besterr = err2
      ! Get the new error modifying parameters.
      do mm = 1, 2
         do nn = 1, 2
            do pp = 1, 2
               azm_opt = azmf + (mm - 1)*180
               if (azm_opt .gt. 360) azm_opt = azm_opt - 360
               dip_opt = dipf*(-1)**nn
               tilt_opt = tilf*(-1)**pp
               call setrotmat(azm_opt, dip_opt, tilt_opt, rotmat_t, revrotmat_t)
               err = 0.0
               do i = 1, 3
                  do j = 1, 3
                     err = err + abs( &
                           (eigveco(i, j))**2 - (rotmat_t(i, j)**2) &
                           )
                  end do
               end do

               if (err .lt. besterr) then
                  besterr = err
                  azmf = azm_opt
                  dipf = dip_opt
                  tilf = tilt_opt
               end if
            end do
         end do
      end do

      ! The following is an extract of the klp program
      ! Set the data extends considering dimensions
      !
      ! Determine the extents of the data provided
      !
      il = int(1.0 + 0.04*real(nd))
      ih = nd - il
      usexdir = .true.
      call dblemodsortem(x, nd, 4, y, z, v, id)
      deltax = x(ih) - x(il)
      if (deltax .le. EPSLON) usexdir = .false.
      useydir = .true.
      call dblemodsortem(y, nd, 4, x, z, v, id)
      deltay = y(ih) - y(il)
      if (deltay .le. EPSLON) useydir = .false.
      usezdir = .true.
      call dblemodsortem(z, nd, 4, x, y, v, id)
      deltaz = z(ih) - z(il)
      if (deltaz .le. EPSLON) usezdir = .false.

      !
      !     The size of the cells in the covariance table is important.  The
      !     calculation is based on an assumption of equal spacing in all
      !     directions.  This is only referential
      !
      xcs = 0.0; ycs = 0.0; zcs = 0.0; dims = 0.0
      if (usexdir) dims = dims + 1
      if (useydir) dims = dims + 1
      if (usezdir) dims = dims + 1
      if (usexdir) xcs = deltax/(real(nd)**(1.0/dims))
      if (useydir) ycs = deltay/(real(nd)**(1.0/dims))
      if (usezdir) zcs = deltaz/(real(nd)**(1.0/dims))
      ! check if the spacing is similar to the spacing of the samples
      ! specially for the direction of the downhole
      if (xcs .gt. avgmin) xcs = avgmin
      if (ycs .gt. avgmin) ycs = avgmin
      if (zcs .gt. zavr) zcs = zavr
      if (usexdir .and. useydir) then
         if ((xcs/ycs) .ge. 0.5 .and. (xcs/ycs) .le. 2.0) then
            avgsiz = 0.5*(xcs + ycs)
            xcs = avgsiz
            ycs = avgsiz
         end if
      end if
      nxc = 1; nyc = 1; nzc = 1
      if (usexdir) nxc = 2*MAXCOVS + 1
      if (useydir) nyc = 2*MAXCOVS + 1
      if (usezdir) nzc = 2*MAXCOVS + 1
      if (debug) then
         write (*, *)
         write (*, '(A,2f18.3)') '  X ', deltax, xcs
         write (*, '(A,2f18.3)') '  Y ', deltay, ycs
         write (*, '(A,2f18.3)') '  Z ', deltaz, zcs
      end if

      !
      ! Components of the principal directions in the major and establish what
      ! directions to optimize
      !
      usemaj = .false.
      usemin = .false.
      useter = .false.
      dxv = deltax
      dyv = 0.0
      dzv = 0.0
      dmaj = lrotmat(1, 1)*dxv + lrotmat(1, 2)*dyv + lrotmat(1, 3)*dzv
      dmin = lrotmat(2, 1)*dxv + lrotmat(2, 2)*dyv + lrotmat(2, 3)*dzv
      dter = lrotmat(3, 1)*dxv + lrotmat(3, 2)*dyv + lrotmat(3, 3)*dzv
      if (usexdir) then
         if (abs(dmaj) .gt. EPSLON) usemaj = .true.
         if (abs(dmin) .gt. EPSLON) usemin = .true.
         if (abs(dter) .gt. EPSLON) useter = .true.
      end if
      dxv = 0.0
      dyv = deltay
      dzv = 0.0
      dmaj = lrotmat(1, 1)*dxv + lrotmat(1, 2)*dyv + lrotmat(1, 3)*dzv
      dmin = lrotmat(2, 1)*dxv + lrotmat(2, 2)*dyv + lrotmat(2, 3)*dzv
      dter = lrotmat(3, 1)*dxv + lrotmat(3, 2)*dyv + lrotmat(3, 3)*dzv
      if (useydir) then
         if (abs(dmaj) .gt. EPSLON) usemaj = .true.
         if (abs(dmin) .gt. EPSLON) usemin = .true.
         if (abs(dter) .gt. EPSLON) useter = .true.
      end if
      dxv = 0.0
      dyv = 0.0
      dzv = deltaz
      dmaj = lrotmat(1, 1)*dxv + lrotmat(1, 2)*dyv + lrotmat(1, 3)*dzv
      dmin = lrotmat(2, 1)*dxv + lrotmat(2, 2)*dyv + lrotmat(2, 3)*dzv
      dter = lrotmat(3, 1)*dxv + lrotmat(3, 2)*dyv + lrotmat(3, 3)*dzv
      if (usezdir) then
         if (abs(dmaj) .gt. EPSLON) usemaj = .true.
         if (abs(dmin) .gt. EPSLON) usemin = .true.
         if (abs(dter) .gt. EPSLON) useter = .true.
      end if

      ! *************************************************************
      ! Calculate experimental directional variograms
      ! *************************************************************

      !
      ! Set variogram calculation parameters based on the extent of the data
      ! in the principal directions
      !
      dxv = 0.0; dyv = 0.0; dzv = 0.0
      dmajmin = 1.0e21; dmajmax = -1.0e21
      dminmin = 1.0e21; dminmax = -1.0e21
      dtermin = 1.0e21; dtermax = -1.0e21
      ! Allocate arrays to storage rotated coordinates
      allocate (dmajco(nd), dminco(nd), dterco(nd))
      do i = 1, nd
         if (usexdir) dxv = x(i)
         if (useydir) dyv = y(i)
         if (usezdir) dzv = z(i)
         dmaj = lrotmat(1, 1)*dxv + lrotmat(1, 2)*dyv + lrotmat(1, 3)*dzv
         dmin = lrotmat(2, 1)*dxv + lrotmat(2, 2)*dyv + lrotmat(2, 3)*dzv
         dter = lrotmat(3, 1)*dxv + lrotmat(3, 2)*dyv + lrotmat(3, 3)*dzv
         ! Fill arrays with the new coordinates
         dmajco(i) = dmaj
         dminco(i) = dmin
         dterco(i) = dter
         if (dmaj .lt. dmajmin) dmajmin = dmaj
         if (dmaj .gt. dmajmax) dmajmax = dmaj
         if (dmin .lt. dminmin) dminmin = dmin
         if (dmin .gt. dminmax) dminmax = dmin
         if (dter .lt. dtermin) dtermin = dter
         if (dter .gt. dtermax) dtermax = dter
      end do
      dmaj_ref = dmajmax - dmajmin
      dmin_ref = dminmax - dminmin
      dter_ref = dtermax - dtermin
      if (debug) write (*, *) ' '
      if (debug) write (*, *) 'Referential maximum search using rotated coordinates:'
      if (debug) write (*, '(3f18.3)') dmaj_ref, dmin_ref, dter_ref
      if (debug) write (*, *) ' '
      !************************************************************************
      ! Lag analysis:
      ! 1. Find average distances in the 3 directions of continuity = lag distance

      ! Allocate all anisotropic distances in one array and fill it
      allocate (xyz_anis(nd, 3))
      xyz_anis(:, 1) = dmajco
      xyz_anis(:, 2) = dminco
      xyz_anis(:, 3) = dterco
      if (debug) write (*, *) 'Calculating anisotropic lag distances'

      call datspac_horz(xyz_anis, nd, idcol, id, avgmin_anish, avgavg_anish, avgmax_anish, &
                        debug)

      ! Run only for 3d data and only if DHID available
      if (dat3d) then
         if (.not. not_dhid) then
            ! Get the averages:
            if (ndir .eq. 3) then
               call samplespac_vert(xyz_anis, nd, idcol, id, zavr_anis, debug, lavr)
            end if
         end if
      end if
      ! Subroutine to find consistent lag spacings:
      ! Define number of lags from: reference_ext/spacing
      nlaghorz = max(deltax, deltay)/(avgmin)
      if (ndir .eq. 3) nlagvert = deltaz/(zavr)
      if (ndir .eq. 2) nlagvert = 1
      if (ndir .eq. 3) then
         zavr_anis = ssavr
         call spac3d(nd, xyz_anis, avgmin, zavr_anis, nlaghorz, nlagvert, dep_type, &
                     majspc, minspc, terspc)
      end if
      ! lag distances:
      dmajlag = max(EPSLON, abs(avgavg_anish))
      dminlag = max(EPSLON, abs(avgmin_anish))
      ! set conditions for 3D and 2D datasets
      if (ndir .eq. 3) then
         dterlag = max(EPSLON, abs(terspc))
      end if
      ! If only 2D data: set generic parameters for tertiary direction
      if (ndir .eq. 2) dterlag = 1
      if (ndir .eq. 2) zavr_anis = 1
      ! 2. Calculate the max distance in each direction
      ! Calculate the half of distance on the direction of the rotated coordinates:
      if (dep_type .eq. 1) then
         max_maj_dist = 0.5*abs(maxval(dmajco) - minval(dmajco))
         max_min_dist = 0.5*abs(maxval(dminco) - minval(dminco))
         if (ndir .eq. 3) then
            max_ter_dist = 0.5*abs(maxval(dterco) - minval(dterco))
         else
            max_ter_dist = 0.0
         end if
      end if
      if (dep_type .eq. 2) then
         max_maj_dist = 0.5*abs(maxval(dmajco) - minval(dmajco))
         max_min_dist = 0.5*abs(maxval(dminco) - minval(dminco))
         if (ndir .eq. 3) then
            max_ter_dist = 0.5*abs(maxval(dterco) - minval(dterco))
         else
            max_ter_dist = 0.0
         end if
      end if

      ! 2.1. If geolocated variogram: set max distance to 1/3 of domain size
      if (geolocvar .eq. 1) then
         if (max_maj_dist .gt. 0.3*max_dist) max_maj_dist = 0.3*max_dist
         if (max_min_dist .gt. 0.3*max_dist) max_min_dist = 0.3*max_dist
         ! Allow more extension for 3D datasets
         if (ndir .eq. 3) then
            if (max_maj_dist .gt. 0.4*max_dist) max_maj_dist = 0.4*max_dist
            if (max_min_dist .gt. 0.4*max_dist) max_min_dist = 0.4*max_dist
            if (max_ter_dist .gt. 0.4*max_dist) max_ter_dist = 0.4*max_dist
         end if
      end if

      ! 2.2. Check if order relations are respected, otherwise enforce them
      if (max_min_dist .gt. max_maj_dist) max_min_dist = max_maj_dist
      if (max_ter_dist .gt. max_min_dist) max_ter_dist = max_min_dist
      dmaj = max_maj_dist
      dmin = max_min_dist
      dter = max_ter_dist

      if (ndir .eq. 2) dter = 0.0
      if (debug) write (*, *) 'Maximum search by direction - Rotated Coord.:'
      if (debug) write (*, '(3f18.3)') dmaj, dmin, dter

      ! 3. Find the number of lags needed per direction and show lag and max dist
      nmaj_lag = int(dmaj/dmajlag) + 1
      nmin_lag = int(dmin/dminlag) + 1
      nter_lag = int(dter/dterlag) + 1
      ! Show results on screen for user:
      if (debug) then
         write (*, *)
         write (*, '(A,3f15.3)') ' Extent ', dmaj, dmin, dter
         write (*, '(A,3f15.3)') ' Lags   ', dmajlag, dminlag, dterlag
         write (*, *) ' '
      end if
      ! 4. Allocate the arrays according to the number of lags+1
      !major direction
      allocate (lvmajnum(nmaj_lag + 1))
      allocate (lvmajcov, lvmajmu, lvmajmuh, lvmajsu, lvmajsuh, &
                lvmajdis, lvmajwt, mold=lvmajnum)
      !minor direction
      allocate (lvminnum(nmin_lag + 1))
      allocate (lvmincov, lvminmu, lvminmuh, lvminsu, lvminsuh, &
                lvmindis, lvminwt, mold=lvminnum)
      ! tertiary direction
      allocate (lvternum(nter_lag + 1))
      allocate (lvtercov, lvtermu, lvtermuh, lvtersu, lvtersuh, &
                lvterdis, lvterwt, mold=lvternum)

      ! **************************************************************
      ! Loop over all pairs of data calculating experimental variogram
      !
      if (dep_type .eq. 1) then
         x = xxx
         y = yyy
         z = zzz
         v = vv
         id = iidd
      end if
      call varexp(x, y, z, v, lrotmat, nd, usemaj, usemin, useter, nmaj_lag, nmin_lag, nter_lag, &
                  dmajlag, dminlag, dterlag, zavr_anis, usexdir, useydir, usezdir, debug, dep_type, wtscurrent, &
                  lvmajnum, lvmajcov, lvmajmu, lvmajmuh, lvmajsu, lvmajsuh, lvmajdis, lvmajwt, &
                  lvminnum, lvmincov, lvminmu, lvminmuh, lvminsu, lvminsuh, lvmindis, lvminwt, &
                  lvternum, lvtercov, lvtermu, lvtermuh, lvtersu, lvtersuh, lvterdis, lvterwt, &
                  anchx(anx), anchy(anx), anchz(anx), anx, geolocvar, strattrans)
      !
      ! Ensure that the min value of each direction accomplish the order
      ! relations, otherwise ignore to model
      if (ndir .eq. 3) then
         if (lvmindis(1) .lt. lvterdis(1)) then
            lvterwt(1) = 0.0
            lvterwt(1) = 0.0
         end if
      end if
      if (lvmajdis(1) .lt. lvmindis(1)) then
         lvminwt(1) = 0.0
         lvmajwt(1) = 0.0
      end if

      !***************************************************************
      !Optimization function to do the autofit
      !**************************************************************
      !
      ! Since the points are already rotated:
      ! set rotation matrices to identity matrix for fitting
      !
      ! define nst to 3 as default:
      nst = 3
      do ist = 1, nst
         do i = 1, 3
            do j = 1, 3
               rotmat(ist, i, j) = 0.0
            end do
         end do
         ! rotate 90 degrees to be in gslib convention
         ! of cova3.
         rotmat(ist, 1, 2) = -1.0
         rotmat(ist, 2, 1) = 1.0
         rotmat(ist, 3, 3) = 1.0
      end do
      !
      ! Make sure variogram is standard and with long enough range
      !
      c0 = 0! this is the nugget effect, it is hardcoded now (could be set by the user)

      cc = 1
      silln = c0

      do ist = 1, nst
         silln = silln + cc(ist)
      end do
      c0 = c0*1.0/silln
      do ist = 1, nst
         cc(ist) = cc(ist)*1.0/silln
      end do
      !
      !   Incorporate varmodel modified to do the autofit
      !
      ! fix it as spherical
      it = 1
      it_(:, :) = 1
      cc_(:, 1) = 1.000000000000000D-002
      cc_(:, 2) = 1.000000000000000D+021

      c0_(1) = 0.000000000000000D+000
      c0_(2) = 0.000000000000000D+000

      azm(:, 1) = -180.0
      azm(:, 2) = 180.0
      dip(:, 1) = -90.0
      dip(:, 2) = 90.0
      tilt(:, 1) = -90.0
      tilt(:, 2) = 90.0

      ahmax(:, 1) = 0.75*minval(lvmajdis)
      ahmax(:, 2) = 4*maxval(lvmajdis)
      ahmin(:, 1) = 0.75*minval(lvmindis)
      ahmin(:, 2) = 4*maxval(lvmindis)
      avert(:, 1) = 0.75*minval(lvterdis)
      avert(:, 2) = 4*maxval(lvterdis)
      ! Allows for more extension for tabular deposits
      if (dep_type .eq. 2) then
         ahmax(:, 1) = 0.75*minval(lvmajdis)
         ahmax(:, 2) = 20*maxval(lvmajdis)
         ahmin(:, 1) = 0.75*minval(lvmindis)
         ahmin(:, 2) = 20*maxval(lvmindis)
         avert(:, 1) = 0.75*minval(lvterdis)
         avert(:, 2) = 4*maxval(lvterdis)
      end if

      if (.not. useter) avert = 0
      sillf(1) = 1
      sillf(2) = 1

      maxiter = 2000000
      rseed = 53783
      minpairs = 10
      npairswt = .false.
      invdistwt = .false.

      fixhmaxvert = .false.
      hmaxvert = 1.0

      fixhminhmax = .false.
      hminhmax = 1.0
      ! add up to create one single table for the variogram model(include 3 dirs)
      nvarmaj = 0; nvarmin = 0; nvarter = 0
      nvarmaj = nmaj_lag
      nvarmin = nmin_lag
      if (useter) nvarter = nter_lag
      nvargs = nvarmaj + nvarmin + nvarter
      ! allocate the needed arrays,
      allocate (varlagdist(nvargs))
      allocate (varvalue(nvargs))
      allocate (varnumpairs(nvargs))
      allocate (varwts(nvargs))
      allocate (varazm(nvargs))
      allocate (vardip(nvargs))

      varlagdist(1:nvarmaj) = lvmajdis
      varlagdist(nvarmaj + 1:nvarmaj + nvarmin) = lvmindis
      if (useter) varlagdist(nvarmaj + nvarmin + 1:nvargs) = lvterdis

      varvalue(1:nvarmaj) = lvmajcov
      varvalue(nvarmaj + 1:nvarmaj + nvarmin) = lvmincov
      if (useter) varvalue(nvarmaj + nvarmin + 1:nvargs) = lvtercov

      varnumpairs(1:nvarmaj) = int(lvmajnum)
      varnumpairs(nvarmaj + 1:nvarmaj + nvarmin) = int(lvminnum)
      if (useter) varnumpairs(nvarmaj + nvarmin + 1:nvargs) = int(lvternum)

      varwts(1:nvarmaj) = lvmajwt
      varwts(nvarmaj + 1:nvarmaj + nvarmin) = lvminwt
      if (useter) varwts(nvarmaj + nvarmin + 1:nvargs) = lvterwt

      ! The next parameters are hard coded to make the code follow the rotmat
      ! convention

      varazm(1:nvarmaj) = 0
      varazm(nvarmaj + 1:nvarmaj + nvarmin) = 90
      if (useter) varazm(nvarmaj + nvarmin + 1:nvargs) = 0
      vardip(1:nvarmaj) = 0
      vardip(nvarmaj + 1:nvarmaj + nvarmin) = 0
      if (useter) vardip(nvarmaj + nvarmin + 1:nvargs) = -90
      ! Optimization function for variogram modelling
      call varmodelfit(nst, c0_, it_, cc_, azm, dip, tilt, ahmax, ahmin, avert, sillf, &
                       maxiter, rseed, minpairs, npairswt, invdistwt, fixhmaxvert, hmaxvert, &
                       fixhminhmax, hminhmax, nvargs, varlagdist, varvalue, varnumpairs, varwts, &
                       varazm, vardip, opt_c0, opt_it, opt_cc, opt_azm, opt_dip, opt_tilt, &
                       opt_ahmax, opt_ahmin, opt_avert)

      if (debug) then
         write (*, *) 'Optimized variogram model'
         write (*, '(1f15.3)') opt_cc(1)
         write (*, '(3f15.3)') opt_ahmax(1), opt_ahmin(1), opt_avert(1)
         write (*, '(1f15.3)') opt_cc(2)
         write (*, '(3f15.3)') opt_ahmax(2), opt_ahmin(2), opt_avert(2)
         write (*, '(1f15.3)') opt_cc(3)
         write (*, '(3f15.3)') opt_ahmax(3), opt_ahmin(3), opt_avert(3)
         write (*, *) ' '
      end if
      !*******************************************************
      ! Evaluate if the variograms have a different relation of anisotropy
      !*******************************************************
      ! Is the user asking for attempting the variogram calculations?
      if (varimp .eq. 1) then
         optmajmin1 = opt_ahmax(3)/opt_ahmin(3)
         optmajter1 = opt_ahmax(3)/opt_avert(3)
         if (optmajmin1 .lt. 1.1) then
            ! Experimental
            ! Save first exp var on temp arrays
            ! get the correct size from the arrays
            size_maj_arr = size(lvmajcov)
            size_min_arr = size(lvmincov)
            size_ter_arr = size(lvtercov)
            !major direction
            allocate (tlvmajnum(size_maj_arr))
            allocate (tlvmajcov, tlvmajmu, tlvmajmuh, tlvmajsu, tlvmajsuh, &
                      tlvmajdis, tlvmajwt, mold=tlvmajnum)
            !minor direction
            allocate (tlvminnum(size_min_arr))
            allocate (tlvmincov, tlvminmu, tlvminmuh, tlvminsu, tlvminsuh, &
                      tlvmindis, tlvminwt, mold=tlvminnum)
            ! tertiary direction
            allocate (tlvternum(size_ter_arr))
            allocate (tlvtercov, tlvtermu, tlvtermuh, tlvtersu, tlvtersuh, &
                      tlvterdis, tlvterwt, mold=tlvternum)
            ! Deallocate all arrays
            deallocate (varlagdist, varvalue, varnumpairs, varwts, varazm, vardip)
            write (*, *) 'Attempting to improve experimental variogram'
            ! Calculate experimental variogram with flexible constrains
            call varexpflx(x, y, z, v, lrotmat, nd, usemaj, usemin, useter, nmaj_lag, nmin_lag, nter_lag, &
                           dmajlag, dminlag, dterlag, zavr_anis, usexdir, useydir, usezdir, debug, dep_type, wtscurrent, &
                           tlvmajnum, tlvmajcov, tlvmajmu, tlvmajmuh, tlvmajsu, tlvmajsuh, tlvmajdis, tlvmajwt, &
                           tlvminnum, tlvmincov, tlvminmu, tlvminmuh, tlvminsu, tlvminsuh, tlvmindis, tlvminwt, &
                           tlvternum, tlvtercov, tlvtermu, tlvtermuh, tlvtersu, tlvtersuh, tlvterdis, tlvterwt, &
                           anchx(anx), anchy(anx), anchz(anx), anx, geolocvar, strattrans)
            ! Calculate the new model, update the sizes of the arrays base
            ! on results from varexpflx
            nvarmaj = 0; nvarmin = 0; nvarter = 0
            nvarmaj = size(tlvmajcov)
            nvarmin = size(tlvmincov)
            if (useter) nvarter = size(tlvtercov)
            nvargs = nvarmaj + nvarmin + nvarter
            ! allocate the needed arrays,
            allocate (varlagdist(nvargs))
            allocate (varvalue(nvargs))
            allocate (varnumpairs(nvargs))
            allocate (varwts(nvargs))
            allocate (varazm(nvargs))
            allocate (vardip(nvargs))

            varlagdist(1:nvarmaj) = tlvmajdis
            varlagdist(nvarmaj + 1:nvarmaj + nvarmin) = tlvmindis
            if (useter) varlagdist(nvarmaj + nvarmin + 1:nvargs) = tlvterdis

            varvalue(1:nvarmaj) = tlvmajcov
            varvalue(nvarmaj + 1:nvarmaj + nvarmin) = tlvmincov
            if (useter) varvalue(nvarmaj + nvarmin + 1:nvargs) = tlvtercov

            varnumpairs(1:nvarmaj) = int(tlvmajnum)
            varnumpairs(nvarmaj + 1:nvarmaj + nvarmin) = int(tlvminnum)
            if (useter) varnumpairs(nvarmaj + nvarmin + 1:nvargs) = int(tlvternum)

            varwts(1:nvarmaj) = tlvmajwt
            varwts(nvarmaj + 1:nvarmaj + nvarmin) = tlvminwt
            if (useter) varwts(nvarmaj + nvarmin + 1:nvargs) = tlvterwt

            ! The next parameters are just hard coded to make the code follow the rotmat
            ! convention

            varazm(1:nvarmaj) = 0
            varazm(nvarmaj + 1:nvarmaj + nvarmin) = 90
            if (useter) varazm(nvarmaj + nvarmin + 1:nvargs) = 0
            vardip(1:nvarmaj) = 0
            vardip(nvarmaj + 1:nvarmaj + nvarmin) = 0
            if (useter) vardip(nvarmaj + nvarmin + 1:nvargs) = -90
            ! Optimization function for variogram modelling
            call varmodelfit(nst, c0_, it_, cc_, azm, dip, tilt, ahmax, ahmin, avert, sillf, &
                             maxiter, rseed, minpairs, npairswt, invdistwt, fixhmaxvert, hmaxvert, &
                             fixhminhmax, hminhmax, nvargs, varlagdist, varvalue, varnumpairs, varwts, &
                             varazm, vardip, opt_c0, opt_it, opt_cc, opt_azm, opt_dip, opt_tilt, &
                             opt_ahmax, opt_ahmin, opt_avert)

            if (debug) then
               write (*, *) 'Optimized variogram model'
               write (*, '(1f15.3)') opt_cc(1)
               write (*, '(3f15.3)') opt_ahmax(1), opt_ahmin(1), opt_avert(1)
               write (*, '(1f15.3)') opt_cc(2)
               write (*, '(3f15.3)') opt_ahmax(2), opt_ahmin(2), opt_avert(2)
               write (*, '(1f15.3)') opt_cc(3)
               write (*, '(3f15.3)') opt_ahmax(3), opt_ahmin(3), opt_avert(3)
            end if
            !
            ! Compare and pick the best variogram.
            !
            optmajmin2 = opt_ahmax(3)/opt_ahmin(3)
            optmajter2 = opt_ahmax(3)/opt_avert(3)
            if (anis .eq. 1) then ! did the user provide anis ratio?
               ! which anisratio is closer
               ! to the ones provided by the user?
               difratio1 = abs(rtmajmin - optmajmin1)
               difratio2 = abs(rtmajmin - optmajmin2)
               difratio3 = abs(rtmajter - optmajter1)
               difratio4 = abs(rtmajter - optmajter2)
               if (difratio1 .lt. difratio2 .and. &
                   difratio3 .lt. difratio4) then
                  write (*, *) 'First model is better'
                  deallocate (varlagdist, varvalue, varnumpairs, varwts, varazm, vardip)
                  nvarmaj = 0; nvarmin = 0; nvarter = 0
                  nvarmaj = nmaj_lag
                  nvarmin = nmin_lag
                  if (useter) nvarter = nter_lag
                  nvargs = nvarmaj + nvarmin + nvarter
                  ! allocate the needed arrays,
                  allocate (varlagdist(nvargs))
                  allocate (varvalue(nvargs))
                  allocate (varnumpairs(nvargs))
                  allocate (varwts(nvargs))
                  allocate (varazm(nvargs))
                  allocate (vardip(nvargs))
                  ! lag distances
                  varlagdist(1:nvarmaj) = lvmajdis
                  varlagdist(nvarmaj + 1:nvarmaj + nvarmin) = lvmindis
                  if (useter) varlagdist(nvarmaj + nvarmin + 1:nvargs) = lvterdis
                  ! variogram values
                  varvalue(1:nvarmaj) = lvmajcov
                  varvalue(nvarmaj + 1:nvarmaj + nvarmin) = lvmincov
                  if (useter) varvalue(nvarmaj + nvarmin + 1:nvargs) = lvtercov
                  ! number of pairs
                  varnumpairs(1:nvarmaj) = int(lvmajnum)
                  varnumpairs(nvarmaj + 1:nvarmaj + nvarmin) = int(lvminnum)
                  if (useter) varnumpairs(nvarmaj + nvarmin + 1:nvargs) = int(lvternum)
                  ! weights
                  varwts(1:nvarmaj) = lvmajwt
                  varwts(nvarmaj + 1:nvarmaj + nvarmin) = lvminwt
                  if (useter) varwts(nvarmaj + nvarmin + 1:nvargs) = lvterwt
               else ! the second model is better
                  write (*, *) 'Second model is better'
                  deallocate (varlagdist, varvalue, varnumpairs, varwts, varazm, vardip)
                  nvarmaj = 0; nvarmin = 0; nvarter = 0
                  nvarmaj = size(tlvmajcov)
                  nvarmin = size(tlvmincov)
                  if (useter) nvarter = size(tlvtercov)
                  nvargs = nvarmaj + nvarmin + nvarter
                  ! allocate the needed arrays,
                  allocate (varlagdist(nvargs))
                  allocate (varvalue(nvargs))
                  allocate (varnumpairs(nvargs))
                  allocate (varwts(nvargs))
                  allocate (varazm(nvargs))
                  allocate (vardip(nvargs))
                  ! lag distances
                  varlagdist(1:nvarmaj) = tlvmajdis
                  varlagdist(nvarmaj + 1:nvarmaj + nvarmin) = tlvmindis
                  if (useter) varlagdist(nvarmaj + nvarmin + 1:nvargs) = tlvterdis
                  ! variogram values
                  varvalue(1:nvarmaj) = tlvmajcov
                  varvalue(nvarmaj + 1:nvarmaj + nvarmin) = tlvmincov
                  if (useter) varvalue(nvarmaj + nvarmin + 1:nvargs) = tlvtercov
                  ! number of pairs
                  varnumpairs(1:nvarmaj) = int(tlvmajnum)
                  varnumpairs(nvarmaj + 1:nvarmaj + nvarmin) = int(tlvminnum)
                  if (useter) varnumpairs(nvarmaj + nvarmin + 1:nvargs) = int(tlvternum)
                  ! weights
                  varwts(1:nvarmaj) = tlvmajwt
                  varwts(nvarmaj + 1:nvarmaj + nvarmin) = tlvminwt
                  if (useter) varwts(nvarmaj + nvarmin + 1:nvargs) = tlvterwt
               end if
            else ! If no information of ratio, compare by default
               if (optmajmin1 .gt. (optmajmin2 + SMALLDBLE)) then!
                  write (*, *) 'First model is better'
                  deallocate (varlagdist, varvalue, varnumpairs, varwts, varazm, vardip)
                  nvarmaj = 0; nvarmin = 0; nvarter = 0
                  nvarmaj = nmaj_lag
                  nvarmin = nmin_lag
                  if (useter) nvarter = nter_lag
                  nvargs = nvarmaj + nvarmin + nvarter
                  ! allocate the needed arrays,
                  allocate (varlagdist(nvargs))
                  allocate (varvalue(nvargs))
                  allocate (varnumpairs(nvargs))
                  allocate (varwts(nvargs))
                  allocate (varazm(nvargs))
                  allocate (vardip(nvargs))
                  ! lag distances
                  varlagdist(1:nvarmaj) = lvmajdis
                  varlagdist(nvarmaj + 1:nvarmaj + nvarmin) = lvmindis
                  if (useter) varlagdist(nvarmaj + nvarmin + 1:nvargs) = lvterdis
                  ! variogram values
                  varvalue(1:nvarmaj) = lvmajcov
                  varvalue(nvarmaj + 1:nvarmaj + nvarmin) = lvmincov
                  if (useter) varvalue(nvarmaj + nvarmin + 1:nvargs) = lvtercov
                  ! number of pairs
                  varnumpairs(1:nvarmaj) = int(lvmajnum)
                  varnumpairs(nvarmaj + 1:nvarmaj + nvarmin) = int(lvminnum)
                  if (useter) varnumpairs(nvarmaj + nvarmin + 1:nvargs) = int(lvternum)
                  ! weights
                  varwts(1:nvarmaj) = lvmajwt
                  varwts(nvarmaj + 1:nvarmaj + nvarmin) = lvminwt
                  if (useter) varwts(nvarmaj + nvarmin + 1:nvargs) = lvterwt
               else
                  write (*, *) 'Second model is better'
                  deallocate (varlagdist, varvalue, varnumpairs, varwts, varazm, vardip)
                  nvarmaj = 0; nvarmin = 0; nvarter = 0
                  nvarmaj = size(tlvmajcov)
                  nvarmin = size(tlvmincov)
                  if (useter) nvarter = size(tlvtercov)
                  nvargs = nvarmaj + nvarmin + nvarter
                  ! allocate the needed arrays
                  allocate (varlagdist(nvargs))
                  allocate (varvalue(nvargs))
                  allocate (varnumpairs(nvargs))
                  allocate (varwts(nvargs))
                  allocate (varazm(nvargs))
                  allocate (vardip(nvargs))
                  ! lag distances
                  varlagdist(1:nvarmaj) = tlvmajdis
                  varlagdist(nvarmaj + 1:nvarmaj + nvarmin) = tlvmindis
                  if (useter) varlagdist(nvarmaj + nvarmin + 1:nvargs) = tlvterdis
                  ! variogram values
                  varvalue(1:nvarmaj) = tlvmajcov
                  varvalue(nvarmaj + 1:nvarmaj + nvarmin) = tlvmincov
                  if (useter) varvalue(nvarmaj + nvarmin + 1:nvargs) = tlvtercov
                  ! number of pairs
                  varnumpairs(1:nvarmaj) = int(tlvmajnum)
                  varnumpairs(nvarmaj + 1:nvarmaj + nvarmin) = int(tlvminnum)
                  if (useter) varnumpairs(nvarmaj + nvarmin + 1:nvargs) = int(tlvternum)
                  ! weights
                  varwts(1:nvarmaj) = tlvmajwt
                  varwts(nvarmaj + 1:nvarmaj + nvarmin) = tlvminwt
                  if (useter) varwts(nvarmaj + nvarmin + 1:nvargs) = tlvterwt
               end if
            end if
            ! Deallocate all temporary arrays
            deallocate (tlvmajnum, tlvmajcov, tlvmajmu, tlvmajmuh, tlvmajsu, &
                        tlvmajsuh, tlvmajdis, tlvmajwt)
            deallocate (tlvminnum, tlvmincov, tlvminmu, tlvminmuh, tlvminsu, &
                        tlvminsuh, tlvmindis, tlvminwt)
            deallocate (tlvternum, tlvtercov, tlvtermu, tlvtermuh, tlvtersu, &
                        tlvtersuh, tlvterdis, tlvterwt)
         end if
      end if
      ! ----------------------------------------------------------------------
      ! Prepare output files
      ! ---------------------------------------------------------------------
      ! For geolocated variograms
      if (geolocvar .eq. 1) then
         dbgfl = dbgnames(anx)
         outfl_exp = expnames(anx)
         outfl_mod = modelnames(anx)
         varmodgs = modgsnames(anx)
         wtsfl = wtsflnames(anx)
      end if
      ! For the general variogram (file names declared by user)
      if (geolocvar .eq. 0 .and. geoloc) then
         dbgfl = trim(adjustl(dbgfl))//".out"
         outfl_exp = trim(adjustl(outfl_exp))//".out"
         outfl_mod = trim(adjustl(outfl_mod))//".out"
         varmodgs = trim(adjustl(varmodgs))//".out"
         wtsfl = trim(adjustl(wtsfl))//".out"
      end if
      !-----------------------------------------------------------
        !! Write the results of the eigenvectors as debug file
      !-----------------------------------------------------------
      i_exp = 9
      open (i_exp, file=dbgfl, status='UNKNOWN')
      write (i_exp, 101)
101   format('Eigenvectors', /, '1', /, 'x y z')
      do i = 1, 3
         do k = 1, 3
            write (i_exp, '(f14.4)') right_eiveco(i, k)
         end do
      end do
      close (i_exp)
      !-----------------------------------------------------------
        !! Write the results of the variogram (not correlogram)
      !-----------------------------------------------------------
      ! Assigning the referential angles
      allocate (varazm_ref(nvargs))
      allocate (vardip_ref(nvargs))
      varazm_ref(1:nvarmaj) = az_major
      varazm_ref(nvarmaj + 1:nvarmaj + nvarmin) = az_minor
      if (ndir .eq. 3) varazm_ref(nvarmaj + nvarmin + 1:nvargs) = az_ter

      vardip_ref(1:nvarmaj) = dip_major
      vardip_ref(nvarmaj + 1:nvarmaj + nvarmin) = dip_minor
      if (ndir .eq. 3) vardip_ref(nvarmaj + nvarmin + 1:nvargs) = dip_ter

      ! Hard coding indexes:
      allocate (varindex(nvargs))

      varindex(1:nvarmaj) = 1
      varindex(nvarmaj + 1:nvarmaj + nvarmin) = 2
      if (ndir .eq. 3) varindex(nvarmaj + nvarmin + 1:nvargs) = 3

      i_exp = 10
      open (i_exp, file=outfl_exp, status='UNKNOWN')
      write (i_exp, 102)
102   format('Experimental variogram', /, '6', /, 'Variogram Index', /, &
             'Lag Distance', /, 'Number of Pairs', /, 'Variogram Value', /, &
             'Inferred Azimuth', /, 'Inferred Dip')

      do i = 1, nvargs
         write (i_exp, '(i1,f16.4,i10,3f14.4)') varindex(i), varlagdist(i), varnumpairs(i), &
            varvalue(i), varazm_ref(i), vardip_ref(i)
      end do
      close (i_exp)

      !***********************************************************
        !! Write the results of the variogram model
      !***********************************************************
      ! Hard code some parameters to make the writing file similar to
      ! varmodel output
      allocate (npoints(ndir))
      allocate (pointsep(ndir))
      if (ndir .eq. 2) then
         npoints = 1000
         pointsep(1) = 1.1*maxval(opt_ahmax)/1000
         pointsep(2) = 1.1*maxval(opt_ahmin)/1000
      end if

      if (ndir .eq. 3) then
         npoints = 1000
         pointsep(1) = 1.1*maxval(opt_ahmax)/1000
         pointsep(2) = 1.1*maxval(opt_ahmin)/1000
         pointsep(3) = 1.1*maxval(opt_avert)/1000
      end if

      allocate (pointazm(ndir))
      allocate (pointdip(ndir))

      pointazm(1) = az_major
      pointazm(2) = az_minor
      if (useter) pointazm(3) = az_ter
      pointdip(1) = dip_major
      pointdip(2) = dip_minor
      if (useter) pointdip(3) = dip_ter

      i_exp = 11
      open (i_exp, file=outfl_mod, status='UNKNOWN')
      write (i_exp, 103)
103   format('Variogram model points', /, '7', /, 'Variogram Index', /, &
             'Lag Distance', /, 'Number of Pairs', /, 'Variogram Value', /, &
             'Variogram Number', /, 'Inferred Azimuth', /, 'Inferred Dip')

      call varmodelsph(nst, opt_c0, opt_it, opt_cc, opt_ahmax, &
                       opt_ahmin, opt_avert, pointsep, varmodelsph_, lagdist_)
      ! hard wired variogram indexes, azm and dip
      varind(1:1000) = 1
      varind(1001:2000) = 2
      if (ndir .eq. 3) varind(2001:3000) = 3

      model_varazm(1:1000) = az_major
      model_varazm(1001:2000) = az_minor
      if (ndir .eq. 3) model_varazm(2001:3000) = az_ter

      model_vardip(1:1000) = dip_major
      model_vardip(1001:2000) = dip_minor
      if (ndir .eq. 3) model_vardip(2001:3000) = dip_major
      if (ndir .eq. 3) then
         do i = 1, 3000
            write (i_exp, '(i6,f18.6,i10,f18.6,i6,2f18.6)') varind(i), lagdist_(i), 1, &
               varmodelsph_(i), 1, model_varazm(i), model_vardip(i)
         end do
      end if
      if (ndir .eq. 2) then
         do i = 1, 2000
            write (i_exp, '(i6,f18.6,i10,f18.6,i6,2f18.6)') varind(i), lagdist_(i), 1, &
               varmodelsph_(i), 1, model_varazm(i), model_vardip(i)
         end do
      end if

      close (i_exp)

      !***********************************************************
        !! Write the results of the variogram model as text file
      !***********************************************************

      i_exp = 12
      open (i_exp, file=varmodgs, status='UNKNOWN')
      ! For 3D variograms, user needs to write the angles in the right dimension
      if (ndir .eq. 3) then
         write (i_exp, '(i3,f10.2,a)') 3, 0, '    -nst, nugget effect'
         do ist = 1, nst
            write (i_exp, '(i3,4f10.2,a)') 1, opt_cc(ist), azmf, dipf, tilf, &
               '    -it,cc,azm,dip,tilt'
            write (i_exp, '(3f14.2,a)') opt_ahmax(ist), opt_ahmin(ist), opt_avert(ist), &
               '    -a_hmax, a_hmin, a_vert (ranges)'
         end do
      end if

      if (ndir .eq. 2) then
         write (i_exp, '(i3,f10.2,a)') 3, 0, '    -nst, nugget effect'
         do ist = 1, nst
            write (i_exp, '(i3,4f10.2,a)') 1, opt_cc(ist), a1, a2, a3, &
               '    -it,cc,azm,dip,tilt'
            write (i_exp, '(3f14.2,a)') opt_ahmax(ist), opt_ahmin(ist), opt_avert(ist), &
               '    -a_hmax, a_hmin, a_vert (ranges)'
         end do
      end if

      close (i_exp)
      !***********************************************************
        !! Write the results of the weights
      !***********************************************************
      i_exp = 13
      open (i_exp, file=wtsfl, status='UNKNOWN')
      fmt = '(I2.2)'
      write (ii, fmt) anx !!ii is the character version of i
      write (i_exp, 105)
105   format('Weights', /, '4', /, 'X', /, 'Y', /, 'Z')
      write (i_exp, '(a,i3)') 'Anchor ', anx
      do i = 1, nd
         write (i_exp, '(4f14.4)') xyz(i, 1), xyz(i, 2), xyz(i, 3), wtsanch(i, anx)
      end do
      close (i_exp)

      !Deallocate all arrays to move to the following anchor if needed
      deallocate (sq_df1, sqdfdist1, dx1, dy1, dz1, wtsomni1)
      deallocate (dmajco, dminco, dterco)
      deallocate (xyz_anis)
      deallocate (lvmajnum, lvmajcov, lvmajmu, lvmajmuh, lvmajsu, lvmajsuh, &
                  lvmajdis, lvmajwt)
      deallocate (lvminnum, lvmincov, lvminmu, lvminmuh, lvminsu, lvminsuh, &
                  lvmindis, lvminwt)
      deallocate (lvternum, lvtercov, lvtermu, lvtermuh, lvtersu, lvtersuh, &
                  lvterdis, lvterwt)
      deallocate (varlagdist, varvalue, varnumpairs, varwts, varazm, vardip)
      deallocate (varazm_ref, vardip_ref, varindex)
      deallocate (npoints, pointsep, pointazm, pointdip)

      deallocate (xx, yy, zz, vvax, wtsax)

      if (ndir .eq. 3 .and. dep_type .eq. 1) deallocate (x_tmp, y_tmp, z_tmp, v_tmp, w_tmp)

   end do
40 continue

   ! After finishing general variogram, calculate the variograms per anchor
   if (geoloc) then
      geolocvar = 1
      nanchor = nanchororig
      write (*, *) ' '
      write (*, *) ' General variogram calculated'
      write (*, *) '  Starting calculations for each anchor'
      write (*, *) '  --------------------------------------  '
      geoloc = .FALSE.
      go to 111
   end if

   write (*, 994) VERSION
994 format(/' Autovar Finished, Version: ', f7.3/)
98 continue
99 continue
end program Autovar

subroutine makepar
   !-----------------------------------------------------------------------
   !
   !                      Write a Parameter File
   !                      **********************
   !
   !
   !
   !-----------------------------------------------------------------------
   lun = 99
   open (lun, file='Autovar.par', status='UNKNOWN')
   write (lun, 10)
10 format('                Parameters for Autovar', /, &
          '                ************************', /, /, &
          'START OF PARAMETERS:')

   write (lun, 11)
11 format('cluster.dat                      ', &
          '-file with data')
   write (lun, 12)
12 format('1  2  0  3   5                   ', &
          '-   columns for X,Y,Z,var,id')
   write (lun, 13)
13 format('-1.0e21   1.0e21                 ', &
          '-   trimming limits')

   write (lun, 15)
15 format('1   2                             ', &
          '-Geolocated variogram? (1=yes), Number of anchors')

   write (lun, 16)
16 format('anchor.out                                ', &
          '-file with anchors')

   write (lun, 17)
17 format('2 3 4                         ', &
          '-X Y Z column of anchor file')

   write (lun, 18)
18 format('1                         ', &
          '-Attemp to improve variogram?')

   write (lun, 19)
19 format('1  2  4                             ', &
          '- Use general anisotropy ratio,maj/min,maj/ter')

   write (lun, 21)
21 format('50   0.5    1.0                  ', &
          '-nx,xmn,xsiz')
   write (lun, 22)
22 format('50   0.5    1.0                  ', &
          '-ny,ymn,ysiz')
   write (lun, 23)
23 format('1    0.5    1.0                  ', &
          '-nz,zmn,zsiz')
   write (lun, 20)
20 format('autovarmod.out                   ', &
          '-file for variogram model output')
   write (lun, 24)
24 format('autovar.out                      ', &
          '-   file for experimental pairs output')
   write (lun, 25)
25 format('debug.out                        ', &
          '-   file for directions')
   write (lun, 26)
26 format('varfit.out                        ', &
          '-   variogram model text GSLIB format')
   write (lun, 27)
27 format('weights.out                        ', &
          '-   weights for each anchor')
   write (lun, 28)
28 format('1                                ', &
          '-   Debug on screen (1=yes)')

   write (lun, 29)
29 format('1                        ', &
          '- Disseminated(1) or stratigraphic(2) deposit type? ')

   close (lun)
   return
end

subroutine chknam(str, len)
   !-----------------------------------------------------------------------
   !
   !                   Check for a Valid File Name
   !                   ***************************
   !
   ! This subroutine takes the character string "str" of length "len" and
   ! removes all leading blanks and blanks out all characters after the
   ! first blank found in the string (leading blanks are removed first).
   !
   !
   !
   !-----------------------------------------------------------------------
   character(len=*), intent(inout) :: str
   integer itrim
   !
   ! Remove leading blanks:
   str = adjustl(str)
   !
   ! find first two blanks and blank out remaining characters:
   itrim = index(str, '   ')
   if (itrim > 0) str(itrim:) = ' '
   !
   ! Look for "-fi"
   itrim = index(str, '-fi')
   if (itrim > 0) str(itrim:) = ' '
   !
   ! Look for "\fi"
   itrim = index(str, '\fi')
   if (itrim > 0) str(itrim:) = ' '
   !
   ! Return with modified file name:
   return
end

