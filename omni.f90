    module alldir

    public ::  omnipair
    contains

    ! Subroutine to calculate the omnidirectional pairs inside the cube
    ! Modified from varcalc.
    subroutine omnipair(ndv,maxdistv,tminv,tmaxv,xv,yv,zv,vv1,&
        geoloc,wts,outputidxv,sqdfdistv,dxv1,dyv1,dzv1,&
        sqdfv,idfv,wtsomniv,debugv)
    ! Variable declaration
    use sortem
    implicit none
    ! input variables
    integer, intent (in) :: ndv,geoloc
    real(kind=8), dimension(ndv),intent(in) :: xv,yv,zv,vv1,wts
    real(kind=8), intent(in) :: maxdistv,tminv,tmaxv
    logical, intent(in) :: debugv
    ! output variables
    real(kind=8), allocatable,intent(out) ::dxv1(:),dyv1(:),dzv1(:),sqdfv(:),sqdfdistv(:),&
        wtsomniv(:)
    integer, allocatable,intent(out) :: idfv(:)
    integer,intent(out) :: outputidxv
    ! internal variables
    real(kind=8) :: h,srch_rnge,xdist,ydist,zdist,ymin,zmin,ymax,zmax,&
        farthest,azmtolrad,diptolrad,dist,maxx_dist,av,ss
    integer, allocatable :: introtindex(:)
    real(kind=8), allocatable :: dblerotindex(:),xyz(:,:)
    integer :: i,j,ntotal,nps,startidx
    logical :: omni
    real(kind=8), parameter :: PI=4*atan(1.0d0)
    real(kind=8), parameter :: DEG2RAD=PI/180d0, RAD2DEG=180d0/PI
    ! End of variable declaration

    srch_rnge = 0.2*maxdistv
    nps = 0 !number of pairs selected
    ! Get the average and variance of the selected points
    av=0;ss=0
    do i = 1, ndv
        av = av + vv1(i)
        ss = ss + vv1(i)*vv1(i)
    end do
    av = av / max(real(ndv),1.0)
    ss =(ss / max(real(ndv),1.0)) - av * av
    ! allocate in convenient way and do calculations
    allocate(xyz(ndv,3))
    xyz(:,1) = xv(:)
    xyz(:,2) = yv(:)
    xyz(:,3) = zv(:)
    do i=1,ndv-1   ! This point, i, is the "tail" variable in GSLIB notation

        do j=1,ndv   ! This point, j, is the "head" variable in GSLIB notation
            ! Check if the points are at the same point on the xv' axis
            xdist = (xyz(j,1)-xyz(i,1))
            ydist = (xyz(j,2)-xyz(i,2))
            zdist = (xyz(j,3)-xyz(i,3))
            ! for finding the right distances.
            dist = sqrt(xdist**2+ydist**2+zdist**2)
            if (dist .le. srch_rnge) then
                nps = nps+1
            else
                continue
            end if
        end do
    end do

    ! Note: the calculation of the pairs come from a modification of varcalc.
    ! Allocate the needed memory
    ntotal = nps

    allocate(sqdfv(ntotal), dxv1(ntotal),dyv1(ntotal), dzv1(ntotal),idfv(ntotal),&
        sqdfdistv(ntotal),wtsomniv(ntotal))

    if (debugv) write(*,*) ' working on selection of pairs '



    ! Initialize index tracking for values
    allocate(dblerotindex(ndv))
    do i=1,ndv
        dblerotindex(i) = i
    end do

    ! Set omnidirectional search
    azmtolrad = 90d0*DEG2RAD
    diptolrad = 90d0*DEG2RAD
    omni = .true.
    ! Sort along X
    call dblemodsortem(xyz(:,1),ndv,3,xyz(:,2),xyz(:,3),dblerotindex)
    ! Track the sorted indices for computing values
    allocate(introtindex(ndv))
    do i=1,ndv
        introtindex(i) = int(dblerotindex(i))
    end do

    !calculate the maximum distance between pairs to
    !Restrict to selected part of the anchor:
    if (geoloc .eq. 1) then
        maxx_dist = 1*sqrt((maxval(xyz(:,1))-minval(xyz(:,1)))**2 + &
            (maxval(xyz(:,2))-minval(xyz(:,2)))**2 + &
            (maxval(xyz(:,3))-minval(xyz(:,3)))**2)
    end if
    if (geoloc .eq. 0) then
        !calculate the maximum distance between pairs to only calculate points inside this range
        !Restrict to 0.2 domain size:
        maxx_dist = 0.2*sqrt((maxval(xyz(:,1))-minval(xyz(:,1)))**2 + &
            (maxval(xyz(:,2))-minval(xyz(:,2)))**2 + &
            (maxval(xyz(:,3))-minval(xyz(:,3)))**2)
    end if

    ! set the search radius as two times the value entered by the user
    ! in the anisotropy ranges.

    !rma2 = 2*rma
    !rmi2 = 2*rmi
    !rte2 = 2*rte

    ! Main variogram loop over all points indexed by i
    outputidxv = 0

    do i=1,ndv-1   ! This point, i, is the "tail" variable in GSLIB notation
        ! Farthest point to consider from the current location
        if (ndv .lt. 10000) then
            farthest = xyz(i,1) + maxx_dist
        else
            farthest = xyz(i,1) + maxx_dist/2
        endif
        ! Pre-calculated bandwidth min/max yv's and zv's
        ymax = xyz(i,2) + maxx_dist
        ymin = xyz(i,2) - maxx_dist
        zmax = xyz(i,3) + maxx_dist
        zmin = xyz(i,3) - maxx_dist
        ! Loop over all possible pairs forward of the current point
        startidx = 1

        do j=startidx,ndv-1   ! This point, j, is the "head" variable in GSLIB notation
            ! If we are too far away, skip checking other pairs
            if (xyz(j,1) .gt. farthest) exit
            ! Check horizontal bandwidths
            if (xyz(j,2) .gt. ymax) cycle
            if (xyz(j,2) .lt. ymin) cycle
            ! Check vertical bandwidths
            if (xyz(j,3) .gt. zmax) cycle
            if (xyz(j,3) .lt. zmin) cycle
            ! Check if the points are at the same point on the xv' axis
            xdist = (xyz(j,1)-xyz(i,1))! I remove the absval to consider direction
            ydist = (xyz(j,2)-xyz(i,2))
            zdist = (xyz(j,3)-xyz(i,3))
            ! for finding the right distances.
            if (i.eq.j) then
                cycle
            end if

            ! Find the distance between the values
            h = sqrt((xyz(i,1)-xyz(j,1))**2 + &
                (xyz(i,2)-xyz(j,2))**2 + &
                (xyz(i,3)-xyz(j,3))**2)

            if (h .le. srch_rnge) then
                ! We are inside the searchradius, so calculate each sqdfv value

                ! Are both values licit?

                if (vv1(introtindex(i)) .le. tminv .or. &
                    vv1(introtindex(i)) .gt. tmaxv) cycle
                if (vv1(introtindex(j)) .le. tminv .or. &
                    vv1(introtindex(j)) .gt. tmaxv) cycle

                ! Not trimming, calculate and add the values on
                outputidxv = outputidxv + 1
                ! Set the first table output
                idfv(outputidxv) = outputidxv
                sqdfdistv(outputidxv) = h
                dxv1(outputidxv) = xdist
                dyv1(outputidxv) = ydist
                dzv1(outputidxv) = zdist
                !wtsomniv(outputidxv) = abs((wts(i) + wts(j))/2)
                wtsomniv(outputidxv) = abs((wts(i) * wts(j)))
                ! Pseudo Correlogram
                sqdfv(outputidxv)=(vv1(introtindex(i))*vv1(introtindex(j))-(av)**2)/ss

            else
                continue
            end if

        end do
    end do
    return

    end subroutine omnipair


    end module alldir