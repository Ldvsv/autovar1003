    module stratco
    public :: stratcot
    contains
    ! Subroutine to calculate the stratigraphic coordinates using the proportional
    ! method. Modified from StratCoord.exe. Parameters are inferred from dataset.
    ! Assumes drillholes are ordered and grouped.

    subroutine stratcot(xyzstr,vstr,ndstr,idcolstr,idstr,zavrstr,debugstr,lavrstr,xsstr,ysstr,zsstr,vsstr)
    ! Variable declaration
    use sortem
    implicit none

    ! Input variables from main
    integer, intent(in) :: ndstr,idcolstr!number of data and idstr column number
    real(kind = 8), dimension(ndstr,3), intent(in):: xyzstr
    real(kind = 8), dimension(ndstr), intent(in):: vstr ! value of variables
    real(kind = 8), dimension(ndstr), intent(in) :: idstr !column with the dhid
    real(kind = 8) :: zavrstr,lavrstr ! Average spacing and average dh lenght
    logical, intent(in):: debugstr
    !Internal variables
    integer, parameter :: istyle=1
    real(kind = 8), dimension(ndstr,3) :: xyzidstr
    real(kind = 8), dimension(ndstr) :: vv
    real(kind = 8), dimension(ndstr) :: iid
    real, parameter :: EPSLON = 1.0e-2, MAXCOL=50
    real(kind=8) :: rnd,zdist,ztemp,lenght_temp,elevmin,elevmax
    real(kind=8), allocatable ::dh(:)
    real(kind=8), allocatable :: zdh(:),zdistdh(:),lenght(:),scoord(:)
    integer :: idh,idd,i,j,k,l,ndh,n,j_acc
    !outputs
    !real(kind= 8), intent(out):: xyzstrat(ndstr,4)
    real(kind=8), allocatable,intent(out) ::xsstr(:),ysstr(:),zsstr(:),vsstr(:)
    !
    ! START
    !
    !drillhole
    if (idcolstr .ne. 0) then
        idh = idcolstr
        idd = 1 !if 1 to use DHID
    else
        idd = 0
    end if
    ! Assign to local variable
    if(idd==1) then
        allocate(dh(ndstr))
        dh = idstr
        ! sort first and count the number of changes to find the number of dh
        call dblemodsortem(dh,ndstr,0)
        ndh=1

        do i=2,ndstr
            if(dh(i)/=dh(i-1)) ndh = ndh + 1
        end do
    end if

    ! Read in the values
    xyzidstr(:,1) = xyzstr(:,1)
    xyzidstr(:,2) = xyzstr(:,2)
    xyzidstr(:,3) = xyzstr(:,3)
    iid = real(idstr)
    vv = vstr
    ! Allocate array for output
    allocate(xsstr(ndstr),ysstr(ndstr),zsstr(ndstr),vsstr(ndstr))
    !
    ! Find the Zstrat coordinates
    ! Grouping by DH

    if(idd==1) then
        idh = 1
        j_acc = 0
        ndh = 1
        do i=2,ndstr
            if (i.lt.j_acc+1) cycle
            j = 0
            if(dh(i)/=dh(i-1)) ndh = ndh + 1
            !
            ! Analyze one dh at the time
            !
            do n = 1,ndstr
                if (dh(n).gt.dh(i)) exit
                if (dh(i) .eq. dh(n)) then
                    j = j+1
                end if
            end do
            allocate(zdh(j))
            allocate(scoord(j))
            j_acc = j_acc +j
            zdh = xyzidstr(j_acc-j+1:j_acc,3)

            ! Find the min/max elevation for this drillhole
            elevmin = minval(zdh)
            elevmax = maxval(zdh)
            elevmin = elevmin - zavrstr/2.0
            elevmax = elevmax + zavrstr/2.0
            ! Calculate the zstrat coord for this drillhole
            do l=1,j
                if(istyle.eq.1) then
                    scoord(l) = (zdh(l)-elevmin) /&
                        (elevmax-elevmin) * lavrstr
                end if

                if(istyle.eq.2) then
                    elevmin = elevmax - zavrstr
                    scoord(l) = (zdh(l)-elevmin) /&
                        (elevmax-elevmin) * lavrstr
                end if

                if(istyle.eq.3) then
                    elevmax = elevmin + zavrstr
                    scoord(l) = (zdh(l)-elevmin) /&
                        (elevmax-elevmin) * lavrstr
                end if
            end do

            ! storage in the array
            zsstr(j_acc-j+1:j_acc) = scoord
            ! save X, and Y to avoid mixing data
            xsstr(j_acc-j+1:j_acc) = xyzidstr(j_acc-j+1:j_acc,1)
            ysstr(j_acc-j+1:j_acc) = xyzidstr(j_acc-j+1:j_acc,2)
            vsstr(j_acc-j+1:j_acc) = vv(j_acc-j+1:j_acc)
            ! deallocate to analyze next drillhole
            deallocate(zdh)
            deallocate(scoord)
        end do
    end if
    !if (debugstr) then
    !    write(*,*) ' Stratigraphic transformation using proportional method '
    !    write(*,*) ''
    !    write(*,*) ''
    !end if
    return
    end subroutine stratcot

    end module stratco