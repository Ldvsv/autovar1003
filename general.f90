    module gen_subs
    public setrotmat, getoctvector,get_octant
    contains
   

    subroutine setrotmat(azm,dip,tilt,forwardrotmat,reverserotmat)
    !-----------------------------------------------------------------------
    !
    !        Sets up the forward rotation matrix and reverse matrix
    !
    ! (c) Jared Deutsch, 2014
    !-----------------------------------------------------------------------
    implicit none
    real(kind=8), parameter :: PI = 4*atan(1.0d0)
    real(kind=8), intent(in) :: azm,dip,tilt
    real(kind=8), dimension(3,3) :: azmrotmat,diprotmat,tiltrotmat
    real(kind=8), dimension(3,3), intent(out) :: forwardrotmat,reverserotmat
    real(kind=8) :: angle

    ! Get the angle to rotate about the Z axis in radians
    angle = (90d0-azm)*PI/180d0
    ! Rotation matrix for azimuth correction
    azmrotmat = 0d0
    azmrotmat(3,3) = 1d0
    azmrotmat(1,1) = cos(angle)
    azmrotmat(1,2) = -1d0*sin(angle)
    azmrotmat(2,1) = sin(angle)
    azmrotmat(2,2) = cos(angle)

    ! Get the angle to rotate about the new X' axis
    angle = -1d0*(dip)*PI/180d0
    ! Rotation matrix for dip correction
    diprotmat = 0d0
    diprotmat(2,2) = 1d0
    diprotmat(1,1) = cos(angle)
    diprotmat(1,3) = sin(angle)
    diprotmat(3,1) = -1d0*sin(angle)
    diprotmat(3,3) = cos(angle)

    ! Get the angle to rotate about the new Y' axis
    angle = -1d0*(tilt)*PI/180d0
    ! Rotation matrix for tilt correction
    tiltrotmat = 0d0
    tiltrotmat(1,1) = 1d0
    tiltrotmat(2,2) = cos(angle)
    tiltrotmat(2,3) = sin(angle)
    tiltrotmat(3,2) = -1d0*sin(angle)
    tiltrotmat(3,3) = cos(angle)

    ! Complete forward rotation matrix is the product of these 2 matrices
    forwardrotmat = matmul(matmul(azmrotmat,diprotmat),tiltrotmat)
    ! Reverse rotation matrix is the transpose of the forward matrix
    ! as these matrices are orthogonal
    reverserotmat = transpose(forwardrotmat)
    end subroutine setrotmat

    subroutine getoctvector(vector, vec_oct)
    ! Subroutine to identify the octant of a vector
    ! The first vector component is the X, the second is the Y and the
    ! third is Z coordinate.
    implicit none
    ! Internal variables:
    integer :: i,j,k
    ! Input variable
    real(kind=8), intent(in) :: vector(3)
    ! Output variable
    integer, intent(out) :: vec_oct
    ! Start
    if (vector(1) .ge. 0.0) then ! 1,2,5,6
        if (vector(2) .ge. 0.0) then ! 1,5
            if (vector(3) .ge. 0.0) then
                vec_oct = 1
            else
                vec_oct = 5
            end if

        else ! 5,6
            if (vector(3) .ge. 0.0) then
                vec_oct = 2
            else
                vec_oct = 6
            end if
        end if
    else ! 3,4,7,8
        if (vector(2) .ge. 0.0) then ! 4,8
            if (vector(3) .ge. 0.0) then
                vec_oct = 4
            else
                vec_oct = 8
            end if

        else ! 3,7
            if (vector(3) .ge. 0.0) then
                vec_oct = 3
            else
                vec_oct = 7
            end if
        end if
    end if

    return

    end subroutine getoctvector

    subroutine get_octant(rotmat, octmaj, octmin, octter,right_rotmat)
    ! Subroutine to ensure that the eigenvectors are consistent with
    ! GSLIB angles.

    ! This is necessary to enforce the correct direction of the eigenvectors
    ! It will work only for 3D and when all the inferred angles are different than
    ! zero! There is no need to test this on 2D.

    ! Internal variables:
    implicit none
    integer :: i,j,k,l,ini_octindex,ref_valid_idx,optoctindex
    integer :: valid_octs(48)
    integer :: opt_octmaj, opt_octmin, opt_octter
    logical :: matcha
    real(kind = 8) :: vecmaj(3),vecmin(3),vecter(3),&
        optvecmaj(3),optvecmin(3),optvecter(3)

    ! Input variable
    real(kind=8), intent(in) :: rotmat(3,3)
    ! Output variable
    integer, intent(out) :: octmaj, octmin, octter
    real(kind = 8), intent(out) :: right_rotmat(3,3)
    !
    ! Start
    !
    ! Declare all the valid combinations for the eigenvectors:
    valid_octs = [132, 142, 143, 154, 183, 184, 213, 214, 243, 251, &
        254 , 261, 314, 321, 324, 361, 362, 372, 421, 431,&
        432 , 472, 473, 483, 512, 541, 542, 574, 581, 584, 612,&
        613 , 623, 651, 652, 681, 723, 724, 734, 752, 762, 763,&
        831 , 834, 841, 863, 873, 874]
    ! Assign the right values to the vectors from the eigenvectors
    vecmaj(1) = rotmat(1,1)
    vecmaj(2) = rotmat(2,1)
    vecmaj(3) = rotmat(3,1)

    vecmin(1) = rotmat(1,2)
    vecmin(2) = rotmat(2,2)
    vecmin(3) = rotmat(3,2)

    vecter(1) = rotmat(1,3)
    vecter(2) = rotmat(2,3)
    vecter(3) = rotmat(3,3)
    ! Get the octant for each vector
    call getoctvector(vecmaj, octmaj)
    call getoctvector(vecmin, octmin)
    call getoctvector(vecter, octter)

    ! Get the index for the octants
    ini_octindex = int(octmaj*100+octmin*10+octter)
    ! Is this a valid set of eigenvectors for the GSLIB angle convention?
    do i =1, 48
        ref_valid_idx = valid_octs(i)
        if (ini_octindex .ne. ref_valid_idx) then
            matcha = .false.
            continue ! test all the values.
        else ! if the index is correct, just exit the subroutine
            matcha = .true.
            go to 101
        end if
    end do

    ! Get all possible combinatios until finding a valid octant index
    do i = 1,2
        do j = 1,2
            do k= 1,2
                ! ensure that major is on the 1,2,5, or 6 octant
                ! this ensures that the angles are correctly identified.
                if ((octmaj .eq. 1) .or. (octmaj .eq. 2) &
                    .or. (octmaj .eq. 5) .or. (octmaj .eq. 6)) then
                    continue
                    optvecmaj = vecmaj! keep value if already there
                else
                    optvecmaj = vecmaj*(-1)**i ! This gens opposite vectors
                end if

                optvecmin = vecmin*(-1)**j
                optvecter = vecter*(-1)**k
                ! Get new octants
                call getoctvector(optvecmaj, opt_octmaj)
                call getoctvector(optvecmin, opt_octmin)
                call getoctvector(optvecter, opt_octter)
                ! Get the index
                optoctindex = int(opt_octmaj*100+opt_octmin*10+&
                    opt_octter)
                ! Is this a valid set of eigenvectors?
                do l =1, 48
                    ref_valid_idx = valid_octs(l)
                    if (optoctindex .ne. ref_valid_idx) then
                        matcha = .false.
                        continue ! test all the values.
                    else ! if the index is correct, save values and
                        ! exit the subroutine
                        matcha = .true.
                        ! Octants
                        octmaj = opt_octmaj
                        octmin = opt_octmin
                        octter = opt_octter
                        ! Vectors
                        vecmaj = optvecmaj
                        vecmin = optvecmin
                        vecter = optvecter
                        go to 101
                    end if
                end do
            end do
        end do
    end do

101 continue
    right_rotmat(:,1) = vecmaj
    right_rotmat(:,2) = vecmin
    right_rotmat(:,3) = vecter
    return
    end subroutine get_octant



    end module gen_subs