    ! varmodel is a powerful, easy to use variogram modeling and fitting program.
    !
    ! Compile as a python DLL with the command:
    !  call f2py -c -m --fcompiler=gnu95 varmodel covasubs.for random.f90 varmodel.f90
    !
    !-----------------------------------------------------------------------------
    ! Copyright (c) 2014, Jared L. Deutsch, Matthew V. Deutsch and Clayton V. Deutsch
    ! All rights reserved.
    !
    ! Redistribution and use in source and binary forms, with or without modification,
    ! are permitted provided that the following conditions are met:
    !
    ! 1. Redistributions of source code must retain the above copyright notice,
    !    this list of conditions and the following disclaimer.
    !
    ! 2. Redistributions in binary form must reproduce the above copyright notice,
    !    this list of conditions and the following disclaimer in the documentation
    !    and/or other materials provided with the distribution.
    !
    ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    ! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    ! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    ! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
    ! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    ! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    ! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    ! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    ! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
    ! THE POSSIBILITY OF SUCH DAMAGE.
    !-----------------------------------------------------------------------------

    subroutine varmodelfit(nstvmf, c0vmf, itvmf, ccvmf, azmvmf, dipvmf, tiltvmf, ahmaxvmf, ahminvmf, avertvmf, sillvmf, &
                           maxitervmf, rseedvmf, minpairsvmf, npairswtvmf, invdistwtvmf, fixhmaxvertanisvmf, hmaxvertanisvmf, &
                           fixhminhmaxanisvmf, hminhmaxanisvmf, nvargsvmf, varlagdistvmf, varvaluevmf, varnumpairsvmf, varwtsvmf, &
                           varazmvmf, vardipvmf, optc0vmf, optitvmf, optccvmf, optazmvmf, optdipvmf, opttiltvmf, &
                           optahmaxvmf, optahminvmf, optavertvmf)
       !-----------------------------------------------------------------------
       !
       !             varmodelfit - Semi-automatic variogram fitting
       !             **********************************************
       !
       ! This program fits 3D directional variograms given a set of constraints
       ! and partially optimized parameters.
       !
       !  Parameters:
       !    - Partially optimized variogram parameters
       !       nstvmf - number of structures (cannot be optimized here, must be set)
       !       c0vmf(2) - min and max allowable values for the nugget effect
       !       itvmf(nstvmf,2) - min and max structure types allowed
       !       ccvmf(nstvmf,2) - min and max allowable values for structure contributions
       !       azmvmf(nstvmf,2) - min and max allowable values for ang1/azimuth
       !       dipvmf(nstvmf,2) - min and max allowable values for ang2/azimuth
       !       tiltvmf(nstvmf,2) - min and max allowable values for ang3/azimuth
       !       ahmaxvmf(nstvmf,2) - min and max allowable values for major range
       !       ahminvmf(nstvmf,2) - min and max allowable values for minor range
       !       avertvmf(nstvmf,2) - min and max allowable values for vertical range
       !       sillvmf(2) - min and max allowable sillvmf value
       !    - Optimization parameters
       !       maxitervmf - maximum number of iterations
       !       rseedvmf - random number seed
       !       minpairsvmf - minimum number of pairs to consider
       !       npairswtvmf - logical indicating whether or not to optimize with # pairs
       !       invdistwtvmf - logical indicating whether or not to use inverse distance weighting
       !       fixhmaxvertanisvmf - logical indicating whether or not to fix hmax/vert anisotropy
       !       hmaxvertanisvmf - hmax/vert anisotropy
       !       fixhminhmaxanisvmf - logical indicating whether or not to fix hmin/hmax anisotropy
       !       hminhmaxanisvmf - hmin/hmax anisotropy
       !    - Experimental variogram values
       !       nvargsvmf - number of experimental variogram points
       !       varlagdistvmf(nvargsvmf) - lag distances
       !       varvaluevmf(nvargsvmf) - variogram value
       !       varnumpairsvmf(nvargsvmf) - number of pairs found for this variogram value
       !       varazmvmf(nvargsvmf) - variogram azimuth
       !       vardipvmf(nvargsvmf) - variogram dips
       !    - Returns:
       !       optc0vmf - nugget effect
       !       optitvmf(nstvmf) - structure types
       !       optccvmf(nstvmf) - values for structure contributions
       !       optazmvmf(nstvmf) - values for ang1/azimuth
       !       optdipvmf(nstvmf) - values for ang2/azimuth
       !       opttiltvmf(nstvmf) - values for ang3/azimuth
       !       optahmaxvmf(nstvmf) - values for major range
       !       optahminvmf(nstvmf) - values for minor range
       !       optavertvmf(nstvmf) - values for vertical range
       !
       ! (C) Jared Deutsch 2014
       !-----------------------------------------------------------------------
       use random ! Mersenne-Twister random number generator
       implicit none
       real, parameter :: EPSLON = 1.0e-6
       real(kind=8), parameter :: BIGDBLE = 1d21, SMALLDBLE = 1d-6, PI = 4*atan(1.0d0)
       real(kind=8), parameter :: DEG2RAD = PI/180d0
       ! Parameters
       integer, intent(in) :: nstvmf, maxitervmf, nvargsvmf, minpairsvmf, rseedvmf
       real(kind=8), dimension(2), intent(in) :: c0vmf, sillvmf
       integer, dimension(nstvmf, 2), intent(in) :: itvmf
       real(kind=8), dimension(nstvmf, 2), intent(in) :: ccvmf, azmvmf, dipvmf, tiltvmf, ahmaxvmf, ahminvmf, avertvmf
       logical, intent(in) :: npairswtvmf, invdistwtvmf, fixhmaxvertanisvmf, fixhminhmaxanisvmf
       real(kind=8), intent(in) :: hmaxvertanisvmf, hminhmaxanisvmf
       integer, dimension(nvargsvmf), intent(in) :: varnumpairsvmf
       real(kind=8), dimension(nvargsvmf), intent(in) :: varlagdistvmf, varvaluevmf, varazmvmf, &
                                                         vardipvmf, varwtsvmf
       ! Return variables
       real(kind=8), intent(out) :: optc0vmf
       integer, dimension(nstvmf), intent(out) :: optitvmf
       real(kind=8), dimension(nstvmf), intent(out) :: optccvmf, optazmvmf, optdipvmf, opttiltvmf, &
                                                       optahmaxvmf, optahminvmf, optavertvmf
       ! Internal variables
       real(kind=8), allocatable, dimension(:) :: vrlagdistvmf, vrvaluevmf, vrazmvmf, vrdipvmf, vrwtsvmf
       real(kind=8), dimension(nvargsvmf) :: varvaluevmf1  ! to set values above the sill as one temporarily and model this
       real(kind=8) :: old_c0vmf, old_sillvmf, old_obj, opt_obj
       integer, dimension(nstvmf) :: old_itvmf
       real(kind=8), dimension(nstvmf) :: old_ccvmf, old_azmvmf, old_dipvmf, old_tiltvmf, &
                                          old_ahmaxvmf, old_ahminvmf, old_avertvmf
       integer :: i, j, iter, ist, ivarg, nvargsused, ipert, npoint
       real(kind=8), dimension(nvargsvmf) :: varwts, varfitvalue, suminvdist
       real(kind=8) :: npairstotal, rpert, randomshift
       real(kind=8) :: maxlagdist(3)
       logical, dimension(nvargsvmf) :: usevar
       logical :: fixc0, fixsill, success
       logical, dimension(nstvmf) :: fixcc
       ! Assign the temporary array to modify
       do i = 1, nvargsvmf
          if (varvaluevmf(i) .le. 1) then
             varvaluevmf1(i) = varvaluevmf(i)
          else
             varvaluevmf1(i) = 1
          end if
       end do
       ! Select only points that has a weigth greater than 0
       ! Calculate the new size of the array and allocate
       npoint = 0
       do i = 1, nvargsvmf
          if (varwtsvmf(i) .gt. 0.) npoint = npoint + 1
       end do
       allocate (vrlagdistvmf(npoint))
       allocate (vrvaluevmf, vrazmvmf, vrdipvmf, vrwtsvmf, mold=vrlagdistvmf)
       ! Fill the array
       j = 0
       do i = 1, nvargsvmf
          if (varwtsvmf(i) .gt. 0.) then
             j = j + 1
             vrlagdistvmf(j) = varlagdistvmf(i)
             vrvaluevmf(j) = varvaluevmf1(i)
             vrazmvmf(j) = varazmvmf(i)
             vrdipvmf(j) = vardipvmf(i)
             vrwtsvmf(j) = varwtsvmf(i)
          end if
       end do

       ! Initialize the random number generator
       call init_genrand(rseedvmf)
       ! Pre-calculate the sum inverse distance and total number of pairs
       npairstotal = 0
       suminvdist = SMALLDBLE
       do ivarg = 1, npoint
          if (vrlagdistvmf(ivarg) .gt. 0d0) then
             ! Only calculate sum inverse distance weight for the same direction
             do j = 1, npoint
                if ((vrazmvmf(ivarg) .eq. vrazmvmf(j)) .and. &
                    (vrdipvmf(ivarg) .eq. vrdipvmf(j)) .and. &
                    (vrlagdistvmf(j) .gt. 0d0)) then
                   suminvdist(ivarg) = suminvdist(ivarg) + 1d0/vrlagdistvmf(j)
                end if
             end do
          end if
          if (varnumpairsvmf(ivarg) .gt. 0) then
             npairstotal = npairstotal + varnumpairsvmf(ivarg)
          end if
       end do
       ! Figure out what variograms we are using to speed up calculations
       nvargsused = 0
       maxlagdist = SMALLDBLE
       do ivarg = 1, npoint
          if (varnumpairsvmf(ivarg) .ge. max(minpairsvmf, 1)) then
             usevar(ivarg) = .true.
             nvargsused = nvargsused + 1
          else
             usevar(ivarg) = .false.
          end if
       end do
       ! Have we got anything to work with?
       if (nvargsused .lt. 2) then
          write (*, *) 'Too few variogram values to even consider optimizing'
          write (*, *) 'Check that the minimum number of pairs specified is reasonable.'
          stop
       end if
       ! Weights for variogram points
       do ivarg = 1, npoint
          ! Initialize for equal weighting
          varwts(ivarg) = 1d0
          ! Inverse distance weighting of points
          if (invdistwtvmf) then
             if (vrlagdistvmf(ivarg) .gt. 0d0) then
                varwts(ivarg) = (varwts(ivarg)*1d0/vrlagdistvmf(ivarg)/suminvdist(ivarg))**2
             end if
          end if
          ! Is this a valid variogram point or are there too few pairs?
          if (varnumpairsvmf(ivarg) .lt. max(minpairsvmf, 1)) then
             varwts(ivarg) = 0
          else
             ! Number of pairs weighting
             if (npairswtvmf) then
                varwts(ivarg) = (varwts(ivarg)*varnumpairsvmf(ivarg)/npairstotal)
             end if
          end if
       end do
       !-----------------------------------------------------------------------
       ! Initial parameter estimates if not fixed
       !-----------------------------------------------------------------------
       ! Sill
       if (sillvmf(1) .ne. sillvmf(2)) then
          ! Estimate a reasonable sillvmf value to start with
          old_sillvmf = 0
          do ivarg = 1, npoint
             if (usevar(ivarg)) then
                old_sillvmf = old_sillvmf + vrvaluevmf(ivarg)
             end if
          end do
          old_sillvmf = old_sillvmf/nvargsused
          ! If this is less than sillvmf(1) and greater than sillvmf(2), pick the closest
          if (old_sillvmf .lt. sillvmf(1)) old_sillvmf = sillvmf(1)
          if (old_sillvmf .gt. sillvmf(2)) old_sillvmf = sillvmf(2)
          fixc0 = .false.
       else
          old_sillvmf = sillvmf(1)
          fixsill = .true.
       end if
       ! Nugget effect
       if (c0vmf(1) .ne. c0vmf(2)) then
          optc0vmf = 0
          if (optc0vmf .lt. c0vmf(1)) optc0vmf = c0vmf(1)
          if (optc0vmf .gt. c0vmf(2)) optc0vmf = c0vmf(2)
          fixc0 = .false.
       else
          optc0vmf = c0vmf(1)
          fixc0 = .true.
       end if
       ! Nested variance contributions
       ! Initial pass to determine fixed structures
       do ist = 1, nstvmf
          if (ccvmf(ist, 1) .ne. ccvmf(ist, 2)) then
             if (.not. fixsill) then
                ! No fixed sillvmf, so set to fraction of variability
                optccvmf(ist) = old_sillvmf/nstvmf
             else
                optccvmf(ist) = 0
             end if
             fixcc(ist) = .false.
          else
             optccvmf(ist) = ccvmf(ist, 1)
             fixcc(ist) = .true.
          end if
       end do
       ! Fix the variance contributions (sills)
       call fixsills(nstvmf, c0vmf, sillvmf, ccvmf, fixc0, fixcc, optc0vmf, optccvmf, success)
       if (.not. success) then
          write (*, *) 'WARNING: Unsuccessful sillvmf optimization on first pass'
          write (*, *) '  Are your constraints consistent and possible?'
       end if
       ! Structure types
       optitvmf = 1
       ! Angles
       optazmvmf = 0
       optdipvmf = 0
       opttiltvmf = 0
       ! Ranges
       do ist = 1, nstvmf
          if (ahmaxvmf(ist, 1) .ne. ahmaxvmf(ist, 2)) then
             ! Find the combination of closest angle/furthest distnace
             do ivarg = 1, npoint
                if ((vrazmvmf(ivarg) .eq. optazmvmf(ist)) .and. &
                    (vrdipvmf(ivarg) .eq. optdipvmf(ist))) then
                   if (vrlagdistvmf(ivarg) .gt. maxlagdist(1)) maxlagdist(1) = vrlagdistvmf(ivarg)
                end if
             end do
             if (maxlagdist(1) .le. SMALLDBLE) then
                do ivarg = 1, npoint
                   if (vrlagdistvmf(ivarg) .gt. maxlagdist(1)) maxlagdist(1) = vrlagdistvmf(ivarg)
                end do
             end if
             optahmaxvmf(ist) = real(ist)/real(nstvmf)*maxlagdist(1)
             if (optahmaxvmf(ist) .lt. ahmaxvmf(ist, 1)) optahmaxvmf(ist) = ahmaxvmf(ist, 1)
             if (optahmaxvmf(ist) .gt. ahmaxvmf(ist, 2)) optahmaxvmf(ist) = ahmaxvmf(ist, 2)
          else
             optahmaxvmf(ist) = ahmaxvmf(ist, 1)
          end if
          if (ahminvmf(ist, 1) .ne. ahminvmf(ist, 2)) then
             ! Find the combination of closest angle/furthest distnace
             do ivarg = 1, npoint
                if ((vrazmvmf(ivarg) + 90d0 .eq. optazmvmf(ist) .or. vrazmvmf(ivarg) - 90d0 .eq. optazmvmf(ist)) .and. &
                    (vrdipvmf(ivarg) .eq. optdipvmf(ist))) then
                   if (vrlagdistvmf(ivarg) .gt. maxlagdist(2)) maxlagdist(2) = vrlagdistvmf(ivarg)
                end if
             end do
             if (maxlagdist(2) .le. SMALLDBLE) then
                do ivarg = 1, npoint
                   if (vrlagdistvmf(ivarg) .gt. maxlagdist(2)) maxlagdist(2) = vrlagdistvmf(ivarg)
                end do
             end if
             optahminvmf(ist) = real(ist)/real(nstvmf)*maxlagdist(2)
             if (optahminvmf(ist) .lt. ahminvmf(ist, 1)) optahminvmf(ist) = ahminvmf(ist, 1)
             if (optahminvmf(ist) .gt. ahminvmf(ist, 2)) optahminvmf(ist) = ahminvmf(ist, 2)
          else
             optahminvmf(ist) = ahminvmf(ist, 1)
          end if
          if (avertvmf(ist, 1) .ne. avertvmf(ist, 2)) then
             ! Find the combination of closest angle/furthest distnace
             do ivarg = 1, npoint
                if ((vrazmvmf(ivarg) .eq. optazmvmf(ist)) .and. &
                    (vrdipvmf(ivarg) + 90d0 .eq. optdipvmf(ist) .or. vrdipvmf(ivarg) - 90d0 .eq. optdipvmf(ist))) then
                   if (vrlagdistvmf(ivarg) .gt. maxlagdist(3)) maxlagdist(3) = vrlagdistvmf(ivarg)
                end if
             end do
             if (maxlagdist(3) .le. SMALLDBLE) then
                do ivarg = 1, npoint
                   if (vrlagdistvmf(ivarg) .gt. maxlagdist(3)) maxlagdist(3) = vrlagdistvmf(ivarg)
                end do
             end if
             optavertvmf(ist) = real(ist)/real(nstvmf)*maxlagdist(3)
             if (optavertvmf(ist) .lt. avertvmf(ist, 1)) optavertvmf(ist) = avertvmf(ist, 1)
             if (optavertvmf(ist) .gt. avertvmf(ist, 2)) optavertvmf(ist) = avertvmf(ist, 2)
          else
             optavertvmf(ist) = avertvmf(ist, 1)
          end if
       end do
       !-----------------------------------------------------------------------
       ! Optimization
       !-----------------------------------------------------------------------
       ! Calculate initial objective function
       call varmodelpts(nstvmf, optc0vmf, optitvmf, optccvmf, optazmvmf, optdipvmf, &
                        opttiltvmf, optahmaxvmf, optahminvmf, optavertvmf, &
                        npoint, vrlagdistvmf, vrazmvmf, vrdipvmf, varfitvalue)
       opt_obj = 0
       do ivarg = 1, npoint
          opt_obj = opt_obj + varwts(ivarg)*(varfitvalue(ivarg) - vrvaluevmf(ivarg))**2
          !opt_obj = opt_obj + varwts(ivarg)*abs(varfitvalue(ivarg)-varvaluevmf(ivarg))
       end do
       write (*, *) ' '
       write (*, '(A,1f18.3)') 'Starting objective value = ', opt_obj
       ! Copy over the initial values
       old_obj = opt_obj
       old_c0vmf = optc0vmf
       old_itvmf = optitvmf
       old_ccvmf = optccvmf
       old_azmvmf = optazmvmf
       old_dipvmf = optdipvmf
       old_tiltvmf = opttiltvmf
       old_ahmaxvmf = optahmaxvmf
       old_ahminvmf = optahminvmf
       old_avertvmf = optavertvmf
       old_sillvmf = optc0vmf + sum(optccvmf)
       ! Structure type, range, sillvmf, angle optimization
       do iter = 1, maxitervmf
          ! Draw a random value to pick what to change
          rpert = grnd()
          ! Nugget effect and variance contribution change
          if (rpert .lt. 0.5d0) then
             ipert = floor(grnd()*(nstvmf) + 1)
             ! Random shift is +/- 7.5% of old sillvmf
             randomshift = 0.15d0*(grnd() - 0.5d0)*old_sillvmf
             optccvmf(ipert) = optccvmf(ipert) + randomshift

             ! Fix the variance contributions (sills)
             call fixsills(nstvmf, c0vmf, sillvmf, ccvmf, fixc0, fixcc, optc0vmf, optccvmf, success)

          else
             ! Range change
             ! Pick a structure number
             ist = floor(grnd()*nstvmf + 1)
             ! Pick a range value - ahmaxvmf,ahminvmf or avertvmf
             ipert = floor(grnd()*3 + 1)
             ! Random shift is +/- 7.5% of maximum range in this direction
             randomshift = 0.15d0*(grnd() - 0.5d0)
             if (ipert .eq. 3) then
                optavertvmf(ist) = optavertvmf(ist) + randomshift*maxval(optavertvmf)
                if (ist .lt. nstvmf) then
                   if (optavertvmf(ist) .gt. optavertvmf(ist + 1)) optavertvmf(ist) = optavertvmf(ist + 1)
                end if
                if (optavertvmf(ist) .lt. avertvmf(ist, 1)) optavertvmf(ist) = avertvmf(ist, 1)
                if (optavertvmf(ist) .gt. avertvmf(ist, 2)) optavertvmf(ist) = avertvmf(ist, 2)
                if (fixhmaxvertanisvmf) then
                   optahmaxvmf(ist) = optavertvmf(ist)*hmaxvertanisvmf
                end if
             elseif (ipert .eq. 2) then
                optahminvmf(ist) = optahminvmf(ist) + randomshift*maxval(optahminvmf)
                if (ist .lt. nstvmf) then
                   if (optahminvmf(ist) .gt. optahminvmf(ist + 1)) optahminvmf(ist) = optahminvmf(ist + 1)
                end if
                if (optahminvmf(ist) .lt. ahminvmf(ist, 1)) optahminvmf(ist) = ahminvmf(ist, 1)
                if (optahminvmf(ist) .gt. ahminvmf(ist, 2)) optahminvmf(ist) = ahminvmf(ist, 2)
                if (fixhminhmaxanisvmf) then
                   optahmaxvmf(ist) = optahminvmf(ist)/hminhmaxanisvmf
                end if
             else
                optahmaxvmf(ist) = optahmaxvmf(ist) + randomshift*maxval(optahmaxvmf)
                if (ist .lt. nstvmf) then
                   if (optahmaxvmf(ist) .gt. optahmaxvmf(ist + 1)) optahmaxvmf(ist) = optahmaxvmf(ist + 1)
                end if
                if (optahmaxvmf(ist) .lt. ahmaxvmf(ist, 1)) optahmaxvmf(ist) = ahmaxvmf(ist, 1)
                if (optahmaxvmf(ist) .gt. ahmaxvmf(ist, 2)) optahmaxvmf(ist) = ahmaxvmf(ist, 2)
                if (fixhmaxvertanisvmf) then
                   optavertvmf(ist) = optahmaxvmf(ist)/hmaxvertanisvmf
                end if
                if (fixhminhmaxanisvmf) then
                   optahminvmf(ist) = optahmaxvmf(ist)*hminhmaxanisvmf
                end if
             end if
             !Enforce order relationships:
             !Correct to have the minor always less than the major
             if (optahminvmf(ist) .gt. optahmaxvmf(ist)) optahminvmf(ist) = optahmaxvmf(ist)
             if (optavertvmf(ist) .gt. optahminvmf(ist)) optavertvmf(ist) = optahminvmf(ist)

          end if
          ! Enforce range order parameters
          do ist = 1, nstvmf
             if (optahminvmf(ist) .gt. optahmaxvmf(ist)) then
                optahminvmf(ist) = optahmaxvmf(ist)
             end if
             if (optavertvmf(ist) .gt. optahmaxvmf(ist)) then
                optavertvmf(ist) = optahmaxvmf(ist)
             end if
             if (optavertvmf(ist) .gt. optahminvmf(ist)) then
                optavertvmf(ist) = optahminvmf(ist)
             end if
          end do

          ! Calculate a new objective value
          call varmodelpts(nstvmf, optc0vmf, optitvmf, optccvmf, optazmvmf, optdipvmf, &
                           opttiltvmf, optahmaxvmf, optahminvmf, optavertvmf, &
                           npoint, vrlagdistvmf, vrazmvmf, vrdipvmf, varfitvalue)
          opt_obj = 0
          do ivarg = 1, npoint
             opt_obj = opt_obj + varwts(ivarg)*(varfitvalue(ivarg) - vrvaluevmf(ivarg))**2
          end do
          ! Accept?
          if (opt_obj .lt. old_obj) then
             ! Accept the iteration
             old_obj = opt_obj
             old_c0vmf = optc0vmf
             old_itvmf = optitvmf
             old_ccvmf = optccvmf
             old_azmvmf = optazmvmf
             old_dipvmf = optdipvmf
             old_tiltvmf = opttiltvmf
             old_ahmaxvmf = optahmaxvmf
             old_ahminvmf = optahminvmf
             old_avertvmf = optavertvmf
             old_sillvmf = optc0vmf + sum(optccvmf)
          else
             ! Reject the iteration
             opt_obj = old_obj
             optc0vmf = old_c0vmf
             optitvmf = old_itvmf
             optccvmf = old_ccvmf
             optazmvmf = old_azmvmf
             optdipvmf = old_dipvmf
             opttiltvmf = old_tiltvmf
             optahmaxvmf = old_ahmaxvmf
             optahminvmf = old_ahminvmf
             optavertvmf = old_avertvmf
          end if
       end do
       write (*, '(A,1f21.3)') 'Final objective value = ', opt_obj
       ! Final check to fix sills before returning
       call fixsills(nstvmf, c0vmf, sillvmf, ccvmf, fixc0, fixcc, optc0vmf, optccvmf, success)
       if (.not. success) then
          write (*, *) 'ERROR: Unsuccessful sillvmf optimization on final pass'
          write (*, *) '  Are your constraints consistent and possible?'
       end if
       !If 2D dataset: fix optavertvmf as zero
       if ((avertvmf(3, 1) .eq. 0) .and. (avertvmf(3, 2) .eq. 0)) optavertvmf = 0
       return
    end subroutine varmodelfit

    subroutine varmodelpts(nstpts, c0pts, itpts, ccpts, azmpts, dippts, tiltpts, ahmaxpts, ahminpts, avertpts, &
                           nvargspts, varlagdistpts, varazmpts, vardippts, varmodelvalspts)
       !-----------------------------------------------------------------------
       ! This subroutine calculates variogram model points at input lags with
       ! provided lag distances, azimuths and dips using the GSLIB cova3
       ! subroutine.
       !
       ! (C) Jared Deutsch 2014
       !-----------------------------------------------------------------------
       implicit none
       real(kind=8), parameter :: SMALLDBLE = 1d-6, PI = 4*atan(1.0d0)
       real(kind=8), parameter :: DEG2RAD = PI/180d0
       ! Parameters
       integer, intent(in) :: nstpts, nvargspts
       real(kind=8), intent(in) :: c0pts
       integer, dimension(nstpts), intent(in) :: itpts
       real(kind=8), dimension(nstpts), intent(in) :: ccpts, azmpts, dippts, tiltpts, ahmaxpts, ahminpts, avertpts
       real(kind=8), dimension(nvargspts), intent(in) :: varlagdistpts, varazmpts, vardippts
       ! Return
       real(kind=8), dimension(nvargspts), intent(out) :: varmodelvalspts
       ! Internal variables
       integer :: ist, ivarg, i, j
       ! GSLIB interface variables
       integer :: MAXROT
       real(kind=8), allocatable :: rotmat(:, :, :)
       real :: cmax, maxcov, cova, x, y, z, c0gslib
       real, dimension(nstpts) :: ccgslib, aagslib, anis1, anis2
       ! GSLIB interface allocation
       MAXROT = nstpts
       allocate (rotmat(MAXROT, 3, 3))
       c0gslib = real(c0pts)
       ! GSLIB anisotropy ratios
       do ist = 1, nstpts
          anis1(ist) = real(ahminpts(ist)/max(ahmaxpts(ist), SMALLDBLE))
          anis2(ist) = real(avertpts(ist)/max(ahmaxpts(ist), SMALLDBLE))
          ccgslib(ist) = real(ccpts(ist))
          aagslib(ist) = real(ahmaxpts(ist))
       end do
       ! GSLIB rotation matrix
       rotmat = 0
       do ist = 1, nstpts
          call setrot_(real(azmpts(ist)), real(dippts(ist)), real(tiltpts(ist)), anis1(ist), anis2(ist), &
                       ist, MAXROT, rotmat)
       end do
       ! Maximum covariance for calculating variogram values from cova3
       call cova33(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, nstpts, nstpts, c0gslib, &
                   itpts, ccgslib(1:nstpts), aagslib, 1, MAXROT, rotmat, cmax, maxcov)
       ! Calculate variogram values
       do ivarg = 1, nvargspts
          x = real(sin(DEG2RAD*varazmpts(ivarg))*cos(DEG2RAD*vardippts(ivarg))*varlagdistpts(ivarg))
          y = real(cos(DEG2RAD*varazmpts(ivarg))*cos(DEG2RAD*vardippts(ivarg))*varlagdistpts(ivarg))
          z = real(sin(DEG2RAD*vardippts(ivarg))*varlagdistpts(ivarg))
          call cova33(0.0, 0.0, 0.0, x, y, z, 1, nstpts, nstpts, c0gslib, &
                      itpts, ccgslib(1:nstpts), aagslib, 1, MAXROT, rotmat, cmax, cova)
          varmodelvalspts(ivarg) = dble(maxcov - cova)
       end do
       return
    end subroutine varmodelpts

    subroutine fixsills(nstfs, c0fs, sillfs, ccfs, fixc0fs, fixccfs, optc0fs, optccfs, successfs)
       !-----------------------------------------------------------------------
       ! This subroutine enforces the constraints imposed by the user on sills
       ! if possible. As usual, garbage in = garbage out. No excessive error
       ! checking is performed, but maybe it should be.
       !
       ! (C) Jared Deutsch 2014
       !-----------------------------------------------------------------------
       implicit none
       integer, parameter :: maxiter = 100
       real(kind=8), parameter :: abstol = 1.0d-6
       integer, intent(in) :: nstfs
       real(kind=8), dimension(2), intent(in) :: c0fs, sillfs
       real(kind=8), dimension(nstfs, 2), intent(in) :: ccfs
       logical, intent(in) :: fixc0fs
       logical, dimension(nstfs), intent(in) :: fixccfs
       logical, intent(out) :: successfs
       real(kind=8), intent(inout) :: optc0fs
       real(kind=8), dimension(nstfs), intent(inout) :: optccfs
       real(kind=8) :: unfixedvar, inputsill, inputunfixed, rtol, closestsill
       integer :: ist, niter
       logical sillsfixed
       ! Counter initialization
       sillsfixed = .false.
       rtol = abstol*sillfs(1)
       niter = 0
       successfs = .true.
       ! Main iterative loop to correct sills
       do while (.not. sillsfixed)
          ! Make sure that constraints are enforced and calculate input sillfs
          inputsill = 0d0
          if (optc0fs .lt. c0fs(1)) optc0fs = c0fs(1)
          if (optc0fs .gt. c0fs(2)) optc0fs = c0fs(2)
          inputsill = inputsill + optc0fs
          do ist = 1, nstfs
             if (optccfs(ist) .lt. ccfs(ist, 1)) optccfs(ist) = ccfs(ist, 1)
             if (optccfs(ist) .gt. ccfs(ist, 2)) optccfs(ist) = ccfs(ist, 2)
             inputsill = inputsill + optccfs(ist)
          end do
          ! Are we done?
          if ((sillfs(1) .eq. sillfs(2)) .and. &
              (abs(inputsill - sillfs(1)) .lt. rtol)) then
             sillsfixed = .true.
          elseif ((inputsill .ge. sillfs(1)) .and. &
                  (inputsill .le. sillfs(2))) then
             sillsfixed = .true.
          else
             ! Not done - try and fix - what sillfs constraint are we closest to?
             if (abs(sillfs(2) - inputsill) .lt. abs(sillfs(1) - inputsill)) then
                closestsill = sillfs(2)
             else
                closestsill = sillfs(1)
             end if
             ! How much variability is not fixed?
             unfixedvar = closestsill
             inputunfixed = 0d0
             if (fixc0fs) then
                unfixedvar = unfixedvar - c0fs(1)
             else
                inputunfixed = inputunfixed + optc0fs
             end if
             do ist = 1, nstfs
                if (fixccfs(ist)) then
                   unfixedvar = unfixedvar - ccfs(ist, 1)
                else
                   inputunfixed = inputunfixed + optccfs(ist)
                end if
             end do
             ! Scale all the variable sills
             if (inputunfixed .gt. 0d0) then
                if (.not. fixc0fs) optc0fs = optc0fs*unfixedvar/inputunfixed
                do ist = 1, nstfs
                   if (.not. fixccfs(ist)) then
                      optccfs(ist) = optccfs(ist)*unfixedvar/inputunfixed
                   end if
                end do
             end if
          end if
          ! Maximum iteration condition
          niter = niter + 1
          if (niter .gt. maxiter) then
             sillsfixed = .true.
             successfs = .false.
          end if
          ! End main iterative loop
       end do
       ! Seriously, double check these constraints again - these dominate
       if (optc0fs .lt. c0fs(1)) optc0fs = c0fs(1)
       if (optc0fs .gt. c0fs(2)) optc0fs = c0fs(2)
       do ist = 1, nstfs
          if (optccfs(ist) .lt. ccfs(ist, 1)) optccfs(ist) = ccfs(ist, 1)
          if (optccfs(ist) .gt. ccfs(ist, 2)) optccfs(ist) = ccfs(ist, 2)
       end do
       return
    end subroutine fixsills

    subroutine varmodelsph(nstsph, optc0sph, optitsph, optccsph, optahmaxsph, &
                           optahminsph, optavertsph, pointsepsph, varmodelsphsph, lagdistsph)
       ! Subroutine to calculate a spherical variogram model for 1000 lag distances
       ! in each direction
       implicit none

       integer, intent(in) :: nstsph
       real(kind=8), dimension(3), intent(in) :: optccsph, optahmaxsph, optahminsph, optavertsph
       integer, dimension(3), intent(in) :: optitsph
       integer :: i
       real(kind=8), dimension(3), intent(in) :: pointsepsph
       real(kind=8) :: optc0sph, hr1, hr2, hr3
       real(kind=8), dimension(1000) :: dist1, dist2, dist3, sph1, sph2, sph3, &
                                        LMC_major, LMC_minor, LMC_ter
       real(kind=8), dimension(3000), intent(out) :: varmodelsphsph, lagdistsph

       ! Model in major direction
       do i = 1, 1000
          dist1(i) = pointsepsph(1)*i
          hr1 = dist1(i)/optahmaxsph(1)
          if (hr1 .lt. 1) then
             sph1(i) = 1.5*(hr1) - 0.5*(hr1)**3
          else
             sph1(i) = 1
          end if

          hr2 = dist1(i)/optahmaxsph(2)
          if (hr2 .lt. 1) then
             sph2(i) = 1.5*(hr2) - 0.5*(hr2)**3
          else
             sph2(i) = 1
          end if

          hr3 = dist1(i)/optahmaxsph(3)
          if (hr3 .lt. 1) then
             sph3(i) = 1.5*(hr3) - 0.5*(hr3)**3
          else
             sph3(i) = 1
          end if

          LMC_major(i) = optccsph(1)*sph1(i) + optccsph(2)*sph2(i) + &
                         optccsph(3)*sph3(i)
       end do

       ! Model in major direction
       do i = 1, 1000
          dist2(i) = pointsepsph(2)*i
          hr1 = dist2(i)/optahminsph(1)
          if (hr1 .lt. 1) then
             sph1(i) = 1.5*(hr1) - 0.5*(hr1)**3
          else
             sph1(i) = 1
          end if

          hr2 = dist2(i)/optahminsph(2)
          if (hr2 .lt. 1) then
             sph2(i) = 1.5*(hr2) - 0.5*(hr2)**3
          else
             sph2(i) = 1
          end if

          hr3 = dist2(i)/optahminsph(3)
          if (hr3 .lt. 1) then
             sph3(i) = 1.5*(hr3) - 0.5*(hr3)**3
          else
             sph3(i) = 1
          end if

          LMC_minor(i) = optccsph(1)*sph1(i) + optccsph(2)*sph2(i) + &
                         optccsph(3)*sph3(i)

       end do

       ! Model in tertiary direction
       do i = 1, 1000
          dist3(i) = pointsepsph(3)*i
          hr1 = dist3(i)/optavertsph(1)
          if (hr1 .lt. 1) then
             sph1(i) = 1.5*(hr1) - 0.5*(hr1)**3
          else
             sph1(i) = 1
          end if

          hr2 = dist3(i)/optavertsph(2)
          if (hr2 .lt. 1) then
             sph2(i) = 1.5*(hr2) - 0.5*(hr2)**3
          else
             sph2(i) = 1
          end if

          hr3 = dist3(i)/optavertsph(3)
          if (hr3 .lt. 1) then
             sph3(i) = 1.5*(hr3) - 0.5*(hr3)**3
          else
             sph3(i) = 1
          end if

          LMC_ter(i) = optccsph(1)*sph1(i) + optccsph(2)*sph2(i) + &
                       optccsph(3)*sph3(i)
       end do

       varmodelsphsph(1:1000) = LMC_major
       varmodelsphsph(1001:2000) = LMC_minor
       varmodelsphsph(2001:3000) = LMC_ter
       lagdistsph(1:1000) = dist1
       lagdistsph(1001:2000) = dist2
       lagdistsph(2001:3000) = dist3

       return
    end subroutine varmodelsph

