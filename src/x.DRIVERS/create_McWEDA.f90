! copyright info:
!
!                             @Copyright 2022
!                           Fireball Committee
! Hong Kong Quantum AI Laboratory, Ltd. - James P. Lewis, Chair
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek
! Arizona State University - Otto F. Sankey

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! California Institute of Technology - Brandon Keith
! Czech Institute of Physics - Prokop Hapala
! Czech Institute of Physics - Vladimír Zobač
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
! Synfuels China Technology Co., Ltd. - Pengju Ren
! Washington University - Pete Fedders
! West Virginia University - Ning Ma and Hao Wang
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman
!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! Program Description
! ===========================================================================
!       This is the main driver for CREATE.

! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1430 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        program create

! /GLOBAL
        use M_welcome 

! /SYSTEM
        use M_species
        use M_atom_functions
        use M_integrals_2c

! /HARRIS
        use M_overlap
        use M_kinetic
        use M_vna_Harris
        use M_vnl
        use M_rho_2c_Harris
        use M_rhoS_2c_Harris
        use M_vxc_Harris
        use M_Coulomb

! /DOGS
       use M_dipole_z
       use M_vna_DOGS
       use M_vxc_DOGS

! /NAC
        use M_Goverlap

! /HARRIS 3C
        use M_bcna_Harris
        use M_rho_3c_Harris
        use M_rhoS_3c_Harris

! /DOGS 3C
        use M_bcna_DOGS

        implicit none

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies, kspecies        ! counters over species

! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end

        character (len = 25) interactions

        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle_2c
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle_3c

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! ===========================================================================
! ---------------------------------------------------------------------------
!                             W E L C O M E !
! ---------------------------------------------------------------------------
! ===========================================================================
        call cpu_time (time_begin)
        open (unit = ilogfile, file = 'output.log', status = 'replace')
        call welcome_create

! Read in create.inp
        write (ilogfile,*) ' Reading in create.inp '
        call read_Fdata_location
        allocate (species_PP (nspecies))
        call read_create

! Read in wavefunction
        write (ilogfile,*)
        write (ilogfile,*) ' Reading wavefunctions and potentials. '
        call read_wavefunctions
        call read_napotentials
        call read_vpp

! Write the info.dat file
        call write_info

! ***************************************************************************
!                T W O - C E N T E R    I N T E R A C T I O N S
! ***************************************************************************
! Initialize the size of the two-center bundles
        write (ilogfile,*)
        write (ilogfile,*) ' Sizing the Fdata_2c arrays. '
        allocate (Fdata_bundle_2c (nspecies, nspecies))
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            ! cut some lengthy notation
            pFdata_bundle_2c=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle_2c%nFdata_cell_2c = 0
          end do
        end do

! Very important!
! The order of the initializers here much match the order that routines
! are called below.
! Initialize two-center Harris routines
        call initialize_overlap
        call initialize_kinetic
        call initialize_vna_Harris
        call initialize_vnl
        call initialize_rho_2c_Harris
        call initialize_rhoS_2c_Harris
        call initialize_vxc_Harris

! Initialize two-center DOGS routines
        call initialize_dipole_z
        call initialize_vna_DOGS
        call initialize_vxc_DOGS

! Initialize interactions for short-range (double-counting) corrections
        call initialize_Coulomb

! Allocate two-center array sizes
        call size_Fdata_2c

! The ordering here much match the ordering of the initializers below.
! Calculate the two-center Harris interactions:
        call overlap
        call kinetic
        call vna_Harris
        call vnl_2c
        call rho_2c_Harris
        call rhoS_2c_Harris
        call vxc_Harris

! Calculate the two-center DOGS interactions:
        call dipole_z
        call vna_DOGS
        call vxc_DOGS

! Calculate interactions for short-range (double-counting) corrections:
        call Coulomb

! Write out the number of interactions
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            ! cut some lengthy notation
            pFdata_bundle_2c=>Fdata_bundle_2c(ispecies, jspecies)

            ! open directory file
            write (interactions,'("/2c.",i2.2,".",i2.2,".dat")')             &
      &       species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
      &           status = 'unknown')
            write (13,*) pFdata_bundle_2c%nFdata_cell_2c
            close (unit = 13)
          end do
        end do

! ***************************************************************************
!             T H R E E - C E N T E R    I N T E R A C T I O N S
! ***************************************************************************
! Initialize the size of the three-center bundles
        write (ilogfile,*)
        write (ilogfile,*)
        write (ilogfile,*) ' Sizing the Fdata_3c arrays. '
        allocate (Fdata_bundle_3c (nspecies, nspecies, nspecies))
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            do kspecies = 1, nspecies
              ! cut some lengthy notation
              pFdata_bundle_3c=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
              pFdata_bundle_3c%nFdata_cell_3c = 0
            end do
          end do
        end do

! Initialize three-center Harris routines
        call initialize_rho_3c_Harris
        call initialize_rhoS_3c_Harris
        call initialize_bcna_Harris

! Initialize three-center DOGS routines
        call initialize_bcna_DOGS

! Allocate two-center array sizes
        call size_Fdata_3c

! Now calculate the three-center Harris interactions:
        call rho_3c_Harris
        call rhoS_3c_Harris
        call bcna_Harris

! Now calculate the three-center DOGS interactions:
        call bcna_DOGS

! Write out the number of interactions
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            do kspecies = 1, nspecies
              ! cut some lengthy notation
              pFdata_bundle_3c=>Fdata_bundle_3c(ispecies, jspecies, kspecies)

              ! open directory file
              write (interactions,'("/3c.",i2.2,".",i2.2,".",i2.2,".dat")')  &
     &          species(ispecies)%nZ, species(jspecies)%nZ, species(kspecies)%nZ
              open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &              status = 'unknown')
              write (13,*) pFdata_bundle_3c%nFdata_cell_3c
              close (unit = 13)
            end do
          end do
        end do

        ! A "destroy" for the Fdata needs to be made here.
        deallocate (Fdata_bundle_2c)
        write (ilogfile,*)
        write (ilogfile,*) ' Create is Done! '

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Program
! ===========================================================================
        call cpu_time (time_end)
        write (ilogfile,*) time_begin
        write (ilogfile,*) time_end

        stop
        end program create
