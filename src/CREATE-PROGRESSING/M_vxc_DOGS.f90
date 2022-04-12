! copyright info:
!
!                             @Copyright 2009
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! Caltech - Brandon Keith
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
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

! M_vxc_Harris
! Program Description
! ===========================================================================
!      This is a module calculating the integrals of two centers for the
! the exchange-correlation interactions.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute, Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! ===========================================================================
! Module Declaration
! ===========================================================================
        module M_vxc_Harris
        use M_atom_functions
        use M_species
        use M_integrals_2c

! Type Declaration
! ===========================================================================
! This type contains the density and derivatives on rho, z grid
        type T_rho_2c_store
          real, pointer :: rho_2c (:, :)   ! The density - function of r and z
          real, pointer :: rhop_2c (:, :)  ! derivative with respect to rho
          real, pointer :: rhopp_2c (:, :) ! second derivative with respect to rho
          real, pointer :: rhoz_2c (:, :)  ! derivative with respect to z
          real, pointer :: rhozz_2c (:, :) ! second derivative with respect to z
          real, pointer :: rhopz_2c (:, :) ! mixed derivative - rho and z
        end type T_rho_2c_store

! the grid that contains the densities
        type T_rho_bundle_2c

          ! the density on grid of d (distance between two centers)
          type (T_rho_2c_store), pointer :: rho_2c_store (:)
        end type T_rho_bundle_2c

! stored densities into bundles based on ispecies and jspecies pair
        type(T_rho_bundle_2c), pointer :: rho_bundle_2c (:, :)

! module procedures
        contains

! ===========================================================================
! initialize_vxc_Harris
! ===========================================================================
! Program Description
! ===========================================================================
!       We need to determine how many interactions belong to each nspecies
! bundle pair. This routine just counts how many total interactions contribute
! to that bundle.  Something like overlap is obviously only 1 interaction
! added, but something like vna needs number of interactions based on the
! number of shells.
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine initialize_vxc_Harris
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species

        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 3
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_vxc_Harris


! ===========================================================================
! vxc_Harris
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine calls the subroutines required to calculate the vna_ontop
! (both left/right cases) and atom cases.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine vxc_Harris
        implicit none

! Parameters and Data Declaration
! ===========================================================================
! None

! ===========================================================================
! Argument Declaration and Description
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
        write (*,*)
        write (*,*) ' ******************************************************* '
        write (*,*) '        E X C H A N G E   C O R R E L A T I O N          '
        write (*,*) '                  I N T E R A C T I O N S                '
        write (*,*) ' ******************************************************* '
        write (*,*)

        write (*,*) ' Building the two center density on grid '
        call rho_2c_store

!        write (*,*) ' Calling ontop left case. '
!        call vna_ontopL_Harris

!        write (*,*)
!        write (*,*) ' Calling ontop right case. '
!        call vna_ontopR_Harris

!        write (*,*)
!        write (*,*) ' Calling atom case. '
!        call vna_atom_Harris

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ==========================================================================
        return
        end subroutine vxc_Harris


! ===========================================================================
! rho_2c_store
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This routine calculates and stores the combined density of two species
! as a function of r, z and d. This is important because later the gradients
! with respect to these variables are calculated.
!
! On output: rho_2c:      The density as a function of species type,
!                         r, z and d.  Output is placed in common block
!                         density located in wavefunctions.inc
!            rhop_2c:     derivative with respect to rho
!            rhopp_2c:    second derivative with respect to rho
!            rhoz_2c:     derivative with respect to z
!            rhozz_2c:    second derivative with respect to rho
!            rhopz_2c:    mixed derivative with respect to rho and z
!
! ====================================================================
! Code written by:
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ====================================================================
!
! Program Declaration
! ====================================================================
	    subroutine rho2c_store
   	    implicit none

! Argument Declaration and Description
! ===========================================================================


! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species
  	    integer iexc
        integer igrid                       !< number of grid points

        real d                              !< distance between the two centers

        real rcutoff1, rcutoff2             !< cutoffs for the two centers

        real zmin, zmax
        real rhomin, rhomax

        type (T_rho_2c_bundle), pointer :: prho_bundle
        type (T_rho_2c_store), pointer :: prho_2c







! Local Parameters and Data Declaration
! ====================================================================
        integer j1(5)
        data j1 /1, 1, 1, 2, 2/

        integer j2(5)
        data j2 /0, -1, 1, -1, 1/








! Allocate Arrays
! ===========================================================================
        allocate (rho_bundle_2c (nspecies, nspecies))

! Procedure
! ===========================================================================
! Set which exchange-correlation we are calculating
        ispecies = 1
      	iexc = PP_species(ispecies)%iexc

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            prho_bundle=>rho_bundle_2c(ispecies, jspecies)

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoffA_max
            dmax = wf(ispecies)%rcutoffA_max + wf(jspecies)%rcutoffA_max
            drr = dmax/float(ndd_vxc - 1)
            d = -drr

            ! Set integration limits
	        rhomin = 0.0d0
    	    rhomax = min(rcutoff1, rcutoff2)

! Loop over grid
            write (*,100) species(ispecies)%nZ, species(jspecies)%nZ
            allocate (prho_bundle%rho_2c_store(ndd_vxc))
            do igrid = 1, ndd_vxc
              prho_2c => prho_bundle%rho_2c_store(igrid)
              d = d + drr

              ! Set integration limits
	          zmin = max(-rcutoff1, d - rcutoff2)
    	      zmax = min(rcutoff1, d + rcutoff2)

              dz = ((rcutoff1 + rcutoff2)/2.0d0)/dfloat(nz_rho_store)
              nnz = int(zmax - zmin)/dz)
              if (mod(nnz,2) .eq. 0) nnz = nnz + 1

              drho = dz
              nnrho = int(rhomax - rhomin)/drho
              if (mod(nnrho,2) .eq. 0) nnrho = nnrho + 1

! Allocate all of the rho_2c arrays
              allocate (prho_2c%rho_2c (nnrho, nnz))
              allocate (prho_2c%rhop_2c (nnrho, nnz))
              allocate (prho_2c%rhopp_2c (nnrho, nnz))
              allocate (prho_2c%rhoz_2c (nnrho, nnz))
              allocate (prho_2c%rhozz_2c (nnrho, nnz))
              allocate (prho_2c%rhopz_2c (nnrho, nnz))








!Set up dqorb and iderorb.iderorb is the shell where
! the charge will be changed in the +-Q derivative stuff, and dqorb is
! the charge amount changed.
	do aspecies = 1, nspecies
    		iderorb(aspecies) = species(aspecies)%nssh
  			dqorb(ispecies) = 0.5d0
		if (species(aspecies)%nssh .eq. 1) dqorb(aspecies) = 0.25d0
		do issh = 1, species(aspecies)%nssh
  			dqint(issh,aspecies) = dqorb(aspecies)/species(aspecies)%nssh
		end do
	end do




! Initialize 1/4pi
    xinv4pi = 1.0d0/(4.0d0*3.141592653589793238462643D0)


! First initialize in (index) array. This is for the charge correction
! calculated later.
        in(1) = ispecies
        in(2) = jspecies

        j1ch = j1(ideriv)
        j1at = in(j1ch)
        jssh = iderorb(j1at)








! **************************************************************************
! Here we loop over z and r computing the sum of the densities
! for species in1 and in2 at each value of d, z and r.
              do iz = 1, nnz_rho_store
                z1 = zmin + (iz - 1)*dz
                z2 = z1 - d

                do irho = 1, nnrho
                  rho = rhomin + (irho - 1)*drho
                  r1 = sqrt(z1**2 + rho**2)
                  r2 = sqrt(z2**2 + rho**2)

! Evaluate the density at this grid point.
                  density = 0.0d0
                  do issh = 1, species(ispecies)%nssh
                    xnocc = species(ispecies)%shell(issh)%xnocc
                    density = density + xnocc*psiofr(r1, ispecies, issh)**2
                  end do

                  do kssh = 1, species(jspecies)%nssh
                    xnocc = species(jspecies)%shell(kssh)%xnocc
                    density = density + xnocc*psiofr(r2, jspecies,kssh)**2
                  end do






! Here the derivative with respect to the charge correction term
! is calculated.
          if (j1ch .eq. 1) r = r1
          if (j1ch .eq. 2) r = r2

          dens = dens +j2(ideriv)*dqorb(j1at)*psiofr(r,j1at,jssh)**2
          rhostore(ispecies, jspecies)%rho2c(irho,iz) = dens*xinv4pi
         end do
        end do

! **************************************************************************
! Now calculate the derivatives
! Only calculate the derivatives if doing GGA exchange-correlation.
              if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6	             &
     &            .or. iexc .eq. 9 .or. iexc .eq. 10) then

! Calculate rhop2c and rhopp2c.
              do iz = 1, nnz
                do irho = 2, nnrho - 1

! First derivative:
           rhostore(ispecies, jspecies)%rhop2c(irho,iz) =	&
           (rhostore(ispecies, jspecies)%rho2c(irho+1,iz) - &
     			rhostore(ispecies, jspecies)%rho2c(irho-1,iz))/(2.0d0*drho)
! Second derivative:
           rhostore(ispecies, jspecies)%rhopp2c(irho,iz) =	&
           (rhostore(ispecies, jspecies)%rho2c(irho+1,iz) - &
     				2.0d0*rhostore(ispecies, jspecies)%rho2c(irho,iz) &
            + rhostore(ispecies, jspecies)%rho2c(irho-1,iz))/(drho**2)
          end do

! Find endpoint values for the derivatives calculated above.
          rhostore(ispecies, jspecies)%rhop2c(1,iz) = 2.0d0* &
          			rhostore(ispecies, jspecies)%rhop2c(2,iz) - &
          					rhostore(ispecies, jspecies)%rhop2c(3,iz)
          rhostore(ispecies, jspecies)%rhop2c(nnrho,iz) = 2.0d0* &
          			rhostore(ispecies, jspecies)%rhop2c(nnrho-1,iz) &
                           - rhostore(ispecies, jspecies)%rhop2c(nnrho-2,iz)

          rhostore(ispecies, jspecies)%rhopp2c(1,iz) = 2.0d0* &
          			rhostore(ispecies, jspecies)%rhopp2c(2,iz) - &
          					rhostore(ispecies, jspecies)%rhopp2c(3,iz)
          rhostore(ispecies, jspecies)%rhopp2c(nnrho,iz) = 2.0d0* &
          			rhostore(ispecies, jspecies)%rhopp2c(nnrho-1,iz) &
                            - rhostore(ispecies, jspecies)%rhopp2c(nnrho-2,iz)
         end do

! Calculate rhoz2c and rhozz2c.
         do irho = 1, nnrho
          do iz = 2, nnz - 1

! First derivative:
           rhostore(ispecies, jspecies)%rhoz2c(irho,iz) =	&
           	(rhostore(ispecies, jspecies)%rho2c(irho,iz+1) - &
           		rhostore(ispecies, jspecies)%rho2c(irho,iz-1))/(2.0d0*dz)

! Second derivative:
           rhostore(ispecies, jspecies)%rhozz2c(irho,iz) =	&
 	          (rhostore(ispecies, jspecies)%rho2c(irho,iz+1) - &
    			 2.0d0*rhostore(ispecies, jspecies)%rho2c(irho,iz) &
            	+ rhostore(ispecies, jspecies)%rho2c(irho,iz-1))/(dz**2)
          end do

! Find enpoint values for the derivatives calculated above.
          rhostore(ispecies, jspecies)%rhoz2c(irho,1) = 2.0d0* &
          		rhostore(ispecies, jspecies)%rhoz2c(irho,2) - &
          			rhostore(ispecies, jspecies)%rhoz2c(irho,3)
          rhostore(ispecies, jspecies)%rhoz2c(irho,nnz) = 2.0d0* &
     	     	rhostore(ispecies, jspecies)%rhoz2c(irho,nnz-1)	&
                      - rhostore(ispecies, jspecies)%rhoz2c(irho,nnz-2)

          rhostore(ispecies, jspecies)%rhozz2c(irho,1) = 2.0d0* &
          		rhostore(ispecies, jspecies)%rhozz2c(irho,2) - &
          			rhostore(ispecies, jspecies)%rhozz2c(irho,3)
          rhostore(ispecies, jspecies)%rhozz2c(irho,nnz) = 2.0d0* &
          		rhostore(ispecies, jspecies)%rhozz2c(irho,nnz-1) &
          			- rhostore(ispecies, jspecies)%rhozz2c(irho,nnz-2)
         end do

! **************************************************************************
! Now calculate the cross term derivatives - rhopz2c.
         do irho = 1, nnrho
          do iz = 2, nnz - 1
           rhostore(ispecies, jspecies)%rhopz2c(irho,iz) = &
          	 (rhostore(ispecies, jspecies)%rhop2c(irho,iz+1) - &
          	 	rhostore(ispecies, jspecies)%rhop2c(irho,iz-1))/(2.0d0*dz)
          end do

! Now calculate the derivatives at the endpoints.
          rhostore(ispecies, jspecies)%rhopz2c(irho,1) = 2.0d0* &
          	rhostore(ispecies, jspecies)%rhopz2c(irho,2) - &
          		rhostore(ispecies, jspecies)%rhopz2c(irho,3)
          rhostore(ispecies, jspecies)%rhopz2c(irho,nnz) = 2.0d0* &
          	rhostore(ispecies, jspecies)%rhopz2c(irho,nnz-1) &
              - rhostore(ispecies, jspecies)%rhopz2c(irho,nnz-2)
         end do
        end if

! Format Statements
! ============================================================================
100     format (2x, ' Evaluating rho_2c_store arrays for nZ = ', i3,         &
     &              ' and nZ = ', i3)

        return
        end subroutine rho2c_store




















! ===========================================================================
! vna_ontopL_Harris
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(1)|psi2>.  Thus V(1) is located at
! the left site of one of the orbitals.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine vna_ontopL_Harris
        implicit none

        include "../include/gridsizes.h"

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species
        integer igrid                       !< number of grid points
        integer index_2c, nME2c_max         !< basically the number of non-zero
        integer isorp, ideriv               !< the number of different types
        integer nFdata_cell_2c              !< indexing of interactions

        real dmax                           !< max distance between two centers
        real drr		            !< distance between mesh points
        real d                              !< distance between the two centers
        real rcutoff1, rcutoff2             !< cutoffs for the two centers

        real zmin, zmax
        real rhomin, rhomax

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 25) filename
        character (len = 25) Fdata_location

        logical skip

! Procedure
! ============================================================================
! Assign values to the unrequired variables for this specific interaction.
        ideriv = 999

! We are doing only Harris here, so set isorp = 0
        isorp = 0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
        	pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

            call make_munu (nFdata_cell_2c, ispecies, jspecies)
            nME2c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME2c_max))

            ! Open ouput file for this species pair
            write (filename, '("/vna_ontopL_", i2.2,".",i2.2,".",i2.2,".dat")')&
     &  	       isorp, species(ispecies)%nZ, species(jspecies)%nZ
            Fdata_location = 'coutput'
            inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
            if (skip) cycle
            open (unit = 11, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown')

            ! Open mu, nu, mvalue file and write out values.
            write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")') &
     &             P_vna_ontopL, species(ispecies)%nZ, species(jspecies)%nZ
            Fdata_location = 'coutput'
            open (unit = 12, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown', position = 'append')

            ! write the mapping - stored in mu, nu, and mvalue
            write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%mvalue_2c(index_2c),                   &
     &                    index_2c = 1, nME2c_max)

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoffA_max
            dmax = na(ispecies)%rcutoffA_max + wf(jspecies)%rcutoffA_max
            drr = dmax/float(ndd_vna - 1)
            d = -drr

            ! Set integration limits
	        rhomin = 0.0d0
    	    rhomax = min(rcutoff1, rcutoff2)

! Loop over grid
            write (*,100) species(ispecies)%nZ, species(jspecies)%nZ
            do igrid = 1, ndd_vna
              d = d + drr

              ! Set integration limits
	          zmin = max(-rcutoff1, d - rcutoff2)
    	      zmax = min(rcutoff1, d + rcutoff2)

              call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies, &
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d, nz_vna, nrho_vna,                &
     &                                   rint_vna_ontopL, phifactor, zmin,   &
     &                                   zmax, rhomin, rhomax, pFdata_cell%fofx)
              ! Write out details.
              write (11,*) (pFdata_cell%fofx(index_2c),                      &
	 &                                       index_2c = 1, nME2c_max)
            end do !igrid
            write (11,*)
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Evaluating vna ontopL integrals for nZ = ', i3,        &
     &              ' and nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vna_ontopL_Harris


! ===========================================================================
! rint_vna_ontopL
! ===========================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute
! Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! Program Declaration
! ===========================================================================
	    function rint_vna_ontopL (itype, ispecies, jspecies, isorp, d, rho,  &
     &                            z1, z2, ideriv, index_2c)
	    implicit none

        real rint_vna_ontopL

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real d                              !< distance between the two centers
        real, intent(in) :: rho
        real, intent(in) :: z1, z2

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy
        integer l1, m1, n1               ! quantum  numbers
        integer l2, m2, n2

        real dummy
        real psi1val, psi2val
        real r1, r2
        real vofr

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        idummy = ideriv
        dummy = d

! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)

! Pick up "passed" data
        n1 =  pFdata_cell%N_mu(index_2c)
        l1 =  pFdata_cell%L_mu(index_2c)
        m1 =  pFdata_cell%M_mu(index_2c)

        n2 =  pFdata_cell%N_nu(index_2c)
        l2 =  pFdata_cell%L_nu(index_2c)
        m2 =  pFdata_cell%M_nu(index_2c)

! Set parameters for actual function
        r1 = sqrt(z1**2 + rho**2)
        r2 = sqrt(z2**2 + rho**2)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.
! Find psi1 value at point r1 as based on above and rho called for by
! adaptive_simpson
        psi1val = psiofr (r1, ispecies, n1)

! Find psi2 value at point r2 as based on above and rho  called for by
! adaptive_simpson
        psi2val = psiofr (r2, jspecies, n2)

! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
        psi1val = rescaled_psi (l1, m1, rho, r1, z1, psi1val)
        psi2val = rescaled_psi (l2, m2, rho, r2, z2, psi2val)
        vofr = vnaofr (r1, ispecies, isorp)

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_vna_ontopL = psi1val*vofr*psi2val*rho

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_vna_ontopL


! ===========================================================================
! vna_ontopR_Harris
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(2)|psi2>.  Thus V(2) is located at
! the right site of the orbitals.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine vna_ontopR_Harris
        implicit none

        include "../include/gridsizes.h"

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species
        integer igrid                       !< number of grid points
        integer index_2c, nME2c_max         !< basically the number of non-zero
        integer isorp, ideriv               !< the number of different types
        integer nFdata_cell_2c              !< indexing of interactions

        real dmax                           !< max distance between two centers
        real drr		            !< distance between mesh points
        real d                              !< distance between the two centers
        real rcutoff1, rcutoff2             !< cutoffs for the two centers

        real zmin, zmax
        real rhomin, rhomax

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 25) filename
        character (len = 25) Fdata_location

        logical skip

! Procedure
! ============================================================================
! Assign values to the unrequired variables for this specific interaction.
        ideriv = 999

! We are doing only Harris here, so set isorp = 0
        isorp = 0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
        	pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

            call make_munu (nFdata_cell_2c, ispecies, jspecies)
            nME2c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME2c_max))

            ! Open ouput file for this species pair
            write (filename, '("/vna_ontopR_", i2.2,".",i2.2,".",i2.2,".dat")')&
     &  	       isorp, species(ispecies)%nZ, species(jspecies)%nZ
            Fdata_location = 'coutput'
            inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
            if (skip) cycle
            open (unit = 11, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown')

            ! Open mu, nu, mvalue file and write out values.
            write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")') &
     &             P_vna_ontopR, species(ispecies)%nZ, species(jspecies)%nZ
            Fdata_location = 'coutput'
            open (unit = 12, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown', position = 'append')

            ! write the mapping - stored in mu, nu, and mvalue
            write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%mvalue_2c(index_2c),                   &
     &                    index_2c = 1, nME2c_max)

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoffA_max
            dmax = wf(ispecies)%rcutoffA_max + na(jspecies)%rcutoffA_max
            drr = dmax/float(ndd_vna - 1)
            d = -drr

            rhomin = 0.0d0
            rhomax = min(rcutoff1, rcutoff2)

! Loop over grid
            write (*,100) species(ispecies)%nZ, species(jspecies)%nZ
            do igrid = 1, ndd_vna
              d = d + drr
              ! Set integration limits
	          zmin = max(-rcutoff1, d - rcutoff2)
    	      zmax = min(rcutoff1, d + rcutoff2)
              call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies, &
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d,  nz_vna, nrho_vna,          &
     &                                   rint_vna_ontopR, phifactor, zmin,   &
     &                                   zmax, rhomin, rhomax, pFdata_cell%fofx)
              ! Write out details.
              write (11,*) (pFdata_cell%fofx(index_2c),                      &
	 &                                       index_2c = 1, nME2c_max)
            end do !igrid
            write (11,*)
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Evaluating vna ontopR integrals for nZ = ', i3,        &
     &              ' and nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vna_ontopR_Harris


! ===========================================================================
! rint_vna_ontopL
! ===========================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute
! Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! Program Declaration
! ===========================================================================
        function rint_vna_ontopR (itype, ispecies, jspecies, isorp, d, rho,  &
      &                           z1, z2, ideriv, index_2c)
        implicit none

        real rint_vna_ontopR

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real d                              !< distance between the two centers
        real, intent(in) :: rho
        real, intent(in) :: z1, z2

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy
        integer l1, m1, n1               ! quantum  numbers
        integer l2, m2, n2

        real dummy
        real psi1val, psi2val
        real r1, r2
        real vofr

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        idummy = ideriv
        dummy = d

! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)

!Pick up "passed" data
        n1 =  pFdata_cell%N_mu(index_2c)
        l1 =  pFdata_cell%L_mu(index_2c)
        m1 =  pFdata_cell%M_mu(index_2c)

        n2 =  pFdata_cell%N_nu(index_2c)
        l2 =  pFdata_cell%L_nu(index_2c)
        m2 =  pFdata_cell%M_nu(index_2c)

!Set parameters for actual function
        r1 = sqrt(z1**2 + rho**2)
        r2 = sqrt(z2**2 + rho**2)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.
! Find psi1 value at point r1 as based on above and rho called for by
! adaptive_simpson
        psi1val = psiofr (r1, ispecies, n1)

!find psi2 value at point r2 as based on above and rho  called for by
! adaptive_simpson
        psi2val = psiofr (r2, jspecies, n2)

! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
        psi1val = rescaled_psi (l1, m1, rho, r1, z1, psi1val)
        psi2val = rescaled_psi (l2, m2, rho, r2, z2, psi2val)
        vofr = vnaofr (r2, jspecies, isorp)

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_vna_ontopR = psi1val*vofr*psi2val*rho

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_vna_ontopR


! ===========================================================================
! vna_atom
! ===========================================================================
! SUbroutine Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V|psi2>; the two wavefunctions psi1, and
! psi2 are located on one site and the potential V is located on the other
! site.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
! This subroutine then writes the results to file.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! Subroutine Declaration
! ===========================================================================
        subroutine vna_atom_Harris
        implicit none

        include "../include/gridsizes.h"

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species
        integer igrid                       !< number of grid points
        integer index_2c, nME2c_max         !< basically the number of non-zero
        integer isorp, ideriv               !< the number of different types
        integer nFdata_cell_2c              !< indexing of interactions

        real dmax                           !< max distance between two centers
        real drr		            !< distance between mesh points
        real d                              !< distance between the two centers
        real rcutoff1, rcutoff2             !< cutoffs for the two centers

        real zmin, zmax
        real rhomin, rhomax

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 25) filename
        character (len = 25) Fdata_location

        logical skip

! Procedure
! ============================================================================
! Assign values to the unrequired variables for this specific interaction.
        ideriv = 999

! We are doing only Harris here, so set isorp = 0
        isorp = 0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
        	pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

            call make_munu_atom (nFdata_cell_2c, ispecies, jspecies)
            nME2c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME2c_max))

            ! Open ouput file for this species pair
            write (filename, '("/vna_atom_", i2.2,".",i2.2,".",i2.2,".dat")')&
     &  	       isorp, species(ispecies)%nZ, species(jspecies)%nZ
            Fdata_location = 'coutput'
            inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
            if (skip) cycle
            open (unit = 11, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown')

            ! Open mu, nu, mvalue file and write out values.
            write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")') &
     &             P_vna_atom, species(ispecies)%nZ, species(jspecies)%nZ
            Fdata_location = 'coutput'
            open (unit = 12, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown', position = 'append')

            ! write the mapping - stored in mu, nu, and mvalue
            write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%mvalue_2c(index_2c),                   &
     &                    index_2c = 1, nME2c_max)

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoffA_max
            dmax = wf(ispecies)%rcutoffA_max + na(jspecies)%rcutoffA_max
            drr = dmax/float(ndd_vna - 1)
            d = -drr

	        rhomin = 0.0d0
    	    rhomax = min(rcutoff1, rcutoff2)

! Loop over grid
            write (*,100) species(ispecies)%nZ, species(jspecies)%nZ
            do igrid = 1, ndd_vna
              d = d + drr

              ! Set integration limits
	          zmin = max(-rcutoff1, d - rcutoff2)
    	      zmax = min(rcutoff1, d + rcutoff2)

              call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies, &
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d, nz_vna, nrho_vna, rint_vna_atom, &
     &                                   phifactor, zmin, zmax, rhomin,      &
     &                                   rhomax, pFdata_cell%fofx)
              ! Write out details.
              write (11,*) (pFdata_cell%fofx(index_2c),                      &
	 &                                       index_2c = 1, nME2c_max)
            end do !igrid
            write (11,*)
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None


! Format Statements
! ===========================================================================
100     format (2x, ' Evaluating vna atom integrals for nZ = ', i3,          &
     &              ' and nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vna_atom_Harris


! ===========================================================================
! rint_vna_atom
! ===========================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute
! Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! Program Declaration
! ===========================================================================
        function rint_vna_atom (itype, ispecies, jspecies, isorp, d, rho, z1,&
      &                         z2, ideriv, index_2c)
        implicit none

        real rint_vna_atom

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real d                              !< distance between the two centers
        real, intent(in) :: rho
        real, intent(in) :: z1, z2

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy
        integer l1, m1, n1               ! quantum  numbers
        integer l2, m2, n2

        real dummy
        real psi1val, psi2val
        real r1, r2
        real vofr

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        idummy = ideriv
        dummy = d

! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)

!Pick up "passed" data
        n1 =  pFdata_cell%N_mu(index_2c)
        l1 =  pFdata_cell%L_mu(index_2c)
        m1 =  pFdata_cell%M_mu(index_2c)

        n2 =  pFdata_cell%N_nu(index_2c)
        l2 =  pFdata_cell%L_nu(index_2c)
        m2 =  pFdata_cell%M_nu(index_2c)

!Set parameters for actual function
        r1 = sqrt(z1**2 + rho**2)
        r2 = sqrt(z2**2 + rho**2)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.
! Find psi1 value at point r1 as based on above and rho called for by
! adaptive_simpson
        psi1val = psiofr (r1, ispecies, n1)

!find psi2 value at point r2 as based on above and rho  called for by
! adaptive_simpson
        psi2val = psiofr (r1, ispecies, n2)

! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
        psi1val = rescaled_psi (l1, m1, rho, r1, z1, psi1val)
        psi2val = rescaled_psi (l2, m2, rho, r1, z1, psi2val)
        vofr = vnaofr (r2, jspecies, isorp)

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_vna_atom = psi1val*vofr*psi2val*rho

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_vna_atom

! End Module
! =============================================================================
        end module
