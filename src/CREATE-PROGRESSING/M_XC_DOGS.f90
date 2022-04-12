! copyright info:
!
!                             @Copyright 2007
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
! M_case4.f90
! Program Description
! ============================================================================
!      This is a module calculating the integrals of two centers for the
!
! ============================================================================
! Code written by:
! Hong Wang
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)

! ============================================================================
! Module Declaration
! ============================================================================
	module M_XC_DOGS
    use M_atom_functions
    use M_species
	use M_integrals_2c
	use M_VNL
	use Exchange_Extra



      implicit none
! Type Declaration
! ============================================================================
! two-center interactions arrays

! To cut down on storage space, we actually change the storage procedure
! from previous FIREBALL code. Not all atoms have the same number of
! interactions or interaction types. Before - we would store things based
! on the maximum number of fdata points, maximum number of interactions types,
! maximum number of matrix elements - so even hydrogen-hydrogen (just ss
! and/or ss*) stored a 4x4 or an 8x8 matrix even when not needed.  This was
! quite inefficient.

! The new approach is to define some Fdata types which store the actual
! Fdata points. The smallest unit storage is called Fdata_cell_2C, containing
! all Fdata for a particular interaction/subinteraction type.
! ===========================================================================
       type T_Rhostore
          	real, pointer :: rho2c(:,:)       !The density as a function of species type,
!                        r, z and d.  Output is placed in common block
!                        density located in wavefunctions.inc
            real, pointer ::  rhop2c(:,:)			!		:     derivative with respect to rho
            real, pointer ::  rhopp2c(:,:) 	   !	second derivative with respect to rho
            real, pointer ::  rhoz2c(:,:)  	   !	derivative with respect to z
            real, pointer ::  rhozz2c(:,:) 	   !	second derivative with respect to rho
            real, pointer ::  rhopz2c(:,:) 	   !	mixed derivative with respect to rho and z

        end type T_Rhostore

        type(T_Rhostore), pointer :: Rhostore (:,:)


	integer, allocatable :: iderorb(:)
	real, allocatable :: dqint(:,:), dqorb(:)
	integer ideriv
	integer, parameter :: ixcgridfactor = 15

! module procedures
        contains
! ===========================================================================
! initialize_xc_DOGS
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
!and
!
! Barry Haycock
! Department of Physics,
! Dublin Institute of Technology,
! Dublin 2.
! +353 1 402 7960
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine initialize_xc_DOGS
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies, ideriv          !< counters for number of species

        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
           do ideriv = 2, 5
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 3
            end do !ideriv
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
        end subroutine initialize_xc_DOGS






! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine calls the subroutines required to calculate the vna_ontop l/r and
! atom cases.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine xc_DOGS
        implicit none

        include "../include/gridsizes.h"
! Parameters and Data Declaration
! ===========================================================================
! None
! ===========================================================================
! Argument Declaration and Description
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None
! Procedure
! ===========================================================================
		write (*,*) 'Calling XC Ontop'
		call xc_ontop_DOGS
		write (*,*) 'Returned from calling XC Ontop'

		write (*,*) 'Calling exc atom'
		call xc_atom_DOGS
		write (*,*) 'Returned from calling XC Atom'

		write (*,*) 'Calling exc correction'
		call xc_correction_DOGS
		write (*,*) 'Returned from calling XC Correction'


! Format Statements
! ===========================================================================
! None
! ==========================================================================
        return
        end subroutine xc_DOGS


! Program Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(1)|psi2>.  Thus V(1) is located at
! the site of one of the orbitals.  The potential V(1) is something like Vxc
! for the exchange correlation potential, Vna for the neutral atom potential,
! or 1 for the overlap term.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
!
! This subroutine then writes the results to file.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2,
! Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960

        subroutine xc_ontop_DOGS
        implicit none

        include "../include/gridsizes.h"


        real dmax  	!max distance between two centers
        real drr, d		!distance between points in the mesh

        integer ispecies, jspecies !For loop over species
        integer iexc, issh, isorp, nFdata_cell_2c
        integer igrid, index_2c, nME2c_max


        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle
        character (len=25) filename
        character (len=25) Fdata_location


        real zmin, zmax
		real rcutoff1, rcutoff2
		real rhomin, rhomax

		integer, parameter :: nrho_points = 1601
		integer, parameter :: nz_points = 2501
! Procedure
! ============================================================================
    	isorp = 999
    	!Allocate some memory
        allocate (rhostore(nspecies, nspecies))
    	allocate (dqorb(nspecies))
    	allocate (iderorb(nspecies))
    	allocate (dqint(4,nspecies))		!this '4' here is very messy!!





!Set up dqorb and iderorb.iderorb is the shell where
! the charge will be changed in the +-Q derivative stuff, and dqorb is
! the charge amount changed.
	do ispecies = 1, nspecies
    		iderorb(ispecies) = species(ispecies)%nssh
   ! 	write (*,*) 'ispecies, iderorb', ispecies, iderorb(ispecies)
		dqorb(ispecies) = 0.5d0
		if (species(ispecies)%nssh .eq. 1) dqorb(ispecies) = 0.25d0
		do issh = 1, species(ispecies)%nssh
  			dqint(issh,ispecies) = dqorb(ispecies)/species(ispecies)%nssh
		end do
	end do




        do ispecies = 1, nspecies
        	do jspecies = 1, nspecies

!allocate some arrays
			allocate (rhostore(ispecies, jspecies)%rho2c(nrho_points, nz_points))
			allocate (rhostore(ispecies, jspecies)%rhop2c(nrho_points, nz_points))
    	    allocate (rhostore(ispecies, jspecies)%rhopp2c(nrho_points, nz_points))
       		allocate (rhostore(ispecies, jspecies)%rhoz2c(nrho_points, nz_points))
        	allocate (rhostore(ispecies, jspecies)%rhozz2c(nrho_points, nz_points))
        	allocate (rhostore(ispecies, jspecies)%rhopz2c(nrho_points, nz_points))
        	iexc = PP_species(ispecies)%iexc
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
        	pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

        	write(*,*) 'Calling Make MuNu For XC_Ontop'
	    	call make_munu(nFdata_cell_2c, ispecies, jspecies)
    		write (*,*) 'Back from calling make MuNu'

            nME2c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME2c_max))
!I've made this a little more complicated than need be. For Harris, ideriv = 1, for DOGS iderive is from 2 to 5, this is actually ix in the orignal code + 1

            	do ideriv =2, 5
! Open ouput file for this species pair
           	write (filename, '("/xc_ontop_"i2.2".",i2.2,".",i2.2,".dat")')          &
     		&     ideriv-1, species(ispecies)%nZ, species(jspecies)%nZ
            	Fdata_location = 'coutput'
            	open (unit = 12, file = trim(Fdata_location)//trim(filename),   &
     		&     status = 'unknown')
!Set up grid loop control constants
     		dmax = species(ispecies)%rcutoffA_max + &
     		& species(jspecies)%rcutoffA_max
        	drr = dmax/ float(ndd_xc-1)
        	d = -drr
!Loop over grid

		do igrid = 1, ndd_xc
        	d = d + drr
			call rho2c_store(ispecies, jspecies, iexc, d) !, ideriv)

 			rcutoff1 = species(ispecies)%rcutoffA_max
			rcutoff2 = species(jspecies)%rcutoffA_max

              ! Set integration limits
	          zmin = max(-rcutoff1, d - rcutoff2)
    	      zmax = min(rcutoff1, d + rcutoff2)

	          rhomin = 0.0d0
    	      rhomax = min(rcutoff1, rcutoff2)

			  call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies, &
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d, dmax, nz_xc, nrho_xc,          &
     &                                   rint_xc_ontop, phifactor,zmin,       &
     &                                   zmax, rhomin, rhomax, pFdata_cell%fofx)

!Write out details.
	        write (12,*) (pFdata_cell%fofx(index_2c), index_2c = 1, nME2c_max)

		end do!igrid
		write (12,*) '   ' !blank line

		end do !jpecies
		end do !ispecies
		end do !ideriv
!de!allocate some arrays

			deallocate (rhostore)
    		deallocate (dqorb)
    		deallocate (iderorb)
    		deallocate (dqint)
      return
 		! Format Statements
! ===========================================================================
2222    format(I4,'  ', f10.7,'   ', f12.4)

! End Subroutine
! =============================================================================
        end subroutine xc_ontop_DOGS



! rint
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2,
! Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! Program Declaration
! ===========================================================================

	function  rint_xc_ontop(itype, ispecies, jspecies, isorp, d, rho, z1, z2, ideriv, index_2c)

! Argument Declaration
! ===========================================================================
	real, intent(in) :: z1, z2, d
	integer, intent(in) :: ispecies, jspecies, isorp, itype, ideriv, index_2c

	integer n1, n2, m1, m2, l1, l2

	real, intent(in) :: rho
	real r1, r2
	real rint_xc_ontop
	real psi1val, psi2val
	real vofr
	real fraction
	integer iexc
	type (T_Fdata_cell_2c), pointer :: pFdata_cell
    type (T_Fdata_bundle_2c), pointer :: pFdata_bundle


	pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
    pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)

! Procedure
! ===========================================================================
!Pick up "passed" data
	n1 =  pFdata_cell%N_mu(index_2c)
    n2 =  pFdata_cell%N_nu(index_2c)
    m1 =  pFdata_cell%M_mu(index_2c)
	m2 =  pFdata_cell%M_nu(index_2c)
	l1 =  pFdata_cell%L_mu(index_2c)
	l2 =  pFdata_cell%L_nu(index_2c)

	!d = z1 - z2

 	!Set parameters for actual function
	r1 = sqrt(z1**2 + rho**2)
	r2 = sqrt(z2**2 + rho**2)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.


! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
!find psi1 value at point r1 as based on above and rho called for by adaptive_simpson
	psi1val = psiofr (r1, ispecies, n1)
!find psi2 value at point r2 as based on above and rho  called for by adaptive_simpson
	psi2val = psiofr (r2, jspecies, n2)
	psi1val = rescaled_psi (l1, m1, rho, r1, z1, psi1val)
	psi2val = rescaled_psi (l2, m2, rho, r2, z2, psi2val)

!In the atom case, V is on the second atom, jspecies.

        fraction = PP_species(ispecies)%fraction
        iexc = PP_species(ispecies)%iexc
        vofr = vxc(rho, z1, iexc, fraction, ispecies, jspecies, d)
!Actual function (Ylm's are calculated above and multilplied after integration
	rint_xc_ontop=psi1val*vofr*psi2val*rho
! Format Statements
! ===========================================================================
!
! End Function
! =============================================================================
	end function rint_xc_ontop

! ============================================================================
! make_munu
! ============================================================================
! Subroutine Description
! ============================================================================
!       This subroutine calculates the following information (for all pairs
! of atoms (in1,in2)):
!
! num_orb (in1): number of orbitals in atom-type in1
! mu (index,in1,in2): the mu-position for each matrix-element (index) between
!                     atom-type in1 and atom-type in2
! nu (index,in1,in2): the nu-position for each matrix-element (index) between
!                     atom-type in1 and atom-type in2

! (on the BOX ( num_orb(in1) x num_orb(in2)))
!
! Atoms 1 and 2 (bondcharge) are along the z-axis; the third atom is in the
! xz-plane. The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            0
!
!   P-shell :           py   pz   px
!                       -1   0    +1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1   0    +1     +2
!
!
! For the "atom" case, the number of non-zero matrix elements
! is dependent only on one atom. Both wavefunctions are located at
! the same site. This does that by "acknowledging" the jspecies call,
! But comopletely ignoring it, except for addressing the N_Nu L_Nu M_Nu's
!
! Subroutine Declaration
! ============================================================================
        subroutine make_munuatom(itype, ispecies, jspecies)
        implicit none

        include 'constants.h'
        include '../include/gridsizes.h'

! Auguments Declaration and Description
! None

! Parameters and Data Declaration
! ============================================================================
! None

! Input
! ============================================================================

! Local Variable Declaration adn Description
! ============================================================================
        integer, intent(in) :: itype

        integer index_2c                ! counter for matrix location - mu, nu
        integer ispecies
        integer jspecies                ! index for looping over the species
        integer issh
        integer jssh                    ! index for looping over shells
        integer mvalue

        integer nME2c_max

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len=25) filename
        character (len=25) Fdata_location
        integer l1, m1, n1, l2, m2, n2

! Allocate Arrays
! ============================================================================
! None

! Procedure
! ============================================================================
! Loop over the pairs of species.  For each species pair, establish what the
! quantum number values for the orbital mu (the left orbital) and nu (the
! right orbital).
            ! cut some lengthy notation
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            nME2c_max = 0

! First, find the maximum number of matrix elements for each species pair.
            do issh = 1, species(ispecies)%nssh
              l1 = species(ispecies)%shell(issh)%lssh
              do jssh = 1, species(ispecies)%nssh
                l2 = species(ispecies)%shell(jssh)%lssh
                do mvalue = -min(l1,l2), min(l1,l2)
                  nME2c_max = nME2c_max + 1
                end do
              end do
            end do

! Now allocate the sizes for mu_2c, nu_2c, and the quantum numbers NLM
! for each mu and nu pair.
! Overlap Interactions:
            pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)
            pFdata_cell%nME = nME2c_max

            allocate (pFdata_cell%mu_2c(nME2c_max))
            allocate (pFdata_cell%nu_2c(nME2c_max))
            allocate (pFdata_cell%mvalue_2c(nME2c_max))

            allocate (pFdata_cell%N_mu(nME2c_max))
            allocate (pFdata_cell%L_mu(nME2c_max))
            allocate (pFdata_cell%M_mu(nME2c_max))

            allocate (pFdata_cell%N_nu(nME2c_max))
            allocate (pFdata_cell%L_nu(nME2c_max))
            allocate (pFdata_cell%M_nu(nME2c_max))

! Set the values for NLM of each mu, nu pair.
            index_2c = 0
            n1 = 0
            do issh = 1, species(ispecies)%nssh
              l1 = species(ispecies)%shell(issh)%lssh
              n1 = n1 + l1 + 1
            n2 = 0
            do jssh = 1, species(ispecies)%nssh
               l2 = species(ispecies)%shell(jssh)%lssh
               n2 = n2 + l2 + 1

               do mvalue = -min(l1,l2), min(l1,l2)

                 index_2c = index_2c + 1
                 pFdata_cell%mu_2c(index_2c) = n1 + mvalue
                 pFdata_cell%nu_2c(index_2c) = n2 + mvalue
                 pFdata_cell%mvalue_2c(index_2c) = 0

                 pFdata_cell%N_mu(index_2c) = issh
                 pFdata_cell%L_mu(index_2c) = l1
                 pFdata_cell%M_mu(index_2c) = mvalue

                 pFdata_cell%N_nu(index_2c) = jssh
                 pFdata_cell%L_nu(index_2c) = l2
                 pFdata_cell%M_nu(index_2c) = mvalue

                end do
                n2 = n2 + l2
             end do
             n1 = n1 + l1
           end do

            ! Open mu, nu, mvalue file and write out values.
            write (filename, '("/",i2.2,"munu_2c.",i2.2,".",i2.2,".dat")')          &
     &        itype, species(ispecies)%nZ, species(jspecies)%nZ
            Fdata_location = 'coutput'
            open (unit = 12, file = trim(Fdata_location)//trim(filename),   &
     &            status = 'unknown', position = 'append')

            ! read the mapping - stored in mu, nu, and mvalue
            write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%mvalue_2c(index_2c),                  &
     &                    index_2c = 1, nME2c_max)


! End Subroutine
! =============================================================================
        return
        end subroutine make_munuAtom

! XC_atom
! Program Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(1)|psi2>.  Thus V(1) is located at
! the site of one of the orbitals.  The potential V(1) is something like Vxc
! for the exchange correlation potential, Vna for the neutral atom potential,
! or 1 for the overlap term.
! This subroutine then writes the results to file.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2,
! Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960

        subroutine xc_atom_DOGS
        implicit none

        include "../include/gridsizes.h"


        real dmax  	!max distance between two centers
        real drr		!distance between points in the mesh

        integer ispecies, jspecies !For loop over species
        integer iexc, issh, isorp
        integer igrid, index_2c, nME2c_max
        integer nmu, lmu, mmu, nnu, lnu, mnu
		integer l1, m1, n1, l2, m2, n2
		integer index_start, nFdata_cell_2c


        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle
        character (len=25) filename
        character (len=25) Fdata_location
        !Simpsons Quadrature Variables.
		real temp
		real dz, d
		integer nnz
		real zmult (1001)
		integer iz
		integer nz
		real zmin, zmax
		real rhomin, rhomax
		real phifact
		real sum

        integer, parameter :: nrho_points = 1601
		integer, parameter :: nz_points = 2501

		real z1
		real rcutoff1, rcutoff2
! Procedure
! ============================================================================
    	isorp = 999
    !Allocate some memory
        allocate (rhostore(nspecies, nspecies))
    	allocate (dqorb(nspecies))
    	allocate (iderorb(nspecies))
    	allocate (dqint(4,nspecies))		!this '4' here is very messy!!



!Set up dqorb and iderorb.iderorb is the shell where
! the charge will be changed in the +-Q derivative stuff, and dqorb is
! the charge amount changed.
	do ispecies = 1, nspecies
    		iderorb(ispecies) = species(ispecies)%nssh
		dqorb(ispecies) = 0.5d0
		if (species(ispecies)%nssh .eq. 1) dqorb(ispecies) = 0.25d0
		do issh = 1, species(ispecies)%nssh
  			dqint(issh,ispecies) = dqorb(ispecies)/species(ispecies)%nssh
		end do
	end do




        do ispecies = 1, nspecies
        	do jspecies = 1, nspecies
!allocate some arrays
			allocate (rhostore(ispecies, jspecies)%rho2c(nrho_points, nz_points))
			allocate (rhostore(ispecies, jspecies)%rhop2c(nrho_points, nz_points))
    	    allocate (rhostore(ispecies, jspecies)%rhopp2c(nrho_points, nz_points))
       		allocate (rhostore(ispecies, jspecies)%rhoz2c(nrho_points, nz_points))
        	allocate (rhostore(ispecies, jspecies)%rhozz2c(nrho_points, nz_points))
        	allocate (rhostore(ispecies, jspecies)%rhopz2c(nrho_points, nz_points))
			iexc = PP_species(ispecies)%iexc
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
        	pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

    		write(*,*) 'Calling Make MuNu for XC_Atom'

   		 	call make_munuAtom(nFdata_cell_2c, ispecies, jspecies)
    		write (*,*) 'Back from calling make MuNu'

           	nME2c_max = pFdata_cell%nME
!I've made this a little more complicated than need be. For Harris, ideriv = 1, for DOGS iderive is from 2 to 5, this is actually ix in the orignal code + 1
       	do ideriv = 2, 5
! Open ouput file for this species pair
       		write (filename, '("/xc_atom_"i2.2".",i2.2,".",i2.2,".dat")')          &
     		&  (ideriv-1), species(ispecies)%nZ, species(jspecies)%nZ
            	Fdata_location = 'coutput'

            	open (unit = 12, file = trim(Fdata_location)//trim(filename),   &
     		&            status = 'unknown')

!Set up grid loop control constants
     		dmax = species(ispecies)%rcutoffA_max + &
     					& species(jspecies)%rcutoffA_max
        	drr = dmax/ float(ndd_xc-1)
        	d = -drr
			allocate (pFdata_cell%fofx(nME2c_max))
			index_start = 1
!-----------------------------------------------------------------------------------------------
! Do the d = 0 case (1C Case)
!-----------------------------------------------------------------------------------------------
		    if (ispecies == jspecies) then
			d = d + drr

			index_start = 2
			do index_2c = 1, nME2c_max
        		!nz = 427
        		nz = 106
        		n1 = pFdata_cell%N_mu(index_2c)
             	l1 = pFdata_cell%L_mu(index_2c)
                m1 = pFdata_cell%M_mu(index_2c)

                n2 = pFdata_cell%N_nu(index_2c)
                l2 = pFdata_cell%L_nu(index_2c)
                m2 = pFdata_cell%M_nu(index_2c)



 				rcutoff1 = species(ispecies)%rcutoffA_max
				rcutoff2 = species(jspecies)%rcutoffA_max



! Initialize the sum to zero
        sum = 0.0d0

! Set integration limits
        zmin = max(-rcutoff1, -rcutoff2)
        zmax = min(rcutoff1, rcutoff2)

        rhomin = 0.0d0
        rhomax = min(rcutoff1, rcutoff2)
! This factor (phifactor) comes from a result of multiplying the two Ylm's together.
! and after the factor of pi is multiplied out after the phi integration.
! For the coulomb case, we need spherically symmetric charge densities,
! so add up all m's squared

        phifact = phifactor(l1,m1,l2,m2)
! Integration is over z (z-axis points from atom 1 to atom 2) and rho (rho is
! radial distance from z-axis).

!-----------------------------------------------------------------------------------
! Non Adaptive Simpson's Setup
!-----------------------------------------------------------------------------------
! Strictly define what the density of the mesh should be.  Make the density of
! the number of points equivalent for all cases. Change the number of points
! to be integrated to be dependent upon the distance between the centers and
! this defined density.
        dz = (rcutoff1 + rcutoff2)/(dfloat(nz)*2.0d0)
        nnz = int((zmax - zmin)/dz)
        if (mod(nnz,2) .eq. 0) nnz = nnz + 1
! Set up Simpson's rule factors. First for the z integration and then for
! the rho integration.
        zmult(1) = dz/3.0d0
        zmult(nnz) = dz/3.0d0
        do iz = 2, nnz - 1, 2
         zmult(iz) = 4.0d0*dz/3.0d0
        end do
        do iz = 3, nnz - 2, 2
         zmult(iz) = 2.0d0*dz/3.0d0
        end do
!-----------------------------------------------------------------------------------
! Non Adaptive Simpson's Setup
!-----------------------------------------------------------------------------------

			call rho2c_store(ispecies, jspecies, iexc, d)
			sum = 0.00
			do iz = 1, nnz
       		  z1 = zmin + dfloat(iz-1)*dz
			  temp = zint_xc_atom1c(z1, rcutoff1, rcutoff2, ispecies, jspecies, l1,m1, n1, m2, l2, n2, isorp)
			  temp = temp * zmult(iz)
			  sum = sum + temp
			end do

		sum = phifact*sum

!Store 1C data for write-out to file with 2C data
        pFdata_cell%fofx(index_2c) = sum

		end do !index
      	write (12,*) (pFdata_cell%fofx(index_2c), index_2c = 1, nME2c_max)
 !     	write (12,*) '   '
 		d = 0.00
 		end if
!-----------------------------------------------------------------------------------------------
! End of 1c- Onwards to 2C
!-----------------------------------------------------------------------------------------------

!Loop over grid
		do igrid = index_start, ndd_xc
        	d = d + drr
			allocate (pFdata_cell%fofx(nME2c_max))
			!write (*,*) 'Calling rho2c_store'
			call rho2c_store(ispecies, jspecies, iexc, d)
			!write(*,*) 'Back from calling rho2c_store'
			rcutoff1 = species(ispecies)%rcutoffA_max
			rcutoff2 = species(jspecies)%rcutoffA_max

              ! Set integration limits
	          zmin = max(-rcutoff1, d - rcutoff2)
    	      zmax = min(rcutoff1, d + rcutoff2)
	          rhomin = 0.0d0
    	      rhomax = min(rcutoff1, rcutoff2)

			call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies, &
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d, dmax, nz_xc, nrho_xc,          &
     &                                   rint_xc_atom, phifactor, zmin,       &
     &                                   zmax, rhomin, rhomax, pFdata_cell%fofx)

     	!Write out details.
	        	write (12,*) (pFdata_cell%fofx(index_2c), index_2c = 1, nME2c_max)

		end do!igrid
		write (12,*) '   ' !blank line


	end do !jpecies
	end do !ispecies
	end do !ideriv

!de!allocate some arrays
			deallocate (rhostore)
    		deallocate (dqorb)
    		deallocate (iderorb)
    		deallocate (dqint)
     return
! Format Statements
! ===========================================================================
2222    format(I4,'  ', f10.7,'   ', f12.4)

! End Subroutine
! =============================================================================
        end subroutine xc_atom_DOGS


! rint
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2,
! Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! Program Declaration
! ===========================================================================

	function  rint_xc_atom(itype, ispecies, jspecies, isorp, d, rho, z1, z2, ideriv, index_2c)
! Argument Declaration
! ===========================================================================
	real, intent(in) :: rho
	real, intent(in) :: z1, z2, d
	integer, intent(in) :: ispecies, jspecies, isorp, itype, ideriv, index_2c
	real r1, r2

	real rint_xc_atom
	real psi1val, psi2val
	real vofr
	integer n1, n2, m1, m2, l1, l2

	real fraction
	integer iexc

	type (T_Fdata_cell_2c), pointer :: pFdata_cell
    type (T_Fdata_bundle_2c), pointer :: pFdata_bundle


	pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
    pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)
! Procedure
! ===========================================================================
!Pick up "passed" data
	n1 =  pFdata_cell%N_mu(index_2c)
    n2 =  pFdata_cell%N_nu(index_2c)
    m1 =  pFdata_cell%M_mu(index_2c)
	m2 =  pFdata_cell%M_nu(index_2c)
	l1 =  pFdata_cell%L_mu(index_2c)
	l2 =  pFdata_cell%L_nu(index_2c)

	!d = z1 - z2
!Set parameters for actual function
	r1 = sqrt(z1**2 + rho**2)
	r2 = sqrt(z2**2 + rho**2)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.


! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
!find psi1 value at point r1 as based on above and rho called for by adaptive_simpson

	psi1val = psiofr (r1, ispecies, n1)
!find psi2 value at point r2 as based on above and rho  called for by adaptive_simpson
	psi2val = psiofr (r1, ispecies, n2)

	psi1val = rescaled_psi (l1, m1, rho, r1, z1, psi1val)
	psi2val = rescaled_psi (l2, m2, rho, r1, z1, psi2val)

!In the atom case, V is on the second atom, jspecies.

        fraction = PP_species(ispecies)%fraction
        iexc = PP_species(ispecies)%iexc
        vofr = dvxc (ispecies, jspecies, rho, z1, r1, iexc, fraction, ideriv, d)
!Actual function (Ylm's are calculated above and multilplied after integration
		rint_xc_atom=psi1val*vofr*psi2val*rho
! Format Statements
! ===========================================================================
!
! End Function
! =============================================================================
	end function rint_xc_atom



! rint
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2,
! Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! Program Declaration
! ===========================================================================

	function  rint_xc_atom1c(rho,ispecies, jspecies, z1, z2, l1,m1, n1, m2, l2, n2, isorp)
! Argument Declaration
! ===========================================================================
	integer, intent (in) :: ispecies, jspecies
	integer, intent (in) :: l1,m1, n1, m2, l2, n2, isorp
	real, intent(in) :: rho, z1, z2
	real r1, r2, d
	real rint_xc_atom1c
	real psi1val, psi2val
	real vofr
	real fraction
	integer iexc

! Procedure
! ===========================================================================
	d = z1 - z2

!Set parameters for actual function
	r1 = sqrt(z1**2 + rho**2)
	r2 = sqrt(z2**2 + rho**2)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.


! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
!find psi1 value at point r1 as based on above and rho called for by adaptive_simpson

	psi1val = psiofr (r1, ispecies, n1)
!find psi2 value at point r2 as based on above and rho  called for by adaptive_simpson
	psi2val = psiofr (r1, ispecies, n2)

	psi1val = rescaled_psi (l1, m1, rho, r1, z1, psi1val)
	psi2val = rescaled_psi (l2, m2, rho, r1, z1, psi2val)

!In the atom case, V is on the second atom, jspecies.

        fraction = PP_species(ispecies)%fraction
        iexc = PP_species(ispecies)%iexc
        vofr = dvxc (ispecies, jspecies, rho, z1, r1, iexc, fraction, ideriv, d)
!Actual function (Ylm's are calculated above and multilplied after integration
		rint_xc_atom1c = psi1val*vofr*psi2val*rho

! Format Statements
! ===========================================================================
!
! End Function
! =============================================================================
	end function rint_xc_atom1c


! zint
! Program Description
! ===========================================================================
! This program contains, for the overlap integral, the actual function.
! It is written ins such a way as the function is called by a single variable
! only, this means that the Adaptive_Simpson subroutine can integrate it.
! All other required variables are global.
! This is the function in z, which contains the function in rho (as explained
! in twocenter_overlap) adaptive simpson is called within zint to integrate over
! rho (rint) in situ.
!
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2,
! Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! Program Declaration
! ===========================================================================
	function zint_xc_atom1c(z, rcutoff1, rcutoff2, ispecies, jspecies, l1,m1, n1, m2, l2, n2, isorp)

! Argument Declaration
! ===========================================================================
	integer, intent (in) :: ispecies, jspecies
	integer, intent (in) :: l1,m1, n1, m2, l2, n2, isorp
	real, intent(in) :: z, rcutoff1, rcutoff2
	real :: integral
	real zint_xc_atom1c

	real rhomin, rhomax
	real z1, z2
!Simpsons Variables.
	real temprho
	real drho
	integer nnrho
	real rhomult (1001)
	integer irho
	real rho
	integer nrho

! Procedure
! ===========================================================================
!Set parameters for integration
	nrho = 106
	!nrho = 427
	z1 = z
	z2 = z1
	rhomin = 0.0d0
	rhomax = min(rcutoff1, rcutoff2)


	integral = 0.00

!-----------------------------------------------------------------------------------
! Non Adaptive Simpson's Setup
!-----------------------------------------------------------------------------------
! Strictly define what the density of the mesh should be.  Make the density of
! the number of points equivalent for all cases. Change the number of points
! to be integrated to be dependent upon the distance between the centers and
! this defined density.
        drho = max(rcutoff1,rcutoff2)/dfloat(nrho)
        nnrho = int((rhomax - rhomin)/drho)
        if (mod(nnrho,2) .eq. 0) nnrho = nnrho + 1

        rhomult(1) = drho/3.0d0
        rhomult(nnrho) = drho/3.0d0
        do irho = 2, nnrho - 1, 2
         rhomult(irho) = 4.0d0*drho/3.0d0
        end do
        do irho = 3, nnrho - 2, 2
         rhomult(irho) = 2.0d0*drho/3.0d0
        end do
!-----------------------------------------------------------------------------------
! Non Adaptive Simpson's Setup
!-----------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
! Actually use the normal Simpson. Set uprho integration
!-----------------------------------------------------------------------------------
			integral = 0.00

			do irho = 1, nnrho
          	  rho = rhomin + dfloat(irho-1)*drho
			  temprho = rint_xc_atom1c(rho,ispecies, jspecies, z1, z2, l1,m1, n1, m2, l2, n2, isorp)
			  temprho = temprho * rhomult(irho)
			  integral = integral + temprho
			end do


	zint_xc_atom1c = integral

! Format Statements
! ===========================================================================
!
! End Function
! =============================================================================
	end function zint_xc_atom1c

! vxc.f
! Program Description
! ===========================================================================
!
! This function computes vxc(n1+n2) with the densities ni of atom
! in_i at the distance r_i from their centers.
!
! On input:   r, z: geometry information for the charge gradient
!
! On output:  vxc: vxc[n1(r1) + n2(r2)]

! We calculate vxc(n1+n2). The catch comes in when we compute
! derivatives. We compute neutral, neutral for ideriv1. For other ideriv's we
! have the following KEY:
!
! (xy) means charge on (1,2). Case 1 (KEY=1),
! neutral neutral corresponds to (00) etc.
! KEY = 1,2,3,4,5 for ideriv=1,2,3,4,5
!
!                             dq Atom 2 axis
!                                   |
!                                   + (0+) KEY=5
!                                   |
!                                   |
!                                   |
!                       KEY=2       |  KEY=1
!                      (-0)         |(00)        (+0) KEY=3
!                     --+-----------O-----------+--
!                                   |           dq Atom 1 axis
!                                   |
!                                   |
!                                   |
!                                   |
!                                   +(0-) KEY=4
!                                   |
!
! ======================================================================
! Original code from Juergen Fritsch

! Code rewritten by:
! James P. Lewis and Richard B. Evans
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ======================================================================
!
! Program Declaration
! ======================================================================
        real function vxc (r, z, iexc, fraction1, ispecies, jspecies, d)
        !use Exchange_Extra
        implicit none



! Argument Declaration and Description
! ======================================================================
! Input
        integer iexc

        real fraction1
        real fraction
        real r
        real z
        integer ispecies, jspecies
        real d

! Local Parameters and Data Declaration
! ======================================================================
!        real abohr
!        parameter (abohr = 0.529177249d0)

! Local Variable Declaration and Description
! ======================================================================
        real dens
        real densp
        real denspp
        real densz
        real denszz
        real denspz
        real dnuxc2c
        real dnuxcs2c
        real exc2c
!        real hartree
        real rin
        real vxc2c

        real rcutoff1, rcutoff2
        real zmin, zmax
        real rhomin, rhomax

        real, pointer :: rho2c(:,:)
        real, pointer :: rhop2c(:,:)
        real, pointer :: rhopp2c(:,:)
        real, pointer :: rhoz2c(:,:)
        real, pointer :: rhozz2c(:,:)
        real, pointer :: rhopz2c(:,:)

        real dzraw, dz, drho
        integer nnz, nnrho
        integer, parameter :: nrho_points = 1601
		integer, parameter :: nz_points = 2501

! Procedure
! ======================================================================
	fraction = fraction1

	rcutoff1 = species(ispecies)%rcutoffA_max
	rcutoff2 = species(jspecies)%rcutoffA_max

        zmin = min(-rcutoff1, d - rcutoff2)
        zmax = max(rcutoff1, d + rcutoff2)

        rhomin = 0.0d0
        rhomax = max(rcutoff1, rcutoff2)



        dens = 0.00
        densp = 0.00
        denspp = 0.00
        densz = 0.00
        denszz = 0.00
        denspz = 0.00
        dnuxc2c = 0.00
        dnuxcs2c = 0.00





        dzraw = min(wf(ispecies)%dr_min,wf(jspecies)%dr_min)	!drr_rho(itype2))

        dz = dzraw*ixcgridfactor

! If the grid size is made too large or too small, then set to defaults.
        if (dz .gt. 0.05d0) dz = 0.05d0
        if (dz .lt. 0.002d0) dz = 0.002d0

        nnz = int((zmax - zmin)/dz) + 1

        drho = dz
        nnrho = int((rhomax - rhomin)/drho) + 1

!		All the extra functions in interpolate 2c
!cut down on some notations.
		rho2c =>  rhostore(ispecies, jspecies)%rho2c
        rhop2c => rhostore(ispecies, jspecies)%rhop2c
        rhopp2c => rhostore(ispecies, jspecies)%rhopp2c
        rhoz2c => rhostore(ispecies, jspecies)%rhoz2c
        rhozz2c => rhostore(ispecies, jspecies)%rhozz2c
        rhopz2c => rhostore(ispecies, jspecies)%rhopz2c



! Interpolate the density and gradients of the density at the given
! point (r, z).
        call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                           nnrho, nrho_points, nnz, nz_points, rho2c,  &
                           dens)
! Only interpolate the derivatives if doing GGA exchange-correlation.
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 &
           .or. iexc .eq. 9 .or. iexc .eq. 10) then
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, rhop2c, &
                            densp)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, &
                            rhopp2c, denspp)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, &
                            rhopz2c, denspz)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points,  &
                            rhoz2c, densz)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, &
                            rhozz2c, denszz)
        end if
! Convert to atomic units
        rin = r/abohr
        dens = dens*abohr**3
        densp = densp*abohr**4
        densz = densz*abohr**4
        denspp = denspp*abohr**5
        denszz = denszz*abohr**5
        denspz = denspz*abohr**5

! Here energy and potential due to exchange and correlation are calculated.
        call get_potxc2c (iexc, fraction, rin, dens, densp, denspp, &
                         densz, denszz, denspz, exc2c, vxc2c, dnuxc2c, &
                         dnuxcs2c)

! Answers are in Hartrees convert to eV.
!        hartree = 14.39975d0/abohr
        vxc = hartree*vxc2c

! Format Statements
! ===========================================================================

        return
        end function vxc

! dvxc.f
! Program Description
! ===========================================================================
!
! This function computes vxc(n1+n2) - vxc(n1) with the densities ni of atom
! in_i at the distance r_i from their centers.
!
! On input:  in1,in2: atomic indices
!            r, z: geometry information for the charge gradient
!            ix: switch for the derivatives
!
! On output:  vxc: vxc[n1(r1) + n2(r2)] - vxc[n1(r1)]

! We calculate vxc(n1+n2) - vxc(n1). The catch comes in when we compute
! derivatives. We compute neutral, neutral for ideriv1. For other ideriv's we
! have the following KEY:
!
! (xy) means charge on (1,2). Case 1 (KEY=1),
! neutral neutral corresponds to (00) etc.
! KEY = 1,2,3,4,5 for ideriv=1,2,3,4,5
!
!                             dq Atom 2 axis
!                                   |
!                                   + (0+) KEY=5
!                                   |
!                                   |
!                                   |
!                       KEY=2       |  KEY=1
!                      (-0)         |(00)        (+0) KEY=3
!                     --+-----------O-----------+--
!                                   |           dq Atom 1 axis
!                                   |
!                                   |
!                                   |
!                                   |
!                                   +(0-) KEY=4
!                                   |
!
! ===========================================================================
! Original code from Juergen Fritsch

! Code rewritten by:
! James P. Lewis and Richard B. Evans
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        real function dvxc (in1, in2, r, z, r1, iexc, fraction1, ix, d)
        implicit none


!        include '../wavefunctions.inc'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc
        integer in1, in2
        integer ix

        real fraction1
        real fraction
        real r
        real r1
        real z, d

! Local Parameters and Data Declaration
! ===========================================================================
 !       real abohr
 !       parameter (abohr = 0.529177249d0)

! Local Variable Declaration and Description
! ===========================================================================
        integer in3

        real rcutoff1, rcutoff2
        real zmin, zmax
        real rhomin, rhomax

        real, pointer :: rho2c(:,:)
        real, pointer :: rhop2c(:,:)
        real, pointer :: rhopp2c(:,:)
        real, pointer :: rhoz2c(:,:)
        real, pointer :: rhozz2c(:,:)
        real, pointer :: rhopz2c(:,:)

        real dens
        real densp
        real denspp
        real densz
        real denszz
        real denspz
        real dexc1c
        real dexc2c
        real dnuxc
        real dnuxcs
        real dnuxc2c
        real dnuxcs2c
        real dvxc1c
        real dvxc2c
!        real hartree
        real rin

        real dzraw, dz, drho
        integer nnz, nnrho
		integer, parameter :: nrho_points = 1601
		integer, parameter :: nz_points = 2501



! Procedure
! ===========================================================================
	fraction = fraction1
	rcutoff1 = species(in1)%rcutoffA_max
	rcutoff2 = species(in2)%rcutoffA_max

        zmin = min(-rcutoff1, d - rcutoff2)
        zmax = max(rcutoff1, d + rcutoff2)

        rhomin = 0.0d0
        rhomax = max(rcutoff1, rcutoff2)



        dens = 0.00
        densp = 0.00
        denspp = 0.00
        densz = 0.00
        denszz = 0.00
        denspz = 0.00
        dnuxc2c = 0.00
        dnuxcs2c = 0.00




!		All the extra functions in interpolate 2c
!cut down on some notations.
		rho2c =>  rhostore(in1, in2)%rho2c
        rhop2c => rhostore(in1, in2)%rhop2c
        rhopp2c => rhostore(in1, in2)%rhopp2c
        rhoz2c => rhostore(in1, in2)%rhoz2c
        rhozz2c => rhostore(in1, in2)%rhozz2c
        rhopz2c => rhostore(in1, in2)%rhopz2c



        dzraw = min(wf(in1)%dr_min,wf(in2)%dr_min)

        dz = dzraw*ixcgridfactor

! If the grid size is made too large or too small, then set to defaults.
        if (dz .gt. 0.05d0) dz = 0.05d0
        if (dz .lt. 0.002d0) dz = 0.002d0

        nnz = int((zmax - zmin)/dz) + 1

        drho = dz
        nnrho = int((rhomax - rhomin)/drho) + 1

! By default in3 = in2
        in3 = in2

! Two-center piece: vxc[n1 + n2(r,z)]
! ***************************************************************************
! Interpolate the density and gradients of the density at the given
! point (r, z).
        call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                           nnrho, nrho_points, nnz, nz_points, rho2c, &
                           dens)


! Only interpolate the derivatives if doing GGA exchange-correlation.
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 &
           .or. iexc .eq. 9 .or. iexc .eq. 10) then
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, rhop2c, &
                            densp)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, &
                            rhopp2c, denspp)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, &
                            rhopz2c, denspz)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, &
                            rhoz2c, densz)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, &
                            rhozz2c, denszz)
        end if

! Convert to atomic units
        rin = r/abohr
        dens = dens*abohr**3
        densp = densp*abohr**4
        densz = densz*abohr**4
        denspp = denspp*abohr**5
        denspz = denspz*abohr**5
        denszz = denszz*abohr**5

! Here energy and potential due to exchange and correlation are calculated.

        call get_potxc2c (iexc, fraction, rin, dens, densp, denspp, &
                         densz, denszz, denspz, dexc2c, dvxc2c, &
                         dnuxc2c, dnuxcs2c)

! One-center piece: vxc[n1(r1)]
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in1
        call density_calc (iexc, ix, 1, in1, in2, in3, r1, drho, &
                          dens, densp, denspp)

        rin = r1/abohr
        dens = dens*abohr**3
        densp = densp*abohr**4
        denspp = denspp*abohr**5


        call get_potxc1c (iexc, fraction, rin, dens, densp, denspp, &
                         dexc1c, dvxc1c, dnuxc, dnuxcs)

! Answers are in Hartrees convert to eV.
 !       hartree = 14.39975d0/abohr
        dvxc = hartree*(dvxc2c - dvxc1c)

! Format Statements
! ===========================================================================

        return
        end function
! rho2c_store.f
! Program Description
! ====================================================================
!
! This routine calculates and stores the combined density of species
! 1 and 2 as a function of r, z and d.
!
! On input:
!            iexc:       The exchange-correlation option used
!            itype1:     The species type of atom 1
!            itype2:     The species type of atom 2
!            rcutoff1:   The radius cutoff of atom 1
!            rcutoff2:   The radius cutoff of atom 2
!
!     By common blocks:
!            nsshxc:     Number of shells for each species
!
! On output: rho2c:      The density as a function of species type,
!                        r, z and d.  Output is placed in common block
!                        density located in wavefunctions.inc
!            rhop2c:     derivative with respect to rho
!            rhopp2c:    second derivative with respect to rho
!            rhoz2c:     derivative with respect to z
!            rhozz2c:    second derivative with respect to rho
!            rhopz2c:    mixed derivative with respect to rho and z
!
! ====================================================================
! Code written by:
! Richard B. Evans
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ====================================================================
!
! Program Declaration
! ====================================================================
	subroutine rho2c_store (ispecies, jspecies, iexc, d)
   	implicit none




! Input Declaration and Description
! ====================================================================
! Input
  	integer iexc
    	integer ispecies
       	integer jspecies, aspecies
    	integer ndd, nz, nrho, nnz, nnrho

    	real d, dz, drho
    	real rcutoff1
    	real rcutoff2
    	real rhomin, rhomax

! Local Parameters and Data Declaration
! ====================================================================
        integer j1(5)
        data j1 /1, 1, 1, 2, 2/

        integer j2(5)
        data j2 /0, -1, 1, -1, 1/

! Local Variable Declaration and Description
! ====================================================================

        integer irho
        integer issh
        integer iz
        integer j1at
        integer j1ch
        integer jssh, kssh

        integer in(2)

        real dens
        real dzraw
        real r
        real r1
        real r2
		real rho
        real xinv4pi
        real z1
        real z2

        real zmin, zmax
		integer, parameter :: nrho_points = 1601
		integer, parameter :: nz_points = 2501

! Procedure
! ====================================================================
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


	ndd = 107
	nz = 106
	nrho = 106
	rcutoff1 = species(ispecies)%rcutoffA_max
	rcutoff2 = species(jspecies)%rcutoffA_max



! First initialize in (index) array. This is for the charge correction
! calculated later.
        in(1) = ispecies
        in(2) = jspecies

        j1ch = j1(ideriv)
        j1at = in(j1ch)
        jssh = iderorb(j1at)

!Fix the endpoints and initialize the increments dz and drho.
        zmin = min(-rcutoff1, d - rcutoff2)
        zmax = max(rcutoff1, d + rcutoff2)

        rhomin = 0.0d0
        rhomax = max(rcutoff1, rcutoff2)

        dzraw = min(wf(ispecies)%dr_min,wf(jspecies)%dr_min)	!drr_rho(itype2))
        dz = dzraw*ixcgridfactor

! If the grid size is made too large or too small, then set to defaults.
        if (dz .gt. 0.05d0) dz = 0.05d0
        if (dz .lt. 0.002d0) dz = 0.002d0

        nnz = int((zmax - zmin)/dz) + 1

        drho = dz
        nnrho = int((rhomax - rhomin)/drho) + 1

! Performs some checks to make sure that the dimensions on rho2c are
! alright.
        if (nnrho .gt. nrho_points) then
         write (*,*) ' In rho2c_store.f, nnrho = ', nnrho
         write (*,*) ' The dimension nrho_points =', nrho_points
         write (*,*) ' needs to be increased. '
         stop 'error in rho2c_store'
        end if

        if (nnz .gt. nz_points) then
         write (*,*) ' In rho2c_store.f, nnz = ', nnz
         write (*,*) ' The dimension nz_points =', nz_points
         write (*,*) ' needs to be increased. '
         stop 'error in rho2c_store'
        end if





! **************************************************************************
! Here we loop over z and r computing the sum of the densities
! for species in1 and in2 at each value of d, z and r.

        do iz = 1, nnz
         z1 = zmin + (iz-1)*dz
         z2 = z1 - d

         do irho = 1, nnrho
          rho = rhomin + (irho-1)*drho
          r1 = sqrt(z1**2 + rho**2)
          r2 = sqrt(z2**2 + rho**2)

          dens = 0.0d0
          do issh = 1, species(ispecies)%nssh
           dens = dens + species(ispecies)%shell(issh)%xnocc* &
           						psiofr(r1, ispecies, issh)**2
          end do
          do kssh = 1, species(jspecies)%nssh
           dens = dens + species(jspecies)%shell(kssh)%xnocc* &
           						psiofr(r2, jspecies,kssh)**2
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
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6	&
           .or. iexc .eq. 9 .or. iexc .eq. 10) then

! **************************************************************************
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
! **************************************************************************
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
! ====================================================================

        return
        end subroutine rho2c_store


! xc_correction
! Program Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(1)|psi2>.  Thus V(1) is located at
! the site of one of the orbitals.  The potential V(1) is something like Vxc
! for the exchange correlation potential, Vna for the neutral atom potential,
! or 1 for the overlap term.
! This subroutine then writes the results to file.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2,
! Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960

        subroutine xc_correction_DOGS
        implicit none

		include "../include/gridsizes.h"

        real dmax  	!max distance between two centers
        real drr		!distance between points in the mesh

        integer ispecies, jspecies !For loop over species
        integer iexc, issh, isorp, nFdata_cell_2c
        integer igrid, index_2c, nME2c_max



        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle
        character (len=25) filename
        character (len=25) Fdata_location

		real d
		real zmin, zmax
		real rhomin, rhomax

		real rcutoff1, rcutoff2
        integer, parameter :: nrho_points = 1601
		integer, parameter :: nz_points = 2501

! Procedure
! ============================================================================
		isorp = 999
!Allocate some memory
        allocate (rhostore(nspecies, nspecies))
    	allocate (dqorb(nspecies))
    	allocate (iderorb(nspecies))
    	allocate (dqint(4,nspecies))		!this '4' here is very messy!!

!Set up dqorb and iderorb.iderorb is the shell where
! the charge will be changed in the +-Q derivative stuff, and dqorb is
! the charge amount changed.
	do ispecies = 1, nspecies
    		iderorb(ispecies) = species(ispecies)%nssh
    	    	dqorb(ispecies) = 0.5d0
		if (species(ispecies)%nssh .eq. 1) dqorb(ispecies) = 0.25d0
		do issh = 1, species(ispecies)%nssh
  		  	dqint(issh,ispecies) = dqorb(ispecies)/species(ispecies)%nssh
		end do
	end do


    	write (*,*) 'No Call to MuNU for XC_CORRECTION'
    	write (*,*) 'This is not a matrix element, but an'
    	write (*,*) 'over-counting correction to the'
    	write (*,*) 'exchange-correlation interaction,'
    	write (*,*) 'so set index_max = 1'


! This is not a matrix element, but an over-counting correction to the
! exchange-correlation interaction, so set index_max = 1




  	do ispecies = 1, nspecies
  		do jspecies = 1, nspecies
!allocate some arrays
				allocate (rhostore(ispecies, jspecies)%rho2c(nrho_points, nz_points))
				allocate (rhostore(ispecies, jspecies)%rhop2c(nrho_points, nz_points))
    	    	allocate (rhostore(ispecies, jspecies)%rhopp2c(nrho_points, nz_points))
       			allocate (rhostore(ispecies, jspecies)%rhoz2c(nrho_points, nz_points))
        		allocate (rhostore(ispecies, jspecies)%rhozz2c(nrho_points, nz_points))
        		allocate (rhostore(ispecies, jspecies)%rhopz2c(nrho_points, nz_points))

				iexc = PP_species(ispecies)%iexc
	            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
    	        pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
        	    nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
        		pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)


        		pFdata_cell%nME = 1
            	nME2c_max = pFdata_cell%nME

            		allocate(pFdata_cell%N_mu(nME2c_max))
            		allocate(pFdata_cell%N_nu(nME2c_max))
            		allocate(pFdata_cell%M_mu(nME2c_max))
	        		allocate(pFdata_cell%M_nu(nME2c_max))
	        		allocate(pFdata_cell%L_mu(nME2c_max))
	        		allocate(pFdata_cell%L_nu(nME2c_max))

            		pFdata_cell%N_mu(nME2c_max) = 1
            		pFdata_cell%N_nu(nME2c_max) = 1
            		pFdata_cell%M_mu(nME2c_max) = 1
	        		pFdata_cell%M_nu(nME2c_max) = 1
	        		pFdata_cell%L_mu(nME2c_max) = 1
	        		pFdata_cell%L_nu(nME2c_max) = 1

            		allocate (pFdata_cell%fofx(nME2c_max))

!I've made this a little more complicated than need be. For Harris, ideriv = 1, for DOGS iderive is from 2 to 5, this is actually 'ix' in the orignal code + 1
            	do ideriv =2, 5
! Open ouput file for this species pair
 		        write (filename, '("/xc_corr_"i2.2".",i2.2,".",i2.2,".dat")')          &
     			&         (ideriv-1), species(ispecies)%nZ, species(jspecies)%nZ
            		Fdata_location = 'coutput'

            		open (unit = 12, file = trim(Fdata_location)//trim(filename),   &
    			 &            status = 'unknown')

!Set up grid loop control constants
     			dmax = species(ispecies)%rcutoffA_max + &
     						& species(jspecies)%rcutoffA_max
        		drr = dmax/ float(ndd_xc-1)
        		d = -drr

!Loop over grid

		do igrid = 1, ndd_xc
        		d = d + drr
			allocate (pFdata_cell%fofx(nME2c_max))

			call rho2c_store(ispecies, jspecies, iexc, d)
	 		rcutoff1 = species(ispecies)%rcutoffA_max
			rcutoff2 = species(jspecies)%rcutoffA_max

		! Set integration limits
	          zmin = max(-rcutoff1, d - rcutoff2)
    	      zmax = min(rcutoff1, d + rcutoff2)
	          rhomin = 0.0d0
    	      rhomax = min(rcutoff1, rcutoff2)

			call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies, &
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d, dmax, nz_xc, nrho_xc,          &
     &                                   rint_xc_correction, twopi,zmin,       &
     &                                   zmax, rhomin, rhomax, pFdata_cell%fofx)


!Write out details.
				write (12,*) (pFdata_cell%fofx(index_2c), index_2c = 1, nME2c_max)

		end do!igrid
			write (12,*) '   ' !blank line

	end do !ideriv
	end do !jpecies
	end do !ispecies

!de!allocate some arrays
			deallocate (rhostore)
    		deallocate (dqorb)
    		deallocate (iderorb)
    		deallocate (dqint)

! Format Statements
! ===========================================================================
2222    format(I4,'  ', f10.7,'   ', f12.4)

! End Subroutine
! =============================================================================
        end subroutine xc_correction_DOGS





! rint
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2,
! Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! Program Declaration
! ===========================================================================

	function  rint_xc_correction(itype, ispecies, jspecies, isorp, d, rho, z1, z2, ideriv, index_2c)

! Argument Declaration
! ===========================================================================
	real, intent(in) :: rho
	real, intent(in) :: z1, z2, d
	integer, intent(in) :: ispecies, jspecies, isorp, itype, ideriv, index_2c
	real r1, r2

	real rint_xc_correction

	real psi1val, psi2val
	real vofr
	integer n1, n2, m1, m2, l1, l2
	real fraction
	integer iexc

	!real d

	type (T_Fdata_cell_2c), pointer :: pFdata_cell
    type (T_Fdata_bundle_2c), pointer :: pFdata_bundle


	pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
    pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)

! Procedure
! ===========================================================================
!Pick up "passed" data
	n1 =  pFdata_cell%N_mu(index_2c)
    n2 =  pFdata_cell%N_nu(index_2c)
    m1 =  pFdata_cell%M_mu(index_2c)
	m2 =  pFdata_cell%M_nu(index_2c)
	l1 =  pFdata_cell%L_mu(index_2c)
	l2 =  pFdata_cell%L_nu(index_2c)

	!d = z1 - z2

 	!Set parameters for actual function
	r1 = sqrt(z1**2 + rho**2)
	r2 = sqrt(z2**2 + rho**2)

        fraction = PP_species(ispecies)%fraction
        iexc = PP_species(ispecies)%iexc				!??? also chech bohr/angs thing, and dnuxc2c
        vofr = dexc (ispecies, jspecies, rho, z1, r1, r2, iexc, fraction, ideriv, d)
		rint_xc_correction=vofr*rho
! Format Statements
! ===========================================================================
!
! End Function
! =============================================================================
	end function rint_xc_correction


! dexc.f
! Program Description
! ===========================================================================
!
! This function computes
!
!    (n1+n2)*(exc(n1+n2) - vxc(n1+n2)) - sum_i ni*(exc(ni) - vxc(n1))
!
! with the densities ni of atom in_i at the distance r_i from their centers.
!
! On input:  in1,in2: atomic indices
!            r, z: geometry information for the charge gradient
!            ix: switch for the derivatives
!
! On output:  dexc

! The catch comes in when we compute derivatives. We compute neutral,
! neutral for ideriv1. For other ideriv's we have the following KEY:
!
! (xy) means charge on (1,2). Case 1 (KEY=1),
! neutral neutral corresponds to (00) etc.
! KEY = 1,2,3,4,5 for ideriv=1,2,3,4,5
!
!                             dq Atom 2 axis
!                                   |
!                                   + (0+) KEY=5
!                                   |
!                                   |
!                                   |
!                       KEY=2       |  KEY=1
!                      (-0)         |(00)        (+0) KEY=3
!                     --+-----------O-----------+--
!                                   |           dq Atom 1 axis
!                                   |
!                                   |
!                                   |
!                                   |
!                                   +(0-) KEY=4
!                                   |
!
! ===========================================================================
! Original code from Juergen Fritsch

! Code rewritten by:
! James P. Lewis and Richard B. Evans
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        real function dexc (in1, in2, r, z, r1, r2, iexc, fraction, ix, d)
        use Exchange_Extra
        implicit none


 !       include '../wavefunctions.inc'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc
        integer in1, in2
        integer ix

        real fraction
        real fraction1
        real r
        real r1
        real r2
        real z, d

! Local Parameters and Data Declaration
! ===========================================================================
!        real abohr
!        parameter (abohr = 0.529177249d0)

! Local Variable Declaration and Description
! ===========================================================================
        integer in3

        real dens
        real densin
        real densp
        real denspin
        real denspp
        real densppin
        real densz
        real denszin
        real denszz
        real denszzin
        real denspz
        real denspzin
        real dexc1c
        real dexc2c
        real dnuxc
        real dnuxcs
        real dnuxc2c
        real dnuxcs2c
        real dvxc1c
        real dvxc2c
!        real hartree
        real rin

        real rcutoff1, rcutoff2
        real zmin, zmax
        real rhomin, rhomax

        real, pointer :: rho2c(:,:)
        real, pointer :: rhop2c(:,:)
        real, pointer :: rhopp2c(:,:)
        real, pointer :: rhoz2c(:,:)
        real, pointer :: rhozz2c(:,:)
        real, pointer :: rhopz2c(:,:)

        real dzraw, dz, drho
        integer nnz, nnrho
		integer, parameter :: nrho_points = 1601
		integer, parameter :: nz_points = 2501

! Procedure
! ===========================================================================

		fraction1 = fraction

		rcutoff1 = species(in1)%rcutoffA_max
		rcutoff2 = species(in2)%rcutoffA_max

        zmin = min(-rcutoff1, d - rcutoff2)
        zmax = max(rcutoff1, d + rcutoff2)

        rhomin = 0.0d0
        rhomax = max(rcutoff1, rcutoff2)



        dens = 0.00
        densp = 0.00
        denspp = 0.00
        densz = 0.00
        denszz = 0.00
        denspz = 0.00
        dnuxc2c = 0.00
        dnuxcs2c = 0.00




!		All the extra functions in interpolate 2c
!cut down on some notations.
		rho2c =>  rhostore(in1, in2)%rho2c
        rhop2c => rhostore(in1, in2)%rhop2c
        rhopp2c => rhostore(in1, in2)%rhopp2c
        rhoz2c => rhostore(in1, in2)%rhoz2c
        rhozz2c => rhostore(in1, in2)%rhozz2c
        rhopz2c => rhostore(in1, in2)%rhopz2c



        dzraw = min(wf(in1)%dr_min,wf(in2)%dr_min)	!drr_rho(itype2))

        dz = dzraw*ixcgridfactor

! If the grid size is made too large or too small, then set to defaults.
        if (dz .gt. 0.05d0) dz = 0.05d0
        if (dz .lt. 0.002d0) dz = 0.002d0

        nnz = int((zmax - zmin)/dz) + 1

        drho = dz
        nnrho = int((rhomax - rhomin)/drho) + 1
! By default set in3 = in2
        in3 = in2

! Two-center piece: [n1 + n2(r,z)]*(exc[n1 + n2(r,z)] - vxc[n1 + n2(r,z)])
! ***************************************************************************
! Interpolate the density and gradients of the density at the given
! point (r, z).
        call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                           nnrho, nrho_points, nnz, nz_points, rho2c, &
                           dens)

! Only interpolate the derivatives if doing GGA exchange-correlation.
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 &
           .or. iexc .eq. 9 .or. iexc .eq. 10) then
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, rhop2c, &
                            densp)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, &
                            rhopp2c, denspp)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, &
                            rhopz2c, denspz)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, &
                            rhoz2c, densz)
         call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                            nnrho, nrho_points, nnz, nz_points, &
                            rhozz2c, denszz)
        end if

! Convert to atomic units
        rin = r/abohr
        densin = dens*abohr**3
        denspin = densp*abohr**4
        denszin = densz*abohr**4
        densppin = denspp*abohr**5
        denspzin = denspz*abohr**5
        denszzin = denszz*abohr**5

! Here energy and potential due to exchange and correlation are calculated.
        call get_potxc2c (iexc, fraction1, rin, densin, denspin, &
                         densppin, denszin, denszzin, denspzin, dexc2c, &
                         dvxc2c, dnuxc2c, dnuxcs2c)

! Answers are in Hartrees convert to eV.
!        hartree = 14.39975d0/abohr
        dexc = hartree*dens*(dexc2c - dvxc2c)

! One-center piece for atom 1: n1*(exc[n1(r1)] - vxc[n1(r1)])
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in1
        call density_calc (iexc, ix, 1, in1, in2, in3, r1, drho, &
                          dens, densp, denspp)

        rin = r1/abohr
        densin = dens*abohr**3
        denspin = densp*abohr**4
        densppin = denspp*abohr**5
        call get_potxc1c (iexc, fraction1, rin, densin, denspin, &
                         densppin, dexc1c, dvxc1c, dnuxc, dnuxcs)

! Answers are in Hartrees convert to eV.
        dexc = dexc - hartree*dens*(dexc1c - dvxc1c)

! One-center piece for atom 2: n2*(exc[n2(r1)] - vxc[n2(r1)])
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in2
        call density_calc (iexc, ix, 2, in1, in2, in3, r2, drho, &
                          dens, densp, denspp)

        rin = r2/abohr
        densin = dens*abohr**3
        denspin = densp*abohr**4
        densppin = denspp*abohr**5
        call get_potxc1c (iexc, fraction1, rin, densin, denspin, &
                         densppin, dexc1c, dvxc1c, dnuxc, dnuxcs)

! Answers are in Hartrees convert to eV.
        dexc = dexc - hartree*dens*(dexc1c - dvxc1c)

! Format Statements
! ===========================================================================

        return
        end function





! density_calc.f
! Program Description
! ======================================================================
!
!  The subroutine density_calc will calculate the density and its
!  first and second derivatives for a given radius, and atom type
!
! ======================================================================
! Code written by:
! Richard B. Evans
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! Salt Lake City, UT 84112-0850
! (801) 581-8606 (office)      email: rbevans@hec.utah.edu
! (801) 581-4353 (fax)
!
! ======================================================================
        subroutine density_calc (iexc, ix, ispec, itype1, itype2,	&
                                itype3, r, dr, dens, ddens, dddens)
        implicit none



! Argument Declaration and Description
! ======================================================================
! Input
        integer iexc
        integer ispec
        integer itype1
        integer itype2
        integer itype3
        integer ix

        real r
        real dr

! Output
        real dens
        real ddens
        real dddens

! Local Parameters and Data Declaration
! ===========================================================================
        integer j1(7)
        data j1 /1, 1, 1, 2, 2, 3, 3/

        integer j2(7)
        data j2 /0, -1, 1, -1, 1, -1, 1/

! Local Variable Delcaration and Descrition
! ======================================================================
        integer issh
        integer itype
        integer j1ch
        integer j1at
        integer jssh

        integer in (3)

        real densdr
        real dens2dr
        real dens_dr
 !       real psiofr
        real xinv4pi

!        external psiofr

! Procedure
! ===========================================================================
! Initialize 1/4pi
        xinv4pi = 1.0d0/(4.0d0*3.141592653589793238462643D0)

        in(1) = itype1
        in(2) = itype2
        in(3) = itype3
        itype = in(ispec)

! Here the density is computed for r, r+dr and r-dr
        dens = 0.0d0
        do issh = 1, species(itype)%nssh
         dens = dens + species(itype)%shell(issh)%xnocc*psiofr(r,itype,issh)**2
        end do

! Here the density with respect to the charge correction term is calculated
! and added to dens, densdr and dens2dr.  The variable switch determines
! whether the correction is for the one, two or three center case.
        j1ch = j1(ix)
        if (j1ch .eq. ispec) then
         j1at = in(j1ch)
         jssh = iderorb(j1at)
         dens = dens + j2(ix)*dqorb(j1at)*psiofr(r,j1at,jssh)**2
        end if


! Only calculate the derivatives if doing GGA exchange-correlation.
! ***************************************************************************
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6	&
           .or. iexc .eq. 9 .or. iexc .eq. 10) then

         densdr = 0.0d0
         dens_dr = 0.0d0
         do issh = 1, species(itype)%nssh
          densdr =	&
          densdr + species(itype)%shell(issh)%xnocc*psiofr(r+dr,itype,issh)**2
          dens_dr =	&
          dens_dr + species(itype)%shell(issh)%xnocc*psiofr(r-dr,itype,issh)**2
         end do

         if (j1ch .eq. ispec) then
          densdr =	&
          densdr + j2(ix)*dqorb(j1at)*psiofr(r+dr,j1at,jssh)**2
          dens_dr =	&
          dens_dr + j2(ix)*dqorb(j1at)*psiofr(r+dr,j1at,jssh)**2
         end if

! Here the first and second derivatives of the density is computed.
         if ((r - dr) .gt. 1.0d-5) then
          ddens = (densdr - dens_dr)/(2.0d0*dr)
          dddens = (densdr - 2.0d0*dens + dens_dr)/dr**2
         else

! At the endpoint do a forward difference. First, we need the point at r+2dr.
          dens2dr = 0.0d0
          do issh = 1, species(itype)%nssh
           dens2dr =	&
           dens2dr + species(itype)%shell(issh)%xnocc*psiofr(r+2.0*dr,itype,issh)**2
          end do

          if (j1ch .eq. ispec) then
           dens2dr = dens2dr	&
           + j2(ix)*dqorb(j1at)*psiofr(r+2.0*dr,j1at,jssh)**2
          end if

          ddens = (densdr - dens)/dr
          dddens = (dens2dr - 2.0d0*densdr + dens)/dr**2
         end if
        end if

! Convert to the correct units
        dens = dens*xinv4pi
        ddens = ddens*xinv4pi
        dddens = dddens*xinv4pi

! Format Statements
! ===========================================================================

        return
        end subroutine







! interpolate2d.f
! Program Description
! ===========================================================================
!       This routine is a two-dimensional interpolater on a 4x4 sub-grid.
!
! ===========================================================================
! Code rewritten by:
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine interpolate2d (rhoin, rhomin, rhomax, drho, zin, &
                                 zmin, zmax, dz, nnrho, nrho_points, &
                                 nnz, nz_points, frho, answer)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer nnrho
        integer nnz
        integer nrho_points
        integer nz_points

        real drho
        real dz
        real rhoin
        real rhomax
        real rhomin
        real zin
        real zmax
        real zmin

! This is the function being interpolated
        real frho (nrho_points, nz_points)

! Output
        real answer

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer imidrho
        integer imidz
        integer j
        integer jj
        integer k

        real e36t
        real f0p3
        real f0p6
        real f1m2
        real f1m3
        real f1p3
        real f1p6
        real ftp
        real gradrho
        real gradtest
        real gradz
        real prho
        real prod
        real pz
        real tp

        real b (0:5)
        real bb (0:5, -2:3)
        real g (-2:3)
        real fun (-1:2, -1:2)

! Procedure
! ===========================================================================
! Check and make sure that the point (r, z) is within the limits of the
! stored density.

        if (rhoin .gt. (rhomax + 1.0d-5) .or. rhoin .lt. rhomin) then
         write (*,*) ' What the heck is going on in interpolate2d.f !'
         write (*,*) ' ************* error !!! ************* '
         write (*,*) ' Input rhoin = ', rhoin
         write (*,*) ' Max. data = ', rhomax, ' Min. data = ', rhomin
        end if

        if (zin .gt. (zmax + 1.0d-5) .or. zin .lt. (zmin - 1.0d-5)) then
         write (*,*) ' What the heck is going on in interpolate2d.f !'
         write (*,*) ' ************* error !!! ************* '
         write (*,*) ' Input zin = ', zin
         write (*,*) ' Max. data = ', zmax, ' Min. data = ', zmin
        end if

! Much of the code in this routine is a dodge to avoid interpolating, even
! though the mission of this subprogram is to interpolate. The dodge can be
! approached at two levels: first, if the gradient is locally very small, we
! can just return. This is CRUCIAL for multirc problems, but seems to help for
! monorc.

! Set-up everything for interpolation.
        imidrho = int((rhoin - rhomin)/drho) + 1
        imidz = int((zin - zmin)/dz) + 1

        if (imidrho .lt. 2) imidrho = 2
        if (imidrho .gt. nnrho) imidrho = nnrho
        if (imidz .lt. 2) imidz = 2
        if (imidz .gt. nnz) imidz = nnz

        prho = (rhoin - rhomin)/drho - dfloat(imidrho - 1)
        pz = (zin - zmin)/dz - dfloat(imidz - 1)

        fun(-1,-1) = frho(imidrho - 1,imidz - 1)
        fun(-1, 0) = frho(imidrho - 1,imidz)
        fun(-1, 1) = frho(imidrho - 1,imidz + 1)
        fun(-1, 2) = frho(imidrho - 1,imidz + 2)

        fun(0,-1) = frho(imidrho,imidz - 1)
        fun(0, 0) = frho(imidrho,imidz)
        fun(0, 1) = frho(imidrho,imidz + 1)
        fun(0, 2) = frho(imidrho,imidz + 2)

        fun(1,-1) = frho(imidrho + 1,imidz - 1)
        fun(1, 0) = frho(imidrho + 1,imidz)
        fun(1, 1) = frho(imidrho + 1,imidz + 1)
        fun(1, 2) = frho(imidrho + 1,imidz + 2)

        fun(2,-1) = frho(imidrho + 2,imidz - 1)
        fun(2, 0) = frho(imidrho + 2,imidz)
        fun(2, 1) = frho(imidrho + 2,imidz + 1)
        fun(2, 2) = frho(imidrho + 2,imidz + 2)

! **************************************************************************
! If the gradient is small, then do quadratic quick bivariate interpolation.
        gradrho = (fun(0,0) - fun(1,0))/drho
        gradz = (fun(0,0) - fun(0,1))/dz

        gradtest = abs(gradrho) + abs(gradz)

! Form the criterion for a quick interpolation. Empirically, I find that
! gradtest < 1.0d-05 is adequate. If you dont want this option, change
! 1.0d-05 in the next line to 0.0d0!
        if (gradtest .lt. 1.0d-05) then
         answer = (1.0d0 - prho - pz)*fun(0,0) + prho*fun(1,0) &
                                              + pz*fun(0,1)
         return
        end if

! **************************************************************************
! Phase III. All else fails. Interpolate carefully. Original pfed
! interpolator with minimal multiplies.
        e36t = 1.0d0/36.0d0

        do k = - 1, 2
         f1m2 = fun(k,-1) + fun(k,-1)
         f1m3 = f1m2 + fun(k,-1)

         f0p3 = fun(k,0) + fun(k,0) + fun(k,0)
         f0p6 = f0p3 + f0p3

         f1p3 = fun(k,1) + fun(k,1) + fun(k,1)
         f1p6 = f1p3 + f1p3

         bb(3,k) = - fun(k,-1) + f0p3 - f1p3 + fun(k,2)
         bb(2,k) = f1m3 - f0p6 + f1p3
         bb(1,k) = - f1m2 - f0p3 + f1p6 - fun(k,2)

         tp = fun(k,0)
         tp = tp + tp + tp
         bb(0,k) = tp + tp

         prod = bb(3,k)*pz
         do j = 1, 2
          jj = 3 - j
          ftp = bb(jj,k)
          prod = (prod + ftp)*pz
         end do
         g(k) = prod + bb(0,k)
        end do

        f1m2 = g(-1) + g(-1)
        f1m3 = f1m2 + g(-1)

        f0p3 = g(0) + g(0) + g(0)
        f0p6 = f0p3 + f0p3

        f1p3 = g(1) + g(1) + g(1)
        f1p6 = f1p3 + f1p3

        b(3) = -g(-1) + f0p3 - f1p3 + g(2)
        b(2) = f1m3 - f0p6 + f1p3
        b(1) = -f1m2 - f0p3 + f1p6 - g(2)
        tp = g(0) + g(0) + g(0)
        b(0) = tp + tp

        prod = b(3)*prho
        do j = 1, 2
         jj = 3 - j
         ftp = b(jj)
         prod = (prod + ftp)*prho
        end do
        prod = prod + b(0)

! Final answer
        answer = e36t*prod

! Format Statements
! ===========================================================================

        return
        end subroutine interpolate2d


end module
