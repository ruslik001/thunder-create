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
       module M_ext_hubbard_DOGS
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
!           						             r, z and d.  Output is placed in common block
!           						             density located in wavefunctions.inc
            real, pointer ::  rhop2c(:,:)			!		:     derivative with respect to rho
            real, pointer ::  rhopp2c(:,:) 	   !	second derivative with respect to rho
            real, pointer ::  rhoz2c(:,:)  	   !	derivative with respect to z
            real, pointer ::  rhozz2c(:,:) 	   !	second derivative with respect to rho
            real, pointer ::  rhopz2c(:,:) 	   !	mixed derivative with respect to rho and z

        end type T_Rhostore

        type(T_Rhostore), pointer :: Rhostore (:,:)



		integer, allocatable :: iderorb(:)
		real, allocatable :: dqint(:,:), dqorb(:)
		integer, parameter :: ixcgridfactor = 15
! module procedures
        contains
! ===========================================================================
! initialize_eh_DOGS
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
        subroutine initialize_eh_DOGS
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
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            end do
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
        end subroutine initialize_eh_DOGS
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
		subroutine rho2c_store (ispecies, jspecies, iexc, d, ideriv)! rcutoff1,
                               !rcutoff2, d, ix)      !ix is ideriv +1, ideriv is actually for dogs, a loop from 1-4, for harris, is zero
   		implicit none





! Input Declaration and Description
! ====================================================================
! Input
	  		integer iexc
	   	 	integer ideriv
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
        	integer, parameter :: nrho_points = 1601
			integer, parameter :: nz_points = 2501


        	real zmin, zmax

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

    	    nnz = nint((zmax - zmin)/dz) + 1

    	    drho = dz
    	    nnrho = nint((rhomax - rhomin)/drho) + 1

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
! Subroutine Declaration
! ============================================================================
        subroutine make_munu14(interaction, ispecies, jspecies)
        implicit none

        include 'constants.h'
        include 'gridsizes.h'

! Auguments Declaration and Description
! None

! Parameters and Data Declaration
! ============================================================================
! None

! Input
! ============================================================================

! Local Variable Declaration adn Description
! ============================================================================
        integer index_2c                ! counter for matrix location - mu, nu
        integer ispecies
        integer jspecies                ! index for looping over the species
        integer issh
        integer jssh                    ! index for looping over shells
        integer mvalue

        integer l1, l2
        integer n1, n2
        integer nME2c_max

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len=25) filename
        character (len=25) Fdata_location

		integer, intent(in) ::  interaction

! Allocate Arrays
! ============================================================================
! None

! Procedure
! ============================================================================

! cut some lengthy notation
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            nME2c_max = 0

! First, find the maximum number of matrix elements for each species pair.
            do issh = 1, species(ispecies)%nssh
              l1 = species(ispecies)%shell(issh)%lssh

              do jssh = 1, species(jspecies)%nssh
                l2 = species(jspecies)%shell(jssh)%lssh
                  nME2c_max = nME2c_max + 1
              end do
            end do

! Now allocate the sizes for mu_2c, nu_2c, and the quantum numbers NLM
! for each mu and nu pair.
! Overlap Interactions:
            pFdata_cell=>pFdata_bundle%Fdata_cell_2c(interaction)
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
            do jssh = 1, species(jspecies)%nssh
               l2 = species(jspecies)%shell(jssh)%lssh
               n2 = n2 + l2 + 1
				 mvalue = 0
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
                n2 = n2 + l2
             end do
             n1 = n1 + l1
           end do

            ! Open mu, nu, mvalue file and write out values.
            write (filename, '("/",i2.2,"munu_2c.",i2.2,".",i2.2,".dat")')          &
     &       interaction, species(ispecies)%nZ, species(jspecies)%nZ

            Fdata_location = 'coutput'
!            read the mapping - stored in mu, nu, and mvalue
            open (unit = 13, file = trim(Fdata_location)//trim(filename),   &
     &            status = 'unknown', position = 'append')

            write (13,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
            write (13,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
            write (13,*) (pFdata_cell%mvalue_2c(index_2c),                  &
     &                    index_2c = 1, nME2c_max)

! End Subroutine
! =============================================================================
        return
        end subroutine make_munu14

! na_atom
! Program Description
! ===========================================================================
! This subroutine simply loops over the species and the grid and calls
! twocenter_overlap to do the grunt work, its a driver, for want of a better
! term, the bookeeper.
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

        subroutine eh_DOGS
        implicit none

		include '../include/gridsizes.h'

        real dmax  	!max distance between two centers
        real drr, d		!distance between points in the mesh

        integer ispecies, jspecies !For loop over species
        integer ndd, iexc, issh
        integer igrid, index_2c, nME2c_max
        integer ideriv, isorp
        real rcutoff1, rcutoff2



        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle
        character (len=25) filename
        character (len=25) Fdata_location

        real zmin, zmax
        real rhomin, rhomax

		integer nFdata_cell_2c

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

            		call make_munu (nFdata_cell_2c, ispecies, jspecies)
            		nME2c_max = pFdata_cell%nME
            		allocate (pFdata_cell%fofx(nME2c_max))


				do ideriv = 2, 5
! Open ouput file for this species pair
		           write (filename, '("/xc_eh_",i2.2,".",i2.2,".",i2.2,".dat")')		&
     &  		        ideriv-1, species(ispecies)%nZ, species(jspecies)%nZ
            		Fdata_location = 'coutput'

            		open (unit = 12, file = trim(Fdata_location)//trim(filename),   &
     &          		  status = 'unknown')

!Set up grid loop control constants
     			dmax = species(ispecies)%rcutoffA_max + &
     								& species(jspecies)%rcutoffA_max
     		   drr = dmax/ float(ndd_eh-1)
        		d = -drr

!Loop over grid
				do igrid = 1, ndd_eh
        			d = d + drr

					call rho2c_store(ispecies, jspecies, iexc, d, ideriv)
					rcutoff1 = species(ispecies)%rcutoffA_max
					rcutoff2 = species(jspecies)%rcutoffA_max
              ! Set integration limits
	          zmin = max(-rcutoff1, d - rcutoff2)
    	      zmax = min(rcutoff1, d + rcutoff2)

	          rhomin = 0.0d0
    	      rhomax = min(rcutoff1, rcutoff2)

				call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies, &
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d, dmax, nz_eh, nrho_eh,  &
     &                                   rint_eh, twopi,zmin,       &
     &                                   zmax, rhomin, rhomax, pFdata_cell%fofx)

!Write out details.
			        write (12,*) (pFdata_cell%fofx(index_2c), index_2c = 1, nME2c_max)

				end do!igrid
				write (12,*) '   ' !blank line

				end do !ideriv
			end do !jpecies
		end do !ispecies
		deallocate (rhostore)


 		! Format Statements
! ===========================================================================
2222    format(I4,'  ', f10.7,'   ', f12.4)

! End Subroutine
! =============================================================================
        end subroutine eh_DOGS





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

	function  rint_eh(itype, ispecies, jspecies, isorp, d, rho, z1, &
     &                         z2, ideriv, index_2c)

! Argument Declaration
! ===========================================================================
! Argument Declaration
! ===========================================================================
    integer, intent (in) :: ispecies, jspecies     ! two centers species
    integer, intent (in) :: itype, isorp, ideriv   ! which interaction
    integer, intent (in) :: index_2c               ! which matrix element

    real d                              !< distance between the two centers
    real, intent(in) :: rho
    real, intent(in) :: z1, z2

	real r1, r2
	real rint_eh
	real psi1val, psi2val
	real vofr
	real fraction
	integer iexc

	 integer l1, m1, n1               ! quantum  numbers
     integer l2, m2, n2

	type (T_Fdata_cell_2c), pointer :: pFdata_cell
    type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
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


! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
!find psi1 value at point r1 as based on above and rho called for by adaptive_simpson
		psi1val = psiofr (r1, ispecies, n1)
		psi1val = psi1val ** 2/(4.0 * Pi)

!find psi2 value at point r2 as based on above and rho  called for by adaptive_simpson
		psi2val = psiofr (r2, jspecies, n2)
		psi2val = psi2val ** 2/(4.0 * Pi)

        fraction = PP_species(ispecies)%fraction
        iexc = PP_species(ispecies)%iexc

        vofr = dpotxc12 (rho, z1, d, ispecies, jspecies, iexc, fraction)


!Actual function (Ylm's are calculated above and multilplied after integration
		rint_eh = psi1val*vofr*psi2val*rho
! Format Statements
! ===========================================================================
!
! End Function
! =============================================================================
	end function rint_eh



! dpotxc12.f
! Program Description
! ===========================================================================
!
! This function computes the extended hubbard interaction: n1*n2*dnuxc(n1+n2)
! with the densities ni of atom in_i at the distance r_i from their centers.
!
! On input:  r, z: geometry information for the charge gradient
!
! On output:  dpotxc12
! ===========================================================================

! Code rewritten by:
! Otto F. Sankey (visiting)
! James P. Lewis
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
        real function dpotxc12 (r, z, d, ispecies, jspecies, iexc, fraction)
        implicit none



! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc, ispecies, jspecies

        real fraction
        real r, d
        real z

! Local Parameters and Data Declaration
! ===========================================================================
        real abohr
        parameter (abohr = 0.529177249d0)

! Local Variable Declaration and Description
! ===========================================================================

        real dens
        real densin
        real densp
        real denspp
        real denspz
        real densz
        real denszz
        real dexc2c
        real dnuxc2c
        real dnuxc2cs
        real dvxc2c
        !real hartree
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

		real dz, drho, dzraw
		integer nnrho, nnz
		integer, parameter :: nrho_points = 1601
		integer, parameter :: nz_points = 2501

! Procedure

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
        dnuxc2c = 0.00




!		All the extra functions in interpolate 2c
!cut down on some notations.
		rho2c =>  rhostore(ispecies, jspecies)%rho2c
        rhop2c => rhostore(ispecies, jspecies)%rhop2c
        rhopp2c => rhostore(ispecies, jspecies)%rhopp2c
        rhoz2c => rhostore(ispecies, jspecies)%rhoz2c
        rhozz2c => rhostore(ispecies, jspecies)%rhozz2c
        rhopz2c => rhostore(ispecies, jspecies)%rhopz2c


        dzraw = min(wf(ispecies)%dr_min,wf(jspecies)%dr_min)	!drr_rho(itype2))

        dz = dzraw*ixcgridfactor



! If the grid size is made too large or too small, then set to defaults.
        if (dz .gt. 0.05d0) dz = 0.05d0
        if (dz .lt. 0.002d0) dz = 0.002d0

        nnz = int((zmax - zmin)/dz) + 1

        drho = dz
        nnrho = int((rhomax - rhomin)/drho) + 1

! ===========================================================================
! Two-center piece: [n1 + n2(r,z)]*(exc[n1 + n2(r,z)] - vxc[n1 + n2(r,z)])
! ***************************************************************************
! Interpolate the density and gradients of the density at the given
! point (r, z).
        call interpolate2d (r, rhomin, rhomax, drho, z, zmin, zmax, dz,	&
                           nnrho, nrho_points, nnz, nz_points, rho2c,	&
                           dens)

! Convert to atomic units
        rin = r/abohr
        densin = dens*abohr**3

! Here energy and potential due to exchange and correlation are calculated.
! OFS+JPL 1999 extended hubbard addition. We add nu to the output list.
! We only use this for interaction=12, and at this time only Ceperly Alder.
! Here we ONLY have Ceper. Alder. (See create.f)
        call get_potxc2c (iexc, fraction, rin, densin, densp, denspp,	&
                         densz, denszz, denspz, dexc2c, dvxc2c,	&
                         dnuxc2c, dnuxc2cs)

! dnuxc2c is energy*volume, so must be in Hartree*abohr**3
!        hartree = 14.39975d0/abohr

        dpotxc12 = hartree*dnuxc2c*abohr**3

! Format Statements
! ===========================================================================

        return
        end function dpotxc12

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
