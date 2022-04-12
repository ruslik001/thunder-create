! copyright info:
!
!                             @Copyright 2009
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
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

! Module Declaration
! ===========================================================================
        module M_xc3c_3c_Harris

		use M_species
		use M_atom_functions
		use M_integrals_3c
		use Exchange_Extra
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



		integer, allocatable :: iderorb(:)
		real, allocatable :: dqint(:,:), dqorb(:)
! module procedures
	contains

! ===========================================================================
! initialize_xc3c_3c_Harris
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
        subroutine initialize_xc_3c_HARRIS
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies, kspecies!< counters for number of species
        integer isorp                       !< counter for shells

        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over species for
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies

! For xc_3c_Harris
            do kspecies = 1, nspecies
              pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
              do isorp = 1, species(kspecies)%nssh
               pFdata_bundle%nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c + 1
              end do
            end do ! kspecies
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
        end subroutine initialize_xc_3c_HARRIS


! ===========================================================================
! evaluate_integral_3c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the general three-center
! matrix elements of the form <psi1|V(1)|psi2> for the spherical density.
! ===========================================================================
        subroutine xc_3c_Harris
        implicit none

        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ibcba, inaba, itheta    ! looping counters
        integer iexc
        integer index_3c, nME3c_max     ! different mu, nu types
        integer iounit                  ! file for writing
        integer isorp                   ! loop over shells
        integer ispmin, ispmax
        integer ispecies, jspecies, kspecies  ! species numbers
        integer nFdata_cell_3c          !< indexing of interactions

        real dbc, dna                  ! distances between centers
        real dbcx, dnax, distance_bc
        real rcutoff1, rcutoff2, rcutoff3   ! atom cutoffs

        real, dimension (P_ntheta) :: ctheta
        real, dimension (P_ntheta) :: ctheta_weights

        integer aspecies, assh

		real, pointer :: qpl (:, :, :)

        character (len = 40) filename
        character (len=25) Fdata_location

		type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Set up dqorb and iderorb.iderorb is the shell where
! the charge will be changed in the +-Q derivative stuff, and dqorb is
! the charge amount changed.
    	allocate (dqorb(nspecies))
    	allocate (iderorb(nspecies))
    	allocate (dqint(4,nspecies))
 		do aspecies = 1, nspecies

    	    iderorb(aspecies) = species(aspecies)%nssh
    	    dqorb(aspecies) = 0.5d0
			if (species(aspecies)%nssh .eq. 1) dqorb(aspecies) = 0.25d0
			do assh = 1, species(aspecies)%nssh
  			  	dqint(assh,aspecies) = dqorb(aspecies)/species(aspecies)%nssh
			end do
		end do


! Initialize the Legendre coefficients
        call gleg (ctheta, ctheta_weights, P_ntheta)

! Loop over all the species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
     		do kspecies = 1, nspecies
              pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
              pFdata_bundle%nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c + 1
              nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c
        	  pFdata_cell=>pFdata_bundle%Fdata_cell_3c(nFdata_cell_3c)

	          call make_munu_3c (nFdata_cell_3c, ispecies, jspecies, kspecies)

			  nME3c_max = pFdata_cell%nME

              write (filename, '("/",i2.2, "_munu_3c.",                      &
     &                               i2.2,".",i2.2,".",i2.2,".dat")')        &
     &               P_xc3c, species(ispecies)%nZ, species(jspecies)%nZ,  &
     &                          species(kspecies)%nZ
              Fdata_location = 'coutput'
!             inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
!             if (skip) cycle
              open (unit = 12, file = trim(Fdata_location)//trim(filename),   &
     &              status = 'unknown', position = 'append')

              ! Write out the mapping - stored in mu, nu, and mvalue
              write (12,*) (pFdata_cell%mu_3c(index_3c), index_3c = 1, nME3c_max)
              write (12,*) (pFdata_cell%nu_3c(index_3c), index_3c = 1, nME3c_max)
              write (12,*) (pFdata_cell%mvalue_3c(index_3c),                  &
     &                      index_3c = 1, nME3c_max)

              ! Set up grid loop control constants
              rcutoff1 = species(ispecies)%rcutoffA_max
              rcutoff2 = species(jspecies)%rcutoffA_max
              rcutoff3 = species(kspecies)%rcutoffA_max
		      dbc = rcutoff1 + rcutoff2
       	      dna = rcutoff3 + max(rcutoff1,rcutoff2)

              ispmin = 1
              ispmax = species(kspecies)%nssh
              allocate (qpl(P_ntheta, nME3c_max, ispmin:(ispmax - ispmin + 1)))

! ----------------------------------------------------------------------------
! Begin the big loops over dbc and dna.'
! ----------------------------------------------------------------------------
! Loop over all bondcharge distances.
              do ibcba = 1, nbc_xc
                dbcx = dfloat(ibcba - 1)*dbc/dfloat(nbc_xc - 1)

! for all bondcharges-- we set b=dbcx/2.
                distance_bc = dbcx/2.0d0

! Loop over all neutral atom distances.
! The distance is measured from the bondcharge center (b=dbcx/2)
                do inaba = 1, nna_xc
                  dnax = dfloat(inaba - 1)*dna/dfloat(nna_xc - 1)
                  call evaluate_integral_3c (nFdata_cell_3c, ispecies,       &
     &                                       jspecies, kspecies, iexc,       &
     &                                       ctheta, ctheta_weights, dbcx,   &
     &                                       dnax, nnr_rho, nntheta_rho,     &
     &                                       psiofr, phiint_xc3c, qpl)


! ----------------------------------------------------------------------------
! qpl's are the answer
! ------------------------------------- ---------------------------------------
! Write the qpl coefficients into the data files each combination in1, in2, in3,
! itheta(=1,ntheta_max), isorp gives an individual file.  The values for the
! different non-zero matrix elements of a given combination are written out
! after the index loop.
! ----------------------------------------------------------------------------
                  iounit = 12
                  do isorp = ispmin, ispmax
                    do itheta = 1, P_ntheta
                      iounit = iounit + 1
                      write (filename, '("/", "xc_3c_", i2.2, "_", i2.2,    &
     &                                   ".", i2.2, ".", i2.2, ".", i2.2,     &
     &                                   ".dat")')                            &
     &                  itheta, isorp, species(ispecies)%nZ,                  &
     &                  species(jspecies)%nZ, species(kspecies)%nZ
                      Fdata_location = 'coutput'

! Read the mapping - stored in mu, nu, and mvalue
                      open (unit = (iounit),                                  &
     &                      file = trim(Fdata_location)//trim(filename),      &
     &                      status = 'unknown', position = 'append')

                      write (iounit,*)                                        &
     &                  (qpl(itheta,index_3c,isorp), index_3c = 1, nME3c_max)
                    end do
                  end do
                end do   ! end of the dna loop
              end do  ! the end of the dbc loop

              ! close files
              iounit = 12
              do itheta = 1, P_ntheta
                do isorp = ispmin, ispmax
                  iounit = iounit + 1
                  close (unit = iounit)
                end do
              end do

    		end do  ! end loop over kspecies
		  end do  ! end loop over jspecies
        end do  ! end loop over ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
	    end subroutine xc_3c_Harris



	subroutine phiint_xc3c(itype, ispecies, jspecies, kspecies, ispmin, &
                                ispmax, iexc, r, ds, zr, r1, r2, rna, avgVmat)

    implicit none
    include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispecies, jspecies, kspecies, iexc, ispmin, ispmax, itype
        real, intent (in) :: r, ds, r2, r1, rna(3), zr

! Output
		real, intent (out) :: avgVmat(:,:)

! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
      integer, parameter :: numInthpP=25
      real hunderedfortieth
      real, parameter :: pi = 4.0d0*atan(1.0d0)
      real phifactor(-3:3)
      integer ip, ix
      !real hundredfortieth
	  type (T_Fdata_cell_3c), pointer :: pFdata_cell
      type (T_Fdata_bundle_3c), pointer :: pFdata_bundle
	  integer nME3c_max
	  real nphiinv, HP, phi
	  real inthpPhi(numInthpP*2+1)
	  real wPhi(numInthpP)
	  integer iphiMax
	  real cphi, sphi, xr, yr, r3
	  real vpot, prod, averagephi
	  integer, allocatable :: mleft(:), mright(:)
	  real psipsi


! Procedure
! ===========================================================================
      hunderedfortieth = 0.007142857142857142857142857142857143d0
      iphiMax = NuminthpP

      avgVmat = 0.00

      nphiinv=1.0d0/dfloat(numInthpP-1)
      HP = PI*nphiinv
 	  !hundredfortieth = 1/140

      inthpPhi(1)=0.0D0
      wPhi(1) = HP*41.0D0/140

      do ip= 2, numInthpP-1
         inthpPhi(ip) = (ip-1.0D0)*HP
         if (mod(ip,6).eq.2) wPhi(ip)=HP*hunderedfortieth*216.0D0
         if (mod(ip,6).eq.3) wPhi(ip)=HP*hunderedfortieth*27.0D0
         if (mod(ip,6).eq.4) wPhi(ip)=HP*hunderedfortieth*272.0D0
         if (mod(ip,6).eq.5) wPhi(ip)=HP*hunderedfortieth*27.0D0
         if (mod(ip,6).eq.0) wPhi(ip)=HP*hunderedfortieth*216.0D0
         if (mod(ip,6).eq.1) wPhi(ip)=HP*hunderedfortieth*82.0D0
      end do
      wPhi(iphiMax)= HP*hunderedfortieth*41.0D0

! JOM-test
!      inthpPhi(iphiMax)= 2.0D0*PI
      	inthpPhi(iphiMax)= PI



! ========================================================
!          Do integral over phi:
! ========================================================
!          The phi factors depend only on m.
			do ip = 1, iPhiMax
             phifactor(0)=1.0d0

             phi = inthpPhi(ip)
             cphi = cos(phi)
             sphi = sin(phi)
!
!            Do the phi integral, with the phifactors.
!            Note: We order the p-orbitals
!            here x,y,z (or pi,pi',sig), NOT z,x,y.
!            Note that px, and xz now are +1. And so on!
!
             phifactor(1)=cphi
             phifactor(-1)=sphi
! d's
             phifactor(2)=cphi*cphi-sphi*sphi
             phifactor(-2)=cphi*sphi

             xr=r*ds*cphi
             yr=r*ds*sphi
!
             r3=sqrt((xr-rna(1))**2+(yr-rna(2))**2+(zr-rna(3))**2)

!			write(*,*) 'hunderedfortieth, r3, phifactor', hunderedfortieth, r3, phifactor
!
! ---------------------------------------------------
!
             do  iX=ispmin,ispmax
				pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
        		pFdata_cell=>pFdata_bundle%Fdata_cell_3c(itype)
            	nME3c_max = pFdata_cell%nME


                vpot=dvxc3c (iexc, r, r2, r3, ispecies, jspecies, kspecies, ix+1)
!              note: dc,ds defined at the beginning of the theta loop
!

			   averagephi = 0.5 /PI
               prod=vpot*wphi(ip)		!baz *averagephi

               allocate (mleft(nME3c_max))
               allocate (mright(nME3c_max))

                 mleft = pFdata_cell%m_Mu(:)

                 mright = pFdata_cell%m_Nu(:)

                 avgVmat(ix,:)=avgVmat(ix,:) +   &
     &                           prod*phifactor(mleft)*phifactor(mright)



  			 deallocate (mleft)
             deallocate (mright)

             end do !IX1

		end do
! ===============================================
!          The end of the phi integral.

		end subroutine phiint_xc3c

! dvxc3c.f
! Program Description
! ===========================================================================
!
!       This subroutine computes vxc(n1+n2+n3) - vxc(n1+n2) with the
! densities ni of atom in_i at the distance ri from their centers.  Due to
! the small contribution of the three center case to the overall energy,
! only the lda level of theory will be used for this calculation.
!
! Input Variables
!
!     iexc:        Desired exchange/correlation functional
!     in1,in2,in3: Atomic indices
!     r1,r2,r3:    Radii for n1, n2, and n3
!     IX:          Derivative w.r.t the j1-th charge using the
!                  j2-th magnitude
!
! Common Variables
!
!     orbocc:      Orbital occupations
!     numorb:      Number of orbitals
!
! Output Variables
!
!     dvxc3c:      vxc(n1+n2+n3) - vxc(n1+n2)
!
!
! ===========================================================================
! Original code from Juergen Fritsch

! Code rewritten by:
! Richard B. Evans
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
        real function dvxc3c (iexc, r1, r2, r3, in1, in2, in3, IX)
        implicit none


! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc
        integer in1
        integer in2
        integer in3
        integer ix

        real r1
        real r2
        real r3

! Local Parameters and Data Declaration
! ===========================================================================
        real abohr
        parameter (abohr = 0.529177249d0)

! Local Variable Declaration and Description
! ===========================================================================

        real dens
        real dens1
        real dens2
        real dens3
        real densp
        real densp1
        real densp2
        real densp3
        real denspp
        real denspp1
        real denspp2
        real denspp3
        real dnuxc
        real dnuxcs
        real exc2c
        real exc3c
        real fraction
        real hartree
        real rin
        real vxc2c
        real vxc3c
        real drho

! Procedure
! ===========================================================================


! This exchange-correlation routine deals only with three-center interactions,
! therefore, we do not do exact exchange for the three-center terms and
! fraction should be initiallized to 1.0d0.
        fraction = 1.0d0

		drho = min(wf(in1)%dr_min, wf(in2)%dr_min, &
													wf(in3)%dr_min)


! Evaluate the spherically averaged wave functions for all orbitals of the
! two atoms
        call density_calc (iexc, ix, 1, in1, in2, in3, r1, drho, &					!Need to work out if drho is based overall or on each.
                          dens1, densp1, denspp1)									!Is it the max overall from all three species/ all species Period?
        call density_calc (iexc, ix, 2, in1, in2, in3, r2, drho, &
                          dens2, densp2, denspp2)
        call density_calc (iexc, ix, 3, in1, in2, in3, r3, drho, &
                          dens3, densp3, denspp3)

! Three-center-piece: dvxc3c[n1(r1) + n2(r2) + n3(r3)]
! Compute the exchange correlation potential for the three-center case
! ***************************************************************************
! The total three-center density is the sum of the three.
! Note that rin, densp, and denspp are not used in the LDA limits.
        dens = dens1 + dens2 + dens3
        densp = densp1 + densp2 + densp3
        denspp = denspp1 + denspp2 + denspp3

        rin = r1/abohr
        dens = dens*abohr**3
        densp = densp*abohr**4
        denspp = denspp*abohr**5
        call get_potxc1c (iexc, fraction, rin, dens, densp, denspp, &
                         exc3c, vxc3c, dnuxc, dnuxcs)


! Two-center-piece: dvxc2c[n1(r1) + n2(r2)]
! Compute the exchange correlation potential for the three-center case
! ***************************************************************************
! The total two-center density is the sum of the two - dens1 + dens2.
! Note that rin, densp, and denspp are not used in the LDA limits.
        dens = dens1 + dens2
        densp = densp1 + densp2
        denspp = denspp1 + denspp2

        rin = r1/abohr
        dens = dens*abohr**3
        densp = densp*abohr**4
        denspp = denspp*abohr**5
        call get_potxc1c (iexc, fraction, rin, dens, densp, denspp, &
                         exc2c, vxc2c, dnuxc, dnuxcs)

! Answers are in Hartrees convert to eV.
        hartree = 14.39975d0/abohr
        dvxc3c = hartree*(vxc3c - vxc2c)

! Format Statements
! ===========================================================================

        return
        end function

! get_potxc1c.f90
! Program Description
! ===========================================================================
!  This program will access the requested exchange and correlation
!  functionals from the following list:
!          1  LDA   Wigner
!          2  LDA   Hedin/Lundqvist
!          3  LDA   Ceperley/Alder Perdew/Zunger (1980)
!          4  GGA   Perdew/Wang (1991)
!          5  GGA   Becke (1988) X, Perdew (1986) C
!          6  GGA   Perdew/Burke/Ernzerhof (1996)
!          7  LDA   Zhao/Parr
!          8  LDA   Ceperley/Alder Perdew/Wang (1991)
!          9  GGA   Becke (1988) X, Lee/Yang/Parr (1988) C
!         10  GGA   Perdew/Wang (1991) X, Lee/Yang/Parr (1988) C
!         11  LSDA  Vosko/Wilk/Nusair (1980)
!         12  GGA   Becke (1988) X, Lee/Yang/Parr (1988) C
!                   with exact exchange
! The numerical value above is assigned to the variable iexc output
! exchange-correlation potential. This program has been modified to account
! for the fact that the density is a sum of two densities at two different
! centers.  Also, the potential is evaluated at one point in space and the
! integrals (matrix elements) are evaluated elsewhere.
!
! input
!    iexc        xc scheme
!    r           radial coordinate
!    rho         sum of atomic density in au
!    rhop        sum of atomic density gradient (with respect to r) in au
!    rhopp       sum of second gradient (with respect to r) in au
!
! output
!    vpxc        xc potential
!    newexc      xc energy
!    dnuxc
!    dnuxcs
! ===========================================================================
! Code written by:
! James P. Lewis
! Eduardo Mendez
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
!
! Program Declaration
! ===========================================================================
        subroutine FORget_potxc1c (iexc, fraction, r, rho, rhop, rhopp, newexc, &
     &                          vpxc, dnuxc, dnuxcs)
!        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iexc

        real, intent(in) :: fraction
        real, intent(inout) :: r
        real, intent(inout) :: rho
        real, intent(in) :: rhop
        real, intent(in) :: rhopp

! Output
        real, intent(out) :: newexc
        real, intent(out) :: vpxc
        real, intent(out) :: dnuxc
        real, intent(out) :: dnuxcs

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ix

        real aln
        real dec
        real dex
        real drvexc
        real ecp
        real ex
        real exc
        real fx
        real fxc
        real rh
        real rs
        real x
        real zeta

        real, dimension (2) :: cpot
        real, dimension (2) :: d
        real, dimension (2) :: dp
        real, dimension (2) :: dpp
        real, dimension (2) :: xpot

! Procedure
! ===========================================================================
! Initialize to zero.
        dnuxc = 0.0d0
        dnuxcs = 0.0d0

! If r is really small, then set to manageably small number.
        if (r .lt. 1.0d-4) r = 1.0d-4

! Rho must be positive, but not too small
        if (rho .lt. 1.0d-8) then
         rho = 0.0d0
         dnuxc = 0.0d0
         newexc = 0.0d0
         vpxc = 0.0d0
         return
        else if (rho .lt. 1.0d-5) then
         rho = 1.0d-5
        end if

! Determine exchange-correlation potentials
! exchange (X) only
        if (iexc .eq. 11) then
         zeta = 0.0d0
         d(1) = rho*0.5*(1 + zeta)
         d(2) = rho*0.5*(1 - zeta)
         call lsdavwn (d, dex, dec, xpot, cpot, dnuxc, dnuxcs)
         newexc = dex + dec
         vpxc = xpot(2) + cpot(2)                ! Holds for zeta = 0.0d0 only

! correlation (C) Wigner
        else if (iexc .eq. 1) then
         rh = rho
         call wigner (rh, ex, fx, exc, fxc)
         vpxc = fxc
         dex = ex
         dec = exc - ex

! C Hedin-Lundqvist
        else if (iexc .eq. 2) then
         rh = rho
         if (rh .ne. 0.0d0) then
          rs = 0.62035049d0*rh**(-1.0d0/3.0d0)
          x = rs/21.0d0
          aln = log(1.0d0 + 1.0d0/x)
          ecp = aln + (x**3*aln - x*x) + x/2 - 1.0d0/3.0d0
          dex = -0.458175d0/rs - 0.0225d0*ecp
          vpxc = -0.6109d0/rs - 0.0225d0*aln
         else
          dex = 0.0d0
          vpxc = 0.0d0
         end if

! XC Ceperley - Alder
! as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)
        else if (iexc .eq. 3) then
         rh = rho
         call ceperley_alder (rh, ex, fx, exc, fxc, drvexc, dnuxc)
         vpxc = fxc
         dex = ex
         dec = exc - ex

! X Perdew, C Perdew, generalized gradient approximation 1992
        else if (iexc .eq. 4) then
         rh = rho
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         call ggaxrad1c (3, r, d, dp, dpp, xpot, dex)
         call ggacrad1c (2, r, d, dp, dpp, cpot, dec)
         vpxc = xpot(1) + cpot(1)

! X Becke, C Perdew, generalized gradient approximation
        else if (iexc .eq. 5) then
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         call ggaxrad1c (2, r, d, dp, dpp, xpot, dex)
         call ggacrad1c (3, r, d, dp, dpp, cpot, dec)
         vpxc = xpot(1) + cpot(1)

! XC Wigner-scaled LDA of PRA 46, R5320 (1992)
        else if (iexc .eq. 7) then
         rh = rho
         call wigscaled(rh, ex, fx, exc, fxc)
         vpxc = fxc
         dex = ex
         dec = exc - ex

! XC Ceperley-Alder in Perdew-Wang parametrization of 1991
        else if (iexc .eq. 8) then
         d = 0.5d0*rho
         dp = 0.0d0
         dpp = 0.0d0
         call ggaxrad1c (1, r, d, dp, dpp, xpot, dex)
         call ggacrad1c (1, r, d, dp, dpp, cpot, dec)
         vpxc = xpot(1) + cpot(1)

! C Lee-Yang-Parr
        else if (iexc .eq. 9 .or. iexc .eq. 10 .or. iexc .eq. 12) then

! X Becke gga by default
         ix = 2

! X Perdew-Wang gga
         if (iexc .eq. 10) ix = 3
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         call ggaxrad1c (ix, r, d, dp, dpp, xpot, dex)
         call ggacrad1c (4, r, d, dp, dpp, cpot, dec)
         if (iexc .ne. 12) then
          vpxc = xpot(1) + cpot(1)
         else
          vpxc = (1.0d0 - fraction)*xpot(1) + cpot(1)
          dex = (1.0d0 - fraction)*dex
         end if

! XC burke-perdew-ernzerhof gga 1996
        else if(iexc .eq. 6) then
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         call ggaxrad1c (5, r, d, dp, dpp, xpot, dex)
         call ggacrad1c (5, r, d, dp, dpp, cpot, dec)
         vpxc = xpot(1) + cpot(1)

! If the improper iexc option was entered then the program will stop.
        else
         write (*,*) ' In get_potxc1c.f90 - '
         write (*,*) ' stop: xc option not implemented', iexc
         stop
        end if

! Calculate the exchange energy by combining the exchange and correlation
! energies and subtracting the exchange/correlation potential energy
! This comment seems to be old fashioned
        newexc = dec + dex

! Format Statements
! ===========================================================================
        return
        end subroutine

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
        end subroutine density_calc



! End Module
! =============================================================================
        end module
