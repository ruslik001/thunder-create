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
!      This is a module with the extra subroutines and functions for exchange,
!	Everything here is called by pot2cxc
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
       module Exchange_Extra
	   implicit none
    	 integer, parameter :: long = selected_real_kind (4,99)
         real, parameter :: abohr = 0.529177249d0
         real, parameter :: abohr15 = 0.3849477153d0
         real, parameter :: beta = 0.9d0
         real, parameter :: eq2 = 14.39975d0
         real, parameter :: Hartree = 14.39975d0/abohr
         real, parameter :: pi = 3.141592653589793238462643d0
         real, parameter :: ryd = 13.6057981d0
         real, parameter :: tolerance = 1.0d-5
! module procedures
        contains
! get_potxc.f90
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
!         11  LSDA  Volko/Wilk/Nusair (1980)
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
!    rhoz        sum of atomic density gradient (with respect to z) in au
!    rhozz       sum of second gradient (with respect to z) in au
!
! output
!    vpxc        xc potential
!    newexc      xc energy
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Eduardo Mendez
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
        subroutine get_potxc2c (iexc, fraction, r, rho, rhop, rhopp, rhoz,  &
     &                          rhozz, rhopz, newexc, vpxc, dnuxc, dnuxcs)
       ! use precision
        implicit none

!! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc

        real, intent(in) :: fraction
        real, intent(inout) :: r
        real, intent(inout) :: rho
        real, intent(in) :: rhop
        real, intent(in) :: rhopp
        real, intent(in) :: rhopz
        real, intent(in) :: rhoz
        real, intent(in) :: rhozz

! Output
        real, intent(out) :: dnuxc
        real, intent(out) :: dnuxcs
        real, intent(out) :: newexc
        real, intent(out) :: vpxc

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
        real rs
        real x
        real zeta

        real, dimension (2) :: cpot
        real, dimension (2) :: d
        real, dimension (2) :: dp
        real, dimension (2) :: dpp
        real, dimension (2) :: dpz
        real, dimension (2) :: dz
        real, dimension (2) :: dzz
        real, dimension (2) :: xpot

! Procedure
! ===========================================================================
! Initialize variables (only used in 3,11)
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

! Determine exchange-correlation potentials lsdavwn
        if (iexc .eq. 11) then
         zeta = 0.0d0
         d(1) = 0.5d0*rho*(1 + zeta)
         d(2) = 0.5d0*rho*(1 - zeta)
         call lsdavwn (d, dex, dec, xpot, cpot, dnuxc, dnuxcs)
         vpxc = xpot(1) + cpot(1)                ! Holds only for zeta = 0 only

! correlation (C) Wigner
        else if (iexc .eq. 1) then
         call wigner (rho, ex, fx, exc, fxc)
         vpxc = fxc
         dex = ex
         dec = exc - ex

! C Hedin-Lundqvist
        else if (iexc .eq. 2) then
         if (rho .ne. 0.0d0) then
          rs = 0.62035049d0*rho**(-1.0d0/3.0d0)
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
         call ceperley_alder (rho, ex, fx, exc, fxc, drvexc, dnuxc)
         vpxc = fxc
         dex = ex
         dec = exc - ex


! X Perdew, C Perdew, generalized gradient approximation 1992
        else if (iexc .eq. 4) then	!Concerned about this, rhop etc is much bigger than (2)- BaZ
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         dz = 0.5d0*rhoz
         dzz = 0.5d0*rhozz
         dpz = 0.5d0*rhopz
         call ggaxrad2c (3, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
         call ggacrad2c (2, r, d, dp, dpp, dz, dzz,      cpot, dec)
         vpxc = xpot(1) + cpot(1)

! X Becke, C Perdew, generalized gradient approximation
        else if (iexc .eq. 5) then
         d = 0.5d0*rho
         dp = 0.5d0*rhop
         dpp = 0.5d0*rhopp
         dz = 0.5d0*rhoz
         dzz = 0.5d0*rhozz
         dpz = 0.5d0*rhopz
         call ggaxrad2c (2, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
         call ggacrad2c (3, r, d, dp, dpp, dz, dzz,      cpot, dec)
         vpxc = xpot(1) + cpot(1)

! XC Wigner-scaled LDA of PRA 46, R5320 (1992)
        else if(iexc .eq. 7) then
         call wigscaled (rho, ex, fx, exc, fxc)
         vpxc = fxc
         dex = ex
         dec = exc - ex

! XC Ceperley-Alder in Perdew-Wang parametrization of 1991
        else if (iexc .eq. 8) then
         d = 0.5d0*rho
         dp = 0.0d0
         dpp = 0.0d0
         dz = 0.0d0
         dzz = 0.0d0
         dpz = 0.0d0
         call ggaxrad2c (1, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
         call ggacrad2c (1, r, d, dp, dpp, dz, dzz,      cpot, dec)
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
         dz = 0.5d0*rhoz
         dzz = 0.5d0*rhozz
         dpz = 0.5d0*rhopz
         call ggaxrad2c (ix, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
         call ggacrad2c (4 , r, d, dp, dpp, dz, dzz,      cpot, dec)
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
         dz = 0.5d0*rhoz
         dzz = 0.5d0*rhozz
         dpz = 0.5d0*rhopz
         call ggaxrad2c (5, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
         call ggacrad2c (5, r, d, dp, dpp, dz, dzz,      cpot, dec)
         vpxc = xpot(1) + cpot(1)

! If the improper iexc option was entered then the program will stop.
        else
         write (*,*) ' In get_potxc2c.f90 - '
         write (*,*) ' stop: xc option not implemented', iexc
         stop
        end if

! Calculate the exchange energy by combining the exchange and correlation
! energies and subtracting the exchange/correlation potential energy
        newexc = dec + dex

! Format Statements
! ===========================================================================

        return
        end subroutine get_potxc2c

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
        subroutine get_potxc1c (iexc, fraction, r, rho, rhop, rhopp, newexc, &
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

! lsdavwn.f90
! Program Description
! ===========================================================================
!       This routine computes the exchange and correlation potenials and
! energies for the Vosko, Wilk, Nusair LSDA functional. Each spin component
! is considered.
!
! See
!      S.H. VOSKO and L. WILK and M. NUSAIR
!      Can. J. Phys., 58, 1200 (1980)
!
! ===========================================================================
! Code written by:
! Eduardo Mendez
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine lsdavwn (rho, ex, ec, xpot, cpot, dnuxc, dnuxcs)
!        use constants
!        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in), dimension (2) :: rho

! Output
        real, intent (out) :: dnuxc
        real, intent (out) :: dnuxcs
        real, intent (out) :: ec
        real, intent (out) :: ex

        real, intent (out), dimension (2) :: cpot
        real, intent (out), dimension (2) :: xpot

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: Ap = 0.0621814d0
        real, parameter :: bp = 3.72744d0
        real, parameter :: cp = 12.9352d0
        real, parameter :: x0p = -0.10498d0

        real, parameter :: Aa = 0.033773728d0
        real, parameter :: ba = 1.13107d0
        real, parameter :: ca = 13.0045d0
        real, parameter :: x0a = -0.00475840d0

        real, parameter :: Af = 0.01554535d0
        real, parameter :: bf = 7.06042d0
        real, parameter :: cf = 18.0578d0
        real, parameter :: x0f = -0.32500d0

        real, parameter :: epsilon = 1.0d-10

! Local Variable Declaration and Description
! ===========================================================================
        real density
        real densitys

        real, dimension (3) :: cdpot
        real, dimension (3) :: xdpot

! spin polarization and derivatives
        real zeta, zp1, zp2, zp1p2, zpp1, zpp2
        real x, xp, xpp
        real g, gp, gpp
        real XXp, XXf, XXa , Qp, Qf, Qa, jp, jf, ja
        real ecP, ecF, ecA, ecPp, ecFp, ecAp, ecPpp, ecFpp, ecApp
        real cte, h, hp, hpp
        real ecpx, ecpz, ecppx, ecppz, ecpxpz
        real d1ec, d2ec, dd1ec, dd2ec, d1d2ec
        real exP, exPp, exPpp
        real expd, expz, exppd, exppz, expdpz
        real d1ex, d2ex, dd1ex, dd2ex, d1d2ex


! Allocate Arrays
! ===========================================================================

! Procedure
! =========================================================================
! Initialize some parameters
        stop ! zpp1 and zpp2 are not set
        density = rho(1) + rho(2)
        densitys = rho(1) - rho(2)
        zeta = densitys/density
        if (density .le. epsilon) then
         zeta = 0.0d0
         ec = 0.0d0
         ex = 0.0d0
         cpot = 0.0d0
         xpot = 0.0d0
         cdpot = 0.0d0
         xdpot = 0.0d0
         return
        end if

! Define simple derivatives
! *************************************************************************
        zp1 = 2.0d0*rho(1)/density**2
        zp2 = -2.0d0*rho(2)/density**2
        zp1p2 = 2.0d0*zeta/density**2

        x = (3.0d0/(4.0d0*pi*density))**(1.0d0/6.0d0)
        xp = - (1.0d0/6.0d0)*x/density
        xpp = 7.0d0*x/(36.0d0*density**2)

        g = (9.0d0/8.0d0)*((1.0d0 + zeta)**(4.0d0/3.0d0)                     &
          + (1.0d0 - zeta)**(4.0d0/3.0d0) - 2.0d0)
        gp = (3.0d0/2.0d0)*((1.0d0 + zeta)**(1.0d0/3.0d0)                    &
           - (1.0d0 - zeta)**(1.0d0/3.0d0))
        gpp = (1.0d0/2.0d0)*((1.0d0 + zeta)**(-2.0d0/3.0d0)                  &
            - (1.0d0 - zeta)**(-2.0d0/3.0d0))

! Intermediate variables
        XXp = x**2.0d0 + bp*x + cp
        XXf = x**2.0d0 + bf*x + cf
        XXa = x**2.0d0 + ba*x + ca
        Qp = (4.0d0*cp - bp*bp)**0.5d0
        Qf = (4.0d0*cf - bf*bf)**0.5d0
        Qa = (4.0d0*ca - ba*ba)**0.5d0
        jp = 2.0d0*log(x - x0p) - log(XXp)                                   &
           + 2.0d0*((2.0d0*x0p + bp)/Qp)*atan(Qp/(2.0d0*x + bp))
        jf = 2.0d0*log(x - x0f) - log(XXf)                                   &
           + 2.0d0*((2.0d0*x0f + bf)/Qf)*atan(Qf/(2.0d0*x + bf))
        ja = 2.0d0*log(x - x0a) - log(XXa)                                   &
           + 2.0d0*((2.0d0*x0a + ba)/Qa)*atan(Qa/(2.0d0*x + ba))

! epsilon derivatives
        ecP = Ap*(2.0d0*log(x) - log(XXp)                                    &
                 + (2.0d0*bp/Qp)*atan(Qp/(2.0d0*x + bp))                    &
                 - (bp*x0p/(x0p*x0p + bp*x0p + cp))*jp)
        ecF = Af*(2.0d0*log(x) - log(XXf)                                    &
                 + (2.0d0*bf/Qp)*atan(Qf/(2.0d0*x + bf))                    &
                 - (bf*x0f/(x0f*x0f + bf*x0f + cf))*jp)
        ecA = Aa*(2.0d0*log(x) - log(XXa)                                    &
                 + (2.0d0*ba/Qa)*atan(Qa/(2.0d0*x + ba))                    &
                 - (ba*x0a/(x0a*x0a + ba*x0a + ca))*ja)

        ecPp = 2.0d0*Ap*cp/(XXp*x) - 2.0d0*Ap*bp*x0p/((x - x0p)*XXp)
        ecFp = 2.0d0*Af*cf/(XXf*x) - 2.0d0*Af*bf*x0f/((x - x0f)*XXf)
        ecAp = 2.0d0*Aa*ca/(XXa*x) - 2.0d0*Aa*ba*x0a/((x - x0a)*XXa)

        ecPpp = - 2.0d0*Ap*cp*(3.0d0*x**2 + 2.0d0*bp*x + cp)/(x*XXp)**2      &
               + 2.0d0*Ap*bp*x0p*((2.0d0*x + bp)*(x - x0p) + XXp)           &
                      /(XXp*(x - x0p))**2
        ecFpp = - 2.0d0*Af*cf*(3.0d0*x**2 + 2.0d0*bf*x + cf)/(x*XXf)**2      &
               + 2.0d0*Af*bf*x0f*((2.0d0*x + bf)*(x - x0f) + XXf)           &
                      /(XXf*(x - x0f))**2
        ecApp = - 2.0d0*Aa*ca*(3.0d0*x**2 + 2.0d0*ba*x + ca)/(x*XXa)**2      &
               + 2.0d0*Aa*ba*x0a*((2.0d0*x + ba)*(x - x0a) + XXa)           &
                      /(XXa*(x - x0a))**2

        cte = 4.0d0/(9.0d0*(2.0d0**(1.0d0/3.0d0) - 1.0d0))

        h = cte*((ecF - ecP)/ecA) - 1.d0
        hp = cte*((ecFp - ecPp)/ecA - (ecF - ecP)*(ecAp/ecA))
        hpp = cte*((ecFpp - ecPpp)/ecA - (ecFp - ecPp)*ecAp/ecA**2           &
                 - (ecF - ecP)*ecApp/ecA - (ecFp - ecPp)*ecAp/ecA           &
                 + (ecF - ecP)*(ecAp/ecA)**2)


! Correlation functional ( and partials to z and x ):
        if (zeta .ne. 0.0d0) then
         ec = ecP + ecA*g*(1 + h*zeta**4)
        else
         ec = ecP
        end if

        ecpx = ecPp + ecAp*g*(1 + h*zeta**4) + ecA*g*hp*zeta**4
        ecpz = ecA*gp*(1.0d0 + h*zeta**4) + ecA*g*h*4*zeta**3

        ecppx = ecPp + ecApp*g*(1.0d0 + h*zeta**4) + 2.0d0*ecAp*g*hp*zeta**4 &
              + ecA*g*hpp*zeta**4
        ecppz = ecA*gpp*(1.0d0 + h*zeta**4) + ecA*gp*h*zeta**3               &
              + ecA*g*h*12.0d0*zeta**2
        ecpxpz = ecAp*gp*(1.0d0 + h*zeta**4) + ecA*gp*hp*zeta**4             &
               + ecAp*g*h*4.0d0*zeta**3 + ecA*g*hp*4.0d0*zeta**3

! Partial derivatives VWN exchanche functional
        d1ec = xp*ecpx + zp1*ecpz
        d2ec = xp*ecpx + zp2*ecpz

! Second partial derivatives
        dd1ec = xp**2*ecpp + 2.0d0*xp*zp1*ecpxpz + xpp*ecpx                  &
              + zp1*zp1*ecppz + zpp1*ecpz
        dd2ec = xp**2*ecpp + 2.0d0*xp*zp2*ecpxpz + xpp*ecpx                  &
              + zp2*zp1*ecppz + zpp2*ecpz
        d1d2ec = xp**2*ecpp+ xp*(zp1 + zp2)*ecpxpz + xpp*ecpx                &
               + zp1*zp2*ecppz + zp1p2*ecpz

! ****************************************************************************
!
!       VNN EXCHANGE FUNCTIONAL
!
! ****************************************************************************
        exP = (-3.0d0/2.0d0)*(3.0d0*density/pi)**(1.0d0/3.0d0)
        exPp = exP/(3.0d0*density)
        exPpp = - 2.0d0*exP/(3.0d0*density)**2

        ex =(1.0d0 + 4.0d0*g/9.0d0)*exP
        expd = ex/(3.0d0*density)
        exppd = -2.0d0*ex/(9.0d0*density**2)
        expz = exP*gp
        exppz = exP*gpp
        expdpz = exPp*gp

        d1ex = expd + zp1*expz
        d2ex = expd + zp2*expz

        dd1ex = exppd + 2.0d0*zp1*expdpz + expd + zp1*zp1*exppz + zpp1*expz
        dd2ex = exppd + 2.0d0*zp2*expdpz + expd + zp2*zp2*exppz + zpp2*expz
        d1d2ex = exppd + (zp1 + zp2)*expdpz + expd + zp1*zp2*exppz + zp1p2*expz

! Functions in Rydberg units - divide by factor of 2 to get Hartree
! ****************************************************************************
        xpot(1) = 0.5d0*(density*d1ex + ex)
        xpot(2) = 0.5d0*(density*d2ex + ex)
        cpot(1) = 0.5d0*(density*d1ec + ec)
        cpot(2) = 0.5d0*(density*d2ec + ec)
        ex = 0.5d0*ex
        ec = 0.5d0*ec

        cdpot(1) = 0.5d0*dd1ec
        cdpot(2) = 0.5d0*d1d2ec
        cdpot(3) = 0.5d0*dd2ec
        xdpot(1) = 0.5d0*dd1ex
        xdpot(2) = 0.5d0*d1d2ex
        xdpot(3) = 0.5d0*dd2ex

        dnuxc = 0.25d0*density*(xdpot(1) + 2.0d0*xdpot(2) + xdpot(3))        &
              + 0.5d0*(d1ec + d2ec) + 0.5d0*(d1ex + d2ex)                   &
              + 0.25d0*density*(cdpot(1) + 2.0d0*cdpot(2) + cdpot(3))

        dnuxcs = 0.25d0*density*(xdpot(1) - 2.0d0*xdpot(2) + xdpot(3))       &
                + 0.5d0*(d1ec - d2ec) + 0.5d0*(d1ex - d2ex)                 &
                + 0.25d0*density*(cdpot(1) - 2.0d0*cdpot(2) + cdpot(3))

        dnuxcs = 0.5d0*(ecA + 4.0d0*exP/9.0d0)/density


! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end subroutine lsdavwn


!***********************************************************************
! Wigner interpolation formula, E. Wigner, Phys. Rev. 46, 1002 (1934).
! Hartree a.u.
! see W.E. Pickett, Comp.Phys.Rep. 9, 115 (1989), pg. 187 : w2 = 7.79
!***********************************************************************
!
      subroutine wigner(rh,ex,fx,exc,fxc)
      implicit none

      real AX,EPS,EX,EXC,FX,FXC,PIF,RH,RS,THRD,W1,W2,X,Y

      data thrd,w1,w2,ax,pif,eps/.333333333333333333d0, &
                              -.44d0, &
                               .779d1, &
                              -.738558766382022447d0, &
                               .620350490899400087d0, &
                              1e-100/

!
      if(rh .lt. eps) then
        ex = 0.d0
        fx = 0.d0
        exc= 0.d0
        fxc= 0.d0
      else
        x  = rh**thrd
        ex = ax * x
        fx = 4.d0*thrd*ex
        rs = pif/x
        y  = 1.d0/(rs + w2)
        exc= ex + w1*y
        fxc= fx + w1*(4.d0*thrd*rs + w2)*y*y
      endif
!
      return
      end subroutine wigner


! ceperley_alder.f90
! Program Description
! ===========================================================================
!       This routine compute the ceperley-alder form of the LDA as
! parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981).  The units
! of this program are in atomic units, so the density but be changed to atomic
! units after input and the final answer converted to eV-Angstrom units.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine ceperley_alder (rho_in, epsx, potx, epsxc, potxc, drvexc, &
                                  dpotxc)
!        use constants
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: rho_in

! Output
        real, intent (out) :: dpotxc
        real, intent (out) :: drvexc
        real, intent (out) :: epsx
        real, intent (out) :: epsxc
        real, intent (out) :: potx
        real, intent (out) :: potxc

! Local Parameters and Data Declaration
! ===========================================================================
!        real, parameter :: tolerance = 1.0d-10

! Local Variable Declaration and Description
! ===========================================================================
        real density
        real densityp
        real densitypp
        real dpotc
        real dpotx
        real depsc
        real ddepsc
        real rho
        real rs
        real rsl
        real sqrs

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Convert density to Angstrom units.
!       rho = rho_in*abohr**3
        rho = rho_in

        if (rho .le. tolerance) then
         epsx = 0.d0
         potx = 0.d0
         epsxc = 0.d0
         potxc = 0.d0
         dpotxc = 0.0d0
         return
        end if

! Initialize some constants related to density.
        rs = 0.62035049d0/rho**(1.0d0/3.0d0)
        if (rho .lt. 0.23873241d0) then
         sqrs = sqrt(rs)
         density = 1.0d0 + 1.0529d0*sqrs + 0.3334d0*rs
         epsxc = -0.4581652d0/rs - 0.1423d0/density
         potxc = epsxc - rs*(0.15273333d0/rs**2                               &
                     + (0.02497128d0/sqrs + 0.01581427d0)/density**2)

         densityp =  1.0529d0/(2.0d0*sqrs) + 0.3334d0
         densitypp = -0.5d0*1.0529d0/(2.0d0*rs*sqrs)
         depsc = 0.1423d0*densityp/(density*density)
         ddepsc = - 2.0d0*0.1423d0*densityp*densityp/density**3              &
               + 0.1423d0*densitypp/(density*density)
        else
         rsl = log(rs)
         epsxc = - 0.4581652d0/rs - 0.0480d0 + 0.0311d0*rsl - 0.0116d0*rs    &
                + 0.002d0*rs*rsl
         potxc = epsxc - rs*(0.15273333d0/rs**2 + 0.01036667d0/rs            &
                     - 0.003866667d0 + 0.00066667d0*(1.0d0 + rsl))
         depsc = 0.0311d0/rs - 0.0116d0 + 0.0020d0*(rsl + 1.0d0)
         ddepsc = -0.0311d0/(rs*rs) + 0.0020d0/rs
        end if

! Exchange-only energy and potential
        epsx = - 0.7385587664d0*rho**(1.0d0/3.0d0)
        potx = 4.0d0/3.0d0*epsx

! Extended hubbard additions.
        drvexc = (potxc - epsxc)/rho

! Now dpotxc; we compute dpot/dn. We use dpot/dn = 2*dexc/dn + n*d2(exc)/dn2.
! Here dexc/dn = drvexc, and
! dpot/dn = (-2/(9*n*n))*ex + 4*rs/(9*n*n)*dec + rs*rs/(9*n*n)*ddec
! Let dpotc = dpot/dn = 4*rs/(9*n*n)*dec + rs*rs/(9*n*n)*ddec
        dpotc = (4.0d0*rs/(9.0d0*rho*rho))*depsc                             &
              + (rs*rs/(9.0d0*rho*rho))*ddepsc
        dpotx = - (2.0d0/(9.0d0*rho*rho))*epsx
        dpotxc = 2.0d0*drvexc + rho*(dpotx + dpotc)

! Convert output to eV-Angstrom units
!       epsx = epsx*Hartree
!       epsxc = epsxc*Hartree
!       potx = potx*Hartree
!       potxc = potxc*Hartree
!       drvexc = drvexc*Hartree*abohr**3
!       dpotxc = dpotxc*Hartree*abohr**3

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end subroutine ceperley_alder


! ggaxrad2c.f90
! Program Description
! ===========================================================================
!
!      This routine calculates the exchange potential and energy density.
! Spherical symmetry is used. LSDA - GGA
!
! input
!    mode = 1    LSDA
!    mode = 2    GGA-X Becke
!    mode = 3    GGA-X Perdew
!    mode = 5    GGA-X Burke-Perdew-Ernzerhof
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992
!
! ===========================================================================
! Code rewritten to FORTRAN 90 by:
! James P. Lewis
! Campus Box 7260
! Department of Biochemistry and Biophysics
! University of North Carolina
! Chapel Hill, NC 27599-7260
! FAX 919-966-2852
! Office telephone 919-966-4644
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine ggaxrad2c (mode, rin, rho, rhop, rhopp, rhoz, rhozz,     &
     &                        rhopz, xpot, xen)
 !       use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mode

        real, intent (in) :: rin

        real, intent (in), dimension (2) :: rho
        real, intent (in), dimension (2) :: rhop
        real, intent (in), dimension (2) :: rhopp
        real, intent (in), dimension (2) :: rhopz
        real, intent (in), dimension (2) :: rhoz
        real, intent (in), dimension (2) :: rhozz

! Output
        real, intent (out) :: xen

        real, intent (out), dimension (2) :: xpot

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: eps = 1.0d-15

! Local Variable Declaration and Description
! ===========================================================================
        integer ispin

        real density
        real densityp
        real densitypp
        real densitypz
        real densityz
        real densityzz
        real ex
        real fermik
        real pi
        real r
        real s
        real u
        real v
        real vx

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize
        pi = 3.141592653589793238462643D0

! If r is really small, then set to manageably small number.
        r = rin
        if (rin .lt. 1.0d-4) r = 1.0d-4

! exchange GGA, loop for up & down spin
        xen = 0.0d0
        do ispin = 1, 2
         if (rho(ispin) .le. eps) then
          xpot(ispin) = 0.0d0
         else
          density = 2.0d0*rho(ispin)
          if (mode .eq. 1) then
           call xlda (density, vx, ex)
          else if (mode .eq. 2 .or. mode .eq. 3 .or. mode .eq. 5) then
           densityp = 2.0d0*rhop(ispin)
           densitypp = 2.0d0*rhopp(ispin)
           densityz = 2.0d0*rhoz(ispin)
           densityzz = 2.0d0*rhozz(ispin)
           densitypz = 2.0d0*rhopz(ispin)
           fermik = 2.0d0*(3.0d0*pi*pi*density)**(1.0d0/3.0d0)

! s = abs(grad d)/(2kf*d)
! u = (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!      >>  grad(abs(grad d) has mixed derivatives ! <<
! v = (laplacian d)/(d*(2*kf)**2)
           s = sqrt(densityp**2 + densityz**2)/(fermik*density)
           u = (densityp**2*densitypp + 2.0d0*densityp*densityz*densitypz    &
     &          + densityz**2*densityzz)/(s*density**3*fermik**4)
           v = (densitypp + densityp/r + densityzz)/(density*fermik*fermik)

           select case (mode)
            case (2)
             call xbecke (density, s, u, v, ex, vx)
            case (3)
             call exch (density, s, u, v, ex, vx)
            case (5)
             call exchpbe (density, s, u, v, 1, 1, ex, vx)
           end select
          else
           stop 'ggaxrad2c : mode improper'
          end if
          xpot(ispin) = vx
          xen = xen + rho(ispin)*ex
         end if
        end do

! energy
        xen = xen/max(rho(1) + rho(2), eps)

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end subroutine ggaxrad2c

!**********************************************************************
!
! correlation potential & energy density
! spherical symmetry
! LSDA - GGA
! Hartree a.u.
!
! input
!    mode = 1    LSDA
!    mode = 2    GGA-C Perdew 91
!    mode = 3    GGA-C Perdew 86
!    mode = 4    GGA-C Lee-Yang-Parr 1988
!    mode = 5    GGA-C Burke-Perdew-Ernzerhof
!    r           radius
!    rh()        spin up/down density
!    rhp()       1st derivative of rh
!    rhpp()      2nd derivative of rh
!
! output
!    zet         spin polarization
!    cpot()       - " -  correlation potential
!    cen         correlation energy density
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992
!**********************************************************************
!
      subroutine ggacrad2c (mode, r, rh, rhp, rhpp, rhz, rhzz, &
                           cpot, cen)
!
      implicit none

      real r, rh, rhp, rhpp, rhz, rhzz, cpot, cen
      real thrd,pi,pisq3,thrd2,crs,eps
      real alfc,d,dp,dp11,dp12,dp22,dpp,dvcdn,dvcup,gks2
      real gks2sq,h,rs,t,uu,vcdn,vcup,vv,ww,zet,ztp
      real fk,sk,g,ec,ecrs,eczet

      integer mode

      dimension rh(2),rhp(2),rhpp(2),cpot(2)
      dimension rhz(2),rhzz(2)
      data thrd,pi /.333333333333333333d0,.314159265358979312d1/
      data pisq3,thrd2 /.296088132032680740d2,.666666666666666666d0/
      data crs,eps /1.91915829267751281d0,1.d-15/
!
! LSDA
      d = rh(1)+rh(2)
      cen=0.d0
      if(d .le. eps) then
        cen=0.d0
        cpot(1)=0.d0
        cpot(2)=0.d0
      else
        if(mode .ne. 4) then
          zet=(rh(1)-rh(2))/d
          fk=(pisq3*d)**thrd
          rs=crs/fk
!
! GGA correction to LSDA
          if(mode .eq. 1) then
            call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
            dvcup=0.d0
            dvcdn=0.d0
            h=0.d0
          else if(mode .eq. 2) then
            call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
            dp=rhp(1)+rhp(2)
            dpp=rhpp(1)+rhpp(2)
            ztp=(rhp(1)-rhp(2)-zet*dp)/d
            sk=2.d0*sqrt(fk/pi)
            g=((1.d0+zet)**thrd2 + (1.d0-zet)**thrd2) /2.d0
            gks2=2.d0*sk*g
            gks2sq=gks2*gks2
            t=abs(dp)/(d*gks2)
            uu=abs(dp)*dpp/(d*d*gks2sq*gks2)
            vv=(dpp +2.d0*dp/r)/(d*gks2sq)
            ww=dp*ztp/(d*gks2sq)
            call corgga(rs,zet,t,uu,vv,ww,h,dvcup,dvcdn, &
                       fk,sk,g,ec,ecrs,eczet)
          else if(mode .eq. 3) then
            call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
            dp  =rhp(1)+rhp(2)
            dpp =rhpp(1)+rhpp(2)
            uu  =abs(dp)*dpp
            vv  =dpp +2.d0*dp/r
            dp11=rhp(1)*rhp(1)
            dp22=rhp(2)*rhp(2)
            dp12=rhp(1)*rhp(2)
            if(rhp(1) .ne. 0.d0 .or. rhp(2) .ne. 0.d0) &
           call corga86(rh(1),rh(2),dp11,dp22,dp12,uu,vv,h,dvcup,dvcdn)
          else if(mode .eq. 5) then
            dp=rhp(1)+rhp(2)
            dpp=rhpp(1)+rhpp(2)
            ztp=(rhp(1)-rhp(2)-zet*dp)/d
            sk=2.d0*sqrt(fk/pi)
            g=((1.d0+zet)**thrd2 + (1.d0-zet)**thrd2) /2.d0
            gks2=2.d0*sk*g
            gks2sq=gks2*gks2
            t=abs(dp)/(d*gks2)
            uu=abs(dp)*dpp/(d*d*gks2sq*gks2)
            vv=(dpp +2.d0*dp/r)/(d*gks2sq)
            ww=dp*ztp/(d*gks2sq)
            call corpbe(rs,zet,t,uu,vv,ww,1,1, &
                          ec,vcup,vcdn,h,dvcup,dvcdn)
          else
            stop 'ggacrad : mode improper'
          endif
          cpot(1)=vcup+dvcup
          cpot(2)=vcdn+dvcdn
          cen=ec+h
        else if(mode .eq. 4) then
          call corlyp2c (.true., r, rh(1), rh(2), rhp(1), rhp(2), &
                        rhpp(1), rhpp(2), rhz(1), rhz(2), rhzz(1), &
                        rhzz(2), cen, cpot(1))
        else
          stop 'ggacrad : mode improper'
        endif
      endif
!
      return
      end subroutine ggacrad2c

! ggaxrad1c.f90
! Program Description
! ===========================================================================
!
!      This routine calculates the exchange potential and energy density.
! Spherical symmetry is used. LSDA - GGA
!
! input
!    mode = 1    LSDA
!    mode = 2    GGA-X Becke
!    mode = 3    GGA-X Perdew
!    mode = 5    GGA-X Burke-Perdew-Ernzerhof
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992
!
! ===========================================================================
! Code rewritten to FORTRAN 90 by:
! James P. Lewis
! Campus Box 7260
! Department of Biochemistry and Biophysics
! University of North Carolina
! Chapel Hill, NC 27599-7260
! FAX 919-966-2852
! Office telephone 919-966-4644
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine ggaxrad1c (mode, rin, rho, rhop, rhopp, xpot, xen)
!        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mode

        real, intent (in) :: rin

        real, intent (in), dimension (2) :: rho
        real, intent (in), dimension (2) :: rhop
        real, intent (in), dimension (2) :: rhopp

! Output
        real, intent (out) :: xen

        real, intent (out), dimension (2) :: xpot

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: eps = 1.0d-15

! Local Variable Declaration and Description
! ===========================================================================
        integer ispin

        real density
        real densityp
        real densitypp
        real ex
        real fermik
        real pi
        real r
        real s
        real u
        real v
        real vx

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize
        pi = 3.141592653589793238462643D0

! If r is really small, then set to manageably small number.
        r = rin
        if (rin .lt. 1.0d-4) r = 1.0d-4

! exchange GGA, loop for up & down spin
        xen = 0.0d0
        do ispin = 1, 2
         if (rho(ispin) .le. eps) then
          xpot(ispin) = 0.0d0
         else
          density = 2.0d0*rho(ispin)
          if (mode .eq. 1) then
           call xlda (density, vx, ex)
          else if (mode .eq. 2 .or. mode .eq. 3 .or. mode .eq. 5) then
           densityp = 2.0d0*rhop(ispin)
           densitypp = 2.0d0*rhopp(ispin)
           fermik = 2.0d0*(3.0d0*pi*pi*density)**(1.0d0/3.0d0)

! s = abs(grad d)/(2kf*d)
! u = (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!      >>  grad(abs(grad d) has mixed derivatives ! <<
! v = (laplacian d)/(d*(2*kf)**2)
           s = abs(densityp)/(fermik*density)
           u = abs(densityp)*densitypp/(density*density*fermik**3)
           v = (densitypp + 2.0d0*densityp/r)/(density*fermik*fermik)

           select case (mode)
            case (2)
             call xbecke (density, s, u, v, ex, vx)
            case (3)
             call exch (density, s, u, v, ex, vx)
            case (5)
             call exchpbe (density, s, u, v, 1, 1, ex, vx)
           end select
          else
           stop 'ggaxrad1c : mode improper'
          end if
          xpot(ispin) = vx
          xen = xen + rho(ispin)*ex
         end if
        end do

! energy
        xen = xen/max(rho(1) + rho(2), eps)

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end subroutine

! ggacrad1c.f
! Program Description
! ===========================================================================
!
!      This routine calculates the correlation potential and energy density.
! Spherical symmetry is used. LSDA - GGA
!
! input
!    mode = 1    LSDA
!    mode = 2    GGA-X Becke
!    mode = 3    GGA-X Perdew
!    mode = 5    GGA-X Burke-Perdew-Ernzerhof
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992
!
! ===========================================================================
! Code rewritten to FORTRAN 90 by:
! James P. Lewis
! Campus Box 7260
! Department of Biochemistry and Biophysics
! University of North Carolina
! Chapel Hill, NC 27599-7260
! FAX 919-966-2852
! Office telephone 919-966-4644
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine ggacrad1c (mode, rin, rho, rhop, rhopp, cpot, cen)
!        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mode

        real, intent (in) :: rin

        real, intent (in), dimension (2) :: rho
        real, intent (in), dimension (2) :: rhop
        real, intent (in), dimension (2) :: rhopp

! Output
        real, intent (out) :: cen

        real, intent (out), dimension (2) :: cpot

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: crs = 1.91915829267751281d0
        real, parameter :: eps = 1.0d-15
        real, parameter :: pi  = 3.14159265358979312d0
        real, parameter :: thrd = 0.333333333333333333d0
        real, parameter :: pisq3 = 29.6088132032680740d0

! Local Variable Declaration and Description
! ===========================================================================
        real alfc
        real density
        real densityp
        real densityp11
        real densityp12
        real densityp22
        real densitypp
        real ec
        real ecrs
        real eczet
        real fermik
        real g
        real gsfermik
        real h
        real r
        real rs
        real sfermik
        real t
        real uu
        real vv
        real ww
        real zet
        real ztp
        real fk
        real sk

        real, dimension (2) :: dvc, vc
        real, dimension (2) :: flip

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! If r is really small, then set to manageably small number.
        r = rin
        if (rin .lt. 1.0d-4) r = 1.0d-4

! LSDA
        density = rho(1) + rho(2)
        densityp = rhop(1) + rhop(2)
        densitypp = rhopp(1) + rhopp(2)
        cen = 0.0d0
        if (density .le. eps) then
         cen = 0.0d0
         cpot(1) = 0.0d0
         cpot(2) = 0.0d0
         return
        end if

        if (mode .ne. 4) then
         zet = (rho(1) - rho(2))/density
         ztp = (rhop(1) - rhop(2) - zet*densityp)/density
         fermik = (3.0d0*pi*pi*density)**(1.0d0/3.0d0)
         rs = crs/fermik
         call corlsd (rs, zet, ec, vc(1), vc(2), ecrs, eczet, alfc)

! GGA correction to LSDA
         select case (mode)
          case (1)
           dvc = 0.0d0
           h = 0.0d0
          case (2)
           sfermik = 2.0d0*sqrt(fermik/pi)
           g = ((1.0d0 + zet)**(2.0d0/3.0d0)                               &
     &        + (1.0d0 - zet)**(2.0d0/3.0d0))/2.0d0
           gsfermik = 2.0d0*sfermik*g
           t = abs(densityp)/(density*gsfermik)
           uu = abs(densityp)*densitypp/(density*density*gsfermik**3)
           vv = (densitypp + 2.0d0*densityp/r)/(density*gsfermik*gsfermik)
           ww = densityp*ztp/(density*gsfermik*gsfermik)
           fk=(pisq3*density)**thrd
           sk=2*sqrt(fk/pi)
           call corgga (rs, zet, t, uu, vv, ww, h, dvc(1), dvc(2),         &
     &                  fk,sk,g,ec,ecrs,eczet)

          case (3)
           uu = abs(densityp)*densitypp
           vv = densitypp + 2.0d0*densityp/r
           densityp11 = rhop(1)*rhop(1)
           densityp22 = rhop(2)*rhop(2)
           densityp12 = rhop(1)*rhop(2)
           if (rhop(1) .ne. 0.0d0 .or. rhop(2) .ne. 0.0d0)then
            call corga86 (rho(1),rho(2), densityp11, densityp22,             &
     &                    densityp12, uu, vv, h, dvc(1), dvc(2))
           end if
          case (5)
           sfermik = 2.0d0*sqrt(fermik/pi)
           g = ((1.0d0 + zet)**(2.0d0/3.0d0)                                 &
     &        + (1.0d0 - zet)**(2.0d0/3.0d0))/2.0d0
           gsfermik = 2.0d0*sfermik*g
           t = abs(densityp)/(density*gsfermik)
           uu = abs(densityp)*densitypp/(density*density*gsfermik**3)
           vv = (densitypp + 2.0d0*densityp/r)/(density*gsfermik*gsfermik)
           ww = densityp*ztp/(density*gsfermik*gsfermik)
           call corpbe (rs, zet, t, uu, vv, ww, 1, 1, ec, vc(1), vc(2), h,   &
     &                  dvc(1),dvc(2))
          end select
          cpot = vc + dvc
          cen = ec + h
        else if (mode .eq. 4) then
         flip = rhopp + 2.0d0*rhop/r
         call corlyp1c (.true., rho(1), rho(2), rhop(1), rhop(2), flip(1),   &
     &                  flip(2), cen, cpot(1), cpot(2))
        else
         stop 'ggacrad1c : mode improper'
        end if

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end subroutine
!
!
!c**********************************************************************
! LDA - scaled Wigner exchange-correlation
! from [Q. Zhao, R.G. Parr, PRA 46, R5320 (1992)] eqs. (5) and (7)
!
! Input
! rh        density
!
! Output
! ex        exchange energy per electron
! fx           ""    potential
!           ... based on continuum-LDA (electron-gas like)
!
! exc       exchange-correlation energy per electron
! fxc          ""        ""      potential
!
! Hartree a.u.
! Martin Fuchs, FHI der MPG, Berlin, 01-1993
!c**********************************************************************
      subroutine wigscaled(rh,ex,fx,exc,fxc)
      implicit none

      real AA,AX,BB,EX,EXC,FX,FXC,RH,THRD,X,XYLN,Y,Z

      data thrd,aa,bb,ax    /.333333333333333333d0, &
                           .93222d0, &
                           .947362d-2, &
                           .738558766382022406d0/
!
      if(rh .le. 0.d0) then
        rh = 0.d0
        ex = 0.d0
        fx = 0.d0
        exc= 0.d0
        fxc= 0.d0
      else
        x = bb*rh**thrd
        y = x/(x + 1.d0)
        z = -aa*x/bb
        xyln = x*log(y)
!
        exc = z * (1.d0 + xyln)
        fxc = thrd*z *(4.d0 + 5.d0*xyln + y)
!
! electron-gas based LDA exchange
!
        ex = ax*z/aa
        fx = 4.d0*thrd*ex
      endif
!
      return
      end subroutine wigscaled
!
!**********************************************************************
!  LDA exchange
!  Hartree a.u.
!***********************************************************************
      subroutine xlda(d,vx,ex)
!
      implicit none

      real AX,D,EPS,EX,THD,THD4,VX

      data ax,thd,thd4,eps/-.738558766382022406d0, &
                          .333333333333333333d0, &
                          .133333333333333333d1, &
                         1e-100/
!
      if(d .le. eps) then
        vx=0.d0
        ex=0.d0
      else
        ex=ax*d**thd
        vx=thd4*ex
      endif
!
      return
      end subroutine xlda
!
!
!  Becke exchange for a spin-unpolarized electronic system
!
!  Gradient-corrected exchange energy based on
!     [A.D. Becke, J.Chem.Phys.96, 2155, 1992].
!  The LSDA energy functional, obtained as E{n+,n-}=(E{2n+}+E{2n-})/2,
!     and the functional derivative formula are given by
!     [J.P. Perdew , PRB 33, 8800, 1986].
!     [J.P. Perdew , PRB 34, 7406, 1986].
!  see also [G.Ortiz ...,PRB 43, 6376 (1991)] eq. (A2)
!
!  Hartree a.u.
!
!  input
!  d            density
!  s            abs(grad d)/(2kf*d)
!  u            (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!           >>  grad(abs(grad d) has mixed derivatives ! <<
!  v            (laplacian d)/(d*(2*kf)**2)
!
!  output
!  ex           exchange energy per electron
!  vx           exchange potential
!
! Martin Fuchs, FHI der MPG, Berlin, 02-1993
!

      subroutine xbecke(d,s,u,v,ex,vx)

      implicit none

      real d,s,u,v,ex,vx,c,b,bb,ax,thrd,thrd4
      real dd1,ddi,f,fac,fs,fss,g,g1,x,y0,y1,y2

      data c / .779555417944150792d1/
      data b / .42d-2/
      data bb/-.451357747124625192d-2/
      data ax/-.738558766382022406d0/
      data thrd,thrd4/.333333333333333333d0,.1333333333333333333d1/

! exchange enhancement factor f
      x  = c*s
      y1 = 1.d0/sqrt(1.d0+x*x)
      y0 = log(x+1.d0/y1)
      y2 = -x*y1*y1*y1
      ddi= 1.d0/(1.d0 + 6.d0*b*x*y0)
      dd1= 6.d0*b*(y0+x*y1)
      g  = 1.d0 - 0.5d0*x*dd1*ddi
      fs = -2.d0*bb*c*c*ddi
      g1 = -3.d0*b*(y0+x*(3.d0*y1+x*y2-dd1*dd1*ddi/(6.d0*b)))
      fss= fs*c*(g1 - g*dd1)*ddi
      fs = fs*g
      f  = 1.d0 - bb*x*x*ddi

! LDA only
      fac= ax*d**thrd

! energy
      ex = fac*f

! potential
      vx = fac*(thrd4*f-(u-thrd4*s*s*s)*fss-v*fs)

      return
      end subroutine xbecke

!***********************************************************************

!
!  Perdew - Wang GGA91 exchange-correlation functional
!
!  J.P. Perdew et.al., Phys.Rev.B 46, 6671 (1992)
!
!  w/out gradients it's the new parametrization of the Ceperley-Alder
!  xc-energy data from
!
!  J.P. Perdew et.al., Phys.Rev.B 45, 13244 (1992)
!
!  pure spin singularity numerically removed
!
! original version by J P Perdew
! modified -- M. Fuchs (cmf), FHI Berlin, 07-1992
!
!
!  gga91 exchange for a spin-unpolarized electronic system
!  input d : density
!  input s:  abs(grad d)/(2*kf*d)
!  input u:  (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!  input v: (laplacian d)/(d*(2*kf)**2)
!  output:  exchange energy per electron (ex) and potential (vx)
!
      subroutine exch(d,s,u,v,ex,vx)
!
      implicit none
      real d,s,u,v,ex,vx,a1,a2,a3,a4,ax,a,b1,thrd,thrd4
      real f,fac,fs,fss,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11
      real s2,s3,s4

      data a1,a2,a3,a4/0.19645d0,0.27430d0,0.15084d0,100.d0/
      data ax,a,b1/-0.7385588d0,7.7956d0,0.004d0/
      data thrd,thrd4/0.333333333333d0,1.33333333333d0/
      fac = ax*d**thrd
      s2 = s*s
      s3 = s2*s
      s4 = s2*s2
      p0 = 1.d0/sqrt(1.d0+a*a*s2)
      p1 = log(a*s+1.d0/p0)
      p2 = exp(-a4*s2)
      p3 = 1.d0/(1.d0+a1*s*p1+b1*s4)
      p4 = 1.d0+a1*s*p1+(a2-a3*p2)*s2
      f = p3*p4
      ex = fac*f
!  local exchange option
!     ex = fac
!  energy done. now the potential:
      p5 = b1*s2-(a2-a3*p2)
      p6 = a1*s*(p1+a*s*p0)
      p7 = 2.d0*(a2-a3*p2)+2.d0*a3*a4*s2*p2-4.d0*b1*s2*f
      fs = p3*(p3*p5*p6+p7)
      p8 = 2.d0*s*(b1-a3*a4*p2)
      p9 = a1*p1+a*a1*s*p0*(3.d0-a*a*s2*p0*p0)
      p10 = 4.d0*a3*a4*s*p2*(2.d0-a4*s2)-8.d0*b1*s*f-4.d0*b1*s3*fs
      p11 = -p3*p3*(a1*p1+a*a1*s*p0+4.d0*b1*s3)
      fss = p3*p3*(p5*p9+p6*p8)+2.d0*p3*p5*p6*p11+p3*p10+p7*p11
      vx = fac*(thrd4*f-(u-thrd4*s3)*fss-v*fs)
!
!  local exchange option:
!     vx = fac*thrd4
      return
      end subroutine exch

!==========================================================================
! c original version by K Burke
! Perdew-Burke-Ernzerhof GGA
! see: J.P. Perdew, K. Burke, M. Ernzerhof, Phys Rev Lett 77, 3865 (1996).
! this collection contains the older PW91 and Becke exchange as well,
! everything not needed for the PBE GGA is commented out w/ "c--"
!
! Martin Fuchs, FHI der MPG, Berlin, 11-1996
!==========================================================================
!
!  --------------------------------------------------------------------
! |WARNING!  PBE is a simplification of PW91, which yields almost      |
! |identical numerical results with simpler formulas from a simpler    |
! |derivation.  If you should find significant DIFFERENCES between     |
! |PBE and PW91 results, please consult kieron@merlin.phy.tulane.edu   |
! |or perdew@mailhost.tcs.tulane.edu.  Thank you.                      |
!  --------------------------------------------------------------------
! Note: Neglects small grad (zeta) contributions to the correlation
! energy.
!
! Programs implement functional in PBE paper, July 1996 version.
!
!----------------------------------------------------------------------
      subroutine exchpbe(rho,S,U,V,lgga,lpot,EX,VX)
!----------------------------------------------------------------------
!  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  K Burke's modification of PW91 codes, May 14, 1996
!  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
!  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
!   (for U,V, see PW86(24))
!  input lgga:  (=0=>don't put in gradient corrections, just LDA)
!  input lpot:  (=0=>don't get potential and don't need U and V)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
!     {\bf 40},  3399  (1989) (E).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Formulas:
!     e_x[unif]=ax*rho^(4/3)  [LDA]
!     ax = -0.75*(3/pi)^(1/3)
!     e_x[PBE]=e_x[unif]*FxPBE(s)
!     FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! uk, ul defined after [a](13)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      implicit none

      real AX,EX,EXUNIF,FS,FSS,FXPBE,P0,RHO,S
      real s2,thrd,thrd4,u,uk,ul,um,v,vx

      integer  LGGA,LPOT

      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct LDA exchange energy density
      exunif = AX*rho**THRD
      if(lgga.eq.0)then
        ex=exunif
        vx=ex*thrd4
        return
      endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct PBE enhancement factor
      S2 = S*S
      P0=1.d0+ul*S2
      FxPBE = 1d0+uk-uk/P0
      EX = exunif*FxPBE
      if(lpot.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  Energy Done, now the potential:
!  find first and second derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
!  Fss=d Fs/ds
      Fs=2.d0*uk*ul/(P0*P0)
      Fss=-4.d0*ul*S*Fs/P0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate potential from [b](24)
      vx = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
      return
      end subroutine exchpbe


!******************************************************************

      subroutine corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
!  uniform-gas correlation of perdew and wang 1991
!  input: seitz radius (rs), relative spin polarization (zet)
!  output: correlation energy per electron (ec), up- and down-spin
!  potentials (vcup,vcdn), derivatives of ec wrt rs (ecrs) & zet (eczet)
!  output: correlation contribution (alfc) to the spin stiffness
      implicit none
      real rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc
      real gam,fzz,thrd,thrd4
      real ALFM,ALFRSM,COMM,EP,EPRS,EU,EURS,F,FZ,Z4
      data gam,fzz/0.5198421d0,1.709921d0/
      data thrd,thrd4/0.333333333333d0,1.333333333333d0/
      f = ((1.d0+zet)**thrd4+(1.d0-zet)**thrd4-2.d0)/gam
      call gcor(0.0310907,0.21370,7.5957,3.5876,1.6382, &
         0.49294,1.00,rs,eu,eurs)
      call gcor(0.01554535,0.20548,14.1189,6.1977,3.3662, &
         0.62517,1.00,rs,ep,eprs)
      call gcor(0.0168869,0.11125,10.357,3.6231,0.88026, &
         0.49671,1.00,rs,alfm,alfrsm)
!  alfm is minus the spin stiffness alfc
      alfc = -alfm
      z4 = zet**4
      ec = eu*(1.d0-f*z4)+ep*f*z4-alfm*f*(1.d0-z4)/fzz
!  energy done. now the potential:
      ecrs = eurs*(1.d0-f*z4)+eprs*f*z4-alfrsm*f*(1.d0-z4)/fzz
      fz = thrd4*((1.d0+zet)**thrd-(1.d0-zet)**thrd)/gam
      eczet = 4.d0*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu &
             -(1.d0-z4)*alfm/fzz)
      comm = ec -rs*ecrs/3.d0-zet*eczet
      vcup = comm + eczet
      vcdn = comm - eczet
      return
      end subroutine corlsd

!********************************************************************

      subroutine gcor(a,a1,b1,b2,b3,b4,p,rs,gg,ggrs)
!
      implicit none
      real a,a1,b1,b2,b3,b4,p,rs,gg,ggrs,p1,q0,rs12,rs32,rsp
      real q1,q2,q3
      p1 = p + 1.d0
      q0 = -2.d0*a*(1.d0+a1*rs)
      rs12 = sqrt(rs)
      rs32 = rs12**3
      rsp = rs**p
      q1 = 2.d0*a*(b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
      q2 = log(1.d0+1.d0/q1)
      gg = q0*q2
      q3 = a*(b1/rs12+2.d0*b2+3.d0*b3*rs12+2.d0*b4*p1*rsp)
      ggrs = -2.d0*a*a1*q2-q0*q3/(q1**2+q1)
      return
      end subroutine gcor
!
!  gga91 correlation
!  input rs: seitz radius
!  input zet: relative spin polarization
!  input t: abs(grad d)/(d*2.*ks*g)
!  input uu: (grad d)*grad(abs(grad d))/(d**2 * (2*ks*g)**3)
!  input vv: (laplacian d)/(d * (2*ks*g)**2)
!  input ww: (grad d)*(grad zet)/(d * (2*ks*g)**2
!  output h: nonlocal part of correlation energy per electron
!  output dvcup,dvcdn:  nonlocal parts of correlation potentials
!

      subroutine corgga(rs,zet,t,uu,vv,ww,h,dvcup,dvcdn, &
                       fk,sk,g,ec,ecrs,eczet)
!
      implicit none
      real rs,zet,t,uu,vv,ww,h,dvcup,dvcdn
      real xnu,cc0,cx,alf,c1,c2,c3,c4,c5,c6,a4,thrdm,thrd2,bet,delt
      real g3,g4,B,B2,BEC,BG,CC,CCRS,COEFF,COMM
      real FAC,FACT0,FACT1,FACT2,FACT3,FACT4,FACT5
      real GZ,H0,H0B,H0BT,H0RS,H0RST,H0T,H0TT,H0Z,H0ZT,H1,H1RS
      real H1RST,H1T,H1TT,H1Z,H1ZT,HRS,HRST,HT,HTT,HZ,HZT
      real PON,PREF,Q4,Q5,Q6,Q7,Q8,Q9,R0,R1,R2,R3,R4
      real RS2,RS3,RSTHRD,T2,T4,T6
      real fk,sk,g,ec,ecrs,eczet

      data xnu,cc0,cx,alf/15.75592d0,0.004235d0,-0.001667212d0,0.09d0/
      data c1,c2,c3,c4/0.002568d0,0.023266d0,7.389d-6,8.723d0/
      data c5,c6,a4/0.472d0,7.389d-2,100.d0/
      data thrdm,thrd2/-0.333333333333d0,0.666666666667d0/
      bet = xnu*cc0
      delt = 2.d0*alf/bet
      g3 = g**3
      g4 = g3*g
      pon = -delt*ec/(g3*bet)
      b = delt/(exp(pon)-1.d0)
      b2 = b*b
      t2 = t*t
      t4 = t2*t2
      t6 = t4*t2
      rs2 = rs*rs
      rs3 = rs2*rs
      q4 = 1.d0+b*t2
      q5 = 1.d0+b*t2+b2*t4
      q6 = c1+c2*rs+c3*rs2
      q7 = 1.d0+c4*rs+c5*rs2+c6*rs3
      cc = -cx + q6/q7
      r0 = (sk/fk)**2
      r1 = a4*r0*g4
      coeff = cc-cc0-3.d0*cx/7.d0
      r2 = xnu*coeff*g3
      r3 = exp(-r1*t2)
      h0 = g3*(bet/delt)*log(1.d0+delt*q4*t2/q5)
      h1 = r3*r2*t2
      h = h0+h1
!  local correlation option:
!     h = 0.0d0
!  energy done. now the potential:
      ccrs = (c2+2*c3*rs)/q7 - q6*(c4+2*c5*rs+3*c6*rs2)/q7**2
      rsthrd = rs/3.d0
      r4 = rsthrd*ccrs/coeff
!
      if(abs(zet) .ge. 1.d0) then
        if(zet .lt. 0.d0) zet=-1.d0+1.d-15
        if(zet .gt. 0.d0) zet=1.d0-1.d-15
      endif
      gz = ((1.d0+zet)**thrdm - (1.d0-zet)**thrdm)/3.d0
      fac = delt/b+1.d0
      bg = -3.d0*b2*ec*fac/(bet*g4)
      bec = b2*fac/(bet*g3)
      q8 = q5*q5+delt*q4*q5*t2
      q9 = 1.d0+2.d0*b*t2
      h0b = -bet*g3*b*t6*(2.d0+b*t2)/q8
      h0rs = -rsthrd*h0b*bec*ecrs
      fact0 = 2.d0*delt-6.d0*b
      fact1 = q5*q9+q4*q9*q9
      h0bt = 2.d0*bet*g3*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
      h0rst = rsthrd*t2*h0bt*bec*ecrs
      h0z = 3.d0*gz*h0/g + h0b*(bg*gz+bec*eczet)
      h0t = 2*bet*g3*q9/q8
      h0zt = 3.d0*gz*h0t/g+h0bt*(bg*gz+bec*eczet)
      fact2 = q4*q5+b*t2*(q4*q9+q5)
      fact3 = 2.d0*b*q5*q9+delt*fact2
      h0tt = 4.d0*bet*g3*t*(2.d0*b/q8-(q9*fact3/q8)/q8)
      h1rs = r3*r2*t2*(-r4+r1*t2/3.d0)
      fact4 = 2.d0-r1*t2
      h1rst = r3*r2*t2*(2.d0*r4*(1.d0-r1*t2)-thrd2*r1*t2*fact4)
      h1z = gz*r3*r2*t2*(3.d0-4.d0*r1*t2)/g
      h1t = 2.d0*r3*r2*(1.d0-r1*t2)
      h1zt = 2.d0*gz*r3*r2*(3.d0-11.d0*r1*t2+4.d0*r1*r1*t4)/g
      h1tt = 4.d0*r3*r2*r1*t*(-2.d0+r1*t2)
      hrs = h0rs+h1rs
      hrst = h0rst+h1rst
      ht = h0t+h1t
      htt = h0tt+h1tt
      hz = h0z+h1z
      hzt = h0zt+h1zt
      comm = h+hrs+hrst+t2*ht/6.d0+7.d0*t2*t*htt/6.d0
      pref = hz-gz*t2*ht/g
      fact5 = gz*(2.d0*ht+t*htt)/g
      comm = comm-pref*zet-uu*htt-vv*ht-ww*(hzt-fact5)
      dvcup = comm + pref
      dvcdn = comm - pref
!  local correlation option:
!     dvcup = 0.0d0
!     dvcdn = 0.0d0
      return
      end subroutine corgga

!***********************************************************************

!  gradient-correction to correlation energy from
!
!  [J.P.Perdew, PRB 33, 8822 (1986) and PRB 34, 7406 (1986)]
!
!  to correlation part of the Becke-Perdew gradient-corrected
!  xc-functional
!
!  input
!  d1,d2 : up/down spindensity
!  dp12..: grad(d1)*grad(d2) {* == vector product}, MUST NOT == 0
!  uu    : (grad d)*grad(abs(grad d)) , d = d1 + d2
!  vv    : laplacian d
!
!  output
!  dec  : correction to correlation energy per electron
!  dvcup :         - "" -            potential maj. spin
!  dvcdn :         - "" -            potential min. spin
!
!  Hartree a.u.
!
!  Martin Fuchs, FHI der MPG, Berlin, 07-1993
!
!
      subroutine corga86(d1,d2,dp11,dp22,dp12,uu,vv,dec,dvcup,dvcdn)
!
      implicit none

      real d1,d2,dp11,dp22,dp12,uu,vv,dec,dvcup,dvcdn
      real a1,a2,a3,a4,a5,a6,a7
      real t13,t23,t43,t53,t76,crt2
      real c,c1,c2,cp,d,d43,dd,dm13,dpnorm,dpnorm2,fi
      real uuu,vvv,www,zzz,zet,zz1,zz2
!
      data a1,a2,a3,a4,a5,a6,a7/ 2.568d-3, &
                                1.443307452d-2, &
                                2.843543831d-6, &
                                5.411317331d0, &
                                1.816419933d-1, &
                                1.763993811d-2, &
                                8.12908d-4/
      data t13,t23,t43,t53,t76/.33333333333333333d0, &
                             .66666666666666667d0, &
                             .13333333333333333d1, &
                             .16666666666666667d1, &
                             .11666666666666667d1/
      data crt2/.1587401052d1/
!
      d    = d1+d2
      dm13 = 1.d0/d**t13
      d43  = d**t43
!
! gradient expansion coefficient
      c1 = a1 + dm13*(a2+dm13*a3)
      c2 = 1.d0 + dm13*(a4+dm13*(a5+dm13*a6))
      c  = 1.667d-3 + c1/c2
!
      dpnorm = sqrt(dp11+dp22+2.d0*dp12)
      dpnorm2= dpnorm*dpnorm
!
      fi = a7*dpnorm/(c*d**t76)
!
! spin interpolation
      zet= (d1-d2)/d
      dd = sqrt(.5d0*((1.d0+zet)**t53 + (1.d0-zet)**t53))
!
! dC(n)/dn
      cp = -t13/(d43*c2*c2) &
         *(c2*(a2+2.d0*a3*dm13)-c1*(a4+dm13*(2.d0*a5+dm13*3.d0*a6)))
!
! spin-independent term
      www=( (fi-1.d0)*(cp/c-t43/d)+(fi-2.d0)*fi*(cp/c+t76/d) )*dpnorm2
      uuu=(3.d0-fi)*fi*uu/dpnorm
      vvv=(fi-2.d0)*vv
!
! spin dependent term
      zzz=crt2*.5d0*t53/(dd*d43)**2*(d1**t23-d2**t23)
      zz1= zzz*((1.d0-fi)*d2*dpnorm2-(2.d0-fi)*d*(dp12+dp22))
      zz2=-zzz*((1.d0-fi)*d1*dpnorm2-(2.d0-fi)*d*(dp12+dp11))
!
      dec=c/(dd*exp(fi)*d43)
!
      dvcup=dec*(www+uuu+vvv+zz1)
      dvcdn=dec*(www+uuu+vvv+zz2)
      dec=dec*dpnorm2/d
!
      return
      end subroutine corga86

!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      subroutine corpbe(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn,&
                       H,DVCUP,DVCDN)
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
!       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
!       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
!       :  UU,VV,WW, only needed for PBE potential
!       : lgga=flag to do gga (0=>LSD only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof,
!     {\sl Generalized gradient approximation made simple}, sub.
!     to Phys. Rev.Lett. May 1996.
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      implicit none

      real ALFM,ALFRSM,B,B2,BEC,BET,BG,COMM,DELT,DVCDN
      real DVCUP,EC,ECRS,ECZET,EP,EPRS,ETA,EU,EURS,F,FAC
      real FACT0,FACT1,FACT2,FACT3,FACT5,FZ,FZZ,G,G3,G4
      real GAM,GAMMA,GZ,H,HB,HBT,HRS,HRST,HT,HTT,HZ,HZT
      real Q4,Q5,Q8,Q9,RS,RSTHRD,RTRS,SIXTHM,T,T2,T4
      real T6,THRD,THRD2,THRD4,THRDM,UU,VCDN,VCUP,VV,WW,Z4,ZET
      real PON,PREF

      integer LGGA,LPOT

! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at
!          |zeta|=1.
      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(GAM=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*GAM))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gamma)
      parameter(eta=1.d-12)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=sqrt(rs)
      call gcor2(0.0310907,0.21370,7.5957,3.5876,1.6382, &
         0.49294,rtrs,EU,EURS)
      call gcor2(0.01554535,0.205480,14.1189,6.1977,3.3662, &
         0.62517,rtRS,EP,EPRS)
      call gcor2(0.0168869,0.111250,10.3570,3.62310,0.880260, &
         0.496710,rtRS,ALFM,ALFRSM)
      Z4 = ZET**4
      F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
             -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      if(lgga.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
      G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(EXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      H = G3*(BET/DELT)*log(1.D0+DELT*Q4*T2/Q5)
      if(lpot.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3.D0
      GZ=(((1.d0+zet)**2+eta)**sixthm- &
     ((1.d0-zet)**2+eta)**sixthm)/3.d0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      hRST = RSTHRD*T2*hBT*BEC*ECRS
      hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2.d0*BET*G3*Q9/Q8
      hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      return
      end subroutine corpbe
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      subroutine gcor2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
! slimmed down version of GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.

      implicit none

      real A,A1,B1,B2,B3,B4,rtrs,GG,GGRS,Q0,Q1,Q2,Q3

      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = log(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      return
      end subroutine gcor2

! corlyp2c.f90
! Program Description
! ========================================================================
!
!  The subroutine corlyp calculates the correlation energy using
!  the Lee, Yang and Parr generalized gradient approximation.
!
!
! Input Variables
!
! tpot .... T  evaluate correlation energy and potential
!           F  evaluate energy only (a posteriori scheme)
! x ....... dummy
! pa ...... spin up density
! pb ...... spin down density
! dpaofr .. 1st partial derivative of spin up with respect to r
! dpaofz .. 1st partial derivative of spin up with respect to z
! d2paofr . 2nd partial derivative of spin up with respect to r
! d2paofz . 2nd partial derivative of spin up with respect to z
! dpbofr .. 1st partial derivative of spin down with respect to r
! dpbofz .. 1st partial derivative of spin down with respect to z
! d2pbofr . 2nd partial derivative of spin down with respect to r
! d2pbofz . 2nd partial derivative of spin down with respect to z
!
! Output Variables
!
! ec ...... correlation energy per electron
! vp ...... correlation potential
!
! ===========================================================================
! Code written by:
! Richard B. Evans and James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine corlyp2c (tpot, r, pa, pb, dpaofr, dpbofr, d2paofr,     &
     &                       d2pbofr, dpaofz, dpbofz, d2paofz, d2pbofz,    &
     &                       ec, vp)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        logical  tpot

        real r
        real pa
        real pb
        real dpaofr
        real dpaofz
        real dpbofr
        real dpbofz
        real d2paofr
        real d2paofz
        real d2pbofr
        real d2pbofz

! Output
        real ec
        real vp

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: aa = 0.04918d0
        real, parameter :: bb = 0.132d0
        real, parameter :: cc = 0.2533d0
        real, parameter :: dd = 0.349d0
        real, parameter :: t13 = 1.0d0/3.0d0
        real, parameter :: t23 = 2.0d0/3.0d0
        real, parameter :: t43 = 4.0d0/3.0d0
        real, parameter :: t53 = 5.0d0/3.0d0
        real, parameter :: t73 = 7.0d0/3.0d0
        real, parameter :: t83 = 8.0d0/3.0d0

! Local Variable Delcaration and Descrition
! ======================================================================
! constants
        real Cf
        real pi
        real dep
        real dp_one
        real ep

! density and derivatives
        real p             ! Sum of spin up and spin down densities
        real dpr           ! Sum of spin up and down partials wrt r
        real dpz           ! Sum of spin up and down partials wrt z
        real d2pr          ! Sum of spin 2nd partials wrt r
        real d2pz          ! Sum of spin 2nd partials wrt z

        real gamma
        real dgammap       ! partial derivative of gamma wrt pa
        real dgammar       ! partial derivative of gamma wrt r
        real dgammaz       ! partial derivative of gamma wrt z
        real d2gammaz      ! 2nd Partial derivative of gamma wrt z
        real d2gammar      ! 2nd Partial derivative of gamma wrt r

        real Fofp
        real dFofp         ! partial derivative of Fofp wrt pa
        real dFofpz        ! partial derivative of Fofp wrt z
        real dFofpr        ! partial derivative of Fofp wrt r
        real d2Fofpr       ! 2nd Partial derivative of F wrt r
        real d2Fofpz       ! 2nd Partial derivative of F wrt z

        real Gofp
        real dGofp         ! partial derivative of Gofp wrt pa
        real dGofpz        ! partial of Gofp wrt z
        real dGofpr        ! partial of Gofp wrt r
        real d2Gofpr       ! 2nd partial of Gofp wrt r
        real d2Gofpz       ! 2nd partial of Gofp wrt z
        real d2Gofp        ! Laplacian of Gofp

        real tW            ! Weizsacker kinetic-energy density for p
        real taW           ! Weizsacker K.E. dens. for spin up
        real tbW           ! Weizsacker K.E. dens. for spin down

! Procedure
! ===========================================================================
! The notation used for our variables is similar to the original
! notation used by Lee, Yang and Parr. (Physical Review B 37 785 (1988))
! Initialize pi
        pi = 3.141592653589793238462643D0

        Cf = 0.3d0*(3.0d0*pi**2)**t23

! Here the option to calculate the potential is evaluated
        if (tpot) then
         p = pa + pb

         dpr = dpaofr + dpbofr
         dpz = dpaofz + dpbofz

         d2pr = d2paofr + d2pbofr
         d2pz = d2paofz + d2pbofz

         ep = p**(-t53)*exp(-cc*p**(-t13))
         dep = cc*t13*p**(-t43) - t53/p

         dp_one = 1.0d0 + dd*p**(-t13)

         gamma = 2.0d0*(1.0d0 - (pa**2 + pb**2)/p**2)
         dgammap = - 4.0d0*pa/p**2 - 2.0d0*(gamma - 2.0d0)/p
         dgammar = - 4.0d0*(pa*dpaofr + pb*dpbofr)/p**2                     &
                  - 2.0d0*dpr*(gamma - 2.0d0)/p
         dgammaz = - 4.0d0*(pa*dpaofz + pb*dpbofz)/p**2                     &
                  - 2.0d0*dpz*(gamma - 2.0d0)/p
         d2gammar =                                                         &
         - 4.0d0*(pa*d2paofr + dpaofr**2 + pb*d2pbofr + dpbofr**2)/p**2    &
         + 8.0d0*(pa*dpaofr + pb*dpbofr)*dpr/p**3                          &
         - 2.0d0*d2pr*(gamma - 2.0d0)/p                                    &
         - 2.0d0*dpr*(p*dgammar - (gamma - 2.0d0)*dpr)/p**2
         d2gammaz =                                                         &
         - 4.0d0*(pa*d2paofz + dpaofz**2 + pb*d2pbofz + dpbofz**2)/p**2    &
         + 8.0d0*(pa*dpaofz + pb*dpbofz)*dpz/p**3                          &
         - 2.0d0*d2pz*(gamma - 2.0d0)/p                                    &
         - 2.0d0*dpz*(p*dgammaz - (gamma - 2.0d0)*dpz)/p**2

         Fofp  = gamma/dp_one
         dFofp = (dgammap + dd*t13*p**(-t43)*Fofp)/dp_one
         dFofpr = (dgammar + dd*t13*p**(-t43)*Fofp*dpr)/dp_one
         dFofpz = (dgammaz + dd*t13*p**(-t43)*Fofp*dpz)/dp_one
         d2Fofpr =                                                          &
         (dp_one*(d2gammar - dd*t13*t43*p**(-t73)*Fofp*dpr**2              &
                  + dd*t13*p**(-t43)*(dFofpr*dpr + Fofp*d2pr))             &
          + dd*t13*p**(-t43)*dpr*(dgammar + dd*t13*p**(-t43)*Fofp*dpr))    &
         /dp_one**2
         d2Fofpz =                                                          &
         (dp_one*(d2gammaz - dd*t13*t43*p**(-t73)*Fofp*dpz**2              &
                  + dd*t13*p**(-t43)*(dFofpz*dpz + Fofp*d2pz))             &
         + dd*t13*p**(-t43)*dpz*(dgammaz + dd*t13*p**(-t43)*Fofp*dpz))     &
         /dp_one**2

         Gofp  = Fofp*ep
         dGofp = dFofp*ep + Gofp*dep
         dGofpr = dFofpr*ep + Gofp*dep*dpr
         dGofpz = dFofpz*ep + Gofp*dep*dpz
         d2Gofpr =                                                          &
         (d2Fofpr + dFofpr*dpr*dep)*ep + (dGofpr*dpr + Gofp*d2pr)*dep      &
          + Gofp*dpr**2*(t53/p**2 - cc*t13*t43*p**(-t73))
         d2Gofpz =                                                          &
         (d2Fofpz + dFofpz*dpz*dep)*ep + (dGofpz*dpz + Gofp*d2pz)*dep      &
          + Gofp*dpz**2*(t53/p**2 - cc*t13*t43*p**(-t73))

         d2Gofp = d2Gofpr + dGofpr/r + d2Gofpz

         vp = - aa*(dFofp*p + Fofp)                                         &
         - 2.0d0**t53*aa*bb*Cf                                             &
                     *(dGofp*(pa**t83 + pb**t83) + t83*Gofp*pa**t53)       &
         - aa*bb/4.0d0*(p*d2Gofp + 4.0d0*(dGofpr*dpr + dGofpz*dpz)         &
                        + 4.0d0*Gofp*(d2pr + dpr/r + d2pz)                 &
                        + dGofp*(p*(d2pr + dpr/r + d2pz)                   &
                        - (dpr**2 + dpz**2)))                              &
         - aa*bb/36.0d0                                                    &
             *(3.0d0*pa*d2Gofp + 4.0d0*(dpaofr*dGofpr + dpaofz*dGofpz)     &
               + 4.0d0*Gofp*(d2paofr + dpaofr/r + d2paofz)                 &
               + 3.0d0*dGofp*(pa*(d2paofr + dpaofr/r + d2paofz)            &
                              + pb*(d2pbofr + dpbofr/r + d2pbofz))         &
               + dGofp*(dpaofr**2 + dpaofz**2 + dpbofr**2 + dpbofz**2))
        else
         vp = 0.0d0
        end if

! correlation energy per electron
        tW = ((dpr**2 + dpz**2)/p - (d2pr + dpr/r + d2pz))/8.0d0
        taW = ((dpaofr**2 + dpaofz**2)/pa                                   &
              - (d2paofr + dpaofr/r + d2paofz))/8.0d0
        tbW = ((dpbofr**2 + dpbofz**2)/pb                                   &
              - (d2pbofr + dpbofr/r + d2pbofz))/8.0d0

        ec = 2.0d0**t23*Cf*(pa**t83 + pb**t83) - p*tW                       &
           + (pa*taW + pb*tbW)/9.0d0                                       &
           + (pa*(d2paofr + dpaofr/r + d2paofz)                            &
              + pb*(d2pbofr + dpbofr/r + d2pbofz))/18.0d0
        ec = - aa*gamma/dp_one*(p + 2.0d0*bb*p**(-t53)*exp(-cc*p**(-t13))*ec)
        ec = ec/p

! Format Statements
! ===========================================================================

        return
        end subroutine corlyp2c

! corlyp1c.f90
! Program Description
! ===========================================================================
!       Lee Yang Parr correlation energy functional one-dimensional densities
! only no provisions taken against division by zero.
!
! See e.g. C. Lee et al. Phys. Rev. B 37 785 (1988)
!
! Hartree A.U.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine corlyp1c (tpot, dp, dm, dp1, dm1, dp2, dm2, ec, vcp0, vcm0)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: dm, dm1, dm2
        real, intent (in) :: dp, dp1, dp2

        logical, intent (in) :: tpot

! Output
        real, intent (out) :: ec
        real, intent (out) :: vcm0
        real, intent (out) :: vcp0

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: aa = 0.04918d0
        real, parameter :: bb = 0.132d0
        real, parameter :: cc = 0.2533d0
        real, parameter :: dd = 0.349d0
        real, parameter :: c5 = 4.55779986d0
        real, parameter :: c6 = 1.0d0/72.0d0
        real, parameter :: c7 = 1.0d0/18.0d0
        real, parameter :: c8 = 0.125d0
        real, parameter :: t13 = 1.0d0/3.0d0
        real, parameter :: t89 = 8.0d0/9.0d0

! Local Variable Declaration and Description
! ===========================================================================
        real c1, c2, c3, c4, c9
        real chf
        real d0xt13, d0xt53
        real d0, d1, d2
        real dmt53, dpt53
        real dxsq
        real ga
        real gafm, gafp
        real gb
        real h
        real h2
        real hf
        real hff
        real sc
        real sc2
        real scf
        real t43, t53, t83
        real yafm, yafp
        real yb, yb1, yb2
        real ybfm, ybfp
        real yy1
        real yz, yz1, yz2
        real z1, z2
        real zfm, zfp
        real zeta

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize some parameters
        c1 = -4.0d0*aa
        c2 = dd
        c3 = 2.0d0*bb
        c4 = cc
        t53 = 5.0d0*t13
        t43 = 4.0d0*t13
        t83 = 2.0d0*t43
        c9 = t43 + t89

        d0 = dp + dm
        dxsq = 1.0d0/(d0*d0)
        d1 = dp1 + dm1
        d2 = dp2 + dm2
        d0xt13 = d0**(-t13)
        d0xt53 = d0xt13*d0xt13/d0
        dpt53 = dp**t53
        dmt53 = dm**t53

! Polarization factor
        zeta = c1*(dp*dm)*dxsq

! Scaling function
        sc = 1.0d0/(1.0d0 + c2*d0xt13)
        h = c3*d0xt53*exp(-c4*d0xt13)

! kinetic energy density expansion
        ga = c5*(dp*dpt53 + dm*dmt53)
        gb = c6*(dp1*dp1 - dp*dp2 + dm1*dm1 - dm*dm2) + c7*(dp*dp2 + dm*dm2) &
           + c8*(d0*d2 - d1*d1)

! Calculate potential
        if (tpot) then
         gafp = t83*c5*dpt53
         gafm = t83*c5*dmt53

         scf = t13*c2*d0xt13/d0*sc*sc
         sc2 = scf*(d2 + 2.0d0*(scf/sc - 2.0d0*t13/d0)*d1*d1)

         chf = t13*(c4*d0xt13 - 5.0d0)/d0
         hf = chf*h
         hff = h*(chf**2 + t13*(5.0d0 - 4.0d0*t13*c4*d0xt13)*dxsq)
         h2 = (hf*d2 + hff*d1*d1)

         zfp = (c1*dm - 2.0d0*zeta*d0)*dxsq
         zfm = (c1*dp - 2.0d0*zeta*d0)*dxsq
         yz = zeta/c1
         yy1 = dp*dm1 + dm*dp1
         yz1 = (yy1 - 2.0d0*yz*d1*d0)*dxsq
         yz2 = (2.0d0*yz*d1*d1 - 2.0d0*(yz1*d1 + yz*d2)*d0                   &
             - 2.0d0*d1*yy1/d0 + (dp*dm2 + 2.0d0*dp1*dm1 + dm*dp2))*dxsq
         z1 = c1*yz1
         z2 = c1*yz2

         yafp = sc*(d0*zfp + zeta) + zeta*d0*scf
         yafm = sc*(d0*zfm + zeta) + zeta*d0*scf

         yb = sc*zeta*h
         ybfp = sc*(h*zfp + zeta*hf) + zeta*h*scf
         ybfm = sc*(h*zfm + zeta*hf) + zeta*h*scf
         yb1 = sc*(h*z1 + zeta*hf*d1) + zeta*h*scf*d1
         yb2 = (sc*hf + h*scf)*d1*z1 + h*sc*z2 + (sc*z1 + zeta*scf*d1)*hf*d1 &
             + zeta*sc*h2 + (zeta*hf*d1 + h*z1)*scf*d1 + zeta*h*sc2

! Collect contributions
         vcp0 = yafp + ybfp*(ga + gb)                                        &
              + yb*(gafp + 2.0d0*c8*(c9*dp2 + 2.0d0*dm2))                   &
              + yb1*2.0d0*c8*(c9*dp1 + 2.0d0*dm1) + yb2*c8*(t43*dp + dm)
         vcm0 = yafm + ybfm*(ga + gb)                                        &
              + yb*(gafm + 2.0d0*c8*(c9*dm2 + 2.0d0*dp2))                   &
              + yb1*2.0d0*c8*(c9*dm1 + 2.0d0*dp1) + yb2*c8*(t43*dm + dp)
        else
         vcp0 = 0.0d0
         vcm0 = 0.0d0
        endif

! Correlation energy per electron
        ec = zeta*sc*(d0 + h*(ga + gb))/d0

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end subroutine



end module
