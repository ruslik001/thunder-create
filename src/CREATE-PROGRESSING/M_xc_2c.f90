! copyright info:
!
!                             @Copyright 2010
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

! M_xc
! Program Description
! ===========================================================================
!      This is a module containing a collection of exchange-correlation
! functionals - energies and potentials along with derivatives of both are
! calculated.
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
! Module Declaration
! ===========================================================================
        module M_xc

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! get_potxc_2c
! ===========================================================================
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
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine get_potxc_2c (iexc, xc_fraction, r, rho, rhop, rhopp,     &
     &                           rhoz, rhozz, rhopz, newexc, vpxc, dnuxc,    &
     &                           dnuxcs)
        implicit none

!! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc

        real, intent(in) :: xc_fraction
        real, intent(inout) :: r        ! radial distance

        ! density and derivatives
        real, intent(inout) :: rho      ! return zero value if rho is small
        real, intent(in) :: rhop
        real, intent(in) :: rhopp
        real, intent(in) :: rhoz
        real, intent(in) :: rhozz
        real, intent(in) :: rhopz

! Output
        real, intent(out) :: dnuxc
        real, intent(out) :: dnuxcs
        real, intent(out) :: newexc
        real, intent(out) :: vpxc

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ix

        real aln
        real exc                     ! exchange-correlation energy
        real ecp                     ! correlation potential
        real ex                      ! exchange potential
        real dec, dex                ! derivative of correlation and exchange
        real drvexc
        real fx
        real fxc
        real rs
        real x
        real zeta

! density and derivatives - spin cases
        real, dimension (2) :: d, dp, dpp, dz, dzz, dpz

! exchange and correlation potentials - for spin
        real, dimension (2) :: xpot
        real, dimension (2) :: cpot

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
          call ggaxrad_2c (3, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
          call ggacrad_2c (2, r, d, dp, dpp, dz, dzz,      cpot, dec)
          vpxc = xpot(1) + cpot(1)

! X Becke, C Perdew, generalized gradient approximation
        else if (iexc .eq. 5) then
          d = 0.5d0*rho
          dp = 0.5d0*rhop
          dpp = 0.5d0*rhopp
          dz = 0.5d0*rhoz
          dzz = 0.5d0*rhozz
          dpz = 0.5d0*rhopz
          call ggaxrad_2c (2, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
          call ggacrad_2c (3, r, d, dp, dpp, dz, dzz,      cpot, dec)
          vpxc = xpot(1) + cpot(1)

! XC Wigner-scaled LDA of PRA 46, R5320 (1992)
        else if (iexc .eq. 7) then
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
          call ggaxrad_2c (1, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
          call ggacrad_2c (1, r, d, dp, dpp, dz, dzz,      cpot, dec)
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
          call ggaxrad_2c (ix, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
          call ggacrad_2c (4 , r, d, dp, dpp, dz, dzz,      cpot, dec)
          if (iexc .ne. 12) then
            vpxc = xpot(1) + cpot(1)
          else
            vpxc = (1.0d0 - xc_fraction)*xpot(1) + cpot(1)
            dex = (1.0d0 - xc_fraction)*dex
          end if

! XC burke-perdew-ernzerhof gga 1996
        else if (iexc .eq. 6) then
          d = 0.5d0*rho
          dp = 0.5d0*rhop
          dpp = 0.5d0*rhopp
          dz = 0.5d0*rhoz
          dzz = 0.5d0*rhozz
          dpz = 0.5d0*rhopz
          call ggaxrad_2c (5, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
          call ggacrad_2c (5, r, d, dp, dpp, dz, dzz,      cpot, dec)
          vpxc = xpot(1) + cpot(1)

! If the improper iexc option was entered then the program will stop.
        else
          write (*,*) ' In get_potxc_2c.f90 - '
          write (*,*) ' stop: xc option not implemented', iexc
          stop
        end if

! Calculate the exchange energy by combining the exchange and correlation
! energies and subtracting the exchange/correlation potential energy
        newexc = dec + dex

! Format Statements
! ===========================================================================
! None

        return
        end subroutine get_potxc_2c


! ===========================================================================
! lsdavwn
! ===========================================================================
! Program Description
! ===========================================================================
!       This routine computes the exchange and correlation potenials and
! energies for the Vosko, Wilk, Nusair LSDA functional. Each spin component
! is considered.
!
! See
!      S.H. Vosko and L. Wilk and M. Nusair
!      Can. J. Phys., 58, 1200 (1980)
!
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
! Program Declaration
! ===========================================================================
        subroutine lsdavwn (rho, ex, ec, xpot, cpot, dnuxc, dnuxcs)
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
        real pi

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! =========================================================================
! Initialize some parameters
! Initialize
        pi = 4.0d0*atan(1.0d0)

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
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine lsdavwn


! ===========================================================================
! wigner
! ===========================================================================
! Program Description
! ===========================================================================
!       This routine computes the exchange and correlation potenials and
! energies for the Wigner interpolation formula, see
!
! E. Wigner, Phys. Rev. 46, 1002 (1934). Hartree a.u.
! and W.E. Pickett, Comp.Phys.Rep. 9, 115 (1989), pg. 187 : w2 = 7.79
!
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
! Program Declaration
! ==========================================================================
        subroutine wigner (rho, ex, fx, exc, fxc)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: rho

! Output
        real, intent (out) :: ex
        real, intent (out) :: exc

        real, intent (out) :: fx
        real, intent (out) :: fxc

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: ax = -0.738558766382022447d0
        real, parameter :: pif = 0.620350490899400087d0
        real, parameter :: w1 = -0.44d0
        real, parameter :: w2 = 7.79d0

        real, parameter :: eps = 1.0d-10

! Local Variable Declaration and Description
! ===========================================================================
        real rs
        real x
        real y

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
        if (rho .lt. eps) then
          ex = 0.0d0
          fx = 0.0d0
          exc = 0.0d0
          fxc = 0.0d0
        else
          x  = rho**(1.0d0/3.0d0)
          ex = ax*x
          fx = (4.0d0/3.0d0)*ex
          rs = pif/x
          y  = 1.0d0/(rs + w2)
          exc = ex + w1*y
          fxc = fx + w1*((4.0d0/3.0d0)*rs + w2)*y**2.0d0
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine wigner


! ===========================================================================
! ceperley_alder
! ===========================================================================
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
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine ceperley_alder (rho_in, epsx, potx, epsxc, potxc, drvexc, &
                                   dpotxc)
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
        real, parameter :: eps = 1.0d-15


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
        rho = rho_in

        if (rho .le. eps) then
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
          potxc = epsxc - rs*(0.15273333d0/rs**2                             &
                        + (0.02497128d0/sqrs + 0.01581427d0)/density**2)

          densityp =  1.0529d0/(2.0d0*sqrs) + 0.3334d0
          densitypp = -0.5d0*1.0529d0/(2.0d0*rs*sqrs)
          depsc = 0.1423d0*densityp/(density*density)
          ddepsc = - 2.0d0*0.1423d0*densityp*densityp/density**3             &
                   + 0.1423d0*densitypp/(density*density)
        else
          rsl = log(rs)
          epsxc = - 0.4581652d0/rs - 0.0480d0 + 0.0311d0*rsl - 0.0116d0*rs   &
                  + 0.002d0*rs*rsl
          potxc = epsxc - rs*(0.15273333d0/rs**2 + 0.01036667d0/rs           &
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

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine ceperley_alder


! ===========================================================================
! ggaxrad_2c
! ===========================================================================
! Program Description
! ===========================================================================
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
! Martin Fuchs, FHI der MPG, Berlin, 07-1992!
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
! Program Declaration
! ===========================================================================
        subroutine ggaxrad_2c (mode, rin, rho, rhop, rhopp, rhoz, rhozz,     &
     &                         rhopz, xpot, xen)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mode

        real, intent (in) :: rin

        real, intent (in), dimension (2) :: rho
        real, intent (in), dimension (2) :: rhop
        real, intent (in), dimension (2) :: rhopp
        real, intent (in), dimension (2) :: rhoz
        real, intent (in), dimension (2) :: rhozz
        real, intent (in), dimension (2) :: rhopz

! Output
        real, intent (out) :: xen

        real, intent (out), dimension (2) :: xpot

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: eps = 1.0d-15

! Local Variable Declaration and Description
! ===========================================================================
        integer ispin

! Value of density and corresponding derivatives at the point r, z
        real density
        real density_p, density_pp
        real density_z, density_zz
        real density_pz

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
! None

! Procedure
! ===========================================================================
! Initialize
        pi = 4.0d0*atan(1.0d0)

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
              density_p = 2.0d0*rhop(ispin)
              density_pp = 2.0d0*rhopp(ispin)
              density_z = 2.0d0*rhoz(ispin)
              density_zz = 2.0d0*rhozz(ispin)
              density_pz = 2.0d0*rhopz(ispin)
              fermik = 2.0d0*(3.0d0*pi*pi*density)**(1.0d0/3.0d0)

! s = abs(grad d)/(2kf*d)
! u = (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!      >>  grad(abs(grad d) has mixed derivatives ! <<
! v = (laplacian d)/(d*(2*kf)**2)
              s = sqrt(density_p**2 + density_z**2)/(fermik*density)
              u = (density_p**2*density_pp                                   &
     &             + 2.0d0*density_p*density_z*density_pz &
     &             + density_z**2*density_zz)/(s*density**3*fermik**4)
              v = (density_pp + density_p/r + density_zz)/(density*fermik*fermik)

              select case (mode)
               case (2)
                 call xbecke (density, s, u, v, ex, vx)
               case (3)
                 call exch (density, s, u, v, ex, vx)
               case (5)
                 call exchpbe (density, s, u, v, 1, 1, ex, vx)
              end select
            else
              stop 'ggaxrad_2c : mode improper'
            end if
            xpot(ispin) = vx
            xen = xen + rho(ispin)*ex
          end if
        end do

! energy
        xen = xen/max(rho(1) + rho(2), eps)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine ggaxrad_2c


! ===========================================================================
! ggacrad_2c
! ===========================================================================
! Program Description
! ===========================================================================
!      This routine calculates the correlation potential and energy density.
! Spherical symmetry is used. LSDA - GGA, Hartree a.u.
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
! Program Declaration
! ===========================================================================
        subroutine ggacrad_2c (mode, rin, rho, rhop, rhopp, rhoz, rhozz,     &
     &                         cpot, cen)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mode

        real, intent (in) :: rin

        real, intent (in), dimension (2) :: rho
        real, intent (in), dimension (2) :: rhop
        real, intent (in), dimension (2) :: rhopp
        real, intent (in), dimension (2) :: rhoz
        real, intent (in), dimension (2) :: rhozz

! Output
        real, intent (out) :: cen

        real, intent (out), dimension (2) :: cpot

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: eps = 1.0d-15

        real, parameter :: pisq3 = 2.96088132032680740d1
        real, parameter :: crs = 1.91915829267751281d0

! Local Variable Declaration and Description
! ===========================================================================
! Value of density and corresponding derivatives at the point r, z
        real density
        real density_p, density_pp
        real density_p11, density_p12, density_p22
        real density_pz

        real pi

! inputs
        real alfc, h, rs, t, uu, vv, ww
        real zet, ztp, fk, sk, g

! answers
        real ec, ecrs, eczet, vcdn, vcup, dvcdn, dvcup


! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize
        pi = 4.0d0*atan(1.0d0)
! LSDA
        density = rho(1) + rho(2)
        cen = 0.0d0

        if (density .le. eps) then
          cen = 0.0d0
          cpot(1) = 0.0d0
          cpot(2) = 0.0d0
        else
          if (mode .ne. 4) then
            zet = (rho(1) - rho(2))/density
            fk = (pisq3*density)**(1.0d0/3.0d0)
            rs = crs/fk
!
! GGA correction to LSDA
            if (mode .eq. 1) then
              call corlsd (rs, zet, ec, vcup, vcdn, ecrs, eczet, alfc)
              dvcup = 0.0d0
              dvcdn = 0.0d0
              h = 0.0d0
            else if (mode .eq. 2) then
              call corlsd (rs, zet, ec, vcup, vcdn, ecrs, eczet, alfc)
              density_p = rhop(1) + rhop(2)
              density_pp = rhopp(1) + rhopp(2)
              ztp = (rhop(1) - rhop(2) - zet*density_p)/density
              sk = 2.0d0*sqrt(fk/pi)
              g = ((1.0d0 + zet)**(2.0d0/3.0d0) + (1.0d0 - zet)**(2.00d0/3.0d0))/2.0d0
              t = abs(density_p)/(density*(2.0d0*sk*g)**2.0d0)
              uu = abs(density_p)*density_pp/(density**2.0d0*(2.0d0*sk*g)**3.0d0)
              vv = (density_pp + 2.0d0*density_p/rin)/(density*(2.0d0*sk*g)**2.0d0)
              ww = density_p*ztp/(density*(2.0d0*sk*g)**2.0d0)
              call corgga (rs, zet, t, uu, vv, ww, h, dvcup, dvcdn, fk, sk,  &
     &                     g, ec, ecrs, eczet)
            else if (mode .eq. 3) then
              call corlsd (rs, zet, ec, vcup, vcdn, ecrs, eczet, alfc)
              density_p  = rhop(1) + rhop(2)
              density_pp = rhopp(1) + rhopp(2)
              uu  = abs(density_p)*density_pp
              vv  = density_pp + 2.0d0*density_p/rin
              density_p11 = rhop(1)*rhop(1)
              density_p22 = rhop(2)*rhop(2)
              density_p12 = rhop(1)*rhop(2)
              if (rhop(1) .ne. 0.d0 .or. rhop(2) .ne. 0.0d0) &
                call corga86 (rho(1), rho(2), density_p11, density_p22,      &
     &                        density_p12, uu, vv, h, dvcup, dvcdn)
              else if (mode .eq. 5) then
              density_p = rhop(1) + rhop(2)
              density_pp = rhopp(1) + rhopp(2)
              ztp = (rhop(1) - rhop(2) - zet*density_p)/density
              sk = 2.0d0*sqrt(fk/pi)
              g = ((1.0d0 + zet)**(2.0d0/3.0d0) + (1.0d0 - zet)**(2.0d0/3.0d0))/2.0d0
              t = abs(density_p)/(density*(2.0d0*sk*g)**2.0d0)
              uu = abs(density_p)*density_pp/(density**2.0d0*(2.0d0*sk*g)**3.0d0)
              vv = (density_pp + 2.d0*density_p/rin)/(density*(2.0d0*sk*g)**2.0d0)
              ww = density_p*ztp/(density*(2.0d0*sk*g)**2.0d0)
              call corpbe (rs, zet, t, uu, vv, ww, 1, 1, ec, vcup, vcdn, h,  &
     &                     dvcup, dvcdn)
            else
              stop 'ggacrad : mode improper'
            end if
            cpot(1) = vcup + dvcup
            cpot(2) = vcdn + dvcdn
            cen = ec + h
          else if (mode .eq. 4) then
            call corlyp_2c (.true., rin, rho(1), rho(2), rhop(1), rhop(2),   &
     &                     rhopp(1), rhopp(2), rhoz(1), rhoz(2), rhozz(1),   &
                           rhozz(2), cen, cpot(1))
          else
            stop 'ggacrad : mode improper'
          end if
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine ggacrad_2c

! End Module
! =============================================================================
        end module
