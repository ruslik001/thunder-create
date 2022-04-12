! Module Declaration
! ===========================================================================
        module M_XC3C

		use M_Species
		use M_atom_functions
		use Exchange_Extra
		use M_Three_Harris

		include 'interactions_2c.h'

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


	contains



 !===========================================================================
     subroutine xc3c

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer ispecies
        integer jspecies
        integer kspecies
        integer iexc         ! type of exchange-correlation
        integer index_max    ! maximum number of matrix elements
        integer ispmax
        integer ispmin
        integer ispnum
        integer nzx



        real dbc  ! maximum bond distance: rc1 + rc2
        real dna  ! maximum neutral atom distance: 2*rc3
        real rcutoff1  ! largest radius of i-th atom (in Ang.)
        real rcutoff2  ! 1, 2, 3 = left, right, neutral atom
        real rcutoff3

        real, dimension (ntheta_max) :: ctheta
        real, dimension (ntheta_max) :: ctheta_weights


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================


        integer, allocatable :: n1(:), l1(:), m1(:), n2(:), l2(:), m2(:)

        integer nME
        integer ibcba
        integer iam
        integer inaba
        integer index
        integer iounit
        integer isorp
        integer itheta
        integer jtheta
        integer nphi2ba
        integer nphi
        integer nrba
        integer nr
        integer ntheta
        integer nthetaba

        integer, dimension (inter_max) :: nabs

        real dbcx
        real dnax
        real distance_bc
        real cost
        real pl
        real plm
        real plmm
        real sint
        real temp

        real, dimension (ntheta_max) :: answer
! Final results after sin(theta) or cos(theta) factor.
        real, allocatable :: ggstore (:,:,:)
! Results from the integrator
        real, allocatable :: gmat(:,:)
        real, allocatable :: qpl(:,:,:)
        real, dimension (3) :: rna

        character (len = 40) filename
        character (len = 12) ftype
        character (len=25) Fdata_location

		type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle


		integer aspecies, assh
! Procedure
! ===========================================================================

	ispmin = 0
	ispmax = 0
	ispnum = 1

	ctheta = 0.0
	ctheta_weights = 0.0
	iexc = 3


	call Make_munu3c(2)
     write(*,*) 'back Munu'
     call gleg (ctheta, ctheta_weights, 5)
	write(*,*) 'gleg'
     do ispecies = 1, nspecies

     	do jspecies = 1, nspecies

     		do kspecies = 1, nspecies






! ====================================================================
! Set up dqorb and iderorb.iderorb is the shell where
! the charge will be changed in the +-Q derivative stuff, and dqorb is
! the charge amount changed.
    	allocate (dqorb(nspecies))
    	allocate (iderorb(nspecies))
    	allocate (dqint(4,nspecies))

		rcutoff1 = species(ispecies)%rcutoffA_max
		rcutoff2 = species(jspecies)%rcutoffA_max
		rcutoff3 = species(kspecies)%rcutoffA_max
		dbc = rcutoff1 + rcutoff2
       	dna = rcutoff3 + max(rcutoff1,rcutoff2)

		do aspecies = 1, nspecies

    	    iderorb(aspecies) = species(aspecies)%nssh
    	    dqorb(aspecies) = 0.5d0
			if (species(aspecies)%nssh .eq. 1) dqorb(aspecies) = 0.25d0
			do assh = 1, species(aspecies)%nssh
  			  	dqint(assh,aspecies) = dqorb(aspecies)/species(aspecies)%nssh
			end do
		end do
		pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_3c(2)
		nME = pFdata_cell%nME

	    allocate (n1(nME))
	    allocate (l1(nME))
	    allocate (m1(nME))

	    allocate (n2(nME))
	    allocate (l2(nME))
	    allocate (m2(nME))

	    allocate (gmat(ispmin:(ispmax-ispmin), nME))

		allocate (ggstore(ispmin:(ispmax-ispmin), nME, ntheta_max))
		allocate (qpl(ntheta_max, nMe, ispmin:(ispmax-ispmin)))
	    n1 = pFdata_cell%N_mu(:)
        l1 = pFdata_cell%L_mu(:)
        m1 = pFdata_cell%M_mu(:)

        n2 = pFdata_cell%N_nu(:)
        l2 = pFdata_cell%L_nu(:)
        m2 = pFdata_cell%M_nu(:)


        nphi2ba = nphi2ba_xc
        nrba = nrba_xc
        nthetaba = nthetaba_xc



        nr = 2*nrba + 1
        ntheta = 2*nthetaba + 1
        nphi = 2*nphi2ba + 1
        if (nr .gt. 1000) stop 'error---must re-dimension vrho'

! Now open the data files:
! what you get looks like this: bcna_01_01.14.06.14.dat
!                           or  xc3c_05_07.06.06.14.dat

         ftype = 'coutput/xc3c'




! ----------------------------------------------------------------------------
! 1. Call gauss-legendre routine to set up the angles.  The Legendre
!    polynomials go from l=0 and up, P_1 has one node, P_2 has two...
!    Sankey's theory has ntheta_max = 4, so here we have, ntheta_max + 1 = 5
! ----------------------------------------------------------------------------
! Set up nabs while we are here
        do index = 1, nME 	!index_max
         nabs(index) = abs(m1(index) - m2(index))
         write(*,*) nabs(index), abs(m1(index)- m2(index)),m1(index), m2(index)
         if (nabs(index) .gt. 2) stop 'WRONG NABS IN BCNA!!!!!'
        end do

! ----------------------------------------------------------------------------
! 2. Begin the big loops over dbc and dna.
! ----------------------------------------------------------------------------
! Loop over all bondcharge distances.
        do ibcba = 1, nbcba
         dbcx = dfloat(ibcba - 1)*dbc/dfloat(nbcba - 1)

! for all bondcharges-- we set b=dbcx/2.
         distance_bc = dbcx/2.d0

! Loop over all neutral atom distances.
! The distance is measured from the bondcharge center (b=dbcx/2)
         do inaba = 1, nnaba
          dnax = dfloat(inaba - 1)*dna/dfloat(nnaba - 1)
! ----------------------------------------------------------------------------
! 3. Since threecenter_integral internaly loops over ispmin to ispmax.
!    The thetas are roots of P_(ntheta_max+1)
! ----------------------------------------------------------------------------

          do itheta = 1, ntheta_max
           cost = ctheta(itheta)
           sint = sqrt(1 - cost*cost)

           rna(1) = sint*dnax
           rna(2) = 0.0d0
! -------------------------------------------------------
!           rna(3)=cost*dnax + distance_bc  ! old staff
! -------------------------------------------------------
! JOM : now the integral is done with the origin at the
!       center of the bondcharge
           rna(3)=cost*dnax               ! new staff
! -------------------------------------------------------
! threecenter_integral computes 3-c-integrals for ispmin...ispmax potentials
! and for all non-zero orbital combinations 1...index_max for one fixed
! dcb, rna configuration.  (index_max = index_max3c(itype1,itype2) )
! THIS IS WHERE ALL THE TIME IS SPENT

           call xc3c_integral(dbcx, rna, nr, ntheta, nphi, ispecies, jspecies, kspecies,   &
     &                               gmat, index_max, iexc,ispmin, ispmax)

! ----------------------------------------------------------------------------
! 4. Correct integrals as either type A or B
! ----------------------------------------------------------------------------
           do index = 1, index_max
            if (nabs(index) .eq. 1)then
             do isorp = ispmin, ispmax
!             type B: V=sin(theta)*Sum(l)* P*Q (nabs=1)
              if (sint .lt. 0.001d0) stop 'sin theta is zero!!!'
              ggstore(isorp,index,itheta) = gmat(isorp,index)/sint
             end do
            else
!            type A: V=Sum(l) P*Q. (nabs=0,2):do nothing
             do isorp = ispmin, ispmax
              ggstore(isorp,index,itheta) = gmat(isorp,index)
             end do
            end if
           end do
          end do


! Now all qpl coefficients are computed for a fixed dbc, dna pair, for all
! potentials  isorp = ispmin, ispmax and for all itheta(...as). The results
! for different matrix elements of one combination of isorp and itheta
! are written sequentially into one file.
! Begin the loop over all possible potentials on the third site (in3).
          do isorp = ispmin, ispmax

! ----------------------------------------------------------------------------
! 5. Now the time has come for the Gauss-Legendre integration.  We are still
!    in the isorp-loop.  Since the Gauss-Legendre integration has to be done
!    separately for each kind of potential and for each kind of orbital
!    combination, we need to loop here once more over the number of
!    non-vanishing three center integrals
! ----------------------------------------------------------------------------
           do index = 1, index_max
! Looping over the thetas, which are the roots of a Pl
            answer = 0.0d0
            do itheta = 1, ntheta_max
             plmm = 1.0d0
             plm = ctheta(itheta)
             temp = ctheta_weights(itheta)*ggstore(isorp,index,itheta)
             if (abs(temp) .lt. 1.0d-10) temp = 0.0d0
             answer(1) = answer(1) + plmm*temp
             answer(2) = answer(2) + plm*temp
             do jtheta = 3, ntheta_max
              pl = (plm*ctheta(itheta)                                  &
     &                 *(2.0d0*jtheta - 3.0d0) - (jtheta - 2.0d0)*plmm) &
     &             /(jtheta - 1.0d0)
              answer(jtheta) = answer(jtheta) + pl*temp
              plmm = plm
              plm = pl
             end do
            end do


! Normalize the coefficient, and write them out to qpl's.
            do itheta = 1, ntheta_max
             qpl(itheta,index,isorp) = answer(itheta)*(2.d0*itheta - 1.d0)*0.5d0
            end do
           end do !  end of the GL loop
          end do !  end of the isorp loop


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
           do itheta = 1, ntheta_max

            iounit = iounit + 1
            write (filename, '("/", "xc3c_", i2.2 ,"_",i2.2,"_",i2.2,".",i2.2,".",i2.2,".dat")')          &
     &       itheta, isorp, species(ispecies)%nZ, species(jspecies)%nZ, species(kspecies)%nZ

			Fdata_location = 'coutput'

!read the mapping - stored in mu, nu, and mvalue
			open (unit = (iounit), file = trim(Fdata_location)//trim(filename),   &
     &            status = 'unknown', position = 'append')

            write (iounit,*) (qpl(itheta,index,isorp), index = 1, index_max)
           end do
          end do
         end do   ! end of the dna loop
        end do! the end of the dbc loop

! Close files
        iounit = 12
        do itheta = 1, ntheta_max
         do isorp = ispmin, ispmax
          iounit = iounit + 1

          close (unit = iounit)
         end do
        end do
    	deallocate (dqorb)
    	deallocate (iderorb)
    	deallocate (dqint)
    	deallocate (gmat)
    	deallocate (ggstore)
    	deallocate (qpl)
    	deallocate (n1)
	    deallocate (l1)
	    deallocate (m1)

	    deallocate (n2)
	    deallocate (l2)
	    deallocate (m2)

    		end do !kspecies
		end do !jspecies
	end do  !ispecies

! Format Statements
! ===========================================================================
100     format (70('='))
110     format ('   itheta:  ',i2,'  isorp:  ',i2)
120     format (a70)
130     format (3(2x, i3))
140     format (2x, f9.4, 2x, i3)
150     format (i5, f9.4, ' <=== ', a4)
200     format (4d18.8)

        return
        end subroutine xc3c

!     ============================================================
       subroutine xc3c_integral(dbcx, rna, nr, ntheta, nphi, ispecies, jspecies, kspecies,   &
     &                               gmat, index_max, iexc,ispmin, ispmax)

!     ============================================================
!
!     Calculates the matrix elements of the neutral atom for the
!     specific configuration of dbcx (dist. between bc) and rna=position of
!     neutral atom.
!
!     ------------------------------------------------------------
!
      implicit none
!

!
!     ------------------------------------------------------------
!
!     Passed variables:
      integer numInthpR, numinthpT, numinthpP
!
! JOM : the following parameters are the grids for integration
! JOM : I think these numbers should be of the
!       type (6*N +1)  (Newton-Cotes 7-point rule)
! so far I recommend (note that symmetry (y,-y) is now used
! in the phi-integral
      parameter (numInthpR=49)
      parameter (numInthpT=49)
      parameter (numInthpP=25)

      integer nr, ntheta
!
      integer ispecies, jspecies, kspecies      ! types of the atoms
      integer iexc             ! exchange-correlation approximation
      integer nrr,nt,nphi      ! number of integration points
!
      integer index_max        ! maximal number of matrix elements
      integer interaction      ! number for the interaction
      integer ispmin,ispmax    ! range for isorp loop
      logical ispher           ! spherical approx.
!
      integer n1 (inter_max)   ! shell number of left  atom
      integer n2 (inter_max)   !                 right atom
      integer l1 (inter_max)   ! angular momentum of that shell
      integer l2 (inter_max)   !
      integer m1 (inter_max)   ! m-value in that shell
      integer m2 (inter_max)   !
!
      real rc1           ! i=1,2,3  for left,right, neutral atm
      real rc2
!
      real h
      real HR,HT,HP
      integer tempR, tempT, tempP

      real dbcx          ! bond charge distance  (A)
      real rna(3)        ! neutral atom location (A)
!
      integer tmpImid,imid
      real xmin
      real xxp
      real psipsi
!
      real, allocatable :: gmat(:,:)

      real psiofr,vpot,vnnaofr,dvxc3c
      external dvxc3c, vnnaofr,psiofr
!
!     ............................................................
!
!     interaction = 1: bcna  neutral atom
!                   2: xc3c  exchange correlation
!                   3: xc3c  exchange correlation (SNXC and OLSXC)
!
!     for the exchange correaltion case,
!
!     gmat(ix,*) stores the derivatives w.r.t charges:
!
!     gmat(1,*) : in1,0    gmat(2,*) : in1,-q1  gmat(3,*) : in1,+q1
!     gmat(4,*) : in2,-q2  gmat(5,*) : in2,+q2
!     gmat(6,*) : in3,-q3  gmat(7,*) : in3,+q3
!
!     ............................................................
!
!     ------------------------------------------------------------
!
!     Internal variables:
!
!

!
      real znormMU(inter_max),znormNU(inter_max),		&
     &       psiAmat(inter_max),psiBmat(inter_max)
!
!
      real thfactMU(0:3,-3:3),    &!  dimensions up to f-orbitals
     &       thfactNU(0:3,-3:3),	&
     &       phifactor(-3:3)
!
      !real avgVmat(0:10,inter_max)
      real fsimp(5000)
!
      real inthpR(numInthpR*2+1), inthpTheta(numInthpT*2+1),	&
     & inthpPhi(numInthpP*2+1)
      real wR(numInthpR*2+1),wTheta(numInthpT*2+1), &
     & wPhi(numInthpP*2+1)
!
      integer   numbphi
      integer irMax, itMax,iPhiMax
      parameter (numbphi=5000)
      real  paramR, paramT, paramP
      real    phiy(numbphi),      &
     &          cphiy(numbphi),sphiy(numbphi)

      real w1, w2
!
      integer nn1,nl1,nm1,nn2,nl2,nm2,i,ir,ix,ix1,it,ip,inm, &
     &        nmax,isorpX, NinthpR, NinthpT, NinthpP
!
! JOM r1
      real  r1
      real  sq3,sq15,pi,dr,dtheta,dphi,r,theta,phi, &
     &        rmin,rmax,dc1,ds1,dc,ds,zr,r2,dc2,ds2,xr,yr,r3,   &
     &        cphi,sphi,simp,simp2,simpson,    &
     &        averagephi,averagetheta,    &
     &        prod2,dsth,stuffmunu,thrd,nrrinv,ntinv,nphiinv

      real, allocatable :: prod(:)

      real sq12,sq4				!,hunderedfortieth

	  type (T_Fdata_cell_3c), pointer :: pFdata_cell
      type (T_Fdata_bundle_3c), pointer :: pFdata_bundle


!
! =====================================================================

!     The size of the matrix is determined by Nsh(in1) and Nsh(in2)
!     We do only those matrix elements that are not zero
!     The number of the non-zero matrix elements is INMAX
!
!     Normalization factors

      sq3=1.73205080756887729352744634150587d0
      sq15=3.87298334620741688517926539978240d0
      sq12=3.46410161513775458705489268301174d0
      sq4 = 2.0d0
      pi= 3.14159265358979323846264338327950d0
      thrd=0.333333333333333333333333333333333d0

      rmin=0.0d0

      ! added my murat manguoglu
 !     hunderedfortieth=0.007142857142857142857142857142857143d0 ! 1/140

! JOM

      rmax = 0.5d0*(rc1+rc2)

      NinthpR = numInthpR
      NinthpT = numInthpT
      NinthpP = numInthpP

      !modified by murat manguoglu
      irMax = NinthpR
      itMax = NinthpT
      iphiMax = NinthpP
      ! end of modification by murat manguoglu

      paramR = 2.0D0  ! 2*d/pi
      paramT = 2.0D0
      paramP = 1.0D0

      nrrinv=1.0d0/dfloat(irMax-1)
      ntinv=1.0d0/dfloat(itMax-1)
      nphiinv=1.0d0/dfloat(iphiMax-1)

      HR = (rmax-rmin)*nrrinv
      HT = PI*ntinv
! JOM-test

      HP = PI*nphiinv
! JOM
! I think I will change the integral in phi from (0,2pi) to (0,pi): latter

! now we choose d = Pi/4

!
!

!end of modification by murat manguoglu

! ===================================================================
! Here is the correct list:
!   m:     -2       -1         0        1         2
!          xy       yz      3z^2-r^2   xz       x^2-y^2
! sq15 * (  1        1        1/sq12    1       1/sq4  )
!
! ====================================================================
		pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies)

        pFdata_cell=>pFdata_bundle%Fdata_cell_3c(2)
		index_max = pFdata_cell%nME



      do inm=1,index_max
        nl1=pFdata_cell%l_Mu(inm)
        nl2=pFdata_cell%l_Nu(inm)
        if(nl1.eq.0)znormMU(inm)=1.0d0
        if(nl2.eq.0)znormNU(inm)=1.0d0
        if(nl1.eq.1)znormMU(inm)=sq3
        if(nl2.eq.1)znormNU(inm)=sq3
        if(nl1.eq.2)znormMU(inm)=sq15
        if(nl2.eq.2)znormNU(inm)=sq15
! First we calculate the m values.
        nm1=pFdata_cell%m_Mu(inm)
        nm2=pFdata_cell%m_nu(inm)
! working on d for mu
        if(nl1.eq.2)then
          if(nm1.eq.0)znormMU(inm)=znormMU(inm)/sq12
          if(nm1.eq.2)znormMU(inm)=znormMU(inm)/sq4
        end if
! working on d for nu
        if(nl2.eq.2)then
          if(nm2.eq.0)znormNU(inm)=znormNU(inm)/sq12
          if(nm2.eq.2)znormNU(inm)=znormNU(inm)/sq4
        end if
 	end do !inm=1,index_max

!
!     isorp=0,1,2,3, for neutral atom, s part, p part, d part.
!     We add 1 because of 0 being the neutral atom.
!

      if(nphi.gt.numbphi)then
        write(*,*)' nphi=',nphi
        write(*,*)' numbphi=',numbphi
        write(*,*)' In threecenter_integral-- redimension numbphi'
        stop 'error in threecenter_integral'
      end if
!
!     Set up some constants.
!     The phi integration does not depend on the r integration.
!     Therefore, it is done outside the r loop.
!
      do ip= 1, iphiMax
        phi=inthpPhi(ip)
        cphiy(ip)=cos(phi)
        sphiy(ip)=sin(phi)
	  end do
! ========================================================

 	    gmat = 0.0

		call Rint(ispmin, ispmax, ispecies, jspecies, kspecies, gmat, dbcx, rna)


!
!     Finally,the normalization factors for the different orbitals.
!     For instance a p orbital is sqrt(3) * x/r * r(r). The sqrt(3)
!     factor (and ones like it) are now included.
!
! JOM also averagetheta*averagephi = 1/(2 pi) is added now
!

!JOM
        allocate (prod(index_max))
        prod = 0.00

        prod=znormMU*znormNU*0.5d0/pi

        do ix=ispmin,ispmax
           gmat(ix,:)=gmat(ix,:)*prod
        end do



!
! ================================================================
!     SUMMARY
! We have computed gmat(isorp,mu,nu). mu, and nu are 1 to 9 for
! sp^3d^5.
! We are in molecular coordinates with z along sigma, x along pi, and
! y along pi'. Many matrix elements og gmat are zero. We have computed them
! anyway, and indeed find they are zero. Just to avoid at a later time,
! any trouble with roundoffs, I will now set those that are supposed to be
! zero, identically to zero.
! Here is the matrix, and the zero's are indicated.
!
!          \ s   x   y   z    3z^2-r^2    x^2-y^2    xz    xy      yz
!
!  s                 0                                     0        0
!
!  x                 0                                     0        0
!
!  y         0   0       0        0          0        0
!
!  z                 0                                     0        0
!
! 3z^2-r^2           0                                     0        0
!
! x^2-y^2            0                                     0        0
!
!  xz                0                                     0        0
!
!  xy        0   0       0        0          0        0
!
!  yz        0   0       0        0          0        0
!
!
       deallocate (prod)
       return
       end subroutine xc3c_integral


	subroutine Rint(ispmin, ispmax, ispecies, jspecies, kspecies, gmat, dbcx, rna)
	implicit none


! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispmin, ispmax
        integer, intent (in) :: ispecies, jspecies, kspecies

        real, intent (in) :: dbcx, rna(3)

! Output

		real, allocatable :: gmat(:,:)

! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
		integer, parameter :: numInthpR=49
		integer ir, irmax, inm, index_max, ix
      	real, allocatable :: thetamat(:,:)
      	real prod
		integer NinthpR				!, irMax
		real nrrinv, rmax, rmin
		real HR, rc1, rc2
		real, allocatable :: wR(:), inthpR(:)
		real r

	  	type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle
! Procedure
! ===========================================================================

      NinthpR = numInthpR
      irMax = NinthpR
      nrrinv=1.0d0/dfloat(irMax-1)
	  rc1 = species(ispecies)%rcutoffA_max			!PLEASE BAZ double-check these.
      rc2 = species(jspecies)%rcutoffA_max

      rmax = 0.5d0*(rc1+rc2)
      rmin = 0.00

	  HR = (rmax-rmin)*nrrinv

      allocate(inthpR(irMax))
      allocate(wR(irMax))

      inthpR(1)=rmin
      wR(1) = HR*hunderedfortieth*41.0D0
      do ir= 2, irMax-1
         inthpR(ir) = rmin+(ir-1.0D0)*HR
           if (mod(ir,6).eq.2) wR(ir)=HR*hunderedfortieth*216.0D0
           if (mod(ir,6).eq.3) wR(ir)=HR*hunderedfortieth*27.0D0
           if (mod(ir,6).eq.4) wR(ir)=HR*hunderedfortieth*272.0D0
           if (mod(ir,6).eq.5) wR(ir)=HR*hunderedfortieth*27.0D0
           if (mod(ir,6).eq.0) wR(ir)=HR*hunderedfortieth*216.0D0
           if (mod(ir,6).eq.1) wR(ir)=HR*hunderedfortieth*82.0D0
      end do

      wR(irMax)= HR*hunderedfortieth*41.0D0
      inthpR(irMax)= rmax

! ========================================================
!                 Do integral over r:
! ========================================================

		pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies)

        pFdata_cell=>pFdata_bundle%Fdata_cell_3c(2)
		index_max = pFdata_cell%nME

		allocate (thetamat((ispmin):(ispmax-ispmin), index_max))

        do ir= 1, irMax
        r=inthpR(ir)
!
!        dorbital addition OFS
!        Fireball2000: Jose y Jandro.
!        Zero out the array for theta integral.
!

             thetamat = 0.0d0


!        nli, nni tells us which shell is used for each atom.
!

! JOM we have to define psiAmat latter
!         do inm=1,index_max
!           nn1=n1(inm)
!           psiAmat(inm)=psiofr(in1,nn1,r)
!! jel-spher
!           if(ispher) psiAmat(inm) = sqrt(psiAmat(inm)**2.0d0)
!         end do
! end JOM



 		call thetaint(r, dbcx, ispmin, ispmax, ispecies, jspecies, kspecies, thetamat, rna)

!       ====================================================
!       The end of the integral over theta.
!       =====================================================

!       now finish off the r integral!

        prod=wR(ir)*r*r

	    gmat = gmat + prod * thetamat


	end do !ir
! ========================================================
	deallocate(thetamat)
	!return
	end subroutine rint







	subroutine thetaint(r, dbcx, ispmin, ispmax, ispecies, jspecies, kspecies, thetamat, rna)
	implicit none


! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispmin, ispmax
        integer, intent (in) :: ispecies, jspecies, kspecies
        real, intent (in) :: r, dbcx, rna(3)

! Output
      real, allocatable :: thetamat(:,:)
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
		integer, parameter :: numInthpT=49
		real ntinv
		real HT
		real inthpTheta(numInthpT*2+1)
		real, allocatable :: psiAmat(:), psiBmat(:)
		real, allocatable :: wTheta(:)
		real, allocatable :: avgVmat(:,:)
		integer it
		integer itMax, NinthpT
		real theta, dc, ds, r1,r2, rc1, rc2
		real zr, dc1, dc2, ds1, ds2
		integer inm, index_max
		integer nn1, nn2,nl1, nl2, nm1, nm2, ix
		integer iexc !iexc needs work Baz
		real dsth, prod, stuffmunu


		!     dimensions up to f-orbitals
        real thfactMU(0:3,-3:3),thfactNU(0:3,-3:3),phifactor(-3:3)


		type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle
! Procedure
! ===========================================================================
! ========================================================

		pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies)

        pFdata_cell=>pFdata_bundle%Fdata_cell_3c(2)
		index_max = pFdata_cell%nME

	 allocate (wTheta(numInthpT))
	 allocate (psiAmat(index_max))
	 allocate (psiBmat(index_max))
	 allocate (avgVmat(ispmin:(ispmax-ispmin), index_max))

     NinthpT = numInthpT
     itMax = NinthpT
     ntinv=1.0d0/dfloat(itMax-1)
     HT = PI*ntinv

      inthpTheta(1)=0.0D0
      wTheta(1) = HT*hunderedfortieth*41.0D0
      do it= 2, itMax-1
         inthpTheta(it) = (it-1.0D0)*HT
           if (mod(it,6).eq.2) wTheta(it)=HT*hunderedfortieth*216.0D0
           if (mod(it,6).eq.3) wTheta(it)=HT*hunderedfortieth*27.0D0
           if (mod(it,6).eq.4) wTheta(it)=HT*hunderedfortieth*272.0D0
           if (mod(it,6).eq.5) wTheta(it)=HT*hunderedfortieth*27.0D0
           if (mod(it,6).eq.0) wTheta(it)=HT*hunderedfortieth*216.0D0
           if (mod(it,6).eq.1) wTheta(it)=HT*hunderedfortieth*82.0D0
      end do
      wTheta(itMax)= HT*hunderedfortieth*41.0D0
      inthpTheta(itMax)= PI

! ========================================================
!              Do integral over theta:
! ========================================================
!
         do it=1, itMax
!
           theta=inthpTheta(it)
!          Theta stuff for atom 1
			!		write(*,*) 'GREPME Theta', theta
!=========================================================

           dc=cos(theta)
           ds=sin(theta)
!

       rc1 = species(ispecies)%rcutoffA_max			!PLEASE BAZ double-check these.
       rc2 = species(jspecies)%rcutoffA_max

	   r1 = r**2 + 0.25d0*(dbcx**2) + r*dbcx*dc
           if (r1 .le. 0.0d0) then
            r1 = 0.0d0
           else
            r1 = sqrt(r1)
           end if
           if(r1.gt.rc1) return      ! outside integration range

!           r2 = r**2 + dbcx**2 - 2*r*dbcx*dc
	   r2 = r**2 + 0.25d0*(dbcx**2) - r*dbcx*dc
           if (r2 .le. 0.0d0) then
            r2 = 0.0d0
           else
            r2 = sqrt(r2)
           end if
           if(r2.gt.rc2) return      ! outside integration range

!                                       ! other theta values might work
!
           do inm=1,index_max
             nn1 = pFdata_cell%N_mu(inm)
             nn2 = pFdata_cell%N_nu(inm)
             psiAmat(inm)=psiofr(r1,ispecies,nn1)
             psiBmat(inm)=psiofr(r2,jspecies,nn2)
! BaZ- As far as I can rememebr, sphere does not count in the interaction == 2 case.
! jel-spher
!             if(ispher) then
!              psiAmat(inm) = sqrt(psiAmat(inm)**2.0d0)
!              psiBmat(inm) = sqrt(psiBmat(inm)**2.0d0)
!              end if
           enddo
!
           zr=r*dc
!          Theta stuff for atom 1.
!          Be careful for r1 very small.
!          Find cos(theta1), sin(theta1).
           if(r1.gt.0.00001)then
             dc1=(zr+0.5d0*dbcx)/r1
             ds1=ds*r/r1
           else
             dc1=1.0d0
             ds1=0.0d0
           end if
!
!          Theta stuff for atom 2.
!          Be careful for r2 very small.
!          Find cos(theta2), sin(theta2).
           if(r2.gt.0.00001)then
             dc2=(zr-0.5d0*dbcx)/r2
             ds2=ds*r/r2
           else
             dc2=1.0d0
             ds2=0.0d0
           end if

!
! -------------------------------------------------------
!          Theta factors for A and B.
! -------------------------------------------------------
!          Use the (l,m) notation. For example 3z^2-1 becomes
!          (2,0), px becomes (1,1), and so on.
!
!          ATOM A .........................................
!
!          S
           thfactMU(0,0)=1.0d0
!
!          P
!          Note: We order the orbitals here x,y,z (or pi,pi',sig)
!
!          watch it: theta factors for d-orbitals should
!                    also contain the m-dependency of the
!                    prefactors znormMU and znormNU.
!                    This is not the case so far.
!
           thfactMU(1,1)=ds1
           thfactMU(1,-1)=ds1
           thfactMU(1,0)=dc1
!
!          D Order of d-orbitals is 3z^2-1, x^2-y^2, xz, xy, yz
!
           thfactMU(2,0)=3.0d0*dc1*dc1-1.0d0
           thfactMU(2,2)=ds1*ds1
           thfactMU(2,1)=ds1*dc1
           thfactMU(2,-2)=ds1*ds1
           thfactMU(2,-1)=ds1*dc1
!
!          ATOM B .............................................
!
!          S
           thfactNU(0,0)=1.0d0
!
!          P
!
           thfactNU(1,1)=ds2
           thfactNU(1,-1)=ds2
           thfactNU(1,0)=dc2
!
!          D Order of d-orbitals is 3z^2-1, x^2-y^2, xz, xy, yz
!
           thfactNU(2,0)=3.0d0*dc2*dc2-1.0d0
           thfactNU(2,2)=ds2*ds2
           thfactNU(2,1)=ds2*dc2
           thfactNU(2,-2)=ds2*ds2
           thfactNU(2,-1)=ds2*dc2
!
! --------------------------------------------------------------
!          Done with theta factors.
! --------------------------------------------------------------

               avgVmat = 0.0d0

!
! ========================================================
!          Do integral over phi:
! ========================================================

!          average over phi (divide by 2*pi)
! JOM : I add this normalization factor at the end, out
!       of the loops
!           averagephi=0.5d0/pi
! HAO : It should be averagephi = 1.d0/pi, since interal of phi is now
!       done in (0, pi).

!          The phi factors depend only on m.

           phifactor(0)=1.0d0

			iexc = 3
			call phiint(zr, r, ds, r1, r2, rna, ispecies, jspecies, kspecies, ispmin, ispmax, iexc, avgVmat)

           dsth=ds
! JOM I add this norm. factor outside loops
!           averagetheta=0.5d0
!           prod=wTheta(it)*averagetheta*dsth
           prod=wTheta(it)*dsth




           do inm=1,index_max
             nl1=pFdata_cell%l_Mu(inm)
             nm1=pFdata_cell%m_Mu(inm)
             nl2=pFdata_cell%l_nu(inm)
             nm2=pFdata_cell%m_nu(inm)
! JOM add psiAmat here

             stuffmunu=prod*thfactMU(nl1,nm1)*	&
     &              thfactNU(nl2,nm2)*psiBmat(inm)*psiAmat(inm)

               thetamat(:,inm)=thetamat(:,inm)+	&
     &                          avgVmat(:,inm)*stuffmunu

	end do !inm

!
!
!
	end do !ith



! ========================================================
	deallocate(avgVmat)
	end subroutine thetaint




	subroutine phiint(zr, r, ds,r1, r2, rna, ispecies, jspecies, kspecies, ispmin, ispmax, iexc, avgVmat)
    implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispecies, jspecies, kspecies, iexc, ispmin, ispmax
        real, intent (in) :: r, ds, r2, r1, rna(3), zr

! Output
		real, allocatable :: avgVmat(:,:)

! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
      integer, parameter :: numInthpP=25
      real phifactor(-3:3)
      integer ip, index, ix, ix1
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


! Procedure
! ===========================================================================




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

!
! ---------------------------------------------------
!
             do  iX=ispmin,ispmax
				pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies)

        		pFdata_cell=>pFdata_bundle%Fdata_cell_3c(2)

            	nME3c_max = pFdata_cell%nME

                IX1=IX+1
                vpot = dvxc3c (iexc, r, r2, r3, ispecies, jspecies, kspecies, IX1)

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


             end do !IX1
             deallocate (mleft)
             deallocate (mright)
		end do
! ===============================================
!          The end of the phi integral.

		end subroutine phiint

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

		drho = min(wf_species(in1)%dr_min, wf_species(in2)%dr_min, &
													wf_species(in3)%dr_min)


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


   end module
