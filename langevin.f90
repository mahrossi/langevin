! *********************************************************
! * An awkwardly simple code for the dynamics of an       *
! * an n-dimensional harmonic oscillator, to demonstrate  *
! * the use of colored-noise Langevin equation in         * 
! * thermostatting.                                       *   
! *                                                       *
! * Code has been kept as modular as possible, and it     *
! * should be relatively straightforward to adapt it for  *
! * other serial MD codes. In parallel codes using domain *
! * decomposition, the additional degrees of freedom      *
! * should "follow" the corresponding atoms.              *
! *                                                       *
! * Feel free to copy, modify, lend, borrow, nuke this    *
! * code, which is licensed under GPLv3 [www.gnu.org]     *
! * e-mail me at michele dot ceriotti at gmail dot com    *
! *********************************************************
! **********************************************************
! Modified in order to read an arbitrary numerical potential
! and only in 1D!!!!
! **********************************************************

program langevin
  use md_tools
  use md_gle
  use splines

  implicit none

  ! * options to be given in input file *
  ! seed              initial seed for PRNG
  ! temp              target temperature in atomic units
  ! dt                timestep
  ! nstep             number of steps to be performed
  ! stride            output every stride frames
  ! ndof              dimensionality of the oscillator
  ! wfile             name of the file to read the hessian 
  !                   decomposition from (see below)
  ! traj              logical, whether to output p,q trajectory
  ! potname           NEEDED: file name that contains the potential in 2 column format
  !                   with header in the form '# ndata'
  ! mass              mass of the particle in atomic units
  ! * thermostatting options * 
  ! wnw               "optimal" frequency for langevin (0. 0-> WN off) 
  ! glew              range shift freq. for GLE (0. 0-> GLE off) 
  ! * a note on thermostats: all the thermostats can be used at once,
  !   even if that makes not much sense. again, this is just a small
  !   test code
  ! * a note on units: mass is set to one, and k_B as well, so that
  !   at equilibrium <p^2>=omega^2 <q^2>=temp
  
  real*8 dum, dt, temp, nhw, wnw, glew, tau
  integer argc, seed, nstep, stride, nchains, nhmts, tstride, ndata, neq, ntraj
  character *256 fname, wfile, prefix, potname, dummy
  namelist /inp/ seed, wfile, dt, temp, nstep, stride, tstride, nchains, & 
                 nhmts, nhw, wnw, glew, potname, mass, neq, ntraj, tau

  real*8 :: q, p, f, wm(1,1), mass, dxsq
  real*8, allocatable :: potydata(:), potxdata(:), splinedy2(:), qt(:)
  real*8 v, k, h ! yes, it's all we need!
  real*8 dt2
#ifdef USELIBS
  integer :: ipiv(1)  ! needed to invert matrix
#endif
  integer irnd, istep, i, j, itraj, m, s

  ! reads command line
  argc = iargc()
  if ( argc .ne. 1 ) then
    write(6,*) '* Call me as: langevin <input>'
    stop
  endif  

  ! reads input file
  tstride=0
  mass=1
  glew=0.d0
  call getarg (1,fname)
  open(101,file=fname)
  read(101,inp)
  close (unit=101)
  irnd=-seed
  dt2=dt*0.5
  kt=temp
  write(*,*) 'Performing simulation with mass', mass
  write(*,*) 'Performing simulation with friction', 2.*wnw
  if(glew .gt. 0.d0) then
     write(*,*) 'Support for colored noise not available'
     stop
  endif

  allocate(qt(0:nstep/tstride))
!  allocate(wtw(ndof,ndof))
!  allocate(wm(ndof,ndof))
!#ifdef USELIBS
!  allocate(ipiv(ndof))
!#endif 

  ! reads potential file and creates the necessary
  ! arrays for interpolation

  OPEN(11, FILE=potname)
  READ(11,*) dummy, ndata
  IF (ALLOCATED(potydata)) DEALLOCATE(potydata)
  IF (ALLOCATED(potxdata)) DEALLOCATE(potxdata)
  IF (ALLOCATED(splinedy2)) DEALLOCATE(splinedy2)
  ALLOCATE(potydata(ndata))
  ALLOCATE(potxdata(ndata))
  ALLOCATE(splinedy2(ndata))

  DO i=1,ndata
    READ(11,*) potxdata(i), potydata(i)
  ENDDO
  CALL spline(potxdata, potydata, ndata, splinedy2)

  qmax=potxdata(ndata)
  ! init random seed
  dum=randu(irnd)

  ! init momenta
  p=rang(irnd)*sqrt(kt*mass)
  ! init coordinates 
  q=0.


  call force(q, potxdata, potydata, ndata, splinedy2, f)  


  ! perform a thermalization with an andersen thermostat
  ! for neq trajectories of same length as the main loop
  do itraj=1, neq
     p=rang(irnd)*sqrt(kt*mass)
     do istep=1,nstep

       !does some thermalization with the simplest possible andersen thermostat
       if (randu(irnd) < dt/tau) then 
         p=rang(irnd)*sqrt(kt*mass)
       endif

       ! hamiltonian step for dynamics
       p=p+f*dt2
       q=q+p*dt/mass
       call force(q, potxdata, potydata, ndata, splinedy2, f)
       p=p+f*dt2
     enddo
  enddo


  !initializes thermostats
  if(wnw .gt. 0.d0) call wn_init(dt2,wnw, mass)

  ! opens file for output
  if (tstride>0) open(102,file='traj-p.out') 
  if (tstride>0) open(103,file='traj-q.out') 
  open(104,file='statis.out') 
  open(105,file='msd.out')
  write(105,*) "# ", nstep/tstride, ntraj

  ! we are already at the main dynamics loop!
  do itraj=1, ntraj
     qt=0.0
     p=rang(irnd)*sqrt(kt*mass)
     langham=0.d0
     do istep=1,nstep

     !calls the active thermostats
     if (wnw .gt. 0.d0 .or. glew .gt. 0.d0) then 
        call kin(p, k, mass)
        langham=langham+k
        if (wnw .gt. 0.d0)  call wn_step(p,irnd)
        call kin(p, k, mass)
        langham=langham-k
      endif

    ! hamiltonian step for dynamics
      p=p+f*dt2
      q=q+p*dt/mass
      call force(q, potxdata, potydata, ndata, splinedy2, f)
      p=p+f*dt2

      ! thermostats, second bit
      if (wnw .gt. 0.d0) then 
        call kin(p, k, mass)
        langham=langham+k
        if (wnw .gt. 0.d0)  call wn_step(p,irnd)
        call kin(p, k, mass)
        langham=langham-k
      endif

    ! computes properties & outputs
      if (mod(istep,stride).eq.0) then
         call pot(q, potxdata, potydata, ndata, splinedy2, v)
         call kin(p, k, mass)
         h=k+v

         if (wnw .gt. 0.d0 .or. glew .gt. 0.d0) then
           h=h+langham
         end if
         if (tstride.gt.0 .and. mod(istep,tstride).eq.0) then
           qt((istep/tstride)-1)=q
           write(102,'(1e17.8 )', advance='NO') istep*dt
           write(103,'(1e17.8 )', advance='NO') istep*dt
           write(102,'( 1e17.8 )', advance='NO') p
           write(103,'( 1e17.8 )', advance='NO') q
           write(102,*) ""
           write(103,*) ""
         endif
         write(104,'(4e17.8)') istep*dt, v, k, h
       endif
     enddo
     ! now calculate MSD

      m = nstep/tstride-1
      do s = 0,m
        dxsq = 0.d0
        do j = 0,m-s
          dxsq = dxsq+(qt(j+s)-qt(j))**2 
        enddo
        dxsq = dxsq/(m-s+1)
        write(105,*) s*tstride*dt,dxsq
      enddo


     write(102, *)" "
     write(103, *)" "
     write(104, *)" "
     write(105, *)" "
  enddo
  if (tstride>0) close(102)
  if (tstride>0) close(103)
  close(104)
  


end program langevin

