module m_neighbour

  implicit none
  integer, parameter :: rp = kind(1.0)

  integer, parameter :: &
       first_alloc = 100, & ! allocate to first_alloc*nat
       batch_size = 20 ! if lists exceed first_alloc, increase in steps of batch_size*nat


  !> single bin
  type t_bin
     integer :: tag
     integer, private :: batchsize = 200
     integer :: n_members=0
     integer :: n_tot=0
     integer, allocatable :: global_idx_of_member(:)
     real(rp), allocatable :: pos_member(:,:)
     real( rp ) :: center(3)
   contains
     procedure :: bin_add
  end type t_bin


  !> neigbour object
  type t_neighbour

     ! total size of arrays
     integer, private :: ntot = 0
     ! current head_idx
     integer, private :: head_idx = 0


     !> flag to indicate if neighbor list is computed or not
     integer :: active = -1

     !> distance cutoff value
     real(rp) :: rcut

     !> number of neighbors of atom i (excluding self)
     integer, pointer, contiguous :: nneig(:)

     !> concatenated list of neigbors (excluding idx)
     integer, pointer, contiguous :: neiglist(:)

     !> concatenated list of vectors of neiglist in pbc (vector for idx is excluded,
     !! but assumed to be placed at [0.0, 0.0, 0.0] )
     real(rp), pointer, contiguous :: veclist(:,:)

     integer, allocatable, private :: partial_sumlist(:)

     integer :: bin_tag=1

     !> list of bins
     type( t_bin ), allocatable :: bins(:,:,:)
   contains
     procedure :: get_list
     procedure :: get_veclist
     procedure, private :: add
     final :: t_neighbour_destroy
  end type t_neighbour

  interface t_neighbour
     procedure :: compute_binned_pbc
  end interface t_neighbour



contains

  subroutine t_neighbour_destroy( self )
    type( t_neighbour ), intent(inout) :: self
    if( associated( self% nneig) )deallocate( self% nneig )
    if( associated( self% neiglist) )deallocate( self% neiglist )
    if( associated( self% veclist) )deallocate( self% veclist )
    if( allocated( self% partial_sumlist))deallocate( self% partial_sumlist)
  end subroutine t_neighbour_destroy



  function get_list( self, idx, list )result(n)
    !! Output the `list` of indices from neiglist, which are neighbours of `idx`.
    !! Return `n` which i sthe size of `list`
    !! The index `idx` is NOT included in the returned array `list`.
    class( t_neighbour ), intent(inout) :: self

    !> input atomic index
    integer, intent(in)                 :: idx

    !> output list of neighbours to atom `idx`
    integer, allocatable, intent(out)   :: list(:)

    !> `n`, size of returned `list` is `(n)`; if `idx` is invalid `n=-1`
    integer :: n

    integer :: i_start, i_end

    if( idx .le. 0 .or. idx .gt. size(self% nneig) ) then
       ! idx out of range
       n = -1
       return
    end if

    ! The list corresponding to idx is a section of total
    ! neiglist, between indices: partial_sum(idx-1) and partial_sum(idx)

    ! find indices where the list for idx starts and ends
    i_end = self% partial_sumlist(idx)
    i_start = i_end - self% nneig(idx)+1

    n = i_end - i_start + 1
    allocate( list(1:n), source=self% neiglist(i_start:i_end))

  end function get_list



  function get_veclist( self, idx, veclist )result(n)
    !! Output vectors from `veclist`, which are neighbors of `idx`.
    !! Return `n` which is the size of `veclist`
    !! The vector of atom `idx` is NOT included, but is placed at (0.0, 0.0, 0.0)
    class( t_neighbour ), intent(inout) :: self

    !> input atomic index
    integer, intent(in)                 :: idx

    !> output list of vectors neighbour to `idx`
    real(rp), allocatable               :: veclist(:,:)

    !> `n`, size of returned `veclist` is `(3, n)`; if `idx` is invalid `n=-1`
    integer :: n

    integer :: i_start, i_end

    if( idx .le. 0 .or. idx .gt. size(self% nneig) ) then
       ! idx out of range
       n = -1
       return
    end if

    ! find indices where the list for idx starts and ends
    i_end = self% partial_sumlist(idx)
    i_start = i_end - self% nneig(idx)+1

    n = i_end - i_start + 1
    allocate( veclist, source=self% veclist(1:3, i_start:i_end) )

  end function get_veclist




  function compute_binned_pbc( nat, ityp, pos, lat, rcut )result(self)
    !! Compute the neighbour list in pbc.
    !! The idea is to actually detect the boundaries of the box, and set properly shifted
    !! fake atoms to the other side of the box, and do this on all sides (hi, and lo).
    !! Then we compute the neighbours by binning, and without pbc.
    !! The edges of the box are detected by transforming the cutoff `rcut` into cristalline coords,
    !! and then each atom that is within that distance of any box boundary is considered as box edge,
    !! and is replicated on the other side of the box (in all applicable directions).
    !! The bins are computed in reciprocal coords.
    !!
    !! call as:
    !!
    !!```f90
    !!  type( t_neighbour ), pointer :: neigh
    !!    ...
    !!  neigh => t_neighbour( nat, ityp, pos, lat, rcut )
    !!    ...
    !!  deallocate( neigh )
    !!```
    implicit none
    !> number of atoms
    integer, intent(in) :: nat

    !> integer atomic types (actually unused)
    integer, intent(in) :: ityp(nat)

    !> atomic positions, shape(3,nat)
    real(rp), intent(in) :: pos(3,nat)

    !> lattice vectors in rows, in units of `pos`
    real(rp), intent(in) :: lat(3,3)

    !> distance cutoff, in units of `pos`
    real(rp), intent(in) :: rcut

    !> t_neighbour class
    type( t_neighbour ), pointer :: self

    integer :: i, j, ii, jj, kk, ll, mm
    real(rp) :: dij, rij(3)
    real(rp) :: invlat(3,3)
    real(rp) :: r(3), rf(3), kadd(3)
    integer, dimension(3) :: ibin, addbin
    integer :: xbin, ybin, zbin
    real(rp) :: rcut2

    integer :: nbins, nbins_check
    real(rp) :: kvec(3), kmax, lo_buffer, hi_buffer
    integer :: jj_lb, kk_lb, ll_lb
    integer :: jj_ub, kk_ub, ll_ub


    ! cutoff square
    rcut2 = rcut*rcut

    allocate( t_neighbour::self )

    ! inverse lattice
    call inverse3x3( lat, invlat)

    ! maximal "displacement vector" in cartesian coords ...
    ! this is an overshoot proportional to sqrt(3.0)*rcut in a cubic box,
    ! but no idea in general.
    ! Is probably not optimal, but good-enough for now. Boxes with weird shapes
    ! could be problematic from this point of view (e.g. some ax similar norm
    ! to rcut, while others much larger).
    kvec(:) = rcut
    ! see how this vector transforms to crist coords
    call cart_to_crist(kvec, lat, invlat)
    ! take the maxval as unified cutoff distance in crist coords over all axes
    ! (this should maybe resolve to i,j,k axes at some point?)
    kmax = maxval(abs(kvec))

    ! write(*,*) "kvec:", kvec
    ! write(*,*) 1.0/kvec(1)
    ! write(*,*) "kmax",kmax

    ! take the floor() for nbins, since fewer bins means each bin is larger in size
    nbins = floor(1.0_rp/kmax)
    ! protection against zero bins
    nbins = max(1,nbins)
    ! nbins = nint(1.0_rp/kmax)
    ! nbins = 64
    ! nbins_check = ceiling(1.0_rp/kmax)

    allocate( self% nneig(1:nat), source=0)
    allocate( self% neiglist(1:first_alloc*nat) )
    allocate( self% veclist(1:3,1:first_alloc*nat) )
    allocate( self% partial_sumlist(1:nat) )
    self% partial_sumlist(1) = 0
    ! save first ntot (alloc size)
    self% ntot = first_alloc*nat


    ! buffer region near low-end
    ! lo_buffer = -0.5+kmax
    lo_buffer = kmax
    ! buffer region near high-end
    ! hi_buffer = 0.5-kmax
    hi_buffer = 1.0_rp-kmax

    ! write(*,*) "lo_buffer =",lo_buffer
    ! write(*,*) "hi_buffer =",hi_buffer
    ! write(*,*) "nbins=",nbins

    ! create bins
    allocate( self% bins(0:nbins+1, 0:nbins+1, 0:nbins+1) )
    do i = 1, nat
       r = pos(:,i)

       ! atom needs to be inside box
       call cart_to_crist( r, lat, invlat )
       ! shift to -0.5:0.5
       call periodic(r)
       ! shift to 0:1
       r = r + [0.5_rp, 0.5_rp, 0.5_rp ]

       ! write(*,*)
       ! write(*,*) " ==> idx i:",i
       ! write(*,*) r
       ! write(*,*) r*4

       ! compute bin of this atom in reciprocal
       ! ibin goes from 1 to nbins
       ibin = [floor(r*nbins)] + [1,1,1]

       ! write(*,"(a,1x,3i4, 3f9.4)") "ibin:",ibin, r

       ! get the cartesian coords of this atom non-shifted
       rf = r
       call crist_to_cart( rf, lat, invlat )
       ! add to bin
       call self% bins( ibin(1), ibin(2), ibin(3) )% bin_add( 1, i, rf )

       ! detect boundaries of the box:
       ! if atom is add the boundary, add a fake atom with identical index to the other side
       ! of the box. The fake atom has its coords shifted by appropriate vector.

       ! single components (boundary planes)
       do ii = 1, 3
          addbin = ibin
          kadd = 0.0_rp
          if( check1( r(ii), kmax, kadd(ii)) ) then
             ! need to add shifted atm
             rf = r + kadd

             ! which bin to add? assume low-end, if kadd=1.0 then it's high-end
             addbin(ii) = 1
             if( kadd(ii) > 0.5_rp ) addbin(ii) = nbins

             call crist_to_cart( rf, lat, invlat )
             call self% bins( addbin(1), addbin(2), addbin(3))% bin_add( 11, i, rf )
          end if
       end do

       ! two components simultaneously (boundary edges)
       do ii = 1, 3
          do jj = ii+1, 3
             addbin = ibin
             kadd = 0.0_rp
             if( check2( r, ii, jj, kmax, kadd ) )then
                ! which bin: high-end where kadd=1.0; low-end where kadd=-1.0
                addbin = merge( nbins, addbin, kadd > 0.5_rp )
                addbin = merge( 1, addbin, kadd < -0.5_rp )
                rf = r + kadd
                call crist_to_cart( rf, lat, invlat )
                call self% bins( addbin(1), addbin(2), addbin(3))% bin_add( 12, i, rf )
             end if
          end do
       end do

       ! three components simultaneously (boundary corners)
       do ii = 1, 3
          do jj = ii+1, 3
             do kk = jj+1, 3
                addbin = ibin
                kadd = 0.0_rp
                if( check3( r, ii, jj, kk, kmax, kadd) ) then
                   ! which bin: high-end where kadd=1.0; low-end where kadd=-1.0
                   addbin = merge( nbins, addbin, kadd > 0.5_rp )
                   addbin = merge( 1, addbin, kadd < -0.5_rp )
                   rf = r + kadd
                   call crist_to_cart( rf, lat, invlat )
                   call self% bins( addbin(1), addbin(2), addbin(3))% bin_add( 13, i, rf )
                end if
             end do
          end do
       end do

    end do

    ! write(*,*) "start"


    ! head_idx is the last index in neiglist
    self% head_idx = 0
    do i = 1, nat

       ! write(*,*) "&&& ", i
       ! write(*,"(a,1x,3f9.4)") "pos",pos(:,i)

       r = pos(:,i)

       ! shift inside box; range 0:1 in crist
       call cart_to_crist( r, lat, invlat )
       call periodic( r )
       r = r + [0.5_rp, 0.5_rp, 0.5_rp ]

       ! find my bin in reciprocal
       ibin = [floor(r*nbins)] + [1,1,1]

       ! write(*,*) r
       ! write(*,"(a,1x,3i4,3f9.4)")"ibin",ibin,r

       ! r into cartesian
       call crist_to_cart( r, lat, invlat )

       xbin = ibin(1)
       ybin = ibin(2)
       zbin = ibin(3)

       !
       ! loop over the bins:
       ! current ibin, and +-1 in each direction (cube around ibin)
       !

       ! lower bound
       jj_lb = -1
       kk_lb = -1
       ll_lb = -1

       ! upper bound
       jj_ub = 1
       kk_ub = 1
       ll_ub = 1

       ! if we are in bin=1 (low-edge), no need to check the -1 bin
       ! if( ibin(1) == 1 ) jj_lb = 0
       ! if( ibin(2) == 1 ) kk_lb = 0
       ! if( ibin(3) == 1 ) ll_lb = 0

       ! if we are in the bin=nbins (high-edge), no need to check the +1 bin
       ! if( ibin(1) == nbins ) jj_ub = 0
       ! if( ibin(2) == nbins ) kk_ub = 0
       ! if( ibin(3) == nbins ) ll_ub = 0

       ! write(*,*) i, "ibin",ibin
       ! write(*,*) "jj_lb", jj_lb, "jj_ub", jj_ub
       ! write(*,*) "kk_lb", kk_lb, "kk_ub", kk_ub
       ! write(*,*) "ll_lb", ll_lb, "ll_ub", ll_ub
       do jj = jj_lb, jj_ub
          do kk = kk_lb, kk_ub
             do ll = ll_lb, ll_ub

                ! set the bin where we take atoms
                xbin = ibin(1) + jj
                ybin = ibin(2) + kk
                zbin = ibin(3) + ll

                ! loop over the atoms of that bin
                ! write(*,"(a,1x,3i4,3x,i0)") " >> ", xbin, ybin, zbin, self% bins(xbin,ybin,zbin)%n_members
                do ii = 1, self% bins( xbin, ybin, zbin )% n_members

                   ! vector between each member of the bin and r(i)
                   rij = self% bins( xbin, ybin, zbin )% pos_member(:,ii) - r

                   ! distance square from i to j
                   dij = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

                   ! write(*,'(i16,1x,3f9.4,5x,f9.4)') ii, self%bins(xbin,ybin,zbin)%pos_member(:,ii), sqrt(dij)
                   if( dij .le. rcut2 ) then
                      ! get the true index of the atom j, and add it to neighbour list of i
                      j = self% bins( xbin, ybin, zbin )% global_idx_of_member(ii)
                      call self% add(i, j, rij )
                   end if

                end do

             end do
          end do
       end do

       ! add current nneig(i) to partial sum list
       if( i .eq. 1) then
          self% partial_sumlist(1) = self% nneig(1)
       else
          self% partial_sumlist(i) = self% partial_sumlist(i-1) + self% nneig(i)
       end if

    end do

  end function compute_binned_pbc




  subroutine bin_add( self, bin_tag, idx, vec )
    !! add atom `idx` and vector `vec` to this bin
    class( t_bin ), intent(inout) :: self
    integer, intent(in) :: bin_tag
    integer, intent(in) :: idx
    real(rp), intent(in) :: vec(3)
    integer, allocatable :: tmp(:)
    real(rp), allocatable :: tmp_v(:,:)
    integer :: oldsize, newsize

    ! write(*,*) "   enter bin_add. n_tot:",self% n_tot, self% n_members, self% tag
    ! write(*,*) bin_tag, vec, idx
    self% tag = bin_tag
    ! if( .not.allocated( self% global_idx_of_member) ) then
    if( self% n_tot == 0 ) then
       allocate( self% pos_member(1:3,1:self%batchsize))
       allocate( self% global_idx_of_member(1:self%batchsize))
       self% n_tot = self%batchsize
       self% n_members = 0
       ! bin_tag = bin_tag + 1
       ! write(*,*) "   allocated new bin. n_tot:", self% n_tot
    end if

    ! realloc
    if( self% n_members + 1 .gt. self% n_tot ) then
       oldsize = self% n_tot
       newsize = oldsize + self%batchsize
       ! write(*,*) " $$$ bin realloc", oldsize, newsize
       ! write(*,*) "  batchsize:",self% batchsize
       allocate( tmp, source=self%global_idx_of_member )
       deallocate( self% global_idx_of_member )
       allocate( self%global_idx_of_member(1:newsize) )
       self%global_idx_of_member(1:oldsize) = tmp(:)
       deallocate(tmp)

       allocate( tmp_v, source=self%pos_member )
       deallocate( self% pos_member )
       allocate( self% pos_member(1:3,1:newsize) )
       self% pos_member(1:3,1:oldsize) = tmp_v(:,:)
       deallocate( tmp_v )

       self% n_tot = newsize
    end if
    self% n_members = self% n_members + 1
    self% global_idx_of_member( self% n_members ) = idx
    self% pos_member(:, self% n_members ) = vec

    ! write(*,*) "exit bin_add", self% n_members, self% n_tot
    ! write(*,*) "added pos", self% pos_member(:, self% n_members)
  end subroutine bin_add


  subroutine add( self, origin_idx, idx, vec )
    !! add index `idx` and vector `vec` to atom `origin_idx`
    implicit none
    class( t_neighbour ), intent(inout) :: self
    integer, intent(in)  :: origin_idx   !! index of origin atom
    integer, intent(in)  :: idx          !! index to add to neiglist
    real(rp), intent(in) :: vec(3)       !! vec to add to veclist

    integer :: newsize, oldsize
    integer, allocatable :: tmp(:)
    real(rp), allocatable :: tmp_v(:,:)

    integer :: i_s

    if( idx .eq. origin_idx ) return
    ! starting index
    i_s = 1
    if( origin_idx > 1) i_s = self% partial_sumlist(origin_idx-1) + 1
    ! atom origin_idx already has neiglist:
    ! write(*,"(*(i5,:,1x))") self% neiglist(i_s:self% head_idx )
    if( any(self% neiglist(i_s:self% head_idx) .eq. idx) )then
       ! this can happen when box is small compared to rcut, we get the same atom from
       ! the periodic image as its own neighbor, is actually ok, no need to return

       ! write(*,*) "idx=", idx, "already in:"
       ! write(*,"(2x,8(i8,:,2x))") self% neiglist(i_s:self%head_idx)
       ! write(*,*) "vec:",vec
       ! return
    end if


    ! increase nneig of origin_idx
    self% nneig( origin_idx ) = self% nneig( origin_idx ) + 1

    ! increase the head_idx
    self% head_idx = self% head_idx + 1

    ! check realloc
    if( self% head_idx .gt. self% ntot ) then
       oldsize = self% ntot
       newsize = oldsize + batch_size*size(self%nneig)
       ! write(*,*) "calling realloc", oldsize, newsize
       ! move alloc to tmp arrays, then
       ! allocate newsize and copy old data from tmp
       allocate( tmp, source=self% neiglist )
       deallocate( self% neiglist )
       ! call move_alloc( self% neiglist, tmp )
       allocate( self% neiglist(1:newsize) )
       self% neiglist(1:oldsize) = tmp(:)
       deallocate( tmp )
       !
       allocate( tmp_v, source=self% veclist )
       deallocate( self% veclist )
       ! call move_alloc( self% veclist, tmp_v )
       allocate( self% veclist(1:3, 1:newsize) )
       self% veclist(1:3,1:oldsize) = tmp_v(:,:)
       deallocate( tmp_v )
       ! increment ntot
       self% ntot = newsize
    end if

    ! add idx to head of neiglist
    self% neiglist( self% head_idx ) = idx

    ! add vec to head of veclist
    self% veclist( 1:3, self% head_idx ) = vec(1:3)

  end subroutine add



  pure subroutine periodic(c)
    !--------------------------------
    ! periodic boundary condition, for 3 dimensional vector input in crist coords.
    !--------------------------------
    implicit none
    real(RP), dimension(3),intent(inout) :: c
    integer :: i

    ! do i = 1, 3
    !    if( c(i) .lt. -0.5 ) c(i) = c(i) + 1.0
    !    if( c(i) .ge. 0.5 ) c(i) = c(i) - 1.0
    ! end do

    c(1) = c(1) - int( (c(1)/ 0.5_rp) )
    c(2) = c(2) - int( (c(2)/ 0.5_rp) )
    c(3) = c(3) - int( (c(3)/ 0.5_rp) )
  end subroutine periodic


  pure subroutine crist_to_cart( rij, lat, invlat )
    ! lat is lattice vectors in rows:
    !
    !  lat = a1 a2 a3
    !        b1 b2 b3
    !        c1 c2 c3
    !
    ! invlat is inverse of lat
    implicit none
    real(rp), intent(inout) :: rij(3)
    real(rp), intent(in) :: lat(3,3)
    real(rp), intent(in) :: invlat(3,3)

    rij = matmul( lat, rij )
  end subroutine crist_to_cart

  pure subroutine cart_to_crist( rij, lat, invlat )
    ! lat is lattice vectors in rows:
    !
    !  lat = a1 a2 a3
    !        b1 b2 b3
    !        c1 c2 c3
    !
    ! invlat is inverse of lat
    implicit none
    real(rp), intent(inout) :: rij(3)
    real(rp), intent(in) :: lat(3,3)
    real(rp), intent(in) :: invlat(3,3)

    rij = matmul( invlat, rij )
  end subroutine cart_to_crist


  pure subroutine inverse3x3( mat, inv )
    implicit none
    real(rp), intent(in) :: mat(3,3)
    real(rp), intent(out) :: inv(3,3)

    real(rp) :: det, invdet

    det = 0.0_rp

    ! calculate the determinant
    det = det + mat(1,1)*mat(2,2)*mat(3,3) &
         + mat(1,2)*mat(2,3)*mat(3,1) &
         + mat(1,3)*mat(2,1)*mat(3,2) &
         - mat(1,3)*mat(2,2)*mat(3,1) &
         - mat(1,2)*mat(2,1)*mat(3,3) &
         - mat(1,1)*mat(2,3)*mat(3,2)
    ! invert the determinant
    invdet = 1/det
    ! calculate the inverse matrix
    inv(1,1) = invdet  * ( mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2) )
    inv(2,1) = -invdet * ( mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1) )
    inv(3,1) = invdet  * ( mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1) )
    inv(1,2) = -invdet * ( mat(1,2)*mat(3,3) - mat(1,3)*mat(3,2) )
    inv(2,2) = invdet  * ( mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1) )
    inv(3,2) = -invdet * ( mat(1,1)*mat(3,2) - mat(1,2)*mat(3,1) )
    inv(1,3) = invdet  * ( mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2) )
    inv(2,3) = -invdet * ( mat(1,1)*mat(2,3) - mat(1,3)*mat(2,1) )
    inv(3,3) = invdet  * ( mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1) )
  end subroutine inverse3x3

  function check1( r, kmax, kadd ) result( add )
    ! check the i component of r
    ! If r(i) is below the low_buffer(=kmax), then we need to add 1.0 on axis i;
    ! if r(i) is above the high_buffer(=1.0-kmax) then we need to add -1.0 on axis i.
    ! return .true. when `kadd` vector has any component different from 0.0
    real(rp), intent(in) :: r
    real(rp), intent(in) :: kmax
    real(rp), intent(inout) :: kadd
    logical :: add
    add = .true.
    if( r .le. kmax ) then
       kadd = 1.0_rp
       return
    elseif( r .gt. 1.0-kmax ) then
       kadd = -1.0_rp
       return
    end if
    kadd = 0.0_rp
    add = .false.
  end function check1

  function check2( r, i, j, kmax, kadd ) result(add)
    ! check components i and j of r simultaneously
    real(rp), intent(in) :: r(3)
    integer, intent(in) :: i, j
    real(rp), intent(in) :: kmax
    real(rp), intent(out) :: kadd(3)
    logical :: add
    logical :: ai, aj
    real(rp) :: k_cp(3)
    add = .false.
    kadd = 0.0_rp
    k_cp = kadd
    ai = check1( r(i), kmax, k_cp(i) )
    aj = check1( r(j), kmax, k_cp(j) )
    if( ai .and. aj ) then
       add = .true.
       kadd = k_cp
    end if
  end function check2

  function check3( r, i, j, k, kmax, kadd )result(add)
    ! check components i, j, k of r simultaneously
    real(rp), intent(in) :: r(3)
    integer, intent(in) :: i, j, k
    real(rp), intent(in) :: kmax
    real(rp), intent(inout) :: kadd(3)
    logical :: add
    logical :: ai, aj, ak
    real(rp) :: k_cp(3)
    add = .false.
    kadd = 0.0_rp
    k_cp = kadd
    ai = check1( r(i), kmax, k_cp(i) )
    aj = check1( r(j), kmax, k_cp(j) )
    ak = check1( r(k), kmax, k_cp(k) )
    if( ai .and. aj .and. ak )then
       add = .true.
       kadd = k_cp
    end if
  end function check3


end module m_neighbour
