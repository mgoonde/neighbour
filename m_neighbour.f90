module m_neighbour

  implicit none

  integer, parameter :: rp = kind(1.0)    ! modify precision as needed

  private
  public :: rp
  public :: t_neighbour

  integer, parameter :: &
       bin_batchsize = 200, &
       first_alloc = 100, & ! allocate to first_alloc*nat
       batch_size = 20 ! if lists exceed first_alloc, increase in steps of batch_size*nat


  !> neigbour object
  type, public :: t_neighbour

     ! total size of arrays
     integer, private :: ntot = 0
     ! current head_idx
     integer, private :: head_idx = 0
     ! dimension of input
     ! integer, private :: ndim = 3

     !> number of neighbors of atom `i` (excluding self)
     integer, pointer, contiguous :: nneig(:)

     !> concatenated list of neigbors (excluding `idx`)
     integer, pointer, contiguous :: neiglist(:)

     !> concatenated list of vectors of neiglist in pbc (vector for idx is exc
     !! but assumed to be placed at `[0.0, 0.0, 0.0]` )
     real(rp), pointer, contiguous :: veclist(:,:)

     !> partial sums of `nneig(:)`, used for retrieving the slices of
     !! arrays `neiglist(:)` and `veclist(:,:)` which correspond to a certain atom
     integer, allocatable :: partial_sumlist(:)

     !> list of bins
     type( t_bin ), allocatable, private :: bins(:,:,:)


     !> |unused| flag to indicate if neighbor list is computed or not
     integer, private :: active = -1

     !> |unused| distance cutoff value
     real(rp), private :: rcut

     !> |unused|
     integer, private :: bin_tag=1
   contains
     procedure :: get_nn
     procedure :: get
     procedure :: expand
     procedure, private :: add
     final :: t_neighbour_destroy
  end type t_neighbour

  !> overload the name `t_neighbour` with the `compute_binned_pbc` functions.
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
  interface t_neighbour
     procedure :: compute_binned_pbc
     procedure :: compute_binned_pbc_one
  end interface t_neighbour



  !> single bin
  type, private :: t_bin
     ! not needed
     integer :: tag=0

     ! (re-)allocation size
     integer, private :: batchsize = bin_batchsize

     ! current number of atoms in bin
     integer :: n_members=0

     ! current allocation size
     integer :: n_tot=0

     ! indices of bin members
     integer, allocatable :: global_idx_of_member(:)

     ! position vectors of bin members
     real(rp), allocatable :: pos_member(:,:)

     ! unused
     real( rp ) :: center(3)
   contains
     procedure :: bin_add
  end type t_bin



contains

  subroutine t_neighbour_destroy( self )
    type( t_neighbour ), intent(inout) :: self
    if( associated( self% nneig) )deallocate( self% nneig )
    if( associated( self% neiglist) )deallocate( self% neiglist )
    if( associated( self% veclist) )deallocate( self% veclist )
    if( allocated( self% partial_sumlist))deallocate( self% partial_sumlist)
  end subroutine t_neighbour_destroy




  function get_nn( self, idx, list, veclist )result(n)
    !! Get the first neighbor shell of `idx` from neighbour list.
    !! Return `n` which is the number of neighbors.
    !! The vector of atom `idx` is NOT included in output.
    class( t_neighbour ), intent(inout) :: self

    !> input atomic index
    integer, intent(in)                 :: idx

    !> output list of neighbours to atom `idx`
    integer, allocatable, intent(out), optional   :: list(:)

    !> output list of vectors neighbour to `idx`, which is assumed at [0.0, 0.0, 0.0]
    real(rp), allocatable, intent(out), optional  :: veclist(:,:)

    !> `n`, number of first neighbors; if `idx` is invalid, or the
    !! neighbour list has not been computed `n=-1`
    integer :: n

    integer :: i_start, i_end

    if(  idx .le. 0 .or. &
         idx .gt. size(self% nneig) .or. &
         .not.allocated(self%partial_sumlist)) then
       ! idx out of range, or neiglist not computed
       n = -1
       return
    end if

    ! find indices where the list for idx starts and ends
    i_end = self% partial_sumlist(idx)
    i_start = i_end - self% nneig(idx)+1

    n = i_end - i_start + 1
    if( present(list)) allocate( list(1:n), source=self% neiglist(i_start:i_end))
    if(present(veclist)) allocate( veclist, source=self% veclist(:, i_start:i_end) )

  end function get_nn


  function get( self, idx, list, veclist, nbond )result(n)
    !! Get the neighbor data of `idx` from neighbour list, up to `nbond`.
    !! Return `n` which is the number of neighbors.
    !! The vector of atom `idx` is NOT included in output.
    class( t_neighbour ), intent(inout) :: self

    !> input atomic index
    integer, intent(in)                 :: idx

    !> output list of neighbours to atom `idx`
    integer, allocatable, intent(out), optional  :: list(:)

    !> output list of vectors neighbour to `idx`
    real(rp), allocatable, intent(out), optional :: veclist(:,:)

    !> how many bond shells to get (default=1). If `nbond>1`, the output list is not sorted
    integer, intent(in), optional :: nbond

    !> `n`, number of first neighbors; if `idx` is invalid, or the
    !! neighbour list has not been computed `n=-1`
    integer :: n

    integer :: nb, i, ndim
    integer, allocatable :: inlist(:)
    real(rp), allocatable :: v_inlist(:,:)

    if(  idx .le. 0 .or. &
         idx .gt. size(self% nneig) .or. &
         .not.allocated(self%partial_sumlist)) then
       ! idx out of range, or neiglist not computed
       n = -1
       return
    end if

    nb = 1
    if(present(nbond))nb=nbond
    if( nb < 1 ) then
       n = 0
       return
    end if


    ! create array with just idx
    allocate(inlist(1:1), source=idx )
    if(present(veclist)) then
       ndim = size(veclist,1)
       allocate(v_inlist(1:ndim,1:1) )
       v_inlist(:,:) = 0.0_rp
    end if

    ! expand it by nb
    if( present(veclist)) then
       n = self% expand( nb, inlist, veclist=v_inlist )
    else
       n = self% expand( nb, inlist )
    end if

    ! check error
    if( n .lt. 0 ) return

    ! first index of `inlist` is idx, remove it for output
    if(present(veclist)) allocate( veclist(1:ndim,1:n-1))
    if( present(list)) allocate(list(1:n-1))
    do i = 1, n-1
       if( present(veclist)) veclist(:,i) = v_inlist(:,i+1)
       if( present(list))  list(i) = inlist(i+1)
    end do
    n = n - 1

  end function get


  function expand( self, nbond, list, veclist )result(n)
    !! expand the list in input by `nbond` number of bond shells.
    implicit none
    class( t_neighbour ), intent(inout) :: self

    !> number of bons shells to expand
    integer, intent(in) :: nbond

    !> on input: list to be expanded; on output: expanded list
    integer, allocatable, intent(inout) :: list(:)

    !> veclist on input; on output modified to include the expansion vectors
    real(rp), allocatable, intent(inout), optional :: veclist(:,:)

    !> number of elements in output list
    integer :: n

    integer :: i, ntot, ncur, idx
    integer, allocatable :: work(:), tmp(:)
    real(rp), allocatable :: vwork(:,:), vtmp(:,:)
    integer, parameter :: batchsize=50
    integer :: ibond, ii, nl_idx, nn
    integer, allocatable :: l_idx(:)
    real(rp), allocatable :: vl_idx(:,:)
    real(rp), allocatable :: origin(:)

    if(.not.allocated(self%partial_sumlist)) then
       ! neiglist not computed
       n = -1
       return
    end if

    if( nbond == 0) then
       n=0
       return
    end if


    n = size(list)

    ntot = n+batchsize
    ! allocate larger list for work
    allocate( work(1:ntot))
    work(1:n) = list
    ! vecs
    if(present(veclist)) then
       allocate( vwork(1:size(veclist,1),1:ntot) )
       vwork(:,1:n) = veclist
       ! origin vec
       allocate(origin(1:size(veclist,1)))
    end if


    ! current size
    ncur = n

    do ibond = 1, nbond
       ! go over current work
       nn = ncur
       do i = 1, nn
          ! this idx
          idx = work(i)
          ! neighlist of this idx
          nl_idx = self% get_nn( idx, list=l_idx )
          ! vec
          if( present(veclist)) then
             nl_idx = self% get_nn( idx, veclist=vl_idx )
             origin = vwork(:,i)
          end if

          if( nl_idx < 0 ) then
             ! error
             n = -1
             return
          end if

          ! add to work
          do ii = 1, nl_idx
             ! skip if already there
             if( any(work(1:ncur) .eq. l_idx(ii)) )cycle
             if( ncur+1 .gt. ntot ) then
                ! realloc
                ntot = ntot+batchsize
                call move_alloc( work, tmp )
                allocate( work(1:ntot) )
                work(1:ncur) = tmp(:)
                deallocate(tmp)
                ! vec
                if( present(veclist)) then
                   call move_alloc( vwork, vtmp )
                   allocate( vwork(1:size(veclist,1),1:ntot) )
                   vwork(:, 1:ncur) = vtmp(:,:)
                   deallocate(vtmp)
                end if
             end if
             ncur = ncur + 1
             work(ncur) = l_idx(ii)
             ! vec
             if( present(veclist)) then
                vwork(:,ncur) = vl_idx(:,ii) + origin
             end if

          end do
          deallocate( l_idx )
          if(present(veclist)) deallocate( vl_idx )
       end do
    end do

    ! set output
    n = ncur
    deallocate(list)
    allocate(list(1:ncur))
    list(:) = work(1:ncur)

    if( present(veclist)) then
       deallocate(veclist)
       allocate(veclist(1:size(vwork,1),1:ncur))
       veclist(:,:) = vwork(:,1:ncur)
       deallocate(origin)
    end if

  end function expand


  function compute_binned_pbc_one( nat, ityp, pos, lat, rcut, sort_by ) result(self)
    !! version where `rcut` is a single value, meaning the same `rcut` for all `ityp` values.
    implicit none
    !> number of atoms
    integer, intent(in) :: nat

    !> integer atomic types
    integer, intent(in) :: ityp(nat)

    !> atomic positions, shape(ndim, nat)
    real(rp), intent(in) :: pos(:,:)

    !> lattice vectors in rows, in units of `pos`
    real(rp), intent(in) :: lat(size(pos,1), size(pos,1))

    !> distance cutoff, in units of `pos`
    real(rp), intent(in) :: rcut

    !> optional string to sort the neighbour list; possible values `"distance"`, `"index"`
    character(*), intent(in), optional :: sort_by

    !> t_neighbour class
    type( t_neighbour ), pointer :: self

    real(rp) :: rrcut(1,1)
    character(:), allocatable :: s_by

    rrcut(:,:) = rcut
    s_by = "none"
    if( present(sort_by))s_by = sort_by
    self => compute_binned_pbc( nat, ityp, pos, lat, rrcut, s_by )
  end function compute_binned_pbc_one



  function compute_binned_pbc( nat, ityp, pos, lat, rcut, sort_by )result(self)
    !! version where `rcut` can have multiple values, depending on ityp value.
    implicit none
    !> number of atoms
    integer, intent(in) :: nat

    !> integer atomic types
    integer, intent(in) :: ityp(nat)

    !> atomic positions, shape(ndim, nat)
    real(rp), intent(in) :: pos(:,:)

    !> lattice vectors in rows, in units of `pos`
    real(rp), intent(in) :: lat(size(pos,1), size(pos,1))

    !> distance cutoff, in units of `pos`. Is a symmetric square matrix with values of cutoff for the
    !! corresponding `ityp` - `ityp`. For example `rcut(1,3) = 1.2` means the cutoff distance between
    !! atoms of type `ityp=1` and `ityp=3` is 1.2. The size of matrix `rcut` is thus related to how many
    !! different atomic types there are in the `ityp` list.
    real(rp), intent(in) :: rcut(:,:)

    !> optional string to sort the neighbour list; possible values `"distance"`, `"index"`
    character(*), intent(in), optional :: sort_by

    !> t_neighbour class
    type( t_neighbour ), pointer :: self

    integer :: i, j, ii, jj, kk, ll, mm
    integer :: typ_i, typ_j
    real(rp) :: dij, max_rcut
    real(rp), dimension(size(pos,1), size(pos,1)) :: invlat
    real(rp), dimension(size(pos,1)) :: r, rf, rij, kvec
    integer, dimension(3) :: ibin, addbin ! keep dimension 3, for lower dim they just contain no data
    integer :: xbin, ybin, zbin
    real(rp) :: rcut2

    integer, dimension(3) :: nbins, nbins_check
    integer :: ndim, idim
    integer :: ntyp
    real(rp), dimension(3) :: kmax, kadd
    real(rp), dimension(3) :: lo_buffer, hi_buffer
    integer :: jj_lb, kk_lb, ll_lb
    integer :: jj_ub, kk_ub, ll_ub
    integer :: i_start, i_end, n
    real(rp), allocatable :: d_o(:,:)

    ndim = size( pos,1 )
    ! self% ndim = ndim

    ! number of expected different values in ityp
    ntyp = size( rcut, 1 )
#ifdef DEBUG
    write(*,*) "got ntyp:",ntyp
#endif

    ! cutoff square
    max_rcut = maxval(rcut)
    rcut2 = max_rcut*max_rcut
#ifdef DEBUG
    write(*,*) "max_rcut=",max_rcut
#endif

    allocate( t_neighbour::self )

    ! inverse lattice
    call inverse( lat, invlat)

    ! maximal "displacement vector" in cartesian coords ...
    ! this is an overshoot proportional to sqrt(3.0)*rcut in a cubic box,
    ! but no idea in general.
    ! Is probably not optimal, but good-enough for now. Boxes with weird shapes
    ! could be problematic from this point of view (e.g. some ax similar norm
    ! to rcut, while others much larger).
    kvec(:) = max_rcut*1.0_rp
    ! see how this vector transforms to crist coords
    call cart_to_crist(kvec, lat, invlat)
    ! take the maxval as unified cutoff distance in crist coords over all axes
    ! kmax = maxval(abs(kvec))
    kmax(:) = maxval(abs(kvec))
    do idim = 1, ndim
       kmax(idim) = abs(kvec(idim))
    end do

    ! write(*,*) "kvec:", kvec
    ! write(*,*) 1.0/kvec(1)
    ! write(*,*) "kmax",kmax

    ! take the floor() for nbins, since fewer bins means each bin is larger in size
    nbins = floor(1.0_rp/kmax)
    ! protection against zero bins
    do idim = 1, ndim
       nbins(idim) = max(1,nbins(idim))
    end do

    ! nbins = nint(1.0_rp/kmax)
    ! nbins = 64
    ! nbins_check = ceiling(1.0_rp/kmax)

    allocate( self% nneig(1:nat), source=0)
    allocate( self% neiglist(1:first_alloc*nat) )
    allocate( self% veclist(1:ndim,1:first_alloc*nat) )
    allocate( self% partial_sumlist(1:nat) )
    self% partial_sumlist(1) = 0
    ! save first ntot (alloc size)
    self% ntot = first_alloc*nat


    ! buffer region near low-end
    ! lo_buffer = -0.5+kmax
    lo_buffer = kmax
    ! buffer region near high-end
    ! hi_buffer = 0.5-kmax
    hi_buffer(:) = 1.0_rp
    do idim = 1, ndim
       hi_buffer(idim) = 1.0_rp-kmax(idim)
    end do

#ifdef DEBUG
    write(*,*) "lo_buffer =",lo_buffer
    write(*,*) "hi_buffer =",hi_buffer
    write(*,*) "nbins=",nbins
    write(*,*) "rcut",rcut
#endif

    ! create bins
    allocate( self% bins(0:nbins(1)+1, 0:nbins(2)+1, 0:nbins(3)+1) )
    ! init bins to zero
    self% bins(:,:,:)% n_members = 0
    self% bins(:,:,:)% n_tot = 0
    self% bins(:,:,:)% tag = 0

    do i = 1, nat
       r = pos(:,i)

       ! atom needs to be inside box
       call cart_to_crist( r, lat, invlat )
       ! shift to -0.5:0.5
       call periodic(r)
       ! shift to 0:1
       r = r + [( 0.5_rp, idim=1,ndim )]

       ! write(*,*)
#ifdef DEBUG
       write(*,*) " ==> idx i:",i
       write(*,*) "pos",pos(:,i)
       write(*,*) "r",r
       ! write(*,*) r*4
#endif

       ! compute bin of this atom in reciprocal
       ! ibin goes from 0 to nbins.
       ! for ndim<3, the ibin keeps dimension 3, but has value 1 for extra dims
       ibin(:) = 1
       do idim = 1, ndim
          ! ibin(idim) = floor( r(idim)*nbins(idim) ) + 1
          ibin(idim) = nint( r(idim)*nbins(idim) )
       end do

#ifdef DEBUG
       write(*,"(a,1x,3i4, 3f9.4)") "ibin:",ibin, r
#endif

       ! get the cartesian coords of this atom non-shifted
       rf = r
       call crist_to_cart( rf, lat, invlat )
       ! add to bin
       call self% bins( ibin(1), ibin(2), ibin(3) )% bin_add( 1, i, rf )

       ! detect boundaries of the box:
       ! if atom is add the boundary, add a fake atom with identical index to the other side
       ! of the box. The fake atom has its coords shifted by appropriate vector.

       ! single components (boundary planes)
       do ii = 1, ndim
          addbin = ibin
          kadd = 0.0_rp
          if( check1( r(ii), kmax(ii), kadd(ii)) ) then
             ! need to add shifted atm
             rf = r + kadd(1:ndim)

             ! which bin to add? assume low-end, if kadd=1.0 then it's high-end
             addbin(ii) = 0
             if( kadd(ii) > 0.5_rp ) addbin(ii) = nbins(ii)

             call crist_to_cart( rf, lat, invlat )
             call self% bins( addbin(1), addbin(2), addbin(3))% bin_add( 11, i, rf )
          end if
       end do

       ! two components simultaneously (boundary edges)
       do ii = 1, ndim
          do jj = ii+1, ndim
             addbin = ibin
             kadd = 0.0_rp
             if( check2( r, ii, jj, kmax, kadd ) )then
                ! which bin: high-end where kadd=1.0; low-end where kadd=-1.0
                addbin = merge( nbins, addbin, kadd > 0.5_rp )
                addbin = merge( 0, addbin, kadd < -0.5_rp )
                rf = r + kadd(1:ndim)
                call crist_to_cart( rf, lat, invlat )
                call self% bins( addbin(1), addbin(2), addbin(3))% bin_add( 12, i, rf )
             end if
          end do
       end do

       ! three components simultaneously (boundary corners)
       do ii = 1, ndim
          do jj = ii+1, ndim
             do kk = jj+1, ndim
                addbin = ibin
                kadd = 0.0_rp
                if( check3( r, ii, jj, kk, kmax, kadd) ) then
                   ! which bin: high-end where kadd=1.0; low-end where kadd=-1.0
                   addbin = merge( nbins, addbin, kadd > 0.5_rp )
                   addbin = merge( 0, addbin, kadd < -0.5_rp )
                   rf = r + kadd(1:ndim)
                   call crist_to_cart( rf, lat, invlat )
                   call self% bins( addbin(1), addbin(2), addbin(3))% bin_add( 13, i, rf )
                end if
             end do
          end do
       end do

    end do

#ifdef DEBUG
    write(*,*) "start"
#endif


#ifdef WRBIN

    mm = 0
    do ii = 0, nbins(1)+1
       do jj = 0, nbins(2)+1
          do kk = 0, nbins(3)+1
             mm = mm + 1
             do i = 1, self% bins( ii,jj,kk)% n_members
                j = self%bins(ii,jj,kk)% global_idx_of_member(i)
                ! write(*,*) ityp(j), pos(:,j), j, mm
                write(*,*) ityp(j), self%bins(ii,jj,kk)% pos_member(:,i), j, ii,jj,kk
             end do
          end do
       end do
    end do
#endif



    ! head_idx is the last index in neiglist
    self% head_idx = 0
    do i = 1, nat

#ifdef DEBUG
       write(*,*) "&&& ", i
       write(*,*) "ityp", ityp(i)
       write(*,"(a,1x,3f9.4)") "pos",pos(:,i)
#endif

       r = pos(:,i)

       ! if rcut is passed as matrix, get the actual atom type
       if( ntyp > 1 ) then
          typ_i = ityp(i)
       else
          ! rcut has same value for all, ignore type
          typ_i = 1
       end if


       ! shift inside box; range 0:1 in crist
       call cart_to_crist( r, lat, invlat )
       call periodic( r )
       ! r = r + [0.5_rp, 0.5_rp, 0.5_rp ]
       r = r + [( 0.5_rp, i=1,ndim )]

       ! find my bin in reciprocal
       ibin(:) = 1
       do idim = 1, ndim
          ! ibin(idim) = floor( r(idim)*nbins(idim) ) + 1
          ibin(idim) = nint( r(idim)*nbins(idim) )
       end do

#ifdef DEBUG
       write(*,*) "r", r
       write(*,"(a,1x,3i4)")"ibin",ibin
#endif

       ! r into cartesian
       call crist_to_cart( r, lat, invlat )

       ! write(*,*) "r cart:",r
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

       ! if we are in bin=0 (low-edge), no need to check the -1 bin
       if( ibin(1) == 0 ) jj_lb = 0
       if( ibin(2) == 0 ) kk_lb = 0
       if( ibin(3) == 0 ) ll_lb = 0

       ! if we are in the bin=nbins+1 (high-edge), no need to check the +1 bin
       if( ibin(1) == nbins(1)+1 ) jj_ub = 0
       if( ibin(2) == nbins(2)+1 ) kk_ub = 0
       if( ibin(3) == nbins(3)+1 ) ll_ub = 0

       ! write(*,*) i, "ibin",ibin
#ifdef DEBUG
       write(*,*) "jj_lb", jj_lb, "jj_ub", jj_ub
       write(*,*) "kk_lb", kk_lb, "kk_ub", kk_ub
       write(*,*) "ll_lb", ll_lb, "ll_ub", ll_ub
#endif
       do jj = jj_lb, jj_ub
          do kk = kk_lb, kk_ub
             do ll = ll_lb, ll_ub

                ! set the bin where we take atoms
                xbin = ibin(1) + jj
                ybin = ibin(2) + kk
                zbin = ibin(3) + ll

                ! loop over the atoms of that bin
#ifdef DEBUG
                write(*,"(a,1x,3i4,3x,i0)") " >> ", xbin, ybin, zbin, self% bins(xbin,ybin,zbin)%n_members
#endif
                do ii = 1, self% bins( xbin, ybin, zbin )% n_members

                   ! vector between each member of the bin and r(i)
                   rij = self% bins( xbin, ybin, zbin )% pos_member(:,ii) - r

                   ! distance square from i to j
                   ! dij = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
                   dij = dot_product(rij,rij)

#ifdef DEBUG
                   j = self% bins( xbin, ybin, zbin )% global_idx_of_member(ii)
                   write(*,'(i16,1x,3f9.4,5x,f9.4,i12,i4)') &
                        ii, &
                        self%bins(xbin,ybin,zbin)%pos_member(:,ii), &
                        sqrt(dij), &
                        j, &
                        ityp(j)
#endif


                   if( dij .le. rcut2 ) then
                      ! get the true index of the atom j, and add it to neighbour list of i
                      j = self% bins( xbin, ybin, zbin )% global_idx_of_member(ii)
                      ! if rcut is passed as matrix, get the actual type of atom
                      if( ntyp > 1 ) then
                         typ_j = ityp(j)
                      else
                         ! rcut has same value for all
                         typ_j = 1
                      end if

                      ! add j based on actual value of rcut, depending on ityp
#ifdef DEBUG
                      write(*,*) "typ_i", typ_i, "typ_j", typ_j, "rcut(i,j)",rcut(typ_i, typ_j )
#endif
                      if( sqrt(dij) .le. rcut(typ_i, typ_j) ) then
                         call self% add(i, j, rij )
                      end if

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


       ! optional sorting
       if( present( sort_by ) ) then
          i_end = self% partial_sumlist(i)
          i_start = i_end - self% nneig(i)+1
          n = i_end - i_start + 1
          select case( sort_by )
          case( "index")
             call bubble_sort_i1d( self% neiglist(i_start:i_end))
          case( "distance" )
             allocate( d_o(1:2,1:n) )
             do ii = 1, n
                ! distance
                d_o(1,ii) = norm2( self% veclist(1:ndim, i_start+ii-1) )
                ! actual index
                d_o(2,ii) = real(i_start+ii-1,rp)
             end do
             call bubble_sort_r2d( n, d_o )
             self% veclist( 1:3, i_start:i_end ) = self% veclist(1:ndim, nint(d_o(2,:)) )
             self% neiglist(i_start:i_end) = self% neiglist(nint(d_o(2,:)))
             deallocate( d_o )
          end select
       end if

    end do

  end function compute_binned_pbc




  subroutine bin_add( self, bin_tag, idx, vec )
    !! add atom `idx` and vector `vec` to this bin
    class( t_bin ), intent(inout) :: self
    integer, intent(in) :: bin_tag
    integer, intent(in) :: idx
    real(rp), intent(in) :: vec(:)
    integer, allocatable :: tmp(:)
    real(rp), allocatable :: tmp_v(:,:)
    integer :: oldsize, newsize
    integer :: ndim
    ndim = size(vec)

    ! write(*,*) "   enter bin_add. n_tot:",self% n_tot
    ! write(*,*) "     n_members:",self% n_members
    ! write(*,*) "     tag:",self% tag
    ! write(*,*) bin_tag, vec, idx
    self% tag = bin_tag
    ! if( .not.allocated( self% global_idx_of_member) ) then
    if( self% n_tot == 0 ) then
       allocate( self% pos_member(1:ndim,1:bin_batchsize))
       allocate( self% global_idx_of_member(1:bin_batchsize))
       self% n_tot = bin_batchsize
       self% n_members = 0
       ! bin_tag = bin_tag + 1
       ! write(*,*) "   allocated new bin. n_tot:", self% n_tot
    end if

    ! realloc
    if( self% n_members + 1 .gt. self% n_tot ) then
       oldsize = self% n_tot
       newsize = oldsize + bin_batchsize
       ! write(*,*) " $$$ bin realloc", oldsize, newsize
       ! write(*,*) "  batchsize:",bin_batchsize
       allocate( tmp, source=self%global_idx_of_member )
       deallocate( self% global_idx_of_member )
       allocate( self%global_idx_of_member(1:newsize) )
       self%global_idx_of_member(1:oldsize) = tmp(:)
       deallocate(tmp)

       allocate( tmp_v, source=self%pos_member )
       deallocate( self% pos_member )
       allocate( self% pos_member(1:ndim,1:newsize) )
       self% pos_member(1:ndim,1:oldsize) = tmp_v(:,:)
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
    real(rp), intent(in) :: vec(:)       !! vec to add to veclist

    integer :: newsize, oldsize
    integer, allocatable :: tmp(:)
    real(rp), allocatable :: tmp_v(:,:)

    integer :: i_s, ndim

    ndim = size(vec)

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
       allocate( self% veclist(1:ndim, 1:newsize) )
       self% veclist(1:ndim,1:oldsize) = tmp_v(:,:)
       deallocate( tmp_v )
       ! increment ntot
       self% ntot = newsize
    end if

    ! add idx to head of neiglist
    self% neiglist( self% head_idx ) = idx

    ! add vec to head of veclist
    self% veclist( 1:ndim, self% head_idx ) = vec(1:ndim)

  end subroutine add



  pure subroutine periodic(c)
    !--------------------------------
    ! periodic boundary condition, for n dimensional vector input in crist coords.
    !--------------------------------
    implicit none
    real(RP), dimension(:),intent(inout) :: c
    integer :: i, ndim

    ndim = size(c)
    do i = 1, ndim
       if( c(i) .lt. -0.5_rp ) c(i) = c(i) + 1.0_rp
       if( c(i) .ge. 0.5_rp ) c(i) = c(i) - 1.0_rp
    end do

    ! c(1) = c(1) - int( (c(1)/ 0.5_rp) )
    ! c(2) = c(2) - int( (c(2)/ 0.5_rp) )
    ! c(3) = c(3) - int( (c(3)/ 0.5_rp) )
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
    real(rp), intent(inout) :: rij(:)
    real(rp), intent(in) :: lat(size(rij), size(rij))
    real(rp), intent(in) :: invlat(size(rij), size(rij))

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
    real(rp), intent(inout) :: rij(:)
    real(rp), intent(in) :: lat(size(rij), size(rij))
    real(rp), intent(in) :: invlat(size(rij), size(rij))

    rij = matmul( invlat, rij )
  end subroutine cart_to_crist

  pure subroutine inverse2x2( mat, inv )
    implicit none
    real(rp), intent(in) :: mat(2,2)
    real(rp), intent(out) :: inv(2,2)

    real(rp) :: det, invdet

    det = 0.0_rp

    ! calculate the determinant
    det = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
    ! inverse determinant
    invdet = 1.0_rp/det
    ! inverse matrix elements
    inv(1,1) = invdet  * mat(2,2)
    inv(1,2) = -invdet * mat(1,2)
    inv(2,1) = -invdet * mat(2,1)
    inv(2,2) = invdet  * mat(1,1)
  end subroutine inverse2x2

  pure subroutine inverse3x3( mat, inv )
    implicit none
    real(rp), intent(in) :: mat(3,3)
    real(rp), intent(out) :: inv(3,3)

    real(rp) :: det, invdet

    ! calculate the determinant
    det =  mat(1,1)*mat(2,2)*mat(3,3) &
         + mat(1,2)*mat(2,3)*mat(3,1) &
         + mat(1,3)*mat(2,1)*mat(3,2) &
         - mat(1,3)*mat(2,2)*mat(3,1) &
         - mat(1,2)*mat(2,1)*mat(3,3) &
         - mat(1,1)*mat(2,3)*mat(3,2)
    ! invert the determinant
    invdet = 1.0_rp/det
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

  pure subroutine inverse( mat, inv )
    implicit none
    real(rp), intent(in) :: mat(:,:)
    real(rp), intent(out) :: inv(size(mat,1), size(mat,1))

    integer :: ndim
    ndim = size(mat,1)
    select case( ndim )
    case( 1 )
       inv(1,1) = 1.0_rp/mat(1,1)
    case( 2 )
       call inverse2x2( mat, inv )
    case( 3 )
       call inverse3x3( mat, inv )
    end select
  end subroutine inverse

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
    real(rp), intent(in) :: r(:)
    integer, intent(in) :: i, j
    real(rp), intent(in) :: kmax(size(r))
    real(rp), intent(out) :: kadd(size(r))
    logical :: add
    logical :: ai, aj
    real(rp) :: k_cp(size(r))
    add = .false.
    kadd = 0.0_rp
    k_cp = kadd
    ai = check1( r(i), kmax(i), k_cp(i) )
    aj = check1( r(j), kmax(j), k_cp(j) )
    if( ai .and. aj ) then
       add = .true.
       kadd = k_cp
    end if
  end function check2

  function check3( r, i, j, k, kmax, kadd )result(add)
    ! check components i, j, k of r simultaneously
    real(rp), intent(in) :: r(:)
    integer, intent(in) :: i, j, k
    real(rp), intent(in) :: kmax(size(r))
    real(rp), intent(inout) :: kadd(size(r))
    logical :: add
    logical :: ai, aj, ak
    real(rp) :: k_cp(size(r))
    add = .false.
    kadd = 0.0_rp
    k_cp = kadd
    ai = check1( r(i), kmax(i), k_cp(i) )
    aj = check1( r(j), kmax(j), k_cp(j) )
    ak = check1( r(k), kmax(k), k_cp(k) )
    if( ai .and. aj .and. ak )then
       add = .true.
       kadd = k_cp
    end if
  end function check3



  ! bubble sort, modified from rosetta code:
  ! https://rosettacode.org/wiki/Sorting_algorithms/Bubble_sort#Fortran
  subroutine bubble_sort_i1d(a)
    integer, intent(inout), dimension(:) :: a
    integer :: temp
    integer :: i, j
    logical :: swapped

    do j = size(a)-1, 1, -1
       swapped = .false.
       do i = 1, j
          if (a(i) > a(i+1)) then
             temp = a(i)
             a(i) = a(i+1)
             a(i+1) = temp
             swapped = .true.
          end if
       end do
       if (.not. swapped) exit
    end do
  end subroutine bubble_sort_i1d
  ! d_o(1, :) = d
  ! d_o(2, :) = o
  subroutine bubble_sort_r2d( dim, a )
    integer, intent(in) :: dim
    real(rp), intent(inout), dimension(2,dim) :: a
    real(rp) :: temp(2)
    integer :: i, j
    logical :: swapped

    do j = dim-1, 1, -1
       swapped = .false.
       do i = 1, j
          if (a(1, i) > a(1, i+1)) then
             temp(1) = a(1,i)
             temp(2) = a(2,i)

             a(1,i) = a(1,i+1)
             a(2,i) = a(2,i+1)

             a(1,i+1) = temp(1)
             a(2,i+1) = temp(2)
             swapped = .true.
          end if
       end do
       if (.not. swapped) exit
    end do
  end subroutine bubble_sort_r2d


end module m_neighbour

