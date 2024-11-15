program main

  use m_neighbour

  implicit none
  integer :: nat
  integer, allocatable :: typ(:)
  real(rp), allocatable, target :: pos(:,:)
  real(rp) :: rcut
  real(rp) :: lat(3,3)
  type( t_neighbour ), pointer :: neigh
  integer, allocatable :: list(:), ityplist(:)
  integer :: n, idx
  integer :: i, j
  real(rp), allocatable :: veclist(:,:)
  character(256) :: line
  integer :: n_begin, n_end, nshells


  ! read xyz
  read(*,*) nat
  read(*,'(a256)') line
  n_begin = index(line, "Lattice=") + 9
  line = line(n_begin:)
  n_end = index(line,'"')-1
  line = line(:n_end)
  read(line, *) lat
  allocate( typ(1:nat) )
  allocate( pos(1:3,1:nat))
  do i = 1, nat
     read(*,*) typ(i), pos(:,i)
  end do

  ! cutoff radius
  rcut = 3.0_rp

  ! launch neighbour list
  ! neigh => t_neighbour( nat, typ, pos, lat, rcut )
  ! neigh => t_neighbour( nat, typ, pos, lat, rcut, sort_by="index" )
  neigh => t_neighbour( nat, typ, pos, lat, rcut, sort_by="distance" )


  ! arbitrary atom index
  idx = 3

  ! get number of first neighbors
  n = neigh% get_nn(idx)
  write(*,*) "number of first neighbors:", n

  ! get list of nearest neighbors
  write(*,*) "get list"
  n = neigh% get_nn( idx, list=list )
  write(*,*) n
  write(*,"(10i8)") list

  deallocate( list )


  ! get list of vectors expanded by `nbond=2`
  ! The list does not include the original idx itself, include it with `include_idx` flag
  write(*,*) "get list and veclist with nbond=2"
  n = neigh% get( idx, list=list, ityplist=ityplist, veclist=veclist, nbond=2, include_idx=.true. )

  ! Note that with `nbond>1` we lose the sorting order of neiglist...
  write(*,*) n
  write(*,*)
  do i = 1, n
     ! write(*,*) typ(list(i)),veclist(:,i), list(i)
     write(*,*) ityplist(i), veclist(:,i), list(i)
  end do
  deallocate( list, ityplist, veclist )


  write(*,*) repeat("-",20)

  ! create an arbitrary list with 4 indices
  allocate(list(1:4))
  list = [1, 34, 897, 2 ]
  ! and their vectors from `pos`
  allocate(veclist(1:3,1:4))
  veclist = pos(:,list)

  ! expand this list by 1 bond shell
  ! The result contains also the original list.
  nshells = 1
  n = neigh% expand( nshells, list, ityplist=ityplist, veclist=veclist )
  write(*,*) "expanded list:"
  write(*,"(10i8)") list
  write(*,*) n
  write(*,*)
  do i = 1, n
     ! write(*,*) typ(list(i)), veclist(:,i), list(i)
     write(*,*) ityplist(i), veclist(:,i), list(i)
  end do

  deallocate(list, veclist, ityplist )

  deallocate( typ, pos )
  deallocate( neigh )

end program
