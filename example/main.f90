program main

  use m_neighbour

  implicit none
  integer :: nat
  integer, allocatable :: typ(:)
  real(rp), allocatable, target :: pos(:,:)
  real(rp) :: rcut
  real(rp) :: lat(3,3)
  type( t_neighbour ), pointer :: neigh
  integer, allocatable :: list(:)
  integer :: n, idx
  integer :: i, j
  real(rp), allocatable :: veclist(:,:)
  character(256) :: line
  integer :: n_begin, n_end


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
  rcut = 6.0_rp

  ! launch neighbour list
  neigh => t_neighbour( nat, typ, pos, lat, rcut )


  ! print list of atom index 3
  idx = 3
  n = neigh% get_list( idx, list )
  write(*,*) n
  write(*,*) list

  ! print vectors of idx=3, sorted by distance
  n = neigh% get_veclist( idx, veclist, sort_by="distance" )
  ! it does not include the idx itself, include it manually
  write(*,*) n+1
  write(*,*)
  ! vector of idx is assumed at zero
  write(*,*) typ(idx), [0.0_rp, 0.0_rp, 0.0_rp]
  do i = 1, n
     write(*,*) typ(list(i)),veclist(:,i)
  end do
  deallocate( list, veclist )


  ! output the neigbor list of all atoms
  ! do i = 1, nat
  !    idx = i
  !    n = neigh% get_list( idx, list, sort_by="index" )
  !    write(*,"(i8)") idx
  !    ! write(*,"(a,i0)") "idx=",idx
  !    write(*,"(10i8)") list
  !    deallocate(list)
  ! end do


  deallocate( typ, pos )
  ! deallocate( list )
  ! deallocate( veclist )
  deallocate( neigh )

end program
