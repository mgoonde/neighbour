pbc neighbour list in fortran, using bins.

Compile:

```bash
gfortran -g -O2 -march=native -ffast-math -funroll-loops -c m_neighbour.f90
```

Use in program:

```f90
program main
  use m_neighbour, only: t_neighbour

  type( t_neighbour ), pointer :: neigh

  integer                  :: nat       ! number of atoms
  integer, allocatable     :: typ(:)    ! allocate to size (1:nat)
  real(wp), allocatable    :: pos(:,:)  ! allocate to size (1:3, 1:nat)
  real(wp), dimension(3,3) :: box       ! lattice vectors in rows
  real(wp)                 :: rcut      ! distance cutoff

  ...

  !
  ! compute the neighbour list with distance cutoff `rcut`:
  neigh => t_neighbour( nat, typ, pos, box, rcut )

  !
  ! get list of neighbors of atom `my_idx`:
  block
    integer :: my_idx, n_list
    integer, allocatable :: list(:)
    n_list = neigh% get_list( my_idx, list )
  end block

  !
  ! get the vectors of neighbours of atom `my_idx`:
  block
    integer :: my_idx, n_list
    real(wp), allocatable :: veclist(:,:)
    n_list = neigh% get_veclist( my_idx, veclist )
  end block

end program
```
