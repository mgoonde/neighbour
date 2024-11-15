# Description:

Periodic Boundary Conditions (pbc) neighbour list in fortran, using bins.
Encapsulated in a derived type [t_neighbour](https://mgoonde.github.io/neighbour/type/t_neighbour.html).



# Compile:

Compile the object with some useful flags:

```bash
gfortran -g -O2 -march=native -ffast-math -funroll-loops -c m_neighbour.f90
```

# Caller program:

To compile: link your program with `m_neighbour.o` compiled as above, and set the `-I` to this directory.
See also the `example` directory.

Pseudocode example program:

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
  ! get list of nearest neighbors of atom `my_idx`:
  block
    integer :: my_idx, n_list
    integer, allocatable :: list(:)
    n_list = neigh% get_nn( my_idx, list=list )
  end block

  !
  ! get the vectors and atomic types of 2nd shell neighbours of atom `my_idx`,
  ! and include the origin `my_idx` in the output:
  block
    integer :: my_idx, n_list
    real(wp), allocatable :: veclist(:,:)
    integer, allocatable :: ityplist(:)
    n_list = neigh% get( my_idx, veclist=veclist, ityplist=ityplist, nbond=2, include_idx=.true. )
  end block

  !
  ! set an arbitrary list of atomic indices, and expand it by 2 shells
  block
    integer :: nshells, n_list
    integer, allocatable :: list(:)
    list=[ 3, 15, 87, 32 ]
    nshells = 2
    n_list = neigh% expand( nshells, list )
  end block

end program
```

The cutoff value `rcut` can be an array where each element is a cutoff value between according pairs of atomic types. For example with 4 atomic types:

```f90
real(wp) :: rcut(4,4)  ! rcut(1,1) is the cutoff between atoms of type 1;
                       ! rcut(3,2) is the cutoff between atoms of type 3 and 2, etc.
                       ! the rcut matrix has to be symmetric.
```
