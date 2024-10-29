#gfortran -g -c m_neighbour.f90
#gfortran -g -o main.x main.f90 m_neighbour.o
gfortran -g -O2 -march=native -ffast-math -funroll-loops -c m_neighbour.f90
gfortran -g -O2 -march=native -ffast-math -funroll-loops -o main.x main.f90 m_neighbour.o
