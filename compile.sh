#gfortran -g -c m_neighbour.f90
gfortran -g -cpp -O2 -march=native -ffast-math -funroll-loops -c m_neighbour.f90
#gfortran -g -cpp -DDEBUG -O2 -march=native -ffast-math -funroll-loops -c m_neighbour.f90
