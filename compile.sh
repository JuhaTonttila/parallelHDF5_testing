mpif90 -c hdio.f90 -I /usr/lib/x86_64-linux-gnu/hdf5/mpich/include/ 
mpif90 *.o -L/usr/lib/x86_64-linux-gnu/hdf5/mpich/lib  -lhdf5_fortran  -lhdf5hl_fortran #-lhdf5_hl -lhdf5