PROGRAM hdio
  USE mpi
  USE hdf5
  USE h5ds
  USE h5lt


  INTEGER, PARAMETER :: nxp=3,nyp=3,nzp=3
  INTEGER, PARAMETER :: nxpg=6,nypg=6!,nzp=3
  REAL :: data1(nzp,nxp,nyp), data2(nzp,nxp,nyp)
  REAL :: xt(nxp), yt(nyp), zt(nzp)

  INTEGER(HID_T)  :: writeprpID   ! writeout property list handle
  INTEGER(HID_T)  :: h5ID_       ! file ID
  
  INTEGER(HID_T) :: ds3d       ! 3d dataspace 

  INTEGER(HID_T) :: dims_(3), dimsx(1), dimsy(1), dimsz(1)  !! dimension lengths
  INTEGER(HID_T) :: counts_(3),countsx_(1), countsy_(1),countsz_(1)
  INTEGER(HID_T) :: offsets_(3), offsetsx_(1), offsetsy_(1),offsetsz_(1)

  INTEGER(HID_T) :: dsid3d(3)  ! dimension scale ids for 3d variables

  INTEGER :: comm_,mpie
  INTEGER :: mpi_size,mpi_rank

  comm_ = MPI_COMM_WORLD
  info = MPI_INFO_NULL

  dims_ = [nzp,nxpg,nypg]
  dimsx = [nxpg]
  dimsy = [nypg]
  dimsz = [nzp]

  CALL MPI_INIT(mpie)
  CALL MPI_COMM_SIZE(comm_,mpi_size,mpie)
  CALL MPI_COMM_RANK(comm_,mpi_rank,mpie)
  
  data1 = 1.*mpi_rank
  data2 = 10.*mpi_rank

  counts_ = [nzp,nxp,nyp]
  countsx_ = [nxp]
  countsy_ = [nyp]
  countsz_ = [nzp]
  offsetsz_ = [0]
  zt = [20.,40.,60.]
  IF (mpi_rank == 0) THEN
    offsets_ = [0,0,0]
    offsetsx_ = [0]
    offsetsy_ = [0]
    yt = [100.,200.,300.]
    xt = [100.,200.,300.]
  ELSE IF (mpi_rank == 1) THEN
    offsets_ = [0,nxp,0]
    offsetsx_ = [nxp]
    offsetsy_ = [nyp]
    yt = [100.,200.,300.]
    xt = [400.,500.,600.]
  ELSE IF (mpi_rank == 2) THEN
    offsets_ = [0,0,nyp]
    offsetsx_ = [0]
    offsetsy_ = [nyp]
    yt = [400.,500.,600.]
    xt = [100.,200.,300.]
  ELSE IF (mpi_rank == 3) THEN
    offsets_ = [0,nxp,nyp]
    offsetsx_ = [nxp]
    offsetsy_ = [nyp]
    yt = [400.,500.,600.]
    xt = [400.,500.,600.]
  END IF



  CALL testkind([1,2,3,4])

  CALL initialize(comm_,h5ID_)


  ! Create the axis variables and set as dimension scales
  CALL add_dimensions(h5ID_,dimsx,"xt")
  CALL add_dimensions(h5ID_,dimsy,"yt")
  CALL add_dimensions(h5ID_,dimsz,"zt")


  ! Create the 3d variable dataspace and the variables
  CALL get_ds_ids(h5ID_,3,["yt","xt","zt"],dsid3d)
  CALL set_varspace(dims_,ds3d)
  CALL add_variable(h5ID_,"data_1",ds3d,dsid3d)
  CALL add_variable(h5ID_,"data_2",ds3d,dsid3d)
  CALL Sclose(ds3d)
  CALL Dclose(dsid3d)


    ! Create the mpio property list
  CALL set_mpio_writeproperty(writeprpID)
  ! Write dimension axes vectors 
  CALL write_axis(h5ID_,"xt",offsetsx_,countsx_,dimsx,writeprpID,xt)
  CALL write_axis(h5ID_,"yt",offsetsy_,countsy_,dimsy,writeprpID,yt)
  CALL write_axis(h5ID_,"zt",offsetsz_,countsz_,dimsz,writeprpID,zt)

  CALL write_var3d(h5ID_,"data_1",offsets_,counts_,dims_,writeprpID,data1)
  CALL write_var3d(h5ID_,"data_2",offsets_,counts_,dims_,writeprpID,data2)

  ! Close the write out property list
  CALL Pclose(writeprpID)

  ! Close file
  CALL Fclose(h5ID_)
  
  ! Close the instance
  CALL terminate()
  
  CALL MPI_FINALIZE(mpie)
  
CONTAINS

  SUBROUTINE initialize(comm,h5ID)

    INTEGER, INTENT(in) :: comm 
    INTEGER(HID_T), INTENT(out) :: h5ID

    INTEGER :: err
    INTEGER(HID_T) filespace

      ! Initialize HDF5 library
    CALL h5open_f(err)

    ! Create file access property list
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, filespace, err )
    ! Set parallel IO in the property list
    CALL h5pset_fapl_mpio_f(filespace, comm, MPI_INFO_NULL, err)
  
    ! Create a new output file (CHECK THE FILE OPEN METHOD)
    CALL h5fcreate_f("testi.h5", H5F_ACC_TRUNC_F, h5ID, err, access_prp=filespace )

    ! Close the file property list instance
    CALL h5pclose_f(filespace, err)

  END SUBROUTINE initialize

  ! --------------------------------

  SUBROUTINE set_varspace(dlen,dspace)
    INTEGER(HID_T), INTENT(in) :: dlen(:)
    INTEGER(HID_T), INTENT(out) :: dspace
    INTEGER :: n, err
    n = SIZE(dlen)
    CALL h5screate_simple_f(n, dlen, dspace, err)
  END SUBROUTINE set_varspace

  ! ---------------------------------

  SUBROUTINE set_mpio_writeproperty(plist)
    INTEGER(HID_T), INTENT(out) :: plist
    INTEGER :: err
    CALL h5Pcreate_f(H5P_DATASET_XFER_F, plist, err)
    CALL h5Pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, err)
  END SUBROUTINE set_mpio_writeproperty

  ! ----------------------------------

  SUBROUTINE add_dimensions(h5ID,dlen,dname)
    INTEGER(HID_T), INTENT(in) :: h5ID
    INTEGER(HID_T), INTENT(in) :: dlen(1)
    CHARACTER(len=*), INTENT(in) :: dname

    INTEGER(HID_T):: did
    INTEGER(HID_T) :: dataspace
    INTEGER :: err

    ! Create dimension
    CALL h5screate_simple_f(1, dlen, dataspace, err)
    CALL h5dcreate_f(h5ID, dname, H5T_NATIVE_REAL, &
                     dataspace, did, err)

    ! make into dimension scale
    CALL H5DSset_scale_f(did,err,dname)

    CALL h5Sclose_f(dataspace,err)
    CALL Dclose([did])

    ! add writing!




  END SUBROUTINE add_dimensions

  ! ----------------------------

  SUBROUTINE add_variable(h5ID,vname,dspace,dsid)
    INTEGER(HID_T), INTENT(in) :: h5ID
    CHARACTER(len=*), INTENT(in) :: vname
    INTEGER(HID_T), INTENT(in) :: dspace
    INTEGER(HID_T), INTENT(in), OPTIONAL :: dsid(:)   ! Dimension scale ids
    INTEGER(HID_T) :: did 
    INTEGER :: err, i
    ! create the dataset
    CALL h5Dcreate_f(h5ID, vname, H5T_NATIVE_REAL, &
                     dspace, did, err)
    IF ( PRESENT(dsid) ) THEN
      DO i = 1,SIZE(dsid)
        CALL h5dsattach_scale_f(did,dsid(i),i,err)
        !CALL H5DSset_label_f(did, i, "aa",err)
      END DO
    END IF
    CALL Dclose([did])
  END SUBROUTINE add_variable

  ! -----------------------------


  SUBROUTINE write_var3d(h5ID,vname,offsets,counts,dims,wplist,data)
    INTEGER(HID_T), INTENT(in) :: h5ID
    CHARACTER(len=*), INTENT(in) :: vname
    INTEGER(HID_T), INTENT(in) :: offsets(3), counts(3), dims(3)
    INTEGER(HID_T), INTENT(in) :: wplist
    REAL, INTENT(in) :: data(:,:,:)

    INTEGER(HID_T) :: memspace,dspace,varID
    INTEGER :: err

    CALL h5Dopen_f(h5ID,vname,varID,err)
    CALL h5Dget_space_f(varID,dspace,err)

    CALL h5Sselect_hyperslab_f(dspace,H5S_SELECT_SET_F,offsets,counts,err)
    CALL h5Screate_simple_f(3,counts,memspace,err)
    CALL h5Dwrite_f(varID,H5T_NATIVE_REAL,data,dims,err, &
                    file_space_id=dspace, mem_space_id=memspace, &
                    xfer_prp=wplist)

    CALL Sclose(memspace)
    CALL Sclose(dspace)
    CALL Dclose([varID])

  END SUBROUTINE write_var3d

  SUBROUTINE write_axis(h5ID,vname,offsets,counts,dims,wplist,data)
    INTEGER(HID_T), INTENT(in) :: h5ID
    CHARACTER(len=*), INTENT(in) :: vname
    INTEGER(HID_T), INTENT(in) :: offsets(1), counts(1), dims(1)
    INTEGER(HID_T), INTENT(in) :: wplist
    REAL, INTENT(in) :: data(:)

    INTEGER(HID_T) :: memspace,dspace,varID
    INTEGER :: err

    CALL h5Dopen_f(h5ID,vname,varID,err)
    CALL h5Dget_space_f(varID,dspace,err)

    CALL h5Sselect_hyperslab_f(dspace,H5S_SELECT_SET_F,offsets,counts,err)
    CALL h5Screate_simple_f(1,counts,memspace,err)
    CALL h5Dwrite_f(varID,H5T_NATIVE_REAL,data,dims,err, &
                    file_space_id=dspace, mem_space_id=memspace, &
                    xfer_prp=wplist)



    CALL Sclose(memspace)
    CALL Sclose(dspace)
    CALL Dclose([varID])

  END SUBROUTINE write_axis

  SUBROUTINE get_ds_ids(h5ID,n,names,ids)
    INTEGER(HID_T), INTENT(in) :: h5ID
    INTEGER, INTENT(in) :: n
    CHARACTER(len=*), INTENT(in) :: names(n)
    INTEGER(HID_T), INTENT(out) :: ids(n)
    INTEGER :: i,err
    ids = 0
    DO i = 1,n
      CALL h5Dopen_f(h5ID,names(i),ids(i),err)
    END DO
  END SUBROUTINE get_ds_ids
  SUBROUTINE Dclose(vid)
    INTEGER(HID_T), INTENT(in) :: vid(:)
    INTEGER :: err,i
    DO i = 1,SIZE(vid)
      CALL h5dclose_f(vid(i),err)
    END DO
  END SUBROUTINE Dclose

  SUBROUTINE Sclose(dspace)
    INTEGER(HID_T), INTENT(in) :: dspace
    INTEGER :: err
    CALL h5sclose_f(dspace,err)
  END SUBROUTINE Sclose

  SUBROUTINE Pclose(prp)
    INTEGER(HID_T), INTENT(in) :: prp
    INTEGER :: err
    CALL h5pclose_f(prp,err)
  END SUBROUTINE Pclose

  SUBROUTINE terminate()
    INTEGER :: err
    CALL h5close_f(err)
  END SUBROUTINE terminate

  SUBROUTINE Fclose(h5ID)
    INTEGER(HID_T), INTENT(in) :: h5ID
    INTEGER :: err
    CALL h5fclose_f(h5ID,err)
  END SUBROUTINE Fclose


  SUBROUTINE testkind(offsets)
    INTEGER, INTENT(in) :: offsets(:)

    INTEGER(HID_T), ALLOCATABLE :: local(:)

    ALLOCATE(local,SOURCE=INT(offsets,KIND=HID_T))

    WRITE(*,*) offsets

  END SUBROUTINE testkind

END PROGRAM hdio
