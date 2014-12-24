        Subroutine StateVarriablenetCDFCreate(WATERSHEDfILE,wsycoordinate,wsxcoordinate,stateNCIDARRAY,nstepday,stateoutSymbol,stateOutUnits,stateOutfilelist)
        use netCDF
        implicit none
        integer n,i,j,LengthOutFile
        integer::nstepday
        parameter(n=4)
        character*200:: stateoutSymbol(n+nstepday), stateoutFile(n+nstepday), stateOutUnits(n+nstepday),stateOutfilelist(n+nstepday)
        integer:: outvar(n+nstepday), outvarpos(n+nstepday), stateNCIDARRAY(n+nstepday)
        character(len=10) :: IntToChar
        integer, parameter :: NDIMS=2 ! x-coord and y-coord
        character*50:: WATERSHEDFILE,NCOutfileName
        character*50:: wsycoordinate,wsxcoordinate
        integer:: ncidout, dimid1, dimid2, vartype, varid
        Integer:: dimlen1,dimlen2,dimlen1OUT,dimlen2OUT,NumOP
        character*200:: statencfolder,statencfoldertemp
        byte, allocatable:: inbyte1(:)
        character, allocatable:: inchar1(:)
        integer*2, allocatable:: inshort1(:)
        integer*4, allocatable:: ininteger1(:)
        real*4, allocatable:: inreal1(:)
        REAL*8, allocatable:: indouble1(:)
        Character*200, allocatable:: OutPointFiles (:)
        byte, allocatable:: inbyte2(:)
        character, allocatable:: inchar2(:)
        integer*2, allocatable:: inshort2(:)
        integer*4, allocatable:: ininteger2(:)
        real*4, allocatable:: inreal2(:)
        real*8, allocatable:: indouble2(:)
        integer::scanlength
        Real*8, allocatable:: dimvalue1(:),dimvalue2(:)
        Real*8, allocatable:: DimValue1Org(:),DimValue2Org(:)
        integer, parameter:: NF90_BYTEs1 = 1, NF90_CHARs1 = 2, NF90_SHORTs1 = 3, NF90_INTs1 = 4, NF90_FLOATs1 = 5, NF90_DOUBLEs1 = 6
        integer, parameter:: NF90_BYTEs2 = 1, NF90_CHARs2 = 2, NF90_SHORTs2 = 3, NF90_INTs2 = 4, NF90_FLOATs2 = 5, NF90_DOUBLEs2 = 6
        character (20) :: DimName1,DimName2
        character (100) :: DimUnit1,DimUnit2

        Integer:: ncid
        integer :: dimids(NDIMS)
        integer::x_dimid,y_dimid,time_dimid,x_varid,y_varid,time_varid
        integer:: Varid4
        character (len = *), parameter :: UNITS = "units"
        character (len = *), parameter :: missing_value = "missing_value"
        Real:: MissingValues,Fillvalues
        double precision:: DimValueX

        call check(nf90_open(WATERSHEDfILE,nf90_nowrite, ncidout))           ! open the netcdf file
        call check(nf90_inq_dimid(ncidout,wsycoordinate,DimID1))
        call check(nf90_inq_dimid(ncidout,wsxcoordinate,DimID2))
        call check(nf90_inquire_dimension(ncidout,DimID1,len=dimlen1))       ! Information about dimensionID 1
        call check(nf90_inquire_dimension(ncidout,DimID2,len=dimlen2))       ! information about dimensionID 2  
        
        allocate(dimvalue1(dimlen1))
        allocate(dimvalue2(dimlen2))
        
        Allocate(inbyte1(dimlen1))
        Allocate(inchar1(dimlen1))
        Allocate(inshort1(dimlen1))
        Allocate(ininteger1(dimlen1))
        Allocate(inreal1(dimlen1))
        Allocate(indouble1(dimlen1))
        Allocate(DimValue1Org(dimlen1))

        Allocate(inbyte2(dimlen2))
        Allocate(inchar2(dimlen2))
        Allocate(inshort2(dimlen2))
        Allocate(ininteger2(dimlen2))
        Allocate(inreal2(dimlen2))
        Allocate(indouble2(dimlen2))
        Allocate(DimValue2Org(dimlen2))
        MissingValues=-9999.0
        Fillvalues=-9999.0    
        call check(nf90_inq_varid(ncidout,wsycoordinate, VarId))        ! information about variableID for a given VariableName
        call check(NF90_inquire_variable(ncidout,VarID,xtype=vartype))
        call check(nf90_get_att(ncidout,Varid, UNITS,DimUnit1)) 

        if(vartype .eq. NF90_BYTEs1)then
            call check(nf90_get_var(ncidout,Varid,inbyte1)) 
            DimValue1Org=REAL(inbyte1,8)
        elseif(vartype .eq. NF90_CHARs1)then
            call check(nf90_get_var(ncidout,Varid,inchar1)) 
            Write (6, *) "Error: A site variable can't be a character in ", WATERSHEDfILE
        elseif(vartype .eq. NF90_SHORTs1)then
            call check(nf90_get_var(ncidout,Varid,inshort1)) 
            DimValue1Org=REAL(inshort1,8)
        elseif(vartype .eq. NF90_INTs1)then
            call check(nf90_get_var(ncidout,Varid,ininteger1)) 
            DimValue1Org=REAL(ininteger1,8)
        elseif(vartype .eq. NF90_FLOATs1)then 
            call check(nf90_get_var(ncidout,Varid,inshort1)) 
            DimValue1Org=REAL(inreal1,8) 
        elseif(vartype .eq. NF90_DOUBLEs1)then
            call check(nf90_get_var(ncidout,Varid,indouble1))  
            DimValue1Org=REAL(indouble1,8) 
        END IF 
        DimValue1=DimValue1Org

        call check(nf90_inq_varid(ncidout,wsxcoordinate, VarId))     ! information about variableID for a given VariableName 
        call check(NF90_inquire_variable(ncidout,VarID,xtype=vartype))
        call check(nf90_get_att(ncidout,Varid, UNITS,DimUnit2)) 

        if(vartype .eq. NF90_BYTEs2)then
            call check(nf90_get_var(ncidout,Varid,inbyte2))
            DimValue2Org=REAL(inbyte2,8)
        elseif(vartype .eq. NF90_CHARs2)then
            call check(nf90_get_var(ncidout,Varid,inchar2))
            Write (6, *) "Error: A site variable can't be a character in ", WATERSHEDfILE
        elseif(vartype .eq. NF90_SHORTs2)then
            call check(nf90_get_var(ncidout,Varid,inshort2))
            DimValue2Org=REAL(inshort2,8)
        elseif(vartype .eq. NF90_INTs2)then
            call check(nf90_get_var(ncidout,Varid,ininteger2))
            DimValue2Org=REAL(ininteger2,8)
        elseif(vartype .eq. NF90_FLOATs2)then
            call check(nf90_get_var(ncidout,Varid,inreal2)) 
            DimValue2Org=REAL(inreal2,8) 
        elseif(vartype .eq. NF90_DOUBLEs2)then
            call check(nf90_get_var(ncidout,Varid,indouble2)) 
            DimValue2Org=REAL(indouble2,8) 
        END IF 
        DimValue2=DimValue2Org
        
        Deallocate(inbyte1)
        Deallocate(inchar1)
        Deallocate(inshort1)
        Deallocate(ininteger1)
        Deallocate(inreal1)
        Deallocate(indouble1)
        Deallocate(DimValue1Org)
        
        Deallocate(inbyte2)
        Deallocate(inchar2)
        Deallocate(inshort2)
        Deallocate(ininteger2)
        Deallocate(inreal2)
        Deallocate(indouble2)
        Deallocate(DimValue2Org)

       ! the symbol table element lengths have been expanded to match
       ! fixed width lengths.  This accomidates cross-compiler
       ! differences and conforms to later Fortran-2003 standards.
       
 

       do i=1,(n+nstepday)
            CALL create_directory("state")
            NCOutfileName="state"//"\"//trim(stateOutfilelist(i))
            CALL lowercase(stateoutSymbol(i),stateoutSymbol(i))
            CALL lowercase(stateOutUnits(i),stateOutUnits(i))
            call check(nf90_create(NCOutfileName,nf90_clobber,ncid))
            call check(nf90_def_dim(ncid,wsycoordinate,dimlen1,y_dimid))
            call check(nf90_def_dim(ncid,wsxcoordinate,dimlen2,x_dimid))
            ! Define the coordinate variables. We will only define coordinate
            ! variables for lat and lon.  Ordinarily we would need to provide
            ! an array of dimension IDs for each variable's dimensions, but
            ! since coordinate variables only have one dimension, we can
            ! simply provide the address of that dimension ID (lat_dimid) and
            ! similarly for (lon_dimid).
            call check(nf90_def_var(ncid,wsycoordinate,NF90_DOUBLE,y_dimid, y_varid))
            call check(nf90_def_var(ncid,wsxcoordinate,NF90_DOUBLE,x_dimid, x_varid))
            ! Assign units attributes to coordinate variables.
            call check(nf90_put_att(ncid,y_varid,UNITS,DimUnit1))
            call check(nf90_put_att(ncid,x_varid,UNITS,DimUnit2))
            ! The dimids array is used to pass the dimids of the dimensions of
            ! the netCDF variables. Both of the netCDF variables we are creating
            ! share the same four dimensions. In Fortran, the unlimited
            ! dimension must come last on the list of dimids.
            dimids = (/y_dimid,x_dimid/)               
            ! Define the netCDF variables for the pressure and temperature data.
            call check( nf90_def_var(ncid,trim(stateOutSymbol(i)), NF90_REAL, dimids,Varid4))
            ! Assign units attributes to the netCDF variables.
            call check(nf90_put_att(ncid,Varid4,UNITS,trim(stateOutUnits(i))))
            call check(nf90_put_att(ncid,Varid4,missing_value,MissingValues)) 
            call check(nf90_put_att(ncid,Varid4,'_FillValue', FillValues))
            call check(nf90_enddef(ncid))
            CALL Check(NF90_PUT_VAR(ncid,y_varid,DimValue1)) ! latitude/y is dimension 1
            CALL Check(NF90_PUT_VAR(ncid,x_varid,DimValue2)) ! longtude/x is dimension 2
            CALL check(nf90_sync(ncid))
            stateNCIDARRAY(i)=ncid
        end do
        
        end Subroutine StateVarriablenetCDFCreate