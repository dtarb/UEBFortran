        !COMPILER-GENERATED INTERFACE MODULE: Mon Oct 08 18:57:05 2012
        MODULE DIRECTORYCREATE__genmod
          INTERFACE 
            SUBROUTINE DIRECTORYCREATE(NREFYR,NREFMO,NREFDAY,DIMLEN1,   &
     &DIMLEN2,DIMNAME1,DIMNAME2,DIMVALUE1,DIMVALUE2,DIMUNIT1,DIMUNIT2,  &
     &NUMOFFILE,OUTCOUNT,OUTVAR,NUMTIMESTEPPERFILE,OUTSAMPLEFILE,       &
     &OUTPUTNCCONTAINER,NCIDARRAY)
              INTEGER(KIND=4) :: OUTCOUNT
              INTEGER(KIND=4) :: NUMOFFILE
              INTEGER(KIND=4) :: DIMLEN1
              INTEGER(KIND=4) :: NREFYR
              INTEGER(KIND=4) :: NREFMO
              INTEGER(KIND=4) :: NREFDAY
              INTEGER(KIND=4) :: DIMLEN2
              CHARACTER(LEN=20) :: DIMNAME1
              CHARACTER(LEN=20) :: DIMNAME2
              REAL(KIND=4) :: DIMVALUE1(DIMLEN1)
              REAL(KIND=4) :: DIMVALUE2(DIMLEN1)
              CHARACTER(LEN=100) :: DIMUNIT1
              CHARACTER(LEN=100) :: DIMUNIT2
              INTEGER(KIND=4) :: OUTVAR(OUTCOUNT)
              INTEGER(KIND=4) :: NUMTIMESTEPPERFILE(NUMOFFILE)
              CHARACTER(LEN=200) :: OUTSAMPLEFILE(OUTCOUNT)
              CHARACTER(LEN=200) :: OUTPUTNCCONTAINER(NUMOFFILE,OUTCOUNT&
     &)
              INTEGER(KIND=4) :: NCIDARRAY(NUMOFFILE,OUTCOUNT)
            END SUBROUTINE DIRECTORYCREATE
          END INTERFACE 
        END MODULE DIRECTORYCREATE__genmod
