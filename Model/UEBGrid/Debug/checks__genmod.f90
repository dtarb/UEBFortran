        !COMPILER-GENERATED INTERFACE MODULE: Mon Oct 08 18:57:07 2012
        MODULE CHECKS__genmod
          INTERFACE 
            SUBROUTINE CHECKS(SVFILE,MAXNUMOFFILE,ISINPUTFROMNC,        &
     &NUMNCFILES,TOTALNC,STATESITEVNAME)
              INTEGER(KIND=4) :: MAXNUMOFFILE
              CHARACTER(LEN=200) :: SVFILE
              INTEGER(KIND=4) :: ISINPUTFROMNC(11)
              INTEGER(KIND=4) :: NUMNCFILES(11)
              INTEGER(KIND=4) :: TOTALNC
              CHARACTER(LEN=200) :: STATESITEVNAME(32)
            END SUBROUTINE CHECKS
          END INTERFACE 
        END MODULE CHECKS__genmod
