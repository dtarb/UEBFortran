        !COMPILER-GENERATED INTERFACE MODULE: Mon Oct 08 18:57:04 2012
        MODULE READSV__genmod
          INTERFACE 
            SUBROUTINE READSV(PARAM,STATEV,SITEV,SVFILE,SLOPE,AZI,LAT,  &
     &SUBTYPE,DIMLEN2,DIMLEN1,ILAT,JLON,DTBAR,TS_LAST,LONGITUDE)
              REAL(KIND=4) :: PARAM(32)
              REAL(KIND=4) :: STATEV(6)
              REAL(KIND=4) :: SITEV(10)
              CHARACTER(LEN=50) :: SVFILE
              REAL(KIND=4) :: SLOPE
              REAL(KIND=4) :: AZI
              REAL(KIND=4) :: LAT
              INTEGER(KIND=4) :: SUBTYPE
              INTEGER(KIND=4) :: DIMLEN2
              INTEGER(KIND=4) :: DIMLEN1
              INTEGER(KIND=4) :: ILAT
              INTEGER(KIND=4) :: JLON
              REAL(KIND=4) :: DTBAR(12)
              REAL(KIND=4) :: TS_LAST
              REAL(KIND=4) :: LONGITUDE
            END SUBROUTINE READSV
          END INTERFACE 
        END MODULE READSV__genmod
