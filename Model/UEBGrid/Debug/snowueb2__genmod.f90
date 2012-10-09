        !COMPILER-GENERATED INTERFACE MODULE: Mon Oct 08 18:57:04 2012
        MODULE SNOWUEB2__genmod
          INTERFACE 
            SUBROUTINE SNOWUEB2(DT,NT,INPUT,SITEV,STATEV,TSPREVDAY,     &
     &TAVEPREVDAY,NSTEPDAY,PARAM,IFLAG,CUMP,CUMES,CUMEC,CUMMR,CUMGM,OUTV&
     &,MTIME,ATFF,CF,OUTARR)
              INTEGER(KIND=4) :: NT
              REAL(KIND=4) :: DT
              REAL(KIND=4) :: INPUT(8,NT)
              REAL(KIND=4) :: SITEV(*)
              REAL(KIND=4) :: STATEV(*)
              REAL(KIND=4) :: TSPREVDAY(*)
              REAL(KIND=4) :: TAVEPREVDAY(*)
              INTEGER(KIND=4) :: NSTEPDAY
              REAL(KIND=4) :: PARAM(*)
              INTEGER(KIND=4) :: IFLAG(*)
              REAL(KIND=4) :: CUMP
              REAL(KIND=4) :: CUMES
              REAL(KIND=4) :: CUMEC
              REAL(KIND=4) :: CUMMR
              REAL(KIND=4) :: CUMGM
              REAL(KIND=4) :: OUTV(*)
              REAL(KIND=4) :: MTIME(*)
              REAL(KIND=4) :: ATFF
              REAL(KIND=4) :: CF
              REAL(KIND=4) :: OUTARR(53)
            END SUBROUTINE SNOWUEB2
          END INTERFACE 
        END MODULE SNOWUEB2__genmod
