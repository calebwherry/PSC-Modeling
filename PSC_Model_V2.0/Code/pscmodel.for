*     PROGRAM PSC_MODEL_E
*
*     Main program used to run the PSCBOX model in the Eulerian-radius-space
*     version, i.e. using fixed size bins.
*     The program file PSCMODEL.FOR must be linked to PSCBOX.FOR holding the
*     microphysical, thermodynamical, and optical subroutines of the PSC box
*     model.
*
*     A few lines of code in this program must be changed, depending on
*     the computer (operative system) to be used and where input/output files
*     are stored. The statements of the program, which have to be changed, are
*     all preceeded by a comment line, starting with 'CMP'.  All instructions in
*     comment lines, starting with 'CMP', should be invoked, before compiling
*     and running the program. Comment lines starting with 'CMP' are only found
*     in this main program file.
*
*     Input/Output files:
*     -------------------
*
*     The program uses a simple self-explaining ASCII input file, stored in
*     the associated input file directory. The input file can have any name
*     with extension '.INP'. Output file names are generated from the applied
*     input file name, all with extension .DAT and stored in the associated
*     output file directory. The input/output directories must be specified
*     in the main program before compiling (cf. CMP comment lines).
*
*     Time units.
*     -----------
*     The model performs all calculations internally in SI units, i.e. time in
*     seconds. However, other units for input/output can be required.
*     The units (Days, Hours, Minutes, Seconds) of the time specifications in
*     the input file and the tempetature table file can be set in the input file.
*     These units are also used for output.
*
*     Computing control:
*     ------------------
*
*     Two numbers can be set in the input file to control the integration
*     of the PSC-model:
*         Maximum integration step size (s)
*
*         Integration start and stop time (Days, Hours, Minutes, Seconds);
*             when the integration has reached the stop time, the program
*             is terminated (in batch mode), or control is transferred to
*             the interactive routine (using the interactive mode on a PC).
*
*     Output control:
*     ---------------
*
*     A number can be set in the input file to control the frequency
*     of output from the PSC-model (Days, Hours, Minutes, Seconds).
*
*     Model structure:
*     ----------------
*
*     The number of layers in the vertical direction is specified in the
*     input file. A maximum of MAXLAY layers can be used. If a negative 
*     number is specified, the absolute value is used as the number of 
*     layers, but no sedimentation of particles between the layers is 
*     calculated. The layers are numbered from the top to the bottom, 
*     and the top layer is no. 1.
*
*     The potential temperature THETA (K) of the bottom layer and the layer
*     thickness dPOT (K) (delta-potential temperature) are specified in
*     the input file which, together with temperature/pressure specification
*     in the input fili, will have the following influence:
*
*     Temperature, pressure, and potential temperature:
*     -------------------------------------------------
*
*     The temperature can be calculated from a function of time as
*
*         TAIR(time) = RAMP-function(time) 
*
*     The RAMP-function is a continuous piece-wise linear function, specified
*     by 4 coordinate pairs: { (ramp time i, temperature i); i=1..4 };
*     i.e. the graph of this function is made up of 3 straight line segments
*     between the 4 coordinate points. For time values less than ramp time 1,
*     the function is constant equal to temperature 1, and for time values
*     greater than ramp time 4, the function is constant equal to
*     temperuture 4. The 4 coordinate pairs are specified in the input file.
*     The same temperatur applies for all layers when using the RAMP-function.
*     Isentropic conditions are assumed and pressures are calculated from
*     temperatures and the specified potenteal temperatures in each layer.
*     If the RAMP temperature function is to be used, the temperature calculation
*     method flag in the input file must be specified as 0 (zero).
*
*     Alternatively, the temperature history may be specified individually in
*     each layer from a time table of the temperature. The table must be stored
*     in the input directory as an ASCII file with 1+NLAYER columns
*     (time,temp1,temp2...temp-nlayer), (Kelvin), using the above speficied time
*     units. Cubic spline interpolation will be used between the table entries.
*     Time values in the table must be strictly increasing.
*     The name of the temperature table file is specified in
*     the input file, residing in the input directory.
*     Isentropic conditions are assumed and pressures are calculated from
*     temperatures and the specified potential temperatures in each layer.
*     If the temperature table is to be used, the temperature calculation method
*     flag in the input file must be specified as 1 (one).
*
*     Finally, both temperature and pressure can be specified for each layer in
*     an input table file, giving 1+(2*NLAYER) columns
*     (time, temp1, pres1, temp2, pres2,...,temp-nlayer,pres-nlayer)(Kelvin,Pa)
*     The temperature calculation method flag must set to 2 (two).
*     In this case the specified potential temperature is ignored, but dPOT is
*     used to specify the top layer thickness.
*
*     If a temperature (pressure) input file is used the ramp temperature
*     information in the input file is ignored. If a ramp temperature function
*     is used the specification of temperature table file name is ignored.
*
*     The temperature calculations as described above may be overlayed
*     by adding a sinosoidal temperature oscillation. Amplitude (K) and period
*     (user units) are specified in the input file. If no oscilations are to be
*     used, specify amplitude AND period equal to zero. To apply a constant
*     temperature correction in all layers, specify the amplitude with sign,
*     and set the period equal to zero.
*
*
*     Initial values:
*     ---------------
*
*     Initial values of water vapor (ppmv) is specified in the input file,
*     one value for each layer.
*
*     Initial values of nitric acid (ppbv) is specified in the input file,
*     one value for each layer. If a positive value is given, this value is used;
*     if zero is specified, an observed vertical profile of nitric acid from
*     LIMS-data is used, (cf. JGR 89,5125,1984), either from  a northern
*     (positive latitude) hemisphere or southern (negative latitude) hemisphere
*     data set.
*
*     Initial values of the total number density of sulfuric aerosols (1/ccm)
*     is specified in the input file (same for all layers). If a positive value
*     is given, this value is used for all levels; if zero is specified, number
*     density is derived in consistency with SAGE I and II aerosol surface
*     climatology (Hitchman et al., JGR 99,20689, 1994).
*
*     Median radius and geometric standard deviation of the initial 
*     lognormal distribution is given in the input file (same for all 
*     altitudes).
*
*
*     Optical calculations:
*     ---------------------
*     Calculation of particle backscatter ratios, extinction coefficients and
*     depolarisation (Mie scattering or T-matrix) can be specified by a logical
*     flag in the input file. Directory of T-matric coeffifients must be
*     specified in the main program before compiling (cf. CMP comment lines).
*     Wavelengths are speficied as parameters in the main programme and
*     must be set before compilation. The number of size bins must be large
*     to obtain reliable backscatter values. IF(OPTICS)-ENDIF blocks in the main
*     program and subroutine DUMP should also be inspected and changed according
*     to the actual optical parameters required.
*
C----------------------------------------------------------------------------
      PROGRAM PSC_MODEL_E
      IMPLICIT NONE
C----------------------------------------------------------------------------
C     Model size:
      INTEGER NBINS,TYPES,NCLASS,MAXLAY,NLAYER,NLAY,NWORK
C
      PARAMETER (
C         Number of particle size bins:
     +        NBINS=60,
C     +        NBINS=600,
C         Number of types of particles comprehended by the model:
     +        TYPES=4,
C         Dimension of work arry:
     +        NWORK=13,
c         Number of size classes in culumated size distributions:
     +        NCLASS=13,
C         Maximum number of layers:
     +        MAXLAY=25)
C----------------------------------------------------------------------------
C     Particle types:
      INTEGER STS,SAT,NAT,ICE
      PARAMETER(
C     Liquid supercooled ternary solution particle (SA or Type 1b PSC)
     +  STS=1,
C     Solid sulfuric acid tetrahydrate (SAT) particle; i.e. no HNO3 or excess ice
     +  SAT=2,
C     Solid nitric acid trihydrate (NAT) particle; i.e. NAT and SAT but no excess ice (Type 1a PSC)
     +  NAT=3,
C     Solid ice particle; i.e. holding excess ice, NAT, and SAT (Type 2 PSC)
     +  ICE=4)
C----------------------------------------------------------------------------
CMP   Directory names:
CMP   Give full directory names. Check slash conventions (\Windows /Unix).
CMP   Include ending slash in directory names.
CMP   Input directory name:
      CHARACTER*(*), PARAMETER :: INPDIR = '.\Simulations\'
CMP   Output directory name:
      CHARACTER*(*), PARAMETER :: OUTDIR = '.\Simulations\'
	!!
	!!Caleb Wherry
	!!Below (old) notation intializes input and output folders to PWD folders.
CMP   Input directory name:
      !!CHARACTER*10 INPDIR
	!!DATA INPDIR/'.\inpdata\'/
CMP   Output directory name:
      !!CHARACTER*10 OUTDIR
      !!DATA OUTDIR/'.\outdata\'/
	!!
C----------------------------------------------------------------------------
CMP   Optical calculations. Inspect PARAMETERS and
CMP   IF(OPTICS)-ENDIF blocks to check calculations.
CMP   Directory (PWD) where T-Matrix expansion coefficients are stored:
      CHARACTER*(*), PARAMETER :: T_MATRIX_DIR = '.\T-matrix\'
	!!
	!!Caleb Wherry
	!!Below (old) notation only reads file from the C: base directory.
	!!New notation (above) reads file from Present Working Directory (PWD).
	!!
C	CHARACTER *14, T_MATRIX_DIR 
C	DATA T_MATRIX_DIR/'\T_matrix\'/
	!!

C     Optical variables:
      LOGICAL OPTICS
C     Wavelengths (m) for backscatter ratio and extinction calculations.
      INTEGER N_WAVE
      PARAMETER(N_WAVE=7)
      REAL WAVELENGTH(N_WAVE)
C     Wavelengths in meters:
C     DATA WAVELENGTH/0.449E-6,0.480E-6,0.532E-6,
C    +                0.685E-6,0.940E-6,1.022E-6,1.064E-6/
      DATA WAVELENGTH/0.354E-6,0.449E-6,0.532E-6,
     +                0.779E-6,1.545E-6,1.022E-6,1.064E-6/
C     DATA WAVELENGTH/0.354E-6,0.440E-6,0.603E-6,
C    +                0.779E-6,0.922E-6,1.020E-6,1.064E-6/
      REAL ANGLE(N_WAVE)
C     Backscatter angles corresponding to each wavelength:
C      DATA ANGLE/180.0,173.0,180.0,
C     +           180.0,173.0,180.0,180.0/
      DATA ANGLE/180.0,180.0,180.0,
     +           180.0,180.0,180.0,180.0/
C     Refractive indices for each particle type and wavelength:
      REAL REFRACTIVE_INDEX(TYPES,N_WAVE)
      DATA REFRACTIVE_INDEX/
c     at wavelength 1: 354 nm
     +1.52,1.42,1.51,1.31,
c     at wavelength 2: 449 nm
     +1.52,1.42,1.51,1.31,
c     at wavelength 3: 532 nm
     +1.55,1.42,1.55,1.31,
c     at wavelength 4: 779 nm
     +1.57,1.42,1.57,1.31,
c     at wavelength 5: 1545 nm (I changed the ice refractive index to match the Warren data set at this wavelength)
     +1.50,1.42,1.46,1.29,
c     at wavelength 6: 1022 nm
     +1.50,1.42,1.60,1.31,
c     at wavelength 7: 1064 nm
     +1.50,1.42,1.51,1.31/
C     Aspect ratios for each particle type:
      REAL ASPECT(TYPES)
      DATA ASPECT/1.0,1.05,1.05,1.55/
C      DATA ASPECT/1.0,1.05,1.05,1.55/
C     Molecular depolarization ratio:
      REAL DEPOL_M
      DATA DEPOL_M/0.014/
      REAL EXTINCTION_CROSS(NBINS,TYPES,N_WAVE),
     +     BACKSCAT_CROSS(NBINS,TYPES,N_WAVE),
     +     DEPOLARIZATION(NBINS,TYPES,N_WAVE),
     +     AEROSOL_BACKSCATTER_RATIO(N_WAVE,MAXLAY),
     +     EXTINCTION(N_WAVE,MAXLAY),
     +     DEPOL(N_WAVE,MAXLAY)
      REAL RAYLGH
      REAL RAY,OPT1(MAXLAY),OPT2(MAXLAY),OPT_INDEX(MAXLAY)
C----------------------------------------------------------------------------
C     Aerosol and PSC variables:
      REAL
     +        ND(NBINS,TYPES,MAXLAY),
     +        MCS(NBINS,TYPES,MAXLAY),
     +        MCN(NBINS,TYPES,MAXLAY),
     +        MCW(NBINS,TYPES,MAXLAY),
     +        SDND(NBINS,TYPES),
     +        SDMCS(NBINS,TYPES),SDMCN(NBINS,TYPES),SDMCW(NBINS,TYPES),
     +        PTSIZE(NBINS,8),WORK(NBINS,NWORK),
     +        WSA(MAXLAY),WNA(MAXLAY)
      REAL    RMIN,RMAX,
     +        SIZEDIST(NCLASS,MAXLAY),RADIUSCLASS(NCLASS)
      LOGICAL PRESENCE(TYPES,MAXLAY),SEDMNT(TYPES)
C----------------------------------------------------------------------------
C     Ambient air state variables:
      REAL
     +        TAIR(MAXLAY),
     +        PAIR(MAXLAY),PRESS(MAXLAY),
     +        MRWV(MAXLAY),MRNA(MAXLAY),
     +        MRWVI(MAXLAY),MRNAI(MAXLAY),
     +        PPWV,PPNA
C----------------------------------------------------------------------------
C     Time variables:
      REAL TIME,
     +        DT,DTMAX,
     +        PTIME,DTOUT,
     +        VTIME,DTIME,DTVTME,DTDTME,
     +        TSTART,TSTOP,
     +        DTEMPMAX,DTEMP,TEMPE,DTMIN
      PARAMETER(
C         Mimimum time step size (sec):
     +        DTMIN=0.1D0,
C         Maximum temperature change (K) per time step:
     +        DTEMPMAX=0.1)
C----------------------------------------------------------------------------
C     Altitude and latitude variables:
      REAL HEIGHT(MAXLAY),POT(MAXLAY),THICKNESS(MAXLAY),
     +     THETA,DPOT,LATITU
C----------------------------------------------------------------------------
C     Temperature perturbation variables:
      REAL TEMP1,TEMP2,TEMP3,TEMP4,RAMPT1,RAMPT2,RAMPT3,RAMPT4,
     +     OMEGA,TAMPLI,DELT,POTTMP
      INTEGER MAXTAB,NPOINT,IOS,TFUNC
C----------------------------------------------------------------------------
C     Maximum of records in the 1-D temperature/pressure table file:
      PARAMETER(MAXTAB=31000)
      REAL TIMEA(MAXTAB),
     +     TEMPA(MAXTAB,MAXLAY),PRESA(MAXTAB,MAXLAY),
     +     Y2TEMP(MAXTAB,MAXLAY),Y2PRES(MAXTAB,MAXLAY),
     +     YWORK(MAXTAB)
      CHARACTER*1 TUNIT
C----------------------------------------------------------------------------
C     Surface pressure (Pa) and scale height (m):
      REAL P0,HSCALE
      PARAMETER(P0=1000.0E2, HSCALE=6.5E3)
C----------------------------------------------------------------------------
C     Physical constants:
      REAL MRNFIX(MAXLAY),MH2O,MHNO3,MH2SO4,MAIR,
     +     CWV,CNA,CSA,RGAS,G,CP,CPG
      PARAMETER(
     +          MH2O=18.0153E-3,
     +          MHNO3=63.01E-3,
     +          MH2SO4=98.08E-3,
     +          MAIR=28.9644E-3,
     +          RGAS=8.31441,
     +          CWV=MAIR/MH2O,
     +          CNA=MAIR/MHNO3,
     +          CSA=MAIR/MH2SO4,
     +          CP=1004.0E0,
     +          G=9.81E0,
     +          CPG=CP/G)
C     Thermodynamic functions:
      REAL ROSNW,PNANAT,PWVICE
C     Mathematical constants:
      REAL PI
      PARAMETER(PI=3.1415926536E0)
C----------------------------------------------------------------------------
C     Auxilary variables:
      REAL SRATIO(2,MAXLAY),
     +     NDTOT1,MEDIAN1,GSTDEV1,
     +     NDTOT2,MEDIAN2,GSTDEV2,
     +     MRTWV(MAXLAY),MRTNA(MAXLAY),MRTSA(MAXLAY),
     +     NPLSA(MAXLAY),NPFSA(MAXLAY),
     +     NPPSC1(MAXLAY),NPPSC2(MAXLAY),
     +     SALSA(MAXLAY),SAFSA(MAXLAY),
     +     SAPSC1(MAXLAY),SAPSC2(MAXLAY),
     +     MRLSA(MAXLAY),MRFSA(MAXLAY),
     +     MRPSC1(MAXLAY),MRPSC2(MAXLAY),
     +     MASSSA(MAXLAY),MASSNA(MAXLAY),MASSWV(MAXLAY),
     +     RHOTER(MAXLAY),RHOAIR,
     +     X0,X1,X2,X3,X4,
     +     CND,CMCS,CMCN,CMCW
      INTEGER MODE
      LOGICAL ONESTP,VDUMP,DDUMP,AGAIN,QUIT,
     +        TOSC,HOURS,MORE,ERROR_MSG,REDUCE
      DATA MORE/.FALSE./
      INTEGER I,J,K,L
      CHARACTER*100 HEADER
      CHARACTER*100 INPFILE
	CHARACTER*100 INPFILE2
      CHARACTER*100 OUTFILE
      CHARACTER*100 INFFILE
      CHARACTER*100 ERRFILE
      CHARACTER*100 TMPFIL
      DATA DTOUT,DTVTME,DTDTME/3*0.0/
      INTEGER DOTPOS
      REAL DWORK(NBINS)
C----------------------------------------------------------------------------
C     Common block variables:
C
      REAL VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
      REAL TAUMAX
      COMMON /TAU/TAUMAX
C
      LOGICAL SEDI
      COMMON /SED/SEDI
C
      INTEGER TEST
      COMMON /CONTRL/TEST
c
      REAL TABA_CORR
      COMMON /CORRECT/TABA_CORR
	
C	Lengths of INPDIR,OUTDIR, INPFILE, INPFILE2
	INTEGER num1,num2,num3,num4
C----------------------------------------------------------------------------
C     External functions:
      REAL PRSURE,ALTITU,RAMP,MIXRNA,TEMP2D
C----------------------------------------------------------------------------
C     Statement functions:
      REAL PTHETA,TEMP,TETA,THETAP,PPP,HHH,DeltaT
      PTHETA(TETA,TEMP)=P0*(TEMP/TETA)**3.5
      THETAP(PPP,TEMP)=TEMP*(P0/PPP)**0.286
      ALTITU(PPP)=-HSCALE*ALOG(PPP/P0)
      PRSURE(HHH)=P0*EXP(-HHH/HSCALE)
C----------------------------------------------------------------------------
C >>> Program starts here:
C	
	!! Caleb Wherry
	!!
	!! Command Prompt Title & Color
	!!
	CALL system('title MicroPhysical Model')
	CALL system('color 0A')
	!!

      TEST=90
      IF(TEST.GT.0) OPEN(TEST,FILE='TEST.DAT')
      TEST=0
 1    CONTINUE
C----------------------------------------------------------------------------
C     Read input file name and construct output file name:
      CALL INPUTFILE(INPDIR,INPFILE,INPFILE2,MORE)
      OUTFILE=INPFILE
      DOTPOS=INDEX(OUTFILE,'.')

C	!!
C	!! Caleb Wherry
	!!
	!! Changed file system structure. Had to do some
	!! concatenation to form the file names.

	num1 = LEN_TRIM(INPDIR)
	num2 = LEN_TRIM(OUTDIR)
	num3 = LEN_TRIM(INPFILE)
	num4 = LEN_TRIM(INPFILE2)
      
	INPFILE=INPDIR(1:num1)//INPFILE2(1:num4)//'\'//INPFILE(1:num3)
      INFFILE=OUTDIR(1:num2)//INPFILE2(1:num4)//'\Output\'//
     +	OUTFILE(1:DOTPOS)//'inf'
      ERRFILE=OUTDIR(1:num2)//INPFILE2(1:num4)//'\Output\'//
     +	OUTFILE(1:DOTPOS)//'err'
      OUTFILE=OUTDIR(1:num2)//INPFILE2(1:num4)//'\Output\'//
     +	OUTFILE(1:DOTPOS)//'dat'

C     Read input file:
      OPEN(1,FILE=INPFILE, STATUS='OLD')
      READ(1,'(A)') HEADER
      READ(1,*) NLAY
      READ(1,*) RMIN,RMAX
      READ(1,*) DTMAX
      READ(1,'(A1)') TUNIT
      READ(1,*) TSTART
      READ(1,*) TSTOP
      READ(1,*) DTOUT,DTVTME,DTDTME
      READ(1,*) THETA,DPOT
      READ(1,*) TFUNC
      READ(1,*) TMPFIL
      READ(1,*) RAMPT1,TEMP1
      READ(1,*) RAMPT2,TEMP2
      READ(1,*) RAMPT3,TEMP3
      READ(1,*) RAMPT4,TEMP4
      READ(1,*) OMEGA,TAMPLI
      READ(1,*) (MRWVI(L),L=1,IABS(NLAY))
      READ(1,*) (MRNFIX(L),L=1,IABS(NLAY))
      READ(1,*) NDTOT1,NDTOT2
      READ(1,*) MEDIAN1,MEDIAN2
      READ(1,*) GSTDEV1,GSTDEV2
      READ(1,*) LATITU
      READ(1,*) OPTICS
      READ(1,*) MODE
      READ(1,*) TABA_CORR
      CLOSE(1)
C----------------------------------------------------------------------------
C     Test input data:
      OPEN(2,FILE=ERRFILE)
      ERROR_MSG=.FALSE.
      IF(RMIN.GT.MEDIAN1.OR.RMAX.LT.MEDIAN1.OR.RMIN.GT.RMAX) THEN
          WRITE(2,*) 'Mismatch between R-min, R-max, median-1 radius'
          ERROR_MSG=.TRUE.
      ENDIF
      IF(RMIN.GT.MEDIAN2.OR.RMAX.LT.MEDIAN2.OR.RMIN.GT.RMAX) THEN
          WRITE(2,*) 'Mismatch between R-min, R-max, median-2 radius'
          ERROR_MSG=.TRUE.
      ENDIF
      IF(THETA.LT.300.0.AND.IABS(TFUNC).NE.2) THEN
          WRITE(2,*) 'Error in specification pot.temp.',
     +               ' and temperature calculation mode.'
          ERROR_MSG=.TRUE.
      ENDIF
      IF(DPOT.LE.0.0) THEN
          WRITE(2,*) 'Wrong layer thickness specification dPOT'
          ERROR_MSG=.TRUE.
      ENDIF
      IF(ERROR_MSG) THEN
          CLOSE(2,STATUS='KEEP')
          WRITE(*,*) 'Error in input file. Cf. file '
          WRITE(*,*) OUTFILE(1:DOTPOS)//'ERR'
          STOP
      ELSE
          CLOSE(2,STATUS='DELETE')
      ENDIF
      MODE=MAX(2,MIN(MODE,3))
C----------------------------------------------------------------------------
C     Set University of Wyoming OPC radius classes:
      CALL OPC_WYOMING(NCLASS,RADIUSCLASS)
C----------------------------------------------------------------------------
C     Read temperature (pressure) table files if needed:
	deltaT = 0.
      IF(TFUNC.GE.1) THEN
          OPEN(1,FILE=INPDIR//TMPFIL)
          IF(TFUNC.EQ.1) THEN
C             Read time,temperature table values:
              DO I=1,MAXTAB
 31                READ(1,*,IOSTAT=IOS)
     +                  TIMEA(I),(TEMPA(I,L),L=1,IABS(NLAY))
					  Do l = 1, IABS(NLAY)
						Tempa(i,l) = tempa(i,l)+ deltaT
	                  Enddo
                   IF(IOS.GT.0) GOTO 31
                   IF(IOS.LT.0) GOTO 22
              ENDDO
          ELSE IF(TFUNC.EQ.2) THEN
C             Read time,temperature,pressure, thickness table values:
              DO  I=1,MAXTAB
 21                READ(1,*,IOSTAT=IOS)
     +                 TIMEA(I),(TEMPA(I,L),PRESA(I,L),
     +                 L=1,IABS(NLAY))
                   IF(IOS.GT.0) GOTO 21
                   IF(IOS.LT.0) GOTO 22
              ENDDO
          ENDIF
 22       CLOSE(1)
          NPOINT=I-1
      ELSE IF(TFUNC.EQ.0) THEN
          TMPFIL='None'
      ENDIF
C----------------------------------------------------------------------------
C     Write input data into file .INF file about the current simulation:
      OPEN(1,FILE=INFFILE)
      WRITE(1,'(A)') HEADER
      WRITE(1,'(A,A)') 'Input file:  ',INPFILE
      WRITE(1,'(A,A)') 'Output file: ',OUTFILE
      WRITE(1,'(A20,I5,A20,I5)') 'Number of size bins: ',NBINS,
     +                           '  Number of layers:    ',IABS(NLAY)
      write(1,'(A)') 'Temperature/pressure file used: '//TMPFIL
      WRITE(1,'(A,F8.4,A1,F8.4)')
     +           'Minimum and maximum initial radius (microns): ',
     +            RMIN,',',RMAX
      IF(MODE.EQ.2) THEN
          WRITE(1,'(A)') 'Particle surface areas calculated.'
      ELSE IF(MODE.EQ.3) THEN
          WRITE(1,'(A)') 'Particle volumes areas calculated.'
      ENDIF
      WRITE(1,'(A,F8.4)') 'Max. integration time step (s):',DTMAX
      WRITE(1,'(A12,A1)') 'Time units: ',TUNIT
      WRITE(1,'(A,F8.1)') 'Integration start time ',TSTART
      WRITE(1,'(A,F8.1)') 'Integration stop time ',TSTOP
      WRITE(1,'(A,F8.4)') 'Time interval between plot/print output ',
     +            DTOUT
      WRITE(1,'(A,F8.4)') 'Time interval between distribution output ',
     +            DTDTME
      WRITE(1,'(A,F8.4)') 'Time interval between vertical output ',
     +            DTVTME
      IF(TFUNC.EQ.0) THEN
          WRITE(1,'(A)')
     +         'Ramp temperature function applied for all layers.'
          WRITE(1,*) 'Ramp time 1  ',RAMPT1,
     +             '  Temperature (K) ',TEMP1
          WRITE(1,*) 'Ramp time 2  ',RAMPT2,
     +             '  Temperature (K) ',TEMP2
          WRITE(1,*) 'Ramp time 3  ',RAMPT3,
     +             '  Temperature (K) ',TEMP3
          WRITE(1,*) 'Ramp time 4  ',RAMPT4,
     +             '  Temperature (K) ',TEMP4
      ELSE IF(TFUNC.GE.1) THEN
          DO L=1,IABS(NLAY)
              WRITE(1,'(A,I4)') 'Layer: ',L
              IF(TFUNC.EQ.1) THEN
                  WRITE(1,'(A)')
     +                 'Temperature table used: ',INPDIR//TMPFIL
                  WRITE(1,'(A,I6)')'Number of data points read: ',
     +                              NPOINT
                  WRITE(1,*)
     +            'T(',TIMEA(1),')=',TEMPA(1,L)
                  WRITE(1,*)
     +            'T(',TIMEA(NPOINT),')=',TEMPA(NPOINT,L)
              ELSE IF(TFUNC.EQ.2) THEN
                  WRITE(1,'(A)')
     +                 'Temperature/pressure table used: '//TMPFIL
                  WRITE(1,'(A,I6)')'Number of data points read: ',
     +                  NPOINT
                  WRITE(1,*)
     +            'T(',TIMEA(1),')=',TEMPA(1,L)
                  WRITE(1,*)
     +            'T(',TIMEA(NPOINT),')=',TEMPA(NPOINT,L)
                  WRITE(1,*)
     +            'P(',TIMEA(1),')=',PRESA(1,L)
                  WRITE(1,*)
     +            'P(',TIMEA(NPOINT),')=',PRESA(NPOINT,L)
              ENDIF
          END DO
      ENDIF
      IF(OMEGA.NE.0.0.AND.TAMPLI.NE.0.0) THEN
          WRITE(1,'(A,F8.4)') 'Temperature sinus osc. period  ',OMEGA
          WRITE(1,'(A,F8.4)') 'Amplitude (K)',TAMPLI
      ELSE IF(TAMPLI.NE.0.0) THEN
         WRITE(1,'(A,F8.4)') 'Temperature correction (K)',TAMPLI
      ELSE
         WRITE(1,'(A)') 'No temperature correction / sinus oscillations'
      ENDIF
      DO L=1,IABS(NLAY)
          WRITE(1,'(A,I4)') 'Layer: ',L
          WRITE(1,'(A,F8.4)')
     +             'Mixing ratio water vapor, initial value (ppmv) ',
     +                  MRWVI(L)
          IF(MRNFIX(L).GT.0.0) THEN
              WRITE(1,'(A,F8.4)')
     +        'Mixing ratio nitric acid vapor, initial value (ppbv) ',
     +                      MRNFIX(L)
          ELSE
              IF(LATITU.GE.0.0) THEN
                  WRITE(1,*)
     +        'LIMS Northern Hemisphere nitric acid mixing ratio used'
              ELSE IF(LATITU.LT.0.0) THEN
                  WRITE(1,*)
     +        'LIMS Southern Hemisphere nitric acid mixing ratio used'
              ENDIF
          ENDIF
      END DO
      IF(NDTOT1.EQ.0) THEN
          WRITE(1,'(A)')
     +                'SAGE LIQUID sulfate aerosol profile; '
          WRITE(1,'(A,F10.3,A,F10.3)')
     +                'Median (um) = ',MEDIAN1, '    Gstdev= ',GSTDEV1
      ELSE
          WRITE(1,'(A,2F10.3)')
     +                'Number density LIQUID sulfate aerosols (1/ccm) ',
     +                NDTOT1,NDTOT2
          WRITE(1,'(A,F8.4,A,F8.4,A,F8.4,A,F8.4)')
     +                'Median1 (um) = ',MEDIAN1,'  Gstdev1= ',GSTDEV1,
     +                '  Median2 (um) = ',MEDIAN2,'  Gstdev2= ',GSTDEV2
      ENDIF
      WRITE(1,'(A,F10.1)') 'Latitude: ',LATITU
      IF(TFUNC.EQ.0.OR.TFUNC.EQ.1) THEN
          WRITE(1,'(A)') 'Isentropic calculations'
          WRITE(1,'(A,F10.3,A,F10.3)')
     +     'Bottom layer potential temperature (K)',
     +     THETA,' and layer thickness (K) ',DPOT
      ELSE IF(TFUNC.EQ.2) THEN
              WRITE(1,'(A)') 'Pressure/temperature table used.'
              WRITE(1,'(A,F10.3)') 'Top layer thickness (K) ',DPOT
      ENDIF
      IF(NLAY.GE.1) THEN
          WRITE(1,'(A)') 'Sedimentation included'
      ELSE
          WRITE(1,'(A)') 'Sedimentation not included'
      ENDIF
      IF(OPTICS) THEN
          WRITE(1,'(A)') 'Optical calculations included'
      ELSE
          WRITE(1,'(A)') 'Optical calculations not included'
      ENDIF
      WRITE(1,'(A,F10.3)')
     +     'Correction factor for homogeneous NAT nucleation ',
     +            TABA_CORR
      CLOSE(1)
C----------------------------------------------------------------------------
C     Initial settings:
      ONESTP=.FALSE.
      AGAIN=.FALSE.
      QUIT=.FALSE.
      VDUMP=DTVTME.GT.0.0
      DDUMP=DTDTME.GT.0.0
      TOSC=OMEGA.NE.0.0.AND.TAMPLI.NE.0.0
      DO L=1,IABS(NLAY)
          MRWVI(L)=MRWVI(L)*1.0E-6
      END DO
      NLAYER=MAX0(1,MIN0(IABS(NLAY),MAXLAY))
      NDTOT1=NDTOT1*1.0E6
      NDTOT2=NDTOT2*1.0E6
      MEDIAN1=MEDIAN1*1.0E-6
      MEDIAN2=MEDIAN2*1.0E-6
C----------------------------------------------------------------------------
C     Convert input data time settings to seconds:
      IF(INDEX('Hh',TUNIT).GT.0) THEN
          DTOUT=DTOUT*3600.0
          DTVTME=DTVTME*3600.0
          DTDTME=DTDTME*3600.0
          TSTART=TSTART*3600.0
          TSTOP=TSTOP*3600.0
          RAMPT1=RAMPT1*3600.0
          RAMPT2=RAMPT2*3600.0
          RAMPT3=RAMPT3*3600.0
          RAMPT4=RAMPT4*3600.0
          IF(TOSC) OMEGA=2.0*PI/(OMEGA*3600.0)
          IF(TFUNC.GE.1) THEN
              DO I=1,NPOINT
                   TIMEA(I)=TIMEA(I)*3600.0
              ENDDO
          ENDIF
      ELSE IF(INDEX('Dd',TUNIT).GT.0) THEN
          DTOUT=DTOUT*3600.0*24.0
          DTVTME=DTVTME*3600.0*24.0
          DTDTME=DTDTME*3600.0*24.0
          TSTART=TSTART*3600.0*24.0
          TSTOP=TSTOP*3600.0*24.0
          RAMPT1=RAMPT1*3600.0*24.0
          RAMPT2=RAMPT2*3600.0*24.0
          RAMPT3=RAMPT3*3600.0*24.0
          RAMPT4=RAMPT4*3600.0*24.0
          IF(TOSC) OMEGA=2.0*PI/(OMEGA*3600.0*24.0)
          IF(TFUNC.GE.1) THEN
              DO I=1,NPOINT
                   TIMEA(I)=TIMEA(I)*3600.0*24.0
              ENDDO
          ENDIF
      ELSE IF(INDEX('Mm',TUNIT).GT.0) THEN
          DTOUT=DTOUT*60.0
          DTVTME=DTVTME*60.0
          DTDTME=DTDTME*60.0
          TSTART=TSTART*60.0
          TSTOP=TSTOP*60.0
          RAMPT1=RAMPT1*60.0
          RAMPT2=RAMPT2*60.0
          RAMPT3=RAMPT3*60.0
          RAMPT4=RAMPT4*60.0
          IF(TOSC) OMEGA=2.0*PI/(OMEGA*60.0)
          IF(TFUNC.GE.1) THEN
              DO I=1,NPOINT
                   TIMEA(I)=TIMEA(I)*60.0
              ENDDO
          ENDIF
      ENDIF
      DO L=1,NLAYER
          IF(TFUNC.GE.1) THEN
              CALL SPLINE(TIMEA,TEMPA(1,L),NPOINT,
     +                    YWORK,Y2TEMP(1,L))
              IF(TFUNC.EQ.2) THEN
                  CALL SPLINE(TIMEA,PRESA(1,L),NPOINT,
     +                              YWORK,Y2PRES(1,L))
              ENDIF
              TSTART=AMAX1(TSTART,TIMEA(1))
              TSTOP=AMIN1(TSTOP,TIMEA(NPOINT))
          ENDIF
      END DO
C----------------------------------------------------------------------------
C     Initial time setting:
C----------------------------------------------------------------------------
      TIME=TSTART
C----------------------------------------------------------------------------
C     Initial settings for each layer:
C----------------------------------------------------------------------------
      DO L=1,NLAYER
C
C         Set initial temperature:
C         ------------------------
          TAIR(L)=0.0
          IF(TFUNC.GE.1) THEN
              CALL SPLINT(TIMEA,TEMPA(1,L),Y2TEMP(1,L),
     +                    NPOINT,TIME,TAIR(L))
          ELSE IF(TFUNC.EQ.0) THEN
              TAIR(L)=RAMP(TIME,
     +                  RAMPT1,TEMP1,
     +                  RAMPT2,TEMP2,
     +                  RAMPT3,TEMP3,
     +                  RAMPT4,TEMP4)
          ENDIF
C
C         Set initial pressure:
C         ---------------------
          IF(TFUNC.EQ.0.OR.TFUNC.EQ.1) THEN
              POT(L)=THETA+(NLAYER-L)*DPOT
              PRESS(L)=PTHETA(POT(L),TAIR(L))
              HEIGHT(L)=ALTITU(PRESS(L))
          ELSE IF(TFUNC.EQ.2) THEN
              CALL SPLINT(TIMEA,PRESA(1,L),Y2PRES(1,L),
     +                        NPOINT,TIME,PRESS(L))
              HEIGHT(L)=ALTITU(PRESS(L))
              POT(L)=THETAP(PRESS(L),TAIR(L))
          ENDIF
          PAIR(L)=PRESS(L)
C
C         Set initial layer thickness:
C         ----------------------------
          RHOAIR=MAIR*PAIR(L)/(RGAS*TAIR(L))
          IF(L.EQ.1) THEN
              THICKNESS(1)=CPG*TAIR(1)*DPOT/POT(1)
          ELSE
              THICKNESS(L)=(PAIR(L)-PAIR(L-1))/(G*RHOAIR)
          ENDIF
C
C         Initial value of water vapor mixing ratio:
C         ------------------------------------------
          MRWV(L)=MRWVI(L)
C
C         Initial values of nitric acid mixing ratio:
C         -------------------------------------------
          IF(MRNFIX(L).GT.0.0) THEN
C             Fixed nitric acid mixing ratio at all levels:
              MRNA(L)=MRNFIX(L)*1.0E-9
          ELSE
C             LIMS profile:
              MRNA(L)=MIXRNA(HEIGHT(L),LATITU)
          ENDIF
          MRNAI(L)=MRNA(L)
C
C         Initial partial pressures:
C         --------------------------
          PPWV=PAIR(L)*MRWV(L)
          PPNA=PAIR(L)*MRNA(L)
C
          CALL PSCBOX_START(
     +        RMIN,RMAX,
     +        NDTOT1,MEDIAN1,GSTDEV1,
     +        NBINS,TYPES,
     +        TAIR(L),PAIR(L),PPWV,PPNA,
     +        ND(1,1,L),
     +        MCS(1,1,L),MCN(1,1,L),MCW(1,1,L),
     +        SDND,SDMCS,SDMCN,SDMCW,
     +        PRESENCE(1,L),SEDMNT,
     +        PTSIZE)
C
          IF(NDTOT2.GT.0.0) THEN
              NDTOT2=NDTOT2*RGAS*TAIR(L)/(MAIR*PAIR(L))
              CALL LGNDST(NBINS,NDTOT2,MEDIAN2,GSTDEV2,PTSIZE,WORK(1,1))
              DO I=1,NBINS
                  ND(I,STS,L)=ND(I,STS,L)+WORK(I,1)
              END DO
          ENDIF
C
          SRATIO(1,L)=PPNA/PNANAT(TAIR(L),PPWV)
          SRATIO(2,L)=PPWV/PWVICE(TAIR(L))
C
C         Calculate total gas phase volume mixing ratios
C         and condensed phase mass mixing ratios:
C         ----------------------------------------------
          CALL MIXCON(NBINS,TYPES,
     +                ND(1,1,L),
     +                MCS(1,1,L),MCN(1,1,L),MCW(1,1,L),
     +                MASSWV(L),MASSNA(L),MASSSA(L))
          MRTWV(L)=MASSWV(L)*CWV+MRWV(L)
          MRTNA(L)=MASSNA(L)*CNA+MRNA(L)
          MRTSA(L)=MASSSA(L)*CSA
C
C         Optical cross section calcutaions:
C         ----------------------------------
          IF(OPTICS) THEN
              DO K=1,N_WAVE
                  DO J=1,TYPES
                      IF(J.EQ.STS) THEN
                          DO I=1,NBINS
                             CALL MIEBCK(WAVELENGTH(K),PTSIZE(I,1),
     +                                   ANGLE(K),
     +                                   REFRACTIVE_INDEX(J,K),1.0E-8,
     +                                   BACKSCAT_CROSS(I,J,K))
                             CALL MIEEXT(WAVELENGTH(K),PTSIZE(I,1),
     +                                   REFRACTIVE_INDEX(J,K),1.0E-8,
     +                                   EXTINCTION_CROSS(I,J,K))
                             DEPOLARIZATION(I,J,K)=0.0
                          END DO
                      ELSE
                          CALL T_MATRIX(T_MATRIX_DIR,NBINS,PTSIZE(1,1),
     +                                  ASPECT(J),
     +                                  REFRACTIVE_INDEX(J,K),
     +                                  WAVELENGTH(K),ANGLE(K),
     +                                  BACKSCAT_CROSS(1,J,K),
     +                                  EXTINCTION_CROSS(1,J,K),
     +                                  DEPOLARIZATION(1,J,K))
                      ENDIF
                  END DO
              END DO
c
              RHOAIR=MAIR*PAIR(L)/(RGAS*TAIR(L))
              DO K=1,N_WAVE
                  RAY=RAYLGH(WAVELENGTH(K),ANGLE(K),PAIR(L),TAIR(L))
                  AEROSOL_BACKSCATTER_RATIO(K,L)=0.0
                  EXTINCTION(K,L)=0.0
                  DEPOL(K,L)=0.0
                  DO J=1,TYPES
                      IF(PRESENCE(J,L)) THEN
                          DO I=1,NBINS
                              AEROSOL_BACKSCATTER_RATIO(K,L)=
     +                            AEROSOL_BACKSCATTER_RATIO(K,L)+
     +                            ND(I,J,L)*BACKSCAT_CROSS(I,J,K)
                              EXTINCTION(K,L)=EXTINCTION(K,L)+
     +                            ND(I,J,L)*EXTINCTION_CROSS(I,J,K)
                              DEPOL(K,L)=DEPOL(K,L)+
     +                            ND(I,J,L)*DEPOLARIZATION(I,J,K)
                          END DO
                      ENDIF
                  END DO
c                 Aerosol depolarisation:
                  DEPOL(K,L)=DEPOL(K,L)/AEROSOL_BACKSCATTER_RATIO(K,L)
                  AEROSOL_BACKSCATTER_RATIO(K,L)=RHOAIR*
     +                    AEROSOL_BACKSCATTER_RATIO(K,L)/RAY
c                 Volume depolarisation:
cc                DEPOL(K,L)=(DEPOL(K,L)*
cc   +                     AEROSOL_BACKSCATTER_RATIO(K,L)+DEPOL_M)/
cc   +                      (AEROSOL_BACKSCATTER_RATIO(K,L)+1.0)
ccc                  X0=AEROSOL_BACKSCATTER_RATIO(K,L)
ccc                  AEROSOL_BACKSCATTER_RATIO(K,L)= 
ccc     +				AEROSOL_BACKSCATTER_RATIO(K,L)
				x0 = RHOAIR*
     +                    AEROSOL_BACKSCATTER_RATIO(K,L)/RAY
                  X1=DEPOL(K,L)
cc                  DEPOL(K,L)=(X1*((X0+1.0)*DEPOL_M+X0)+DEPOL_M)/
cc    +                       (DEPOL_M*X0+X0+1.0+X1)
                  DEPOL(K,L)=DEPOL(K,L)*100.0
                  EXTINCTION(K,L)=EXTINCTION(K,L)*RHOAIR
              END DO
          ELSE
              DO K=1,N_WAVE
                  AEROSOL_BACKSCATTER_RATIO(K,L)=1.0
                  DEPOL(K,L)=1.0
                  EXTINCTION(K,L)=1.0
              ENDDO
          ENDIF
      ENDDO
C----------------------------------------------------------------------------
C     Output time settings:
C----------------------------------------------------------------------------
      IF(ONESTP) THEN
          PTIME=TIME+DTMAX
      ELSE
          PTIME=TIME+DTOUT
      ENDIF
      IF(DDUMP) THEN
          DTIME=TIME+DTDTME
      ELSE
          DTIME=TSTOP
      ENDIF
      IF(VDUMP) THEN
          VTIME=TIME+DTVTME
      ELSE
          VTIME=TSTOP
      ENDIF
C----------------------------------------------------------------------------
C     Initial integration step size:
C----------------------------------------------------------------------------
      DT=AMIN1(DTMAX,PTIME-TIME)
C----------------------------------------------------------------------------
C     Calculation of gross number densities:
C----------------------------------------------------------------------------
      CALL GROSS(
     +        NBINS,TYPES,NCLASS,MAXLAY,NLAYER,MODE,
     +        TAIR,PAIR,PRESENCE,
     +        ND,
     +        MCS,MCN,MCW,
     +        MRWV,MRNA,
     +        PTSIZE,
     +        MRTWV,MRTNA,MRTSA,
     +        NPLSA,NPFSA,NPPSC1,NPPSC2,
     +        SALSA,SAFSA,SAPSC1,SAPSC2,
     +        MRLSA,MRFSA,MRPSC1,MRPSC2,
     +        MASSSA,MASSNA,MASSWV,
     +        WSA,WNA,
     +        SIZEDIST,RADIUSCLASS)
C----------------------------------------------------------------------------
C     Initial display:
C     Interactive routine:
C
      CALL IA(NBINS,TYPES,MAXLAY,NLAYER,MODE,
     +         TIME,TSTOP,TUNIT,TAIR,PAIR,
     +         ND,MRWV,MRNA,
     +         SRATIO,PRESENCE,ONESTP,AGAIN,QUIT,VDUMP,DDUMP,
     +         PTSIZE,MRTWV,MRTNA,NPLSA,NPFSA,NPPSC1,NPPSC2,
     +         SALSA,SAFSA,SAPSC1,SAPSC2,MRLSA,MRFSA,MRPSC1,MRPSC2,
     +         WSA,WNA,DWORK,WORK(1,1),HEADER,INPFILE)
      IF(QUIT) STOP
C----------------------------------------------------------------------------
C     Initial data output on file:
      IF(DTOUT.LE.TSTOP) CALL DUMP(OUTFILE,
     +        NBINS,TYPES,NCLASS,MAXLAY,NLAYER,MODE,
     +        TIME,TUNIT,
     +        TAIR,
     +        HEIGHT,PAIR,
     +        ND,
     +        MRWV,MRNA,MRWVI,MRNAI,
     +        SRATIO,
     +        AGAIN,QUIT,VDUMP,DDUMP,
     +        PTSIZE,
     +        MRTWV,MRTNA,MRTSA,WSA,WNA,
     +        NPLSA,NPFSA,NPPSC1,NPPSC2,
     +        SALSA,SAFSA,SAPSC1,SAPSC2,
     +        MRLSA,MRFSA,MRPSC1,MRPSC2,
     +        N_WAVE,AEROSOL_BACKSCATTER_RATIO,EXTINCTION,DEPOL,
     +        SIZEDIST,WORK(1,1))
C----------------------------------------------------------------------------
C
C >>> Main loop for column integration starts here:
C
C----------------------------------------------------------------------------
 2    CONTINUE
C----------------------------------------------------------------------------
C     Set sedimentation calculation flag:
      SEDI=NLAY.GE.1
C----------------------------------------------------------------------------
C     Initialize sedimentation flow at top layer:
      IF(SEDI) THEN
          DO J=1,TYPES
              DO I=1,NBINS
                  SDND(I,J)=0.0
                  SDMCS(I,J)=0.0
                  SDMCN(I,J)=0.0
                  SDMCW(I,J)=0.0
              ENDDO
              SEDMNT(J)=.FALSE.
          ENDDO
      ENDIF
C----------------------------------------------------------------------------
C     Calculate temperature, presssure, and mixing raios for each layer:
C----------------------------------------------------------------------------
      DO L=1,NLAYER
C
C         Set temperature:
C         ----------------
          TAIR(L)=0.0
          IF(TFUNC.GE.1) THEN
C             Cubic spline interpolation for temperature in 1-D table:
              CALL SPLINT(TIMEA,TEMPA(1,L),Y2TEMP(1,L),
     +                    NPOINT,TIME,TAIR(L))
          ELSE IF (TFUNC.EQ.0) THEN
              TAIR(L)=RAMP(TIME,
     +             RAMPT1,TEMP1,
     +             RAMPT2,TEMP2,
     +             RAMPT3,TEMP3,
     +             RAMPT4,TEMP4)
          ENDIF
C
C         Temparature oscillation:
C         ------------------------
          DELT=TAMPLI
          IF(TOSC) DELT=TAMPLI*SIN(OMEGA*(TIME-TSTART))
C
C         Set pressure and add temperature oscillations:
C         ----------------------------------------------
          IF(TFUNC.EQ.0.OR.TFUNC.EQ.1) THEN
C              Isentropic pressure:
               TAIR(L)=TAIR(L)+DELT
               PAIR(L)=PTHETA(POT(L),TAIR(L))
               HEIGHT(L)=ALTITU(PAIR(L))
          ELSE IF(TFUNC.EQ.2) THEN
C              Cubic spline interpolation for pressure:
               CALL SPLINT(TIMEA,PRESA(1,L),Y2PRES(1,L),NPOINT,
     +                     TIME,PAIR(L))
               TAIR(L)=TAIR(L)+DELT
               POT(L)=THETAP(PAIR(L),TAIR(L))
               HEIGHT(L)=ALTITU(PAIR(L))
          ENDIF
C
C         Set layer thickness:
C         --------------------
          RHOAIR=MAIR*PAIR(L)/(RGAS*TAIR(L))
          IF(L.EQ.1) THEN
              THICKNESS(1)=CPG*TAIR(1)*DPOT/POT(1)
          ELSE
              THICKNESS(L)=(PAIR(L)-PAIR(L-1))/(G*RHOAIR)
          ENDIF
      ENDDO
C----------------------------------------------------------------------------
c     Calculate time step:
C----------------------------------------------------------------------------
      DT=AMAX1(DTMIN,AMIN1(DTMAX,PTIME-TIME))
      REDUCE=.TRUE.
      DO WHILE(REDUCE.AND.DT.GT.DTMIN)
          REDUCE=.FALSE.
          DO L=1,NLAYER
              IF(TFUNC.GE.1) THEN
                  CALL SPLINT(TIMEA,TEMPA(1,L),Y2TEMP(1,L),
     +                        NPOINT,TIME+DT,TEMPE)
C              ELSE IF(TFUNC.EQ.-2) THEN
C                  TEMPE=TEMP2D(TIME+DT,TSTART,L,LATITU,
C     +                  INPDIR//TMPFIL,TUNIT)
              ELSE IF (TFUNC.EQ.0) THEN
                  TEMPE=RAMP(TIME+DT,
     +             RAMPT1,TEMP1,
     +             RAMPT2,TEMP2,
     +             RAMPT3,TEMP3,
     +             RAMPT4,TEMP4)
              ENDIF
              DELT=TAMPLI
              IF(TOSC) DELT=TAMPLI*SIN(OMEGA*(TIME+DT-TSTART))
              TEMPE=TEMPE+DELT
              DTEMP=ABS(TEMPE-TAIR(L))
              REDUCE=REDUCE.OR.DTEMP.GT.DTEMPMAX
          ENDDO
          IF(REDUCE) DT=DT/2.0
      ENDDO
C
C----------------------------------------------------------------------------
C     Calculate PSC particle size distribution for each layer:
C----------------------------------------------------------------------------
      DO L=1,NLAYER
C
C         Sedimentation flags:
C         --------------------
          IF(.NOT.SEDI) THEN
              DO J=1,TYPES
                 SEDMNT(J)=.FALSE.
              END DO
          ENDIF
C
C         Set partial  pressures:
C         -----------------------
          PPWV=PAIR(L)*MRWV(L)
          PPNA=PAIR(L)*MRNA(L)
c
          CALL PSCBOX(
     +        DT,
     +        THICKNESS(L),
     +        NBINS,TYPES,NWORK,
     +        TAIR(L),PAIR(L),PPWV,PPNA,
     +        ND(1,1,L),
     +        MCS(1,1,L),MCN(1,1,L),MCW(1,1,L),
     +        SDND,SDMCS,SDMCN,SDMCW,
     +        PRESENCE(1,L),SEDMNT,
     +        PTSIZE(1,1),WORK)
C
          MRWV(L)=PPWV/PAIR(L)
          MRNA(L)=PPNA/PAIR(L)
C
C         Store common block variables:
C         -----------------------------
          SRATIO(1,L)=PPNA/PNANAT(TAIR(L),PPWV)
          SRATIO(2,L)=PPWV/PWVICE(TAIR(L))

C----------------------------------------------------------------------------
C     Column loop ends here:
      ENDDO
C----------------------------------------------------------------------------
C     Actual time:
      TIME=TIME+DT
C----------------------------------------------------------------------------
C     Has integration stop time been reached ?
      QUIT=TIME.GE.TSTOP
C----------------------------------------------------------------------------
      IF(QUIT.OR.TIME.GE.AMIN1(PTIME,DTIME,VTIME)) THEN
C         Time for output etc.:
C----------------------------------------------------------------------------
C         Calculation of gross number densities:
          CALL GROSS(
     +        NBINS,TYPES,NCLASS,MAXLAY,NLAYER,MODE,
     +        TAIR,PAIR,PRESENCE,
     +        ND,
     +        MCS,MCN,MCW,
     +        MRWV,MRNA,
     +        PTSIZE,
     +        MRTWV,MRTNA,MRTSA,
     +        NPLSA,NPFSA,NPPSC1,NPPSC2,
     +        SALSA,SAFSA,SAPSC1,SAPSC2,
     +        MRLSA,MRFSA,MRPSC1,MRPSC2,
     +        MASSSA,MASSNA,MASSWV,
     +        WSA,WNA,
     +        SIZEDIST,RADIUSCLASS)
C----------------------------------------------------------------------------
C         Interactive routine:
          CALL IA(NBINS,TYPES,MAXLAY,NLAYER,MODE,
     +         TIME,TSTOP,TUNIT,TAIR,PAIR,
     +         ND,MRWV,MRNA,
     +         SRATIO,PRESENCE,ONESTP,AGAIN,QUIT,VDUMP,DDUMP,
     +         PTSIZE,MRTWV,MRTNA,NPLSA,NPFSA,NPPSC1,NPPSC2,
     +         SALSA,SAFSA,SAPSC1,SAPSC2,MRLSA,MRFSA,MRPSC1,MRPSC2,
     +         WSA,WNA,DWORK,WORK(1,1),HEADER,INPFILE)
C----------------------------------------------------------------------------
C         Optical calcutaions:
          IF(OPTICS) THEN
              DO L=1,NLAYER

C	!!
cmp   !! Niels modification to account for STS refractive indices
	!!

                  J=STS
                  DO K=1,N_WAVE
                     CALL STS_INDEX(WAVELENGTH(K),TAIR(L),WSA(L),WNA(L),
     +                              REFRACTIVE_INDEX(J,K))
                     DO I=1,NBINS
                         CALL MIEBCK(WAVELENGTH(K),PTSIZE(I,1),
     +                                   ANGLE(K),
     +                                   REFRACTIVE_INDEX(J,K),1.0E-8,
     +                                   BACKSCAT_CROSS(I,J,K))
                         CALL MIEEXT(WAVELENGTH(K),PTSIZE(I,1),
     +                                   REFRACTIVE_INDEX(J,K),1.0E-8,
     +                                   EXTINCTION_CROSS(I,J,K))
                         DEPOLARIZATION(I,J,K)=0.0
                     END DO
                  END DO

	!!
cmp	!! 	End of Niels modifications
C	!!

                  RHOAIR=MAIR*PAIR(L)/(RGAS*TAIR(L))
                  DO K=1,N_WAVE
                      RAY=RAYLGH(WAVELENGTH(K),ANGLE(K),PAIR(L),TAIR(L))
                      AEROSOL_BACKSCATTER_RATIO(K,L)=0.0
                      EXTINCTION(K,L)=0.0
                      DEPOL(K,L)=0.0
                      DO J=1,TYPES
                          IF(PRESENCE(J,L)) THEN
                              DO I=1,NBINS
                                  AEROSOL_BACKSCATTER_RATIO(K,L)=
     +                                AEROSOL_BACKSCATTER_RATIO(K,L)+
     +                                ND(I,J,L)*BACKSCAT_CROSS(I,J,K)
                                  EXTINCTION(K,L)=EXTINCTION(K,L)+
     +                                ND(I,J,L)*EXTINCTION_CROSS(I,J,K)
                                  DEPOL(K,L)=DEPOL(K,L)+
     +                                ND(I,J,L)*DEPOLARIZATION(I,J,K)
                              END DO
                          ENDIF
                      END DO
                      DEPOL(K,L)=DEPOL(K,L)/
     +                           AEROSOL_BACKSCATTER_RATIO(K,L)
                      AEROSOL_BACKSCATTER_RATIO(K,L)=RHOAIR*
     +                        AEROSOL_BACKSCATTER_RATIO(K,L)/RAY
C
C I uncommented the next three lines and commented out the following three lines. I assume that this switches between 
C Volume and Aerosol Depolarization
C 
cc                      DEPOL(K,L)=(DEPOL(K,L)*
cc     +                         AEROSOL_BACKSCATTER_RATIO(K,L)+DEPOL_M)/
cc     +                          (AEROSOL_BACKSCATTER_RATIO(K,L)+1.0)
cc                      X0=AEROSOL_BACKSCATTER_RATIO(K,L)
cc                      X1=DEPOL(K,L)
cc                      DEPOL(K,L)=(X1*((X0+1.0)*DEPOL_M+X0)+DEPOL_M)/
cc     +                           (DEPOL_M*X0+X0+1.0+X1)
                      DEPOL(K,L)=DEPOL(K,L)*100.0
                      EXTINCTION(K,L)=EXTINCTION(K,L)*RHOAIR
                  END DO
              ENDDO
          ENDIF
C----------------------------------------------------------------------------
C         Test and adjust time for vertical column output:
          IF(TIME.GE.VTIME.AND.DTVTME.GT.0.0) THEN
              VDUMP=.TRUE.
              VTIME=AMIN1(VTIME+DTVTME,TSTOP)
          ENDIF
C----------------------------------------------------------------------------
C         Test and adjust time for particle distribution output:
          IF((TIME.GE.DTIME.AND.DTDTME.GT.0.0).OR.QUIT) THEN
              DDUMP=.TRUE.
              DTIME=AMIN1(DTIME+DTDTME,TSTOP)
          ENDIF
C----------------------------------------------------------------------------
C         Output data on file:
          CALL DUMP(OUTFILE,
     +        NBINS,TYPES,NCLASS,MAXLAY,NLAYER,MODE,
     +        TIME,TUNIT,
     +        TAIR,
     +        HEIGHT,PAIR,
     +        ND,
     +        MRWV,MRNA,MRWVI,MRNAI,
     +        SRATIO,
     +        AGAIN,QUIT,VDUMP,DDUMP,
     +        PTSIZE,
     +        MRTWV,MRTNA,MRTSA,WSA,WNA,
     +        NPLSA,NPFSA,NPPSC1,NPPSC2,
     +        SALSA,SAFSA,SAPSC1,SAPSC2,
     +        MRLSA,MRFSA,MRPSC1,MRPSC2,
     +        N_WAVE,AEROSOL_BACKSCATTER_RATIO,EXTINCTION,DEPOL,
     +        SIZEDIST,WORK(1,1))
C----------------------------------------------------------------------------
C         Output column data on file:
CC        CALL COLUMN(
CC   +        TIME,TUNIT,
CC   +        NBINS,TYPES,NLAYER,
CC   +        TAIR,
CC   +        PAIR,
CC   +        MRNA,MRWV,
CC   +        ND,
CC   +        MCS,MCN,MCW,
CC   +        SDND,SDMCS,SDMCN,SDMCW,
CC   +        THICKNESS,
CC   +        CND,CMCS,CMCN,CMCW)
C----------------------------------------------------------------------------
C         Start all over again ?
          IF(AGAIN) GOTO 1
C----------------------------------------------------------------------------
C         Update next print/plot time:
          IF(ONESTP) THEN
              PTIME=PTIME+DTMAX
          ELSE
              PTIME=PTIME+DTOUT
          ENDIF
C----------------------------------------------------------------------------
      ENDIF
C----------------------------------------------------------------------------
C     Calculate the next time step:
      IF(.NOT.QUIT) GOTO 2
C----------------------------------------------------------------------------
      IF(QUIT.AND.MORE) THEN
          AGAIN=.TRUE.
          GOTO 1
      ENDIF
C----------------------------------------------------------------------------
C     End of main program:
      STOP
 1000 FORMAT(10(1PE9.2,1X),0P)
 1010 FORMAT(1X,8(A9))
      END
C
C******************************************************************************
C     SUBROUTINE GROSS
C
C     This subroutine calculates integral parameters of size distributions:
C
C     Total volume mixing ratios of
C          nitric acid (MRTNA)
C          water (MRTWV)
C          sulfuric acid (MRTSA)
C     Number densities (#/m**3) of
C          liquid type 1b PSC particles (NPLSA)
C          frozen sulfate aerosols (NPFSA)
C          type 1a PSC (NPPSC1)
C          type 2 PSC (NPPSC2)
C     Volume density (m**3/m**3) or surface area density  (m**2/m**3) of
C          liquid type 1b PSC particles (SALSA)
C          frozen sulfate aerosols (SAFSA)
C          type 1a PSC (SAPSC1)
C          type 2 PSC (SAPSC2)
C     Mean radius (m) of
C          liquid type 1b PSC particles (MRLSA)
C          frozen sulfate aerosols (MRFSA)
C          type 1a PSC (MRPSC1)
C          type 2 PSC (MRPSC2)
C     Condensed phase mass mixing ratios of
C          nitric acid (MASSNA)
C          water (MASSWV)
C          sulfuric acid (MASSSA)
C******************************************************************************
C
      SUBROUTINE GROSS(
     +        NBINS,TYPES,NCLASS,MAXLAY,NLAYER,MODE,
     +        TAIR,PAIR,PRESENCE,
     +        ND,
     +        MCS,MCN,MCW,
     +        MRWV,MRNA,
     +        PTSIZE,
     +        MRTWV,MRTNA,MRTSA,
     +        NPLSA,NPFSA,NPPSC1,NPPSC2,
     +        SALSA,SAFSA,SAPSC1,SAPSC2,
     +        MRLSA,MRFSA,MRPSC1,MRPSC2,
     +        MASSSA,MASSNA,MASSWV,
     +        WSA,WNA,
     +        SIZEDIST,RADIUSCLASS)
C
      IMPLICIT NONE
C     Model size:
      INTEGER NBINS,TYPES,MAXLAY,NLAYER,MODE,I,NCLASS
C
      REAL
     +        TAIR(MAXLAY),PAIR(MAXLAY),
     +        ND(NBINS,TYPES,MAXLAY),
     +        MCS(NBINS,TYPES,MAXLAY),
     +        MCN(NBINS,TYPES,MAXLAY),
     +        MCW(NBINS,TYPES,MAXLAY),
     +        MRWV(MAXLAY),MRNA(MAXLAY),
     +        PTSIZE(NBINS,8),
     +        MRTWV(MAXLAY),MRTNA(MAXLAY),MRTSA(MAXLAY),
     +        NPLSA(MAXLAY),NPFSA(MAXLAY),
     +        NPPSC1(MAXLAY),NPPSC2(MAXLAY),
     +        SALSA(MAXLAY),SAFSA(MAXLAY),
     +        SAPSC1(MAXLAY),SAPSC2(MAXLAY),
     +        MRLSA(MAXLAY),MRFSA(MAXLAY),
     +        MRPSC1(MAXLAY),MRPSC2(MAXLAY),
     +        MASSSA(MAXLAY),MASSNA(MAXLAY),MASSWV(MAXLAY),
     +        WSA(MAXLAY),WNA(MAXLAY),
     +        SIZEDIST(NCLASS,MAXLAY),RADIUSCLASS(NCLASS)
C
      LOGICAL PRESENCE(TYPES,MAXLAY)
C     Particle types:
      INTEGER STS,SAT,NAT,ICE
      PARAMETER(
C     Liquid supercooled ternary solution particle (SA or Type 1b PSC)
     +  STS=1,
C     Solid sulfuric acid tetrahydrate (SAT) particle; i.e. no HNO3 or excess ice
     +  SAT=2,
C     Solid nitric acid trihydrate (NAT) particle; i.e. NAT and SAT but no excess ice (Type 1a PSC)
     +  NAT=3,
C     Solid ice particle; i.e. holding excess ice, NAT, and SAT (Type 2 PSC)
     +  ICE=4)
C
      REAL MAIR,MH2O,MHNO3,MH2SO4,C4,CWV,CNA,CSA,RGAS,RHOAIR,X0
      PARAMETER(
     +          RGAS=8.31441,
     +          MH2O=18.0153E-3,
     +          MHNO3=63.01E-3,
     +          MH2SO4=98.08E-3,
     +          MAIR=28.9644E-3,
     +          C4=MAIR/RGAS,
     +          CWV=MAIR/MH2O,
     +          CNA=MAIR/MHNO3,
     +          CSA=MAIR/MH2SO4)
      INTEGER J,L
C
      REAL VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
      DO L=1,NLAYER
C----------------------------------------------------------------------------
          RHOAIR=C4*PAIR(L)/TAIR(L)
C----------------------------------------------------------------------------
C         Calculate total number density, mean radius and surface area
C         of distributions:
          NPLSA(L)=0.0
          MRLSA(L)=0.0
          SALSA(L)=0.0
          IF(PRESENCE(STS,L)) THEN
              DO I=1,NBINS
                  NPLSA(L)=NPLSA(L)+ND(I,STS,L)
                  SALSA(L)=SALSA(L)+PTSIZE(I,MODE)*ND(I,STS,L)
                  MRLSA(L)=MRLSA(L)+ND(I,STS,L)*PTSIZE(I,1)
              ENDDO
              IF(NPLSA(L).GT.ZERO) THEN
                  MRLSA(L)=MRLSA(L)/NPLSA(L)
              ELSE
                  MRLSA(L)=0.0
              ENDIF
          ENDIF
          NPFSA(L)=0.0
          MRFSA(L)=0.0
          SAFSA(L)=0.0
          IF(PRESENCE(SAT,L)) THEN
              DO I=1,NBINS
                  NPFSA(L)=NPFSA(L)+ND(I,SAT,L)
                  SAFSA(L)=SAFSA(L)+PTSIZE(I,MODE)*ND(I,SAT,L)
                  MRFSA(L)=MRFSA(L)+ND(I,SAT,L)*PTSIZE(I,1)
              ENDDO
              IF(NPFSA(L).GT.ZERO) THEN
                  MRFSA(L)=MRFSA(L)/NPFSA(L)
              ELSE
                  MRFSA(L)=0.0
              ENDIF
          ENDIF
          NPPSC1(L)=0.0
          MRPSC1(L)=0.0
          SAPSC1(L)=0.0
          IF(PRESENCE(NAT,L)) THEN
              DO I=1,NBINS
                  NPPSC1(L)=NPPSC1(L)+ND(I,NAT,L)
                  SAPSC1(L)=SAPSC1(L)+PTSIZE(I,MODE)*ND(I,NAT,L)
                  MRPSC1(L)=MRPSC1(L)+ND(I,NAT,L)*PTSIZE(I,1)
              ENDDO
              IF(NPPSC1(L).GT.ZERO) THEN
                  MRPSC1(L)=MRPSC1(L)/NPPSC1(L)
              ELSE
                  MRPSC1(L)=0.0
              ENDIF
          ENDIF
          NPPSC2(L)=0.0
          MRPSC2(L)=0.0
          SAPSC2(L)=0.0
          IF(PRESENCE(ICE,L)) THEN
              DO I=1,NBINS
                  NPPSC2(L)=NPPSC2(L)+ND(I,ICE,L)
                  SAPSC2(L)=SAPSC2(L)+PTSIZE(I,MODE)*ND(I,ICE,L)
                  MRPSC2(L)=MRPSC2(L)+ND(I,ICE,L)*PTSIZE(I,1)
              ENDDO
              IF(NPPSC2(L).GT.ZERO) THEN
                  MRPSC2(L)=MRPSC2(L)/NPPSC2(L)
              ELSE
                  MRPSC2(L)=0.0
              ENDIF
          ENDIF
          NPLSA(L)=NPLSA(L)*RHOAIR
          NPFSA(L)=NPFSA(L)*RHOAIR
          NPPSC1(L)=NPPSC1(L)*RHOAIR
          NPPSC2(L)=NPPSC2(L)*RHOAIR
          SALSA(L)=SALSA(L)*RHOAIR
          SAFSA(L)=SAFSA(L)*RHOAIR
          SAPSC1(L)=SAPSC1(L)*RHOAIR
          SAPSC2(L)=SAPSC2(L)*RHOAIR
          IF(NPLSA(L).LT.1.0E-2) MRLSA(L)=0.0
          IF(NPFSA(L).LT.1.0E-2) MRFSA(L)=0.0
          IF(NPPSC1(L).LT.1.0E-2) MRPSC1(L)=0.0
          IF(NPPSC2(L).LT.1.0E-2) MRPSC2(L)=0.0
C----------------------------------------------------------------------------
C         Calculate total gas phase volume mixing ratios
C         and condensed phase mass mixing ratios:
          CALL MIXCON(NBINS,TYPES,
     +                ND(1,1,L),
     +                MCS(1,1,L),MCN(1,1,L),MCW(1,1,L),
     +                MASSWV(L),MASSNA(L),MASSSA(L))
          MRTWV(L)=MASSWV(L)*CWV+MRWV(L)
          MRTNA(L)=MASSNA(L)*CNA+MRNA(L)
          MRTSA(L)=MASSSA(L)*CSA
C----------------------------------------------------------------------------
C         Calculate acid weight fractions:
          X0=MASSWV(L)+MASSNA(L)+MASSSA(L)
          WSA(L)=MASSSA(L)/X0
          WNA(L)=MASSNA(L)/X0
C----------------------------------------------------------------------------
C         Calculate cumulated size distribution n(r>ri):
          DO J=1,NCLASS
              SIZEDIST(J,L)=0.0
          END DO
          DO I=NBINS,1,-1
             DO J=1,NCLASS
                 IF(PTSIZE(I,1).GE.RADIUSCLASS(J))
     +               SIZEDIST(J,L)=SIZEDIST(J,L)+
     +               ND(I,STS,L)+ND(I,SAT,L)+ND(I,NAT,L)+ND(I,ICE,L)
             END DO
          END DO
          DO J=1,NCLASS
              SIZEDIST(J,L)=SIZEDIST(J,L)*RHOAIR
          END DO
      ENDDO
      RETURN
      END
C
C******************************************************************************
C     SUBROUTINE DUMP
C     This subroutine dumps model results on disk files.
C******************************************************************************
*
*
*     Output on files:
*     ----------------
*
*
*     Three types of ASCII output files can be created by the program:
*
*         1:    Historical output files
*         2:    Distribution output file
*         3:    Vertical profile output file
*
*     The output file names are generated from the input file name (e.g.
*     XXXXX.INP) and stored in the specified output file directory:
*
*     HISTORICAL OUTPUT FILE is always created for each layer, named
*     XXXXX-Hyy.DAT, where yy is the layer number. The historical output files
*     hold a number of records; one record for each output time, as
*     specified in the input file. Each record holds the following values:
*
*     A:  Time (Days, Hours, Minutes, Seconds; as specified in the input file)
*     B:  Temperature  (K)
*     C:  Total number density of frozen sulfate particles (1/ccm)
*     D:  Total number density of PSC 1 particles (1/ccm)
*     E:  Total number density of PSC 2 particles (1/ccm)
*     F:  Saturation ratio of nitric acid over NAT
*     G:  Saturation ratio of water vapor over ice
*     H:  Mixing ratio of water vapor (ppmv)
*     I:  Mixing ratio of nitric acid vapor (ppbv)
*     J:  Total (gas phase and condensed phase) mix. ratio, water (ppmv)
*     K:  Total (gas phase and condensed phase) mix. ratio, nitric acid (ppbv)
*     L:  Ratio of total to initial mixing ratio of water (%)
*     M:  Ratio of total to initial mixing ratio of nitric acid (%)
*     N:  Total volume of frozen sulfate particles (micron**3/ccm)
*     O:  Total volume of PSC 1a particles (micron**3/ccm)
*     P:  Total volume of PSC 2 particles (micron**3/ccm)
*     Q:  Total volume of all particles (micron**3/ccm)
*     R:  Total number density of PSC 1b (STS) particles (1/ccm)
*     S:  Total volume of PSC 1b (STS) particles (micron**3/ccm)
*     T:  Sulfuric acid weight fraction liquid PSC 1b aerosols [0;1]
*     U:  Mean radius liquid PSC 1b aerosols (micron)
*     V:  Median radius frozen sulfate aerosols (micron)
*     X:  Median radius PSC 1 (micron)
*     Y:  Median radius PSC 2 (micron)
*     Z:  Nitric acid weight fraction liquid PSC 1b aerosols [0;1]
*     AA: NAT condensation temperature (K)
*     BB: Ice frost point temperature (K)
*     CC: Air pressure (hPa)
*     DD: Total mixing ration of sulfuric acid (ppb)
*     EE: Optical parameter 1
*     FF: Optical parameter 2
*     GG: Optical parameter 3
*     HH: Moleratio
*     II-UU: Cumulated size distribution
*     (AEROSOL_BACKSCATTER_RATIO(J,L),J=1,N_WAVE),!      (47-53)
*     (EXTINCTION(J,L),J=1,N_WAVE),               !      (54-60)
*     (DEPOL(J,L),J=1,N_WAVE)                     !      (61-67)
*
*
*     DISTRIBUTION OUTPUT FILES are created, if a non-zero distribution-
*     dump-frequency is specified in the input file, or if a command
*     is given interactively for this (PC). The names of the distribution
*     output files are XXXXX-Dyy.DAT, where yy is the layer number.
*     The distribution output files holds a number of blocks of records.
*     A block is written to each distribution output file at each distribution
*     output time. The first record in a block holds the time (Seconds, Days or
*     Hours as specified in the input file). Then follows a number of records;
*     one record for each size bin giving the following values:
*
*     A:  Particle radius (micron)
*     B:  Differential concentration of liquid sulfate particles (1/micron ccm)
*     C:  Differential concentration of frozen sulfate particles (1/micron ccm)
*     D:  Differential concentration of PSC 1 particles (1/micron ccm)
*     E:  Differential concentration of PSC 2 particles (1/micron ccm)
*     F:  Differential concentration of all particles (1/icronm ccm)
*
*
*     A VERTICAL PROFILE OUTPUT FILE is created, if a non-zero vertical-dump-
*     frequency is specified in the input file, or if a command is given
*     interactively for this (PC). The name of the vertical profile output
*     files is XXXXX-V.DAT. The vertical profile output file holds a number of
*     blocks of records. A block is written to the vertical output
*     file at each vertical output time. The first record in a block holds
*     the time (Days or Hours; as specified in the input file). Then follows
*     a number of records; one record for each layer, starting from the
*     BOTTOM layer. Each record holds the following values:
*
*     A:  Altitude (km)
*     B:  Pressure (hPa)
*     C:  Temperature  (K)
*     D:  Total number density of frozen sulfate particles (1/ccm)
*     E:  Total number density of PSC 1 particles (1/ccm)
*     F:  Total number density of PSC 2 particles (1/ccm)
*     G:  Saturation ratio of nitric acid over NAT
*     H:  Saturation ratio of water vapor over ice
*     I:  Mixing ratio of water vapor (ppmv)
*     J:  Mixing ratio of nitric acid vapor (ppbv)
*     K:  Total (gas phase and condensed phase) mix. ratio, water (ppmv)
*     L:  Total (gas phase and condensed phase) mix. ratio, nitric acid (ppbv)
*     M:  Ratio of total to initial mixing ratio of water (%)
*     N:  Ratio of total to initial mixing ratio of nitric acid (%)
*     O:  Total surface area density of frozen sulfate particles (æm**2/ccm)
*     P:  Total surface area density of PSC 1 particles (æm**2/ccm)
*     Q:  Total surface area density of PSC 2 particles (æm**2/ccm)
*     R:  Total surface area density of all particles (æm**2/ccm)
*     S:  Total number density of liquid sulfate particles (1/ccm)
*     T:  Total surface area density of liquid sulfate particles (æm**2/ccm)
*     U:  Sulfuric acid weight fraction in LIQUID sulfate aerosols [0;1]
*     V:  Median radius liquid sulfate aerosols (æm)
*     X:  Median radius frozen sulfate aerosols (æm)
*     Y:  Median radius PSC 1 (æm)
*     Z:  Median radius PSC 2 (æm)
*     AA: Nitric acid weight fraction in LIQUID sulfate aerosols [0;1]
*     AB: NAT condensation temperature (K)
*     AC: Ice frost point temperature (K)
*     AD: Air pressure (hPa)
*     AE: Total mixing ration of sulfuric acid (ppb)
*     AF: Optical parameter 1
*     AG: Optical parameter 2
*     AH: Optical parameter 3
*     AI: Moleratio
*
*
C
      SUBROUTINE DUMP(OUTFILE,
     +        NBINS,TYPES,NCLASS,MAXLAY,NLAYER,MODE,
     +        TIME,TUNIT,
     +        TAIR,
     +        HEIGHT,PAIR,
     +        ND,
     +        MRWV,MRNA,MRWVI,MRNAI,
     +        SRATIO,
     +        AGAIN,QUIT,VDUMP,DDUMP,
     +        PTSIZE,
     +        MRTWV,MRTNA,MRTSA,WSA,WNA,
     +        NPLSA,NPFSA,NPPSC1,NPPSC2,
     +        SALSA,SAFSA,SAPSC1,SAPSC2,
     +        MRLSA,MRFSA,MRPSC1,MRPSC2,
     +        N_WAVE,AEROSOL_BACKSCATTER_RATIO,EXTINCTION,DEPOL,
     +        SIZEDIST,WORK)
C
C     Model size:
      IMPLICIT NONE
      CHARACTER*(*) OUTFILE
      INTEGER NBINS,TYPES,NCLASS,MAXLAY,NLAYER,DOT,N_WAVE
C
      REAL
     +        TIME,
     +        TAIR(MAXLAY),
     +        HEIGHT(MAXLAY),PAIR(MAXLAY),
     +        ND(NBINS,TYPES,MAXLAY),
     +        MRWV(MAXLAY),MRNA(MAXLAY),MRWVI(MAXLAY),MRNAI(MAXLAY),
     +        SRATIO(2,MAXLAY),
     +        PTSIZE(NBINS,8),
     +        MRTWV(MAXLAY),MRTNA(MAXLAY),MRTSA(MAXLAY),
     +        WSA(MAXLAY),WNA(MAXLAY),
     +        NPLSA(MAXLAY),NPFSA(MAXLAY),NPPSC1(MAXLAY),NPPSC2(MAXLAY),
     +        SALSA(MAXLAY),SAFSA(MAXLAY),SAPSC1(MAXLAY),SAPSC2(MAXLAY),
     +        MRLSA(MAXLAY),MRFSA(MAXLAY),MRPSC1(MAXLAY),MRPSC2(MAXLAY),
     +        AEROSOL_BACKSCATTER_RATIO(N_WAVE,MAXLAY),
     +        EXTINCTION(N_WAVE,MAXLAY),DEPOL(N_WAVE,MAXLAY),
     +        SIZEDIST(NCLASS,MAXLAY),WORK(NBINS)
      CHARACTER*1 TUNIT
      LOGICAL AGAIN,QUIT,VDUMP,DDUMP
C     Particle types:
      INTEGER STS,SAT,NAT,ICE
      PARAMETER(
C     Liquid supercooled ternary solution particle (SA or Type 1b PSC)
     +  STS=1,
C     Solid sulfuric acid tetrahydrate (SAT) particle; i.e. no HNO3 or excess ice
     +  SAT=2,
C     Solid nitric acid trihydrate (NAT) particle; i.e. NAT and SAT but no excess ice (Type 1a PSC)
     +  NAT=3,
C     Solid ice particle; i.e. holding excess ice, NAT, and SAT (Type 2 PSC)
     +  ICE=4)
C
      REAL 
     +     SATOT,TT,TN,TC,TCNAT,TCICE,CONV,MOLERATIO,
     +     COLOR_1,COLOR_2,COLOR_3
      INTEGER L,MODE,I,J
      CHARACTER*80 FILNAM
      LOGICAL FIRST,VERT,DIST
      SAVE FIRST,CONV
      REAL MHNO3,MH2O,MAIR,RGAS,C4,RHOAIR
      PARAMETER(
     +          MH2O=18.0153E-3,
     +          MHNO3=63.01E-3,
     +          MAIR=28.9644E-3,
     +          RGAS=8.31441,
     +          C4=MAIR/RGAS)
C
      REAL TEMP,THETAP,PPP,P0
      PARAMETER(P0=1000.0E2)
      THETAP(PPP,TEMP)=TEMP*(P0/PPP)**0.286
C----------------------------------------------------------------------------
      DATA FIRST/.TRUE./
C----------------------------------------------------------------------------
      IF(FIRST) THEN
          FIRST=.FALSE.
          DOT=INDEX(OUTFILE(2:LEN(OUTFILE)),'.')
          FILNAM=OUTFILE(1:DOT)//'-HXX.DAT'
          DO L=1,NLAYER
              WRITE(FILNAM(DOT+3:DOT+4),'(I2.2)') L
              OPEN(10+L,FILE=FILNAM)
          END DO
          VERT=.FALSE.
          DIST=.FALSE.
          IF(MODE.EQ.2) THEN
              CONV=1.0E6
          ELSE IF(MODE.EQ.3) THEN
              CONV=1.0E12
          ENDIF
cc          OPEN(91,FILE='E354.DAT')
cc          OPEN(92,FILE='E440.DAT')
cc          OPEN(93,FILE='E603.DAT')
cc          OPEN(94,FILE='E779.DAT')
cc          OPEN(95,FILE='E922.DAT')
cc          OPEN(96,FILE='E1020.DAT')
      ENDIF
C----------------------------------------------------------------------------
      IF(AGAIN) FIRST=.TRUE.
C----------------------------------------------------------------------------
      TT=TIME
      IF(INDEX('Hh',TUNIT).GT.0) THEN
          TT=TIME/3600.0
      ELSE IF(INDEX('Dd',TUNIT).GT.0) THEN
          TT=TIME/(24.0*3600.0)
      ELSE IF(INDEX('Mm',TUNIT).GT.0) THEN
          TT=TIME/(60.0)
      ENDIF
C----------------------------------------------------------------------------
C     Write historical data to file:
      DO L=1,NLAYER
          SATOT=SALSA(L)+SAFSA(L)+SAPSC1(L)+SAPSC2(L)
          TN=TCNAT(MRNA(L)*PAIR(L),MRWV(L)*PAIR(L))
          TC=TCICE(MRWV(L)*PAIR(L))
          MOLERATIO=1.0E29
CMP   Definition of 3 optical color indices to output:
CC        COLOR_1=AEROSOL_BACKSCATTER_RATIO(5,L)/
CC   +            AEROSOL_BACKSCATTER_RATIO(2,L)
CC        COLOR_2=AEROSOL_BACKSCATTER_RATIO(3,L)/
CC   +            AEROSOL_BACKSCATTER_RATIO(7,L)
CC        COLOR_3=EXTINCTION(6,L)/EXTINCTION(1,L)
          COLOR_1=AEROSOL_BACKSCATTER_RATIO(5,L)
          COLOR_2=AEROSOL_BACKSCATTER_RATIO(2,L)
          COLOR_3=AEROSOL_BACKSCATTER_RATIO(5,L)/
     +            AEROSOL_BACKSCATTER_RATIO(2,L)
cc          COLOR_1=EXTINCTION(6,L)*1.0E3
cc          COLOR_2=EXTINCTION(1,L)*1.0E3
cc          COLOR_3=EXTINCTION(6,L)/EXTINCTION(1,L)
          IF(WNA(L).GT.1.0E-7)
     +         MOLERATIO=(1.0-WSA(L)-WNA(L))*MHNO3/(WNA(L)*MH2O)
          WRITE(10+L,'(70(1PE12.5),0P)')
     +       TT,TAIR(L),                                 !A,B  (1,2)
     +       AMAX1(1.0E-30,NPFSA(L)*1.0E-6),             !C    (3)
     +       AMAX1(1.0E-30,NPPSC1(L)*1.0E-6),            !D    (4)
     +       AMAX1(1.0E-30,NPPSC2(L)*1.0E-6),            !E    (5)
     +       SRATIO(1,L),SRATIO(2,L),                    !F,G  (6,7)
     +       MRWV(L)*1.0E6,MRNA(L)*1.0E9,                !H,I  (8,9)
     +       MRTWV(L)*1.0E6,MRTNA(L)*1.0E9,              !J,K  (10,11)
     +       MRTWV(L)/MRWVI(L),MRTNA(L)/MRNAI(L),           !L,M  (12,13)
     +       AMAX1(1.0E-30,CONV*SAFSA(L)),               !N    (14)
     +       AMAX1(1.0E-30,CONV*SAPSC1(L)),              !O    (15)
     +       AMAX1(1.0E-30,CONV*SAPSC2(L)),              !P    (16)
     +       AMAX1(1.0E-30,CONV*SATOT),                  !Q    (17)
     +       AMAX1(1.0E-30,NPLSA(L)*1.0E-6),             !R    (18)
     +       AMAX1(1.0E-30,CONV*SALSA(L)),               !S    (19)
     +       WSA(L),MRLSA(L)*1.0E6,MRFSA(L)*1.0E6,       !T,U,V(20,21,22)
     +       MRPSC1(L)*1.0E6,MRPSC2(L)*1.0E6,            !W,X  (23,24)
     +       WNA(L),TN,TC,PAIR(L)/100.0,MRTSA(L)*1.0E9,  !Y,Z,AA,BB,CC (25,26,27,28,29)
     +       COLOR_1,COLOR_2,COLOR_3,MOLERATIO,          !DD,EE,FF,GG  (30,31,32,33)
     +       (SIZEDIST(J,L)*1.0E-6,J=1,NCLASS),          !HH-TT (34-46)
     +       (AEROSOL_BACKSCATTER_RATIO(J,L),J=1,N_WAVE),!      (47-53)
     +       (EXTINCTION(J,L)*1.0E3,J=1,N_WAVE),         !      (54-60)
     +       (DEPOL(J,L),J=1,N_WAVE),                    !      (61-67)
     +        THETAP(PAIR(L),TAIR(L)),HEIGHT(L)/1000.0          !      (68-69)
          IF(QUIT) CLOSE(10+L,STATUS='KEEP')
          IF(AGAIN) CLOSE(10+L,STATUS='DELETE')
      ENDDO
cc      write(91,'(46(1PE12.5,1x),0P)') tt,EXTINCTION(1,1)*1.0E3
cc      write(92,'(46(1PE12.5,1x),0P)') tt,EXTINCTION(2,1)*1.0E3
cc      write(93,'(46(1PE12.5,1x),0P)') tt,EXTINCTION(3,1)*1.0E3
cc      write(94,'(46(1PE12.5,1x),0P)') tt,EXTINCTION(4,1)*1.0E3
cc      write(95,'(46(1PE12.5,1x),0P)') tt,EXTINCTION(5,1)*1.0E3
cc      write(96,'(46(1PE12.5,1x),0P)') tt,EXTINCTION(6,1)*1.0E3
C----------------------------------------------------------------------------
      IF(DDUMP) THEN
C         Dump size distributions on file:
          DDUMP=.FALSE.
C         Generic name of distribution output file:
          DO L=1,NLAYER
              IF(.NOT.DIST) THEN
                  DOT=INDEX(OUTFILE(2:LEN(OUTFILE)),'.')
                  FILNAM=OUTFILE(1:DOT)//'-DXX.DAT'
                  WRITE(FILNAM(DOT+3:DOT+4),'(I2.2)') L
                  OPEN(40+L,FILE=FILNAM)
              ENDIF
              I=NBINS
              WORK(I)=(ND(I,STS,L)+ND(I,SAT,L)+ND(I,NAT,L)+
     +               ND(I,ICE,L))
              DO I=NBINS-1,1,-1
                  WORK(I)=WORK(I+1)+
     +                  (ND(I,STS,L)+ND(I,SAT,L)+ND(I,NAT,L)+
     +                   ND(I,ICE,L))
              END DO
              RHOAIR=C4*PAIR(L)/TAIR(L)
              WRITE(40+L,'(1PE12.5)') TT
              DO I=1,NBINS
                   WRITE(40+L,'(10(1PE12.5),0P)')
     +PTSIZE(I,1)*1.0E6,
CC   +AMAX1(1.0E-30,ND(I,STS,L)*RHOAIR*1.0E-6/PTSIZE(I,6)),
CC   +AMAX1(1.0E-30,ND(I,SAT,L)*RHOAIR*1.0E-6/PTSIZE(I,6)),
CC   +AMAX1(1.0E-30,ND(I,NAT,L)*RHOAIR*1.0E-6/PTSIZE(I,6)),
CC   +AMAX1(1.0E-30,ND(I,ICE,L)*RHOAIR*1.0E-6/PTSIZE(I,6)),
CC   +AMAX1(1.0E-30,(ND(I,STS,L)+ND(I,SAT,L)+ND(I,NAT,L)+
CC   +               ND(I,ICE,L))*RHOAIR*1.0E-6/PTSIZE(I,6))
     +AMAX1(1.0E-30,PTSIZE(I,1)*ND(I,STS,L)*RHOAIR/PTSIZE(I,6)),
     +AMAX1(1.0E-30,PTSIZE(I,1)*ND(I,SAT,L)*RHOAIR/PTSIZE(I,6)),
     +AMAX1(1.0E-30,PTSIZE(I,1)*ND(I,NAT,L)*RHOAIR/PTSIZE(I,6)),
     +AMAX1(1.0E-30,PTSIZE(I,1)*ND(I,ICE,L)*RHOAIR/PTSIZE(I,6)),
     +AMAX1(1.0E-30,PTSIZE(I,1)*(ND(I,STS,L)+ND(I,SAT,L)+ND(I,NAT,L)+
     +               ND(I,ICE,L))*RHOAIR/PTSIZE(I,6)),
     +AMAX1(1.0E-30,WORK(I)*RHOAIR*1.0E-6)
              ENDDO
          ENDDO
          DIST=.TRUE.
      ENDIF
      IF(QUIT.AND.DIST) THEN
          DO L=1,NLAYER
             CLOSE(40+L,STATUS='KEEP')
          ENDDO
      ENDIF
      IF(AGAIN.AND.DIST) THEN
          DO L=1,NLAYER
             CLOSE(40+L,STATUS='DELETE')
          ENDDO
      ENDIF
C----------------------------------------------------------------------------
      IF(VDUMP) THEN
C         Dump vertical profile:
          VDUMP=.FALSE.
C         Name of vertical profile output file:
          IF(.NOT.VERT) THEN
              VERT=.TRUE.
              DOT=INDEX(OUTFILE(2:LEN(OUTFILE)),'.')
              FILNAM=OUTFILE(1:DOT)//'-V.DAT'
              OPEN(10,FILE=FILNAM)
          ENDIF
CC          WRITE(10,'(1PE12.5)') TT
          DO L=1,NLAYER
              SATOT=SALSA(L)+SAFSA(L)+SAPSC1(L)+SAPSC2(L)
              TN=TCNAT(MRNA(L)*PAIR(L),MRWV(L)*PAIR(L))
              TC=TCICE(MRWV(L)*PAIR(L))
CMP   Definition of 3 optical color indices to output:
CC            COLOR_1=AEROSOL_BACKSCATTER_RATIO(5,L)/
CC   +                AEROSOL_BACKSCATTER_RATIO(2,L)
CC            COLOR_2=AEROSOL_BACKSCATTER_RATIO(3,L)/
CC   +                AEROSOL_BACKSCATTER_RATIO(7,L)
CC            COLOR_3=EXTINCTION(6,L)/EXTINCTION(1,L)
          COLOR_1=AEROSOL_BACKSCATTER_RATIO(3,L)
          COLOR_2=AEROSOL_BACKSCATTER_RATIO(2,L)
          COLOR_3=DEPOL(3,L)
              IF(WNA(L).GT.1.0E-7)
     +             MOLERATIO=(1.0-WSA(L)-WNA(L))*MHNO3/(WNA(L)*MH2O)
              WRITE(10,'(70(1PE12.5),0P)')
     +                  TT,HEIGHT(L)/1000.0,PAIR(L)/100.0,TAIR(L),
     +                  AMAX1(1.0E-30,NPFSA(L)*1.0E-6),
     +                  AMAX1(1.0E-30,NPPSC1(L)*1.0E-6),
     +                  AMAX1(1.0E-30,NPPSC2(L)*1.0E-6),
     +                  SRATIO(1,L),SRATIO(2,L),
     +                  MRWV(L)*1.0E6,MRNA(L)*1.0E9,
     +                  MRTWV(L)*1.0E6,MRTNA(L)*1.0E9,
     +                  MRTWV(L)/MRWVI(L),MRTNA(L)/MRNAI(L),
     +                  AMAX1(1.0E-30,CONV*SAFSA(L)),
     +                  AMAX1(1.0E-30,CONV*SAPSC1(L)),
     +                  AMAX1(1.0E-30,CONV*SAPSC2(L)),
     +                  AMAX1(1.0E-30,CONV*SATOT),
     +                  AMAX1(1.0E-30,NPLSA(L)*1.0E-6),
     +                  AMAX1(1.0E-30,CONV*SALSA(L)),
     +                  WSA(L),MRLSA(L)*1.0E6,MRFSA(L)*1.0E6,
     +                  MRPSC1(L)*1.0E6,MRPSC2(L)*1.0E6,
     +                  WNA(L),TN,TC,PAIR(L)/100.0,MRTSA(L)*1.0E9,
     +                  COLOR_1,COLOR_2,COLOR_3,MOLERATIO,
     +                  (SIZEDIST(J,L)*1.0E-6,J=1,NCLASS),
     +                  (AEROSOL_BACKSCATTER_RATIO(J,L),J=1,N_WAVE),
     +                  (EXTINCTION(J,L)*1.0E3,J=1,N_WAVE),
     +                  (DEPOL(J,L),J=1,N_WAVE),
     +                  THETAP(PAIR(L),TAIR(L)),HEIGHT(L)/1000.0
          ENDDO
      ENDIF
      IF(QUIT.AND.VERT) CLOSE(10,STATUS='KEEP')
      IF(AGAIN.AND.VERT) CLOSE(10,STATUS='DELETE')
C----------------------------------------------------------------------------
      RETURN
      END
C
C
C*****************************************************************************
      REAL FUNCTION MIXRNA(HEIGHT,LATITU)
C*****************************************************************************
C     This function calculated the initial mixing ratio of HNO3
C     as function of height using a LIMS-profile.
C     Gille & Russel, JGR 89,5125,1984.
C
C     Input/output:
C
C         HEIGHT:  (m)            Height
C         LATITU:  (latitude)     Northern H (positive); Southern H (neg.)
C         MRNA:                   Mixing ratio of HNO3
C
      IMPLICIT NONE
      REAL    A0,A1,A2,A3,A4,A5,A6,
     +        B0,B1,B2,B3,B4,B5,
     +        H,HH,MIX,HEIGHT,LATITU
      DATA A0,A1,A2,A3,A4,A5,A6/
     #-1.06884E+03, 2.91968E+02,-3.26783E+01, 1.91549E+00,-6.16650E-02,
     # 1.03130E-03,-7.00183E-06/                                                 
      DATA B0,B1,B2,B3,B4,B5/
     # 3.14998E+02,-7.09302E+01, 6.13750E+00,-2.52636E-01, 4.98268E-03,
     #-3.79818E-05/                                                              
C
      HH=HEIGHT/1000.0
      H=AMAX1(15.1,AMIN1(HH,30.0))
      IF(LATITU.GE.0) THEN
          MIX=1.0E-9*((((((A6*H+A5)*H+A4)*H+A3)*H+A2)*H+A1)*H+A0)
      ELSE
          MIX=1.0E-9*(((((B5*H+B4)*H+B3)*H+B2)*H+B1)*H+B0)
      ENDIF
      IF(HH.LT.15.1) THEN
          MIX=AMAX1(MIX-(15.1-HH)*(MIX-1.0E-9)/2.5,1.0E-11)
      ELSE IF(HH.GT.30.0) THEN
          MIX=AMAX1(MIX-(HH-30.0)*(MIX-1.0E-9)/4.0,1.0E-11)
      ENDIF
      MIXRNA=MIX
      RETURN
      END
C*****************************************************************************
      REAL FUNCTION LIMSNA(PAIR,LATITU)
C*****************************************************************************
C     This function calculated the initial mixing ratio of HNO3
C     as function of pressure altitude height using a LIMS-profile.
C     Gille & Russel, JGR 89,5125,1984.
C
C     Input/output:
C
C         PAIR:    (Pa)           Pressure altitude
C         LATITU:  (latitude)     Northern H (positive); Southern H (neg.)
C         MRNA:                   Mixing ratio of HNO3
C
      IMPLICIT NONE
      REAL HSCALE
C     LIMS data fitted with HSCALE=6.5 km. 
C     This parameter should not be changed:
      PARAMETER(HSCALE=6.5)
      REAL    A0,A1,A2,A3,A4,A5,A6,
     +        B0,B1,B2,B3,B4,B5,
     +        H,HH,MIX,LATITU,PAIR
      DATA A0,A1,A2,A3,A4,A5,A6/
     #-1.06884E+03, 2.91968E+02,-3.26783E+01, 1.91549E+00,-6.16650E-02,          
     # 1.03130E-03,-7.00183E-06/                                                 
      DATA B0,B1,B2,B3,B4,B5/
     # 3.14998E+02,-7.09302E+01, 6.13750E+00,-2.52636E-01, 4.98268E-03,          
     #-3.79818E-05/                                                              
C
      HH=HSCALE*ALOG(1013.25E2/PAIR)
      H=AMAX1(15.1,AMIN1(HH,30.0))
      IF(LATITU.GE.0) THEN
          MIX=1.0E-9*((((((A6*H+A5)*H+A4)*H+A3)*H+A2)*H+A1)*H+A0)
      ELSE
          MIX=1.0E-9*(((((B5*H+B4)*H+B3)*H+B2)*H+B1)*H+B0)
      ENDIF
      IF(HH.LT.15.1) THEN
          MIX=AMAX1(MIX-(15.1-HH)*(MIX-1.0E-9)/2.5,1.0E-11)
      ELSE IF(HH.GT.30.0) THEN
          MIX=AMAX1(MIX-(HH-30.0)*(MIX-1.0E-9)/4.0,1.0E-11)
      ENDIF
      LIMSNA=MIX
      RETURN
      END

C*****************************************************************************
C     RAMP FUNCTION
C*****************************************************************************
C
C
C
      FUNCTION RAMP(T,T1,Y1,T2,Y2,T3,Y3,T4,Y4)
C     FUNCTION RAMP(T,T1,Y1,T2,Y2,T3,Y3,T4,Y4)
C     ----------------------------------------
C     RAMP FUNCTION
C
C     ^ Y
C     !
C     !
C     !
C     !
C     !              Y2           Y3
C     !               ooooooooooooo
C     !              o             o
C     !             o                o
C     !            o                  o
C     !           o                    oooooooooooooooo
C     !  ooooooooo Y1                  Y4
C     !
C     !____________________________________________________________>
C                T1   T2          T3   T4
C                                                                T
C
C     RESTRICTIONS: T4 >= T3 >= T2 >= T1
C
C     =================================================================
      IMPLICIT NONE
      REAL RAMP,T,T1,T2,T3,T4,Y1,Y2,Y3,Y4,SUBRAM
C
      IF((T4.LT.T3).OR.(T3.LT.T2).OR.(T2.LT.T1)) THEN
         WRITE(*,*) 'INPUT ERROR TO FUNCTION RAMP'
      ELSE IF(T2-T1.GT.0.0) THEN
         IF(T.LE.T1) THEN
            RAMP=Y1
         ELSE IF(T.GT.T1.AND.T.LT.T2) THEN
            RAMP=Y1+(Y2-Y1)*(T-T1)/(T2-T1)
         ELSE
            RAMP=SUBRAM(T,T2,Y2,T3,Y3,T4,Y4)
         ENDIF
      ELSE
         IF(T.LT.T1) THEN
            RAMP=Y1
         ELSE
            RAMP=SUBRAM(T,T2,Y2,T3,Y3,T4,Y4)
         ENDIF
      ENDIF
      RETURN
      END
C
C
      FUNCTION SUBRAM(T,T2,Y2,T3,Y3,T4,Y4)
      IMPLICIT NONE
      REAL SUBRAM,T,T2,Y2,T3,Y3,T4,Y4,SUBSUB
      IF(T3-T2.GT.0.0) THEN
         IF(T.LT.T3) THEN
            SUBRAM=Y2+(Y3-Y2)*(T-T2)/(T3-T2)
         ELSE
            SUBRAM=SUBSUB(T,T3,Y3,T4,Y4)
         ENDIF
      ELSE
         IF(T.EQ.T2) THEN
            SUBRAM=AMAX1(Y2,Y3)
         ELSE
            SUBRAM=SUBSUB(T,T3,Y3,T4,Y4)
         ENDIF
      ENDIF
      RETURN
      END
C
C
      FUNCTION SUBSUB(T,T3,Y3,T4,Y4)
      IMPLICIT NONE
      REAL SUBSUB,T,T3,Y3,T4,Y4
      IF(T4-T3.GT.0.0) THEN
         IF(T.LE.T4) THEN
            SUBSUB=Y3+(Y4-Y3)*(T-T3)/(T4-T3)
         ELSE
            SUBSUB=Y4
         ENDIF
      ELSE
         IF(T.EQ.T3) THEN
            SUBSUB=AMAX1(Y3,Y4)
         ELSE
            SUBSUB=Y4
         ENDIF
      ENDIF
      RETURN
      END
C
C******************************************************************************
      REAL FUNCTION TEMP2D(TIME,TSTART,LAYER,LATITU,TMPFIL,TUNIT)
C******************************************************************************
C     Ad hoc function used to generate temperatures from 2-D temperature
C     table (i.e. temperature as function of time and altitude layer)
C     Input: TIME :  model time (s)
C            TSTART: model start time (s)
C                    (dummy if temperature oscilation not needed)
C            LAYER:  model layer number
C            LATITU: latitude 
C                    (dummy if temeprature oscilation not needed)
C     Output:        temperature in layer no. LAYER at time TIME
C
C
C
C     The 2-D temperature calculations are based on an input table file, 
C     given in the input file, holding the following records:
C     number of time settings in the table (NTIMES)
C     number of altitude layers (NLAYER)
C     for each time setting the following records:
C             time (units corresponding to input file pscinput.inp)
C     followed by NLAYER records holding (top to bottom)
C             temperature (K)
C
C     Temperatures for each layer are calculated by spline interpolation
C     at times between table time settings. Latitude dependent
C     temperature oscillations are overlayed if logical OSC is set
C     to .TRUE. on compilation.
C
      IMPLICIT NONE
      INTEGER MAXTIM,MAXLAY,NLAYER,NTIMES,I,L,LAYER
      PARAMETER(MAXLAY=15,MAXTIM=150)
C
      REAL TIMEA(MAXTIM),TEMPA(MAXTIM,MAXLAY),
     +        YWORK(MAXTIM),Y2TEMP(MAXTIM,MAXLAY),TIME,TAIR,
     +        LATITU,PERIOD,PI,AMP,TSTART,DAY,TT
      CHARACTER *(*) TMPFIL,TUNIT
C
      LOGICAL FIRST,OSC
      SAVE FIRST,TIMEA,TEMPA,OSC
C
      PARAMETER(PI=3.1415926536,
C     Oscilation period (s):
     +PERIOD=15.0*24.0*3600.0)
      DATA FIRST/.TRUE./
      DATA OSC/.FALSE./
      DATA AMP/2.0/
C
      IF(FIRST) THEN
          FIRST=.FALSE.
          OPEN(1,FILE=TMPFIL)
          READ(1,*) NTIMES
          READ(1,*) NLAYER
          IF(NLAYER.GT.MAXLAY) THEN
              WRITE(*,*) 'MAXLAY parameter in subroutine TEMP2D wrong'
              STOP
          ENDIF
          DO I=1,NTIMES
              READ(1,*) TIMEA(I)
              DO L=1,NLAYER
                   READ(1,*) TEMPA(I,L)
              ENDDO
          ENDDO
          CLOSE(1)
C
          IF(INDEX('Hh',TUNIT).GT.0) THEN
              DO I=1,NTIMES
                  TIMEA(I)=TIMEA(I)*3600.0
              END DO
          ELSE IF(INDEX('Dd',TUNIT).GT.0) THEN
              DO I=1,NTIMES
                  TIMEA(I)=TIMEA(I)*3600.0*24.0
              END DO
          ELSE IF(INDEX('Mm',TUNIT).GT.0) THEN
              DO I=1,NTIMES
                  TIMEA(I)=TIMEA(I)*60.0
              END DO
          ENDIF
C
          DO L=1,NLAYER
              CALL SPLINE(TIMEA,TEMPA(1,L),NTIMES,YWORK,Y2TEMP(1,L))
          ENDDO
      ENDIF
      TT=AMAX1(TIMEA(1),AMIN1(TIME,TIMEA(NTIMES)))
      CALL SPLINT(TIMEA,TEMPA(1,LAYER),Y2TEMP(1,LAYER),
     +            NTIMES,TT,TAIR)
      IF(OSC) THEN
          DAY=TIME/(24.0*3600.0)
          IF(DAY.LT.0.0) DAY=DAY+365.0
          CALL WAVEAM(INT(DAY),LATITU,AMP)
          TAIR=TAIR+AMP*SIN((TIME-TSTART)*2.0*PI/PERIOD)
      ENDIF
      TEMP2D=TAIR
      RETURN
      END
C
      SUBROUTINE WAVEAM(M,P,A)
C     SUBROUTINE USED TO CALCULATE TEMPERATURE WAVE AMPLITUDE
C     FOR NACR 2-D MODEL TEMPERATURES.
C
C     INPUT:
C     M:  DAY-NUMBER (INTEGER)
C     P:  LATITUDE (POSITIV ON NH, NEGATIVE ON SH, REAL)
C     OUTPUT: TEMPERATURE AMPLITUDE (K)
C
      IMPLICIT NONE
      INTEGER M, N
      REAL    A
      REAL    PI, DTR, TMIN, AN, AS
      REAL    D, P, YY, Y, T, TT
C
      PI=4.*ATAN(1.)
      DTR=PI/180.
      TMIN=0.05
      AN=10.
      AS=6.
C
      IF (M.LT.172) N=365+M-172
      IF (M.GE.172) N=M-172
      D=FLOAT(N)
C          P=FLOAT((L-18)*5)
          YY=ABS(P)
          IF (YY.LT.30.) THEN
              Y=0.
          ELSE
              Y=(YY-60.)*90.*DTR/30.
              Y=COS(Y)*COS(Y)
          ENDIF
          T=TMIN
          TT=0.
          IF (P.LT.0.) THEN
              IF (N.LE.10) THEN
                  TT=(D+365.-314.)*90.*DTR/61.
              ELSE 
                  IF (N.LE.132) THEN
                      TT=(D-71.)*90.*DTR/61.
                  ELSE
                      IF (N.GE.253) THEN
                          TT=(D-314.)*90.*DTR/61.
                      ENDIF
                  ENDIF
              ENDIF
              TT=COS(TT)*COS(TT)
              IF (N.LE.10.OR.N.GT.314) THEN
                  T=.25*(1.+TT)
              ELSE
                  IF (N.LE.71.AND.N.GT.10) THEN
                      T=0.25+0.75*TT
                  ELSE
                      IF (N.LE.132.AND.N.GT.71) THEN
                          T=TMIN+(1.-TMIN)*TT
                      ELSE
                          IF (N.GE.253.AND.N.LE.314) THEN
                              T=TMIN+(.5-TMIN)*TT
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
          ELSE
              IF (N.GE.86.AND.N.LE.298) THEN
                  TT=(D-192.)*90.*DTR/106.
                  T=TMIN+(1.-TMIN)*COS(TT)*COS(TT)
              ENDIF
          ENDIF
          A=AN*Y*T
          IF (P.LT.0.) A=AS*Y*T
C          WRITE(*,*) 'DAY= ',M,' LAT= ',L,' AMP= ',A
      RETURN
      END
C
C**************************************************************************
      SUBROUTINE COLUMN(
C**************************************************************************
     +        TIME,TUNIT,
     +        NBINS,TYPES,NLAYER,
     +        TAIR,
     +        PAIR,
     +        MRNA,MRWV,
     +        ND,
     +        MCS,MCN,MCW,
     +        SDND,SDMCS,SDMCN,SDMCW,
     +        THICKNESS,
     +        CND,CMCS,CMCN,CMCW)
      INTEGER NBINS,TYPES,NLAYER
      REAL    TIME,
     +        TAIR(NLAYER),
     +        PAIR(NLAYER),
     +        MRNA(NLAYER),MRWV(NLAYER),
     +        ND(NBINS,TYPES,NLAYER),
     +        MCS(NBINS,TYPES,NLAYER),
     +        MCN(NBINS,TYPES,NLAYER),
     +        MCW(NBINS,TYPES,NLAYER),
     +        SDND(NBINS,TYPES),
     +        SDMCS(NBINS,TYPES),SDMCN(NBINS,TYPES),SDMCW(NBINS,TYPES),
     +        THICKNESS(NLAYER),
     +        CND,
     +        CMCS,
     +        CMCN,
     +        CMCW
      CHARACTER*1 TUNIT
      INTEGER I,J,L
      REAL MH2O,MHNO3,MAIR,RGAS,C3,RHOAIR,X1,TT
      PARAMETER(
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
C       Molar weight of nitric acid (kg/mole)
     +          MHNO3=63.01E-3,
C       Molar weight of air (kg/mole)
     +          MAIR=28.9644E-3,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441,
     +          C3=MAIR/RGAS)
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      SAVE FIRST
      IF(FIRST) THEN
          FIRST=.FALSE.
          OPEN(9,FILE='COLUMN.DAT')
      ENDIF
C----------------------------------------------------------------------------
      TT=TIME
      IF(INDEX('Hh',TUNIT).GT.0) THEN
          TT=TIME/3600.0
      ELSE IF(INDEX('Dd',TUNIT).GT.0) THEN
          TT=TIME/(24.0*3600.0)
      ELSE IF(INDEX('Mm',TUNIT).GT.0) THEN
          TT=TIME/(60.0)
      ENDIF
      CND=0.0
      CMCS=0.0
      CMCN=0.0
      CMCW=0.0
      DO L=1,NLAYER
          RHOAIR=C3*PAIR(L)/TAIR(L)
          DO J=1,TYPES
              DO I=1,NBINS
                  X1=ND(I,J,L)*THICKNESS(L)*RHOAIR
                  CND=CND+X1
                  CMCS=CMCS+MCS(I,J,L)*X1
                  CMCN=CMCN+MCN(I,J,L)*X1
                  CMCW=CMCW+MCW(I,J,L)*X1
              END DO
          END DO
          CMCN=CMCN+
     +         ((MRNA(L)*PAIR(L)*MHNO3*THICKNESS(L))/(RGAS*TAIR(L)))
          CMCW=CMCW+
     +         ((MRWV(L)*PAIR(L)*MH2O*THICKNESS(L))/(RGAS*TAIR(L)))
      END DO

      RHOAIR=C3*PAIR(NLAYER)/TAIR(NLAYER)
      DO J=1,TYPES
          DO I=1,NBINS
              X1=SDND(I,J)*THICKNESS(NLAYER)
              CND=CND+X1
              CMCS=CMCS+SDMCS(I,J)*X1
              CMCN=CMCN+SDMCN(I,J)*X1
              CMCW=CMCW+SDMCW(I,J)*X1
          ENDDO
      END DO
      WRITE(9,'(10(1PE12.5,1X),0P)') TT,CND,CMCS,CMCN,CMCW
      RETURN
      END
C**************************************************************************
CMP  The two following dummy subroutines
C                  SUBROUTINE INPUTFILE and SUBROUTINE IA
C     must be linked to PSCMODEL.FOR for non-interactive batch mode runs
C     using standard fortran90, i.e. the statements in the subroutines must be
C     activated by removing any leading 'CC' if the statements are commented out.
C**************************************************************************
C     SUBROUTINE INPUTFILE(INPDIR,INPFILE,MORE)
C**************************************************************************
        SUBROUTINE INPUTFILE(INPDIR,INPFILE,INPFILE2,MORE)
C
C     Subroutine used to to provide input file name(s) to the PSCMODEL
C     running in batch mode.
C     Depending on the value of the logical flag SCHEDULE the input file
C     name can be given manually as input from the standard input unit
C     (SCHEDULE=.FALSE.), or as sequence of input files (i.e a sequence
C     of simulations can be performed (SCHEDULE=.TRUE.). In this case
C     an ASCII file with name pscmodel.run must reside in the input
C     directory holding the file names of the input files to be processed.
C     Standard fortran 90.
      IMPLICIT NONE
      CHARACTER *(*) INPDIR,INPFILE, INPFILE2
	CHARACTER(len=32) ARGUMENT
      LOGICAL MORE,SCHEDULE,FIRST
      DATA FIRST/.TRUE./
      INTEGER MAXRUNS,NRUN,IRUN
      DATA NRUN,IRUN/0,0/
      PARAMETER(MAXRUNS=100)
      CHARACTER*80 SCH_FILE(MAXRUNS)
      CHARACTER*80 LINE
      SAVE NRUN,IRUN,SCH_FILE,FIRST
      PARAMETER (SCHEDULE=.FALSE.)
      IF(SCHEDULE) THEN
          IF(FIRST) THEN
              WRITE(*,*) 'PSCMODEL schedule. Input files:'
              OPEN(1,FILE=INPDIR//'pscmodel.run')
  1           READ(1,'(A)',END=100) LINE
              NRUN=NRUN+1
              SCH_FILE(NRUN)=LINE
              WRITE(*,*) SCH_FILE(NRUN)
              IF(NRUN.LT.MAXRUNS) GOTO 1
  100         CONTINUE
              CLOSE(1)
              FIRST=.FALSE.
          ENDIF
      ELSE
		!!
		!! Caleb Wherry
	    !! Input *.inp file as a command line argument.
		!! Used when running model with GUI application.
		!! Comment out the following "WRITE" and "READ" lines when using this with a GUI.
		!!
		!
		CALL GETARG(1,ARGUMENT)
		READ(ARGUMENT,*) INPFILE

		CALL GETARG(2,ARGUMENT)
		READ(ARGUMENT,*) INPFILE2
C		!!
C		!! Used to explicitly input the *.inp file you want to use.
C		INPFILE = 'ramp_1_f_s.inp'
		!!
		!!
		!WRITE(6, '("Enter input file name(.inp): ")', ADVANCE='NO')
		!READ(5,'(A)') INPFILE	
		!
			 
          MORE=.FALSE.
          RETURN
      ENDIF
      IRUN=IRUN+1
      INPFILE=SCH_FILE(IRUN)
      MORE=IRUN.LT.NRUN
      RETURN
      END

*************************************************************************
      SUBROUTINE IA(
*************************************************************************
     +        NBINS,TYPES,MAXLAY,NLAYER,MODE,
     +        TIME,TSTOP,TUNIT,
     +        TAIR,PAIR,
     +        ND,
     +        MRWV,MRNA,
     +        SRATIO,
     +        PRESENCE,
     +        ONESTP,AGAIN,QUIT,VDUMP,DDUMP,
     +        PTSIZE,
     +        MRTWV,MRTNA,
     +        NPLSA,NPFSA,NPPSC1,NPPSC2,
     +        SALSA,SAFSA,SAPSC1,SAPSC2,
     +        MDLSA,MDFSA,
     +        MDPSC1,MDPSC2,
     +        WSA,WNA,
     +        X,Y,HEADER,INPFILE)
C     Dummy subroutine giving the simulation time and % of simulation
C     duration performed during runs in batch mode. Standard fortran 90.
      IMPLICIT NONE
      INTEGER NBINS,TYPES,MAXLAY,NLAYER,DOT,MODE

      REAL
     +        TIME,TSTOP,
     +        TAIR(MAXLAY),PAIR(MAXLAY),
     +        ND(NBINS,TYPES,MAXLAY),
     +        MRWV(MAXLAY),MRNA(MAXLAY),
     +        SRATIO(2,MAXLAY),
     +        PTSIZE(NBINS,8),
     +        MRTWV(MAXLAY),MRTNA(MAXLAY),
     +        NPLSA(MAXLAY),NPFSA(MAXLAY),
     +        NPPSC1(MAXLAY),NPPSC2(MAXLAY),
     +        SALSA(MAXLAY),SAFSA(MAXLAY),
     +        SAPSC1(MAXLAY),SAPSC2(MAXLAY),
     +        MDLSA(MAXLAY),MDFSA(MAXLAY),
     +        MDPSC1(MAXLAY),MDPSC2(MAXLAY),
     +        WSA(MAXLAY),WNA(MAXLAY),
     +        X(NBINS),Y(NBINS)
      REAL TT,TS,TSTART,DT
      LOGICAL PRESENCE(TYPES,MAXLAY)
      LOGICAL ONESTP,AGAIN,QUIT,VDUMP,DDUMP
      LOGICAL FIRST
      CHARACTER*1 TUNIT
      CHARACTER*(*) HEADER,INPFILE
      DATA FIRST/.TRUE./
      SAVE FIRST,TSTART,DT,DOT
      IF(INDEX('Hh',TUNIT).GT.0) THEN
          TT=TIME/3600.0
          TS=TSTOP/3600.0
      ELSE IF(INDEX('Dd',TUNIT).GT.0) THEN
          TT=TIME/(3600.0*24.0)
          TS=TSTOP/(3600.0*24.0)
      ELSE IF(INDEX('Mm',TUNIT).GT.0) THEN
          TT=TIME/60.0
          TS=TSTOP/60.0
      ELSE
          TT=TIME
          TS=TSTOP
      ENDIF
      IF(FIRST) THEN
          FIRST=.FALSE.
          WRITE(*,*) INPFILE
          WRITE(*,*) HEADER
          TSTART=TT
          DT=TS-TSTART
          DOT=INDEX(INPFILE(2:),'.')+1
      ENDIF
      WRITE(*,'(1X,A,A1,A,F8.2,F8.2,A,A)')
     +     'Time (',TUNIT,'): ',TT,
     +      100.0*(TT-TSTART)/DT, ' % to end. ',INPFILE(1:DOT+3)
      IF(AGAIN) FIRST=.TRUE.
      RETURN
      END
