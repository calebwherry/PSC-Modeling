******************************************************************************
*     SUBROUTINE PSCBOX(
******************************************************************************
*             DTIME,
*             THKNES,
*             NBINS,TYPES,NWORK,
*             TAIR,PAIR,PPWV,PPNA,
*             ND,
*             MCS,MCN,MCW,
*             SDND,SDMCS,SDMCN,SDMCW,
*             PRESENCE,SEDMNT,
*             PTSIZE,WORK)
*
*     Declaration of input/output variables:
*     INTEGER NBINS,TYPES,NWORK
*     REAL(KIND=4)
*             DTIME,
*             THKNES,
*             TAIR,PAIR,PPWV,PPNA,
*             ND(NBINS,TYPES),
*             MCS(NBINS,TYPES),MCN(NBINS,TYPES),MCW(NBINS,TYPES),
*             SDND(NBINS,TYPES),
*             SDMCS(NBINS,TYPES),SDMCN(NBINS,TYPES),SDMCW(NBINS,TYPES),
*             PTSIZE(NBINS,8),WORK(NBINS,NWORK)
*     LOGICAL PRESENCE(TYPES),SEDMNT(TYPES)
*
*     Input:
*         DTIME:   (s)            Integration time step
*         THKNES:  (m)            Vertical thickness of layer
*         NBINS:                  Number of particle radii bins
*         TYPES:                  Number of particle types comprehended by the model
*         NWORK:                  2nd dimension of work array
*         TAIR:    (K)            Ambient air temperature
*         PAIR:    (Pa)           Ambient air pressure
*
*     Input/output:
*         PPWV:    (Pa)           Ambient partial pressure of water vapor
*         PPNA:    (Pa)           Ambient partial pressure of nitric acid
*         ND:      (par/kg air)   Number density of particles per kg of air
*         MCS:     (kg/particle)  Mass of condensed sulfuric acid per particle
*         MCN:     (kg/particle)  Mass of condensed nitric acid per particle
*         MCW:     (kg/particle)  Mass of condensed water per particle
*         SDND:    (par/m**3)     Sedimentation flow of particles in time interval DTIME
*         SDMCS:   (kg/particle)  Mass of condensed sulfuric acid in sedimenting particles
*         SDMCN:   (kg/particle)  Mass of condensed nitric acid in sedimenting particles
*         SDMCW:   (kg/particle)  Mass of condensed water in sedimenting particles
*         PRESENCE:(logical)      Indicator of presence of different types of particles
*         SEDMNT:  (logical)      Indicator of sedimentation of different types of particles
*
*     Work array:
*         PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume, etc.
*                                 (not to be changed)
*         WORK:                   Work array (allowed to be changed)
*
******************************************************************************
*
*     This subroutine calculates the particle size distributions and chemical
*     compositions in an ensemble of liquid supercooled ternary solution (STS)
*     stratospheric particles (type 1b PSC at low temperatures), frozen stratospheric
*     sulfate aerosols, and polar stratospheric clouds (PSC) of type 1a and 2 in a
*     SINGLE point in space (box) and time. Individual PSC-box models can be placed in
*     grid-points of a larger model, or can be used for air parcel trajectory
*     calculations, or can be stacked to form a column model. The placement of the
*     box-model within a larger model (geographical/vertical) is referred to as the
*     location of the box-model.
*
*     At each location the PSC-box model takes as input the ambient air state
*     variables: temperature (TAIR), pressure (PAIR), partial pressure of water
*     vapor (PPWV), and partial pressure of nitric acid vapor (PPNA).
*
*     The particle size distribution of each particle type is divided into a number
*     of bins (NBINS) on a geometrically increasing volume scale. It is required,
*     that NBINS > 2. NBINS must be the same for all locations of the PSC-box model.
*
*     The number density (ND - particles per unit air mass), i.e. the
*     the particle size distributions, are stored in the following arrays:
*     ND(i,1), i=1,...NBINS    Liquid sulfate aerosols which, at low temperatures,
*                                 take up nitric acid and water, turning into
*                                 supercooled ternary solution Type 1b PSC
*     ND(i,2), i=1,...NBINS    Frozen sulfuric acid tetrahydrate particles (SAT)
*     ND(i,3), i=1,...NBINS    Nitric acid trihydrate (NAT) type 1a PSC
*     ND(i,4), i=1,...NBINS    Water ice type 2 PSC
*     ND-array element no. i holds the particles per unit air mass in the particle
*     radius interval spanned by bin no. i.
*
*     The subroutine calculates the differential size distribition ND(i,1)
*     of liquid ternary sulfuric-nitric acid solution aerosol particles,
*     assuming the ternary aerosols to be in equilibrium with the water and nitric
*     acid vapor at an arbitrary ambient atmospheric state, specified as input
*     (TAIR, PPWV, and PPNA).
*
*     For each type (j=1,2,3,4) of particle, the condensed masses are stored as follows:
*        condensed sulfuric acid mass per particle, MCS(i,j), i=1,...NBINS, j=1,2,3,4
*        condensed nitric acid mass per particle,   MCN(i,j), i=1,...NBINS, j=1,2,3,4
*        condensed water per particle,              MCW(i,j), i=1,...NBINS, j=1,2,3,4
*
*     Currently, the phase transition between liquid/solid PSCs is unknown.
*     The model applies homogeneous volume dependent freezing of ice in HNO3/H2O
*     supercooled solutions  according to Koop et al. (2000). Homogeneous freezing
*     of NAT (NAD) in supercooled solution according to Tabazadeh et al. (2002)
*     can also be simulated if the common block variable TABA_COR > 0.0, set in calling program.
*     Upon freezing of a STS particle all the nitric acid together with water is
*     deposited as nitric acid trihydrate (NAT) and all sulfuric acid together with water
*     as sulfuric acid tetrahydrate (SAT). Any unbound water in hydrates are assumed to form
*     ice. Upon heating an ice PSC the ice evaporates, turning into a type 1a NAT PSC.
*     Upon further heating the NAT will evaporate, whereby the frozen sulfate aerosol
*     core is released (SAT particle).
*
*     Any frozen SAT particles are assumed to melt when the air temperature
*     is above the SAT melting temperature. SAT melting (deliquescence) upon cooling
*     (Koop and Carslaw, 1996) can also be simulated if the logical parameter SAT_DELI
*     is be set to .TRUE. When melting a SAT paticle the number density is transferred
*     from the frozen to the liquid particle category and the equilibrium composition
*     and ternary density is recalculated.
*
*     The SAT particle may also act as nucleation center for type 1a NAT PSC formation
*     upon subsequent cooling, simulating preactivated SAT (Zhang et al., 1996).
*
*     HNO3 condensation on ice is included if logical parameter HNO3_ON_ICE
*     is set to .TRUE.
*
*     The logical input/output variables PRESENCE(j), (j=1,2,3,4) can be used
*     in the calling program to test whether particles of a given type are present.
*
*     The number density of particles (ND), and the masses of condensed substances
*     (MCS, MCN, MCW)) constitute the integration variables of the model.
*     When the subroutine is called, these variables must hold the
*     values, calculated by the routine in the previous time step (or the initial
*     values, cf. below). Upon exit these variables will hold the new values at a
*     time DTIME later. The new values are calculated by a first order explicit
*     Euler expression within the routine. The ambient air state variables (TAIR,
*     PAIR) are assumed to be constant during the time interval DTIME. In order to
*     perform calculations at different geographical locations and altitudes (grid
*     points), these variables must be defined at each point.
*
*     Nitric acid and water are taken up or released from/to the air during
*     condensation/evaporation. On exit new values of gas phase partial pressures
*     of nitric acid (PPNA) and water (PPWV) are recalculated.
*
*     The balance calculations of particle number density and masses of condensed
*     substance use as input the sedimentation flow of particles (SDNA).
*     The variables SDMCS, SDMCN, AND SDMCW hold the masses of condensed H2SO4, HNO3,
*     and H2O of the particles falling in from the layer above. The routine will return
*     as output, in the same variables, the flow to be used in the box-calculations
*     in the layer below.
*     The logical array SEDMNT(j), (j=1,2,3,4) indicate, if there is a
*     sedimentation flow of particles of type j from the layer above on entry,
*     and indicate, in the same variable, if there is a sedimentation
*     flow to the layer below on exit.
*
*     Thus, for the calculations in a column, the subroutine should be called
*     in a sequence from the top-layer to the bottom-layer. At the top layer,
*     the sedimentation input variables should be set equal to zero; the logical
*     sedimentation indicator (SEDMNT) should be set  to .FALSE., and the calling
*     program should not change the values of  these variables between calls in
*     the column sequence. The sedimentation variables will hold the fall out of
*     particles in the bottom layer after the last call of the subroutine in the
*     column sequence. The same array for storing the sedimentation flows and
*     logical indicators can be used at different locations, if the calculations
*     are performed column-wise.
*
*     A positive vertical extent (THKNES) of the PSC-box (vertical distance
*     between layers) must be given as input to the subroutine. This extent
*     need not be the same, neither at all locations, nor at all times.
*
*     The subroutine uses adjustable array dimensions of all arrays.
*     The parameter NBINS and the arrays in the argument list of the
*     subroutine must be set and dimensioned (cf. above) in the calling
*     program.
*
*     The minimum and maximum particle radius are given as parameters to
*     SUBROUTINE SETBIN, which must be called once before any call to PSCBOX
*     (cf. below). Subroutine SETBIN will store the values of the particle radius,
*     surface area, and volume of each bin in the array PTSIZE. This array PTSIZE
*     must be used at all locations, and the array must not be changed between the
*     calls of subroutine PSCBOX.
*
*     Subroutine PSCBOX_START can be used to initialise PSCBOX, including a call
*     to subroutine SETBIN. PSCBOX_START can be used to calculate initial values
*     of particle number density and chemical composition of a liquid stratospheric
*     sulfuric acid particle ensemble with an initial log-normal size distribution
*     (parameters given to the suboutine). Under warmer stratospheric conditions
*     with no PSCs this distribution could  normally be used as initial values of
*     the PSC-box model. In addition PSCBOX_START will initialise sedimentation
*     variables.
*
*     The total particle surface area density of a given particle type can
*     be calculated by SUBROUTINE SURFCE. The condensed phase mixing ratios 
*     of sulfuric acid, nitric acid, and water can be calculated by SUBROUTINE MIXCON.
*     Parameters of log-normal distribution fits to individual particle
*     distribitions can be calculated by SUBROUTINE LGNPAR (cf. below).
*
*     The array WORK can be be used by the calling program unit between calls
*     of the subroutine, and the same array WORK can be used at all locations.
*
*     The PSC-box model uses a number of named common blocks, holding values
*     which can be also be utilised in the calling program unit:
*
*
*     COMMON /TAU/TAUMAX
*     holding the following REAL value:
*     TAUMAX:      Inverse of minimum time constant encountered in the 
*                  microphysical processes (1/s);
*                  Subroutine PSCBOX will integrate the state variables
*                  forward in time and return after DTIME seconds. During
*                  DTIME the air temperature and pressure are assumed to
*                  be constant. PSCBOX will try to integrate the state
*                  variables in one integration step of DTIME seconds;
*                  however, if the microphysical minimum time constants
*                  encountered are smaller than DTIME, the actual
*                  integration step will be reduced to 1/TAUMAX in order
*                  to make a stable integration; in this case PSCBOX will
*                  use more than one step to integrate DTIME seconds forward.
*                  A minimum time step DTMIN, defined as a parameter below,
*                  is used.
*
*     COMMON /SED/SEDI
*     holding the following LOGICAL value:
*     SEDI:        A logical flag to indecate, if particle sedimentation
*                  calculations are to be performed by PSCBOX. This flag
*                  will be set to .TRUE. initially by subroutine SETBIN;
*                  however, this value can be changed to .FALSE. by the
*                  calling program unit (after the call to SETBIN), if
*                  sedimentation calculations are not wanted.
*
*
*     COMMON /CONTRL/TEST
*     holding the following INTEGER value:
*     TEST:        The value of TEST indicates, if values of local
*                  variables in various sub-programmes, used in the PSC-box
*                  model, should be printed to a file. This facility is
*                  mainly used to test the program. If TEST is equal to
*                  zero, no printing will take place. TEST is initially
*                  set to zero by subroutine SETBIN. If this detailed
*                  (voluminous) output is needed, TEST must be set to
*                  a positive value ( >= 50 ) in the calling program
*                  unit, and a formatted file with unit number TEST must
*                  be opened from the calling program unit. Values of
*                  local variables will then be printed on file unit TEST
*                  as long as TEST is non-zero. The file must also be closed
*                  from the calling program unit. TEST must be >= 60 for
*                  very detailed output, and 50 <= TEST < 60 for moderate
*                  detailed output.
*
*
*     COMMON /BINS/VR,ZERO,D1,LOGVR
*     holding the following REAL values:
*     VR:          VRATIO
*     ZERO:        Smallest significant number density of particles
*                  in a single bin (1/kg air)
*     D1:          Parameter, related to VRATIO, used in log-normal dist.
*     LOGVR:       ALOG(VR)
*
*     Values in common /BINS/ are set by subroutine SETBIN, and must
*     not be changed by the calling program.
*
*
*     Multiplication factor for Tabazadeh freezing rates (if =0, no freezing):
*     REAL TABA_CORR
*     COMMON /CORRECT/TABA_CORR
*
*     All calculations are performed in SI units.
*
*     The entire PSC-box model is written in ANSI fortran 90.
*
******************************************************************************
      SUBROUTINE PSCBOX(
     +        DTIME,
     +        THKNES,
     +        NBINS,TYPES,NWORK,
     +        TAIR,PAIR,PPWV,PPNA,
     +        ND,
     +        MCS,MCN,MCW,
     +        SDND,SDMCS,SDMCN,SDMCW,
     +        PRESENCE,SEDMNT,
     +        PTSIZE,WORK)
C----------------------------------------------------------------------------
      IMPLICIT NONE
C----------------------------------------------------------------------------
C     Particle types:
      INTEGER STS,SAT,NAT,ICE
      PARAMETER(
C     Liquid supercooled binary/ternary solution particle (SA or Type 1b PSC)
     +  STS=1,
C     Solid sulfuric acid tetrahydrate (SAT) particle; i.e. no HNO3 or excess ice
     +  SAT=2,
C     Solid nitric acid trihydrate (NAT) particle; i.e. NAT and SAT but no excess ice (Type 1a PSC)
     +  NAT=3,
C     Solid ice particle; i.e. holding excess ice, NAT, and SAT (Type 2 PSC)
     +  ICE=4)
C----------------------------------------------------------------------------
C     Input/output variables:
      INTEGER NBINS,TYPES,NWORK
      REAL(KIND=4)
     +        DTIME,
     +        THKNES,
     +        TAIR,PAIR,PPWV,PPNA,
     +        ND(NBINS,TYPES),
     +        MCS(NBINS,TYPES),MCN(NBINS,TYPES),MCW(NBINS,TYPES),
     +        SDND(NBINS,TYPES),
     +        SDMCS(NBINS,TYPES),SDMCN(NBINS,TYPES),SDMCW(NBINS,TYPES),
     +        PTSIZE(NBINS,8),WORK(NBINS,NWORK)
      LOGICAL PRESENCE(TYPES),SEDMNT(TYPES)
C----------------------------------------------------------------------------
C     Common block variables:
      REAL(KIND=4) VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
      REAL(KIND=4) TAUMAX
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
C----------------------------------------------------------------------------
C     Phase transitions and condensation:
c     SAT deliquescence mode; SAT melting upon cooling (Koop and Carslaw, 1996):
      LOGICAL SAT_DELI,HNO3_ON_ICE
      PARAMETER(SAT_DELI=.TRUE.)
C     HNO3 condensation on ice:
      PARAMETER(HNO3_ON_ICE=.false.)
C----------------------------------------------------------------------------
C     Ice freezing mode:
      INTEGER FREEZE_MODE
      PARAMETER(
C     Homogeneous volume dependent freezing of ice in H2SO4/HNO3/H2O solution:
     +    FREEZE_MODE=1)
C     Parameterized volume dependent freezing of ice around T-ice
C     +    FREEZE_MODE=2)
C----------------------------------------------------------------------------
C     Physical constants:
      REAL(KIND=4)
     +     MH2O,MHNO3,MH2SO4,MAIR,FSAT,FNAT,
     +     RONAT,ROSAT,LNA,RGAS,STICE,SMIN,S_SAT_STS
      PARAMETER(
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
C       Molar weight of nitric acid (kg/mole)
     +          MHNO3=63.01E-3,
C       Molar weight of sulfuric acid (kg/mole)
     +          MH2SO4=98.08E-3,
C       Molar weight of air (kg/mole)
     +          MAIR=28.9644E-3,
C       Hydrate weight fractions:
     +          FSAT=4.0E0*MH2O/MH2SO4,FNAT=3.0E0*MH2O/MHNO3,
C       Density of nitric acid trihydrate (kg/m**3) (Taesler et al. 1975)
     +          RONAT=1.621E3,
C       Density of sulfuric acid tetrahydrate (kg/m**3)
C                             (Acta Cryst. B28,1692,1972)
     +          ROSAT=1.589E3,
C       Heat of evaporation of nitric acid / STS (J/kg)
     +          LNA=8.52E+5,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441,
C       Surface tension ice/air (N/m) (Pruppacher & Klett, 1980, p. 121)
     +          STICE=105.0E-3,
C       Min. HNO3 saturation for SAT diss./STS formation (Koop Carslaw 1996)
     +          S_SAT_STS=15.0E0,
C       Minimum saturation ratio for nucleation of type 1a and type 2 PSC:
     +          SMIN=1.1)
C----------------------------------------------------------------------------
C     External thermodynamic functions needed:
      REAL(KIND=4) FPLAIR,VISAIR,DFWVA,CDTAIR,PNANAT,PWVICE,
     +             LSWV,ROSNW,ROICE,TMSAT,STNAT
C----------------------------------------------------------------------------
C     Integration step size constants:
      REAL(KIND=4) DTMIN,DTMIN2
      PARAMETER(
C       Minumum integration stepsize (s):
     +          DTMIN=0.1,
     +          DTMIN2=DTMIN/2.0)
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      INTEGER I,J,LIMIT_NAT,LIMIT_ICE,TYPE
      LOGICAL NUCLEA_NAT_SAT,NUCLEA_ICE_NAT,NUCLEA_ICE_SAT,SEDIMENT
      REAL(KIND=4)
     +     RHO(4),RHOAIR,MFPAIR,DYNVIS,DWMAIR,KAIR,
     +     PSNNAT,SNANAT,PSWICE,SWVICE,LSW,
     +     FRZRAT_SNW_NAT,FRZRAT_SNW_ICE,
     +     WSA,WNA,
     +     X0,X1,X2,
     +     C1,C2,C3,C4,C5,CWV,CNA,
     +     PPPWV,PPPNA,FFWV,FFNA,MIXWV,MIXNA,
     +     MSA1,MSA2,MMRSA,
     +     TIME,DT,TAU,
     +     MN1,MN2,MW1,MW2,
     +     FREEZE
C
      PARAMETER(
     +        C1=MHNO3/(3.0*MH2O),
     +        C2=(3.0*MH2O+MHNO3)/MHNO3,
     +        C3=MAIR/RGAS,
     +        C4=(1.0+3.0*MH2O/MHNO3)/RONAT,
     +        C5=(1.0+4.0*MH2O/MH2SO4)/ROSAT,
     +        CWV=MAIR/MH2O,
     +        CNA=MAIR/MHNO3)
C----------------------------------------------------------------------------
      PPPNA=PPNA
      PPPWV=PPWV
C----------------------------------------------------------------------------
C     Thermodynamic calculations:
C----------------------------------------------------------------------------
C     Air:
      RHOAIR=C3*PAIR/TAIR
      MFPAIR=FPLAIR(TAIR,PAIR)
      DYNVIS=VISAIR(TAIR)
      DWMAIR=DFWVA(TAIR,PAIR)
      KAIR=CDTAIR(TAIR)
C     Vapor:
      PSNNAT=PNANAT(TAIR,PPPWV)
      SNANAT=PPPNA/PSNNAT
C     Condensed phase:
C     Calculate mass mixing ratio of condensed sulfuric acid,
C     acid weight fractions in liquid particles, and densities:
      IF(PRESENCE(STS)) THEN
          MMRSA=0.0
          X0=0.0
          X1=0.0
          DO I=1,NBINS
              MMRSA=MMRSA+ND(I,STS)*MCS(I,STS)
              X0=X0+(MCS(I,STS)+MCN(I,STS)+MCW(I,STS))*ND(I,STS)
              X1=X1+ND(I,STS)*MCN(I,STS)
          ENDDO
          WSA=MMRSA/X0
          WNA=X1/X0
      ELSE
          CALL WGTTER(TAIR,PPPWV,PPPNA,WSA,WNA)
      ENDIF
      RHO(STS)=ROSNW(TAIR,WSA,WNA)
      RHO(SAT)=ROSAT
      RHO(NAT)=RONAT
      RHO(ICE)=ROICE(TAIR)
      LSW=LSWV(TAIR)
C----------------------------------------------------------------------------
C     Sedimentation bin-flow:
C----------------------------------------------------------------------------
      DO TYPE=1,TYPES
          SEDIMENT=PRESENCE(TYPE).AND.SEDI
          IF(SEDIMENT) THEN
C             Fall velocity::
              CALL TMVELO(NBINS,TYPE,TAIR,PAIR,RHOAIR,
     +                    MFPAIR,DYNVIS,RHO(TYPE),ND(1,TYPE),
     +                    PTSIZE,WORK(1,1),WORK(1,2))
C             Particle loss terms due to sedimentation to layer below.
              DO I=1,NBINS
                  X1=AMIN1(ND(I,TYPE)*DTIME*WORK(I,1)/THKNES,ND(I,TYPE))
                  ND(I,TYPE)=ND(I,TYPE)-X1
                  WORK(I,1)=X1*RHOAIR
                  WORK(I,3)=MCS(I,TYPE)
                  WORK(I,4)=MCN(I,TYPE)
                  WORK(I,5)=MCW(I,TYPE)
              ENDDO
              PRESENCE(TYPE)=.FALSE.
              DO I=1,NBINS
                  PRESENCE(TYPE)=PRESENCE(TYPE).OR.ND(I,TYPE).GT.ZERO
              END DO
          ENDIF
C         Production of particles by inflow from layer above:
          IF(SEDMNT(TYPE).AND.SEDI) THEN
              PRESENCE(TYPE)=.FALSE.
              DO I=1,NBINS
                  X0=ND(I,TYPE)
                  SDND(I,TYPE)=SDND(I,TYPE)/RHOAIR
                  ND(I,TYPE)=ND(I,TYPE)+SDND(I,TYPE)
                  IF(ND(I,TYPE).GT.ZERO) THEN
                      X1=X0/ND(I,TYPE)
                      X2=SDND(I,TYPE)/ND(I,TYPE)
                      MCS(I,TYPE)=MCS(I,TYPE)*X1+SDMCS(I,TYPE)*X2
                      MCN(I,TYPE)=MCN(I,TYPE)*X1+SDMCN(I,TYPE)*X2
                      MCW(I,TYPE)=MCW(I,TYPE)*X1+SDMCW(I,TYPE)*X2
                  ENDIF
                  PRESENCE(TYPE)=PRESENCE(TYPE).OR.ND(I,TYPE).GT.ZERO
              ENDDO
          ENDIF
C         Sedimentation output flow variabels to layer below:
          SEDMNT(TYPE)=SEDIMENT
          IF(SEDMNT(TYPE)) THEN
              DO I=1,NBINS
                  SDND(I,TYPE)=WORK(I,1)
                  SDMCS(I,TYPE)=WORK(I,3)
                  SDMCN(I,TYPE)=WORK(I,4)
                  SDMCW(I,TYPE)=WORK(I,5)
              ENDDO
          ENDIF
      ENDDO
C----------------------------------------------------------------------------
C     Melting of any frozen particles above SAT melting temperature:
C----------------------------------------------------------------------------
      IF(PRESENCE(SAT)) THEN
C         Melting condition: TAIR > melting temperature of SAT:
          IF(TAIR.GT.TMSAT(PPPWV)) THEN
              PRESENCE(STS)=.FALSE.
              PRESENCE(SAT)=.FALSE.
              MMRSA=0.0
              DO I=1,NBINS
                  X0=ND(I,STS)
                  ND(I,STS)=ND(I,STS)+ND(I,SAT)
                  IF(ND(I,STS).GT.ZERO) THEN
                      X1=X0/ND(I,STS)
                      X2=ND(I,SAT)/ND(I,STS)
                      MCS(I,STS)=MCS(I,STS)*X1+MCS(I,SAT)*X2
                      MCN(I,STS)=MCN(I,STS)*X1+MCN(I,SAT)*X2
                      MCW(I,STS)=MCW(I,STS)*X1+MCW(I,SAT)*X2
                  ENDIF
                  MMRSA=MMRSA+ND(I,STS)*MCS(I,STS)
                  ND(I,SAT)=0.0
                  MCS(I,SAT)=0.0
                  MCN(I,SAT)=0.0
                  MCW(I,SAT)=0.0
C                  PRESENCE(STS)=PRESENCE(STS).OR.ND(I,STS).GT.ZERO
                  PRESENCE(STS)=.TRUE.
              ENDDO
              CALL WGTTER(TAIR,PPPWV,PPPNA,WSA,WNA)
              RHO(STS)=ROSNW(TAIR,WSA,WNA)
          ENDIF
      ENDIF
C----------------------------------------------------------------------------
C     Melting upon cooling of SAT particles; changing into type 1b PSC:
C     (Koop and Carslaw, 1996)
C----------------------------------------------------------------------------
      IF(SAT_DELI.AND.PRESENCE(SAT)) THEN
          IF(SNANAT.GT.S_SAT_STS) THEN
              PRESENCE(STS)=.FALSE.
              PRESENCE(SAT)=.FALSE.
              DO I=1,NBINS
                  X0=ND(I,STS)
                  ND(I,STS)=ND(I,STS)+ND(I,SAT)
                  IF(ND(I,STS).GT.ZERO) THEN
                      X1=X0/ND(I,STS)
                      X2=ND(I,SAT)/ND(I,STS)
                      MCS(I,STS)=MCS(I,STS)*X1+MCS(I,SAT)*X2
                      MCN(I,STS)=MCN(I,STS)*X1+MCN(I,SAT)*X2
                      MCW(I,STS)=MCW(I,STS)*X1+MCW(I,SAT)*X2
                  ENDIF
                  ND(I,SAT)=0.0
                  MCS(I,SAT)=0.0
                  MCN(I,SAT)=0.0
                  MCW(I,SAT)=0.0
                  PRESENCE(STS)=PRESENCE(STS).OR.ND(I,STS).GT.ZERO
              ENDDO
          ENDIF
      ENDIF
C----------------------------------------------------------------------------
C     Growth of liquid particles and uptake of HNO3:
C----------------------------------------------------------------------------
      IF(PRESENCE(STS)) THEN
C         Calculate equilibrium composiotion:
          MSA1=WSA*RHO(STS)
          CALL CMPTER(TAIR,PAIR,MMRSA,PPPWV,PPPNA,WSA,WNA)
          RHO(STS)=ROSNW(TAIR,WSA,WNA)
          MSA2=WSA*RHO(STS)
          X1=ALOG(MSA1/MSA2)/LOGVR
          IF(X1.NE.1.0) THEN
C             Redistribute particles in size bins:
              DO I=1,NBINS
                  WORK(I,1)=I+X1
                  WORK(I,3)=REAL(I)
                  ND(I,STS)=ALOG(AMAX1(ND(I,STS),1.0E-20))
              ENDDO
              CALL SPLINE(WORK(1,1),ND(1,STS),NBINS,WORK(1,2),
     +                    WORK(1,4))
              DO I=1,NBINS
                  WORK(I,4)=AMAX1(-5.0,AMIN1(WORK(I,4),5.0))
              ENDDO
              CALL SPLVEC(NBINS,WORK(1,1),ND(1,STS),WORK(1,4),
     +                    X1,WORK(1,3),WORK(1,2))
              DO I=1,NBINS
                  ND(I,STS)=EXP(WORK(I,2))
              ENDDO
              PRESENCE(STS)=.FALSE.
              DO I=1,NBINS
                  PRESENCE(STS)=PRESENCE(STS).OR.ND(I,STS).GT.ZERO
                  IF(ND(I,STS).GT.0.0) THEN
                      MCS(I,STS)=WSA*RHO(STS)*PTSIZE(I,3)
                      MCN(I,STS)=WNA*RHO(STS)*PTSIZE(I,3)
                      MCW(I,STS)=(1.0-WSA-WNA)*RHO(STS)*PTSIZE(I,3)
                  ELSE
                      MCS(I,STS)=0.0
                      MCN(I,STS)=0.0
                      MCW(I,STS)=0.0
                  ENDIF
              END DO
          ENDIF
          IF(TEST.GE.60) THEN
              WRITE(TEST,*) 'TESTING PSCBOX- LIQUID PART. REDIST.'
              WRITE(TEST,1010) 'MMRSA','MSA1','MSA2',
     +                         'X1','WSA','WNA','PPPNA','RHO(STS)'
              WRITE(TEST,1000) MMRSA,MSA1,MSA2,X1,WSA,WNA,PPPNA,RHO(STS)
          ENDIF

      ENDIF
C---------------------------------------------------------------------------------
C     Fast nucleation and condensation/evaporation processes:
C---------------------------------------------------------------------------------
c     Initial settings for time loop:
      TIME=0.0
c     Gasphase:
      MIXWV=PPPWV/(PAIR*CWV)
      MIXNA=PPPNA/(PAIR*CNA)
      PSNNAT=PNANAT(TAIR,PPPWV)
      SNANAT=PPPNA/PSNNAT
      PSWICE=PWVICE(TAIR)
      SWVICE=PPPWV/PSWICE
c     Condensed mass:
      MN1=0.0
      MW1=0.0
      DO J=1,TYPES
          DO I=1,NBINS
              MN1=MN1+MCN(I,J)*ND(I,J)
              MW1=MW1+MCW(I,J)*ND(I,J)
          END DO
      END DO
C     NAT particle fall velocity and ventilation factors:
      CALL TMVELO(NBINS,NAT,TAIR,PAIR,RHOAIR,
     +                    MFPAIR,DYNVIS,RHO(NAT),ND(1,NAT),
     +                    PTSIZE,WORK(1,3),WORK(1,1))
C     Ice particle fall velocity and ventilation factors:
      CALL TMVELO(NBINS,ICE,TAIR,PAIR,RHOAIR,
     +                    MFPAIR,DYNVIS,RHO(ICE),ND(1,ICE),
     +                    PTSIZE,WORK(1,3),WORK(1,2))
C---------------------------------------------------------------------------------
c     Time loop starts here:
C---------------------------------------------------------------------------------
      TAUMAX=1.0/DTIME
      DO WHILE(DTIME-TIME.GT.DTMIN2)
c         Initialise work array (not ventilation factors)
          DO J=3,NWORK
              DO I=1,NBINS
                  WORK(I,J)=0.0
              ENDDO
          ENDDO
C---------------------------------------------------------------------------------
C         Homogeneous freezing rate of NAD/NAT in H2SO4/HNO3/H2O solution:
C---------------------------------------------------------------------------------
          FRZRAT_SNW_NAT=-1.0
          IF(PRESENCE(STS).AND.SNANAT.GT.SMIN.AND.
     +       TABA_CORR.GT.0.0D0) THEN
CC            CALL HOMFRZ_SNW_NAT(TAIR,PPPWV,WSA,WNA,FRZRAT_SNW_NAT)
CC            CALL HOMFRZ_SNW_NAD(TAIR,PPPWV,WSA,WNA,FRZRAT_SNW_NAT)
              CALL SURF_FRZ_SNW_NAD(TAIR,WSA,WNA,FRZRAT_SNW_NAT)
CC            CALL SURF_FRZ_SNW_NAT(TAIR,WSA,WNA,FRZRAT_SNW_NAT)
              FRZRAT_SNW_NAT=TABA_CORR*FRZRAT_SNW_NAT
c              CALL VOL_FRZ_SNW_NAD(TAIR,PPPWV,PPPNA,WSA,WNA,
c     +                             FRZRAT_SNW_NAT)
              IF (FRZRAT_SNW_NAT.GT.0.0) THEN
                  DO I=1,NBINS
                      IF(ND(I,STS).GT.ZERO)
c     +                   WORK(I,3)=FRZRAT_SNW_NAT*PTSIZE(I,3)
     +                   WORK(I,3)=FRZRAT_SNW_NAT*PTSIZE(I,2)
                  END DO
              END IF
          ENDIF
C-----------------------------------------------------------------------------
C         Homogeneous freezing rate of ice in H2SO4/HNO3/H2O solution:
C-----------------------------------------------------------------------------
          FRZRAT_SNW_ICE=-1.0
          IF(PRESENCE(STS).AND.SWVICE.GT.SMIN) THEN
              IF(FREEZE_MODE.EQ.1) THEN
                  CALL HOMFRZ_SNW_IC(TAIR,PAIR,PPPWV,FRZRAT_SNW_ICE)
              ELSE IF(FREEZE_MODE.EQ.2) THEN
                  FRZRAT_SNW_ICE=FREEZE(TAIR,PPPWV)
              ENDIF
              IF (FRZRAT_SNW_ICE.GT.0.0) THEN
                  DO I=1,NBINS
                      IF(ND(I,STS).GT.ZERO)
     +                 WORK(I,4)=FRZRAT_SNW_ICE*PTSIZE(I,3)
                  END DO
              ENDIF
          ENDIF
C----------------------------------------------------------------------------
C         Heterogeneous nucleation rate of NAT upon SAT:
C----------------------------------------------------------------------------
          NUCLEA_NAT_SAT=.FALSE.
          IF(SNANAT.GT.SMIN.AND.PRESENCE(SAT)) THEN
             CALL HETNUC(3,NBINS,TAIR,PPPNA,SNANAT,RHO(NAT),STNAT(TAIR),
     +                    ND(1,SAT),PTSIZE,WORK(1,5),NUCLEA_NAT_SAT)
          ENDIF
C----------------------------------------------------------------------------
C         Heterogeneous nucleation rate of ice upon NAT:
C----------------------------------------------------------------------------
          NUCLEA_ICE_NAT=.FALSE.
          IF(SWVICE.GT.SMIN.AND.PRESENCE(NAT)) THEN
              CALL HETNUC(1,NBINS,TAIR,PPPWV,SWVICE,RHO(ICE),STICE,
     +                    ND(1,NAT),PTSIZE,WORK(1,6),NUCLEA_ICE_NAT)
          ENDIF
C----------------------------------------------------------------------------
C         Heterogeneous nucleation rate of ice upon SAT:
C----------------------------------------------------------------------------
          NUCLEA_ICE_SAT=.FALSE.
          IF(SWVICE.GT.SMIN.AND.PRESENCE(SAT)) THEN
              CALL HETNUC(1,NBINS,TAIR,PPPWV,SWVICE,RHO(ICE),STICE,
     +                    ND(1,SAT),PTSIZE,WORK(1,7),NUCLEA_ICE_SAT)
          ENDIF
C----------------------------------------------------------------------------
C         Potential nitric acid vapor flux to/from
C         a single NAT PSC particle due to condensation/evaporation:
C----------------------------------------------------------------------------
          LIMIT_NAT=-1
          IF(PRESENCE(NAT)) THEN
              CALL CONDEN(3,NAT,TAIR,PAIR,MFPAIR,KAIR,DWMAIR,
     +                    SNANAT,PSNNAT,RHO(NAT),STNAT(TAIR),LNA,
     +                    NBINS,ND(1,NAT),PTSIZE(1,1),WORK(1,1),
     +                    WORK(1,8),LIMIT_NAT)
              IF(LIMIT_NAT.GE.0) THEN
                  DO I=1,NBINS
                      WORK(I,9)=WORK(I,8)/C1
                      WORK(I,12)=WORK(I,8)*C2/(RHO(NAT)*PTSIZE(I,3))
                  ENDDO
              ENDIF
          ENDIF
C----------------------------------------------------------------------------
C         Potential nitric acid and water vapor flux to/from
C         a single ice PSC particle due to condensation/evaporation:
C----------------------------------------------------------------------------
          LIMIT_ICE=-1
          IF(PRESENCE(ICE)) THEN
              IF(HNO3_ON_ICE) THEN
                  CALL CONDEN(3,NAT,TAIR,PAIR,MFPAIR,KAIR,DWMAIR,
     +                        SNANAT,PSNNAT,RHO(NAT),STNAT(TAIR),LNA,
     +                        NBINS,ND(1,ICE),PTSIZE(1,1),WORK(1,2),
     +                        WORK(1,10),LIMIT_ICE)
              ENDIF
              CALL CONDEN(1,ICE,TAIR,PAIR,MFPAIR,KAIR,DWMAIR,
     +                SWVICE,PSWICE,RHO(ICE),STICE,LSW,
     +                NBINS,ND(1,ICE),PTSIZE(1,1),WORK(1,2),
     +                WORK(1,11),LIMIT_ICE)
              IF(LIMIT_ICE.GE.0) THEN
                  IF(HNO3_ON_ICE) THEN
                      DO I=1,NBINS
                          WORK(I,11)=WORK(I,11)+(WORK(I,10)/C1)
                          WORK(I,13)=
     +                           (WORK(I,10)*C2/(RHO(NAT)*PTSIZE(I,3)))+
     +                           (WORK(I,11)/(RHO(ICE)*PTSIZE(I,3)))
                      ENDDO
                  ELSE
                      DO I=1,NBINS
                          WORK(I,13)=
     +                           (WORK(I,11)/(RHO(ICE)*PTSIZE(I,3)))
                      ENDDO
                  ENDIF
              ENDIF
          ENDIF
C---------------------------------------------------------------------------------
C         Integration step size:
C---------------------------------------------------------------------------------
          TAU=1.0/DTIME
          DO I=1,NBINS
              TAU=AMAX1(TAU,
     +              ABS(WORK(I,3)),
     +              ABS(WORK(I,4)),
     +              ABS(WORK(I,5)),
     +              ABS(WORK(I,6)),
     +              ABS(WORK(I,7)),
     +              ABS(WORK(I,12)),
     +              ABS(WORK(I,13)))
          ENDDO
          DT=AMIN1(AMAX1(1.0/TAU,DTMIN),DTIME-TIME)
          TAUMAX=AMAX1(TAU,TAUMAX)
C-------------------------------------------------------------------------------------
C         Homogeneous volume dependent freezing of NAD/NAT in H2SO4/HNO3/H2O solution:
C-------------------------------------------------------------------------------------
          IF(FRZRAT_SNW_NAT.GT.0.0) THEN
              PRESENCE(STS)=.FALSE.
              PRESENCE(NAT)=.FALSE.
              DO I=1,NBINS
                  IF(ND(I,STS).GT.ZERO) THEN
                      X2=AMIN1(ND(I,STS)*WORK(I,3)*DT,ND(I,STS))
                      IF(X2.GT.ZERO) then
                          X1=ND(I,NAT)
                          ND(I,STS)=AMAX1(ND(I,STS)-X2,0.0)
                          ND(I,NAT)=ND(I,NAT)+X2
                          IF(ND(I,NAT).GT.ZERO) THEN
                              X1=X1/ND(I,NAT)
                              X2=X2/ND(I,NAT)
                              MCS(I,NAT)=MCS(I,NAT)*X1+MCS(I,STS)*X2
                              MCN(I,NAT)=MCN(I,NAT)*X1+MCN(I,STS)*X2
                              MCW(I,NAT)=MCW(I,NAT)*X1+MCW(I,STS)*X2
                          ENDIF
                      ENDIF
                  ENDIF
                  PRESENCE(STS)=PRESENCE(STS).OR.ND(I,STS).GT.ZERO
                  PRESENCE(NAT)=PRESENCE(NAT).OR.ND(I,NAT).GT.ZERO
              ENDDO
          ENDIF
C---------------------------------------------------------------------------------
C         Homogeneous volume dependent freezing of ice in H2SO4/HNO3/H2O solution:
C---------------------------------------------------------------------------------
          IF(FRZRAT_SNW_ICE.GT.0.0) THEN
              PRESENCE(STS)=.FALSE.
              PRESENCE(ICE)=.FALSE.
              DO I=1,NBINS
                  IF(ND(I,STS).GT.ZERO) THEN
                      X2=AMIN1(ND(I,STS)*WORK(I,4)*DT,ND(I,STS))
                      IF(X2.GT.ZERO) THEN
                          X1=ND(I,ICE)
                          ND(I,STS)=AMAX1(ND(I,STS)-X2,0.0)
                          ND(I,ICE)=ND(I,ICE)+X2
                          IF(ND(I,ICE).GT.ZERO) THEN
                               X1=X1/ND(I,ICE)
                               X2=X2/ND(I,ICE)
                               MCS(I,ICE)=MCS(I,ICE)*X1+MCS(I,STS)*X2
                               MCN(I,ICE)=MCN(I,ICE)*X1+MCN(I,STS)*X2
                               MCW(I,ICE)=MCW(I,ICE)*X1+MCW(I,STS)*X2
                          ENDIF
                      ENDIF
                  ENDIF
                  PRESENCE(STS)=PRESENCE(STS).OR.ND(I,STS).GT.ZERO
                  PRESENCE(ICE)=PRESENCE(ICE).OR.ND(I,ICE).GT.ZERO
              ENDDO
          ENDIF
C----------------------------------------------------------------------------
C         Heterogeneous nucleation of NAT upon SAT:
C----------------------------------------------------------------------------
          IF(NUCLEA_NAT_SAT) THEN
C             Loss of SAT particles by NUCLEATION to NAT PSC:
              PRESENCE(SAT)=.FALSE.
              DO I=1,NBINS
                  X2=AMIN1(WORK(I,5)*ND(I,SAT)*DT,ND(I,SAT))
                  IF(X2.GT.0.0) THEN
                       ND(I,SAT)=AMAX1(0.0,ND(I,SAT)-X2)
C                      Production of NAT PSC by NUCLEATION:
                       X1=ND(I,NAT)
                       ND(I,NAT)=ND(I,NAT)+X2
                       IF(ND(I,NAT).GT.ZERO) THEN
                           X1=X1/ND(I,NAT)
                           X2=X2/ND(I,NAT)
                           MCS(I,NAT)=MCS(I,NAT)*X1+MCS(I,SAT)*X2
                           MCN(I,NAT)=MCN(I,NAT)*X1+MCN(I,SAT)*X2
                           MCW(I,NAT)=MCW(I,NAT)*X1+MCW(I,SAT)*X2
                           PRESENCE(NAT)=.TRUE.
                       ENDIF
                  ENDIF
                  PRESENCE(SAT)=PRESENCE(SAT).OR.ND(I,SAT).GT.ZERO
              ENDDO
          ENDIF
C----------------------------------------------------------------------------
C         Heterogeneous nucleation of ice upon NAT:
C----------------------------------------------------------------------------
          IF(NUCLEA_ICE_NAT) THEN
C             Loss of NAT particles by NUCLEATION to ice PSC:
              PRESENCE(NAT)=.FALSE.
              DO I=1,NBINS
                  X2=AMIN1(WORK(I,6)*ND(I,NAT)*DT,ND(I,NAT))
                  IF(X2.GT.ZERO) THEN
                       ND(I,NAT)=AMAX1(0.0,ND(I,NAT)-X2)
C                      Production of ice PSC by NUCLEATION:
                       X1=ND(I,ICE)
                       ND(I,ICE)=ND(I,ICE)+X2
                       IF(ND(I,ICE).GT.ZERO) THEN
                           X1=X1/ND(I,ICE)
                           X2=X2/ND(I,ICE)
                           MCS(I,ICE)=MCS(I,ICE)*X1+MCS(I,NAT)*X2
                           MCN(I,ICE)=MCN(I,ICE)*X1+MCN(I,NAT)*X2
                           MCW(I,ICE)=MCW(I,ICE)*X1+MCW(I,NAT)*X2
                           PRESENCE(ICE)=.TRUE.
                       ENDIF
                  ENDIF
                  PRESENCE(NAT)=PRESENCE(NAT).OR.ND(I,NAT).GT.ZERO
              ENDDO
          ENDIF
C----------------------------------------------------------------------------
C         Heterogeneous nucleation of ice upon SAT:
C----------------------------------------------------------------------------
          IF(NUCLEA_ICE_SAT) THEN
C             Loss of SAT particles by NUCLEATION to ice PSC:
              PRESENCE(SAT)=.FALSE.
              DO I=1,NBINS
                  X2=AMIN1(WORK(I,7)*ND(I,SAT)*DT,ND(I,SAT))
                  IF(X2.GT.ZERO) THEN
                       ND(I,SAT)=AMAX1(0.0,ND(I,SAT)-X2)
C                      Production of ice PSC by NUCLEATION:
                       X1=ND(I,ICE)
                       ND(I,ICE)=ND(I,ICE)+X2
                       IF(ND(I,ICE).GT.ZERO) THEN
                           X1=X1/ND(I,ICE)
                           X2=X2/ND(I,ICE)
                           MCS(I,ICE)=MCS(I,ICE)*X1+MCS(I,SAT)*X2
                           MCN(I,ICE)=MCN(I,ICE)*X1+MCN(I,SAT)*X2
                           MCW(I,ICE)=MCW(I,ICE)*X1+MCW(I,SAT)*X2
                           PRESENCE(ICE)=.TRUE.
                       ENDIF
                  ENDIF
                  PRESENCE(SAT)=PRESENCE(SAT).OR.ND(I,SAT).GT.ZERO
              ENDDO
          ENDIF
C----------------------------------------------------------------------------
C         Condensation/evaporation bin flows of NAT particles:
C----------------------------------------------------------------------------
          IF(LIMIT_NAT.GE.0) THEN
C             Mass condensation/evaporation:
              DO I=1,NBINS
                  MCN(I,NAT)=AMAX1(MCN(I,NAT)+WORK(I,8)*DT,0.0)
                  MCW(I,NAT)=AMAX1(MCW(I,NAT)+WORK(I,9)*DT,0.0)
              END DO
C             Growing NAT particles bin flow:
              DO I=LIMIT_NAT+1,NBINS
                  X0=AMIN1(WORK(I,12)*ND(I,NAT)*DT,ND(I,NAT))
                  ND(I,NAT)=ND(I,NAT)-X0
                  WORK(I,3)=MCS(I,NAT)
                  WORK(I,4)=MCN(I,NAT)
                  WORK(I,5)=MCW(I,NAT)
                  WORK(I,6)=X0
              ENDDO
              DO I=LIMIT_NAT+2,NBINS
                  X0=ND(I,NAT)
                  ND(I,NAT)=X0+WORK(I-1,6)
                  IF(ND(I,NAT).GT.ZERO) THEN
                     X1=X0/ND(I,NAT)
                     X2=WORK(I-1,6)/ND(I,NAT)
                     MCS(I,NAT)=AMAX1(MCS(I,NAT)*X1+WORK(I-1,3)*X2,0.0)
                     MCN(I,NAT)=AMAX1(MCN(I,NAT)*X1+WORK(I-1,4)*X2,0.0)
                     MCW(I,NAT)=AMAX1(MCW(I,NAT)*X1+WORK(I-1,5)*X2,0.0)
                  ENDIF
              ENDDO
C             Evaporating NAT particles bin flow:
              DO I=LIMIT_NAT,1,-1
                  X0=AMIN1(-WORK(I,12)*ND(I,NAT)*DT,ND(I,NAT))
                  ND(I,NAT)=ND(I,NAT)-X0
                  WORK(I,3)=MCS(I,NAT)
                  WORK(I,4)=MCN(I,NAT)
                  WORK(I,5)=MCW(I,NAT)
                  WORK(I,6)=X0
              ENDDO
              DO I=LIMIT_NAT-1,1,-1
                  X0=ND(I,NAT)
                  ND(I,NAT)=X0+WORK(I+1,6)
                  IF(ND(I,NAT).GT.ZERO) THEN
                     X1=X0/ND(I,NAT)
                     X2=WORK(I+1,6)/ND(I,NAT)
                     MCS(I,NAT)=AMAX1(MCS(I,NAT)*X1+WORK(I+1,3)*X2,0.0)
                     MCN(I,NAT)=AMAX1(MCN(I,NAT)*X1+WORK(I+1,4)*X2,0.0)
                     MCW(I,NAT)=AMAX1(MCW(I,NAT)*X1+WORK(I+1,5)*X2,0.0)
                  ENDIF
C                 Core return from NAT to SAT particles:
                  X0=MCS(I,NAT)*(1.0+FSAT)/ROSAT
                  IF(X0.GT.PTSIZE(I,3).AND.ND(I,NAT).GT.ZERO) THEN
                       X0=ND(I,SAT)
                       X1=ND(I,NAT)
                       ND(I,NAT)=0.0
                       ND(I,SAT)=ND(I,SAT)+X1
                       IF(ND(I,SAT).GT.ZERO) THEN
                           X0=X0/ND(I,SAT)
                           X1=X1/ND(I,SAT)
                           MCS(I,SAT)=AMAX1(MCS(I,SAT)*X0+MCS(I,NAT)*X1,
     +                                      0.0)
                           MCW(I,SAT)=MCS(I,SAT)*FSAT
                           PRESENCE(SAT)=.TRUE.
                       ENDIF
                       MCS(I,NAT)=0.0
                       MCN(I,NAT)=0.0
                       MCW(I,NAT)=0.0
                  ENDIF
              ENDDO
              IF(TEST.GT.0) THEN
                  WRITE(TEST,*)
                  WRITE(TEST,*) 'TESTING PSCBOX-NAT'
                  WRITE(TEST,1010) 'TAUMAX','TAU','DT','DTIME','TIME'
                  WRITE(TEST,1000) TAUMAX,TAU,DT,DTIME,TIME
                  WRITE(TEST,1010) 'PPPWV','PPPNA','FFWV','FFNA',
     +                             'MIXWV','MIXNA','SNANAT','TAIR'
                  WRITE(TEST,1000) PPPWV,PPPNA,FFWV,FFNA,
     +                             MIXWV,MIXNA,SNANAT,TAIR
                  WRITE(TEST,*) '(ND(I,NAT),I=1,NBINS)'
                  WRITE(TEST,1000) (ND(I,NAT),I=1,NBINS)
                  WRITE(TEST,*) '(MCS(I,NAT),I=1,NBINS)'
                  WRITE(TEST,1000) (MCS(I,NAT),I=1,NBINS)
                  WRITE(TEST,*) '(MCN(I,NAT),I=1,NBINS)'
                  WRITE(TEST,1000) (MCN(I,NAT),I=1,NBINS)
                  WRITE(TEST,*) '(MCW(I,NAT),I=1,NBINS)'
                  WRITE(TEST,1000) (MCW(I,NAT),I=1,NBINS)
              ENDIF
          ENDIF
C----------------------------------------------------------------------------
C         Condensation/evaporation bin flows of ice particles:
C----------------------------------------------------------------------------
          IF(LIMIT_ICE.GE.0) THEN
              DO I=1,NBINS
                  MCN(I,ICE)=AMAX1(MCN(I,ICE)+WORK(I,10)*DT,0.0)
                  MCW(I,ICE)=AMAX1(MCW(I,ICE)+WORK(I,11)*DT,0.0)
              END DO
C             Growing PSC 2 particles bin flow:
              DO I=LIMIT_ICE+1,NBINS
                  X0=AMIN1(WORK(I,13)*ND(I,ICE)*DT,ND(I,ICE))
                  ND(I,ICE)=ND(I,ICE)-X0
                  WORK(I,3)=MCS(I,ICE)
                  WORK(I,4)=MCN(I,ICE)
                  WORK(I,5)=MCW(I,ICE)
                  WORK(I,6)=X0
              ENDDO
              DO I=LIMIT_ICE+2,NBINS
                  X0=ND(I,ICE)
                  ND(I,ICE)=X0+WORK(I-1,6)
                  IF(ND(I,ICE).GT.ZERO) THEN
                     X1=X0/ND(I,ICE)
                     X2=WORK(I-1,6)/ND(I,ICE)
                     MCS(I,ICE)=AMAX1(MCS(I,ICE)*X1+WORK(I-1,3)*X2,0.0)
                     MCN(I,ICE)=AMAX1(MCN(I,ICE)*X1+WORK(I-1,4)*X2,0.0)
                     MCW(I,ICE)=AMAX1(MCW(I,ICE)*X1+WORK(I-1,5)*X2,0.0)
                  ENDIF
              ENDDO
C             PSC 2 evaporating particles bin flow:
              DO I=LIMIT_ICE,1,-1
                  X0=AMIN1(-WORK(I,13)*ND(I,ICE)*DT,ND(I,ICE))
                  ND(I,ICE)=ND(I,ICE)-X0
                  WORK(I,3)=MCS(I,ICE)
                  WORK(I,4)=MCN(I,ICE)
                  WORK(I,5)=MCW(I,ICE)
                  WORK(I,6)=X0
              ENDDO
              DO I=LIMIT_ICE-1,1,-1
                  X0=ND(I,ICE)
                  ND(I,ICE)=X0+WORK(I+1,6)
                  IF(ND(I,ICE).GT.ZERO) THEN
                     X1=X0/ND(I,ICE)
                     X2=WORK(I+1,6)/ND(I,ICE)
                     MCS(I,ICE)=AMAX1(MCS(I,ICE)*X1+WORK(I+1,3)*X2,0.0)
                     MCN(I,ICE)=AMAX1(MCN(I,ICE)*X1+WORK(I+1,4)*X2,0.0)
                     MCW(I,ICE)=AMAX1(MCW(I,ICE)*X1+WORK(I+1,5)*X2,0.0)
                  ENDIF
C                 Core return from ice to NAT particles:
                  X1=(MCS(I,ICE)*(1.0+FSAT)/ROSAT)+
     +               (MCN(I,ICE)*(1.0+FNAT)/RHO(NAT))
                  IF(X1.GT.PTSIZE(I,3).AND.ND(I,ICE).GT.ZERO) THEN
                       X0=ND(I,NAT)
                       X1=ND(I,ICE)
                       ND(I,ICE)=0.0
                       ND(I,NAT)=ND(I,NAT)+X1
                       IF(ND(I,NAT).GT.ZERO) THEN
                           X0=X0/ND(I,NAT)
                           X1=X1/ND(I,NAT)
                           MCS(I,NAT)=AMAX1(MCS(I,NAT)*X0+MCS(I,ICE)*X1,
     +                                      0.0)
                           MCN(I,NAT)=AMAX1(MCN(I,NAT)*X0+MCN(I,ICE)*X1,
     +                                      0.0)
                           MCW(I,NAT)=MCN(I,NAT)*FNAT
                           PRESENCE(NAT)=.TRUE.
                       ENDIF
                       MCS(I,ICE)=0.0
                       MCN(I,ICE)=0.0
                       MCW(I,ICE)=0.0
                  ENDIF
              ENDDO
              IF(TEST.GT.0) THEN
                  WRITE(TEST,*)
                  WRITE(TEST,*) 'TESTING PSCBOX-ICE'
                  WRITE(TEST,1010) 'TAUMAX','TAU','DT','DTIME','TIME'
                  WRITE(TEST,1000) TAUMAX,TAU,DT,DTIME,TIME
                  WRITE(TEST,1010) 'PPPWV','PPPNA','FFWV','FFNA',
     +                             'MIXWV','MIXNA','SWVICE','TAIR'
                  WRITE(TEST,1000) PPPWV,PPPNA,FFWV,FFNA,
     +                             MIXWV,MIXNA,SWVICE,TAIR
                  WRITE(TEST,*) '(ND(I,ICE),I=1,NBINS)'
                  WRITE(TEST,1000) (ND(I,ICE),I=1,NBINS)
                  WRITE(TEST,*) '(MCN(I,ICE),I=1,NBINS)'
                  WRITE(TEST,1000) (MCN(I,ICE),I=1,NBINS)
                  WRITE(TEST,*) '(MCW(I,ICE),I=1,NBINS)'
                  WRITE(TEST,1000) (MCW(I,ICE),I=1,NBINS)
              ENDIF
          ENDIF
C----------------------------------------------------------------------------
C         Mass balance in gas phase:
C----------------------------------------------------------------------------
          MN2=0.0
          MW2=0.0
          DO J=1,TYPES
              DO I=1,NBINS
                  MN2=MN2+MCN(I,J)*ND(I,J)
                  MW2=MW2+MCW(I,J)*ND(I,J)
              END DO
          END DO
          FFNA=MN2-MN1
          FFWV=MW2-MW1
          MN1=MN2
          MW1=MW2
          MIXWV=AMAX1(MIXWV-FFWV,0.0)
          PPPWV=PAIR*CWV*MIXWV
          MIXNA=AMAX1(MIXNA-FFNA,0.0)
          PPPNA=PAIR*CNA*MIXNA
          PSNNAT=PNANAT(TAIR,PPPWV)
          SNANAT=PPPNA/PSNNAT
          PSWICE=PWVICE(TAIR)
          SWVICE=PPPWV/PSWICE
C----------------------------------------------------------------------------
C         Update time:
C----------------------------------------------------------------------------
          TIME=TIME+DT
      ENDDO
C----------------------------------------------------------------------------
      PPNA=PPPNA
      PPWV=PPPWV
      DO TYPE=1,TYPES
          PRESENCE(TYPE)=.FALSE.
          DO I=1,NBINS
              PRESENCE(TYPE)=PRESENCE(TYPE).OR.ND(I,TYPE).GT.ZERO
          ENDDO
          IF(.NOT.(PRESENCE(TYPE))) THEN
              DO I=1,NBINS
                  ND(I,TYPE)=0.0
                  MCS(I,TYPE)=0.0
                  MCN(I,TYPE)=0.0
                  MCW(I,TYPE)=0.0
              END DO
          ENDIF
      END DO
C----------------------------------------------------------------------------
      RETURN
C----------------------------------------------------------------------------
 1000 FORMAT(1X,8(1PE9.2),0P)
 1010 FORMAT(1X,8(A9))
      END
C
******************************************************************************
*     SUBROUTINE PSCBOX_START(
******************************************************************************
*             RMIN,RMAX,
*             NDTOT,MEDIAN,GSTDEV,
*             MRSA,
*             NBINS,TYPES,
*             TAIR,PAIR,PPWV,PPNA,
*             ND,
*             MCS,MCN,MCW,
*             SDND,SDMCS,SDMCN,SDMCW,
*             PRESENCE,SEDMNT,
*             PTSIZE)
*
*     Declaration of input/output variables:
*     INTEGER NBINS,TYPES,NWORK
*     REAL(KIND=4)
*             RMIN,RMAX,
*             NDTOT,MEDIAN,GSTDEV,
*             TAIR,PAIR,PPWV,PPNA,
*             ND(NBINS,TYPES),
*             MCS(NBINS,TYPES),MCN(NBINS,TYPES),MCW(NBINS,TYPES),
*             SDND(NBINS,TYPES),
*             SDMCS(NBINS,TYPES),SDMCN(NBINS,TYPES),SDMCW(NBINS,TYPES),
*             PTSIZE(NBINS,8))
*     LOGICAL PRESENCE(TYPES),SEDMNT(TYPES)
*
*     Input:
*         RMIN:    (microns)      Minimum particle radius
*         RMAX:    (microns)      Minimum particle radius  in initial state
*         NDTOT:   (par/m**3)     Total number density of particles
*         MEDIAN:  (m)            Lognormal median particle radius
*         GSTDEV:                 Lognormal geometric standard deviation
*         NBINS:                  Number of particle radii bins
*         TYPES:                  Number of particle types comprehended by the model
*         TAIR:    (K)            Ambient air temperature
*         PAIR:    (Pa)           Ambient air pressure
*         PPWV:    (Pa)           Ambient partial pressure of water vapor
*         PPNA:    (Pa)           Ambient partial pressure of nitric acid
*
*     Output:
*         ND:      (par/kg air)   Number density of particles per kg of air
*         MCS:     (kg/particle)  Mass of condensed sulfuric acid per particle
*         MCN:     (kg/particle)  Mass of condensed nitric acid per particle
*         MCW:     (kg/particle)  Mass of condensed water per particle
*         SDND:    (par/m**3)     Sedimentation flow of particles in time interval DTIME
*         SDMCS:   (kg/particle)  Mass of condensed sulfuric acid in sedimenting particles
*         SDMCN:   (kg/particle)  Mass of condensed nitric acid in sedimenting particles
*         SDMCW:   (kg/particle)  Mass of condensed water in sedimenting particles
*         PRESENCE:(logical)      Indicator of presence of different types of particles
*         SEDMNT:  (logical)      Indicator of sedimentation of different types of particles
*         PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume, etc.
*                                 (not to be changed)
*
******************************************************************************
*     Subroutine PSCBOX_START can be used to initialise PSCBOX, including a call
*     to subroutine SETBIN. PSCBOX_START can be used to calculate initial values
*     of particle number density and chemical composition of a liquid stratospheric
*     sulfuric acid particle ensemble with an initial log-normal size distribution
*     (parameters given to the suboutine). Under warmer stratospheric conditions
*     with no PSCs this distribution could  normally be used as initial values of
*     the PSC-box model. In addition PSCBOX_START will initialise sedimentation
*     variables.
*
******************************************************************************
      SUBROUTINE PSCBOX_START(
     +        RMIN,RMAX,
     +        NDTOT,MEDIAN,GSTDEV,
     +        NBINS,TYPES,
     +        TAIR,PAIR,PPWV,PPNA,
     +        ND,
     +        MCS,MCN,MCW,
     +        SDND,SDMCS,SDMCN,SDMCW,
     +        PRESENCE,SEDMNT,
     +        PTSIZE)
C----------------------------------------------------------------------------
      IMPLICIT NONE
C----------------------------------------------------------------------------
C     Particle types:
      INTEGER STS,SAT,NAT,ICE
      PARAMETER(
C     Liquid supercooled binary/ternary solution particle (SA or Type 1b PSC)
     +  STS=1,
C     Solid sulfuric acid tetrahydrate (SAT) particle; i.e. no HNO3 or excess ice
     +  SAT=2,
C     Solid nitric acid trihydrate (NAT) particle; i.e. NAT and SAT but no excess ice (Type 1a PSC)
     +  NAT=3,
C     Solid ice particle; i.e. holding excess ice, NAT, and SAT (Type 2 PSC)
     +  ICE=4)
C----------------------------------------------------------------------------
C     Input/output variables:
      INTEGER NBINS,TYPES,NWORK
      REAL(KIND=4)
     +        RMIN,RMAX,
     +        NDTOT,MEDIAN,GSTDEV,NTOT,
     +        TAIR,PAIR,PPWV,PPNA,
     +        ND(NBINS,TYPES),
     +        MCS(NBINS,TYPES),MCN(NBINS,TYPES),MCW(NBINS,TYPES),
     +        SDND(NBINS,TYPES),
     +        SDMCS(NBINS,TYPES),SDMCN(NBINS,TYPES),SDMCW(NBINS,TYPES),
     +        PTSIZE(NBINS,8)
      LOGICAL PRESENCE(TYPES),SEDMNT(TYPES)
C
      REAL WNA,WSA,RHOTER,ROSNW
      INTEGER I,J
C
C----------------------------------------------------------------------------
C     Physical constants:
      REAL(KIND=4) MAIR,RGAS
      PARAMETER(
C       Molar weight of dry air (kg/mole)
     +          MAIR=28.9644E-3,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441)
C----------------------------------------------------------------------------
C     Calculate particle bin radius, area, volume etc.:
      CALL SETBIN(NBINS,RMIN,RMAX,PTSIZE)
C----------------------------------------------------------------------------
C     Conversion to number of particles per kg air:
      NTOT=NDTOT*RGAS*TAIR/(MAIR*PAIR)
C
C     Calculate log-normal distribution:
      CALL LGNDST(NBINS,NTOT,MEDIAN,GSTDEV,PTSIZE,ND(1,STS))
      PRESENCE(STS)=.TRUE.
C
C     Initial density and H2SO4 weight fraction of sulfuric acid solution:
      CALL WGTTER(TAIR,PPWV,PPNA,WSA,WNA)
      RHOTER=ROSNW(TAIR,WSA,WNA)
C
C     Initial condensed mass in liquid particles:
      DO I=1,NBINS
          MCS(I,STS)=RHOTER*WSA*PTSIZE(I,3)
          MCN(I,STS)=RHOTER*WNA*PTSIZE(I,3)
          MCW(I,STS)=RHOTER*(1.0-WSA-WNA)*PTSIZE(I,3)
      END DO
C
C     Initialize SAT and PSC number densities and condesed mass to zero:
C     ------------------------------------------------------------------
      DO J=SAT,ICE
           DO I=1,NBINS
                ND(I,J)=0.0
                MCS(I,J)=0.0
                MCN(I,J)=0.0
                MCW(I,J)=0.0
           ENDDO
           PRESENCE(J)=.FALSE.
      ENDDO
      DO J=1,TYPES
          DO I=1,NBINS
            SDND(I,J)=0.0
            SDMCS(I,J)=0.0
            SDMCN(I,J)=0.0
            SDMCW(I,J)=0.0
            SEDMNT(J)=.FALSE.
          END DO
      ENDDO
      RETURN
      END
******************************************************************************
*     SUBROUTINE SETBIN(NBINS,RMIN,RMAX,PTSIZE)
******************************************************************************
*
*     Subroutine SETBIN must be called once before any call to
*     subroutines PSCBOX, or INITSP.
*
*     This subroutine calculates the radius, surface area, and volume of
*     each particle size bin on a geometrically increasing volume scale.
*     Further the lower and upper volume bin borders are calculated, and
*     the bin radius-width in æm.
*     These values are stored in array PTSIZE to be used througout the
*     PSC box model, and posssibly also in the calling program unit. The
*     values in PTSIZE must not be changed by the calling program unit.
*     The content of array PTSIZE is as follows:
*         Radius (m):                            PTSIZE(I,1)
*         Particle surface area (m**2):          PTSIZE(I,2)
*         Particle volume (m**3):                PTSIZE(I,3)
*         Lower bin border (m**3):               PTSIZE(I,4)
*         Upper bin border (m**3):               PTSIZE(I,5)
*         Bin radius-width (microns):            PTSIZE(I,6)
*         ln(r):                                 PTSIZE(I,7)
*         (ln(r))**2:                            PTSIZE(I,8)
*
*
*
*     Input:
*     INTEGER NBINS
*     REAL(KIND=4)  RMIN,RMAX
*
*     Output:
*     REAL PTSIZE(NBINS,8)
*
*     Input:
*         NBINS:                  Number of radii groups
*         RMIN:    (microns)      Minimum particle radius
*         RMAX:    (microns)      Minimum particle radius  in initial state
*
*     Output:
*         PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume, etc.
*
*
      SUBROUTINE SETBIN(NBINS,RMIN,RMAX,PTSIZE)
      IMPLICIT NONE
      INTEGER NBINS
      REAL(KIND=4) RMIN,RMAX,VRATIO,PTSIZE(NBINS,8)
C
      REAL(KIND=4) NDMIN
      PARAMETER(
C         Smallest significant number of particles per kg air in each bin:
     +    NDMIN = 1.0E-2)
C
C     Mathematical constants:
      REAL(KIND=4) PI
      PARAMETER(PI=3.1415926536)
C
C     Auxiliary local variables:
      INTEGER I
C
C     Common block variables:
      REAL(KIND=4) VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
      LOGICAL SEDI
      COMMON /SED/SEDI
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C----------------------------------------------------------------------------
      VRATIO=EXP(ALOG(RMAX/RMIN)/(REAL(NBINS-1)/3.0))
      D1=((2.0/(VRATIO+1.0))**(1.0/3.0))*(VRATIO**(1.0/3.0)-1.0)
      LOGVR=ALOG(VRATIO)
      DO 100 I=1,NBINS,1
C         Radius (m):
          PTSIZE(I,1)=RMIN*1.0E-6*VRATIO**((I-1.0)/3.0)
C         Particle surface area (m**2):
          PTSIZE(I,2)=4.0*PI*PTSIZE(I,1)**2
C         Particle volume (m**3):
          PTSIZE(I,3)=4.0*PI*(PTSIZE(I,1)**3)/3.0
C         Lower bin border (m**3):
          PTSIZE(I,4)=PTSIZE(I,3)*2.0/(VRATIO+1.0)
C         Upper bin border (m**3):
          PTSIZE(I,5)=PTSIZE(I,4)*VRATIO
C         Bin radius-width (microns):
          PTSIZE(I,6)=PTSIZE(I,1)*D1*1.0E6
C         ln(r):
          PTSIZE(I,7)=ALOG(PTSIZE(I,1))
C         (ln(r))**2:
          PTSIZE(I,8)=PTSIZE(I,7)**2
 100  CONTINUE
C----------------------------------------------------------------------------
      VR=VRATIO
      ZERO=NDMIN
      D1=D1/SQRT(2.0*PI)
      SEDI=.TRUE.
C----------------------------------------------------------------------------
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING SETBIN','    NBINS = ',NBINS
      WRITE(TEST,1010) 'RMIN','VR','ZERO','D1'
      WRITE(TEST,1000) RMIN,VR,ZERO,D1
      WRITE(TEST,*) '(PTSIZE(I,1),I=1,NBINS)'
      WRITE(TEST,1000) (PTSIZE(I,1),I=1,NBINS)
 1000 FORMAT(1X,8(1PE9.2),0P)
 1010 FORMAT(1X,8(A9))
C----------------------------------------------------------------------------
      RETURN
      END
C
*******************************************************************************
*      SUBROUTINE MIXING_FROM_NUMBER_DENSITY(TAIR,PAIR,NDTOT,PPWV,PPNA,
*                                            MEDIAN,GSTDEV,MRSA)
******************************************************************************
*
*     This subroutine calculates a sulfuric acid volume mixing ratio consistent with
*     the total number density of sulfate aerosol, assuming a given lognormal
*     size distribution.
*
*     Declaration of input/output variables:
*
*
*
*      REAL(KIND=4) TAIR,MRSA,PPWV,PPNA,MEDIAN,GSTDEV,NDTOT
*
*     Input:
*
*         TAIR:   (K)            Air temperature
*         PAIR:   (Pa)           Air pressure
*         NDTOT:  (par/m**3  )   Total number density of particles
*         PPWV:   (Pa)           Water vapor partial pressure
*         PPNA:   (Pa)           Nitric acid mixing ratio
*         MEDIAN: (m)            Median particle radius
*         GSTDEV:                Geometric standard deveation
*     Output:
*         MRSA:                  Volume mixing ratio of sulfuric acid
*
      SUBROUTINE MIXING_FROM_NUMBER_DENSITY(TAIR,PAIR,NDTOT,PPWV,PPNA,
     +                                       MEDIAN,GSTDEV,MRSA)
c
      REAL(KIND=4) TAIR,PAIR,MRSA,PPWV,PPNA,MEDIAN,GSTDEV,NDTOT
      REAL(KIND=4) WSA,WNA,RHOTER,ROSNW,MMRSA,ND,F,RHOAIR

      REAL RGAS,MH2SO4,MAIR,PI,C1
      PARAMETER(
     +          RGAS=8.31441,
     +          MH2SO4=98.08E-3,
     +          MAIR=28.9644E-3,
     +          PI=3.1415926536E0,
     +          C1=4.0E0*PI/3.0E0)
      RHOAIR=MAIR*PAIR/(RGAS*TAIR)
      ND=NDTOT/RHOAIR
      CALL WGTTER(TAIR,PPWV,PPNA,WSA,WNA)
      RHOTER=ROSNW(TAIR,WSA,WNA)
      F=C1*(MEDIAN**3)*EXP(4.5E0*(ALOG(GSTDEV))**2)
      MMRSA=ND*WSA*RHOTER*F
      MRSA=MMRSA/(MH2SO4/MAIR)
      RETURN
      END

*******************************************************************************
*      SUBROUTINE NUMBER_DENSITY_FROM_MIXING(TAIR,PAIR,MRSA,PPWV,PPNA,
*                                            MEDIAN,GSTDEV,NDTOT)
******************************************************************************
*
*     This subroutine calculates the total number density of sulfate aerosol
*     consistent with a sulfuric acid volume mixing ratio, assuming a given lognormal
*     size distribution.
*
*     Declaration of input/output variables:
*
*
*
*      REAL(KIND=4) TAIR,MRSA,PPWV,PPNA,MEDIAN,GSTDEV,NDTOT
*
*     Input:
*
*         TAIR:   (K)            Air temperature
*         PAIR:   (Pa)           Air pressure
*         MRSA:                  Volume mixing ratio of sulfuric acid
*         PPWV:   (Pa)           Water vapor partial pressure
*         PPNA:   (Pa)           Nitric acid mixing ratio
*         MEDIAN: (m)            Median particle radius
*         GSTDEV:                Geometric standard deveation
*     Output:
*         NDTOT:   (par/m**3  )   Total number density of particles
*
      SUBROUTINE NUMBER_DENSITY_FROM_MIXING(TAIR,PAIR,MRSA,PPWV,PPNA,
     +                                      MEDIAN,GSTDEV,NDTOT)
c
      REAL(KIND=4) TAIR,PAIR,MRSA,PPWV,PPNA,MEDIAN,GSTDEV,NDTOT
      REAL(KIND=4) WSA,WNA,RHOTER,ROSNW,MMRSA,ND,F,RHOAIR

      REAL RGAS,MH2SO4,MAIR,PI,C1
      PARAMETER(
     +          RGAS=8.31441,
     +          MH2SO4=98.08E-3,
     +          MAIR=28.9644E-3,
     +          PI=3.1415926536E0,
     +          C1=4.0E0*PI/3.0E0)
      MMRSA=MRSA*MH2SO4/MAIR
      CALL WGTTER(TAIR,PPWV,PPNA,WSA,WNA)
      RHOTER=ROSNW(TAIR,WSA,WNA)
      RHOAIR=MAIR*PAIR/(RGAS*TAIR)
      F=C1*(MEDIAN**3)*EXP(4.5E0*(ALOG(GSTDEV))**2)
      ND=MMRSA/(WSA*RHOTER*F)
      NDTOT=ND*RHOAIR
      RETURN
      END
******************************************************************************
*     SUBROUTINE LGNDST(NBINS,NTOTAL,MEDIAN,GSTDEV,PTSIZE,ND)
******************************************************************************
*
*     This subroutine calculates a log-normal distribution with
*     total number density of particles NTOTAL, median radius MEDIAN,
*     and geometric standard deveation GSTDEV. Subroutine SETBIN must
*     have been called prior to the use of this subroutine.
*
*     Input/output variables:
*     INTEGER NBINS
*     REAL(KIND=4)    NTOTAL,MEDIAN,GSTDEV,ND(NBINS)
*
*     Input:
*         NBINS:                  Number of particle radii bins
*         NTOTAL:  (par/kg air)   Total number density of particles
*         MEDIAN:  (m)            Median particle radius
*         GSTDEV:                 Geometric standard deveation
*         PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume
*
*     Output:
*         ND:      (par/kg air)   Number density of particles
*
      SUBROUTINE LGNDST(NBINS,NTOTAL,MEDIAN,GSTDEV,PTSIZE,ND)
C
      IMPLICIT NONE
      INTEGER NBINS
      REAL(KIND=4)    NTOTAL,MEDIAN,GSTDEV,PTSIZE(NBINS,8),ND(NBINS)
C
C
C     Common block variables:
      REAL(KIND=4) VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
C     Auxiliary local variables:
      INTEGER I
      REAL(KIND=4) D2,D3,D4
C----------------------------------------------------------------------------
      D3=ALOG(GSTDEV)
      D2=NTOTAL*D1/D3
      D3=-0.5/(D3*D3)
      D4=ALOG(MEDIAN)
C
      DO 100 I=1,NBINS
          ND(I)=D2*EXP(D3*(PTSIZE(I,7)-D4)**2)
 100  CONTINUE
C----------------------------------------------------------------------------
      RETURN
      END
C
******************************************************************************
*     SUBROUTINE LGNPAR(NBINS,PTSIZE,ND,NTOTAL,MEDIAN,GSTDEV)
******************************************************************************
*
*     This subroutine fits a log-normal distribution to a given
*     particle size distribution ND, and calculates the total number 
*     density of particles NTOTAL, median radius MEDIAN, and geometric 
*     standard deveation GSTDEV. Subroutine SETBIN must have been called 
*     prior to the use of this subroutine.
*
*     Input/output variables:
*     INTEGER NBINS
*     REAL(KIND=4)    NTOTAL,MEDIAN,GSTDEV,ND(NBINS),PTSIZE(NBINS,8)
*
*     Input:
*         NBINS:                  Number of particle radii bins
*         PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume
*         ND:      (par/kg air)   Number density of particles
*
*     Output:
*         NTOTAL:  (par/kg air)   Total number density of particles
*         MEDIAN:  (m)            Median particle radius
*         GSTDEV:                 Geometric standard deveation
*
      SUBROUTINE LGNPAR(NBINS,PTSIZE,ND,NTOTAL,MEDIAN,GSTDEV)
C
      IMPLICIT NONE
      INTEGER NBINS
      REAL(KIND=4)    NTOTAL,MEDIAN,GSTDEV,PTSIZE(NBINS,8),ND(NBINS)
C
C
C     Common block variables:
      REAL(KIND=4) VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
C     Auxiliary local variables:
      INTEGER I
      REAL(KIND=4) ALFA,BETA
C----------------------------------------------------------------------------
      NTOTAL=0.0
      ALFA=0.0
      BETA=0.0
      DO 100 I=1,NBINS
          NTOTAL=NTOTAL+ND(I)
          ALFA=ALFA+ND(I)*PTSIZE(I,7)
          BETA=BETA+ND(I)*PTSIZE(I,8)
 100  CONTINUE
C
      IF(NTOTAL.LE.ZERO) THEN
          MEDIAN=0.0
          GSTDEV=0.0
      ELSE
          ALFA=ALFA/NTOTAL
          BETA=BETA/NTOTAL-ALFA*ALFA
          MEDIAN=EXP(ALFA)
          IF(BETA.GT.0.0) THEN
              GSTDEV=EXP(SQRT(BETA))
          ELSE
              GSTDEV=0.0
          ENDIF
      ENDIF
C----------------------------------------------------------------------------
      RETURN
      END
C
*******************************************************************************
*     SUBROUTINE DSTPAR(NBINS,PTSIZE,ND,NTOTAL,MEANRD,SURFCE)
******************************************************************************
*
*     This subroutine calculates the total number density of particles 
*     NTOTAL, mean radius MEANRD, and surface area density of a differential 
*     particle size distribution ND.
*     Subroutine SETBIN must have been called prior to the use of this 
*     subroutine.
*
*     Input/output variables:
*     INTEGER NBINS
*     REAL(KIND=4)    NTOTAL,MEANRD,SURFCE,ND(NBINS),PTSIZE(NBINS,8)
*
*     Input:
*         NBINS:                  Number of particle radii bins
*         PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume
*         ND:      (par/kg air)   Number density of particles
*
*     Output:
*         NTOTAL:  (par/kg air)   Total number density of particles
*         MEANRD:  (m)            Mean particle radius
*         SURFCE:  (m**2/kg air)  Particle surface density
*
      SUBROUTINE DSTPAR(NBINS,PTSIZE,ND,NTOTAL,MEANRD,SURFCE)
C
      IMPLICIT NONE
      INTEGER NBINS
      REAL(KIND=4)    NTOTAL,MEANRD,SURFCE,PTSIZE(NBINS,8),ND(NBINS)
C
C
C     Common block variables:
      REAL(KIND=4) VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
C     Auxiliary local variables:
      INTEGER I
      REAL(KIND=4) POWER,X,PI,C1
      PARAMETER(POWER=1.0/3.0)
      PARAMETER(PI=3.1415926536)
      PARAMETER(C1=4.0*PI/3.0)
C----------------------------------------------------------------------------
      NTOTAL=0.0
      SURFCE=0.0
      X=0.0
      DO 100 I=1,NBINS
          NTOTAL=NTOTAL+ND(I)
          SURFCE=SURFCE+PTSIZE(I,2)*ND(I)
          X=X+ND(I)*PTSIZE(I,3)
 100  CONTINUE
C
      IF(NTOTAL.LE.ZERO) THEN
          MEANRD=0.0
          NTOTAL=0.0
          SURFCE=0.0
      ELSE
          MEANRD=((X/NTOTAL)**POWER)/C1
      ENDIF
C----------------------------------------------------------------------------
      RETURN
      END
C
C
*****************************************************************************
*     SUBROUTINE SURFCE(NBINS,ND,PTSIZE,S)
******************************************************************************
*
*     This subroutine calculates the total surface area density of an
*     ensemble of particles of a given type, i.e. particle surface area
*     per kg of air.
*
*     Input/output variables:
*     INTEGER NBINS
*     REAL(KIND=4) TAIR,PAIR,ND(NBINS),PTSIZE(NBINS,8),S
*
*     Input:
*         NBINS:                  Number of particle radii groups
*         ND:      (par/kg air)   Number density of particles
*         PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume
*
*     Output:
*         S:       (m**2/kg air)  Particle surface density
*
      SUBROUTINE SURFCE(NBINS,ND,PTSIZE,S)
C
      IMPLICIT NONE
      INTEGER NBINS
      REAL(KIND=4) ND(NBINS),PTSIZE(NBINS,8),S
C
C     Auxiliary local variables:
      INTEGER I
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C----------------------------------------------------------------------------
      S=0.0
      DO 100 I=1,NBINS
          S=S+PTSIZE(I,2)*ND(I)
 100  CONTINUE
C
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING SURFCE','   S = ',S
C----------------------------------------------------------------------------
      RETURN
      END
C
*****************************************************************************
*     SUBROUTINE VOLUME(NBINS,ND,PTSIZE,V)
******************************************************************************
*
*     This subroutine calculates the total volume density of an
*     ensemble of particles of a given type, i.e. particle volume
*     per kg of air.
*
*     Input/output variables:
*     INTEGER NBINS
*     REAL(KIND=4) TAIR,PAIR,ND(NBINS),PTSIZE(NBINS,8),S
*
*     Input:
*         NBINS:                  Number of particle radii groups
*         ND:      (par/kg air)   Number density of particles
*         PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume
*
*     Output:
*         V:       (m**3/kg air)  Particle surface density
*
      SUBROUTINE VOLUME(NBINS,ND,PTSIZE,V)
C
      IMPLICIT NONE
      INTEGER NBINS
      REAL(KIND=4) ND(NBINS),PTSIZE(NBINS,8),V
C
C     Auxiliary local variables:
      INTEGER I
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C----------------------------------------------------------------------------
      V=0.0
      DO 100 I=1,NBINS
          V=V+PTSIZE(I,3)*ND(I)
 100  CONTINUE
C
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING SURFCE','   V = ',V
C----------------------------------------------------------------------------
      RETURN
      END
C
*******************************************************************************
*      SUBROUTINE MIXCON(NBINS,TYPES,ND,MCS,MCN,MCW,MRCWV,MRCNA,MRCSA)
******************************************************************************
*
*     This subroutine calculates the total condensed phase mass mixing ratio
*     of sulfuric acid, nitric acid, and water in the liquid and solid particles.
*
*     Declaration of input/output variables:
*
*     INTEGER NBINS,TYPES
*
*      REAL(KIND=4)
*             ND(NBINS,TYPES),
*             MCS(NBINS,TYPES),MCN(NBINS,TYPES),MCW(NBINS,TYPES),
*             MRCWV,MRCNA,MRCSA
*
*     Input:
*         NBINS:                  Number of particle radii bins
*         TYPES:                  Number of particle types comprehended by the model
*         ND:      (par/kg air)   Number density of particles
*         MCS:     (kg/particle)  Mass of condensed sulfuric acid per particle
*         MCN:     (kg/particle)  Mass of condensed nitric acid per particle
*         MCW:     (kg/particle)  Mass of condensed water per particle
*
*     Output:
*         MRCWV:   (kg/kg air)    Mass mixing ratio of condensed water
*         MRCNA:   (kg/kg air)    Mass mixing ratio of condensed nitric acid
*         MRCSA:   (kg/kg air)    Mass mixing ratio of condensed sulfuric acid
*
C
      SUBROUTINE MIXCON(
     +        NBINS,TYPES,
     +        ND,MCS,MCN,MCW,
     +        MRCWV,MRCNA,MRCSA)
C
      IMPLICIT NONE
      INTEGER NBINS,TYPES
C
      REAL(KIND=4)
     +        ND(NBINS,TYPES),
     +        MCS(NBINS,TYPES),MCN(NBINS,TYPES),MCW(NBINS,TYPES),
     +        MRCWV,MRCNA,MRCSA
C
C     Common block variables:
      REAL(KIND=4) VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C
C     Auxiliary local variables:
      INTEGER I,J
C----------------------------------------------------------------------------
      MRCSA=0.0
      MRCNA=0.0
      MRCWV=0.0
      DO J=1,TYPES
          DO I=1,NBINS
              MRCSA=MRCSA+MCS(I,J)*ND(I,J)
              MRCNA=MRCNA+MCN(I,J)*ND(I,J)
              MRCWV=MRCWV+MCW(I,J)*ND(I,J)
          END DO
      END DO
C
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING MIXCON   MRCWV = ',MRCWV,
     +              '  MRCNA = ',MRCNA,' MRCSA = ',MRCSA
C----------------------------------------------------------------------------
      RETURN
      END
C
*******************************************************************************
*      SUBROUTINE COMPOSITION(NBINS,ND,MCS,MCN,MCW,WSA,WNA)
******************************************************************************
*
*     This subroutine calculates the weight fraction
*     of sulfuric acid and nitric acid in a given particle type.
*
*     Declaration of input/output variables:
*
*     INTEGER NBINS,TYPES
*
*      REAL(KIND=4)
*             ND(NBINS),
*             MCS(NBINS),MCN(NBINS),MCW(NBINS),
*             WSA,WNA
*
*     Input:
*         NBINS:                  Number of particle radii bins
*         ND:      (par/kg air)   Number density of particles
*         MCS:     (kg/particle)  Mass of condensed sulfuric acid per particle
*         MCN:     (kg/particle)  Mass of condensed nitric acid per particle
*         MCW:     (kg/particle)  Mass of condensed water per particle
*
*     Output:
*         WSA:                    Weight fraction of sulfuric acid
*         WNA:                    Weight fraction of nitric acid
*
C
      SUBROUTINE COMPOSITION(
     +        NBINS,
     +        ND,MCS,MCN,MCW,
     +        WSA,WNA)
C
      IMPLICIT NONE
      INTEGER NBINS,TYPES
C
      REAL(KIND=4)
     +        ND(NBINS),
     +        MCS(NBINS),MCN(NBINS),MCW(NBINS),
     +        WSA,WNA,WWV
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C
C     Auxiliary local variables:
      INTEGER I
C----------------------------------------------------------------------------
      WSA=0.0
      WNA=0.0
      WWV=0.0
      DO I=1,NBINS
          WSA=WSA+MCS(I)*ND(I)
          WNA=WNA+MCN(I)*ND(I)
          WWV=WWV+MCW(I)*ND(I)
      END DO
      WWV=WSA+WNA+WWV
      IF(WWV.GT.0.0) THEN
          WSA=WSA/WWV
          WNA=WNA/WWV
      ELSE
          WSA=0.0
          WNA=0.0
      ENDIF
C----------------------------------------------------------------------------
      RETURN
      END
C
*******************************************************************************
*     SUBROUTINE PNDTOT(NBINS,ND,NDTOT)
******************************************************************************
*
*     This subroutine calculates the total number density of an
*     ensemble of particles of a given type, i.e. number of particles
*     per kg of air.
*
*     Declaration of input/output variables:
*
*     INTEGER NBINS
*
*     REAL(KIND=4)
*             NDPSC(NBINS),ND
*
*     Input:
*         NBINS:                  Number of particle radii bins
*         ND:      (par/kg air)   Number density of particles
*
*     Output:
*         NDTOT:   (par/kg air)   Total number density of particles
*
      SUBROUTINE PNDTOT(NBINS,ND,NDTOT)
C
      IMPLICIT NONE
      INTEGER NBINS
      REAL(KIND=4) ND(NBINS),NDTOT
C
C     Auxiliary local variables:
      INTEGER I
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C----------------------------------------------------------------------------
      NDTOT=0.0
      DO 100 I=1,NBINS
          NDTOT=NDTOT+ND(I)
 100  CONTINUE
C
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING PNDTOT   NDTOT = ',NDTOT
C----------------------------------------------------------------------------
      RETURN
      END
C
******************************************************************************
*     SUBROUTINE OPC_WYOMING(NCLASS,RADIUSCLASS)
******************************************************************************
*
*     Radii corresponding to U. Wyoming optical counter size classes:
*     Used for calculations of cumulated size distributions, i.e.
*     number concentrations of particles with radii greater than radiusclass(i)
*
*
*     Output:
*     In array RADIUSCLASS is stored limiting radii in cululated size distributions.
*
*
*
      SUBROUTINE OPC_WYOMING(NCLASS,RADIUSCLASS)
      IMPLICIT NONE
      INTEGER NCLASS
      REAL(KIND=4) RADIUSCLASS(NCLASS)
C----------------------------------------------------------------------------
c     Radii corresponding to U. Wyoming optical counter size classes:
C     Used for calculations of cumulated size distributions, i.e.
c     number concentrations of particles with radii greater than radiusclass(i)
      RADIUSCLASS(1)=0.0
      RADIUSCLASS(2)=0.15E-6
      RADIUSCLASS(3)=0.25E-6
      RADIUSCLASS(4)=0.30E-6
      RADIUSCLASS(5)=0.50E-6
      RADIUSCLASS(6)=0.75E-6
      RADIUSCLASS(7)=1.08E-6
      RADIUSCLASS(8)=1.25E-6
      RADIUSCLASS(9)=1.75E-6
      RADIUSCLASS(10)=2.50E-6
      RADIUSCLASS(11)=3.50E-6
      RADIUSCLASS(12)=5.00E-6
      RADIUSCLASS(13)=10.0E-6
C----------------------------------------------------------------------------
      RETURN
      END
C

******************************************************************************
*     SUBROUTINE HETNUC(VAPOR,NBINS,TAIR,PP,SAT,ROBULK,ST,PND,PTSIZE,NCRATE)
******************************************************************************
*
*     This subroutine calculates the heterogeneous nucleation rate.
*     Cf. Pruppacher & Klett (1980), p. 236.
*
*     The nucleation rate is only calculated for bins actually holding
*     particles to nucleate, i.e. where the particle number density is
*     non-zero, otherwise the nucleation rate is set to zero.
*
*     Each type of nucleating vapor is specified by a specific vapor number.
*     The vapor numbers are defined as:
*
*     1:  Water vapor
*     2:  Sulfuric acid vapor
*     3:  Nitric acid vapor
*
*     Declaration of input/output variables:
*
*     INTEGER VAPOR,NBINS
*     REAL(KIND=4) 
*         TAIR,PP,S,ROBULK,ST,PND(NBINS),PTSIZE(NBINS,3),NCRATE(NBINS)
*
*     Input:
*
*         VAPOR:                  Vapor number
*         NBINS:                  Number of particel radius groups
*         TAIR:    (K)            Ambient air temperature
*         PP:      (Pa)           Partial pressure of nucleating vapor
*         SAT:                    Saturation ration
*         ROBULK:  (kg/m**3)      Bulk density of condensing substance
*         ST:      (J/m**2)       Bulk surface tension of particles
*         PND:     (par/m**3)     Number density of particles to be nucleated
*         PTSIZE:  (m,m**2,m**3)  Particle radius, surface, volume
*
*     Output:
*
*         NCRATE:  (germs/par s)  Heterogeneous nucleation rate (pr particle)
*
      SUBROUTINE HETNUC(VAPOR,NBINS,TAIR,PP,SAT,ROBULK,ST,
     +                  PND,PTSIZE,NCRATE,ANYONE)
C
      IMPLICIT NONE
      INTEGER VAPOR,NBINS
      REAL(KIND=4) 
     +   TAIR,PP,SAT,ROBULK,ST,PND(NBINS),PTSIZE(NBINS,3),NCRATE(NBINS)
      LOGICAL ANYONE
      INTEGER I
C----------------------------------------------------------------------------
C     Common block variables:
      REAL(KIND=4)  VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C
C----------------------------------------------------------------------------
C     Computational constants:
      REAL(KIND=4)  XMIN
C     XMIX:   Smallest number which give ALOG(XMIN) > 0
      PARAMETER(XMIN=1.0E0+1.0E-7)
C----------------------------------------------------------------------------
C     Mathematical constants:
      DOUBLE PRECISION PI
      PARAMETER(PI=3.1415926536E0)
C
C     Physical constants:
      DOUBLE PRECISION MICE,MNAT,CM,GAMMA,MH2O,MHNO3,RGAS,KBOLTZ
      PARAMETER(
C       Compatibility parameter ice (Toon et al. 1989)
     +          MICE=0.95E0,
C       Compatibility parameter NAT (Wofsy et al. 1990)
CC     +          MNAT=0.95E0,
C       Compatibility parameter NAT (Iraci et al. 1995)
     +          MNAT=0.76E0,
C       Kinetic coefficient (1/m**2) (Toon et al. 1989)
CC     +          CM=3.0E+19,
     +          CM=2.0E+23,
C       Correction factor for surface tension of ice (Toon et al.1989):
     +          GAMMA=3.0E-10,
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
C       Molar weight of nitric acid (kg/mole)
     +          MHNO3=63.01E-3,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0,
C       Boltzmanns constant (J/K)
     +          KBOLTZ=1.380662E-23)
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      DOUBLE PRECISION
     +        C1,C2ICE,C2NAT,C3ICE,C3NAT,C4,C5,C6,
     +        B1,B2,B3,B4,B5,B6,B7,
     +        AG,X,FI,F,DELTAF,S,G
C
      PARAMETER(
     +        C1=4.0E0*PI/3.0E0,
     +        C2ICE=2.0E0*MICE,
     +        C2NAT=2.0E0*MNAT,
     +        C3ICE=3.0E0*MICE,
     +        C3NAT=3.0E0*MNAT,
     +        C4=2.0E0/RGAS,
     +        C5=CM/(4.0E0*PI*KBOLTZ),
     +        C6=2.0E0*KBOLTZ/(PI*RGAS))
C----------------------------------------------------------------------------
      S=AMAX1(SAT,XMIN)
      IF(VAPOR.EQ.1) THEN
          AG=C4*MH2O*ST/(TAIR*ROBULK*DLOG(S))
          B6=C5*SQRT(C6*MH2O*ST)*PP/(TAIR*ROBULK)
          B1=C1*ST*AG*AG
          B5=KBOLTZ*TAIR
          ANYONE=.FALSE.
          DO I=1,NBINS
              IF(PND(I).LE.ZERO) THEN
                   NCRATE(I)=0.0E0
              ELSE
                   ANYONE=.TRUE.
                   G=1.0E0+GAMMA/PTSIZE(I,1)
                   X=PTSIZE(I,1)/(AG*G)
                   FI=SQRT(1.0E0+X*(X-C2ICE))
                   B2=(X-MICE)/FI
                   B3=((1.0E0-MICE*X)/FI)**3
                   B4=X*X
                   F=0.5E0*(1.0E0+B3+B4*(X*(2.0E0+(B2*B2-3.0E0)*B2)+
     +                      C3ICE*(B2-1.0E0)))
                   DELTAF=B1*F*G**3
                   B7=B6*PTSIZE(I,2)*EXP(-DELTAF/B5)
                   NCRATE(I)=B7
              ENDIF
          ENDDO
      ELSE
          AG=C4*MHNO3*ST/(TAIR*ROBULK*DLOG(S))
          B6=C5*SQRT(C6*MHNO3*ST)*PP/(TAIR*ROBULK)
          B1=C1*ST*AG*AG
          B5=KBOLTZ*TAIR
          ANYONE=.FALSE.
          DO I=1,NBINS
              IF(PND(I).LE.ZERO) THEN
                   NCRATE(I)=0.0E0
              ELSE
                   ANYONE=.TRUE.
                   X=PTSIZE(I,1)/AG
                   FI=SQRT(1.0E0+X*(X-C2NAT))
                   B2=(X-MNAT)/FI
                   B3=((1.0E0-MNAT*X)/FI)**3
                   B4=X*X
                   F=0.5E0*(1.0E0+B3+B4*(X*(2.0E0+(B2*B2-3.0E0)*B2)+
     +                      C3NAT*(B2-1.0E0)))
                   DELTAF=B1*F
                   B7=B6*PTSIZE(I,2)*EXP(-DELTAF/B5)
                   NCRATE(I)=B7
              ENDIF
          ENDDO
      ENDIF
C----------------------------------------------------------------------------
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING HETNUC','   VAPOR: ',VAPOR
      WRITE(TEST,1010) 'TAIR','PP','S','ROBULK','ST','AG'
      WRITE(TEST,1000) TAIR,PP,S,ROBULK,ST,AG
      WRITE(TEST,1010) 'B1','B2','B3','B4','B5','B6'
      WRITE(TEST,1000) B1,B2,B3,B4,B5,B6
      WRITE(TEST,1010) 'AG','X','FI','F','DELTAF'
      WRITE(TEST,1000) AG,X,FI,F,DELTAF
      WRITE(TEST,*) 'ANYNOE ',ANYONE
      WRITE(TEST,*) '(NCRATE(I),I=1,NBINS)'
      WRITE(TEST,1000) (NCRATE(I),I=1,NBINS)
 1000 FORMAT(10(1PE9.2,1X),0P)
 1010 FORMAT(1X,8(A9))
C
C----------------------------------------------------------------------------
      RETURN
      END
C
******************************************************************************
*     SUBROUTINE CONDEN_NON_EQUIL(VAPOR,PHASE,PAR_TYPE,TAIR,PAIR,MFPAIR,KAIR,DWMAIR,
*                       S,PSAT,ROBULK,ST,LS,NBINS,RADIUS,VENTIL,MF)
******************************************************************************
*
*     This subroutine calculates the potential mass flux of vapor
*     due to condensation/evatoration to/from single particles of
*     the specified radii at specified ambient temperature,
*     ambient pressure, vapor partial pressures and saturation pressure.
*
*     The mass flux is only calculated, if the bin actually holds any
*     particles,specified by PAR_TYPE(I).GE.1 (i.e. bins holding no
*     particles are specified with PAR_TYPE(I).LT.0
*     Destinction is made between condensation/evaporation onto liquid or solid
*     phase, assuming a capacity parameter = 1 for liquid particles
*     and a capacity parameter = 1.61 for solid particles.
*     The input parameter PHASE specifies if condensation/evaporation the
*     mass flux is to be calculated for liquid (PHASE.EQ.1) or solid (PHASE.GT.1)
*     particles. The physical phase of the particles in each bin is specified
*     by the PAR_TYPE parameter as:
*     1:   Liquid phase
*     >1:  Solid phase
*     The mass flux is only calculated for bins where PHASE.EQ.PAR_TYPE(i)
*
*     The particles are assumed to have a radius dependent composition;
*     i.e. the saturation ratio S, density ROBULK, surface tension ST, and
*     ventilation factor depend on radius; otherwise USE subroutine CONDEN
*
*     Each type of condensing vapor is specified by a specific vapor number.
*     The vapor numbers are defined as:
*
*     1:  Water vapor
*     2:  Sulfuric acid vapor
*     3:  Nitric acid vapor
*
*
*     Input/output variables:
*     INTEGER VAPOR,NBINS,PAR_TYPE(NBINS),PHASE
*     REAL(KIND=4)  TAIR,PAIR,MFPAIR,KAIR,DWMAIR,
*          S,PSAT,ROBULK,ST,LS,ND(NBINS),RADIUS(NBINS),MF(NBINS)
*
*     Input:
*         VAPOR:                  Vapor number (1:H2O; 2:H2SO4; 3:HNO3)
*         PHASE:                  Flag indicating if condensation/evaporation
*                                 mass flux is to be calculated for liquid (1)
*                                 or solid (>1) particles
*         PAR_TYPE:                   Physical phase of particles in bin i
*                                   ( 1: liquid; >1: solid; <0: no particles)
*         TAIR:    (K)            Ambient air temperature
*         PAIR:    (Pa)           Ambient air pressure
*         MFPAIR:  (m)            Mean free path of air molecules
*         KAIR:    (J/m s K)      Thermal conductivity of air
*         DWMAIR:  (m**2/s)       Diffusivity of water molecules in air
*         S:                      Saturation ration
*         PSAT:    (Pa)           Saturation pressure of specified vapor
*                                 over plane surface
*         ROBULK:  (kg/m**3)      Bulk density of particles
*         ST:      (J/m**2)       Bulk surface tension of particles
*         LS:      (J/kg)         Heat of sublimation/evaporation
*                                 of specified vapor
*         NBINS:                  Number of particel radius groups
*         ND:     (par/m**3)      Number density of particles
*         RADIUS:  (m)            Particle radii
*         VENTIL:                 Condensation ventilation factors
*
*     Output:
*         MF:      (kg/s)         Potential mass flux of specified vapor
*                                 to single particle
*
      SUBROUTINE CONDEN_NON_EQUIL(VAPOR,PHASE,PAR_TYPE,TAIR,PAIR,
     +                  MFPAIR,KAIR,DWMAIR,
     +                  S,PSAT,ROBULK,ST,LS,NBINS,RADIUS,VENTIL,
     +                  MF)
C
      IMPLICIT NONE
      INTEGER VAPOR,NBINS,PAR_TYPE(NBINS),PHASE
      REAL(KIND=4)  TAIR,PAIR,MFPAIR,KAIR,DWMAIR,
     +        S(NBINS),PSAT(NBINS),ROBULK(NBINS),ST(NBINS),
     +        LS,RADIUS(NBINS),VENTIL(NBINS),
     +        MF(NBINS)
C----------------------------------------------------------------------------
C     Common block variables:
      REAL(KIND=4)  VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C----------------------------------------------------------------------------
C     Mathematical constants:
      REAL(KIND=4)  PI
      PARAMETER(PI=3.1415926536E0)
C
C     Physical constants:
      REAL(KIND=4)  MH2O,MH2SO4,MHNO3,MAIR,MYWV,MYSA,MYNA,THETAN,THETAS,
     +     RGAS,CPAIR,
     +     AH2O,AH2SO4,AHNO3,ALFAT,
     +     CAPALIQ,CAPASOL
      PARAMETER(
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
C       Molar weight of sulfuric acid (kg/mole)
     +          MH2SO4=98.08E-3,
C       Molar weight of nitric acid (kg/mole)
     +          MHNO3=63.01E-3,
C       Molar weight of dry air (kg/mole)
     +          MAIR=28.9644E-3,
C       Reduction factor of mean free path for water vapor molecules
     +          MYWV=0.820E0,
C       Reduction factor of mean free path for sulfuric acid vapor molecules
     +          MYSA=0.696E0,
C       Reduction factor of mean free path for nitric acid vapor molecules
     +          MYNA=0.857E0,
C       Reduction factor of diffusion of nitric acid vapor molecules
cc   +          THETAN=0.5355E0,
     +          THETAN=0.559E0,
C       Reduction factor of diffusion of sulfuric acid vapor molecules
     +          THETAS=0.364E0,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0,
C       Specific heat capacity of dry air (J/(mole K))
     +          CPAIR=7.0E0*RGAS/2.0E0,
C       Sticking coefficient H2O:
     +          AH2O=1.0E0,
C       Sticking coefficient H2SO4:
     +          AH2SO4=1.0E0,
C       Sticking coefficient HNO3:
     +          AHNO3=1.0E0,
cc   +          AHNO3=0.3E0,
C       Thermal accomodation coefficient:
     +          ALFAT=1.0E0,
C       Liquid particle "capacity" factor:
     +          CAPALIQ=1.0E0,
C       Solid particle "capacity" factor:
cc     +          CAPASOL=1.61E0)
     +          CAPASOL=1.11E0)
C
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      REAL(KIND=4) 
     +    C1WV,C1SA,C1NA,
     +    C2WV,C2SA,C2NA,
     +    C3,
     +    C4,
     +    C5,
     +    B1,B1LIQ,B1SOL,
     +    B2,
     +    B3,B3LIQ,B3SOL,
     +    MFP,
     +    B4,
     +    KND,KNT,
     +    A1,A2,A3,
     +    ALFAWV,ALFASA,ALFANA,ALFA,KN,ALFATT,CAPACI
      INTEGER I,J
C
      PARAMETER(
     +    C1WV=4.0E0*PI*MH2O/RGAS,
     +    C1SA=4.0E0*PI*MH2SO4/RGAS,
     +    C1NA=4.0E0*PI*MHNO3/RGAS,
     +    C2WV=2.0E0*MH2O/RGAS,
     +    C2SA=2.0E0*MH2SO4/RGAS,
     +    C2NA=2.0E0*MHNO3/RGAS,
     +    C3=4.0E0*PI*RGAS,
     +    C4=8.0E0*RGAS/(PI*MAIR),
     +    C5=3.0E0*RGAS/(CPAIR-0.5E0*RGAS),
     +    ALFAWV=1.33E0*(1.0E0-AH2O)/AH2O,
     +    ALFASA=1.33E0*(1.0E0-AH2SO4)/AH2SO4,
     +    ALFANA=1.33E0*(1.0E0-AHNO3)/AHNO3,
     +    ALFATT=1.33E0*(1.0E0-ALFAT)/ALFAT
     +        )
C----------------------------------------------------------------------------
C     Internal (statement) functions:
      REAL(KIND=4)  LAMBDA,KC
C
C     LAMBDA:      
C     KC:
C     Functions used to correct for gas-kinetic transport effects in
C     the continuum regime: 
C     Correction functions for gas-kinetic effect in continuum regime
C     are used, (Fuchs & Sutugin, 1971, Toon et al. 1989).
C     Corrections due to sticking coefficients and thermal accomodation
C     coefficients may be taken into consideration.
C     Note that some condensational and thermal accomodation coefficients
C     have been set equal to 1, and the correction factor for particle 
C     shape set equal to 1. 
C
      LAMBDA(KN,ALFA)=(1.33E0+0.71E0/KN)/(1.0E0+1.0E0/KN)+ALFA
      KC(KN,ALFA,CAPACI)=1.0E0/(1.0E0+LAMBDA(KN,ALFA)*CAPACI*KN)
C----------------------------------------------------------------------------
C     Thermodynamical calculations:
      B4=C5*KAIR*TAIR/(PAIR*SQRT(C4*TAIR))
      IF(VAPOR.EQ.1) THEN
C         Water vapor condensing:
          B1=C1WV*DWMAIR/TAIR
          B1LIQ=B1*CAPALIQ
          B1SOL=B1*CAPASOL
          B2=C2WV/TAIR
          B3=(LS*MH2O-RGAS*TAIR)*LS/(C3*TAIR*TAIR*KAIR)
          B3LIQ=B3/CAPALIQ
          B3SOL=B3/CAPASOL
          MFP=MYWV*MFPAIR
          ALFA=ALFAWV
      ELSE IF(VAPOR.EQ.2) THEN
C         Sulfuric acid condensing:
          B1=C1SA*DWMAIR*THETAS/TAIR
          B1LIQ=B1*CAPALIQ
          B1SOL=B1*CAPASOL
          B2=C2SA/TAIR
          B3=(LS*MH2SO4-RGAS*TAIR)*LS/(C3*TAIR*TAIR*KAIR)
          B3LIQ=B3/CAPALIQ
          B3SOL=B3/CAPASOL
          MFP=MYSA*MFPAIR
          ALFA=ALFASA
      ELSE IF(VAPOR.EQ.3) THEN
C         Nitric acid condensing:
          B1=C1NA*DWMAIR*THETAN/TAIR
          B1LIQ=B1*CAPALIQ
          B1SOL=B1*CAPASOL
          B2=C2NA/TAIR
          B3=(LS*MHNO3-RGAS*TAIR)*LS/(C3*TAIR*TAIR*KAIR)
          B3LIQ=B3/CAPALIQ
          B3SOL=B3/CAPASOL
          MFP=MYNA*MFPAIR
          ALFA=ALFANA
      ENDIF
C----------------------------------------------------------------------------
      DO I=1,NBINS,1
          IF(PAR_TYPE(I).LT.0) THEN
             MF(I)=0.0E0
          ELSE IF((PHASE.EQ.1.AND.PAR_TYPE(I).EQ.1).OR.
     +            (PHASE.GT.1.AND.PAR_TYPE(I).GT.1)) THEN
              IF(PAR_TYPE(I).EQ.1) THEN
                  CAPACI=CAPALIQ
                  B1=B1LIQ
                  B3=B3LIQ
              ELSE
                 CAPACI=CAPASOL
                 B1=B1SOL
                 B3=B3SOL
              ENDIF
              KND=MFP/RADIUS(I)
              KNT=B4/RADIUS(I)
              A1=B1*PSAT(I)*RADIUS(I)*VENTIL(I)*KC(KND,ALFA,CAPACI)
              A2=EXP(B2*ST(I)/(RADIUS(I)*ROBULK(I)))
              A3=B3/(KC(KNT,ALFATT,CAPACI)*RADIUS(I)*VENTIL(I))
              MF(I)=A1*(S(I)-A2)/(1.0E0+A1*A2*A3)
          ELSE
              MF(I)=0.0E0
          ENDIF
      ENDDO
C
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING CONDEN_NON_EQUIL','   VAPOR: ',VAPOR
      WRITE(TEST,1010) 'TAIR','PAIR'
      WRITE(TEST,1000) TAIR,PAIR
      WRITE(TEST,1010) 'B1','B2','B3','B4','MFP'
      WRITE(TEST,1000) B1,B2,B3,B4,MFP
CC    WRITE(TEST,*) '(VENTIL(I),I=1,NBINS)'
CC    DO I=0,NBINS/10-1
CC    WRITE(TEST,1000) (VENTIL(I*10+J),J=1,10)
CC    ENDDO
CC    WRITE(TEST,*) '(ROBULK(I),I=1,NBINS)'
CC    DO I=0,NBINS/10-1
CC    WRITE(TEST,1000) (ROBULK(I*10+J),J=1,10)
CC    ENDDO
CC    WRITE(TEST,*) '(ST(I),I=1,NBINS)'
CC    DO I=0,NBINS/10-1
CC    WRITE(TEST,1000) (ST(I*10+J),J=1,10)
CC    ENDDO
      WRITE(TEST,*) '(PSAT(I),I=1,NBINS)'
      DO I=0,NBINS/10-1
      WRITE(TEST,1000) (PSAT(I*10+J),J=1,10)
      ENDDO
      WRITE(TEST,*) '(S(I),I=1,NBINS)'
      DO I=0,NBINS/10-1
      WRITE(TEST,1000) (S(I*10+J),J=1,10)
      ENDDO
      WRITE(TEST,*) '(MF(I),I=1,NBINS)'
      DO I=0,NBINS/10-1
      WRITE(TEST,1000) (MF(I*10+J),J=1,10)
      ENDDO
 1000 FORMAT(10(1PE9.2,1X),0P)
 1010 FORMAT(1X,10(A9))
C----------------------------------------------------------------------------
      RETURN
      END
C
C
******************************************************************************
*     SUBROUTINE CONDEN_NON_EQUIL_COMPO(VAPOR,PHASE,TAIR,PAIR,MFPAIR,KAIR,DWMAIR,
*                       S,PSAT,ROBULK,ST,LS,NBINS,RADIUS,VENTIL,MF)
******************************************************************************
*
*     This subroutine calculates the potential mass flux of vapor
*     due to condensation/evatoration to/from single particles of
*     the specified radii at specified ambient temperature,
*     ambient pressure, vapor partial pressures and saturation pressure.
*
*     The mass flux is only calculated, if the bin actually holds any
*     particles,specified by PAR_TYPE(I).GE.1 (i.e. bins holding no
*     particles are specified with PAR_TYPE(I).LT.0
*     Destinction is made between condensation/evaporation onto liquid or solid
*     phase, assuming a capacity parameter = 1 for liquid particles
*     and a capacity parameter = 1.61 for solid particles.
*     The input parameter PHASE specifies if condensation/evaporation the
*     mass flux is to be calculated for liquid (PHASE.EQ.1) or solid (PHASE.GT.1)
*     particles. The physical phase of the particles in each bin is specified
*     by the PAR_TYPE parameter as:
*     1:   Liquid phase
*     >1:  Solid phase
*     The mass flux is only calculated for bins where PHASE.EQ.PAR_TYPE(i)
*
*     The particles are assumed to have a radius dependent composition;
*     i.e. the saturation ratio S, density ROBULK, surface tension ST, and
*     ventilation factor depend on radius; otherwise USE subroutine CONDEN
*
*     Each type of condensing vapor is specified by a specific vapor number.
*     The vapor numbers are defined as:
*
*     1:  Water vapor
*     2:  Sulfuric acid vapor
*     3:  Nitric acid vapor
*
*
*     Input/output variables:
*     INTEGER VAPOR,NBINS,PAR_TYPE(NBINS),PHASE
*     REAL(KIND=4)  TAIR,PAIR,MFPAIR,KAIR,DWMAIR,
*          S,PSAT,ROBULK,ST,LS,ND(NBINS),RADIUS(NBINS),MF(NBINS)
*
*     Input:
*         VAPOR:                  Vapor number (1:H2O; 2:H2SO4; 3:HNO3)
*         PHASE:                  Flag indicating if condensation/evaporation
*                                 mass flux is to be calculated for liquid (1)
*                                 or solid (>1) particles
*         PAR_TYPE:                   Physical phase of particles in bin i
*                                   ( 1: liquid; >1: solid; <0: no particles)
*         TAIR:    (K)            Ambient air temperature
*         PAIR:    (Pa)           Ambient air pressure
*         MFPAIR:  (m)            Mean free path of air molecules
*         KAIR:    (J/m s K)      Thermal conductivity of air
*         DWMAIR:  (m**2/s)       Diffusivity of water molecules in air
*         S:                      Saturation ration
*         PSAT:    (Pa)           Saturation pressure of specified vapor
*                                 over plane surface
*         ROBULK:  (kg/m**3)      Bulk density of particles
*         ST:      (J/m**2)       Bulk surface tension of particles
*         LS:      (J/kg)         Heat of sublimation/evaporation
*                                 of specified vapor
*         NBINS:                  Number of particel radius groups
*         ND:     (par/m**3)      Number density of particles
*         RADIUS:  (m)            Particle radii
*         VENTIL:                 Condensation ventilation factors
*
*     Output:
*         MF:      (kg/s)         Potential mass flux of specified vapor
*                                 to single particle
*
      SUBROUTINE CONDEN_NON_EQUIL_COMPO(VAPOR,PHASE,TAIR,PAIR,
     +                  MFPAIR,KAIR,DWMAIR,
     +                  ND,S,PSAT,ROBULK,ST,LS,NBINS,RADIUS,VENTIL,
     +                  MF,LIMIT)
C
      IMPLICIT NONE
      INTEGER VAPOR,NBINS,PHASE,LIMIT
      LOGICAL ANYONE
      REAL(KIND=4)  TAIR,PAIR,MFPAIR,KAIR,DWMAIR,ND(NBINS),
     +        S(NBINS),PSAT(NBINS),ROBULK(NBINS),ST(NBINS),
     +        LS,RADIUS(NBINS),VENTIL(NBINS),
     +        MF(NBINS)
C----------------------------------------------------------------------------
C     Common block variables:
      REAL(KIND=4)  VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C----------------------------------------------------------------------------
C     Mathematical constants:
      REAL(KIND=4)  PI
      PARAMETER(PI=3.1415926536E0)
C
C     Physical constants:
      REAL(KIND=4)  MH2O,MH2SO4,MHNO3,MAIR,MYWV,MYSA,MYNA,THETAN,THETAS,
     +     RGAS,CPAIR,
     +     AH2O,AH2SO4,AHNO3,ALFAT,
     +     CAPALIQ,CAPASOL
      PARAMETER(
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
C       Molar weight of sulfuric acid (kg/mole)
     +          MH2SO4=98.08E-3,
C       Molar weight of nitric acid (kg/mole)
     +          MHNO3=63.01E-3,
C       Molar weight of dry air (kg/mole)
     +          MAIR=28.9644E-3,
C       Reduction factor of mean free path for water vapor molecules
     +          MYWV=0.820E0,
C       Reduction factor of mean free path for sulfuric acid vapor molecules
     +          MYSA=0.696E0,
C       Reduction factor of mean free path for nitric acid vapor molecules
     +          MYNA=0.857E0,
C       Reduction factor of diffusion of nitric acid vapor molecules
cc   +          THETAN=0.5355E0,
     +          THETAN=0.559E0,
C       Reduction factor of diffusion of sulfuric acid vapor molecules
     +          THETAS=0.364E0,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0,
C       Specific heat capacity of dry air (J/(mole K))
     +          CPAIR=7.0E0*RGAS/2.0E0,
C       Sticking coefficient H2O:
     +          AH2O=1.0E0,
C       Sticking coefficient H2SO4:
     +          AH2SO4=1.0E0,
C       Sticking coefficient HNO3:
     +          AHNO3=1.0E0,
cc   +          AHNO3=0.3E0,
C       Thermal accomodation coefficient:
     +          ALFAT=1.0E0,
C       Liquid particle "capacity" factor:
     +          CAPALIQ=1.0E0,
C       Solid particle "capacity" factor:
cc     +          CAPASOL=1.61E0)
     +          CAPASOL=1.11E0)
C
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      REAL(KIND=4) 
     +    C1WV,C1SA,C1NA,
     +    C2WV,C2SA,C2NA,
     +    C3,
     +    C4,
     +    C5,
     +    B1,B1LIQ,B1SOL,
     +    B2,
     +    B3,B3LIQ,B3SOL,
     +    MFP,
     +    B4,
     +    KND,KNT,
     +    A1,A2,A3,
     +    ALFAWV,ALFASA,ALFANA,ALFA,KN,ALFATT,CAPACI
      INTEGER I,J
C
      PARAMETER(
     +    C1WV=4.0E0*PI*MH2O/RGAS,
     +    C1SA=4.0E0*PI*MH2SO4/RGAS,
     +    C1NA=4.0E0*PI*MHNO3/RGAS,
     +    C2WV=2.0E0*MH2O/RGAS,
     +    C2SA=2.0E0*MH2SO4/RGAS,
     +    C2NA=2.0E0*MHNO3/RGAS,
     +    C3=4.0E0*PI*RGAS,
     +    C4=8.0E0*RGAS/(PI*MAIR),
     +    C5=3.0E0*RGAS/(CPAIR-0.5E0*RGAS),
     +    ALFAWV=1.33E0*(1.0E0-AH2O)/AH2O,
     +    ALFASA=1.33E0*(1.0E0-AH2SO4)/AH2SO4,
     +    ALFANA=1.33E0*(1.0E0-AHNO3)/AHNO3,
     +    ALFATT=1.33E0*(1.0E0-ALFAT)/ALFAT
     +        )
C----------------------------------------------------------------------------
C     Internal (statement) functions:
      REAL(KIND=4)  LAMBDA,KC
C
C     LAMBDA:      
C     KC:
C     Functions used to correct for gas-kinetic transport effects in
C     the continuum regime: 
C     Correction functions for gas-kinetic effect in continuum regime
C     are used, (Fuchs & Sutugin, 1971, Toon et al. 1989).
C     Corrections due to sticking coefficients and thermal accomodation
C     coefficients may be taken into consideration.
C     Note that some condensational and thermal accomodation coefficients
C     have been set equal to 1, and the correction factor for particle 
C     shape set equal to 1. 
C
      LAMBDA(KN,ALFA)=(1.33E0+0.71E0/KN)/(1.0E0+1.0E0/KN)+ALFA
      KC(KN,ALFA,CAPACI)=1.0E0/(1.0E0+LAMBDA(KN,ALFA)*CAPACI*KN)
C----------------------------------------------------------------------------
C     Thermodynamical calculations:
      B4=C5*KAIR*TAIR/(PAIR*SQRT(C4*TAIR))
      IF(VAPOR.EQ.1) THEN
C         Water vapor condensing:
          B1=C1WV*DWMAIR/TAIR
          B1LIQ=B1*CAPALIQ
          B1SOL=B1*CAPASOL
          B2=C2WV/TAIR
          B3=(LS*MH2O-RGAS*TAIR)*LS/(C3*TAIR*TAIR*KAIR)
          B3LIQ=B3/CAPALIQ
          B3SOL=B3/CAPASOL
          MFP=MYWV*MFPAIR
          ALFA=ALFAWV
      ELSE IF(VAPOR.EQ.2) THEN
C         Sulfuric acid condensing:
          B1=C1SA*DWMAIR*THETAS/TAIR
          B1LIQ=B1*CAPALIQ
          B1SOL=B1*CAPASOL
          B2=C2SA/TAIR
          B3=(LS*MH2SO4-RGAS*TAIR)*LS/(C3*TAIR*TAIR*KAIR)
          B3LIQ=B3/CAPALIQ
          B3SOL=B3/CAPASOL
          MFP=MYSA*MFPAIR
          ALFA=ALFASA
      ELSE IF(VAPOR.EQ.3) THEN
C         Nitric acid condensing:
          B1=C1NA*DWMAIR*THETAN/TAIR
          B1LIQ=B1*CAPALIQ
          B1SOL=B1*CAPASOL
          B2=C2NA/TAIR
          B3=(LS*MHNO3-RGAS*TAIR)*LS/(C3*TAIR*TAIR*KAIR)
          B3LIQ=B3/CAPALIQ
          B3SOL=B3/CAPASOL
          MFP=MYNA*MFPAIR
          ALFA=ALFANA
      ENDIF
      IF(PHASE.EQ.1) THEN
          CAPACI=CAPALIQ
          B1=B1LIQ
          B3=B3LIQ
      ELSE
         CAPACI=CAPASOL
         B1=B1SOL
         B3=B3SOL
      ENDIF
C----------------------------------------------------------------------------
      LIMIT=0
      ANYONE=.FALSE.
      DO I=1,NBINS,1
          IF(ND(I).LE.ZERO) THEN
             MF(I)=0.0E0
          ELSE
              ANYONE=.TRUE.
              KND=MFP/RADIUS(I)
              KNT=B4/RADIUS(I)
              A1=B1*PSAT(I)*RADIUS(I)*VENTIL(I)*KC(KND,ALFA,CAPACI)
              A2=EXP(B2*ST(I)/(RADIUS(I)*ROBULK(I)))
              A3=B3/(KC(KNT,ALFATT,CAPACI)*RADIUS(I)*VENTIL(I))
              MF(I)=A1*(S(I)-A2)/(1.0E0+A1*A2*A3)
          ENDIF
          IF(MF(I).LT.0.0) LIMIT=I
      ENDDO
      IF(.NOT.ANYONE) LIMIT=-1
C
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING CONDEN_NON_EQUIL','   VAPOR: ',VAPOR
      WRITE(TEST,1010) 'TAIR','PAIR'
      WRITE(TEST,1000) TAIR,PAIR
      WRITE(TEST,1010) 'B1','B2','B3','B4','MFP'
      WRITE(TEST,1000) B1,B2,B3,B4,MFP
CC    WRITE(TEST,*) '(VENTIL(I),I=1,NBINS)'
CC    DO I=0,NBINS/10-1
CC    WRITE(TEST,1000) (VENTIL(I*10+J),J=1,10)
CC    ENDDO
CC    WRITE(TEST,*) '(ROBULK(I),I=1,NBINS)'
CC    DO I=0,NBINS/10-1
CC    WRITE(TEST,1000) (ROBULK(I*10+J),J=1,10)
CC    ENDDO
CC    WRITE(TEST,*) '(ST(I),I=1,NBINS)'
CC    DO I=0,NBINS/10-1
CC    WRITE(TEST,1000) (ST(I*10+J),J=1,10)
CC    ENDDO
      WRITE(TEST,*) '(PSAT(I),I=1,NBINS)'
      DO I=0,NBINS/10-1
      WRITE(TEST,1000) (PSAT(I*10+J),J=1,10)
      ENDDO
      WRITE(TEST,*) '(S(I),I=1,NBINS)'
      DO I=0,NBINS/10-1
      WRITE(TEST,1000) (S(I*10+J),J=1,10)
      ENDDO
      WRITE(TEST,*) '(MF(I),I=1,NBINS)'
      DO I=0,NBINS/10-1
      WRITE(TEST,1000) (MF(I*10+J),J=1,10)
      ENDDO
 1000 FORMAT(10(1PE9.2,1X),0P)
 1010 FORMAT(1X,10(A9))
C----------------------------------------------------------------------------
      RETURN
      END
C
******************************************************************************
*     SUBROUTINE CONDEN(VAPOR,PHASE,TAIR,PAIR,MFPAIR,KAIR,DWMAIR,
*                       S,PSAT,ROBULK,ST,LS,NBINS,ND,RADIUS,VENTIL,MF,L)
******************************************************************************
*
*     This subroutine calculates the potential mass flux of vapor
*     due to condensation/evatoration to/from single solid particles of
*     the specified radii at specified ambient temperature,
*     ambient pressure, vapor partial pressures and saturation pressure.
*     The mass flux is only calculated, if the bin actually holds any
*     particles; otherwise the output is set to zero.
*     Destinction is made between condensation/evaporation onto liquid or solid
*     phase, assuming a capacity parameter = 1 for liquid particles
*     and a capacity parameter = 1.61 for solid particles.
*     The input parameter PHASE specifies if condensation/evaporation the
*     mass flux is to be calculated for liquid (PHASE.EQ.1) or solid (PHASE.GT.1)
*     particles. The physical phase of the particles in each bin is specified
*     by the PAR_TYPE parameter as:
*     1:   Liquid phase
*     >1:  Solid phase
*
*     Each type of condensing vapor is specified by a specific vapor number.
*     The vapor numbers are defined as:
*
*     1:  Water vapor
*     2:  Sulfuric acid vapor
*     3:  Nitric acid vapor
*
*
*     Input/output variables:
*     INTEGER VAPOR,NBINS
*     REAL(KIND=4)  TAIR,PAIR,MFPAIR,KAIR,DWMAIR,
*          S,PSAT,ROBULK,ST,LS,ND(NBINS),RADIUS(NBINS),MF(NBINS)
*
*     Input:
*         VAPOR:                  Vapor number
*         PHASE:                  Flag indicating if condensation/evaporation
*                                 mass flux is to be calculated for liquid (1)
*                                 or solid (>1) particles
*         TAIR:    (K)            Ambient air temperature
*         PAIR:    (Pa)           Ambient air pressure
*         MFPAIR:  (m)            Mean free path of air molecules
*         KAIR:    (J/m s K)      Thermal conductivity of air
*         DWMAIR:  (m**2/s)       Diffusivity of water molecules in air
*         S:                      Saturation ration
*         PSAT:    (Pa)           Saturation pressure of specified vapor
*                                 over plane surface
*         ROBULK:  (kg/m**3)      Bulk density of particles
*         ST:      (J/m**2)       Bulk surface tension of particles
*         LS:      (J/kg)         Heat of sublimation/evaporation
*                                 of specified vapor
*         NBINS:                  Number of particel radius groups
*         ND:     (par/m**3)      Number density of particles
*         RADIUS:  (m)            Particle radii
*         VENTIL:                 Condensation ventilation factors
*
*     Output:
*         MF:      (kg/s)         Potential mass flux of specified vapor
*                                 to single particle
*         L:                      Highest bin number of evaporating part.
*
      SUBROUTINE CONDEN(VAPOR,PHASE,TAIR,PAIR,MFPAIR,KAIR,DWMAIR,
     +                  S,PSAT,ROBULK,ST,LS,NBINS,ND,RADIUS,VENTIL,
     +                  MF,LIMIT)
C
      IMPLICIT NONE
      INTEGER VAPOR,NBINS,LIMIT,PHASE
      REAL(KIND=4)  TAIR,PAIR,MFPAIR,KAIR,DWMAIR,
     +        S,PSAT,ROBULK,ST,LS,ND(NBINS),RADIUS(NBINS),VENTIL(NBINS),
     +        MF(NBINS)
C----------------------------------------------------------------------------
C     Common block variables:
      REAL(KIND=4)  VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C----------------------------------------------------------------------------
C     Mathematical constants:
      REAL(KIND=4)  PI
      PARAMETER(PI=3.1415926536E0)
C
C     Physical constants:
      REAL(KIND=4)  MH2O,MH2SO4,MHNO3,MAIR,MYWV,MYSA,MYNA,THETAN,THETAS,
     +     RGAS,CPAIR,
     +     AH2O,AH2SO4,AHNO3,ALFAT,
     +     CAPALIQ,CAPASOL
      PARAMETER(
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
C       Molar weight of sulfuric acid (kg/mole)
     +          MH2SO4=98.08E-3,
C       Molar weight of nitric acid (kg/mole)
     +          MHNO3=63.01E-3,
C       Molar weight of dry air (kg/mole)
     +          MAIR=28.9644E-3,
C       Reduction factor of mean free path for water vapor molecules
     +          MYWV=0.820E0,
C       Reduction factor of mean free path for sulfuric acid vapor molecules
     +          MYSA=0.696E0,
C       Reduction factor of mean free path for nitric acid vapor molecules
     +          MYNA=0.857E0,
C       Reduction factor of diffusion of nitric acid vapor molecules
cc     +          THETAN=0.5355E0,
     +          THETAN=0.559E0,
C       Reduction factor of diffusion of sulfuric acid vapor molecules
     +          THETAS=0.364E0,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0,
C       Specific heat capacity of dry air (J/(mole K))
     +          CPAIR=7.0E0*RGAS/2.0E0,
C       Sticking coefficient H2O:
     +          AH2O=1.0E0,
C       Sticking coefficient H2SO4:
     +          AH2SO4=1.0E0,
C       Sticking coefficient HNO3:
     +          AHNO3=1.0E0,
cc   +          AHNO3=0.3E0,
C       Thermal accomodation coefficient:
     +          ALFAT=1.0E0,
C       Liquid particle "capacity" factor:
     +          CAPALIQ=1.0E0,
C       Solid particle "capacity" factor:
cc     +          CAPASOL=1.61E0)
     +          CAPASOL=1.11E0)
C
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      REAL(KIND=4) 
     +    C1WV,C1SA,C1NA,
     +    C2WV,C2SA,C2NA,
     +    C3,
     +    C4,
     +    C5,
     +    B1,B1LIQ,B1SOL,
     +    B2,
     +    B3,B3LIQ,B3SOL,
     +    B4,
     +    MFP,
     +    KND,KNT,
     +    A1,A2,A3,
     +    ALFAWV,ALFASA,ALFANA,ALFA,KN,ALFATT,CAPACI
      INTEGER I
      LOGICAL ANYONE
C
      PARAMETER(
     +    C1WV=4.0E0*PI*MH2O/RGAS,
     +    C1SA=4.0E0*PI*MH2SO4/RGAS,
     +    C1NA=4.0E0*PI*MHNO3/RGAS,
     +    C2WV=2.0E0*MH2O/RGAS,
     +    C2SA=2.0E0*MH2SO4/RGAS,
     +    C2NA=2.0E0*MHNO3/RGAS,
     +    C3=4.0E0*PI*RGAS,
     +    C4=8.0E0*RGAS/(PI*MAIR),
     +    C5=3.0E0*RGAS/(CPAIR-0.5E0*RGAS),
     +    ALFAWV=1.33E0*(1.0E0-AH2O)/AH2O,
     +    ALFASA=1.33E0*(1.0E0-AH2SO4)/AH2SO4,
     +    ALFANA=1.33E0*(1.0E0-AHNO3)/AHNO3,
     +    ALFATT=1.33E0*(1.0E0-ALFAT)/ALFAT
     +        )
C----------------------------------------------------------------------------
C     Internal (statement) functions:
      REAL(KIND=4)  LAMBDA,KC
C
C     LAMBDA:      
C     KC:
C     Functions used to correct for gas-kinetic transport effects in
C     the continuum regime: 
C     Correction functions for gas-kinetic effect in continuum regime
C     are used, (Fuchs & Sutugin, 1971, Toon et al. 1989).
C     Corrections due to sticking coefficients and thermal accomodation
C     coefficients may be taken into consideration.
C     Note that some condensational and thermal accomodation coefficients
C     have been set equal to 1, and the correction factor for particle 
C     shape set equal to 1. 
C
      LAMBDA(KN,ALFA)=(1.33E0+0.71E0/KN)/(1.0E0+1.0E0/KN)+ALFA
      KC(KN,ALFA,CAPACI)=1.0E0/(1.0E0+LAMBDA(KN,ALFA)*CAPACI*KN)
C----------------------------------------------------------------------------
C     Thermodynamical calculations:
      B4=C5*KAIR*TAIR/(PAIR*SQRT(C4*TAIR))
      IF(VAPOR.EQ.1) THEN
C         Water vapor condensing:
          B1=C1WV*DWMAIR*PSAT/TAIR
          B1LIQ=B1*CAPALIQ
          B1SOL=B1*CAPASOL
          B2=C2WV*ST/(TAIR*ROBULK)
          B3=(LS*MH2O-RGAS*TAIR)*LS/(C3*TAIR*TAIR*KAIR)
          B3LIQ=B3/CAPALIQ
          B3SOL=B3/CAPASOL
          MFP=MYWV*MFPAIR
          ALFA=ALFAWV
      ELSE IF(VAPOR.EQ.2) THEN
C         Sulfuric acid condensing:
          B1=C1SA*DWMAIR*THETAS*PSAT/TAIR
          B1LIQ=B1*CAPALIQ
          B1SOL=B1*CAPASOL
          B2=C2SA*ST/(TAIR*ROBULK)
          B3=(LS*MH2SO4-RGAS*TAIR)*LS/(C3*TAIR*TAIR*KAIR)
          B3LIQ=B3/CAPALIQ
          B3SOL=B3/CAPASOL
          MFP=MYSA*MFPAIR
          ALFA=ALFASA
      ELSE IF(VAPOR.EQ.3) THEN
C         Nitric acid condensing:
          B1=C1NA*DWMAIR*THETAN*PSAT/TAIR
          B1LIQ=B1*CAPALIQ
          B1SOL=B1*CAPASOL
          B2=C2NA*ST/(TAIR*ROBULK)
          B3=(LS*MHNO3-RGAS*TAIR)*LS/(C3*TAIR*TAIR*KAIR)
          B3LIQ=B3/CAPALIQ
          B3SOL=B3/CAPASOL
          MFP=MYNA*MFPAIR
          ALFA=ALFANA
      ENDIF
      IF(PHASE.EQ.1) THEN
          CAPACI=CAPALIQ
          B1=B1LIQ
          B3=B3LIQ
      ELSE
         CAPACI=CAPASOL
         B1=B1SOL
         B3=B3SOL
      ENDIF
C----------------------------------------------------------------------------
      LIMIT=0
      ANYONE=.FALSE.
      DO I=1,NBINS,1
          IF(ND(I).LE.ZERO) THEN
             MF(I)=0.0E0
          ELSE
              ANYONE=.TRUE.
              KND=MFP/RADIUS(I)
              KNT=B4/RADIUS(I)
              A1=B1*RADIUS(I)*VENTIL(I)*KC(KND,ALFA,CAPACI)
              A2=EXP(B2/RADIUS(I))
              A3=B3/(KC(KNT,ALFATT,CAPACI)*RADIUS(I)*VENTIL(I))
              MF(I)=A1*(S-A2)/(1.0E0+A1*A2*A3)
              IF(MF(I).LT.0.0) LIMIT=I
          ENDIF
      ENDDO
      IF(.NOT.ANYONE) LIMIT=-1
C
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING CONDEN','   VAPOR: ',VAPOR,' PHASE: ',PHASE
      WRITE(TEST,1010) 'TAIR','PAIR','S','PSAT','ROBULK','ST','LS',
     +                 'LIMIT'
      WRITE(TEST,1000) TAIR,PAIR,S,PSAT,ROBULK,ST,LS,DBLE(LIMIT)
      WRITE(TEST,1010) 'B1','B2','B3','MFP'
      WRITE(TEST,1000) B1,B2,B3,MFP
      WRITE(TEST,1010) 'B4','KAIR','KND','KNT','A1','A2','A3'
      WRITE(TEST,1000) B4,KAIR,KND,KNT,A1,A2,A3
      WRITE(TEST,*) '(MF(I),I=1,NBINS)'
      WRITE(TEST,1000) (MF(I),I=1,NBINS)
 1000 FORMAT(10(1PE9.2,1X),0P)
 1010 FORMAT(1X,8(A9))
C----------------------------------------------------------------------------
      RETURN
      END
C
******************************************************************************
*     SUBROUTINE TMVELO(NBINS,PHASE,TAIR,PAIR,RHOAIR,MFPAIR,DYNVIS,
*                       ROBULK,PND,PTSIZE,VELOCI,VENTIL)
******************************************************************************
*
*     This subroutine calculates the terminal fall velocity of
*     particles in air, cf. Pruppacher & Klett (1980), p. 323, 337
*     Fuchs (1964), par 7,8,12
*
*     The velocity is only calculated for bins actually holding
*     particles to fall, i.e. where the particle number density is
*     non-zero, otherwise the velocity is set to zero.
*     Particle radii < 500 microns. All particles are assumed to possess
*     a common average density ROBULK.
*
*     Declaration of input/output variables:
*
*     INTEGER NBINS,PHASE
*     REAL(KIND=4) 
*         TAIR,PAIR,RHOAIR,MFPAIR,DYNVIS,
*         ROBULK,PND(NBINS),PTSIZE(NBINS,3),VELOCI(NBINS),VENTIL(NBINS)
*
*     Input:
*
*         NBINS:                  Number of particel radius groups
*         PHASE:                  Physical phase (1: liquid; >1: solid)
*         TAIR:    (K)            Ambient air temperature
*         PAIR:    (Pa)           Ambient air pressure
*         RHOAIR:  (kg/m**3)      Ambient air density
*         MFPAIR:  (m)            Mean free path of air molecules
*         DYNVIS:  (kg/m s)       Dynamic viscosity of air
*         ROBULK:  (kg/m**3)      Bulk density of particles
*         PND:     (par/m**3)     Number density of particles
*         PTSIZE:  (m,m**2,m**3)  Particle radius, surface, volume
*
*     Output:
*
*         VELOCI:  (m/s)          Terminal fall velocity
*         VENTIL:                 Condensation ventilation factors
*
      SUBROUTINE TMVELO(NBINS,PHASE,TAIR,PAIR,RHOAIR,MFPAIR,DYNVIS,
     +                  ROBULK,PND,PTSIZE,VELOCI,VENTIL)
C
      IMPLICIT NONE
      INTEGER NBINS,PHASE
      REAL(KIND=4) 
     +        TAIR,PAIR,RHOAIR,MFPAIR,DYNVIS,
     +        ROBULK,PND(NBINS),PTSIZE(NBINS,3),
     +        VELOCI(NBINS),VENTIL(NBINS)
C----------------------------------------------------------------------------
C     Common block variables:
      REAL(KIND=4)  VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C----------------------------------------------------------------------------
C     Mathematical constants:
      REAL(KIND=4)  PI
      PARAMETER(PI=3.1415926536E0)
C
C     Physical constants:
      REAL(KIND=4)  G,KAPPA,NRELIM
      PARAMETER(
C       Gravitational acceleration (m/s**2)
     +          G=9.8E0,
C       Dynamic shape factor for solid particles (Fuchs, 1964)
     +          KAPPA=1.12E0,
C       Lower Reynolds number for Beard scheme calculation:
     +          NRELIM=0.9E-2)    
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      REAL(KIND=4) 
     +        C2S,C2L,C3S,C3L,
     +        BS0,BS1,BS2,BS3,
     +        BL0,BL1,BL2,BL3,BL4,BL5,BL6,
     +        B5,B6,B7,
     +        ASPECT,
     +        ALFA0,BETA0,ALFA1,BETA1,ALFA2,BETA2,ALFA3,BETA3,
     +        REYNLD,X,KN,SCMIDT
      INTEGER I
C
      PARAMETER(
     +        C2S=2.0E0*G/(4.0E0*PI*9.0E0*KAPPA),
     +        C2L=2.0E0*G/(4.0E0*PI*9.0E0),
c                scmidt=(0.71)**(1/3)
     +        SCMIDT=0.892112141E0,
     +        ASPECT=3.0E0,
     +        C3S=2.0E0*G/ASPECT,
     +        C3L=32.0E0*G/(4.0E0*PI),
     +        ALFA0=4.48E-2,BETA0=-1.35E0,
     +        ALFA1=-9.66E-3,BETA1=1.01E0,
     +        ALFA2=-1.81E-3,BETA2=-4.39E-2,
     +        ALFA3=3.26E-4,BETA3=2.04E-5,
     +        BS0=ALFA0*ASPECT+BETA0,
     +        BS1=ALFA1*ASPECT+BETA1,
     +        BS2=ALFA2*ASPECT+BETA2,
     +        BS3=ALFA3*ASPECT+BETA3,
     +        BL0=-3.18657E0,BL1=0.992696E0,BL2=-0.153193E-2,
     +        BL3=-0.987059E-3,BL4=-0.578878E-3,
     +        BL5=0.855176E-4,BL6=-0.327815E-5)
C
C----------------------------------------------------------------------------
C     Internal (statement) functions:
      REAL(KIND=4)  MILKAN,BEARD_SOLID,BEARD_LIQUID
C     MILKAN:       Millikan correction factor to Stokes velocity
C                   (Fuchs 1964, par. 12 )
C     BEARD_SOLID:  Function to calculate Reynolds number from Best number
C                   (P-K, 10-134, 10-136)
C     BEARD_LIQUID: Function to calculate Reynolds number from Best number
C                   (P-K, 10-108, 10-111)
C
      MILKAN(KN)=1.0E0+KN*(1.246E0+0.42E0*EXP(-0.87E0/KN))
      BEARD_SOLID(X)=EXP(BS0+X*(BS1+X*(BS2+X*BS3)))
      BEARD_LIQUID(X)=
     +      EXP(BL0+X*(BL1+X*(BL2+X*(BL3+X*(BL4+X*(BL5+X*BL6))))))
C----------------------------------------------------------------------------
      IF(PHASE.EQ.1) THEN
C         Liquid particles:
          B5=C2L*(ROBULK-RHOAIR)/DYNVIS
          B6=C3L*ROBULK*RHOAIR/(DYNVIS*DYNVIS)
          B7=DYNVIS/(2.0E0*RHOAIR)
          DO I=1,NBINS
              IF(PND(I).LE.ZERO) THEN
                  VELOCI(I)=0.0E0
                  VENTIL(I)=1.0E0
              ELSE
                  REYNLD=ALOG(B6*PTSIZE(I,3))
                  REYNLD=BEARD_LIQUID(REYNLD)
                  IF(REYNLD.LT.NRELIM) THEN
                       VELOCI(I)=B5*PTSIZE(I,2)*
     +                           MILKAN(MFPAIR/PTSIZE(I,1))
                       VENTIL(I)=1.0E0
                  ELSE
                       VELOCI(I)=B7*REYNLD/PTSIZE(I,1)
                       X=SCMIDT*SQRT(REYNLD)
C                      Ventilation factors; P-K, 13-57, 13-58
                       IF(X.LT.1.4) THEN
                            VENTIL(I)=1.0E0+0.108E0*X*X
                       ELSE
                            VENTIL(I)=0.78E0+0.308E0*X
                       ENDIF
                  ENDIF
              ENDIF
          ENDDO
      ELSE IF(PHASE.GT.1) THEN
C         Solid particles:
          B5=C2S*(ROBULK-RHOAIR)/DYNVIS
          B6=C3S*ROBULK*RHOAIR/(DYNVIS*DYNVIS)
          DO I=1,NBINS
              IF(PND(I).LE.ZERO) THEN
                  VELOCI(I)=0.0E0
                  VENTIL(I)=1.0E0
              ELSE
                  REYNLD=ALOG(B6*PTSIZE(I,3))
                  REYNLD=BEARD_SOLID(REYNLD)
                  IF(REYNLD.LT.NRELIM) THEN
                      VELOCI(I)=B5*PTSIZE(I,2)*
     +                          MILKAN(MFPAIR/PTSIZE(I,1))
                      VENTIL(I)=1.0E0
                  ELSE
                      VELOCI(I)=DYNVIS*REYNLD/(2.0E0*RHOAIR*PTSIZE(I,1))
                      X=SCMIDT*SQRT(REYNLD)
C                     Ventilation factors; P-K, 13-87
                      IF(X.LT.1.0E0) THEN
                           VENTIL(I)=1.0E0+0.14E0*X*X
                      ELSE
                           VENTIL(I)=0.86E0+0.28E0*X
                      ENDIF
                  ENDIF
              ENDIF
          ENDDO
      ENDIF
C----------------------------------------------------------------------------
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING TMVELO'
      WRITE(TEST,1010) 'TAIR','PAIR','RHOAIR','DYNVIS','MFPAIR'
      WRITE(TEST,1000) TAIR,PAIR,RHOAIR,DYNVIS,MFPAIR
      WRITE(TEST,*) '(VELOCI(I),I=1,NBINS)'
      WRITE(TEST,1000) (VELOCI(I),I=1,NBINS)
      WRITE(TEST,*) '(VENTIL(I),I=1,NBINS)'
      WRITE(TEST,1000) (VENTIL(I),I=1,NBINS)
 1000 FORMAT(10(1PE9.2,1X),0P)
 1010 FORMAT(1X,8(A9))
C----------------------------------------------------------------------------
      RETURN
      END
C
******************************************************************************
*     SUBROUTINE COAGUL(TAIR,PAIR,ROBULK,NBINS,PTSIZE,
*                       PND,PPND,LPND,
*                       WORK,AKERNL)
******************************************************************************
*
*     This subroutine calculates the production/loss terms of particle 
*     number density in each size bin due to coagulation processes.
*     
*     Input/output variables:                          
*     INTEGER NBINS
*     REAL(KIND=4)  TAIR,PAIR,ROBULK,PTSIZE(NBINS,3),
*          PND(NBINS),PPND(NBINS),LPNP(NBINS),
*          WORK(3*NBINS),AKERNL(NBINS,NBINS)
*
*     Input:
*         TAIR:    (K)            Ambient air temperature
*         PAIR:    (Pa)           Ambient air pressure
*         ROBULK:  (kg/m**3)      Bulk density of particles
*         NBINS:                  Number of radii groups 
*         PTSIZE:  (m,m**2,m**3)  Particle radii,surface,volume
*         PND:     (par/m**3)     Particle number density 
*
*     Output:
*         PPND:    (par/s m**3)   Production term, particle number density
*         LPND:    (1/s)          Loss term, particle number density
*
*     Work array:
*         WORK:                   Work array
*         AKERNL:  (m**3/s par)   Volume corrected coagulation kernel matrix
*
*
      SUBROUTINE COAGUL(TAIR,PAIR,ROBULK,NBINS,PTSIZE,
     +                  PND,PPND,LPND,
     +                  WORK,AKERNL)
      IMPLICIT NONE
      INTEGER NBINS
      REAL(KIND=4)  TAIR,PAIR,ROBULK,PTSIZE(NBINS,3),
     +     PND(NBINS),PPND(NBINS),LPND(NBINS),
     +     WORK(3*NBINS),AKERNL(NBINS,NBINS)
C
C     Auxiliary local variables:
      INTEGER I,J,K,L,ISTOP
      REAL(KIND=4)     X,XP
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C----------------------------------------------------------------------------
      CALL KERNEL(TAIR,PAIR,NBINS,PTSIZE,ROBULK,
     +            WORK(1),WORK(NBINS+1),WORK(2*NBINS+1),AKERNL)
C----------------------------------------------------------------------------
      XP=0.0E0
      DO I=1,NBINS,1
          K=2*I
          L=K-1
          ISTOP=I-1
          WORK(L)=0.0E0
          WORK(K)=XP
          XP=0.0E0
          DO J=1,NBINS,1
              X=AKERNL(J,I)*PND(J)
              WORK(L)=WORK(L)+X
              IF(J.LE.ISTOP) XP=XP+X
          ENDDO
          LPND(I)=WORK(L)
          IF(I.EQ.1) THEN
              PPND(I)=0.0E0
          ELSE
              J=I-1
              PPND(I)=PND(J)*(WORK(K)+0.5E0*AKERNL(J,J)*PND(J))
          ENDIF
      ENDDO
C
C     Loss term in upper bin is set equal to zero to avoid
C     that particles leave upper bin:
      LPND(NBINS)=0.0E0
C
C
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING COAGUL'
      WRITE(TEST,*) '(PPND(I),I=1,NBINS)'
      WRITE(TEST,1000) (PPND(I),I=1,NBINS)
      WRITE(TEST,*) '(LPND(I),I=1,NBINS)'
      WRITE(TEST,1000) (LPND(I),I=1,NBINS)
 1000 FORMAT(10(1PE9.2,1X),0P)
C----------------------------------------------------------------------------
      RETURN
      END
C
******************************************************************************
*     SUBROUTINE KERNEL(TAIR,PAIR,NBINS,PTSIZE,ROBULK,D,G,DEL,AKERNL)
******************************************************************************
*     
*     This subroutine calculates the volume corrected coagulation
*     kernel matrix.
*
*     This matrix is defined as
*
*                            Vj
*                     Kj,i ------                     j < i
*                          (f-1)Vi  
*
*                             f
*          Aj,i =     Kj,j ------                     j = i
*                          2(f-1)
*
*
*                     Kj,i                            j > i
*
*     where Kj,i is the classical coagulation kernel.
*
*
*     Input/output variables:                          
*     INTEGER NBINS
*     REAL(KIND=4)  TAIR,PAIR,PTSIZE(NBINS,3),ROBULK,
*          D(NBINS),G(NBINS),DEL(NBINS),AKERNL(NBINS,NBINS)
*
*     Input:
*         TAIR:    (K)            Ambient air temperature
*         PAIR:    (Pa)           Ambient air pressure
*         NBINS:                  Number of radii groups 
*         PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume
*         ROBULK:  (kg/m**3)      Bulk density of particles
*         D,G,DEL:                Work arrays, each of dimension NBINS
*
*     Output:
*         AKERNL:  (m**3/s par)   Volume corrected coagulation kernel matrix
*
      SUBROUTINE KERNEL(TAIR,PAIR,NBINS,PTSIZE,ROBULK,D,G,DEL,AKERNL)
C
      IMPLICIT NONE
      INTEGER NBINS
      REAL(KIND=4)  TAIR,PAIR,PTSIZE(NBINS,3),ROBULK,
     +     D(NBINS),G(NBINS),DEL(NBINS),AKERNL(NBINS,NBINS)
C
C     Common block variables:
      REAL(KIND=4)  VR,ZERO,D1,LOGVR
      COMMON /BINS/VR,ZERO,D1,LOGVR
C
C     Physical constants:
      REAL(KIND=4)  KBOLTZ
      PARAMETER(
C       Boltzmanns constant (J/K)
     +        KBOLTZ=1.380662E-23) 
C
C     Mathematical constants:
      REAL(KIND=4)  PI
      PARAMETER(PI=3.1415926536E0)
C
C     External thermodynamic functions needed:
      REAL(KIND=4)  FPLAIR,VISAIR
C     FPLAIR:      Mean free path of air molecules
C     VISAIR:      Viscosity of air
C
C     Internal (statement) functions:
      REAL(KIND=4)  SLIPCO
C     SLIPCO: Cunningham slip-flow correction
C             (Pruppacher & Klett 12-16; Hamill et al. 1977)
C
C     Auxiliary local variables:
      REAL(KIND=4)  
     +    DD,GG,DELDEL,RR,LB,
     +    B1,B2,KN,LAIR,
     +    C1,C2,C3,C4,
     +    VRM1
      INTEGER I,J
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C
      PARAMETER(
     +        C1=KBOLTZ/(6.0E0*PI),
     +        C2=8.0E0*KBOLTZ/PI,
     +        C3=8.0E0/PI,
     +        C4=2.0E0*PI)
C----------------------------------------------------------------------------
      SLIPCO(KN)=1.0E0+KN*(1.246E0+0.42E0*EXP(-0.87E0/KN))
C
      B1=C1*TAIR/VISAIR(TAIR)
      B2=C2*TAIR/ROBULK
      LAIR=FPLAIR(TAIR,PAIR)
      VRM1=VR-1.0E0
C
      DO I=1,NBINS,1
          KN=LAIR/PTSIZE(I,1)
          D(I)=B1*SLIPCO(KN)/PTSIZE(I,1)
          G(I)=SQRT(B2/PTSIZE(I,3))
          LB=C3*D(I)/G(I)
          DEL(I)=((2.0E0*PTSIZE(I,1)+LB)**3-
     +             ((PTSIZE(I,2)/PI)+LB*LB)**1.5E0)/
     +             (6.0E0*PTSIZE(I,1)*LB)-2.0E0*PTSIZE(I,1)
      ENDDO
C
      DO I=1,NBINS,1
          DO J=I,NBINS,1
              RR=PTSIZE(I,1)+PTSIZE(J,1)
              DD=D(I)+D(J)
              GG=SQRT(G(I)**2+G(J)**2)
              DELDEL=SQRT(DEL(I)**2+DEL(J)**2)
              AKERNL(J,I)=C4*RR*DD/((RR/(RR+DELDEL))+(4.0E0*DD/(GG*RR)))
              IF(J.GT.I) THEN
                   AKERNL(I,J)=AKERNL(J,I)*PTSIZE(I,3)/
     +                                     (VRM1*PTSIZE(J,3))
              ELSE IF(J.EQ.I) THEN
                   AKERNL(I,J)=0.5E0*AKERNL(J,I)*VR/VRM1
              ENDIF
          ENDDO
      ENDDO
C
      IF(TEST.LT.60.0) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING KERNEL'
      WRITE(TEST,1010) 'TAIR','ROBULK','DD','GG','DELDEL','RR','LB'
      WRITE(TEST,1000) TAIR,ROBULK,DD,GG,DELDEL,RR,LB
      WRITE(TEST,1010) 'B1','B2','KN','LAIR'
      WRITE(TEST,1000) B1,B2,KN,LAIR
      WRITE(TEST,*) 'AKERNEL'
      DO J=1,NBINS
          WRITE(TEST,1000) (AKERNL(J,I),I=1,NBINS)
      ENDDO
      WRITE(TEST,*) '(D(I),I=1,NBINS)'
      WRITE(TEST,1000) (D(I),I=1,NBINS)
      WRITE(TEST,*) '(G(I),I=1,NBINS)'
      WRITE(TEST,1000) (G(I),I=1,NBINS)
      WRITE(TEST,*) '(DEL(I),I=1,NBINS)'
      WRITE(TEST,1000) (DEL(I),I=1,NBINS)
 1000 FORMAT(10(1PE9.2,1X),0P)
 1010 FORMAT(1X,8(A9))
C----------------------------------------------------------------------------
      RETURN
      END
C
******************************************************************************
*     SUBROUTINE CMPTER(TAIR,PAIR,MMRSA,PPWV,PPNA,WSA,WNA)
******************************************************************************
*
*     This subroutine calculates the equilibrium composition of 
*     stratospheric liquid ternary aewrosols (sulfuric/nitric-acid,water)
*     in equilibrium with the ambient partial pressure of nitric acid 
*     and water vapor.
*
*     The aerosol sulfate mass mixing ratio (MMRSA, i.e. kg. of H2SO4,
*     dissolved in the aerosol liquid phase, per kg. air) is assumed 
*     constant. On entry the aerosol with a liquid composition, given
*     by the nitric and sulfuric acid weight fractions (WNA and WSA), is
*     possibly not in equilibrium with the ambient nitric acid partial
*     pressure (PPNA) at the ambient temperature TAIR. On exit the new
*     equilibrium composition (WNA and WSA) and the new partial nitric 
*     acid pressure is calculated. The difference in nitric acid partial 
*     pressure on entry and exit corresponds to the uptake or release 
*     of HNO3 from the condensed phase.
*
*     The exit-values of WNA and WSA are assumed to be used as input
*     to the subroutine in a subsequent call to calculate the composition
*     of the same aerosol particles at different ambient conditions.
*
*     Subroutine WGTTER can be used to calculate the initial values
*     of WNA and WSA. Either the constant mass mixing ration of dissolved
*     H2SO4 is known a priori or it can be calculated from the the
*     size distribition of the aerosol in the initial condition in the
*     following way: calculate the aerosol particle volume per kg. air (V),
*     then MMRSA = WSA*RHOTER*V, where WSA is the initial H2SO4 weight
*     fraction, and RHOTER is the initial density of the condensed phase,
*     calculated from FUNCTION ROSNW, using the initial composition values.
*
*     Vapor pressures: Source: Luo et al.: GRL 22,247,1995
*
*     The equilibrium compositions (WSA,WNA) are calculated by solving
*             PPNA = f(TAIR,WNA,WSA) and PPWV = g(TAIR,WNA,WSA)
*     where f and g are nitric acid and water vapor saturation pressures
*     over the ternary solution. Equations are solved by Newton-Ralphson
*     iteration.
*
*     The Kelvin effect is not included.
*
*     Range of validity: 1.0E-4 < WNA;WSA < 0.7
*
*
*     Input/output variables:
*     REAL(KIND=4)  TAIR,PAIR,MMRSA,PPWV,PPNA,WSA,WNA
*
*     Input:       
*         TAIR:    K         Temperature of ambient air 
*         PAIR:    Pa        Ambient air pressure
*         MMRSA:   kg/kg air Mass mixing ratio of dissolved H2SO4
*         PPWV:    Pa        Partial pressure of ambient water vapor 
*
*     Input/Output:
*         PPNA:    Pa        Partial pressure of ambient nitric acid vapor 
*         WSA:               Mass fraction of sulfuric acid.
*         WNA:               Mass fraction of nitric acid.
*
      SUBROUTINE CMPTER(TAIR,PAIR,MMRSA,PPWV,PPNA,WSA,WNA)
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PAIR,MMRSA,PPWV,PPNA,WSA,WNA
      REAL(KIND=4)  TOLF,TOLG,TOLX,TOLXY,XOLD,YOLD,
     +        WSGUES,WNGUES,
     +        XMIN,XMAX,YMIN,YMAX,DX,DY,DDX,C1,C2,XX,YY,
     +        MAIR,MHNO3
      PARAMETER(TOLG=0.0001E0,TOLX=1.0E-5,TOLXY=1.0E-4,
     +          WNGUES=0.3E0,WSGUES=0.20E0,
     +          XMIN=1.0E-3,XMAX=0.7E0,
     +          YMIN=1.0E-3,YMAX=0.9E0,
     +          DX=0.0001E0,DY=0.0001E0,
     +          MAIR=28.9644E-3,
     +          MHNO3=63.01E-3)
      REAL(KIND=4)  X,Y,F,G,DFDX,DGDX,DGDY,DELX,DELY,
     +        ALFA,BETA,
     +        LOGPW,PNASNW,PWVSNW,X0G,Y0G,PSNA
      INTEGER ITERAT,MAXITERAT
      PARAMETER(MAXITERAT=100)
CC      LOGICAL INTERA
CC      PARAMETER(INTERA=.TRUE.)
CC      PARAMETER(INTERA=.FALSE.)
C
      ITERAT=0
      LOGPW=ALOG(AMAX1(PPWV,1.0E-20))
      G=PWVSNW(TAIR,YMIN,XMIN)-LOGPW
C
C     Temperature below range of validity:
      IF(G.LE.TOLG) THEN
CC          IF(INTERA) WRITE(*,*) 'A-ITERAT: ',ITERAT
cc          WRITE(70,*) 'exit=1'
          RETURN
      ENDIF
C
      C1=PAIR*MAIR*MMRSA/MHNO3
      C2=WNA/WSA
C
CC      IF(INTERA) ITERAT=0
C
C     Calculate g(y)=0 for x=xmin:
      X=XMIN
      Y=WSGUES
      YOLD=1.0
  2   CONTINUE
      ITERAT=ITERAT+1
      G=PWVSNW(TAIR,Y,X)-LOGPW
      IF(ABS(G).GT.TOLG.AND.ABS(YOLD-Y).GT.TOLXY) THEN
          DGDY=(PWVSNW(TAIR,Y+DY,X)-LOGPW-G)/DY
          DELY=-G/DGDY
CC         IF(INTERA.AND.ITERAT.GT.500) THEN
CC           WRITE(*,*) TAIR,PAIR,MMRSA,PPWV,PPNA
CC             WRITE(*,*) 'X=0 Point on g-line: X,Y',X,Y
CC             WRITE(*,*) 'G-DELY',G,DELY
CC             WRITE(*,*) 'ITERAT= ',ITERAT
CC             PAUSE
CC         ENDIF
          IF((Y.LE.YMIN.AND.DELY.LT.0.0E0).OR.
     +       (Y.GE.YMAX.AND.DELY.GT.0.0E0)) THEN
              Y0G=Y
          ELSE
              YOLD=Y
              Y=AMAX1(YMIN,AMIN1(Y+DELY,YMAX))
              IF(ITERAT.LT.MAXITERAT) THEN
                  GOTO 2
              ELSE
                  Y0G=YOLD
                  ITERAT=0
              ENDIF
          ENDIF
      ELSE
CC         IF(INTERA.AND.ITERAT.GT.500) THEN
CC             WRITE(*,*) TAIR,PAIR,MMRSA,PPWV,PPNA
CC             WRITE(*,*) 'Y0G Point',X,Y
CC             WRITE(*,*) 'G-DELY',G,DELY
CC             WRITE(*,*) 'ITERAT= ',ITERAT
CC             PAUSE
CC         ENDIF
          Y0G=Y
      ENDIF
C
      PSNA=EXP(PNASNW(TAIR,Y0G,XMIN))
      F=PSNA-PPNA+C1*(XMIN/Y0G-C2)
      IF(F.GE.0.0E0) THEN
C         High temperatures, binary H2SO4/H2O solution:
          WSA=Y0G
          WNA=XMIN
          PPNA=AMAX1(PPNA-C1*(WNA/WSA-C2),0.0)
CC          IF(INTERA) WRITE(*,*) 'B-ITERAT: ',ITERAT
cc          WRITE(70,*) 'exit=2'
          RETURN
      ENDIF
C
C     Calculate g(x)=0 for y=ymin:
      Y=YMIN
      X=WNGUES
      XOLD=1.0
  1   CONTINUE
      ITERAT=ITERAT+1
      G=PWVSNW(TAIR,Y,X)-LOGPW
CC      IF(INTERA) ITERAT=ITERAT+1
      IF(ABS(G).GT.TOLG.AND.ABS(X-XOLD).GT.TOLXY) THEN
          DGDX=(PWVSNW(TAIR,Y,X+DX)-LOGPW-G)/DX
          DELX=-G/DGDX
CC         IF(INTERA.AND.ITERAT.GT.500) THEN
CC             WRITE(*,*) TAIR,PAIR,MMRSA,PPWV,PPNA
CC             WRITE(*,*) 'Y=0 Point on g-line: X,Y',X,Y
CC             WRITE(*,*) 'G-DELX',G,DELX
CC             WRITE(*,*) 'ITERAT= ',ITERAT
CC             PAUSE
CC         ENDIF
          IF((X.LE.XMIN.AND.DELX.LT.0.0E0).OR.
     +       (X.GE.XMAX.AND.DELX.GT.0.0E0)) THEN
              X0G=X
          ELSE
              XOLD=X
              X=AMAX1(XMIN,AMIN1(X+DELX,XMAX))
              IF(ITERAT.LT.MAXITERAT) THEN
                  GOTO 1
              ELSE
                  X0G=XOLD
                  ITERAT=0
              ENDIF
          ENDIF
      ELSE
CC         IF(INTERA.AND.ITERAT.GT.500) THEN
CC             WRITE(*,*) TAIR,PAIR,MMRSA,PPWV,PPNA
CC             WRITE(*,*) 'X0G Point',X,Y
CC             WRITE(*,*) 'G-DELX',G,DELX
CC             WRITE(*,*) 'ITERAT= ',ITERAT
CC             PAUSE
CC         ENDIF
          X0G=X
      ENDIF 
C
C
      PSNA=EXP(PNASNW(TAIR,YMIN,X0G))
      F=PSNA-PPNA+C1*(X0G/YMIN-C2)
      IF(F.LE.0.0E0) THEN
C         Low temperatures, binary HNO3/H2O solution:
          WSA=YMIN
          WNA=X0G
          PPNA=AMAX1(PPNA-C1*(WNA/WSA-C2),0.0)
CC          IF(INTERA) WRITE(*,*) 'C-ITERAT: ',ITERAT
cc          WRITE(70,*) 'exit=3'
          RETURN
      ENDIF
C
C     Linear expression for g:
      ALFA=(Y0G-YMIN)/(XMIN-X0G)
      BETA=Y0G-ALFA*XMIN
C
CC      X0G=X0G-2.0E0*DX
C
C     ------------------------
C     open(2,file='TESTTER5.DAT')
C     IF(INTERA) THEN
C         DO X=XMIN,X0G,0.001E0
C         Y=ALFA*X+BETA
C         PSNA=EXP(PNASNW(TAIR,Y,X))
C         TOLF=0.01E0*PSNA
C         F=PSNA-PPNA+C1*(X/Y-C2)
C         WRITE(2,*) X,F,PSNA,C1*(X/Y-C2),PPNA
C         ENDDO
C         CLOSE(2)
C     ENDIF
C---------------------------------
C     Point on f-zero line:
      X=AMIN1(WNA,X0G-4.0E0*DX)
      XOLD=1.0
  13  CONTINUE
      ITERAT=ITERAT+1
      Y=ALFA*X+BETA
      PSNA=EXP(PNASNW(TAIR,Y,X))
      F=PSNA-PPNA+C1*(X/Y-C2)
      TOLF=0.01E0*PSNA
CC      IF(INTERA) ITERAT=ITERAT+1
      IF(ABS(F).GT.TOLF.AND.ABS(X-XOLD).GT.TOLXY) THEN
          DDX=SIGN(DX,Y+ALFA*DX)
          XX=X+DDX
          YY=ALFA*XX+BETA
          PSNA=EXP(PNASNW(TAIR,YY,XX))
          DFDX=(PSNA-PPNA+C1*(XX/YY-C2)-F)/DDX
          DELX=-F/DFDX
CC         IF(INTERA.AND.ITERAT.GT.500) THEN
CC             WRITE(*,*) TAIR,PAIR,MMRSA,PPWV,PPNA
CC             WRITE(*,*) 'New point on f-zero line:'
CC             WRITE(*,*) 'X-Y',X,Y
CC             WRITE(*,*) 'F',F
CC             WRITE(*,*) 'DFDX: ',DFDX,'  DELX:',DELX
CC             WRITE(*,*) 'ITERAT= ',ITERAT
CC             PAUSE
CC         ENDIF
          IF(X.LE.XMIN.AND.DELX.LT.0.0E0) THEN
              WNA=X
              WSA=Y0G
              PPNA=AMAX1(PPNA-C1*(WNA/WSA-C2),0.0)
CC              IF(INTERA) WRITE(*,*) 'D-ITERAT: ',ITERAT
cc              WRITE(70,*) 'exit=4'
              RETURN
          ENDIF
          XOLD=X
          X=AMAX1(XMIN,AMIN1(X+DELX,X0G))
          IF(ABS(DELX).LT.TOLX) THEN
              WNA=X
              WSA=ALFA*X+BETA
              PPNA=AMAX1(PPNA-C1*(WNA/WSA-C2),0.0)
CC              IF(INTERA) WRITE(*,*) 'F-ITERAT: ',ITERAT
cc              WRITE(70,*) 'exit=5'
              RETURN
          ELSE
              IF(ITERAT.LT.MAXITERAT) THEN
                  GOTO 13
              ELSE
                  WNA=X
                  WSA=ALFA*X+BETA
                  PPNA=AMAX1(PPNA-C1*(WNA/WSA-C2),0.0)
cc                WRITE(70,*) 'exit=6'
                  RETURN
              ENDIF
          ENDIF
      ELSE
          WNA=X
          WSA=Y
          PPNA=AMAX1(PPNA-C1*(X/Y-C2),0.0)
CC          IF(INTERA) WRITE(*,*) 'G-ITERAT: ',ITERAT
cc          WRITE(70,*) 'exit=7'
          RETURN
      ENDIF
      END
******************************************************************************
*     SUBROUTINE WGTTER(TAIR,PPWV,PPNA,WSA,WNA)
******************************************************************************
*
*     This subroutine calculates the equilibrium composition of 
*     stratospheric liquid ternary solution (sulfuric/nitric-acid,water)
*     in equilibrium with the ambient partial nitric acid and water vapor.
*
*     Source: Luo et al.: GRL 22,247,1995
*
*     The equilibrium compositions (WSA,WNA) are calculated by solving
*             PPNA = f(TAIR,WNA,WSA) and PPWV = g(TAIR,WNA,WSA)
*     where f and g are nitric acid and water vapor saturation pressures
*     over the ternary solution. Equations are solved by Newton-Ralphson
*     iteration.
*
*     The Kelvin effect is NOT included.
*
*     Range of validity: 1.0E-4 < WNA;WSA < 0.7
*
*
*     Input/output variables:
*     REAL(KIND=4)  TAIR,PPWV,PPNA,WSA,WNA
*
*     Input:       
*         TAIR:    K         Temperature of ambient air 
*         PPWV:    Pa        Partial pressure of ambient water vapor 
*         PPNA:    Pa        Partial pressure of ambient nitric acid vapor 
*
*     Output:
*         WSA:               Mass fraction of sulfuric acid.
*         WNA:               Mass fraction of nitric acid.
*
      SUBROUTINE WGTTER(TAIR,PPWV,PPNA,WSA,WNA)
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PPWV,PPNA,WSA,WNA
      REAL(KIND=4)  TOLF,TOLG,TOLXY,
     +        WSGUES,WNGUES,XOLD,YOLD,
     +        XMIN,XMAX,YMIN,YMAX,DX,DY,DDX
      PARAMETER(TOLF=0.001E0,TOLG=0.0001E0,TOLXY=1.0E-4,
     +          WNGUES=0.3E0,WSGUES=0.20E0,
     +          XMIN=1.0E-3,XMAX=0.7E0,
     +          YMIN=1.0E-3,YMAX=0.9E0,
     +          DX=0.0001E0,DY=0.0001E0)
      REAL(KIND=4)  X,Y,F,G,DFDX,DGDX,DGDY,DELX,DELY,
     +        ALFA,BETA,
     +        LOGPW,LOGPN,PNASNW,PWVSNW,X0G,Y0G
      INTEGER ITERAT,MAXITERAT
      PARAMETER(MAXITERAT=100)
CC    LOGICAL INTERA
CC      PARAMETER(INTERA=.TRUE.)
CC    PARAMETER(INTERA=.FALSE.)
C
C
      ITERAT=0
      LOGPW=ALOG(AMAX1(PPWV,1.0E-20))
      G=PWVSNW(TAIR,YMIN,XMIN)-LOGPW
      IF(G.LE.TOLG) THEN
          WSA=YMIN
          WNA=YMIN
          RETURN
      ENDIF
C
      LOGPN=ALOG(AMAX1(PPNA,1.0E-20))
CC    ITERAT=0

C
C     Calculate g(x)=0 for y=ymin:
      Y=YMIN
      X=WNGUES
      XOLD=1.0
  1   CONTINUE
      ITERAT=ITERAT+1
      G=PWVSNW(TAIR,Y,X)-LOGPW
CC    ITERAT=ITERAT+1
      IF(ABS(G).GT.TOLG.AND.ABS(XOLD-X).GT.TOLXY) THEN
          DGDX=(PWVSNW(TAIR,Y,X+DX)-LOGPW-G)/DX
          DELX=-G/DGDX
CC        IF(INTERA) THEN
CC            WRITE(*,*) 'Y=0 Point on g-line: X,Y',X,Y
CC            WRITE(*,*) 'G-DELX',G,DELX
CC            WRITE(*,*) 'ITERAT= ',ITERAT
CC            PAUSE
CC        ENDIF
          IF((X.LE.XMIN.AND.DELX.LT.0.0E0).OR.
     +       (X.GE.XMAX.AND.DELX.GT.0.0E0)) THEN
              X0G=X
          ELSE
              XOLD=X
              X=AMAX1(XMIN,AMIN1(X+DELX,XMAX))
              IF(ITERAT.LT.MAXITERAT) THEN
                  GOTO 1
              ELSE
                  X0G=XOLD
                  ITERAT=0
              ENDIF
          ENDIF
      ELSE
CC        IF(INTERA) THEN
CC            F=PNASNW(TAIR,Y,X)-LOGPN
CC            WRITE(*,*) 'X0G Point',X,Y
CC            WRITE(*,*) 'G-DELX',G,DELX,'  F: ',F
CC            WRITE(*,*) 'ITERAT= ',ITERAT
CC            PAUSE
CC        ENDIF
          X0G=X
      ENDIF 
C
      F=PNASNW(TAIR,YMIN,X0G)-LOGPN
      IF(F.LE.TOLF) THEN
C         Binary HNO3/H2O solution:
          WNA=X0G
          WSA=YMIN
          RETURN
      ENDIF
C
C     Calculate g(y)=0 for x=xmin:
      X=XMIN
      Y=WSGUES
      YOLD=1.0
  2   CONTINUE
      ITERAT=ITERAT+1
      G=PWVSNW(TAIR,Y,X)-LOGPW
CC    ITERAT=ITERAT+1
      IF(ABS(G).GT.TOLG.AND.ABS(YOLD-Y).GT.TOLXY) THEN
          DGDY=(PWVSNW(TAIR,Y+DY,X)-LOGPW-G)/DY
          DELY=-G/DGDY
CC        IF(INTERA) THEN
CC            WRITE(*,*) 'X=0 Point on g-line: X,Y',X,Y
CC            WRITE(*,*) 'G-DELY',G,DELY
CC            WRITE(*,*) 'ITERAT= ',ITERAT
CC            PAUSE
CC        ENDIF
          IF((Y.LE.YMIN.AND.DELY.LT.0.0E0).OR.
     +       (Y.GE.YMAX.AND.DELY.GT.0.0E0)) THEN
              Y0G=Y
          ELSE
              YOLD=Y
              Y=AMAX1(YMIN,AMIN1(Y+DELY,YMAX))
              IF(ITERAT.LT.MAXITERAT) THEN
                  GOTO 2
              ELSE
                  Y0G=YOLD
                  ITERAT=0
              ENDIF
          ENDIF
      ELSE
CC        IF(INTERA) THEN
CC            WRITE(*,*) 'Y0G Point',X,Y
CC            WRITE(*,*) 'G-DELY',G,DELY
CC            WRITE(*,*) 'ITERAT= ',ITERAT
CC            PAUSE
CC        ENDIF
          Y0G=Y
      ENDIF
C
      F=PNASNW(TAIR,Y0G,XMIN)-LOGPN
      IF(F.GE.0.0E0) THEN
C         Binary H2SO4/H2O solution:
          WNA=XMIN
          WSA=Y0G
          RETURN
      ENDIF
C
C     Calculate crossing of f and g:
      X=X0G/2.0E0
      Y=Y0G/2.0E0
C
C     Linear expression for g:
      ALFA=(Y0G-YMIN)/(XMIN-X0G)
      BETA=Y0G-ALFA*XMIN
C
C     Point on f-zero line:
      XOLD=1.0
  13  CONTINUE
      ITERAT=ITERAT+1
      F=PNASNW(TAIR,Y,X)-LOGPN
CC    ITERAT=ITERAT+1
      IF(ABS(F).GT.TOLF.AND.ABS(XOLD-X).GT.TOLXY) THEN
          DDX=SIGN(DX,Y+ALFA*DX)
          DFDX=(PNASNW(TAIR,Y+ALFA*DDX,X+DDX)-LOGPN-F)/DDX
          DELX=-F/DFDX
CC        IF(INTERA) THEN
CC            G=PWVSNW(TAIR,Y,X)-LOGPW
CC            WRITE(*,*) 'New point on f-zero line:'
CC            WRITE(*,*) 'X-Y',X,Y
CC            WRITE(*,*) 'F-G',F,G
CC            WRITE(*,*) 'DFDX: ',DFDX,'  DELX:',DELX
CC            WRITE(*,*) 'ITERAT= ',ITERAT
CC            PAUSE
CC        ENDIF
          XOLD=X
          X=AMAX1(XMIN,AMIN1(X+DELX,X0G))
C         Assume y to be constrained by linear relationship on g-zero line:
          Y=AMAX1(YMIN,AMIN1(ALFA*X+BETA,YMAX))
          IF(ITERAT.LT.MAXITERAT) THEN
              GOTO 13
          ELSE
              WNA=XOLD
              WSA=Y
              RETURN
          ENDIF
      ELSE
          WNA=X
          WSA=Y
          RETURN
      ENDIF
      END
******************************************************************************
*     SUBROUTINE WGTTER_KELVIN(RADIUS,TAIR,PPWV,PPNA,WSA,WNA,ROBULK)
******************************************************************************
*
*     This subroutine calculates the equilibrium composition of a
*     stratospheric liquid ternary solution (sulfuric/nitric-acid,water)
*     droplet in equilibrium with the ambient partial nitric acid and
*     water vapor.
*
*     Source: Luo et al.: GRL 22,247,1995
*
*     The equilibrium compositions (WSA,WNA) are calculated by solving
*             PPNA = f(TAIR,WNA,WSA) and PPWV = g(TAIR,WNA,WSA)
*     where f and g are nitric acid and water vapor saturation pressures
*     over the ternary solution. Equations are solved by Newton-Ralphson
*     iteration.
*
*     The Kelvin effect is included.
*
*     Range of validity: 1.0E-4 < WNA;WSA < 0.7
*
*
*     Input/output variables:
*     REAL(KIND=4)  TAIR,PPWV,PPNA,WSA,WNA
*
*     Input:       
*         RADIUS:  m         Particle radius
*         TAIR:    K         Temperature of ambient air
*         PPWV:    Pa        Partial pressure of ambient water vapor 
*         PPNA:    Pa        Partial pressure of ambient nitric acid vapor 
*
*     Output:
*         WSA:               Mass fraction of sulfuric acid.
*         WNA:               Mass fraction of nitric acid.
*         ROBULK   kg/m**3   Density
*
      SUBROUTINE WGTTER_KELVIN(RADIUS,TAIR,PPWV,PPNA,WSA,WNA,ROBULK)
      IMPLICIT NONE
      REAL(KIND=4)  RADIUS,TAIR,PPWV,PPNA,WSA,WNA,ROBULK
      REAL(KIND=4)  TOLF,TOLG,TOLXY,
     +        WSGUES,WNGUES,XOLD,YOLD,
     +        XMIN,XMAX,YMIN,YMAX,DX,DY,DDX
      PARAMETER(TOLF=0.0001E0,TOLG=0.0001E0,TOLXY=1.0E-4,
     +          WNGUES=0.3E0,WSGUES=0.20E0,
     +        XMIN=1.0E-3,XMAX=0.7E0,
     +        YMIN=1.0E-3,YMAX=0.9E0,
     +        DX=0.0001E0,DY=0.0001E0)
      REAL(KIND=4)  X,Y,F,G,DFDX,DGDX,DGDY,DELX,DELY,
     +        ALFA,BETA,
     +        LOGPW,LOGPN,PNASNW,PWVSNW,
     +        STSNW,ROSNW,RO,
     +        X0G,Y0G
C
C     Physical constants:
      REAL(KIND=4)  RGAS,MH2O,MHNO3,CW,CN,CCW,CCN
      PARAMETER(
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
C       Molar weight of nitric acid (kg/mole)
     +          MHNO3=63.01E-3,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0,
     +          CW=2.0E0*MH2O/RGAS,
     +          CN=2.0E0*MHNO3/RGAS)
C
      INTEGER ITERAT,MAXITERAT
      PARAMETER(MAXITERAT=100)
CC    LOGICAL INTERA
CC      PARAMETER(INTERA=.TRUE.)
CC    PARAMETER(INTERA=.FALSE.)
C
C
      ITERAT=0
      CCW=CW/(TAIR*RADIUS)
      CCN=CN/(TAIR*RADIUS)
      LOGPW=ALOG(AMAX1(PPWV,1.0E-20))
      RO=ROSNW(TAIR,YMIN,XMIN)
      G=PWVSNW(TAIR,YMIN,XMIN)
     +  +CCW*STSNW(TAIR,YMIN,XMIN)/RO
     +  -LOGPW
      IF(G.LE.TOLG) THEN
          WSA=YMIN
          WNA=YMIN
          ROBULK=RO
          RETURN
      ENDIF
C
      LOGPN=ALOG(AMAX1(PPNA,1.0E-20))
CC    ITERAT=0

C
C     Calculate g(x)=0 for y=ymin:
      Y=YMIN
      X=WNGUES
      XOLD=1.0
  1   CONTINUE
      ITERAT=ITERAT+1
      RO=ROSNW(TAIR,Y,X)
      G=PWVSNW(TAIR,Y,X)
     +  +CCW*STSNW(TAIR,Y,X)/RO
     +  -LOGPW
CC    ITERAT=ITERAT+1
      IF(ABS(G).GT.TOLG.AND.ABS(XOLD-X).GT.TOLXY) THEN
          RO=ROSNW(TAIR,Y,X+DX)
          DGDX=(PWVSNW(TAIR,Y,X+DX)
     +          +CCW*STSNW(TAIR,Y,X+DX)/RO
     +          -LOGPW-G)/DX
          DELX=-G/DGDX
CC        IF(INTERA) THEN
CC            WRITE(*,*) 'Y=0 Point on g-line: X,Y',X,Y
CC            WRITE(*,*) 'G-DELX',G,DELX
CC            WRITE(*,*) 'ITERAT= ',ITERAT
CC            PAUSE
CC        ENDIF
          IF((X.LE.XMIN.AND.DELX.LT.0.0E0).OR.
     +       (X.GE.XMAX.AND.DELX.GT.0.0E0)) THEN
              X0G=X
          ELSE
              XOLD=X
              X=AMAX1(XMIN,AMIN1(X+DELX,XMAX))
              IF(ITERAT.LT.MAXITERAT) THEN
                  GOTO 1
              ELSE
                  X0G=XOLD
                  ITERAT=0
              ENDIF
          ENDIF
      ELSE
CC        IF(INTERA) THEN
CC            RO=ROSNW(TAIR,Y,X)
CC            F=PNASNW(TAIR,Y,X)
CC   +          +CCN*STSNW(TAIR,Y,X)/RO
CC   +          -LOGPN
CC            WRITE(*,*) 'X0G Point',X,Y
CC            WRITE(*,*) 'G-DELX',G,DELX,'  F: ',F
CC            WRITE(*,*) 'ITERAT= ',ITERAT
CC            PAUSE
CC        ENDIF
          X0G=X
      ENDIF 
C
      RO=ROSNW(TAIR,YMIN,X0G)
      F=PNASNW(TAIR,YMIN,X0G)
     +  +CCN*STSNW(TAIR,YMIN,X0G)/RO
     +  -LOGPN
      IF(F.LE.TOLF) THEN
C         Binary HNO3/H2O solution:
          WNA=X0G
          WSA=YMIN
          ROBULK=RO
          RETURN
      ENDIF
C
C     Calculate g(y)=0 for x=xmin:
      X=XMIN
      Y=WSGUES
      YOLD=1.0
  2   CONTINUE
      ITERAT=ITERAT+1
      RO=ROSNW(TAIR,Y,X)
      G=PWVSNW(TAIR,Y,X)
     +  +CCW*STSNW(TAIR,Y,X)/RO
     +  -LOGPW
CC    ITERAT=ITERAT+1
      IF(ABS(G).GT.TOLG.AND.ABS(YOLD-Y).GT.TOLXY) THEN
          RO=ROSNW(TAIR,Y+DY,X)
          DGDY=(PWVSNW(TAIR,Y+DY,X)
     +          +CCW*STSNW(TAIR,Y+DY,X)/RO
     +          -LOGPW-G)/DY
          DELY=-G/DGDY
CC        IF(INTERA) THEN
CC            WRITE(*,*) 'X=0 Point on g-line: X,Y',X,Y
CC            WRITE(*,*) 'G-DELY',G,DELY
CC            WRITE(*,*) 'ITERAT= ',ITERAT
CC            PAUSE
CC        ENDIF
          IF((Y.LE.YMIN.AND.DELY.LT.0.0E0).OR.
     +       (Y.GE.YMAX.AND.DELY.GT.0.0E0)) THEN
              Y0G=Y
          ELSE
              YOLD=Y
              Y=AMAX1(YMIN,AMIN1(Y+DELY,YMAX))
              IF(ITERAT.LT.MAXITERAT) THEN
                  GOTO 2
              ELSE
                  Y0G=YOLD
                  ITERAT=0
              ENDIF
          ENDIF
      ELSE
CC        IF(INTERA) THEN
CC            WRITE(*,*) 'Y0G Point',X,Y
CC            WRITE(*,*) 'G-DELY',G,DELY
CC            WRITE(*,*) 'ITERAT= ',ITERAT
CC            PAUSE
CC        ENDIF
          Y0G=Y
      ENDIF
C
      RO=ROSNW(TAIR,Y0G,XMIN)
      F=PNASNW(TAIR,Y0G,XMIN)
     +  +CCN*STSNW(TAIR,Y0G,XMIN)/RO
     +  -LOGPN
      IF(F.GE.0.0E0) THEN
C         Binary H2SO4/H2O solution:
          WNA=XMIN
          WSA=Y0G
          ROBULK=RO
          RETURN
      ENDIF
C
C     Calculate crossing of f and g:
      X=X0G/2.0E0
      Y=Y0G/2.0E0
C
C     Linear expression for g:
      ALFA=(Y0G-YMIN)/(XMIN-X0G)
      BETA=Y0G-ALFA*XMIN
C
C     Point on f-zero line:
      XOLD=1.0
  13  CONTINUE
      ITERAT=ITERAT+1
      RO=ROSNW(TAIR,Y,X)
      F=PNASNW(TAIR,Y,X)
     + +CCN*STSNW(TAIR,Y,X)/RO
     + -LOGPN
CC    ITERAT=ITERAT+1
      IF(ABS(F).GT.TOLF.AND.ABS(XOLD-X).GT.TOLXY) THEN
          DDX=SIGN(DX,Y+ALFA*DX)
          RO=ROSNW(TAIR,Y+ALFA*DDX,X+DDX)
          DFDX=(PNASNW(TAIR,Y+ALFA*DDX,X+DDX)
     +          +CCN*STSNW(TAIR,Y+ALFA*DDX,X+DDX)/RO
     +          -LOGPN-F)/DDX
          DELX=-F/DFDX
CC        IF(INTERA) THEN
CC            RO=ROSNW(TAIR,Y,X)
CC            G=PWVSNW(TAIR,Y,X)
CC   +          +CCW*STSNW(TAIR,Y,X)/RO
CC   +          -LOGPW
CC            WRITE(*,*) 'New point on f-zero line:'
CC            WRITE(*,*) 'X-Y',X,Y
CC            WRITE(*,*) 'F-G',F,G
CC            WRITE(*,*) 'DFDX: ',DFDX,'  DELX:',DELX
CC            WRITE(*,*) 'ITERAT= ',ITERAT
CC            PAUSE
CC        ENDIF
          XOLD=X
          X=AMAX1(XMIN,AMIN1(X+DELX,X0G))
C         Assume y to be constrained by linear relationship on g-zero line:
          Y=AMAX1(YMIN,AMIN1(ALFA*X+BETA,YMAX))
          IF(ITERAT.LT.MAXITERAT) THEN
              GOTO 13
          ELSE
              WNA=XOLD
              WSA=Y
              ROBULK=RO
              RETURN
          ENDIF
      ELSE
          WNA=X
          WSA=Y
          ROBULK=RO
          RETURN
      ENDIF
      END
******************************************************************************
*     SUBROUTINE WS_WN_REL(RADIUS,TAIR,PPWV,ALFA,BETA)
******************************************************************************
*
*     This subroutine calculates coefficients, ALFA and BETA, of a linear
*     relation
*                    WNA = ALFA*WSA+BETA
*     between nitric (WNA) and sulfuric (WSA) acid weight fractions of
*     stratospheric liquid ternary solution (sulfuric/nitric-acid,water)
*     particles with a given radius in equilibrium with the ambient partial
*     water vapor. If supersaturation with respect to pure water,
*     ALFA=BETA=0 on exit.
*
*     The Kelvin effect is included.
*
*     Source: Luo et al.: GRL 22,247,1995
*
*     Range of validity: 1.0E-4 < WNA;WSA < 0.7
*
*
*     Input/output variables:
*     REAL(KIND=4)  RADIUS,TAIR,PPWV,PPNA,WSA,WNA
*
*     Input:       
*         RADIUS:  m         Particle radius
*         TAIR:    K         Temperature of ambient air
*         PPWV:    Pa        Partial pressure of ambient water vapor 
*         PPNA:    Pa        Partial pressure of ambient nitric acid vapor 
*
*     Output:
*         ALFA:
*         BETA:
*
      SUBROUTINE WS_WN_REL(RADIUS,TAIR,PPWV,ALFA,BETA)
      IMPLICIT NONE
      REAL(KIND=4)  RADIUS,TAIR,PPWV,ALFA,BETA
      REAL(KIND=4)  TOLG,TOLXY,
     +        WSGUES,WNGUES,XOLD,YOLD,
     +        XMIN,XMAX,YMIN,YMAX,DX,DY
      PARAMETER(TOLG=0.0001E0,TOLXY=1.0E-4,
     +          WNGUES=0.3E0,WSGUES=0.20E0,
     +        XMIN=1.0E-3,XMAX=0.7E0,
     +        YMIN=1.0E-3,YMAX=0.9E0,
     +        DX=0.0001E0,DY=0.0001E0)
      REAL(KIND=4)  X,Y,G,DGDX,DGDY,DELX,DELY,
     +        LOGPW,PWVSNW,
     +        STSNW,ROSNW,
     +        X0G,Y0G
      INTEGER ITERAT,MAXITERAT
      PARAMETER(MAXITERAT=100)
C
C     Physical constants:
      REAL(KIND=4)  RGAS,MH2O,CW,CCW
      PARAMETER(
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0,
     +          CW=2.0*MH2O/RGAS)
C
C
      ITERAT=0
      CCW=CW/(TAIR*RADIUS)
      LOGPW=ALOG(AMAX1(PPWV,1.0E-20))

      G=PWVSNW(TAIR,YMIN,XMIN)
     +  +CCW*STSNW(TAIR,YMIN,XMIN)/ROSNW(TAIR,YMIN,XMIN)
     +  -LOGPW
      IF(G.LE.TOLG) THEN
          ALFA=0.0E0
          BETA=0.0E0
          RETURN
      ENDIF
C
C     Calculate g(x)=0 for y=ymin:
      Y=YMIN
      X=WNGUES
      XOLD=1.0
  1   CONTINUE
      ITERAT=ITERAT+1
      G=PWVSNW(TAIR,Y,X)
     +  +CCW*STSNW(TAIR,Y,X)/ROSNW(TAIR,Y,X)
     +  -LOGPW
      IF(ABS(G).GT.TOLG.AND.ABS(XOLD-X).GT.TOLXY) THEN
          DGDX=(PWVSNW(TAIR,Y,X+DX)
     +          +CCW*STSNW(TAIR,Y,X+DX)/ROSNW(TAIR,Y,X+DX)
     +          -LOGPW-G)/DX
          DELX=-G/DGDX
          IF((X.LE.XMIN.AND.DELX.LT.0.0E0).OR.
     +       (X.GE.XMAX.AND.DELX.GT.0.0E0)) THEN
              X0G=X
          ELSE
              XOLD=X
              X=AMAX1(XMIN,AMIN1(X+DELX,XMAX))
              IF(ITERAT.LT.MAXITERAT) THEN
                  GOTO 1
              ELSE
                  X0G=XOLD
                  ITERAT=0
              ENDIF
          ENDIF
      ELSE
          X0G=X
      ENDIF 
C
C     Calculate g(y)=0 for x=xmin:
      X=XMIN
      Y=WSGUES
      YOLD=1.0
  2   CONTINUE
      ITERAT=ITERAT+1
      G=PWVSNW(TAIR,Y,X)
     +  +CCW*STSNW(TAIR,Y,X)/ROSNW(TAIR,Y,X)
     +  -LOGPW
      IF(ABS(G).GT.TOLG.AND.ABS(YOLD-Y).GT.TOLXY) THEN
          DGDY=(PWVSNW(TAIR,Y+DY,X)
     +          +CCW*STSNW(TAIR,Y+DY,X)/ROSNW(TAIR,Y+DY,X)
     +          -LOGPW-G)/DY
          DELY=-G/DGDY
          IF((Y.LE.YMIN.AND.DELY.LT.0.0E0).OR.
     +       (Y.GE.YMAX.AND.DELY.GT.0.0E0)) THEN
              Y0G=Y
          ELSE
              YOLD=Y
              Y=AMAX1(YMIN,AMIN1(Y+DELY,YMAX))
              IF(ITERAT.LT.MAXITERAT) THEN
                  GOTO 2
              ELSE
                  Y0G=YOLD
                  ITERAT=0
              ENDIF
          ENDIF
      ELSE
          Y0G=Y
      ENDIF
C
C     Linear expression for g:
      ALFA=(Y0G-YMIN)/(XMIN-X0G)
      BETA=Y0G-ALFA*XMIN
C
      RETURN
      END
******************************************************************************
*     SUBROUTINE WGTTER(TAIR,PPWV,PPNA,WSA,WNA)
******************************************************************************
*
*     This subroutine calculates the equilibrium composition of 
*     stratospheric liquid ternary solution (sulfuric/nitric-acid,water)
*     in equilibrium with the ambient partial nitric acid and water vapor.
*
*     Source: Beyer et al.: GRL 21,871,1994
*
*     The equilibrium compositions (WSA,WNA) are calculated by solving
*             PPNA = f(TAIR,WNA,WSA) and PPWV = g(TAIR,WNA,WSA)
*     where f and g are nitric acid and water vapor saturation pressures
*     over the ternary solution. Equations are solved by Newton-Ralphson
*     iteration.
*
*     The Kelvin effect is not included.
*
*     If WSA > 0.85 the subroutine returns WSA = 0.85 and sets WNA = 0.0
*     If WNA > 0.50 the subroutine returns WNA = 0.50 and sets WSA = 0.0
*     hence by testing if WSA=0 or WNA=0 the range of validity can be
*     checked.
*
*
*
*     Input/output variables:
*     REAL(KIND=4)  TAIR,PPWV,PPNA,WSA,WNA
*
*     Input:       
*         TAIR:    K         Temperature of ambient air 
*         PPWV:    Pa        Partial pressure of ambient water vapor 
*         PPNA:    Pa        Partial pressure of ambient nitric acid vapor 
*
*     Output:
*         WSA:               Mass fraction of sulfuric acid. [< 0.85 ]
*         WNA:               Mass fraction of nitric acid. [0.0;0.50]
C
CC    SUBROUTINE WGTTER(TAIR,PPWV,PPNA,WSA,WNA)
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,PPWV,PPNA,WSA,WNA
CC    DOUBLE PRECISION A(0:2),B(0:3)
CC    DOUBLE PRECISION AN(0:2),BN(0:2),AW(0:3),BW(0:3)
CC    DATA AN/12.88E0,-4.465E0,-0.569E0/
CC    DATA BN/-2733.0E0,-685.5E0,189.9E0/
CC    DATA AW/12.09E0,-10.17E0,10.55E0,0.1187E0/
CC    DATA BW/-3202.0E0,3466.0E0,-4224.0E0,-59.74E0/
CC    DOUBLE PRECISION C1,C2,TOLF,TOLG,
CC   +        WSGUES,WNGUES,
CC   +        XMIN,XMAX,YMIN,YMAX,XYMAX
CC    PARAMETER(C1=133.32237E0,C2=0.434294482E0,
CC   +          TOLF=0.001E0,TOLG=0.001E0,
CC   +          WNGUES=0.3E0,WSGUES=0.15E0,
CC   +          XMIN=1.0E-10,XMAX=0.5E0,
CC   +          YMIN=0.01E0,YMAX=0.85E0,XYMAX=1.0E0-1.0E-10)
CC    DOUBLE PRECISION TINV,X,Y,F,G,DFDX,DFDY,DGDX,DGDY,DELX,DELY,
CC   +                 EXPY0,EXPY,Y7,EXPY7,XPLUSY,ALFA,BETA,
CC   +                 LOGPW,LOGPN
CC    INTEGER I,ITERAT,MAXITE
CC    PARAMETER(MAXITE=20)
CC    SAVE AW,BW,AN,BN
CC    LOGICAL INTERA
CC      PARAMETER(INTERA=.TRUE.)
CC    PARAMETER(INTERA=.FALSE.)
C
CC    TINV=1.0E0/TAIR
CC    DO 100 I=0,2
CC        A(I)=AN(I)+BN(I)*TINV
CC100 CONTINUE
CC    DO 200 I=0,3
CC        B(I)=AW(I)+BW(I)*TINV
CC200 CONTINUE
C
CC    LOGPW=ALOG10(PPWV/C1)
CC    LOGPN=ALOG10(PPNA/C1)
C
CC    EXPY0=EXP(1.0E0)-1.0E0
C
CC    X=WNGUES
CC    Y=WSGUES
CC    ITERAT=0
C
CC2   CONTINUE
CC    IF(INTERA) WRITE(*,'(A)') ' POINT 2'
C     Point on g-zero line:
CC    IF(X.GT.0.0E0) THEN
CC        Y7=(1.0E0-Y)**7
CC        EXPY=EXP((1.0E0-Y)*Y7)
CC        EXPY7=8.0E0*EXPY*Y7
CC        XPLUSY=X+Y
CC        G=-(ALOG10(1.0E0-XPLUSY)+B(0)+(B(2)*XPLUSY+B(1))*XPLUSY-
CC   +            B(3)*(EXPY-1.0E0)-LOGPW)
CC        IF(ABS(G).GT.TOLG) THEN
CC            DGDX=-C2/ABS(1.0E0-XPLUSY)+B(1)+2.0E0*B(2)*XPLUSY
CC            DGDY=DGDX+B(3)*EXPY7
CC            DELY=G/DGDY
CC            IF((Y.LE.YMIN.AND.DELY.LT.0.0).OR.X.GE.XMAX) THEN
CC                Y=0.0E0
CC                GOTO 3
CC            ELSE IF(Y.GE.YMAX.AND.DELY.GT.0.0) THEN
CC                X=0.0E0
CC                GOTO 2
CC            ENDIF
CC            Y=AMAX1(YMIN,AMIN1(Y+DELY,YMAX))
CC            IF(X+Y.GT.XYMAX) X=XYMAX-Y
CC            ITERAT=ITERAT+1
CC            IF(ITERAT.GT.MAXITE) THEN
CC                Y=0.0E0
CC                GOTO 3
CC            ENDIF
CC            GOTO 2
CC        ELSE
C             Linear equation for g-zero line:
CC            DGDX=-C2/ABS(1.0E0-XPLUSY)+B(1)+2.0E0*B(2)*XPLUSY
CC            DGDY=DGDX+B(3)*EXPY7
CC            BETA=(DGDX*X+DGDY*Y)/DGDY
CC            ALFA=-DGDX/DGDY
CC            IF(ALFA.LT.-1.0E0) THEN
CC                Y=0.0E0
CC                GOTO 3
CC            ENDIF
CC        ENDIF
C
CC        F=-(ALOG10(X)+A(0)+A(1)*(1.0E0-XPLUSY)-
CC   +            A(2)*(EXPY-1.0E0)-LOGPN)
C
CC        IF(INTERA) THEN
CC            WRITE(*,*) 'Point on g-zero line:'
CC            WRITE(*,*) 'X-Y',X,Y
CC            WRITE(*,*) 'G:',G
CC            WRITE(*,*) 'ALFA BETA',ALFA,BETA
CC            PAUSE
CC        ENDIF
CC    ELSE IF(X.EQ.0.0E0) THEN
C         Perform ordinary Newton Ralphson iteration with X=0:
CC4       CONTINUE
CC        IF(INTERA) THEN
CC            WRITE(*,'(A)') ' POINT 4'
CC            WRITE(*,*) 'Point on g-ITERATION:'
CC            WRITE(*,*) 'Y',X,Y
CC            WRITE(*,*) 'G:',G
CC            WRITE(*,*) 'DGDY-DELY',DGDY,DELY
CC            PAUSE
CC        ENDIF
CC        Y7=(1.0E0-Y)**7
CC        EXPY=EXP((1.0E0-Y)*Y7)
CC        EXPY7=8.0E0*EXPY*Y7
CC        G=-(ALOG10(1.0E0-Y)+B(0)+(B(2)*Y+B(1))*Y-
CC   +            B(3)*(EXPY-1.0E0)-LOGPW)
CC        IF(ABS(G).LT.TOLG.OR.(Y.GE.0.95E0.AND.DELY.GT.0.0E0)) THEN
CC            WSA=SNGL(Y)
CC            WNA=0.0
CC            RETURN
CC        ELSE
CC            DGDY=-C2/ABS(1.0E0-Y)+B(1)+2.0E0*B(2)*Y+B(3)*EXPY7
CC            DELY=G/DGDY
CC            Y=AMAX1(YMIN,AMIN1(Y+DELY,0.95E0))
CC            GOTO 4
CC        ENDIF
CC    ENDIF
C
C     Point on f-zero line:
CC3   CONTINUE
CC    IF(INTERA) WRITE(*,'(A)') ' POINT 3'
CC    IF(Y.GT.0.0E0) THEN
CC        DFDX=-A(1)+C2/ABS(X)
CC        DFDY=-A(1)+A(2)*EXPY7
CC        DFDX=DFDX+DFDY*ALFA
CC        DELX=F/DFDX
C
CC        IF(INTERA) THEN
CC            WRITE(*,*) 'New point on f-zero line:'
CC            WRITE(*,*) 'X-Y',X,Y
CC            WRITE(*,*) 'F-G',F,G
CC            WRITE(*,*) 'DELX:',DELX
CC            PAUSE
CC        ENDIF
C
CC        IF(X.GE.XMAX) THEN
CC            WNA=SNGL(X)
CC            WSA=0.0
CC            RETURN
CC        ELSE IF(X.LE.XMIN.AND.DELX.LT.0.0) THEN
CC            WNA=0.0
CC            WSA=SNGL(Y)
CC            RETURN
CC        ENDIF
CC        X=AMAX1(XMIN,AMIN1(X+DELX,XMAX))
C         Assume y to be constrained by linear relationship on g-zero line:
CC        Y=ALFA*X+BETA
C
CC        Y7=(1.0E0-Y)**7
CC        EXPY=EXP((1.0E0-Y)*Y7)
CC        EXPY7=8.0E0*EXPY*Y7
CC        XPLUSY=X+Y
CC        F=-(ALOG10(X)+A(0)+A(1)*(1.0E0-XPLUSY)-
CC   +            A(2)*(EXPY-1.0E0)-LOGPN)
CC        G=-(ALOG10(1.0E0-XPLUSY)+B(0)+(B(2)*XPLUSY+B(1))*XPLUSY-
CC   +            B(3)*(EXPY-1.0E0)-LOGPW)
CC        IF(ABS(G).GT.TOLG) THEN
C             Find new equation for g-zero line:
CC            IF(Y.GT.WSGUES) THEN
CC                Y=XYMAX-X
CC            ELSE
CC                Y=WSGUES
CC            ENDIF
CC            ITERAT=0
CC            GOTO 2
CC        ELSE IF(ABS(F).GT.TOLF) THEN
C             Proceed along g-zero line:
CC            GOTO 3
CC        ELSE
CC            WNA=SNGL(X)
CC            WSA=SNGL(Y)
CC            RETURN
CC        ENDIF
CC    ELSE IF(Y.EQ.0.0E0) THEN
C         Perform ordinary Newton Ralphson iteration with Y=0:
CC5       CONTINUE
CC        IF(INTERA) THEN
CC            WRITE(*,'(A)') ' POINT 5'
CC            WRITE(*,*) 'New point on f-ITERATION:'
CC            WRITE(*,*) 'X',X,Y
CC            WRITE(*,*) 'F',F
CC            WRITE(*,*) 'DFDX-DELX',DFDX,DELX
CC            PAUSE
CC        ENDIF
CC        F=-(ALOG10(X)+A(0)+A(1)*(1.0E0-X)-A(2)*EXPY0-LOGPN)
CC        IF(ABS(F).LT.TOLF.OR.(X.GE.1.0.AND.DELX.GT.0.0E0)) THEN
CC            WNA=SNGL(X)
CC            WSA=0.0
CC            RETURN
CC        ELSE
CC            DFDX=-A(1)+C2/ABS(X)
CC            DELX=F/DFDX
CC            X=AMAX1(XMIN,AMIN1(X+DELX,1.0E0))
CC            GOTO 5
CC        ENDIF
CC    ENDIF
C
CC    END
******************************************************************************
*     SUBROUTINE WGTTER(TAIR,PPWV,PPNA,WSA,WNA)
******************************************************************************
*
*     This subroutine calculates the equilibrium composition of 
*     stratospheric liquid ternary solution (sulfuric/nitric-acid,water)
*     in equilibrium with the ambient partial nitric acid and water vapor.
*
*     Source: Zhang et al.: J. Phys. Chem. 97, 8541, 1993.
*
*     The equilibrium compositions (WSA,WNA) are calculated by solving
*             PPNA = f(TAIR,WNA,WSA) and PPWV = g(TAIR,WNA,WSA)
*     where f and g are nitric acid and sulfuric acid saturation pressures
*     over the ternary solution. Equations are solved by Newton-Ralphson
*     iteration.
*
*     The Kelvin effect is not included.
*
*     Expressions of f and g by Zhang et al. are only valid in limited
*     composions ranges: 0.35 < WSA < 0.75, and WNA < 0.15.
*     If WSA > 0.75 the subroutine returns WSA = 0.75 and sets WNA = 0.0
*     If WNA > 0.15 the subroutine returns WNA = 0.15 and sets WSA = 0.0
*
*
*
*     Input/output variables:
*     REAL(KIND=4)  TAIR,PPWV,PPNA,WSA,WNA
*
*     Input:       
*         TAIR:    K         Temperature of ambient air 
*         PPWV:    Pa        Partial pressure of ambient water vapor 
*         PPNA:    Pa        Partial pressure of ambient nitric acid vapor 
*
*     Output:
*         WSA:               Mass fraction of sulfuric acid. [> 0.35 ]
*         WNA:               Mass fraction of nitric acid. [0.0;0.15]
C
CC    SUBROUTINE WGTTER(TAIR,PPWV,PPNA,WSA,WNA)
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,PPWV,PPNA,WSA,WNA
CC    DOUBLE PRECISION LOGPW,LOGPN
CC    DOUBLE PRECISION B(6),AW(6),BW(6)
CC    DOUBLE PRECISION A(3),AN(3),BN(3)
CC    DOUBLE PRECISION C
CC    DOUBLE PRECISION X,Y,XMIN,XMAX,XACC,YMIN,YMAX,TMIN,TMAX,TINV
CC    DOUBLE PRECISION RTSAFE
CC    EXTERNAL FUNCD1
CC    INTEGER I
CC    DATA AW/-5.10E0,-1.34E0,15.80E0,-25.98E0,-13.24E0,-0.60E0/
CC    DATA BW/3.59E+3,-0.31E+3,-7.16E+3,9.99E+3,4.85E+3,0.22E+3/
CC    DATA AN/9.65E0,2.92E0,-3.60E0/
CC    DATA BN/-3.69E+3,1.04E+3,2.42E+3/
CC    PARAMETER(TMIN=190.0E0,TMAX=230.0E0)
CC      PARAMETER(XMIN=1.0E-10,XMAX=0.50E0,XACC=1.0E-9,
CC    PARAMETER(XMIN=1.0E-10,XMAX=0.15E0,XACC=1.0E-9,
CC     +          YMIN=0.0E0,YMAX=0.75E0)
CC   +          YMIN=0.35E0,YMAX=0.75E0)
CC    PARAMETER(C=133.32237E0)
CC    SAVE AW,BW,AN,BN
CC    COMMON /SANAW/A,B,LOGPW,LOGPN,Y
C
CC    TINV=1.0E0/AMAX1(TMIN,AMIN1(TAIR,TMAX))
CC    DO 100 I=1,3
CC        A(I)=AN(I)+BN(I)*TINV
CC100 CONTINUE
CC    DO 200 I=1,6
CC        B(I)=AW(I)+BW(I)*TINV
CC200 CONTINUE
C
CC    LOGPW=ALOG10(PPWV/C)
CC    LOGPN=ALOG10(PPNA/C)
C
CC    X=RTSAFE(FUNCD1,XMIN,XMAX,XACC)
CC    IF(Y.GE.SNGL(YMAX)) THEN
CC        WSA=SNGL(YMAX)
CC        WNA=0.0
CC    ELSE
CC        IF(X.LT.XMAX) THEN
CC            WNA=SNGL(X)
CC            WSA=AMIN1(SNGL(Y),SNGL(YMAX))
CC        ELSE
CC            WNA=SNGL(XMAX)
CC            WSA=0.0
CC        ENDIF
CC    ENDIF
C
CC    RETURN
CC    END
C
CC    SUBROUTINE FUNCD1(X,G,GDOT)
CC    IMPLICIT NONE
CC    DOUBLE PRECISION X,G,GDOT
CC    DOUBLE PRECISION A(3),B(6),LOGPW,LOGPN,Y,YDOT,LOGY,C
CC    PARAMETER(C=0.434294482E0)
CC    COMMON /SANAW/A,B,LOGPW,LOGPN,Y
C
CC    Y=AMAX1((LOGPN-A(1)-ALOG10(X)-A(3)*X)/A(2),1.0E-10)
cc      Y=(LOGPN-A(1)-ALOG10(X)-A(3)*X)/A(2)
CC    YDOT=-(C/X+A(3))/A(2)
CC    LOGY=ALOG10(Y)
C
CC    G=((B(6)*LOGY+B(5))*LOGY+B(4))*LOGY+B(3)*Y+B(2)*X+B(1)-LOGPW
CC    GDOT=(((3.0E0*B(6)*LOGY+2.0E0*B(5))*LOGY+B(4))*C/Y+B(3))*YDOT+B(2)
C
C      WRITE(*,*) SNGL(A(1)+A(2)*Y+ALOG10(X)+A(3)*X-LOGPN),
C     +           SNGL(G),SNGL(X),SNGL(Y),SNGL(YDOT)
C
C
CC    RETURN
CC    END
C
******************************************************************************
*     SUBROUTINE WGTFSA(RADIUS,TAIR,PPWV,WSAS,RHOSAS,MSA)
******************************************************************************
*
*     This subroutine calculates the acid mass fraction, density, and
*     mass of sulfuric acid in a single aerosol droplet of a specified 
*     radius in equilibrium with ambient water vapor partial pressure 
*     and temperature. The Kelvin curvature effect is neglegted.
*
*     Source: Tabazadeh et al., GRL 24,1931, 1997.
*
*     Input/output variables:
*     REAL(KIND=4)  RADIUS,TAIR,PPWV,WSAS,RHOSAS,MSA
*
*     Input:       
*         RADIUS:  m         Radius of aerosol droplet
*         TAIR:    K         Temperature of ambient air 
*         PPWV:    Pa        Partial pressure of ambient water vapor 
*
*     Output:
*         WSAS:              mass fraction of sulfuric acid. [0.0;1]
*         RHOSAS:  kg/m**3   Density of sulfuric acid solution droplet
*         MSA:     kg        Mass of sulfuric acid in droplet
*
      SUBROUTINE WGTFSA(RADIUS,TAIR,PPWV,WSAS,RHOSAS,MSA)
C
      IMPLICIT NONE
      REAL(KIND=4)  RADIUS,TAIR,PPWV,WSAS,RHOSAS,MSA
C
C     Mathematical constants:
      REAL(KIND=4)  PI
      PARAMETER(PI=3.1415926536E0)
C
C     External functions needed:
      REAL(KIND=4)  PWVWL,ROSAS
C     PWVWL:       Water vapor pressure over pure liquid water
C     ROSAS:       Density of sulfuric acid solution

C     Physical constants:
      REAL(KIND=4)  MH2SO4
      PARAMETER(
C     Molar weight of sulfuric acid (kg/mole)
     +          MH2SO4=98.08E-3)
C
C     Auxiliary local variables:
      DOUBLE PRECISION AW,Y1,Y2,
     +      A11,A21,A31,A12,A22,A32,
     +      B11,B21,B31,B12,B22,B32,
     +      C11,C21,C31,C12,C22,C32,
     +      X11,X21,X31,X12,X22,X32
      SAVE  AW,Y1,Y2,
     +      A11,A21,A31,A12,A22,A32,
     +      B11,B21,B31,B12,B22,B32,
     +      C11,C21,C31,C12,C22,C32,
     +      X11,X21,X31,X12,X22,X32
      REAL(KIND=4)  C2,MS,Y
C
      PARAMETER(
     +        C2=4.0E0*PI/3.0E0)
      DATA A11,A21,A31,A12,A22,A32/
     +    12.372089320E0,11.820654354E0,-180.06541028E0,
     +    13.455394705E0,12.891938068E0,-176.95814097E0/
      DATA B11,B21,B31,B12,B22,B32/
     +    -30.490657554E0,-4.8073063730E0,-93.317846778E0,
     +    -34.285174607E0,-6.4261237757E0,-90.469744201E0/
      DATA C11,C21,C31,C12,C22,C32/
     +    -2.1133114241E0,-5.1727540348E0,273.88132245E0,
     +    -1.7620073078E0,-4.9005471319E0,267.45509988E0/
      DATA X11,X21,X31,X12,X22,X32/
     +    -0.1612551611E0,-0.20786404244E0,-0.38601102592E0,
     +    -0.1921312255E0,-0.23233847708E0,-0.36257048154E0/
C
      AW=PPWV/PWVWL(TAIR)
CC    IF(AW.LT.0.01E0) THEN
C       High sulfuric acid concentration;
C       Use Gmitro & Vermeulen vapor pressures:
CC      CALL WGTGV(RADIUS,TAIR,PPWV,WSAS,RHOSAS,MSA)
CC      RETURN
      IF(AW.GE.1.0E0) THEN
        MS=0.0
      ELSE
cc      IF(AW.GE.0.01E0.AND.AW.LT.0.05E0) THEN
        IF(AW.LT.0.05E0) THEN
            Y1=A11*AW**X11+B11*AW+C11
            Y2=A12*AW**X12+B12*AW+C12
        ELSE IF(AW.GE.0.05E0.AND.AW.LT.0.85E0) THEN
            Y1=A21*AW**X21+B21*AW+C21
            Y2=A22*AW**X22+B22*AW+C22
        ELSE IF(AW.GE.0.85E0) THEN
            Y1=A31*AW**X31+B31*AW+C31
            Y2=A32*AW**X32+B32*AW+C32
        ENDIF
        MS=Y1+(TAIR-190.0E0)*(Y2-Y1)/70.0E0
      ENDIF
      Y=MS*MH2SO4
      WSAS=Y/(1.0+Y)
      RHOSAS=ROSAS(TAIR,WSAS)
      MSA=C2*WSAS*RHOSAS*RADIUS**3
      RETURN
      END
******************************************************************************
*     SUBROUTINE WGTGV(RADIUS,TAIR,PPWV,WSAS,RHOSAS,MSA)
******************************************************************************
*
*     This subroutine calculates the acid mass fraction, density, and
*     mass of sulfuric acid in a single aerosol droplet of a specified 
*     radius in equilibrium with ambient water vapor partial pressure 
*     and temperature.
*
*     The calculation is performed by iteration of
*        ln(PPWV) - [(2Mh2o sigma)/(R T r rho) - ln(ph2osa)] = 0
*     using the secant method. Vapor pressures by Gmitro and Vermeulen
*     (PWVSAS_GV) are used.  
*
*     Input/output variables:
*     REAL(KIND=4)  RADIUS,TAIR,PPWV,WSAS,RHOSAS,MSA
*
*     Input:       
*         RADIUS:  m         Radius of aerosol droplet
*         TAIR:    K         Temperature of ambient air 
*         PPWV:    Pa        Partial pressure of ambient water vapor 
*
*     Output:
*         WSAS:              mass fraction of sulfuric acid. [0.1;1]
*         RHOSAS:  kg/m**3   Density of sulfuric acid solution droplet
*         MSA:     kg        Mass of sulfuric acid in droplet
*
      SUBROUTINE WGTGV(RADIUS,TAIR,PPWV,WSAS,RHOSAS,MSA)
C
      IMPLICIT NONE
      REAL(KIND=4)  RADIUS,TAIR,PPWV,WSAS,RHOSAS,MSA
C
C     Physical constants:
      REAL(KIND=4)  MH2O, RGAS
      PARAMETER(
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0)
C
C     Mathematical constants:
      REAL(KIND=4)  PI
      PARAMETER(PI=3.1415926536E0)
C
C     External functions needed:
      REAL(KIND=4)  PWVSAS_GV,STSAS,ROSAS
C     PWVSAS_GV:      Natural logaritm of water vapor pressure over
C                  sulfuric acid solution
C     STSAS:       Surface tension of sulfuric acid solution
C     ROSAS:       Density of sulfuric acid solution
C
C     Auxiliary local variables:
      REAL(KIND=4)  DELW,DELLP,C1,C2,W0,W1,W2,F0,F1,WGUESS,LPPWV,RO
      INTEGER ITERAT,MAXITE
      REAL(KIND=4)  WMIN
      PARAMETER(
C         Minimum H2SO4 weight fraction:
     +    WMIN=0.1E0,
C         Relative error on iterated weight fraction:
     +        DELW=0.001E0,
C         Relative error on iterated ln(pressure):
     +        DELLP=0.0001E0,
C         Guess of sulfuric acid mass fraction:
     +        WGUESS=0.7E0,
C         Maximum iteration number:
     +        MAXITE=20)
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C
      PARAMETER(
     +        C1=2.0E0*MH2O/RGAS,
     +        C2=4.0E0*PI/3.0E0)
C
C----------------------------------------------------------------------------
      W0=WGUESS
      LPPWV=ALOG(PPWV)
      RO=ROSAS(TAIR,W0)
      F0=LPPWV-C1*STSAS(TAIR,W0)/(TAIR*RADIUS*RO)-PWVSAS_GV(TAIR,W0)
      W1=W0*1.01E0
      ITERAT=0
C----------------------------------------------------------------------------
10    RO=ROSAS(TAIR,W1)
      F1=LPPWV-C1*STSAS(TAIR,W1)/(TAIR*RADIUS*RO)-PWVSAS_GV(TAIR,W1)
      IF(ABS(F1-F0).LT.DELLP) THEN
          WSAS=W1
          RHOSAS=RO
          MSA=C2*WSAS*RHOSAS*RADIUS**3
      ELSE
          W2=AMAX1(0.0E0,AMIN1((F1*W0-F0*W1)/(F1-F0),1.0E0))
          ITERAT=ITERAT+1
          IF(ABS(W2-W1).LT.DELW*ABS(W2).OR.ABS(F1).LT.DELLP.OR.
     +             ITERAT.GT.MAXITE) THEN
              WSAS=W2
              RHOSAS=RO
              MSA=C2*WSAS*RHOSAS*RADIUS**3
          ELSE
              W0=W1
              W1=W2
              F0=F1
              GOTO 10
          ENDIF
      ENDIF
      IF(WSAS.LT.WMIN) THEN
          WSAS=WMIN
          RHOSAS=ROSAS(TAIR,WMIN)
      ENDIF
C----------------------------------------------------------------------------
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING WGTFSA'
      WRITE(TEST,1010) 'RADIUS','TAIR','PPWV','WSAS','RHOSAS',
     +                 'PSAT','S','ITERAT'
      WRITE(TEST,1000) RADIUS,TAIR,PPWV,WSAS,RHOSAS,EXP(PWVSAS_GV(TAIR,
     +        WSAS)),PPWV/EXP(PWVSAS_GV(TAIR,WSAS)),ITERAT
 1000 FORMAT(10(1PE9.2,1X),0P)
 1010 FORMAT(1X,8(A9))
C----------------------------------------------------------------------------
      RETURN 
      END
******************************************************************************
*     SUBROUTINE WSASAS(TAIR,PPWV,WSAS)
******************************************************************************
*
*     This subroutine calculates the sulfuric acid mass fraction in 
*     sulfuric acid solution with a plane surface in equilibrium with 
*     the ambient water vapor partial pressure. 
*
*     Source: D.R. Hanson, A.R. Ravishankara & S. Solomon, JGR 99,3615,1994
*             Numerical fit to data Steele and Hamill (1981)
*
*     Input/output variables:
*     REAL(KIND=4)  TAIR,PPWV,WSAS
*
*     Input:       
*         TAIR:    K         Temperature of ambient air [TAIR.GT. 190 K]
*         PPWV:    Pa        Partial pressure of ambient water vapor 
*
*     Output:
*         WSAS:              mass fraction of sulfuric acid. [0.1;1]
*
      SUBROUTINE WSASAS(TAIR,PPWV,WSAS)
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PPWV,WSAS
      REAL(KIND=4) 
     +     A0,A1,A2,
     +     B0,B1,B2,
     +     LNPPWV,T
      PARAMETER(
     +     A0=-14.458E0,A1=0.19988E0,A2=0.62456E0,
     +     B0=3565.0E0,B1=44.777E0,B2=1.3204E0)
C
      LNPPWV=ALOG(PPWV/100.0E0)
      T=AMAX1(TAIR,190.0E0)
C
      WSAS=((A0+A2*LNPPWV)*T+B0)/(B1+B2*LNPPWV-A1*T)
      WSAS=WSAS/100.0E0
C
      RETURN
      END
******************************************************************************
*     SUBROUTINE RADSAS(TAIR,PPWV,MSA,RADIUS,WSAS,RHOSAS)
******************************************************************************
*
*     This subroutine calculates the radius, acid mass fraction, and density,
*     of a single aerosol droplet of a specified acid mass in equilibrium 
*     with ambient water vapor partial pressure and temperature.
*
*     The calculation is performed by iteration using the secant method.
*
*     Input/output variables:
*     REAL(KIND=4)  TAIR,PPWV,MSA,RADIUS,WSAS,RHOSAS
*
*     Input:       
*         TAIR:    K         Temperature of ambient air 
*         PPWV:    Pa        Partial pressure of ambient water vapor 
*         MSA:     kg        Mass of sulfuric acid in droplet
*                            (suggested value 10**(-18) kg)
*
*     Output:
*         RADIUS:  m         Radius of aerosol droplet
*         WSAS:              Mass fraction of sulfuric acid. [0.1;1]
*         RHOSAS:  kg/m**3   Density of sulfuric acid solution droplet
*
      SUBROUTINE RADSAS(TAIR,PPWV,MSA,RADIUS,WSAS,RHOSAS)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PPWV,MSA,RADIUS,WSAS,RHOSAS
C
C     Physical constants:
      REAL(KIND=4)  MH2O, RGAS
      PARAMETER(
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0)
C
C     Mathematical constants:
      REAL(KIND=4)  PI
      PARAMETER(PI=3.1415926536E0)
C
C     External functions needed:
      REAL(KIND=4)  PWVSAS,STSAS,ROSAS
C     PWVSAS:      Natural logaritm of water vapor pressure over
C                  sulfuric acid solution
C     STSAS:       Surface tension of sulfuric acid solution
C     ROSAS:       Density of sulfuric acid solution
C
C     Auxiliary local variables:
      REAL(KIND=4)  DELW,DELLP,C1,C2,ONETRD,W0,W1,W2,F0,F1,WGUESS,LPPWV,
     +        RO,RAD
      INTEGER ITERAT,MAXITE
      REAL(KIND=4)  WMIN
      PARAMETER(
C         Minimum H2SO4 weight fraction:
     +        WMIN=0.1E0,
C         Relative error on iterated weight fraction:
     +        DELW=0.001E0,
C         Relative error on iterated ln(pressure):
     +        DELLP=0.0001E0,
C         Guess of sulfuric acid mass fraction:
     +        WGUESS=0.7E0,
C         Maximum iteration number:
     +        MAXITE=20)
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C
      PARAMETER(
     +        C1=2.0E0*MH2O/RGAS,
     +        C2=4.0E0*PI/3.0E0,
     +        ONETRD=1.0E0/3.0E0)
C----------------------------------------------------------------------------
      W0=WGUESS
      LPPWV=ALOG(PPWV)
      RO=ROSAS(TAIR,W0)
      RAD=(MSA/(C2*RO*W0))**ONETRD
      F0=LPPWV-C1*STSAS(TAIR,W0)/(TAIR*RAD*RO)-PWVSAS(TAIR,W0)
      W1=W0*1.01E0
      ITERAT=0
C----------------------------------------------------------------------------
 10   RO=ROSAS(TAIR,W1)
      RAD=(MSA/(C2*RO*W1))**ONETRD
      F1=LPPWV-C1*STSAS(TAIR,W1)/(TAIR*RAD*RO)-PWVSAS(TAIR,W1)
      IF(F1.EQ.F0) THEN
          WSAS=W1
          RHOSAS=RO
          RADIUS=RAD
      ELSE
          W2=AMAX1(WMIN,AMIN1((F1*W0-F0*W1)/(F1-F0),1.0E0))
          ITERAT=ITERAT+1
          IF(ABS(W2-W1).LT.DELW*W2.OR.ABS(F1).LT.DELLP.OR.
     +             ITERAT.GT.MAXITE) THEN
              WSAS=W2
              RHOSAS=ROSAS(TAIR,W2)
              RADIUS=(MSA/(C2*RHOSAS*W2))**ONETRD
          ELSE
              W0=W1
              W1=W2
              F0=F1
              GOTO 10
          ENDIF
      ENDIF
C----------------------------------------------------------------------------
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING RADSAS'
      WRITE(TEST,1010) 'RADIUS','TAIR','PPWV','WSAS','RHOSAS',
     +                 'PSAT','S','ITERAT'
      WRITE(TEST,1000) RADIUS,TAIR,PPWV,WSAS,RHOSAS,EXP(PWVSAS(TAIR,
     +        WSAS)),PPWV/EXP(PWVSAS(TAIR,WSAS)),ITERAT
 1000 FORMAT(10(1PE9.2,1X),0P)
 1010 FORMAT(1X,8(A9))
      RETURN 
      END
******************************************************************************
*     SUBROUTINE WGTFNA(RADIUS,TAIR,PPWV,WNAS,RHONAS,MNA)
******************************************************************************
*
*     This subroutine calculates the acid mass fraction, density, and
*     mass of nitric acid in a single aerosol droplet of a specified 
*     radius in equilibrium with ambient water vapor partial pressure 
*     and temperature. The Kelvin curvature effect is neglegted.
*
*     Source: Tabazadeh et al., GRL 24,2007, 1997.
*
*     Input/output variables:
*     REAL(KIND=4)  RADIUS,TAIR,PPWV,WNAS,RHONAS,MNA
*
*     Input:       
*         RADIUS:  m         Radius of aerosol droplet
*         TAIR:    K         Temperature of ambient air 
*         PPWV:    Pa        Partial pressure of ambient water vapor 
*
*     Output:
*         WNAS:              Mass fraction of nitric acid. [0.0;1]
*         RHONAS:  kg/m**3   Density of nitric acid solution droplet
*         MNA:     kg        Mass of nitric acid in droplet
*
      SUBROUTINE WGTFNA(RADIUS,TAIR,PPWV,WNAS,RHONAS,MNA)
C
      IMPLICIT NONE
      REAL(KIND=4)  RADIUS,TAIR,PPWV,WNAS,RHONAS,MNA
C
C     Mathematical constants:
      REAL(KIND=4)  PI
      PARAMETER(PI=3.1415926536)
C
C     External functions needed:
      REAL(KIND=4)  PWVWL,RONAS
C     PWVWL:       Water vapor pressure over pure liquid water
C     RONAS:       Density of nitric acid solution

C     Physical constants:
      REAL(KIND=4)  MHNO3
      PARAMETER(
C       Molar weight of nitric acid (kg/mole)
     +          MHNO3=63.01E-3)
C
C     Auxiliary local variables:
      DOUBLE PRECISION AW,Y1,Y2,
     +      A11,A21,A31,A12,A22,A32,
     +      B11,B21,B31,B12,B22,B32,
     +      C11,C21,C31,C12,C22,C32,
     +      X11,X21,X31,X12,X22,X32
      SAVE  AW,Y1,Y2,
     +      A11,A21,A31,A12,A22,A32,
     +      B11,B21,B31,B12,B22,B32,
     +      C11,C21,C31,C12,C22,C32,
     +      X11,X21,X31,X12,X22,X32
      REAL(KIND=4)  C2,MS,Y
C
      PARAMETER(
     +        C2=4.0E0*PI/3.0E0)
      DATA A11,A21,A31,A12,A22,A32/
     +    5.6811966385E0,7.9415477295E0,-208.46136893E0,
     +    9.0405763286E0,9.5033394917E0,-224.06176695E0/
      DATA B11,B21,B31,B12,B22,B32/
     +    178.44655060E0,-10.675959795E0,-59.746717399E0,
     +    5.3246079021E1,-13.366846037E0,-31.163742477E0/
      DATA C11,C21,C31,C12,C22,C32/
     +    -1.9944014406E0,5.7733271636E0,268.48564500E0,
     +    4.6674594268E1,4.5796103358E0,255.42273520E0/
      DATA X11,X21,X31,X12,X22,X32/
     +    -0.52781515793E0,-0.40853067200E0,-0.13840491342E0,
     +    -0.54424836165E0,-0.53655348172E0,-0.030723152773E0/
C
      AW=PPWV/PWVWL(TAIR)
      IF(AW.GE.1.0E0) THEN
        MS=0.0
      ELSE
        IF(AW.LT.0.05E0) THEN
            Y1=A11*AW**X11+B11*AW+C11
            Y2=A12*AW**X12+B12*AW+C12
        ELSE IF(AW.GE.0.05E0.AND.AW.LT.0.6E0) THEN
            Y1=A21*AW**X21+B21*AW+C21
            Y2=A22*AW**X22+B22*AW+C22
        ELSE IF(AW.GE.0.6E0) THEN
            Y1=A31*AW**X31+B31*AW+C31
            Y2=A32*AW**X32+B32*AW+C32
        ENDIF
        MS=Y1+(TAIR-190.0E0)*(Y2-Y1)/70.0E0
      ENDIF
      Y=MS*MHNO3
      WNAS=Y/(1.0E0+Y)
      RHONAS=RONAS(TAIR,WNAS)
      MNA=C2*WNAS*RHONAS*RADIUS**3
      RETURN
      END
******************************************************************************
*     REAL FUNCTION HOMFZR(SOLID,TAIR,PSOLID,PLIQID,DIFACT,RHOSLD,RHOLQD,ST)
******************************************************************************
*
*     This subroutine calculates the homogeneous nucleation freezing rate.
*
*
*     Declaration of input/output variables:
*
*     INTEGER SOLID
*     REAL(KIND=4)  TAIR,PSOLID,PLIQID,RHOSLD,RHOLQD,ST,DIFACT,FZRATE
*
*     Input: 
*            
*         SOLID:                  Solid indicatior number
*                                 (Ice: SOLID=1, SAT: SOLID=2)
*         TAIR:    (K)            Ambient air temperature
*         PSOLID:  (ln(Pa))       Natural logaritm of saturation pressure 
*                                     over solid phase
*         PLIQID:  (ln(Pa))       Natural logaritm of saturation pressure 
*                                     over liquid phase
*         DIFACT:  (J)            Activation energy for self-diffusion
*         RHOSLD:  (kg/m**3)      Solid phase density
*         RHOLQD:  (kg/m**3)      Liquid phase density
*         ST:      (J/m**2)       Surface tension solid/liquid
*
*     Output:
*
*         HOMFZR:  (1/m**3 s)     Homogeneous freezing rate
*
C
      REAL FUNCTION HOMFZR(SOLID,TAIR,PSOLID,PLIQID,
     +                     DIFACT,RHOSLD,RHOLQD,ST)
C
      IMPLICIT NONE
      INTEGER SOLID
      REAL(KIND=4)  TAIR,PSOLID,PLIQID,DIFACT,RHOSLD,RHOLQD,ST
C----------------------------------------------------------------------------
C     Mathematical constants:
      DOUBLE PRECISION PI
      PARAMETER(PI=3.1415926536E0)
C
C     Physical constants:
      DOUBLE PRECISION KBOLTZ,RGAS,H,MH2SO4,MH2O
      PARAMETER(
C       Boltzmanns constant (J/K)
     +          KBOLTZ=1.380662E-23,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0,
C       Planck's constant (J s)
     +          H=6.626176E-34,
C       Molar weight of sulfuric acid (kg/mole)
     +          MH2SO4=98.08E-3,
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3)
C----------------------------------------------------------------------------
C     Minimum saturation:
      REAL(KIND=4)  SMIN
      PARAMETER(SMIN=1.0E-20)
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      DOUBLE PRECISION PREFAC,X,C1,C2,C3,C1SA,C2SA,C3SA,C1W,C2W,C3W,DG
      PARAMETER(
     +    C1SA=16.0E0*PI*MH2SO4*MH2SO4/(3.0E0*RGAS*RGAS),
     +    C2SA=RGAS/(MH2SO4*H),
     +    C3SA=2.0E0*MH2SO4/RGAS,
     +    C1W=16.0E0*PI*MH2O*MH2O/(3.0E0*RGAS*RGAS),
     +    C2W=RGAS/(MH2O*H),
     +    C3W=2.0E0*MH2O/RGAS)
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C----------------------------------------------------------------------------
      IF(PLIQID-PSOLID.LE.SMIN) THEN
          HOMFZR=0.0E0
      ELSE
          IF(SOLID.EQ.1) THEN
              C1=C1W
              C2=C2W
          ELSE
              C1=C1SA
              C2=C2SA
          ENDIF
          DG=C1*(ST**3)/((RHOSLD*TAIR*(PLIQID-PSOLID))**2)
          PREFAC=C2*RHOLQD*TAIR
          X=-(DIFACT+DG)/(KBOLTZ*TAIR)
          HOMFZR=PREFAC*EXP(X)
      ENDIF
C----------------------------------------------------------------------------
      IF(TEST.LT.60) RETURN
C
      IF(SOLID.EQ.1) THEN
          C3=C3W
      ELSE
          C3=C3SA
      ENDIF
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING HOMFZR'
      WRITE(TEST,1010) 'SOLID','TAIR','RHOSLD','RHOLQD','ST','PSOLID',
     +                 'PLIQID','S'
      WRITE(TEST,1000) SOLID,TAIR,RHOSLD,RHOLQD,ST,PSOLID,
     +                 PLIQID,PLIQID-PSOLID
      WRITE(TEST,1010) 'DG','DIFACT','GERMRAD','PREFAC','X','HOMFZR'
      WRITE(TEST,1000) DG,DIFACT,C3*ST/(RHOSLD*TAIR*(PLIQID-PSOLID)),
     +                  PREFAC,X,HOMFZR
 1000 FORMAT(10(1PE9.2,1X),0P)
 1010 FORMAT(1X,8(A9))
      RETURN
      END
******************************************************************************
*     REAL FUNCTION FRZICE(TAIR,DIFACT,RHOSLD,RHOLQD,ST,AW)
******************************************************************************
*
*     This subroutine calculates the homogeneous nucleation freezing rate
*     of ice in supercooled water or solution.
*     Source: Pruppacher & Klett, 1978, p. 178 and 152.
*     Source: Tabazadeh et al., JGR 102,23845,1997
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,DIFACT,RHOSLD,RHOLQD,ST,AW
*
*     Input: 
*            
*         TAIR:    (K)            Ambient air temperature
*         DIFACT:  (J)            Diffusion activation energy
*         RHOSLD:  (kg/m**3)      Solid phase density
*         RHOLQD:  (kg/m**3)      Liquid phase density
*         ST:      (J/m**2)       Surface tension solid/liquid
*         AW:                     Activity of water (=1 for pure water)
*     Output:
*
*         FRZICE:  (1/m**3 s)     Homogeneous freezing rate
*
C
      REAL FUNCTION FRZICE(TAIR,DIFACT,RHOSLD,RHOLQD,ST,AW)
      IMPLICIT NONE
C
      REAL(KIND=4)  TAIR,DIFACT,RHOSLD,RHOLQD,ST,AW
C----------------------------------------------------------------------------
C     Mathematical constants:
      DOUBLE PRECISION PI
      PARAMETER(PI=3.1415926536E0)
C
C     Physical constants:
      DOUBLE PRECISION KBOLTZ,RGAS,H,MH2O,NC,T0
      PARAMETER(
C       Boltzmanns constant (J/K)
     +          KBOLTZ=1.380662D-23,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0,
C       Planck's constant (J s)
     +          H=6.626176D-34,
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153D-3,
C       Prefactor (1/m**2):
c       Source: Pruppacher and Klett (1980), p. 178:
     +          NC=5.85D+18,
C       Melting temperature (K):
     +          T0=273.15D0)
C----------------------------------------------------------------------------
C     Minimum saturation:
      REAL(KIND=4)  SMIN
      PARAMETER(SMIN=1.0E-20)
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      DOUBLE PRECISION PREFAC,X,C3
      DOUBLE PRECISION C1,C2,DG,AG,KT,ROSLDM,LMELTM,TM,RHOL,RHOS,T,SIG,A
      DOUBLE PRECISION R0,R1,R2,R3,L0,L1,L2,L3,T1,T2,T3
C     Fitting coefficients for ice density and latent heat of melting,
c     Source: Pruppacher and Klett (1980), p. 86, 89:
      PARAMETER(
     +    R0=1.0E+3,R1=0.9167D0,R2=-1.75D-4/2.D0,R3=-5.0D-7/3.0D0,
     +    L0=4.1868D+3*MH2O,L1=79.7D0,L2=0.485D0/2.D0,L3=-2.5D-3/3.D0)
      PARAMETER(
     +    C1=2.0D0*MH2O,
     +    C2=4.0D0*PI/3.0D0,
     +    C3=2.0D0*NC/H)
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C----------------------------------------------------------------------------
      IF(T0-TAIR.LE.SMIN) THEN
          FRZICE=0.0E0
          RETURN
      ELSE
          RHOS=RHOSLD
          RHOL=RHOLQD
          SIG=ST
          A=AW
          T=TAIR
          T1=T-T0
          T2=T1*T1
          T3=T2*T1
          TM=(T+T0)/2.0D0
          ROSLDM =R0*(R1*T1+R2*T2+R3*T3)/T1
          LMELTM =L0*(L1*T1+L2*T2+L3*T3)/T1
          AG=ROSLDM*(LMELTM*DLOG(T0/T)+
     +       RGAS*TM*DLOG(DMAX1(1.0D-3,DMIN1(A,1.0D0))))
          IF(AG.LE.0.0E0) THEN
              FRZICE=0.0E0
              RETURN
          ENDIF
          AG=C1*SIG/AG
          DG=C2*SIG*AG*AG
          KT=KBOLTZ*T
          PREFAC=C3*DSQRT(SIG*KT)*RHOL/RHOS
          X=-(DIFACT+DG)/KT
          X=DMAX1(-100.0D0,X)
          FRZICE=PREFAC*EXP(X)
      ENDIF
C----------------------------------------------------------------------------
      IF(TEST.LT.60) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING FRZICE'
      WRITE(TEST,1010) 'TAIR','RHOSLD','RHOLQD','ST','AW'
      WRITE(TEST,1000) TAIR,RHOSLD,RHOLQD,ST,AW
      WRITE(TEST,1010) 'DG','DIFACT','GERMRAD','PREFAC','X','FRZICE'
      WRITE(TEST,1000) DG,DIFACT,AG,PREFAC,X,FRZICE
 1000 FORMAT(10(1PE9.2,1X),0P)
 1010 FORMAT(1X,8(A9))
      RETURN
      END
******************************************************************************
*     SUBROUTINE HOMFRZ_WL_IC(TAIR,J)
******************************************************************************
*
*     Subroutine used to calculate the homogeneous freezing rate of
*     water ice in supercooled pure water.
*
*     Source: Jensen et al., JGR 99, 10421, 1994
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,J
*
*     Input: 
*         TAIR:    (K)            Ambient air temperature
*
*     Output:
*         J:       (1/m**3 s)     Homogeneous freezing rate
*
      SUBROUTINE HOMFRZ_WL_IC(TAIR,J)
      IMPLICIT NONE
C
      REAL(KIND=4)  TAIR,J
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      REAL(KIND=4)  DIFACT,RHOSLD,RHOLQD,ST,AW
C
C     Thermodynamic functions:
      REAL(KIND=4)  DIFEWL,ROICE,ROWL,STWLIC,FRZICE
C
      DIFACT=DIFEWL(TAIR)
      RHOSLD=ROICE(TAIR)
      RHOLQD=ROWL(TAIR)
      ST=STWLIC(TAIR)
      AW=1.0E0
      J=FRZICE(TAIR,DIFACT,RHOSLD,RHOLQD,ST,AW)
      RETURN
      END
******************************************************************************
*     SUBROUTINE HOMFRZ_SA_IC(TAIR,WSA,J)
******************************************************************************
*
*     Subroutine used to calculate the homogeneous freezing rate of
*     water ice in supercooled sulfuric acid solution.
*
*     Source: Tabazadeh et al., JGR 102,23845,1997
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,WSA,J
*
*     Input: 
*         TAIR:    (K)            Ambient air temperature
*         WSA:                    Mass fraction of sulfuric acid. [0.1;1]
*         PPWV:    (Pa)           Partial pressure of ambient water vapor
*
*     Output:
*         J:       (1/m**3 s)     Homogeneous freezing rate
*
      SUBROUTINE HOMFRZ_SA_IC(TAIR,WSA,PPWV,J)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WSA,PPWV,J
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      REAL(KIND=4)  DIFACT,RHOSLD,RHOLQD,ST,AW
C
C     Thermodynamic functions:
      REAL(KIND=4)  DIFEAT,ROICE,ROSAS,STSAIC,FRZICE,PWVWL
C----------------------------------------------------------------------------
      ST=STSAIC(TAIR,WSA)
C     Diffusion activation energy from Tabazadeh et al. GRL submitted 1999:
      DIFACT=DIFEAT(TAIR)
      RHOSLD=ROICE(TAIR)
      RHOLQD=ROSAS(TAIR,WSA)
      AW=PPWV/PWVWL(TAIR)
      J=FRZICE(TAIR,DIFACT,RHOSLD,RHOLQD,ST,AW)
C
      RETURN
      END
******************************************************************************
*     SUBROUTINE HOMFRZ_SNW_IC(TAIR,PAIR,PPWV,J)
******************************************************************************
*
*     Subroutine used to calculate the homogeneous freezing rate of
*     water ice in supercooled ternary solution.
*
*     Source: Koop et al., Nature 406, 611-614, 2000
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,PAIR,PPWV,J
*
*     Input: 
*         TAIR:    (K)            Ambient air temperature
*         PAIR:    (Pa)           Ambient air pressure
*         PPWV:    (Pa)           Partial pressure of ambient water vapor
*
*     Output:
*         J:       (1/m**3 s)     Homogeneous freezing rate
*
      SUBROUTINE HOMFRZ_SNW_IC(TAIR,PAIR,PPWV,J)
      REAL TAIR,PAIR,PPWV,J,LOGJ,PWVWL
      REAL(KIND=8) DAW,AW,AWI,DMY,VI,VW0,LNT,RGAS,KT0,DKT0DP,
     +       KTI,DKTIDP,P,IV,T,P2,P3
      PARAMETER(
     +          KT0=1.6D0,
     +          DKT0DP=-8.8D0,
     +          KTI=0.22D0,
     +          DKTIDP=-0.17D0,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441D0)
      T=DBLE(TAIR)
      LNT=DLOG(T)
      P=1.0D-9*DBLE(PAIR)
      P2=P*P
      P3=P2*P/6.0D0
      VW0=-230.76D0-0.1478D0*T+(4099.2D0/T)+48.8341D0*LNT
      VI=19.43D0-2.2D-3*T+1.08D-5*T*T
      IV=VW0*(P-0.5D0*KT0*P2-DKT0DP*P3)-VI*(P-0.5D0*KTI*P2-DKTIDP*P3)
      IV=IV*1.0D-3
      DMY=210368.0D0+131.438D0*T-(3.32373D6/T)-41729.1D0*LNT
      AWI=DEXP(DMY/(RGAS*T))
      AW=DBLE(PPWV)/DBLE(PWVWL(TAIR))
      DAW=AW*DEXP(IV/(RGAS*T))-AWI
      LOGJ=REAL(((29180.0D0*DAW-26924.0D0)*DAW+8502.0D0)*DAW-906.7D0)
      J=1.0E6*(10.0E0**(LOGJ))
      RETURN
      END SUBROUTINE
******************************************************************************
*     SUBROUTINE HOMFRZ_SNW_IC(TAIR,PPWV,J)
******************************************************************************
*
*     Subroutine used to calculate the homogeneous freezing rate of
*     water ice in supercooled ternary solution.
*     PRELIMINARY homogeneour freezing model for STS, assuming the
*     the freezing rate at a specific temperature and water vapor partial
*     to be equal to homogeneous freezing of sulfuric acid solution
*     at the same conditions.
*
*     Source: Tabazadeh et al., GRL 27, 1111, 2000
*     Source: Chang et al. J. Phys. Chem. 103, 2673, 1999
*
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,WSA,J
*
*     Input: 
*         TAIR:    (K)            Ambient air temperature
*         PPWV:    (Pa)           Partial pressure of ambient water vapor
*
*     Output:
*         J:       (1/m**3 s)     Homogeneous freezing rate
*
CC    SUBROUTINE HOMFRZ_SNW_IC(TAIR,PPWV,J)
C
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,PPWV,J
C----------------------------------------------------------------------------
C     Auxiliary local variables:
CC    REAL(KIND=4)  DIFACT,RHOSLD,RHOLQD,ST,AW,WSA,MSA,RHOSAS
C
C     Thermodynamic functions:
CC    REAL(KIND=4)  DIFEAT,ROICE,STSAIC,FRZICE,PWVWL
CC    REAL STICE
CC    PARAMETER(
CC   +          STICE=105.0E-3)
C----------------------------------------------------------------------------
CC    CALL WGTFSA(1.0E-6,TAIR,PPWV,WSA,RHOSAS,MSA)
CC    ST=STSAIC(TAIR,WSA)
C     Diffusion activation energy from Tabazadeh et al. GRL submitted 1999:
CC    DIFACT=DIFEAT(TAIR)
CC    RHOSLD=ROICE(TAIR)
CC    RHOLQD=RHOSAS
CC    AW=PPWV/PWVWL(TAIR)
CC    J=FRZICE(TAIR,DIFACT,RHOSLD,RHOLQD,ST,AW)
C
CC    RETURN
CC    END
******************************************************************************
*     SUBROUTINE HOMFRZ_SNW_NAT(TAIR,PPWV,WSA,WNA,J)
******************************************************************************
*
*     Subroutine used to calculate the homogeneous freezing rate of
*     nitric acid trihydrate (NAT) in supercooled ternary solution.
*
*     Source: Tabazadeh et al., Science  291, 2591,2001
*     Source: Salcedo et al., J. Phys. Chem. 105, 1433, 2001
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,PPWV,WSA,WNA,J
*
*     Input: 
*         TAIR:    (K)            Ambient air temperature
*         PPWV:    (Pa)           Partial pressure of ambient water vapor
*         WSA:                    Weight fraction of H2SO4  ]0.0; 0.7[
*         WNA:                    Weight fraction of HNO3   ]0;0.7[
*
*     Output:
*         J:       (1/m**3 s)     Homogeneous freezing rate
*
      SUBROUTINE HOMFRZ_SNW_NAT(TAIR,PPWV,WSA,WNA,J)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PPWV,WSA,WNA,J
C----------------------------------------------------------------------------
      REAL(KIND=4) RGAS
      PARAMETER(
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0)
C     Auxiliary local variables:
      REAL(KIND=8) C,PRE,X
      REAL(KIND=4) G,S
      PARAMETER(C=2.2128D39)
C
C     Thermodynamic functions:
      REAL(KIND=4)  PNASNW,PNANAT
C----------------------------------------------------------------------------
      PRE=C*DSQRT(DBLE(TAIR))
      S=EXP(PNASNW(TAIR,WSA,WNA))/PNANAT(TAIR,PPWV)
      IF(S.LT.1.0) THEN
          J=0.0
      ELSE
          G=(30.9E0-S*0.14E0)*4.1868E3
          X=PRE*DEXP(DBLE(-G/(RGAS*TAIR)))
          J=REAL(X)
      ENDIF
C
      RETURN
      END
******************************************************************************
*     SUBROUTINE HOMFRZ_SNW_NAD(TAIR,PPWV,WSA,WNA,J)
******************************************************************************
*
*     Subroutine used to calculate the homogeneous freezing rate of
*     nitric acid dihydrate (NAD) in supercooled ternary solution.
*
*     Source: Tabazadeh et al., Science 291, 2591, 2001
*     Source: Salcedo et al., J. Phys. Chem. 105, 1433, 2001
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,PPWV,WSA,WNA,J
*
*     Input: 
*         TAIR:    (K)            Ambient air temperature
*         PPWV:    (Pa)           Partial pressure of ambient water vapor
*         WSA:                    Weight fraction of H2SO4  ]0.0; 0.7[
*         WNA:                    Weight fraction of HNO3   ]0;0.7[
*
*     Output:
*         J:       (1/m**3 s)     Homogeneous freezing rate
*
      SUBROUTINE HOMFRZ_SNW_NAD(TAIR,PPWV,WSA,WNA,J)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PPWV,WSA,WNA,J
C----------------------------------------------------------------------------
      REAL(KIND=4) RGAS
      PARAMETER(
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0)
C     Auxiliary local variables:
      REAL(KIND=8) C,PRE,X
      REAL(KIND=4) G,S
      PARAMETER(C=2.7168D39)
C
C     Thermodynamic functions:
      REAL(KIND=4)  PNASNW,PNANAD
C----------------------------------------------------------------------------
      PRE=C*DSQRT(DBLE(TAIR))
      S=EXP(PNASNW(TAIR,WSA,WNA))/PNANAD(TAIR,PPWV)
      IF(S.LT.1.0) THEN
          J=0.0
      ELSE
          G=(28.8E0-S*0.37E0)*4.1868E3
          X=PRE*DEXP(DBLE(-G/(RGAS*TAIR)))
          J=REAL(X)
      ENDIF
C
      RETURN
      END
******************************************************************************
*     SUBROUTINE VOL_FRZ_SNW_NAD(TAIR,PPWV,PPNA,WSA,WNA,JNAD)
******************************************************************************
*
*     Subroutine used to calculate the volume freezing rate of
*     nitric acid dihydrate (NAD) in supercooled ternary solution.
*
*     Source: Möhler et al., ACPD, 2006
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,WSA,WNA,J
*
*     Input: 
*         TAIR:    (K)            Ambient air temperature
*         PPWV:    (Pa)           Partial pressure of water vapor
*         PPNA:    (Pa)           Partial pressure of nitric acid
*         WSA:                    Weight fraction of H2SO4  ]0.0; 0.7[
*         WNA:                    Weight fraction of HNO3   ]0;0.7[
*
*     Output:
*         JNAD:    (1/m**3 s)     Homogeneous freezing rate
*
      SUBROUTINE VOL_FRZ_SNW_NAD(TAIR,PPWV,PPNA,WSA,WNA,JNAD)
C
C
      IMPLICIT NONE
      INTEGER SOLID
      REAL(KIND=4)  TAIR,PPWV,PPNA,WSA,WNA,JNAD
C----------------------------------------------------------------------------
C     Physical constants:
      DOUBLE PRECISION RGAS,H,MHNO3,E
      PARAMETER(
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441D0,
C       Planck's constant (J s)
     +          H=6.626176D-34,
C       Molar weight of nitric acid (kg/mole)
     +          MHNO3=63.01D-3,
C       Conversion factor kcal to J:
     +          E=4.1868D3)
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      DOUBLE PRECISION PREFAC,C,A,B,T,G,S
      PARAMETER(
     +    C=RGAS/(MHNO3*H),
     +    A=2.5D6*E)
C     Thermodynamic functions
      REAL ROSNW,PNANAD
      B(T)=E*(11.2D0-0.1D0*(T-192.0D0))
C----------------------------------------------------------------------------
      S=PPNA/PNANAD(TAIR,PPWV)
      IF(S.LT.1.0D0) THEN
          JNAD=0.0
      ELSE
          PREFAC=C*WNA*ROSNW(TAIR,WSA,WNA)*TAIR
          G=B(DBLE(TAIR))+(A/((DBLE(TAIR)*DLOG(S))**2))
          JNAD=PREFAC*DEXP(-G/(RGAS*TAIR))
      ENDIF
      RETURN
      END
******************************************************************************
*     SUBROUTINE HOMFRZ_NA_NAT(TAIR,WNA,J)
******************************************************************************
*
*     Subroutine used to calculate the homogeneous freezing rate of
*     nitric acid trihydrate (NAT) in supercooled nitric acid solution.
*
*     Source: Salcedo et al., J. Phys. Chem. 105, 1433, 2001
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,WNA,J
*
*     Input: 
*         TAIR:    (K)            Ambient air temperature [160;180K]
*         WNA:                    Weight fraction of HNO3 [0.43; 0.54]
*
*     Output:
*         J:       (1/m**3 s)     Homogeneous freezing rate
*
      SUBROUTINE HOMFRZ_NA_NAT(TAIR,WNA,J)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WNA,J
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      REAL(KIND=8) T,W,JNAT,D(0:4),E(0:4,0:4)
      INTEGER I
C
      DATA E/
     +14917832.065D0,-1242036.2312D0,38606.21085D0,-530.50600968D0,
     +                                           2.7207720872D0,
     +-347829.46943D0,28981.349854D0,-901.79202539D0,12.40911321D0,
     +                                          -0.0637435797D0,
     +3043.506491D0,-253.73474541D0,7.9021536695D0,-0.1088619546D0,
     +                                           5.5994878573D-4,
     +-11.843155128D0,0.9877929804D0,-3.0784815843D-2,4.2450030293D-4,
     +                                           -2.185910865D-6,
     +0.0172914141D0,-1.4426745735D-3,4.4986494412D-5,-6.2081624349D-7,
     +                                           3.19979661D-9/
      SAVE D
      T=DBLE(TAIR)
      W=DBLE(AMIN1(WNA*100.0E0,54.0E0))
C      W=DBLE(WNA*100.0E0)

      IF(T.LT.160.0D0.OR.T.GT.180.0D0.OR.W.LT.43.0) THEN
          J=0.0
      ELSE
          DO I=0,4
              D(I)=(((E(4,I)*W+E(3,I))*W+E(2,I))*W+E(1,I))*W+E(0,I)
          END DO
          JNAT=(((D(4)*T+D(3))*T+D(2))*T+D(1))*T+D(0)
          J=EXP(REAL(JNAT))*1.0E6
      ENDIF
C----------------------------------------------------------------------------
C
      RETURN
      END
******************************************************************************
*     SUBROUTINE HOMFRZ_NA_NAD(TAIR,WNA,J)
******************************************************************************
*
*     Subroutine used to calculate the homogeneous freezing rate of
*     nitric acid dihydrate (NAD) in supercooled nitric acid solution.
*
*     Source: Salcedo et al., J. Phys. Chem. 105, 1433, 2001
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,WNA,J
*
*     Input: 
*         TAIR:    (K)            Ambient air temperature [170;200K]
*         WNA:                    Weight fraction of HNO3 [0.50; 0.64]
*
*     Output:
*         J:       (1/m**3 s)     Homogeneous freezing rate
*
      SUBROUTINE HOMFRZ_NA_NAD(TAIR,WNA,J)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WNA,J
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      REAL(KIND=8) T,W,JNAD,A(0:4),B(0:4,0:4)
      INTEGER I
C
      DATA B/
     +10945616.89D0,-781356.44985D0,20791.74099D0,-244.28026567D0,
     +                                           1.070082528D0,
     +-234214.93924D0,16745.051036D0,-446.41516646D0,5.2562745299D0,
     +                                          -0.0230795817D0,
     +1880.29960210D0,-134.60840154D0,3.5944181746D0,-0.042401854878D0,
     +                                           1.865602733D-4,
     +-6.7114145235D0,0.48101016303D0,-0.012862469834D0,1.5198294564D-4,
     +                                           -6.69888144D-7,
     +8.985811625D-3,-6.4464816758D-4,1.725951771D-5,-2.0423328593D-7,
     +                                           9.016014451D-10/
      SAVE B
      T=DBLE(TAIR)
      W=DBLE(AMIN1(WNA*100.0E0,64.0E0))
C      W=DBLE(WNA*100.0E0)

      IF(T.LT.170.0D0.OR.T.GT.200.0D0.OR.W.LT.50.0) THEN
          J=0.0
      ELSE
          DO I=0,4
              A(I)=(((B(4,I)*W+B(3,I))*W+B(2,I))*W+B(1,I))*W+B(0,I)
          END DO
          JNAD=(((A(4)*T+A(3))*T+A(2))*T+A(1))*T+A(0)
          J=EXP(REAL(JNAD))*1.0E6
      ENDIF
C----------------------------------------------------------------------------
C
      RETURN
      END
******************************************************************************
*     SUBROUTINE HOMFRZ_NA_IC(TAIR,WNA,PPWV,J)
******************************************************************************
*
*     Subroutine used to calculate the homogeneous freezing rate of
*     water ice in supercooled nitric acid solution.
*
*     Source: Tabazadeh et al., GRL 24, 2007, 1997
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,WSA,J
*
*     Input: 
*         TAIR:    (K)        Ambient air temperature
*         WNA:                Mass fraction of nitric acid. [0.1;1]
*         PPWV:    (Pa)       Partial pressure of ambient water vapor
*
*     Output:
*         J:       (1/m**3 s) Homogeneous freezing rate
*
      SUBROUTINE HOMFRZ_NA_IC(TAIR,WNA,PPWV,J)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WNA,PPWV,J
C----------------------------------------------------------------------------
C     Auxiliary local variables:
      REAL(KIND=4)  DIFACT,RHOSLD,RHOLQD,ST,AW
C
C     Thermodynamic functions:
      REAL(KIND=4)  DIFEAT,ROICE,RONAS,STNAIC,FRZICE,PWVWL
C----------------------------------------------------------------------------
      ST=STNAIC(TAIR,WNA)
C     Diffusion activation energy from Tabazadeh et al. JGR 102,23845,1997:
      DIFACT=DIFEAT(TAIR)
      RHOSLD=ROICE(TAIR)
      RHOLQD=RONAS(TAIR,WNA)
      AW=PPWV/PWVWL(TAIR)
      J=FRZICE(TAIR,DIFACT,RHOSLD,RHOLQD,ST,AW)
C
      RETURN
      END
******************************************************************************
*     SUBROUTINE SURF_FRZ_SNW_NAD(TAIR,WSA,WNA,JNAD)
******************************************************************************
*
*     Subroutine used to calculate the surface freezing rate of
*     nitric acid dihydrate (NAD) in supercooled ternary solution.
*
*     Source: Tabazadeh et al., J. Phys. Chem. A. 106, 10238, 2002
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,WSA,WNA,J
*
*     Input: 
*         TAIR:    (K)            Ambient air temperature
*         WSA:                    Weight fraction of H2SO4  ]0.0; 0.7[
*         WNA:                    Weight fraction of HNO3   ]0;0.7[
*
*     Output:
*         JNAD:    (1/m**2 s)     Homogeneous freezing rate
*
      SUBROUTINE SURF_FRZ_SNW_NAD(TAIR,WSA,WNA,JNAD)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WSA,WNA,JNAD
C----------------------------------------------------------------------------
      REAL(KIND=8)NS,KBOLTZ,RGAS,H
      PARAMETER(
C       Number density of molecular surface sites (1/m**2)
     +          NS=1.0D19,
C       Boltzmanns constant (J/K)
     +          KBOLTZ=1.380662D-23,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441D0,
C       Planck's constant (J s)
     +          H=6.626176D-34)
      REAL(KIND=4) MH2SO4,MHNO3,MH2O,XX
      PARAMETER(
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
C       Molar weight of sulfuric acid (kg/mole)
     +          MH2SO4=98.08E-3,
C       Molar weight of nitric acid (kg/mole)
     +          MHNO3=63.01E-3)
C     Auxiliary local variables:
      REAL(KIND=8) C,PRE,X,G,T
      PARAMETER(C=NS*KBOLTZ/H)
C----------------------------------------------------------------------------
      XX=(WSA/MH2SO4)+(WNA/MHNO3)+((1.0-WSA-WNA)/MH2O)
      X=DBLE((WNA/MHNO3)/XX)
c      X=X*0.82D0
      T=DBLE(TAIR)
      PRE=C*X*T
      G=4.1868D3*(11.5593D0+0.0804214D0*T-(71.5133D0-0.256724D0*T)*X)
      X=PRE*DEXP(-G/(RGAS*T))
      JNAD=REAL(X)
C
      RETURN
      END
******************************************************************************
*     SUBROUTINE SURF_FRZ_SNW_NAT(TAIR,WSA,WNA,JNAT)
******************************************************************************
*
*     Subroutine used to calculate the surface freezing rate of
*     nitric acid trihydrate (NAT) in supercooled ternary solution.
*
*     Source: Tabazadeh et al., J. Phys. Chem. A. 106, 10238, 2002
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  TAIR,WSA,WNA,J
*
*     Input: 
*         TAIR:    (K)            Ambient air temperature
*         WSA:                    Weight fraction of H2SO4  ]0.0; 0.7[
*         WNA:                    Weight fraction of HNO3   ]0;0.7[
*
*     Output:
*         JNAT:    (1/m**2 s)     Homogeneous freezing rate
*
      SUBROUTINE SURF_FRZ_SNW_NAT(TAIR,WSA,WNA,JNAT)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WSA,WNA,JNAT
C----------------------------------------------------------------------------
      REAL(KIND=8)NS,KBOLTZ,RGAS,H
      PARAMETER(
C       Number density of molecular surface sites (1/m**2)
     +          NS=1.0D19,
C       Boltzmanns constant (J/K)
     +          KBOLTZ=1.380662D-23,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441D0,
C       Planck's constant (J s)
     +          H=6.626176D-34)
      REAL(KIND=4) MH2SO4,MHNO3,MH2O,XX
      PARAMETER(
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
C       Molar weight of sulfuric acid (kg/mole)
     +          MH2SO4=98.08E-3,
C       Molar weight of nitric acid (kg/mole)
     +          MHNO3=63.01E-3)
C     Auxiliary local variables:
      REAL(KIND=8) C,PRE,X,G,GNAD,GNAT,T,XHNO3,TEMP
      PARAMETER(C=NS*KBOLTZ/H)
C----------------------------------------------------------------------------
      GNAD(XHNO3,TEMP)=(11.5593D0+0.0804214D0*TEMP-
     +                 (71.5133D0-0.256724D0*TEMP)*XHNO3)
      GNAT(TEMP)=-45.2429D0+0.364844D0*TEMP
      XX=(WSA/MH2SO4)+(WNA/MHNO3)+((1.0-WSA-WNA)/MH2O)
      X=DBLE((WNA/MHNO3)/XX)
      T=DBLE(TAIR)
      PRE=C*X*T
      G=GNAD(X,T)*GNAT(T)/GNAD(0.246D0,T)
      X=PRE*DEXP(-G/(RGAS*T))
      JNAT=REAL(X)
C
      RETURN
      END
******************************************************************************
*     REAL FUNCTION FREEZE(TAIR,PPWV)
******************************************************************************
*
*     Function used to calculate a volume freezing rate according to 
*             J(Tair) = B*EXP(a*(T0-Tair))
*     cf. Pruppacker & Klett (1980), p. 275.
*     The number of unfrozen particles N is calculated in the calling
*     program unit from the differential equation
*             dN = - N * J * V * dt
*     where V is the particle volume; i.e. volume dependent freezing
*     is assumed. The constant of proportionality B is given by
*             B = ln(2)/(Tm1*Vd)
*     where Tm1 is the time required to freeze half a population of particles
*     of volume Vd at temperature T0. The sensitivity a is calculated from 
*     the time Tm2 required to freeze half a population of of particles of 
*     volume Vd at a temperature T0+dT as
*             a= ln(Tm2/Tm1)/dT
*     No freezing is assumed to take place for Tair > T0+dT.
*     Values of Tm1, Tm2, Vd, and dT are specified as parameters in the
*     subroutine. T0 is assumed at the ice frost point temperature minus DT.
*
*     Declaration of input/output variables:
*
*     REAL(KIND=4)  T0,TAIR
*
*     Input: 
*
*         TAIR:    (K)            Ambient air temperature
*         PPWV:    (hPa)          Ambient water vapor partial pressure
*
*     Output:
*
*         FREEZE:  (1/m**3 s)     Volume freezing rate
*
*
      REAL FUNCTION FREEZE(TAIR,PPWV)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PPWV
      DOUBLE PRECISION  A,B,PI,RAD,VD,TM1,TM2,DT,T0
      REAL(KIND=4)  TCICE
      LOGICAL FIRST
C
      PARAMETER(
     +        PI=3.1415926536D0,
C         Radius of particles with volume Vd (m):
     +        RAD=1.0D-6,
     +        VD=4.0D0*PI*RAD*RAD*RAD/3.0D0,
C         Time to freeze 1/2 population at T0 (sec):
     +        TM1=1.0D0,
C         Time to freeze 1/2 population at T0+DT (sec):
     +        TM2=100.0D0*24.0D0*3600.0D0,
C         Temperature interval (K) (no freezing above T0+DT):
     +        DT=8.7D0)
      DATA FIRST/.TRUE./
      SAVE FIRST,A,B
C
      IF(FIRST) THEN
          FIRST=.FALSE.
          A=DLOG(TM2/TM1)/DT
          B=DLOG(2.0D0)/(TM1*VD)
      ENDIF
C
      T0=DBLE(TCICE(PPWV))-0.5D0*DT
      IF(TAIR.GT.T0+DT) THEN
          FREEZE=-1.0E0
      ELSE
          FREEZE=REAL(B*DEXP(A*(T0-TAIR)))
      ENDIF
      RETURN
      END

*****************************************************************************
*     SUBROUTINE THERMO(TAIR,PAIR,PPWV,PPNA)
*****************************************************************************
*
*     This subroutine calculates various thermodynamical parameters,
*     all depending on the state of the ambient air, specified by the
*     air temperature (TAIR), air pressure (PAIR), partial pressure of
*     water vapor (PPWV), and partial pressure of nitric acid vapor (PPNA).
*
*     All calculated thermodynamic variables are stored in the named
*     common block /TERMO/, and holds the following values:
*
*     COMMON /TERMO/SNANAT,PSNNAT,SWVICE,PSWICE,STNAT,STICE,
*                   RHONAT,RHOICE,RHOAIR,MFPAIR,DYNVIS,KAIR,DWMAIR,LWV
*     where
*
*     SNANAT:      Saturation ratio of nitric acid over NAT
*     PSNNAT:      Saturation pressure of nitric acid over NAT (Pa)
*     SWVICE:      Saturation ratio of water vapor over ice
*     PSWICE:      Saturation pressure of water vapor over ice (Pa)
*     STNAT:       Surface tension of NAT/air (N/m)
*     STICE:       Surface tension of ice/air (N/m)
*     RHONAT:      Density of NAT (kg/m**3)
*     RHOICE:      Density of ice (kg/m**3)
*     RHOAIR:      Density of air (kg/m**3)
*     MFPAIR:      Mean free path of air molecules (m)
*     DYNVIS:      Dynamic viscosity of air (kg/m s)
*     KAIR:        Thermal conductivity of air (J/m s K)
*     DWMAIR:      Diffusivity of water molecules in air (m**2/s)
*     LSNA:        Heat of sublimation of nitric acid on NAT (J/kg)
*     LWV:         Heat of sublimation of water vapor (J/kg)
*
*     Input:  TAIR         Temperature (K)
*             PAIR:    (Pa)           Ambient air pressure
*             PPWV:    (Pa)           Partial water vapor pressure 
*             PPNA:    (Pa)           Partial pressure of nitric acid 
*
      SUBROUTINE THERMO(TAIR,PAIR,PPWV,PPNA)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PAIR,PPWV,PPNA
C
C     External thermodynamic functions needed:
      REAL(KIND=4)  PNANAT,PWVICE,STNAT,ROICE,
     +     CDTAIR,FPLAIR,DFWVA,VISAIR,LSWV
C     PNANAT:      Nitric acid pressure over nitric acid trihydrate
C     PWVICE:      Water vapor pressure over ice
C     STNAT:       Surface tension of nitric acid trihydrate
C     ROICE:       Density of ice
C     CDTAIR:      Thermal conductivity of air
C     FPLAIR:      Mean free path of air molecules
C     DFWVA:       Diffusion coefficient of water vapor in air
C     VISAIR:      Viscosity of air
C     LSWV:        Heat of sublimation of water vapor 
C
C     Physical constants:
      REAL(KIND=4)  MAIR,RGAS,RONAT,LNA,STICE
      PARAMETER(
C       Molar weight of dry air (kg/mole)
     +          MAIR=28.9644E-3,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0,
C       Density of nitric acid trihydrate (kg/m**3) (Taesler et al. 1975)
     +          RONAT=1.621E+3,
C       Density of nitric acid trihydrate (kg/m**3) (Drdla & Turco 1990)
CC     +          RONAT=1.35E03,
C       Heat of sublimation of nitric acid on NAT (J/kg)
C               (Hanson & Mauersberger 1988b, Wofsy et al., 1990)
     +          LNA=1.86E+6,
C       Surface tension ice/air (N/m) (Pruppacher & Klett, 1980, p. 121)
     +          STICE=105.0E-3)
C
C     Common block variables:
      REAL(KIND=4)     SNANAT,PSNNAT,SWVICE,PSWICE,SFTNAT,SFTICE,
     +        RHONAT,RHOICE,RHOAIR,MFPAIR,DYNVIS,KAIR,DWMAIR,LSNA,LWV
      COMMON /TERMO/SNANAT,PSNNAT,SWVICE,PSWICE,SFTNAT,SFTICE,
     +        RHONAT,RHOICE,RHOAIR,MFPAIR,DYNVIS,KAIR,DWMAIR,LSNA,LWV
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C
C
C     Auxiliary local variables:
      REAL(KIND=4)  C1
      PARAMETER(C1=MAIR/RGAS)
C
      PSNNAT=PNANAT(TAIR,PPWV)
      SNANAT=PPNA/PSNNAT
      PSWICE=PWVICE(TAIR)
      SWVICE=PPWV/PSWICE
      SFTNAT=STNAT(TAIR)
      SFTICE=STICE
      RHONAT=RONAT
      RHOICE=ROICE(TAIR)
      MFPAIR=FPLAIR(TAIR,PAIR)
      RHOAIR=C1*PAIR/TAIR
      DYNVIS=VISAIR(TAIR)
      DWMAIR=DFWVA(TAIR,PAIR)
      KAIR=CDTAIR(TAIR)
      LSNA=LNA
      LWV=LSWV(TAIR)
C
      IF(TEST.LT.50) RETURN
C
      WRITE(TEST,*)
      WRITE(TEST,*) 'TESTING THERMO'
      WRITE(TEST,1010) 'TAIR','PAIR','PPWV','PPNA'
      WRITE(TEST,1000) TAIR,PAIR,PPWV,PPNA
      WRITE(TEST,1010) 'SNANAT','PSNNAT','SWVICE','PSWICE',
     +                 'STNAT','STICE','RHOICE','MFPAIR'
      WRITE(TEST,1000) SNANAT,PSNNAT,SWVICE,PSWICE,SFTNAT,SFTICE,RHOICE,
     +        MFPAIR
      WRITE(TEST,1010) 'RHONAT','RHOAIR','DYNVIS','DWMAIR','KAIR','LWV'
      WRITE(TEST,1000) RHONAT,RHOAIR,DYNVIS,DWMAIR,KAIR,LWV
 1000 FORMAT(1X,8(1PE9.2),0P)
 1010 FORMAT(1X,8(A9))
C
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION FPLAIR(TAIR,PAIR)
      REAL FUNCTION FPLAIR(TAIR,PAIR)
*****************************************************************************
*
*     Molecular mean free path of air molecules.
*
*     Source: Prupacher & Klett:Microphysics of clouds and precipitation,
*                               (1980), 10-106, p.323.
*     Formula used:             FPLAIR=2.281238E-5*TAIR/PAIR
*
*     Input:  TAIR:     Air temperature (K)
*             PAIR:     Air pressure (Pa)
*     Output:           Molecular mean free path (m)
*
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PAIR
      FPLAIR=(2.281238E-5)*TAIR/PAIR
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION CDTAIR(TAIR)
      REAL FUNCTION CDTAIR(TAIR)
*****************************************************************************
*
*     Thermal conduvtivity of air.
*
*     Source: Prupacher & Klett:Microphysics of clouds and precipitation,
*                               (1980), p 418, 13-16
*     Formula used:             CDTAIR=4.381276E-3+7.117560E-5*TAIR
*
*     Input:  TAIR:     Air temperature (K)
*     Output:           Thermal conductivity of air (J/(m sec K))
*
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR
      CDTAIR=4.381276E-3+7.117560E-5*TAIR
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION VISAIR(TAIR)                                                 *
      REAL FUNCTION VISAIR(TAIR)
*****************************************************************************
*
*     Dynamic viscosity of air.
*
*     Source: Prupacher & Klett:Microphysics of clouds and precipitation,
*                       p. 323.
*
*     Input:  TAIR:   Temperature (K)
*     Output: VISAIR: Dynamic viscosity of air  (kg/(m s))
*
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,T
C
      T=TAIR-273.15
      VISAIR=1.0E-5*((-1.2E-5*T+0.0049E0)*T+1.718E0)
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION ROMAIR(TAIR,PAIR,PPWV)                                        *
      REAL FUNCTION ROMAIR(TAIR,PAIR,PPWV)
*****************************************************************************

*
*     Density of moist air.
*     Source: R.J.List: Smithsonian Meteorological Tables, 6.ed.,1951,
*                       p. 295, 347.
*     Formula used:     ROMAIR= (PAIR+PPWV(EPSILON-1))/TR
*
*     Input:
*               TAIR:   ambient temperature (K)  
*               PAIR:   ambient pressure (Pa)    
*               PPWV:   ambient partial pressure of water vapor (Pa)
*     Output:
*               ROMAIR: density of moist air (kg/m**3)
*
C
C       Input/output variables:
        IMPLICIT NONE
        REAL(KIND=4)  TAIR,PAIR,PPWV
C
C       Physical constants:
        REAL(KIND=4)  MH2O,MAIR,RGAS,RAIR
        PARAMETER(
C       Molar mass og water vapor (kg/mole)
     +          MH2O=18.0160E-3,
C       Molar mass of air (kg/mole)
     +          MAIR=28.966E-3,
C       Universal gas constant (J/(mole K))
     +          RGAS=8.31441E0,
C       Gas constant of dry air (J/(kg K))
     +          RAIR=RGAS/MAIR)
C
C       Auxiliary local variables:
        REAL(KIND=4)  EPSI
        PARAMETER (EPSI=MH2O/MAIR-1.0E0)
C
        ROMAIR=(PAIR+PPWV*EPSI)/(RAIR*TAIR)
        RETURN
        END
C
C
*****************************************************************************
*     REAL FUNCTION CPMAIR(PAIR,PPWV)
      REAL FUNCTION CPMAIR(PAIR,PPWV)
*****************************************************************************

*
*     Specific heat capacity at constant pressure of moist air.
*     Source: R.J.List: Smithsonian Meteorological Tables, 6.ed., 
*                       p. 339, 347.
*     Formula used:     CPMAIR= 7/2 RA(1+8/7(R/EPSILON))
*
*     Input:
*               PAIR:    Ambient pressure (Pa)
*               PPWV:    Ambient partial pressure of water vapor (Pa)
*     Output:
*               CPMAIR:  Specific heat capacity  (J/(kg K)).
*
C
C       Input/output variables:
        IMPLICIT NONE
        REAL(KIND=4)  PAIR,PPWV
C
C       Physical constants:
        REAL(KIND=4)  MAIR,RGAS,RAIR
        PARAMETER(
C       Molar mass of air (kg/mole)
     +          MAIR=28.966E-3,
C       Universal gas constant (J/(mol K))
     +          RGAS=8.31441E0,
C       Gas constant of dry air (J/(kg K))
     +          RAIR=RGAS/MAIR)
C
C       Auxiliary local variables:
        REAL(KIND=4)  C1,C2
        PARAMETER (C1=7.0E0*RAIR/2.0E0, C2=4.0E0*RAIR)
C
        CPMAIR=C1+C2*PPWV/(PAIR-PPWV)
        RETURN
        END
*****************************************************************************
*     REAL FUNCTION DFWVA(TAIR,PAIR)                                                 *
      REAL FUNCTION DFWVA(TAIR,PAIR)
*****************************************************************************
*
*     Diffusivity of water vapor in air.
*
*     Source: Prupacher & Klett:Microphysics of clouds and precipitation,
*                               (1980), 13-3, p. 413
*
*     The relation D = E0 (T/T0)**n (P0/P); n=1.94 has been used.
*
*     Input:  TAIR: Temperature (K);  Range:  [180    ,273    ]
*             PAIR: Pressure (Pa)
*     Output: Diffusivity of water vapor in air  (m**2/sec)
*
C
      IMPLICIT NONE
      REAL(KIND=4)  PAIR,TAIR
      REAL(KIND=4)  E0,D1,P0,T0
      PARAMETER(E0=0.211E-4,P0=1.01325E+5,T0=273.15E0,D1=E0*P0)
C
      DFWVA=D1*((TAIR/T0)**1.94E0)/PAIR
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION LEWV(TAIR)
      REAL FUNCTION LEWV(TAIR)
*****************************************************************************
*
*     Latent heat of evaporation of water vapor.
*
*     Source: R.J.List: Smithsonian Meteorological Tables, 6.ed.,1951,
*                       p. 343. Lin. extrap.
*
*     Input:  TAIR:     Air temperature (K)
*     Output:           Heat of evaporation  (J/kg)
*
*     Table values of L, fitted with a 1.order polynomium.
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR
      LEWV=-2.65857E+3*TAIR+3.22427E+6
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION LSWV(TAIR)                                                    *
      REAL FUNCTION LSWV(TAIR)
*****************************************************************************
*
*     Latent heat of sublimation of water vapor.
*
*     Source: R.J.List: Smithsonian Meteorological Tables, 6.ed.,1951,
*                       p. 343.
*     Table values of L, fitted with a 2.order polynomium.
*
*     Input:  TAIR:     Air temperature (K)
*     Output: Heat of sublimation  (J/kg)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,A0,A1,A2
      DATA A0,A1,A2/
     # 2.63335D+06, 1.72542D+03,-3.62100D+00/
      LSWV=(A2*TAIR+A1)*TAIR+A0
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION LMICE(TAIR)                                                    *
      REAL FUNCTION LMICE(TAIR)
*****************************************************************************
*
*     Latent heat of melting of ice
*
*     Source: Prupacher & Klett:Microphysics of clouds and precipitation,
*                               (1980), 4-85, p. 89
*     Table values of L, fitted with a 2nd order polynomium.
*
*     Input:  TAIR:     Air temperature (K)
*     Output:           Heat of melting  (J/mol)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,TT
C
C     Physical constants:
      REAL(KIND=4)  MH2O,C1
      PARAMETER(
C       Molar weight of water (kg/mole)
     +          MH2O=18.0153E-3,
     +          C1=4.1868E+3*MH2O)
      TT=TAIR-273.15E0
      LMICE=(79.7E0+(0.485E0-2.5E-3*TT)*TT)*C1
C      LMICE=6005.2356+18.2719*TT-0.06354*TT*TT
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION ROWL(TAIR)                                              *
      REAL FUNCTION ROWL(TAIR)
*****************************************************************************
*
*     Density of subcooled water.
*
*     Pruppacher & Klett, Microphysics of clouds and precipitation
*                         (1997), p. 87
*     Input:  TAIR: Temperature (K)
*     Output:       Density  (kg/m**3)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,T,A0,A1,A2,A3,A4,A5,A6
      DATA A0,A1,A2,A3,A4,A5,A6/
     +     0.99986E0,6.690E-5,-8.486E-6,
     +     1.518E-7,-6.9984E-9,-3.6449E-10,-7.497E-12/
      T=AMAX1(-45.0E0,TAIR-273.15E0)
      ROWL=1.0E+3*((((((A6*T+A5)*T+A4)*T+A3)*T+A2)*T+A1)*T+A0)
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION ROICE(TAIR)                                             *
      REAL FUNCTION ROICE(TAIR)
*****************************************************************************
*
*     Density of ice.
*
*     Source: Prupacher & Klett:Microphysics of clouds and precipitation,
*                               (1980), p. 86
*
*     Input:  TAIR:  Temperature (K)
*     Output:        Density  (kg/m**3)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,A0,A1,A2,TC
      PARAMETER(A0=0.9167E+3,A1=-1.75E-1,A2=-5.0E-4)
      TC=TAIR-273.15E0
      ROICE=A0+(A1+A2*TC)*TC
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION ROSAS(TAIR,WSA)                                                 *
      REAL FUNCTION ROSAS(TAIR,WSA)
*****************************************************************************
*
*     Density of liquid sulfuric acid solution.
*
*     Source: John H.Perry (ed.):Chemical Engineers' Handbook,
*                                McGraw-Hill, New York 1963, p. 3-79 & 3-80
*
*     The original data set in temp. range 0 C to 20 C and weight pct.
*     0 to 100 % has been fitted with a polynomium of two variables
*     of order 5 in W and lineary in T. Fit quality better than 0.5 %
*
*     Input:  TAIR: Temperature  (K)
*             WSA:  Weight fraction of H2SO4  [0;1] 
*     Output: Density of sulfuric acid solution  (kg/m**3)
*
C
      IMPLICIT NONE
      INTEGER I
      REAL(KIND=4)  TAIR,WSA
      REAL(KIND=4)  C(6)
      REAL(KIND=4)  A(6)
      REAL(KIND=4)  B(6)
      REAL(KIND=4)  D(6)
      DATA (A(I),I=1,6)/
     # 1.00190D+03, 5.50496D+02, 1.54093D+03,-4.89219D+03, 7.56555D+03,
     #-3.92739D+03/
      DATA (B(I),I=1,6)/
     # 1.98378D+01, 1.02256D+03,-1.48665D+03,-7.24651D+02, 3.68348D+03,
     #-2.22159D+03/
      DATA (D(I),I=1,6)/
     #-6.97011E-02,-3.59886D+00, 5.24992D+00, 2.54047D+00,-1.29355D+01,
     # 7.80553D+00/
C
      DO I=1,6
              C(I)=A(I)+B(I)+D(I)*TAIR
      ENDDO
      ROSAS=C(1)+WSA*(C(2)+WSA*(C(3)+WSA*(C(4)+WSA*(C(5)+WSA*C(6)))))
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION RONAS(TAIR,WNA)                                                 *
      REAL FUNCTION RONAS(TAIR,WNA)
*****************************************************************************
*
*     Density of liquid nitric acid solution.
*
*     Source: Granzhan and Laktionova, Russian J. Phys. Chem. 49, 1448, 1975
*
*
*     Input:  TAIR: Temperature  (K)
*             WNA:  Weight fraction of HNO3  [0;1]
*     Output: Density of nitric acid solution  (kg/m**3)
*
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WNA
      REAL(KIND=4)  TC
C
      TC=TAIR-293.15E0
      RONAS=(0.9982E0-0.396E-3*TC+
     +   ((-0.3708E0*WNA+0.3278E0)*WNA+0.5452E0-0.1485E-2*TC)*WNA)*1.0E3
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION ROSNW(TAIR,WSA,WNA)                                                 *
CC    REAL FUNCTION ROSNW(TAIR,WSA,WNA)
*****************************************************************************
*
*     Density of liquid TERNARY sulfuric/nitric acid solution.
*
*
*     Input:  TAIR: Temperature  (K)
*             WSA:  Weight fraction of H2SO4  [0;1] 
*             WNA:  Weight fraction of HNO3  [0;1] 
*     Output: Density of sulfuric/nitric acid solution  (kg/m**3)
*
C
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,WSA,WNA
CC    REAL(KIND=4)  RN,RS,MN,MS,ROWL,ROSAS,RONAS,MHNO3,MH2SO4,WN,WS
CC    PARAMETER(
C       Molar weight of nitric acid (kg/mole)
CC   +          MHNO3=63.01E-3,
C       Molar weight of sulfuric acid (kg/mole)
CC   +          MH2SO4=98.08E-3)
C
CC    IF(WSA.EQ.0.0E0.AND.WNA.EQ.0.0E0) THEN
CC        ROSNW=ROWL(TAIR)
CC    ELSE
CC        MN=WNA/(MHNO3*(1-WNA-WSA))
CC        MS=WSA/(MH2SO4*(1-WNA-WSA))
CC        WN=MHNO3*MN/(1.0E0+MHNO3*MN)
CC        WS=MH2SO4*MS/(1.0E0+MH2SO4*MS)
CC        RN=RONAS(TAIR,WN)
CC        RS=ROSAS(TAIR,WS)
CC        ROSNW=1.0E0/((1.0E0/RS)*MS/(MN+MS)+(1.0E0/RN)*MN/(MN+MS))
CC    ENDIF
C
CC    RETURN
CC    END
c
*****************************************************************************
*     REAL FUNCTION ROSNW(TAIR,WSA,WNA)                                     *
CC    REAL FUNCTION ROSNW(TAIR,WSA,WNA)
*****************************************************************************
*
*     Density of liquid TERNARY sulfuric/nitric acid solution.
*
*     Source: Carslaw et al., GRL 22,1877,1995
*
*     Input:  TAIR: Temperature  (K)
*             WSA:  Weight fraction of H2SO4  [0;1] 
*             WNA:  Weight fraction of HNO3  [0;1] 
*     Output: Density of sulfuric/nitric acid solution  (kg/m**3)
*
CC    IMPLICIT NONE
CC    REAL TAIR,WSA,WNA,MN,MS,RHOS,RHON
CC    REAL MHNO3,MH2SO4
CC    PARAMETER(
C       Molar weight of nitric acid (kg/mole)
CC   +          MHNO3=63.01E-3,
C       Molar weight of sulfuric acid (kg/mole)
CC   +          MH2SO4=98.08E-3)
CC    MN=WNA/(MHNO3*(1-WNA-WSA))
CC    MS=WSA/(MH2SO4*(1-WNA-WSA))
CC    RHOS=1000+(123.64-5.6E-4*TAIR**2)*MS-
CC   +        (29.54-1.81E-4*TAIR**2)*MS**1.5+
CC   +        (2.343-1.487E-3*TAIR-1.324E-5*TAIR**2)*MS**2
CC    RHON=1000+(85.11-5.04E-4*TAIR**2)*MN-
CC   +        (18.96-1.427E-4*TAIR**2)*MN**1.5+
CC   +        (1.458-1.198E-3*TAIR-9.703E-6*TAIR**2)*MN**2
CC    ROSNW=1.0/((1.0/RHOS)*MS/(MN+MS)+(1.0/RHON)*MN/(MN+MS))
CC    RETURN
CC    END
C    Density of ternary solution in g/cm3
*****************************************************************************
*     REAL FUNCTION ROSNW(TAIR,WSA,WNA)                                                 *
CC    REAL FUNCTION ROSNW(TAIR,WSA,WNA)
*****************************************************************************
*
*     Density of liquid TERNARY sulfuric/nitric acid solution.
*
*     Source: Luo et al. GRL 23, 3707,1996,
*
*     Input:  TAIR: Temperature  (K)
*             WSA:  Weight fraction of H2SO4  [0;1] 
*             WNA:  Weight fraction of HNO3  [0;1] 
*     Output: Density of sulfuric/nitric acid solution  (kg/m**3)
*
CC    IMPLICIT NONE
CC    REAL TAIR,WSA,WNA,W,T,WTH,V1,VS,VN,VMCAL
CC    REAL X(22)
CC    DATA X /2.393284E-02,-4.359335E-05,7.961181E-08,0.0,
CC   * -0.198716351,1 .39564574E-03,-2.020633E-06,
CC   * 0.51684706,-3.0539E-03,4.505475E-06,-0.30119511,
CC   * 1.840408E-03,-2 .7221253742E-06,-0.113316741169,
CC   * 8.47763E-04,-1 .22336185E-06,0.3455282,-2.2111E-03,
CC   * 3.503768245E-06,-0.2315332 ,1.60074E-03,-2.5827835E-06/
CC    W=WSA+WNA
CC    T=TAIR
CC    WTH=1.-W
CC    V1=X(1)+X(2)*T+X(3)*T**2+X(4)*T**3
CC    VS=X(5)+X(6)*T+X(7)*T**2+(X(8)+X(9)*T+X(10)*T**2)*W
CC   * +(X(11)+X(12)*T+X(13)*T**2)*W*W
CC    VN=X(14)+X(15)*T+X(16)*T**2+(X(17)+X(18)*T+X(19)*T**2) *W
CC   * +(X(20)+X(21)*T+X(22)*T**2)*W*W
CC    VMCAL=WTH/18.016D0*V1+VS*WSA/98.08D0+VN*WNA/63.016D0
CC    ROSNW=1/VMCAL
CC    RETURN
CC    END
*****************************************************************************
*     REAL FUNCTION ROSNW(TAIR,WSA,WNA)                                                 *
      REAL FUNCTION ROSNW(TAIR,WSA,WNA)
*****************************************************************************
*
*     Density of liquid TERNARY sulfuric/nitric acid solution.
*
*     Source: Martin et al. GRL 27, 197, 2000
*
*     Input:  TAIR: Temperature  (K)
*             WSA:  Weight fraction of H2SO4  [0;1] 
*             WNA:  Weight fraction of HNO3  [0;1] 
*     Output: Density of sulfuric/nitric acid solution  (kg/m**3)
      IMPLICIT NONE
      REAL TAIR,WSA,WNA,X,Y,T293,T253,A0,A1,A2,A3,A4,R293,R253,DT,
     +R00,R01,R02,R03,
     +R10,R11,R12,R13,
     +R20,R21,R22,
     +R30,R31,
     +R40,
     +S00,S01,S02,S03,
     +S10,S11,S12,S13,
     +S20,S21,S22,
     +S30,S31,
     +S40
      DATA
     +R00,R01,R02,R03,
     +R10,R11,R12,R13,
     +R20,R21,R22,
     +R30,R31,
     +R40/
     + 0.9982E00, 5.8487E-3,-2.3873E-5, 7.9815E-7,
     + 7.9119E-3, 2.9369E-5, 1.8949E-6,-3.6905E-8,
     +-7.6431E-5, 2.8093E-6,-6.4247E-8,
     + 2.2885E-6,-2.1422E-8,
     +-1.4651E-8/
      DATA
     +S00,S01,S02,S03,
     +S10,S11,S12,S13,
     +S20,S21,S22,
     +S30,S31,
     +S40/
     + 1.0015E00, 9.7509E-3,-1.8340E-4,2.9113E-6,
     + 9.6589E-3,-3.9433E-5,3.8149E-6,-6.3144E-8,
     +-1.1562E-4, 4.3442E-6,-7.0749E-8,
     + 2.6848E-6,-3.6871E-8,
     +-1.6015E-8/
      PARAMETER(T293=293.0,T253=253.0,DT=T293-T253)
C
      X=100.0*WSA
      Y=100.0*AMIN1(WNA,0.45)
C
      A0=R00+Y*(R01+Y*(R02+R03*Y))
      A1=R10+Y*(R11+Y*(R12+R13*Y))
      A2=R20+Y*(R21+Y*R22)
      A3=R30+Y*R31
      A4=R40
      R293=A0+X*(A1+X*(A2+X*(A3+X*A4)))
C
      A0=S00+Y*(S01+Y*(S02+S03*Y))
      A1=S10+Y*(S11+Y*(S12+S13*Y))
      A2=S20+Y*(S21+Y*S22)
      A3=S30+Y*S31
      A4=S40
      R253=A0+X*(A1+X*(A2+X*(A3+X*A4)))
C
      ROSNW=1000.0*(R253+(R293-R253)*(TAIR-T253)/DT)
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION TMLTSA(WSAS)
      REAL FUNCTION TMLTSA(WSAS)
*****************************************************************************
*
*     Solid/liquid equilibrium temperature of sulfuric acid tetrahydrate
*     and ice as a function of H2SO4 weight fraction.
*
*     Source: Gable et al. J. Amer. Chem. Soc. 72, 1445, 1950.
*
*     Input:      WSA:                Weight fraction of H2SO4  [0;1]
*     Output:     TMLTSA      (K)     Equilibrium temperature
*
      IMPLICIT NONE
      REAL(KIND=4)  WSAS,
     +        A0,A1,A2,A3,A4,
     +        B0,B1,B2,B3,B4
      PARAMETER(
     +        A0=2.72832E+2,A1=-2.83140E0,A2=-6.79090E+2,A3=3.30407E+3,
     +        A4=-7.64069E+3,
     +        B0=-7.01591E+2,B1=6.10969E+3,B2=-1.61399E+4,B3=2.05988E+4,
     +        B4=-1.04882E+4)
C
      IF(WSAS.LE.0.3712E0) THEN
          TMLTSA=(((A4*WSAS+A3)*WSAS+A2)*WSAS+A1)*WSAS+A0
      ELSE
          TMLTSA=(((B4*WSAS+B3)*WSAS+B2)*WSAS+B1)*WSAS+B0
      ENDIF
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION TMSAT(PPWV)
      REAL FUNCTION TMSAT(PPWV)
*****************************************************************************
*     Melting temperature of sulfuric acid tetrahydrate (SAT).
*     Source: Tabazadeh et al., GRL 21,1619,1994.
*
*     Input/output:
*
*         PPWV:    (Pa)           Water vapor partial pressure
*         TMSAT:   (K)            SAT melting temperature 
*
      IMPLICIT NONE
      REAL(KIND=4)  PPWV
C
      TMSAT=3236.0E0/(11.502E0-ALOG10(PPWV/133.32239E0))
      RETURN
      END
C****************************************************************************
CC      REAL FUNCTION STICE(TAIR)
CC      REAL FUNCTION STICE(TAIR)
*****************************************************************************
*
*     Surface tension of ice/vapor.
*
*     Source: O.B.Toon et al.: JGR 94,11359, 1989
*      (NB: According to Pruppacher & Klett, (1980), p. 121 
*      STICE = 105.0E-3 N/m (constant))
*
*     Input:  TAIR:  Temperature (K)
*     Output:        Surface tension  (J/m**2)
*
CC      IMPLICIT NONE
CC      REAL(KIND=4)  TAIR
CC      STICE=0.141E0-0.15E-3*TAIR
CC      RETURN
CC      END
*****************************************************************************
*     REAL FUNCTION STWL(TAIR)                                                    *
      REAL FUNCTION STWL(TAIR)
*****************************************************************************
*
*     Surface tension of water/vapor.
*
*     Source: L.G.Kachurin & V.G. Morachevskii:
*             Kinetics og phase transitions of water in the atmosphere
*             Israel program for scientific translations
*             Jerusalem 1966.
*
*     Input:  TAIR:  Temperature (K)
*     Output:        Surface tension  (J/m**2)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR
      STWL=0.1092E0-0.1227E-3*TAIR
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION STWLIC(TAIR)                                                    *
      REAL FUNCTION STWLIC(TAIR)
*****************************************************************************
*
*     Surface tension of ice/water.
*
*     Source: Prupacher & Klett:Microphysics of clouds and precipitation,
*                               (1980), p. 130
*
*     Input:  TAIR:  Temperature (K)
*     Output:        Surface tension  (J/m**2)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,C1,C2
      PARAMETER(C1=28.5E-3-273.15E0*0.25E-3, C2=0.25E-3)
      STWLIC=C1+C2*TAIR
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION STNAT(TAIR)                                                 *
      REAL FUNCTION STNAT(TAIR)
*****************************************************************************
*
*     Surface tension of nitric acid trihydrate/vapor.
*
*     Source: Drdla & Turco: J.Atm.Chem. 12,319,1991
*
*     Input:  TAIR: Temperature (K)
*     Output: Surface tension of nitric acid trihydrate (J/m**2)
*
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR
      REAL(KIND=4)  C1,C2
      PARAMETER(C1=112.0E-3, C2=-0.155E-3)
C
      STNAT=C1+C2*TAIR
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION STSAIC(TAIR,WSA)                                                 *
      REAL FUNCTION STSAIC(TAIR,WSA)
*****************************************************************************
*
*     Surface tension of sulfuric acid solution/water ice.
*
*     Application of Antonoff's rule: STSAIC=ABS(STSAS(TAIR,WSA)-STICE)    
*     Source: Sabinina & Terpugov: Z. Phys. Chem. A173 ,237, 1935.
*
*
*     Input:  TAIR:   Temperature (K)
*             WSA:    Weight fraction of H2SO4  [0;1]
*     Output: STSAIC: Surface tension of sulfuric acid solution/ice (N/m)
*
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WSA,STSAS,STICE
      PARAMETER(
C       Surface tension ice/air (N/m) (Pruppacher & Klett, 1980, p. 121)
     +          STICE=105.0E-3)
      STSAIC=ABS(STSAS(TAIR,WSA)-STICE)
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION STSAS(TAIR,WSA)                                         *
CC    REAL FUNCTION STSAS(TAIR,WSA)
*****************************************************************************
*
*     Surface tension of sulfuric acid solution/vapor.
*
*     Source: Tabazadeh et al. JGR, 102,23845,1997
*             Sabinina & Terpugov: Z. Phys. Chem. A173 ,237, 1935.
*
*
*     Input:  TAIR: Temperature (K)
*             WSA:  Weight fraction of H2SO4  [0;1]
*     Output: Surface tension of sulfuric acid solution (N/m)
*
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,WSA,W
CC    W=WSA*100.0E0
CC    STSAS=1.0E-3*(142.35E0-0.96525E0*W-TAIR*(0.22954E0-0.0033948E0*W))
CC    RETURN
CC    END
*****************************************************************************
*     REAL FUNCTION STSAS(TAIR,WSA)                                         *
      REAL FUNCTION STSAS(TAIR,WSA)
*****************************************************************************
*
*     Surface tension of sulfuric acid solution/vapor.
*
*     Source: Tabazadeh et al. submitted,1999
*             Myhre et al., J. Chem. Eng. Data 43,617,1998.
*
*
*     Input:  TAIR: Temperature (K)
*             WSA:  Weight fraction of H2SO4  [0;1]
*     Output: Surface tension of sulfuric acid solution (N/m)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WSA
      DOUBLE PRECISION W,T,S180,S220,S260
C
      W=DBLE(WSA)*100.0D0
      T=DBLE(TAIR)
      IF(W.LT.40.0D0) THEN
          T=DMAX1(180.0D0,DMIN1(T,260.0D0))
          S220=(((((8.969257061D-7*W-1.145573827D-4)*W+5.415260617D-3)
     +          *W-1.050692123D-1)*W+5.312072092D-1)*W+82.01197792D0)
          IF(T.LE.220.0D0) THEN
            S180=(((((1.736789787D-6*W-1.912224154D-4)*W+7.485866933D-3)
     +          *W-1.103647657D-1)*W+9.541966318D-2)*W+85.75507114D0)
            STSAS=REAL(1.0D-3*(S220+(5.5D0-0.025D0*T)*(S180-S220)))
          ELSE IF(T.GT.220.0D0) THEN
            S260=(((((2.095358048D-7*W-2.384669516D-5)*W+8.87979880D-4)
     +          *W-9.682499074D-3)*W-6.9631232740D-3)*W+77.40682664D0)
            STSAS=REAL(1.0D-3*(S260+(6.5D0-0.025D0*T)*(S220-S260)))
          ENDIF
      ELSE
          STSAS=1.0E-3*
     +        REAL(142.35D0-0.96525D0*W-TAIR*(0.22954D0-0.0033948D0*W))
      ENDIF
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION STSAST(TAIR)                                            *
CC    REAL FUNCTION STSAST(TAIR)
*****************************************************************************
*
*     Surface tension of sulfuric acid tetrahydrate/liquid solution.
*
*     Source: Luo et al., Ber. Bunsenges. Phys. Chem. 96,334,1992
*
*
*     Input:  TAIR: Temperature (K)
*     Output: Surface tension (N/m)
*
C
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR
CC    REAL(KIND=4)  B,NA,N,X,C1,C2,C3,A1,A2,A3,T0
C
C     Physical constants:
CC    PARAMETER(
C         Number of outgoing bonds per molecule:
CC   +        B=6.12E0,
C         Number of bonds per molecule:
CC   +        N=20.0E0,
C         Average lattice length (m):
CC   +        X=7.0848E-10,
C         Avogadros number (molD-1)
CC   +        NA=6.022045E+23,
C         Conversion, cal to Joule:
CC   +        C2=4.1868E0,
C         Coefficients for calculation of fusion heat of SAT:
CC   +        A1=7323.87E0, A2=73.1492E0, A3=0.062128E0, T0=244.88E0)
CC    PARAMETER(
CC   +        C1=0.5E0*B/(NA*N*X*X),
CC   +        C3=C1*C2)
C
CC    STSAST=C3*(A1+(A2-A3*(TAIR+T0))*(TAIR-T0))
CC    RETURN
CC    END
*****************************************************************************
*     REAL FUNCTION STSAST(TAIR,WSA)                                            
      REAL FUNCTION STSAST(TAIR,WSA)
*****************************************************************************
*
*     Surface tension of sulfuric acid tetrahydrate/liquid solution.
*
*     Source: Luo et al., GRL,1993
*
*
*     Input:  TAIR: Temperature (K)
*             WSA:  Weight fraction of H2SO4  [0;1]
*     Output: Surface tension (N/m)
*
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WSA,W
      W=WSA*100.0E0
      STSAST=(-59.41E0+0.6974E0*W+(0.3322E0-2.1855E-3*W)*TAIR)*1.0E-3
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION STNAS(TAIR,WNA)                                                 *
      REAL FUNCTION STNAS(TAIR,WNA)
*****************************************************************************
*
*     Surface tension of nitric acid solution/vapor.
*
*     Source: Granzhan and Laktionova, Russian J. Phys. Chem. 49, 1448, 1975
*
*
*     Input:  TAIR: Temperature (K)
*             WNA:  Weight fraction of HNO3  [0;1]
*     Output: Surface tension of nitric acid solution (N/m)
*
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WNA
C
      STNAS=(72.75E0-(5.57E0+24.55E0*WNA)*WNA-
     +       0.166E0*(TAIR-293.15E0))*1.0E-3
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION STSNW(TAIR,WSA,WNA)                                                 *
CC    REAL FUNCTION STSNW(TAIR,WSA,WNA)
*****************************************************************************
*
*     Surface tension of liquid TERNARY sulfuric/nitric acid solution.
*
*
*     Input:  TAIR: Temperature  (K)
*             WSA:  Weight fraction of H2SO4  [0;1] 
*             WNA:  Weight fraction of HNO3  [0;1] 
*     Output: Surface tension of sulfuric/nitric acid solution /air  (N/m)
*
C
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,WSA,WNA
CC    REAL(KIND=4)  SN,SS,MN,MS,STWL,STSAS,STNAS,MHNO3,MH2SO4,WN,WS
CC    PARAMETER(
C       Molar weight of nitric acid (kg/mole)
CC   +          MHNO3=63.01E-3,
C       Molar weight of sulfuric acid (kg/mole)
CC   +          MH2SO4=98.08E-3)
C
CC    IF(WSA.EQ.0.0.AND.WNA.EQ.0.0) THEN
CC        STSNW=STWL(TAIR)
CC    ELSE
CC        MN=WNA/(MHNO3*(1-WNA-WSA))
CC        MS=WSA/(MH2SO4*(1-WNA-WSA))
CC        WN=MHNO3*MN/(1.0E0+MHNO3*MN)
CC        WS=MH2SO4*MS/(1.0E0+MH2SO4*MS)
CC        SN=STNAS(TAIR,WN)
CC        SS=STSAS(TAIR,WS)
CC        STSNW=1.0E0/((1.0E0/SS)*MS/(MN+MS)+(1.0E0/SN)*MN/(MN+MS))
CC    ENDIF
C
CC    RETURN
CC    END
*****************************************************************************
*     REAL FUNCTION STSNW(TAIR,WSA,WNA)                                                 *
      REAL FUNCTION STSNW(TAIR,WSA,WNA)
*****************************************************************************
*
*     Surface tension of liquid TERNARY sulfuric/nitric acid solution.
*     Source: Martin et al. GRL 27, 197, 2000
*
*     Input:  TAIR: Temperature  (K)
*             WSA:  Weight fraction of H2SO4  [0;1] 
*             WNA:  Weight fraction of HNO3  [0;1] 
*     Output: Surface tension of sulfuric/nitric acid solution /air  (N/m)
*
      IMPLICIT NONE
      REAL TAIR,WSA,WNA,X,Y,YP,T293,T253,A0,A1,A2,A3,A4,A5,S293,S253,DT,
     +R00,R01,R02,R03,R04,R05,
     +R10,R14,
     +R20,R23,
     +R30,R32,
     +R40,R41,
     +R50,
     +R0,R1,R2,
     +S00,S01,S02,S03,S04,S05,
     +S10,S14,
     +S20,S23,
     +S30,S32,
     +S40,S41,
     +S50,
     +S0,S1,S2
      DATA
     +R00,R01,R02,R03,R04,R05,
     +R10,R14,
     +R20,R23,
     +R30,R32,
     +R40,R41,
     +R50,
     +R0,R1,R2/
     + 7.201E7,3.893E6,9.736E4,-1.832E3,1.282E1,1.076E-1,
     + 1.073E6,-2.811E-1,
     + 8.580E3, 3.358E-1,
     +-1.036E2,-1.866E-1,
     + 2.270E0, 4.895E-3,
     +-2.333E-2,
     +9.949E2,7.321E0,2.917E1/
      DATA
     +S00,S01,S02,S03,S04,S05,
     +S10,S14,
     +S20,S23,
     +S30,S32,
     +S40,S41,
     +S50,
     +S0,S1,S2/
     + 3.601E7,5.894E5,5.397E4,-2.994E3,6.919E1,-5.648E-1,
     + 8.127E5, 1.001E-2,
     + 2.194E4,-9.681E-2,
     +-4.554E2, 8.724E-2,
     + 7.115E0,-3.648E-2,
     +-4.483E-2,
     +6.726E2,9.692E0,8.276E0/
      PARAMETER(T293=293.0,T253=253.0,DT=T293-T253)
C
      X=100.0*WSA
      Y=100.0*AMIN1(WNA,0.45)
C
      A0=R00+Y*(R01+Y*(R02+Y*(R03+Y*(R04+Y*R05))))
      A5=R50
      YP=Y
      A4=R40+R41*YP
      YP=YP*Y
      A3=R30+R32*YP
      YP=YP*Y
      A2=R20+R23*YP
      YP=YP*Y
      A1=R10+R14*YP
      S293=(A0+X*(A1+X*(A2+X*(A3+X*(A4+X*A5)))))/((R0+R1*X+R2*Y)**2)
C
      A0=S00+Y*(S01+Y*(S02+Y*(S03+Y*(S04+Y*S05))))
      A5=S50
      YP=Y
      A4=S40+S41*YP
      YP=YP*Y
      A3=S30+S32*YP
      YP=YP*Y
      A2=S20+S23*YP
      YP=YP*Y
      A1=S10+S14*YP
      S253=(A0+X*(A1+X*(A2+X*(A3+X*(A4+X*A5)))))/((S0+S1*X+S2*Y)**2)
C
      STSNW=1.0E-3*(S253+(S293-S253)*(TAIR-T253)/DT)
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION STNAIC(TAIR,WNA)                                                 *
      REAL FUNCTION STNAIC(TAIR,WNA)
*****************************************************************************
*
*     Surface tension of nitric acid solution/water ice.
*
*     Application of Antonoff's rule: STSAIC=ABS(STNAS(TAIR,WNA)-STICE)    
*     Source: Granzhan and Laktionova, Russian J. Phys. Chem. 49, 1448, 1975
*
*
*     Input:  TAIR:   Temperature (K)
*             WNA:    Weight fraction of HNO3  [0;1]
*     Output: STNAIC: Surface tension of nitric acid solution/ice (N/m)
*
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WNA,STNAS,STICE
      PARAMETER(
C       Surface tension ice/air (N/m) (Pruppacher & Klett, 1980, p. 121)
     +          STICE=105.0E-3)
      STNAIC=ABS(STNAS(TAIR,WNA)-STICE)
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION STNAIC(TAIR,WNA)                                                 *
CC    REAL FUNCTION STNAIC(TAIR,WNA)
*****************************************************************************
*
*     Surface tension of nitric acid solution/water ice.
*     Source: Tabazadeh et al., GRL 24, 2007, 1997.
*     Application of Antonoff's rule: STSAIC=ABS(STNAS(TAIR,WNA)-STICE)    
*     Source: Granzhan and Laktionova, Russian J. Phys. Chem. 49, 1448, 1975
*
*
*     Input:  TAIR:   Temperature (K)
*             WNA:    Weight fraction of HNO3  [0;1]
*     Output: STNAIC: Surface tension of nitric acid solution/ice (N/m)
*
C
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,WNA,WN
CC    WN=WNA*100.0
CC    STNAIC=1.0E-3*(32.25+(0.002455*WN+0.0557)*WN+0.166*(TAIR-293.15))
CC    RETURN
CC    END
******************************************************************************
*     REAL FUNCTION VISSAS(TAIR,WSAS)
      REAL FUNCTION VISSAS(TAIR,WSAS)
*****************************************************************************
*
*     Dynamic viscosity of sulfuric acid liquid solution.
*
*     Source:   Williams & Long, J. Phys. Chem. 99, 3748, 1995 
*     
*
*     Input:  TAIR:  Temperature (K)
*             WSAS:  Weight fraction of H2SO4  
*
*     Output: Dynamic viscosity of sulfuric acid solution ( kg m-1 s-1 )
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WSAS
      REAL(KIND=4)  WP,A,B,N,T0
      PARAMETER(
C       Parameters from Williams and Long fits to viscosity:
     +          N=-1.43,B=425.0E0)
C
      WP=WSAS*100.0E0
      A=279.4E0+(0.358E0*WP-8.8E0)*WP
      T0=203.0E0+(0.0287E0*WP-2.63E0)*WP
      VISSAS=1.0E-3*A*(TAIR**N)*EXP(B/(TAIR-T0))
      RETURN
      END          
******************************************************************************
*     REAL FUNCTION VISSAS(TAIR,WSAS)
CC      REAL FUNCTION VISSAS(TAIR,WSAS)
*****************************************************************************
*
*     Dynamic viscosity of sulfuric acid liquid solution.
*
*     Source: 
*     Experimental data from 
*     K. Schaefer: Landolt-Bornstein, Serie II, Band 5. Teil, Springer 1969
*     have been fitted with functional expression
*            V = A T**n EXP(B/(Tair-T0))
*     cf. Eicher & Zwolinski, J. Phys. Chem. 75, 2016, 1971.
*     
*
*     Input:  TAIR: Temperature (K)
*             WSA:  Weight fraction of H2SO4  [0.37;0.69]
*     Note:   TAIR should be equilibrium temperature, corresponding to WSAS
*
*     Output: Dynamic viscosity of sulfuric acid solution ( kg m-1 s-1 )
*
CC      IMPLICIT NONE
CC      REAL(KIND=4)  TAIR,WSAS
CC      INTEGER NWEIGH,I
CC      PARAMETER(NWEIGH=33)
CC      REAL(KIND=4)  WEIGHT(NWEIGH),A(NWEIGH),N(NWEIGH),B(NWEIGH),T0(NWEIGH)
CC      REAL(KIND=4)  WORK(NWEIGH),A2(NWEIGH),N2(NWEIGH),B2(NWEIGH),T02(NWEIGH)
CC      REAL(KIND=4)  AA,NN,BB,TT0,W
CC      LOGICAL FIRST
C
CC      DATA (WEIGHT(I),I=1,NWEIGH)/
CC     +0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,
CC     +0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,
CC     +0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,
CC     +0.64,0.65,0.66,0.67,0.68,0.69/
CC      DATA (A(I),I=1,NWEIGH)/
CC     +-7.601333, -7.155090, -6.687569, -6.130699, -5.854452,
CC     +-5.047421, -4.255189, -3.121073, -1.927385, -0.313149,
CC     + 1.391153,  3.444294,  5.734247,  7.785439, 10.492723,
CC     +11.912110, 13.527665, 14.154839, 14.746930, 14.572778,
CC     +14.378593, 12.271519, 10.498063,  6.992785,  1.952391,
CC     +-3.786626, -9.959085,-16.583591,-22.466058,-27.052743,
CC     +-28.752641,-29.450797,-29.793597/
CC      DATA (N(I),I=1,NWEIGH)/
CC     + 0.8064693, 0.7394031, 0.6688089, 0.5834553, 0.5432935,
CC     + 0.4208850, 0.3007976, 0.1252286,-0.0598020,-0.3129303,
CC     +-0.5803307,-0.9039132,-1.2649952,-1.5838053,-2.0202958,
CC     +-2.2456562,-2.5031907,-2.5966581,-2.6859410,-2.6526778,
CC     +-2.6189834,-2.2810861,-2.0009938,-1.4405595,-0.6317888,
CC     + 0.2922928, 1.2919280, 2.3616670, 3.3124277, 4.0554005,
CC     + 4.3383722, 4.4608408, 4.5186140/
CC      DATA (B(I),I=1,NWEIGH)/
CC     +591.27376,586.11335,580.97608,575.23241,573.18849,
CC     +558.13985,543.42012,524.83788,505.49332,481.39790,
CC     +456.20943,426.87318,393.98142,359.17093,334.90806,
CC     +322.61725,308.60965,299.93797,293.08724,296.44784,
CC     +302.70727,339.69831,377.81086,436.58333,514.48840,
CC     +597.20183,676.39386,766.41693,845.59082,905.56668,
CC     +917.79973,915.17615,920.22190/
CC      DATA (T0(I),I=1,NWEIGH)/
CC     +142.30595,142.28857,142.27225,142.25480,142.24938,
CC     +143.18639,144.17291,145.19825,146.25613,147.34045,
CC     +148.43930,149.59745,150.93628,153.09946,151.93749,
CC     +150.79941,149.51516,149.48057,149.45524,149.46466,
CC     +149.48120,147.77524,145.95390,144.52354,143.39261,
CC     +143.08890,143.83096,143.94587,144.26027,144.90170,
CC     +146.66404,148.64822,149.70634/
CC      DATA FIRST/.TRUE./
CC      SAVE WEIGHT,A,N,B,T0,AA,NN,BB,TT0,FIRST
C
CC      IF(FIRST) THEN
CC          FIRST=.FALSE.
CC          CALL SPLINE(WEIGHT,A,NWEIGH,WORK,A2)
CC          CALL SPLINE(WEIGHT,N,NWEIGH,WORK,N2)
CC          CALL SPLINE(WEIGHT,B,NWEIGH,WORK,B2)
CC          CALL SPLINE(WEIGHT,T0,NWEIGH,WORK,T02)
CC      ENDIF
CC      W=AMAX1(WEIGHT(1),AMIN1(WSAS,WEIGHT(NWEIGH)))
CC      CALL SPLINT(WEIGHT,A,A2,NWEIGH,W,AA)
CC      CALL SPLINT(WEIGHT,N,N2,NWEIGH,W,NN)
CC      CALL SPLINT(WEIGHT,B,B2,NWEIGH,W,BB)
CC      CALL SPLINT(WEIGHT,T0,T02,NWEIGH,W,TT0)
C
CC      VISSAS=1.0E-3*EXP(AA+NN*ALOG(TAIR)+BB/(TAIR-TT0))
CC      RETURN
CC      END
*******************************************************************************
      REAL FUNCTION DIFESA(TAIR,WSAS)
*     REAL FUNCTION DIFESA(TAIR,WSAS)
*******************************************************************************
*
*     Activation energy for self-diffusion of H2SO4 in 
*     sulfuric acid solution
*
*     Source: 
*
*     The activation energy for self-diffusion is calculated from
*     the method of Luo et al. GRL,1993.
*
*     Experimental viscosity data from 
*     K. Schaefer: Landolt-Bornstein, Serie II, Band 5. Teil, Springer 1969
*     have been fitted with functional expression
*            V = A T**n EXP(B/(Tair-T0))
*     cf. Eicher & Zwolinski, J. Phys. Chem. 75, 2016, 1971.
*
*
*     Input:  TAIR: Temperature (K)
*     Note:   TAIR  must be the equilibrium temperature, 
*                   corresponding to WSAS
*             WSA:  Weight fraction of H2SO4  [0.37;0.69]
*
*     Output: 
*     Activation energy for self-diffusion of H2SO4 in sulfuric acid solution 
*                   (J/molecule)
*
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WSAS
      REAL(KIND=4)  NAVOGA,KCALJ,DGDIS,DH2OB,HH2O,HH2SO4,T25,
     +     C1,C2
      PARAMETER(
C       Avogadro's number (1/mole)
     +          NAVOGA=6.022045E+23,
C       Conversion factor, kcal to J:
CC     +          KCALJ=4.1868E+3,
     +          KCALJ=4.1840E+3,
C       Displacement component of H2SO4 diffusion energy (kcal/mole)
     +          DGDIS=2.07E0,
C       H2O bonding partial energy (kcal/mole):
CC     +          DH2OB=DH2O-DH2OD,
     +          DH2OB=2.7998E0,
C       Partial vaporization energy of H2O at 25 C (kcal/mole)
     +          HH2O=10.52E0,
C       Partial vaporization energy of H2SO4 at 25 C (kcal/mole)
     +          HH2SO4=18.90E0,
C
     +          T25=25.0E0+273.15E0,
     +          C1=DH2OB/HH2O, C2=KCALJ/NAVOGA)
C
      REAL(KIND=4)  LSAS,F,DH2SO4
C
      INTEGER NWEIGH,I
      PARAMETER(NWEIGH=33)
      REAL(KIND=4)  WEIGHT(NWEIGH),N(NWEIGH),B(NWEIGH),T0(NWEIGH)
      REAL(KIND=4)  WORK(NWEIGH),N2(NWEIGH),B2(NWEIGH),T02(NWEIGH)
      REAL(KIND=4)  NN,BB,TT0,WW
      COMMON/FF/NN,BB,TT0,WW,F,DH2SO4
      LOGICAL FIRST
C
      DATA (WEIGHT(I),I=1,NWEIGH)/
     +        0.37,0.38,0.39,
     +        0.40,0.41,0.42,0.43,0.44,0.45,
     +        0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,
     +        0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,
     +        0.64,0.65,0.66,0.67,0.68,0.69/
      DATA (N(I),I=1,NWEIGH)/
     + 0.8064693, 0.7394031, 0.6688089, 0.5834553, 0.5432935,
     + 0.4208850, 0.3007976, 0.1252286,-0.0598020,-0.3129303,
     +-0.5803307,-0.9039132,-1.2649952,-1.5838053,-2.0202958,
     +-2.2456562,-2.5031907,-2.5966581,-2.6859410,-2.6526778,
     +-2.6189834,-2.2810861,-2.0009938,-1.4405595,-0.6317888,
     + 0.2922928, 1.2919280, 2.3616670, 3.3124277, 4.0554005,
     + 4.3383722, 4.4608408, 4.5186140/
      DATA (B(I),I=1,NWEIGH)/
     +591.27376,586.11335,580.97608,575.23241,573.18849,
     +558.13985,543.42012,524.83788,505.49332,481.39790,
     +456.20943,426.87318,393.98142,359.17093,334.90806,
     +322.61725,308.60965,299.93797,293.08724,296.44784,
     +302.70727,339.69831,377.81086,436.58333,514.48840,
     +597.20183,676.39386,766.41693,845.59082,905.56668,
     +917.79973,915.17615,920.22190/
      DATA (T0(I),I=1,NWEIGH)/
     +142.30595,142.28857,142.27225,142.25480,142.24938,
     +143.18639,144.17291,145.19825,146.25613,147.34045,
     +148.43930,149.59745,150.93628,153.09946,151.93749,
     +150.79941,149.51516,149.48057,149.45524,149.46466,
     +149.48120,147.77524,145.95390,144.52354,143.39261,
     +143.08890,143.83096,143.94587,144.26027,144.90170,
     +146.66404,148.64822,149.70634/
      DATA FIRST/.TRUE./
      SAVE FIRST,WEIGHT,N,N2,B,B2,T0,T02
C
      INTEGER TEST
      COMMON /CONTRL/TEST
C
      IF(FIRST) THEN
          FIRST=.FALSE.
          CALL SPLINE(WEIGHT,N,NWEIGH,WORK,N2)
          CALL SPLINE(WEIGHT,B,NWEIGH,WORK,B2)
          CALL SPLINE(WEIGHT,T0,NWEIGH,WORK,T02)
      ENDIF
C
      WW=AMAX1(WEIGHT(1),AMIN1(WSAS,WEIGHT(NWEIGH)))
C
      CALL SPLINT(WEIGHT,N,N2,NWEIGH,WW,NN)
      CALL SPLINT(WEIGHT,B,B2,NWEIGH,WW,BB)
      CALL SPLINT(WEIGHT,T0,T02,NWEIGH,WW,TT0)
      F=(TAIR*(BB*TAIR/((TAIR-TT0)**2)+1.0E0-NN))/
     +        (T25*(BB*T25/((T25-TT0)**2)+1.0E0-NN))
C
      DH2SO4=HH2SO4+LSAS(WW)*1.0E-3
      DIFESA=(C1*DH2SO4+DGDIS)*F*C2
C
      RETURN
      END
*******************************************************************************
      REAL FUNCTION DIFEWL(TAIR)
*     REAL FUNCTION DIFEWL(TAIR)
*******************************************************************************
*
*     Activation energy for H2O self-diffusion in water.
*     Source: Prupacher & Klett:Microphysics of clouds and precipitation,
*                               (1980), p. 178)
*             Pruppacher, JAS 52, 1924, 1995.
*         OR:
*             Jensen et al. GRL 18,1857,1991
*             Jensen et al. JGR 99, 10421, 1994:
*
*     Input:  TAIR:       Temperature (K)
*
*     Output: DIFEWL:     Activation energy for self-diffusion in water (J)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR
C     Pruppacher & Klett expression and data:  NOT TESTED !!!
CC    INTEGER NDATA
CC    PARAMETER(NDATA=13)
CC    REAL(KIND=4)  T,C1,C2,TTAB(NDATA),DTAB(NDATA),YWORK(NDATA),Y2(NDATA),X
CC    PARAMETER(C2=6.9525E-21,C1=C2*5.5E0)
CC    LOGICAL FIRST
CC    DATA FIRST/.TRUE./
CC    DATA TTAB/242.15,241.15,240.15,239.15,238.15,237.15,236.15,235.15,
CC   +          233.15,232.15,231.15,230.15,229.15/
CC    DATA DTAB/
CC   +          10.00,9.60,9.20,8.60,7.90,7.75,7.65,7.45,
CC   +          7.25,6.85,6.00,5.15,4.45/
CC    SAVE FIRST,TTAB,DTAB,Y2
C
CC    IF(FIRST) THEN
CC        FIRST=.FALSE.
CC        CALL SPLINE(TTAB,DTAB,NDATA,YWORK,Y2)
CC    ENDIF
C
CC    IF(TAIR.GT.242.15E0) THEN
CC        T=TAIR-273.15E0
CC        DIFEWL=C1*EXP((-1.330E-2+(2.74E-4+1.085E-6*T)*T)*T)
CC    ELSE
CC        T=AMAX1(229.15E0,TAIR)
CC        CALL SPLINT(TTAB,DTAB,Y2,NDATA,T,X)
CC        DIFEWL=C2*X
CC    ENDIF
C
C     Expression from Jensen et al. GRL 18,1857,1991:
CC    IF(TAIR.GT.228.15E0) THEN
CC        DIFEWL=(-71.18E-13+0.3400E-13E0*TAIR)*1.0E-7
CC    ELSE
CC        DIFEWL=6.391E-20
CC    ENDIF
C
C     Expression from Jensen et al. JGR 99, 10421, 1994:
      DIFEWL=1.29E-19*EXP(0.05E0*(AMAX1(TAIR,195.0E0)-237.15E0))
      RETURN
      END
*******************************************************************************
      REAL FUNCTION DIFEWS(TAIR,WSAS)
*     REAL FUNCTION DIFEWS(TAIR,WSAS)
*******************************************************************************
*
*     Activation energy for H2O self-diffusion in sulfuric acid solution.
*     Experimental viscosity data from  
*     Williams & Long, J. Phys. Chem. 99, 3748, 1995 
*     have been been used together with functional expression
*            DIFEWS = kT(BT/(T-T0)+1-n)
*
*     Input:  TAIR:   Temperature (K)
*             WSAS:   Sulfuric acid weight fraction
*     Output: DIFEWS: Activation energy for self-diffusion of H2O 
*                     in sulfuric acid solution (J)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WSAS
      REAL(KIND=4)  WP,N,B,T0,KBOLTZ
      PARAMETER(
C       Parameters from Williams and Long fits to viscosity:
     +          N=-1.43E0,B=425.0E0,
C       Boltzmanns constant (J/K)
     +          KBOLTZ=1.380662E-23)
C
      WP=WSAS*100.0E0
      T0=203.0E0+(0.0287E0*WP-2.63E0)*WP
      DIFEWS=KBOLTZ*TAIR*((B*TAIR/((TAIR-T0)**2))+1.0E0-N)
      RETURN
      END  
*******************************************************************************
      REAL FUNCTION DIFEAT(TAIR)
*     REAL FUNCTION DIFEAT(TAIR)
*******************************************************************************
*
*     Activation energy for H2O self-diffusion in sulfuric acid solution.
*
*     Input:  TAIR:   Temperature (K)
*
*     Output: DIFEWS: Activation energy for self-diffusion of H2O 
*                     in sulfuric acid solution (J)
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR
      DOUBLE PRECISION T,X
*     Source: Tabazadeh et al., JGR,102, 23845,1997
      DOUBLE PRECISION A0,A1,A2,A3,A4,A5,A6
CC    DATA A0,A1,A2,A3,A4,A5,A6/
CC   +1.1992475202E+4,-5.2298196745E+2,8.232842046E0,-6.4173902638E-2,
CC   +2.6889968134E-4,-5.8279763451E-7,5.1479983159E-10/
CC    SAVE A0,A1,A2,A3,A4,A5,A6
C
*     Source: Tabazadeh et al., 27, 1111, 2000
      DATA A0,A1,A2,A3,A4,A5,A6/
     +-17459.516183D0,458.45827551D0,-4.8492831317D0,0.026003658878D0,
     +-7.1991577798D-5,8.9049094618D-8,-2.4932257419D-11/
      SAVE A0,A1,A2,A3,A4,A5,A6
      DOUBLE PRECISION B0,B1,B2,B3,B4,B5,B6
      DATA B0,B1,B2,B3,B4,B5,B6/
     +104525.93058D0,-1103.7644651D0,1.070332702D0,0.017386254322D0,
     +-1.5506854268D-6,-3.2661912497D-8,6.467954459D-10/
      SAVE B0,B1,B2,B3,B4,B5,B6
C
      T=DBLE(TAIR)
      IF(T.LE.220.0D0) THEN
          X=(((((A6*T+A5)*T+A4)*T+A3)*T+A2)*T+A1)*T+A0
      ELSE IF(T.GT.220.0D0) THEN
          X=(((((B6*T+B5)*T+B4)*T+B3)*T+B2)*T+B1)*T+B0
      ENDIF
      DIFEAT=1.0E-20*REAL(X)
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION PWVWL(TAIR)
      REAL FUNCTION PWVWL(TAIR)
*****************************************************************************
*
*     Saturation pressure of water vapor over plane water surface.
*
*     Source: Tabazadeh et al. GRL 24, 1931, 1997
*
*     Input:  Temperature (K);   Range:  [185,260]
*     Output: Sat. pressure (Pa)
*
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR
      DOUBLE PRECISION C1,C2,C3,C4,TD,P
      DATA C1,C2,C3,C4/
     +    18.452406985D0,-3505.1578807D0,-330918.55082D0,12725068.262D0/

      TD=1.0D0/DBLE(TAIR)
      P=C1+(C2+(C3+C4*TD)*TD)*TD
      PWVWL=100.0E0*EXP(P)

      RETURN
      END
*****************************************************************************
*     REAL FUNCTION PWVWL(TAIR)
cc     REAL FUNCTION PWVWL(TAIR)
*****************************************************************************
*
*     Saturation pressure of water vapor over plane water surface.
*
*     Source: F.W. Murray: J. Appl. Meteor. 6, 203,1967
*
*     Input:  Temperature (K);   Range:  [175    ,265    ]
*     Output: Sat. pressure (Pa)
*
cc
CC      IMPLICIT NONE
cc      REAL(KIND=4)  TAIR,XT
cc      REAL(KIND=4)  C1,C2,C3,C4,C5,C6,C7,TS
cc      DATA C1,C2,C3,C4,C5,C6,C7/
cc     +        7.95357242E12,-18.1972839,5.02808,-70242.1852,-26.1205253,
cc     +        58.0691913,-8.03945282/
cc      DATA TS/373.16/
cc
cc      XT=TS/TAIR
cc      PWVWL=C1*EXP(C2*XT+C3*ALOG(XT)+C4*EXP(C5/XT)+C6*EXP(C7*XT))
cc
cc      RETURN
cc      END
*****************************************************************************
*     REAL FUNCTION PWVICE(TAIR)
      REAL FUNCTION PWVICE(TAIR)
*****************************************************************************
*
*     Saturation pressure of water vapor over plane water surface of ice.
*
*     Source: J. Marti and K. Mauersberger, GRL. 20,363,1993
*     (Source: G. Jancso et al.: J. Phys. Chem. 74,2984,1970)
*     (Source: Tabazadeh et al., JGR 102,23845,1997
*
*     Input:  Temperature (K);   Range:  [173    ,273    ]
*     Output: Sat. pressure (Pa)
*
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR
C
C     Marti's and Mauersberger's expression:
      REAL(KIND=4)  A,B
      PARAMETER(A=-2663.5, B=12.537)
C
      PWVICE=10.0**(A/TAIR+B)
C
C     Jancso's expression:
CC    REAL(KIND=4)  C1,C2,C3,C4,C5,C6
CC    DATA C1,C2,C3,C4,C5,C6/-2481.604,3.5721988,-3.097203E-3,
CC   +                       -1.7649E-7,1.901973,133.32237/
C
CC    PWVICE=C6*10.0**(C1/TAIR+C2*ALOG10(TAIR)+(C3+C4*TAIR)*TAIR+C5)
C
C     Tabazadeh's expression
CC    PWVICE=100.0*EXP(24.313-(6146.8/TAIR))
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION PNAMON(TAIR)
      REAL FUNCTION PNAMON(TAIR)
*****************************************************************************
*     Nitric acid vapor pressure over nitric acid mono/trihydrate
*     Source: Hanson & Mauersberger: J. Phys. Chem. 92, 6167, 1988.
*
*     Input:  TAIR         Temperature (K)
*     Output: PNAMON       Nitric acid pressure (Pa)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR
      PNAMON=133.32237E0*10.0**(13.622E0-3561.3E0/TAIR)
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION PWVMON(TAIR)
      REAL FUNCTION PWVMON(TAIR)
*****************************************************************************
*     Water vapor vapor pressure over nitric acid mono/trihydrate
*     Source: Hanson & Mauersberger: J. Phys. Chem. 92, 6167, 1988.
*
*     Input:  TAIR         Temperature (K)
*     Output: PWVMON       Water vapor pressure (Pa)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR
      PWVMON=133.32237E0*10.0**(10.049E0-2819.2E0/TAIR)
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION PNATRI(TAIR)
      REAL FUNCTION PNATRI(TAIR)
*****************************************************************************
*     Nitric acid vapor pressure over nitric acid ice/trihydrate
*     Source: Hanson & Mauersberger: J. Phys. Chem. 92, 6167, 1988.
*
*     Input:  TAIR         Temperature (K)
*     Output: PNATRI       Nitric acid pressure (Pa)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR
      PNATRI=133.32237E0*10.0**(12.298E0-3968.0E0/TAIR)
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION PWVTRI(TAIR)
      REAL FUNCTION PWVTRI(TAIR)
*****************************************************************************
*     Water vapor vapor pressure over nitric acid ice/trihydrate
*     Source: Hanson & Mauersberger: J. Phys. Chem. 92, 6167, 1988.
*
*     Input:  TAIR         Temperature (K)
*     Output: PWVTRI       Water vapor pressure (Pa)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR
      PWVTRI=133.32237E0*10.0**(10.431E0-2668.7E0/TAIR)
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION PNANAT(TAIR,PWV)
      REAL FUNCTION PNANAT(TAIR,PWV)
*****************************************************************************
*     Nitric acid vapor vapor pressure over nitric acid trihydrate/solid
*     solution.
*     Source: Hanson & Mauersberger: Geophys. Res. Lett. 15, 855, 1988.
*
*     Input:  TAIR         Temperature (K)
*             PWV          Water vapor saturation pressure (Pa)
*     Output: PNANAT       Nitric acid vapor pressure (Pa)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PWV
      PNANAT=133.32237*10.0**(
     +       AMAX1(12.298E0-3968.0E0/TAIR,
     +       AMIN1((-2.7836E0-TAIR*8.8E-4)*(ALOG10(PWV)- 2.124903090E0)+
     +                     (38.9855E0-(11397.0E0/TAIR)+(9.179E-3*TAIR)),
     +             13.622E0-3561.3E0/TAIR)))
      RETURN
      END
*
*****************************************************************************
*     REAL FUNCTION PNANAT(TAIR,PWV)
CC    REAL FUNCTION PNANAT(TAIR,PWV)
*****************************************************************************
*     Nitric acid vapor vapor pressure over nitric acid trihydrate.
*     Source: Worsnop et al. Science 259, 71, 1993.
*
*     Input:  TAIR         Temperature (K)
*             PWV          Water vapor saturation pressure (Pa)
*     Output: PNANAT       Nitric acid vapor pressure (Pa)
*
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,PWV
CC    REAL(KIND=4)  N,LN760,DH0_R,DS0_R,TORR_TO_PA,NLN
CC    PARAMETER(
CC   +     N=3.0,
CC   +     DH0_R=28520.0,
CC   +     DS0_R=78.7,
CC   +     LN760=6.63331843,
CC   +     TORR_TO_PA=133.32237,
CC   +     NLN=(N+1.0)*LN760)
CC    PNANAT=TORR_TO_PA*
CC   +       EXP((-DH0_R/TAIR)+
CC   +           DS0_R+NLN-N*ALOG(PWV/TORR_TO_PA))
CC    RETURN
CC    END
*
*****************************************************************************
*     REAL FUNCTION PNANAD(TAIR,PWV)
      REAL FUNCTION PNANAD(TAIR,PWV)
*****************************************************************************
*     Nitric acid vapor vapor pressure over nitric acid dihydrate.
*     Source: Worsnop et al. Science 259, 71, 1993.
*
*     Input:  TAIR         Temperature (K)
*             PWV          Water vapor saturation pressure (Pa)
*     Output: PNANAD       Nitric acid vapor pressure (Pa)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PWV
      REAL(KIND=4)  N,LN760,DH0_R,DS0_R,TORR_TO_PA,NLN
      PARAMETER(
     +     N=2.0,
     +     DH0_R=22630.0,
     +     DS0_R=64.9,
     +     LN760=6.63331843,
     +     TORR_TO_PA=133.32237,
     +     NLN=(N+1.0)*LN760)
      PNANAD=TORR_TO_PA*
     +       EXP((-DH0_R/TAIR)+
     +           DS0_R+NLN-N*ALOG(PWV/TORR_TO_PA))
      RETURN
      END
*
*****************************************************************************
*     REAL FUNCTION PNANAS(TAIR,PPWV)
      REAL FUNCTION PNANAS(TAIR,PPWV)
*****************************************************************************
*     Nitric acid vapor vapor pressure over nitric acid liquid solution.
*     The expression logP = A-B/T have been extrapolated linearly.
*     Source: Hanson: GRL 17, 421, 1990.
*
*     Input:  TAIR         Temperature (K)
*             PPWV:        Partial pressure of water vapor (Pa) 
*     Output: PNANAS       Nitric acid vapor pressure (Pa)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PPWV
      REAL(KIND=4)  PWV,TORPAS,X,
     +        WALFAA,WBETAA,WALFAB,WBETAB,
     +        NALFAA,NBETAA,NALFAB,NBETAB,
     +        XMIN,XMAX,FUNC,F,FMID,X1,X2,DX,XMID,XX,PW,T,
     +        ALFAA,BETAA,ALFAB,BETAB,XACC
      INTEGER J,JMAX
C
      PARAMETER(
C         Conversion factor Torr to Pascal:
     +    TORPAS=133.32237E0)
C
      DATA WALFAA,WBETAA,WALFAB,WBETAB
     +        /2.51714,8.55148,1.27143E+3,2.13643E+3/
      DATA NALFAA,NBETAA,NALFAB,NBETAB
     +        /-8.48571,12.6133,-4.28857E+3,4.23910E+3/
      DATA XMIN,XMAX,XACC/0.0,1.0,0.001/
      DATA JMAX/40/
C
      FUNC(XX,PW,T,ALFAA,BETAA,ALFAB,BETAB)=
     +        PW-((ALFAA*XX+BETAA)-(ALFAB*XX+BETAB)/T)
C
      PWV=ALOG10(PPWV/TORPAS)
C
      X1=XMIN
      X2=XMAX
C
      FMID=FUNC(X2,PWV,TAIR,WALFAA,WBETAA,WALFAB,WBETAB)
      F=FUNC(X1,PWV,TAIR,WALFAA,WBETAA,WALFAB,WBETAB)
      IF(F.LT.0.)THEN
        X=X1
        DX=X2-X1
      ELSE
        X=X2
        DX=X1-X2
      ENDIF
      DO J=1,JMAX
        DX=DX*.5
        XMID=X+DX
        FMID=FUNC(XMID,PWV,TAIR,WALFAA,WBETAA,WALFAB,WBETAB)
        IF(FMID.LE.0.) X=XMID
        IF(ABS(DX).LT.XACC.OR.FMID.EQ.0.) GOTO 12
      ENDDO
12    CONTINUE
      PNANAS=TORPAS*10.0**((NALFAA*X+NBETAA)-(NALFAB*X+NBETAB)/TAIR)
      RETURN
C
      END
*****************************************************************************
*     REAL FUNCTION PWVSAS(TAIR,WSA)                                        *
CC    REAL FUNCTION PWVSAS(TAIR,WSA)
*****************************************************************************
*
*     Natural logaritm of saturated water vapor pressure over plane
*     sulfuric acid solution.
*
*     Source: Jaecker-Voirol et al.: JGR 95,11857,1990,
*             W.F.Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.
*
*     The formula of Jaecker-Voirol et al. for saturation pressure
*     is used and calculated from partial molal properties given by 
*     Giauque et al.
*     The partial molal properties as function of sulfuric acid
*     mass fraction have been fitted with composit polynomial
*     expressions.
*     
*     
*
*     Input:  TAIR: Temperature (K)
*             WSA:  Weight fraction of H2SO4  [0;0.95] 
*     Output: Natural logaritm of water vapor pressure 
*             over sulfuric acid solution   ( ln(Pa) )
*
*
*     External functions needed for calculation of partial molal 
*     properties of pure components at 25 C as function of W:
CC    REAL(KIND=4)  CPH2O,FFH2O,LH2O,PWVWL
*     CPH2O:  Partial molal heat capacity of sulfuric acid solution.
*     FFH2O:  Partial molal free energy of sulfuric acid solution.
*     LH2O:   Partial molal enthalpy of sulfuric acid
*
C 
C 
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,WSA
CC    REAL(KIND=4)  DCP,F,L,LGPH2O
CC    REAL(KIND=4)  RGAS,T0,T0INV,C1,CPH2O0
CC    PARAMETER(
C     Gas constant (cal/(deg mole)):
CC   +        RGAS=1.98726,
CC   +        T0=298.15,
CC   +        T0INV=1.0/T0,
CC   +        CPH2O0=17.996,
CC   +        C1=1.0/(RGAS*T0))
C 
C 
CC    DCP=CPH2O0-CPH2O(WSA)
CC    F=FFH2O(WSA)
CC    L=LH2O(WSA)
CC    LGPH2O=ALOG(PWVWL(TAIR))
C
CC    PWVSAS=(-F*C1)+
CC   +        (-L-DCP*T0)*(1.0/TAIR-T0INV)/RGAS+
CC   +        DCP*ALOG(T0/TAIR)/RGAS+
CC   +        LGPH2O
CC    RETURN 
CC    END
C 
*****************************************************************************
*     REAL FUNCTION PWVSAS_GV(TAIR,WSA)                                        
      REAL FUNCTION PWVSAS_GV(TAIR,WSA)
*****************************************************************************
*
*     Natural logaritm of saturated water vapor pressure over plane
*     sulfuric acid solution.
*
*     Source: J.I.Gmitro & T.Vermeulen: A.I.Ch.E.J.  10,740,1964.
*             W.F.Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.
*
*     The formula of Gmitro & Vermeulen for saturation pressure
*     is used:
*                 ln(p) = A ln(298/T) + B/T + C + DT
*     with values of A,B,C and D given by Gmitro & Vermeulen,
*     and calculated from partial molal properties given by Giauque et al.
*     
*     
*
*     Input:  TAIR: Temperature (K)
*             WSA:  Weight fraction of H2SO4  [0;1] 
*     Output: Natural logaritm of water vapor pressure 
*             over sulfuric acid solution   ( ln(Pa) )
*
*
*     External functions needed for calculation of partial molal 
*     properties of pure components at 25 C as function of W.
      IMPLICIT NONE
      REAL(KIND=4)  CPH2O,ALH2O,FFH2O,LH2O
*     CPH2O:  Partial molal heat capacity of sulfuric acid solution.
*     ALH2O:  Temparature derivative of CPH2O
*     FFH2O:  Partial molal free energy of sulfuric acid solution.
*     LH2O:   Partial molal enthalpy of sulfuric acid
*
C
C
      REAL(KIND=4)  TAIR,WSA
      REAL(KIND=4)  ADOT,BDOT,CDOT,DDOT
      REAL(KIND=4)  RGAS,MMHGPA
      REAL(KIND=4)  K1,K2
      REAL(KIND=4)  A,B,C,D,CP,L,F,ALFA
C     Physical constants given by Gmitro & Vermeulen:
      PARAMETER(
     +        ADOT=-3.67340,
     +        BDOT=-4143.5,
     +        CDOT=10.24353,
     +        DDOT=0.618943E-3)
      PARAMETER(
C     Gas constant (cal/(deg mole)):
     +        RGAS=1.98726,
C     Natural logarith of conversion factor between atm. and Pa:     
     +     MMHGPA=11.52608845, 
     +     K1=298.15,
     +     K2=K1*K1/2.0)
C
C
      CP=CPH2O(WSA)
      F=-FFH2O(WSA)
      L=-LH2O(WSA)
      ALFA=ALH2O(WSA)
C
      A=ADOT+(CP-K1*ALFA)/RGAS
      B=BDOT+(L-K1*CP+K2*ALFA)/RGAS
      C=CDOT+(CP+(F-L)/K1)/RGAS
      D=DDOT-ALFA/(2.0E0*RGAS)
C
C     WRITE(*,*) 'TAIR= ',TAIR,'  WSA= ',WSA
C     WRITE(*,*) 'CPH2O(WSA)= ',CP
C     WRITE(*,*) 'ALFAH2O(WSA)= ',ALFA
C     WRITE(*,*) 'FFH2O(WSA)= ',F
C     WRITE(*,*) 'LH2O(WSA)= ',L
C
      PWVSAS_GV=A*ALOG(K1/TAIR)+B/TAIR+C+D*TAIR+MMHGPA
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION PWVSAS(TAIR,WSA)                                        
CC    REAL FUNCTION PWVSAS(TAIR,WSA)
*****************************************************************************
*
*     Natural logaritm of saturated water vapor pressure over plane
*     sulfuric acid solution.
*
*     Source: D.R. Hanson, A.R. Ravishankara & S. Solomon, JGR 99,3615,1994
*             Numerical fit to data Steele and Hamill (1981)
*
*     
*     
*
*     Input:  TAIR: Temperature (K)
*             WSA:  Weight fraction of H2SO4  [0;1] 
*     Output: Natural logaritm of water vapor pressure 
*             over sulfuric acid solution   ( ln(Pa) )
*
*
*
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,WSA,WW
CC    REAL(KIND=4)  C1,
CC   +     A0,A1,A2,
CC   +     B0,B1,B2
CC    PARAMETER(
CC   +     A0=-14.458,A1=0.19988,A2=1.62456,
CC   +     B0=3565.0,B1=-44.777,B2=-1.3204)
C
C     Natural log. of conversion factor between mbar and Pa:
CC    PARAMETER(C1=4.60517019)
C
CC    WW=WSA*100.0
C
CC    PWVSAS=(A0+A1*WW+(B0+B1*WW)/TAIR)/(1.0E0-A2-B2*WW/TAIR)+C1
CC    RETURN
CC    END
*****************************************************************************
*     REAL FUNCTION PWVSAS(TAIR,WSA)                                        
CC    REAL FUNCTION PWVSAS(TAIR,WSA)
*****************************************************************************
*
*     Natural logaritm of saturated water vapor pressure over plane
*     sulfuric acid solution.
*
*     Source: Zhang et al. J. Phys. Chem. 97, 7351, 1993
*     
*     
*
*     Input:  TAIR: Temperature (K)
*             WSA:  Weight fraction of H2SO4  [0;1] 
*     Output: Natural logaritm of water vapor pressure 
*             over sulfuric acid solution   ( ln(Pa) )
*
*
*
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,WSA
CC    REAL(KIND=4)  A(4),AA(4),BB(4),C1,C2,TINV
CC    INTEGER I
CC    DATA AA/9.60,4.04,-16.82,16.03/
CC    DATA BB/-2.44E+3,-0.77E+3,3.19E+3,-4.20E+3/
CC    PARAMETER(C1=2.12490303,C2=2.30258509)
CC    SAVE AA,BB
C
CC    TINV=1.0E0/TAIR
CC    DO 100 I=1,4
CC        A(I)=AA(I)+BB(I)*TINV
CC100 CONTINUE
CC    PWVSAS=(((A(4)*WSA+A(3))*WSA+A(2))*WSA+A(1)+C1)*C2
CC    RETURN 
CC    END
*****************************************************************************
*     REAL FUNCTION PWVSAS(TAIR,WSA)                                        
      REAL FUNCTION PWVSAS(TAIR,WSA)
*****************************************************************************
*
*     Natural logaritm of saturated water vapor pressure over plane
*     sulfuric acid solution.
*
*     Source: Tabazadeh et al. GRL, 24,1931, 1997 
*     
*     
*
*     Input:  TAIR: Temperature (K)
*             WSA:  Weight fraction of H2SO4  [0;0.75] 
*     Output: Natural logaritm of water vapor pressure 
*             over sulfuric acid solution   ( ln(Pa) )
*
*
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WSA
      INTEGER NWEI
      PARAMETER(NWEI=15)
      REAL(KIND=4)  WEIGHT(NWEI),AA(NWEI),BB(NWEI),CC(NWEI),
     +     Y2AA(NWEI),Y2BB(NWEI),Y2CC(NWEI),YWORK(NWEI)
      REAL(KIND=4)  A,B,C,DT,LN100,W
      PARAMETER(LN100=4.60517019)
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      DATA WEIGHT/
     +    0.10,0.15,0.20,0.25,0.30,
     +    0.35,0.40,0.45,0.50,0.55,
     +    0.60,0.65,0.70,0.75,0.80/
      DATA AA/
     +    19.726,19.747,19.761,19.794,19.883,
     +    20.078,20.379,20.637,20.682,20.555,
     +    20.405,20.383,20.585,21.169,21.808/
      DATA BB/
     +    -4364.8,-4390.9,-4414.7,-4451.1,-4519.2,
     +    -4644.0,-4828.5,-5011.5,-5121.3,-5177.6,
     +    -5252.1,-5422.4,-5743.8,-6310.6,-6985.9/
      DATA CC/
     +    -147620.,-144690.,-142940.,-140870.,-136500.,
     +    -127240.,-112550., -98811., -94033., -96984.,
     +    -100840., -97966., -83701., -48396., -12170./
      SAVE FIRST,WEIGHT,AA,BB,CC,Y2AA,Y2BB,Y2CC
C
      IF(FIRST) THEN
          FIRST=.FALSE.
          CALL SPLINE(WEIGHT,AA,NWEI,YWORK,Y2AA)
          CALL SPLINE(WEIGHT,BB,NWEI,YWORK,Y2BB)
          CALL SPLINE(WEIGHT,CC,NWEI,YWORK,Y2CC)
      ENDIF
      W=AMIN1(WSA,0.75E0)
      CALL SPLINT(WEIGHT,AA,Y2AA,NWEI,W,A)
      CALL SPLINT(WEIGHT,BB,Y2BB,NWEI,W,B)
      CALL SPLINT(WEIGHT,CC,Y2CC,NWEI,W,C)
      DT=1.0E0/TAIR
      PWVSAS=(C*DT+B)*DT+A+LN100
      RETURN
      END
*****************************************************************************
*     REAL FUNCTION PSASAS(TAIR,WSA)                                        
      REAL FUNCTION PSASAS(TAIR,WSA)
*****************************************************************************
*
*     Natural logaritm of saturated sulfuric acid vapor pressure 
*     over plane sulfuric acid solution.
*
*     Source: G.P. Ayers et al: GRL 7,433,1980
*             W.F.Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.
*
*     The formula of Ayers for saturation pressure is used:
*                 ln(p) = -10156/T + 16.259 + (my-my0)/8.3143 T
*     with values of (my-my0) given by Giauque et al.
*     
*     
*
*     Input:  TAIR: Temperature  (K)
*             WSA:  Weight fraction of H2SO4  [0;1]
*     Output: Natural logaritm of H2SO4 vapor pressure 
*             over sulfuric acid solution    ( ln(Pa) )
*
*
*     External functions needed for calculation of partial molal 
*     properties of pure components at 25 C as function of W.
      IMPLICIT NONE
      REAL(KIND=4)  FFSAS
*     FFSAS:  Partial molal free energy of sulfuric acid solution.
*
C
C
      REAL(KIND=4)  TAIR,WSA
      REAL(KIND=4)  ADOT,BDOT,CDOT
      REAL(KIND=4)  MMHGPA
      REAL(KIND=4)  F
C     Physical constants given by Ayers:
      PARAMETER(
     +        ADOT=-10156.0,
     +        BDOT=16.259,
     +        CDOT=4.1868/8.31441)
      PARAMETER(
C     Natural logarith of conversion factor between atm. and Pa:     
     +     MMHGPA=11.52608845)
C
      F=FFSAS(WSA)
C
C     WRITE(*,*) 'FFSAS(WSA)= ',F
C
      PSASAS=(ADOT-F*CDOT)/TAIR+BDOT+MMHGPA
      RETURN 
      END
C
*****************************************************************************
*     REAL FUNCTION PSASAT(TAIR,PPWV)                                        
      REAL FUNCTION PSASAT(TAIR,PPWV)
*****************************************************************************
*
*     Natural logaritm of saturated sulfuric acid vapor pressure 
*     over sulfuric acid tetrahydrate (SAT)
*
*     Source: J.I.Gmitro & T.Vermeulen: A.I.Ch.E.J.  10,740,1964.
*             G.P. Ayers et al: GRL 7,433,1980
*             W.F.Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.
*             C.M. Gable et al.: J. Aer. Chem. Soc. 72,1445,1950.
*
*     Input:  TAIR:   Temperature  (K) [TAIR < 244.8 K]
*             PPWV:   Partial pressure of water vapor (Pa)
*     Output: Natural logaritm of H2SO4 vapor pressure  ( ln(Pa) )
*
C
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,PPWV
C
C     Coefficients for slope and intercepts:
C     (Coefficinets produced by separate programme)
      REAL(KIND=4)  A,B,M
      PARAMETER(A=-4.65956,B=2.50899E-3)
      REAL(KIND=4)  AW,AS,BW,BS,P0W,P0S
      REAL(KIND=4)  PWVSAS,PSASAS
C
      INTEGER MAXN,I
      PARAMETER(MAXN=21)
      REAL(KIND=4)  T(MAXN),W(MAXN),LNPW(MAXN),LNPS(MAXN),
     +     Y2W(MAXN),Y2S(MAXN),YWORK(MAXN)
      LOGICAL FIRST
C
C     Frost point data from Gable et al.:
      DATA (W(I),T(I),I=1,MAXN)/
     +0.3755, 200.05,0.3831, 201.76,0.3858, 203.24,0.3888, 204.74,
     +0.4008, 209.88,0.3956, 207.86,0.4119, 214.56,0.4169, 216.43,
     +0.4264, 219.67,0.4349, 222.20,0.4441, 224.68,0.4620, 229.50,
     +0.4721, 231.66,0.4775, 233.29,0.4947, 236.52,0.5082, 238.95,
     +0.5308, 242.18,0.5416, 243.26,0.5537, 244.05,0.5636, 244.53,
     +0.5764, 244.79/
      DATA FIRST/.TRUE./
      SAVE FIRST,T,LNPW,LNPS,Y2W,Y2S,AW,BW,AS,BS
C
      IF(FIRST) THEN
          FIRST=.FALSE.
          DO I=1,MAXN
              LNPW(I)=PWVSAS(T(I),W(I))
              LNPS(I)=PSASAS(T(I),W(I))
          ENDDO
          CALL SPLINE(T,LNPW,MAXN,YWORK,Y2W)
          CALL SPLINE(T,LNPS,MAXN,YWORK,Y2S)
C
          AW=(LNPW(3)-LNPW(2))/(T(3)-T(2))
          BW=LNPW(2)-AW*T(2)
          AS=(LNPS(3)-LNPS(2))/(T(3)-T(2))
          BS=LNPS(2)-AS*T(2)
      ENDIF
C
      IF(TAIR.GE.T(1).AND.TAIR.LE.T(MAXN)) THEN
          CALL SPLINT(T,LNPW,Y2W,MAXN,TAIR,P0W)
          CALL SPLINT(T,LNPS,Y2S,MAXN,TAIR,P0S)
      ELSE IF(TAIR.LT.T(1)) THEN
          P0W=AW*TAIR+BW
          P0S=AS*TAIR+BS
      ENDIF
C
      M=A+B*TAIR
C
      PSASAT=P0S+M*(ALOG(PPWV)-P0W)
      RETURN
      END
C
*****************************************************************************
      REAL FUNCTION PNASNW(TAIR,WSA,WNA)
*     REAL FUNCTION PNASNW(TAIR,WSA,WNA)
*****************************************************************************
*
*     Natural logaritm of saturated nitric acid vapor pressure over plane
*     TERNARY sulfuric/nitric acid solution (SNW).
*
*     Source: Luo et al.: GRL 22,247,1995
*     
*
*     Input:  TAIR: Temperature (K) [185; 235 k]
*             WSA:  Weight fraction of H2SO4  ]0.0; 0.7[
*             WNA:  Weight fraction of HNO3   ]0;0.7[ 
*     Output: Natural logaritm of nitric acid vapor pressure 
*             over ternary sulfuric/nitric acid solution   ( ln(Pa) )
*
*
*
C
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WSA,WNA
CC      REAL(KIND=4)  WMIN,WMAX
      REAL(KIND=4)  MBPAS,A,B,SQW1,SQW2,W1,W2
      REAL(KIND=4)  A0,A1,A2,A3,A4,A5,A6,A7,B0,B1,B2,B3,B4,B5,B6,B7
C
CC      PARAMETER(WMIN=1.0E-10,WMAX=0.7E0-WMIN)
      DATA A0,A1,A2,A3,A4,A5,A6,A7/
     +22.74,29.0,-0.457,-5.03,1.447,-29.72,-13.9,6.1/ 
      DATA B0,B1,B2,B3,B4,B5,B6,B7/
     +-7689.8,-2896.1,2859.8,-274.2,-389.5,7281.0,6475.0,801.0/
      DATA MBPAS/4.60517019/
C
CC    W2=AMAX1(WMIN,AMIN1(WSA,WMAX))
CC    W1=AMAX1(0.7E0-W2,AMIN1(WNA,WMAX))
      W1=WNA
      W2=WSA
      SQW1=SQRT(W1)
      SQW2=SQRT(W2)
      A=A0+A1*W1+A2*W2+A3*SQW1+A4*SQW2+A5*W1*W1+A6*W1*W2+A7*W2*W2
      B=B0+B1*W1+B2*W2+B3*SQW1+B4*SQW2+B5*W1*W1+B6*W1*W2+B7*W2*W2
C
      PNASNW=A+(B/TAIR)+ALOG(W1*(W1+0.09E0*W2))+MBPAS
      RETURN
      END
C
C
*****************************************************************************
      REAL FUNCTION PWVSNW(TAIR,WSA,WNA)
*     REAL FUNCTION PWVSNW(TAIR,WSA,WNA)
*****************************************************************************
*
*     Natural logaritm of saturated water vapor pressure over plane
*     TERNARY sulfuric/nitric acid solution (SNW).
*
*     Source: Luo et al.: GRL 22,247,1995
*     
*     
*
*     Input:  TAIR: Temperature (K) [185; 235 k]
*             WSA:  Weight fraction of H2SO4  ]0.0; 0.7[
*             WNA:  Weight fraction of HNO3   ]0;0.7[ 
*     Output: Natural logaritm of water vapor pressure 
*             over ternary sulfuric/nitric acid solution   ( ln(Pa) )
*
*     
*     
*
*     
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,WSA,WNA
CC      REAL(KIND=4)  WMIN,WMAX
      REAL(KIND=4)  WH,MBPAS,W1,W2
CC      PARAMETER(WMIN=1.0E-10,WMAX=0.7E0-WMIN)
      DATA MBPAS/4.60517019/
C
CC    W2=AMAX1(WMIN,AMIN1(WSA,WMAX))
CC    W1=AMAX1(0.7E0-W2,AMIN1(WNA,WMAX))
      W1=WNA
      W2=WSA
      WH=W1+1.4408E0*W2
C
      PWVSNW=23.306E0-4.5261E0*W1-5.3465E0*W2+
     + WH*(7.451E0*W1+12.0E0*W2)-WH*WH*(4.0E0*W1+8.19E0*W2)+
     + (-5814.0E0+1033.0E0*W1+928.9E0*W2-
     +  WH*(2309.0E0*W1+1876.7E0*W2))/TAIR+
     + MBPAS
      RETURN
      END
C
*****************************************************************************
*     REAL FUNCTION PNASNW(TAIR,WSA,WNA)                                        
CC    REAL FUNCTION PNASNW(TAIR,WSA,WNA)
*****************************************************************************
*
*     Natural logaritm of saturated nitric acid vapor pressure over plane
*     TERNARY sulfuric/nitric acid solution (SNW).
*
*     Source: Zhang et al. J. Phys. Chem. 97, 8541, 1993
*             Beyer et al. GRL 21,871,1994
*     
*
*     Input:  TAIR: Temperature (K) 
*             WSA:  Weight fraction of H2SO4  ]0.1; 1.0[
*             WNA:  Weight fraction of HNO3   ]0;0.45[ 
*     Output: Natural logaritm of nitric acid vapor pressure 
*             over ternary sulfuric/nitric acid solution   ( ln(Pa) )
*
*     If WSA or WNA are outside the above ranges the function returns 0.
*
*
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,WSA,WNA
CC    DOUBLE PRECISION A(0:2),AA(0:2),BB(0:2),C1,C2,TINV,WSAS,WNAS
CC    INTEGER I
CC    DATA AA/12.88E0,-4.465E0,-0.569E0/
CC    DATA BB/-2733.0E0,-685.5E0,189.9E0/
CC    PARAMETER(C1=2.12490303E0,C2=2.30258509E0)
CC    SAVE AA,BB
C
C      IF(WSA.LT.0.35.OR.WSA.GT.0.75.OR.WNA.GT.0.15.OR.WNA.LE.0.0) THEN
C
CC    WSAS=WSA
CC    WNAS=AMAX1(WNA,1.0E-20)
CC    TINV=1.0E0/TAIR
CC    DO 100 I=0,2
CC        A(I)=AA(I)+BB(I)*TINV
CC100 CONTINUE
CC    PNASNW=SNGL((ALOG10(WNAS)+A(0)+A(1)*(1.0E0-WNAS-WSAS)-
CC   +             A(2)*(EXP((1-WSAS)**8)-1.0E0)+C1)*C2)
CC    RETURN
CC    END
*****************************************************************************
*     REAL FUNCTION PWVSNW(TAIR,WSA,WNA)                                        
CC    REAL FUNCTION PWVSNW(TAIR,WSA,WNA)
*****************************************************************************
*
*     Natural logaritm of saturated water vapor pressure over plane
*     TERNARY sulfuric/nitric acid solution (SNW).
*
*     Source: Zhang et al. J. Phys. Chem. 97, 8541, 1993
*             Beyer et al. GRL 21,871,1994
*     
*     
*
*     Input:  TAIR: Temperature (K)
*             WSA:  Weight fraction of H2SO4  ]0.1; 1.0[
*             WNA:  Weight fraction of HNO3   ]0;0.45[ 
*     Output: Natural logaritm of water vapor pressure 
*             over ternary sulfuric/nitric acid solution   ( ln(Pa) )
*
*     If WSA or WNA are outside the above ranges the function returns 0.
*     
*     
*
*     
*
CC    IMPLICIT NONE
CC    REAL(KIND=4)  TAIR,WSA,WNA
CC    DOUBLE PRECISION B(0:3),AA(0:3),BB(0:3),C1,C2,
CC   +                 TINV,XPLUSY,WSAS,WNAS
CC    INTEGER I
CC    DATA AA/12.09E0,-10.17E0,10.55E0,0.1187E0/
CC    DATA BB/-3202.0E0,3466.0E0,-4224.0E0,-59.74E0/
CC    PARAMETER(C1=2.12490303E0,C2=2.30258509E0)
CC    SAVE AA,BB
C
C      IF(WSA.LT.0.35.OR.WSA.GT.0.75.OR.WNA.GT.0.15.OR.WNA.LE.0.0) THEN
C          PWVSNW=0.0
C          RETURN
C      ENDIF
CC    IF(WSA+WNA.GT.0.999999) THEN
CC         PWVSNW=0.0
CC         RETURN
CC    ENDIF
C
CC    WSAS=WSA
CC    WNAS=WNA
CC    TINV=1.0E0/TAIR
CC    XPLUSY=WNAS+WSAS
CC    DO 100 I=0,3
CC        B(I)=AA(I)+BB(I)*TINV
CC100 CONTINUE
CC    PWVSNW=SNGL((ALOG10(1.0E0-XPLUSY)+B(0)+(B(2)*XPLUSY+B(1))*XPLUSY-
CC   +             B(3)*(EXP((1.0E0-WSAS)**8)-1.0E0)+C1)*C2)
CC    RETURN
CC    END
C
*******************************************************************************
      REAL FUNCTION PNANSS(TAIR,MRNATO,MRWVTO)
*     REAL FUNCTION PNANSS(TAIR,MRNATO,MRWVTO)
*******************************************************************************
*
*     Nitric acid vapor pressure over HNO3/H2O SOLID solution.
*
*     Source: Tabazadeh & Toon, JGR, 1995
*     
*     
*
*     Input:  TAIR: Temperature (K)
*             MRNATO: Total nitric acid volume mixing ratio 
*             MRWVTO: Total water vapor volume mixing ratio 
*     Output: Nitric acid vapor pressure over HNO3/H2O SOLID solution (Pa)
*
      IMPLICIT NONE
      REAL(KIND=4)  TAIR,MRNATO,MRWVTO
      REAL(KIND=4)  MRNA,MRWV,FW,A,B,A0,A1,A2,B0,B1,B2,F0,F1,
     +              T,TMIN,P,DPDT
      PARAMETER(TMIN=190.5+2.0)
      DATA A0,A1,A2/7988.2E-7,-84.386E-7,0.22286E-7/
      DATA B0,B1,B2/1556.7E-7,-15.772E-7,0.04E-7/
      DATA F0,F1/0.1027,-18.1976/
      SAVE A0,A1,A2,B0,B1,B2,F0,F1
C
      T=AMAX1(TMIN,TAIR)
      MRNA=MRNATO*1.0E+9
      MRWV=MRWVTO*1.0E+6
      FW=AMAX1(0.0E0,AMIN1(MRWV*(F0*T+F1)/MRNA,1.0E0))
      A=(A2*T+A1)*T+A0
      B=(B2*T+B1)*T+B0
      P=(A+B*FW)*133.32237E0
      IF(TAIR.GE.TMIN) THEN
          PNANSS=P
      ELSE
          P=ALOG(P)
          T=T+1.0E0
          FW=AMAX1(0.0E0,AMIN1(MRWV*(F0*T+F1)/MRNA,1.0E0))
          A=(A2*T+A1)*T+A0
          B=(B2*T+B1)*T+B0
          DPDT=ALOG((A+B*FW)*133.32237E0)-P
          PNANSS=EXP(P+DPDT*(TAIR-TMIN))
      ENDIF
      RETURN
      END
C
*******************************************************************************
*     REAL FUNCTION CPH2O(W)
      REAL FUNCTION CPH2O(W)
*******************************************************************************
*
*     Relative partial molal heat capacity of water (cal/(deg mole) in 
*     sulfuric acid solution, as a function of H2SO4 weight fraction [0;1],
*     calculated by cubic spline fitting.
*
*     Source: Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.
*
      IMPLICIT NONE
      INTEGER NPOINT,I
      PARAMETER(NPOINT=109)
      REAL(KIND=4)  W,WTAB(NPOINT),CPHTAB(NPOINT),
     +              Y2(NPOINT),YWORK(NPOINT),CPH
      LOGICAL FIRST
      DATA (WTAB(I),I=1,NPOINT)/
     +0.00000,0.08932,0.09819,0.10792,0.11980,0.13461,0.15360,0.16525,
     +0.17882,0.19482,0.21397,0.23728,0.26629,0.27999,0.29517,0.31209,
     +0.33107,0.35251,0.36430,0.37691,0.39043,0.40495,0.42059,0.43749,
     +0.44646,0.45580,0.46555,0.47572,0.48634,0.49745,0.50908,0.52126,
     +0.53405,0.54747,0.56159,0.57646,0.58263,0.58893,0.59537,0.60195,
     +0.60868,0.61557,0.62261,0.62981,0.63718,0.64472,0.65245,0.66037,
     +0.66847,0.67678,0.68530,0.69404,0.70300,0.71220,0.72164,0.73133,
     +0.73628,0.74129,0.74637,0.75152,0.75675,0.76204,0.76741,0.77286,
     +0.77839,0.78399,0.78968,0.79545,0.80130,0.80724,0.81327,0.81939,
     +0.82560,0.83191,0.83832,0.84482,0.85143,0.85814,0.86495,0.87188,
     +0.87892,0.88607,0.89334,0.90073,0.90824,0.91588,0.92365,0.93156,
     +0.93959,0.94777,0.95610,0.96457,0.97319,0.98196,0.99090,0.99270,
     +0.99452,0.99634,0.99725,0.99817,0.99835,0.99853,0.99872,0.99890,
     +0.99908,0.99927,0.99945,0.99963,0.99982/
      DATA (CPHTAB(I),I=1,NPOINT)/
     + 17.996, 17.896, 17.875, 17.858, 17.840, 17.820, 17.800, 17.791,
     + 17.783, 17.777, 17.771, 17.769, 17.806, 17.891, 18.057, 18.248,
     + 18.429, 18.567, 18.613, 18.640, 18.660, 18.660, 18.642, 18.592,
     + 18.544, 18.468, 18.348, 18.187, 17.995, 17.782, 17.562, 17.352,
     + 17.162, 16.993, 16.829, 16.657, 16.581, 16.497, 16.405, 16.302,
     + 16.186, 16.053, 15.901, 15.730, 15.540, 15.329, 15.101, 14.853,
     + 14.586, 14.296, 13.980, 13.638, 13.274, 12.896, 12.507, 12.111,
     + 11.911, 11.711, 11.514, 11.320, 11.130, 10.940, 10.760, 10.570,
     + 10.390, 10.200, 10.000, 9.8400, 9.7600, 9.7900, 9.9500, 10.310,
     + 10.950, 11.960, 13.370, 15.060, 16.860, 18.550, 20.000, 21.170,
     + 22.030, 22.570, 22.800, 22.750, 22.420, 21.850, 21.120, 20.280,
     + 19.360, 18.350, 17.220, 15.940, 14.490, 12.840, 10.800, 9.8000,
     + 7.8000, 3.8000,0.20000,-5.4000,-7.0000,-8.8000,-10.900,-13.500,
     +-17.000,-22.000,-29.000,-40.000,-59.000/
      DATA FIRST/.TRUE./
      SAVE FIRST,WTAB,CPHTAB,Y2
C
      IF(FIRST) THEN
          FIRST=.FALSE.
          CALL SPLINE(WTAB,CPHTAB,NPOINT,YWORK,Y2)
      ENDIF
      CALL SPLINT(WTAB,CPHTAB,Y2,NPOINT,W,CPH)
      CPH2O=CPH
      RETURN
      END
C
*******************************************************************************
      REAL FUNCTION FFH2O(W)
*     REAL FUNCTION FFH2O(W)
*******************************************************************************
*
*     Relative partial molal free energy water (cal/mole) in 
*     sulfuric acid solution, as a function of H2SO4 weight fraction [0;1],
*     calculated by cubic spline fitting.
*
*     Source: Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.
*
      IMPLICIT NONE
      INTEGER NPOINT,I
      PARAMETER(NPOINT=110)
      REAL(KIND=4)  W,WTAB(NPOINT),FFTAB(NPOINT),
     +              Y2(NPOINT),YWORK(NPOINT),FF
      LOGICAL FIRST
      DATA (WTAB(I),I=1,NPOINT)/
     +0.00000,0.08932,0.09819,0.10792,0.11980,0.13461,0.15360,0.16525,
     +0.17882,0.19482,0.21397,0.23728,0.26629,0.27999,0.29517,0.31209,
     +0.33107,0.35251,0.36430,0.37691,0.39043,0.40495,0.42059,0.43749,
     +0.44646,0.45580,0.46555,0.47572,0.48634,0.49745,0.50908,0.52126,
     +0.53405,0.54747,0.56159,0.57646,0.58263,0.58893,0.59537,0.60195,
     +0.60868,0.61557,0.62261,0.62981,0.63718,0.64472,0.65245,0.66037,
     +0.66847,0.67678,0.68530,0.69404,0.70300,0.71220,0.72164,0.73133,
     +0.73628,0.74129,0.74637,0.75152,0.75675,0.76204,0.76741,0.77286,
     +0.77839,0.78399,0.78968,0.79545,0.80130,0.80724,0.81327,0.81939,
     +0.82560,0.83191,0.83832,0.84482,0.85143,0.85814,0.86495,0.87188,
     +0.87892,0.88607,0.89334,0.90073,0.90824,0.91588,0.92365,0.93156,
     +0.93959,0.94777,0.95610,0.96457,0.97319,0.98196,0.99090,0.99270,
     +0.99452,0.99634,0.99725,0.99817,0.99835,0.99853,0.99872,0.99890,
     +0.99908,0.99927,0.99945,0.99963,0.99982, 1.0000/
      DATA (FFTAB(I),I=1,NPOINT)/
     +0.00000, 22.840, 25.810, 29.250, 33.790, 39.970, 48.690, 54.560,
     + 61.990, 71.790, 85.040, 103.70, 130.70, 145.20, 163.00, 184.50,
     + 211.50, 245.60, 266.40, 290.10, 317.40, 349.00, 385.60, 428.40,
     + 452.50, 478.80, 507.50, 538.80, 573.30, 611.60, 653.70, 700.50,
     + 752.60, 810.60, 875.60, 948.60, 980.60, 1014.3, 1049.7, 1087.1,
     + 1126.7, 1168.7, 1213.5, 1261.2, 1312.0, 1366.2, 1424.3, 1486.0,
     + 1551.8, 1622.3, 1697.8, 1778.5, 1864.9, 1956.8, 2055.8, 2162.0,
     + 2218.0, 2276.0, 2337.0, 2400.0, 2466.0, 2535.0, 2607.0, 2682.0,
     + 2760.0, 2842.0, 2928.0, 3018.0, 3111.0, 3209.0, 3311.0, 3417.0,
     + 3527.0, 3640.0, 3757.0, 3878.0, 4002.0, 4130.0, 4262.0, 4397.0,
     + 4535.0, 4678.0, 4824.0, 4973.0, 5128.0, 5287.0, 5454.0, 5630.0,
     + 5820.0, 6031.0, 6268.0, 6541.0, 6873.0, 7318.0, 8054.0, 8284.0,
     + 8579.0, 8997.0, 9295.0, 9720.0, 9831.0, 9954.0, 10092., 10248.,
     + 10423., 10618., 10838., 11099., 11460., 12014./
      DATA FIRST/.TRUE./
      SAVE FIRST,WTAB,FFTAB,Y2
C
      IF(FIRST) THEN
          FIRST=.FALSE.
          CALL SPLINE(WTAB,FFTAB,NPOINT,YWORK,Y2)
      ENDIF
      CALL SPLINT(WTAB,FFTAB,Y2,NPOINT,W,FF)
      FFH2O=FF
      RETURN
      END
C
*******************************************************************************
      REAL FUNCTION LH2O(W)
*     REAL FUNCTION LH2O(W)
*******************************************************************************
*
*     Relative partial molal heat content of water (cal/mole) in 
*     sulfuric acid solution, as a function of H2SO4 weight fraction [0;1],
*     calculated by cubic spline fitting.
*
*     Source: Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.
*
      IMPLICIT NONE
      INTEGER NPOINT,I
      PARAMETER(NPOINT=110)
      REAL(KIND=4)  W,WTAB(NPOINT),LTAB(NPOINT),
     +              Y2(NPOINT),YWORK(NPOINT),L
      LOGICAL FIRST
      DATA (WTAB(I),I=1,NPOINT)/
     +0.00000,0.08932,0.09819,0.10792,0.11980,0.13461,0.15360,0.16525,
     +0.17882,0.19482,0.21397,0.23728,0.26629,0.27999,0.29517,0.31209,
     +0.33107,0.35251,0.36430,0.37691,0.39043,0.40495,0.42059,0.43749,
     +0.44646,0.45580,0.46555,0.47572,0.48634,0.49745,0.50908,0.52126,
     +0.53405,0.54747,0.56159,0.57646,0.58263,0.58893,0.59537,0.60195,
     +0.60868,0.61557,0.62261,0.62981,0.63718,0.64472,0.65245,0.66037,
     +0.66847,0.67678,0.68530,0.69404,0.70300,0.71220,0.72164,0.73133,
     +0.73628,0.74129,0.74637,0.75152,0.75675,0.76204,0.76741,0.77286,
     +0.77839,0.78399,0.78968,0.79545,0.80130,0.80724,0.81327,0.81939,
     +0.82560,0.83191,0.83832,0.84482,0.85143,0.85814,0.86495,0.87188,
     +0.87892,0.88607,0.89334,0.90073,0.90824,0.91588,0.92365,0.93156,
     +0.93959,0.94777,0.95610,0.96457,0.97319,0.98196,0.99090,0.99270,
     +0.99452,0.99634,0.99725,0.99817,0.99835,0.99853,0.99872,0.99890,
     +0.99908,0.99927,0.99945,0.99963,0.99982, 1.0000/
      DATA (LTAB(I),I=1,NPOINT)/
     +0.00000, 5.2900, 6.1000, 7.1800, 8.7800, 11.210, 15.290, 18.680,
     + 23.700, 31.180, 42.500, 59.900, 89.200, 106.70, 128.60, 156.00,
     + 190.40, 233.80, 260.10, 290.00, 324.00, 362.50, 406.50, 456.10,
     + 483.20, 512.40, 543.60, 577.40, 613.80, 653.50, 696.70, 744.50,
     + 797.20, 855.80, 921.70, 995.70, 1028.1, 1062.3, 1098.3, 1136.4,
     + 1176.7, 1219.3, 1264.7, 1313.0, 1364.3, 1418.9, 1477.3, 1539.9,
     + 1607.2, 1679.7, 1757.9, 1842.7, 1934.8, 2035.4, 2145.5, 2267.0,
     + 2332.0, 2401.0, 2473.0, 2550.0, 2631.0, 2716.0, 2807.0, 2904.0,
     + 3007.0, 3118.0, 3238.0, 3367.0, 3507.0, 3657.0, 3821.0, 3997.0,
     + 4186.0, 4387.0, 4599.0, 4819.0, 5039.0, 5258.0, 5476.0, 5694.0,
     + 5906.0, 6103.0, 6275.0, 6434.0, 6592.0, 6743.0, 6880.0, 7008.0,
     + 7133.0, 7255.0, 7376.0, 7497.0, 7618.0, 7739.0, 7855.0, 7876.0,
     + 7905.0, 7985.0, 8110.0, 8415.0, 8515.0, 8655.0, 8835.0, 9125.0,
     + 9575.0, 10325., 11575., 13500., 15200., 16125./
      DATA FIRST/.TRUE./
      SAVE FIRST,WTAB,LTAB,Y2
C
      IF(FIRST) THEN
          FIRST=.FALSE.
          CALL SPLINE(WTAB,LTAB,NPOINT,YWORK,Y2)
      ENDIF
      CALL SPLINT(WTAB,LTAB,Y2,NPOINT,W,L)
      LH2O=L
      RETURN
      END
C
*******************************************************************************
      REAL FUNCTION LSAS(W)
*     REAL FUNCTION LSAS(W)
*******************************************************************************
*
*     Relative partial molal heat content of H2SO4 (cal/mole) in 
*     sulfuric acid solution, as a function of H2SO4 weight fraction [0;1],
*     calculated by cubic spline fitting.
*
*     Source: Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.
*
      IMPLICIT NONE
      INTEGER NPOINT,I
      PARAMETER(NPOINT=109)
      REAL(KIND=4)  W,WTAB(NPOINT),LTAB(NPOINT),
     +              Y2(NPOINT),YWORK(NPOINT),L
      LOGICAL FIRST
      DATA (WTAB(I),I=1,NPOINT)/
     +0.08932,0.09819,0.10792,0.11980,0.13461,0.15360,0.16525,
     +0.17882,0.19482,0.21397,0.23728,0.26629,0.27999,0.29517,0.31209,
     +0.33107,0.35251,0.36430,0.37691,0.39043,0.40495,0.42059,0.43749,
     +0.44646,0.45580,0.46555,0.47572,0.48634,0.49745,0.50908,0.52126,
     +0.53405,0.54747,0.56159,0.57646,0.58263,0.58893,0.59537,0.60195,
     +0.60868,0.61557,0.62261,0.62981,0.63718,0.64472,0.65245,0.66037,
     +0.66847,0.67678,0.68530,0.69404,0.70300,0.71220,0.72164,0.73133,
     +0.73628,0.74129,0.74637,0.75152,0.75675,0.76204,0.76741,0.77286,
     +0.77839,0.78399,0.78968,0.79545,0.80130,0.80724,0.81327,0.81939,
     +0.82560,0.83191,0.83832,0.84482,0.85143,0.85814,0.86495,0.87188,
     +0.87892,0.88607,0.89334,0.90073,0.90824,0.91588,0.92365,0.93156,
     +0.93959,0.94777,0.95610,0.96457,0.97319,0.98196,0.99090,0.99270,
     +0.99452,0.99634,0.99725,0.99817,0.99835,0.99853,0.99872,0.99890,
     +0.99908,0.99927,0.99945,0.99963,0.99982, 1.0000/
      DATA (LTAB(I),I=1,NPOINT)/
     + 17138.0,17088.0,17037.0,16970.0,16879.0,16748.0,16650.0,
     + 16519.0,16342.0,16102.0,15778.0,15305.0,15052.0,14757.0,14415.0,
     + 14020.0,13565.0,13308.0,13032.0,12735.0,12417.0,12077.0,11617.0,
     + 11531.0,11338.0,11139.0,10932.0,10718.0,10495.0,10262.0,10018.0,
     + 9761.0, 9490.0, 9202.0, 8897.0, 8769.0, 8637.0, 8502.0, 8363.0,
     + 8220.0, 8073.0, 7921.0, 7764.0, 7602.0, 7436.0, 7264.0, 7085.0,
     + 6900.0, 6708.0, 6509.0, 6301.0, 6085.0, 5859.0, 5622.0, 5373.0,
     + 5244.0, 5111.0, 4975.0, 4836.0, 4692.0, 4545.0, 4393.0, 4236.0,
     + 4073.0, 3903.0, 3726.0, 3542.0, 3350.0, 3151.0, 2943.0, 2727.0,
     + 2505.0, 2279.0, 2051.0, 1826.0, 1611.0, 1409.0, 1218.0, 1038.0,
     +  874.0,  731.0,  614.0,  515.0,  424.0,  345.0,  280.0,  225.4,
     +  178.5,  138.9,  105.6,   78.4,   57.2,   42.1,   33.3,   32.2,
     +   31.0,   29.2,   27.1,   23.4,   22.4,   21.2,   19.9,   18.0,
     +   15.5,   12.2,   7.85,   3.10,   0.52, 0.0000/
      DATA FIRST/.TRUE./
      SAVE FIRST,WTAB,LTAB,Y2
C
      IF(FIRST) THEN
          FIRST=.FALSE.
          CALL SPLINE(WTAB,LTAB,NPOINT,YWORK,Y2)
      ENDIF
      CALL SPLINT(WTAB,LTAB,Y2,NPOINT,W,L)
      LSAS=L
      RETURN
      END
C
*******************************************************************************
      REAL FUNCTION ALH2O(W)
*     REAL FUNCTION ALH2O(W)
*******************************************************************************
*
*     Relative partial molal temperature derivative of heat capacity (water) 
*     in sulfuric acid solution, (cal/deg**2), calculated by 
*     cubic spline fitting.
*
*     Source: Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.
*
      IMPLICIT NONE
      INTEGER NPOINT,I
      PARAMETER(NPOINT=96)
      REAL(KIND=4)  W,WTAB(NPOINT),ATAB(NPOINT),
     +              Y2(NPOINT),YWORK(NPOINT),A
      LOGICAL FIRST
      DATA (WTAB(I),I=1,NPOINT)/
     +0.29517,0.31209,
     +0.33107,0.35251,0.36430,0.37691,0.39043,0.40495,0.42059,0.43749,
     +0.44646,0.45580,0.46555,0.47572,0.48634,0.49745,0.50908,0.52126,
     +0.53405,0.54747,0.56159,0.57646,0.58263,0.58893,0.59537,0.60195,
     +0.60868,0.61557,0.62261,0.62981,0.63718,0.64472,0.65245,0.66037,
     +0.66847,0.67678,0.68530,0.69404,0.70300,0.71220,0.72164,0.73133,
     +0.73628,0.74129,0.74637,0.75152,0.75675,0.76204,0.76741,0.77286,
     +0.77839,0.78399,0.78968,0.79545,0.80130,0.80724,0.81327,0.81939,
     +0.82560,0.83191,0.83832,0.84482,0.85143,0.85814,0.86495,0.87188,
     +0.87892,0.88607,0.89334,0.90073,0.90824,0.91588,0.92365,0.93156,
     +0.93959,0.94777,0.95610,0.96457,0.97319,0.98196,0.99090,0.99270,
     +0.99452,0.99634,0.99725,0.99817,0.99835,0.99853,0.99872,0.99890,
     +0.99908,0.99927,0.99945,0.99963,0.99982, 1.0000/
      DATA (ATAB(I),I=1,NPOINT)/
     + 0.0190, 0.0182, 0.0180, 0.0177, 0.0174, 0.0169, 0.0167, 0.0164,
     + 0.0172, 0.0212, 0.0239, 0.0264, 0.0276, 0.0273, 0.0259, 0.0238,
     + 0.0213, 0.0190, 0.0170, 0.0155, 0.0143, 0.0133, 0.0129, 0.0124,
     + 0.0120, 0.0114, 0.0106, 0.0097, 0.0084, 0.0067, 0.0047, 0.0024,
     +-0.0002,-0.0031,-0.0063,-0.0097,-0.0136,-0.0178,-0.0221,-0.0263,
     +-0.0303,-0.0340,-0.0352,-0.0360,-0.0362,-0.0356,-0.0343,-0.0321,
     +-0.0290,-0.0251,-0.0201,-0.0137,-0.0058, 0.0033, 0.0136, 0.0254,
     + 0.0388, 0.0550, 0.0738, 0.0962, 0.1198, 0.1300, 0.1208, 0.0790,
     + 0.0348, 0.0058,-0.0102,-0.0211,-0.0292,-0.0350,-0.0390,-0.0418,
     +-0.0432,-0.0436,-0.0429,-0.0411,-0.0384,-0.0346,-0.0292,-0.0220,
     +-0.0130,-0.0110,-0.0080,-0.0060,-0.0040,-0.0030,-0.0030,-0.0020,
     +-0.0020,-0.0020,-0.0020,-0.0010,-0.0010, 0.0000, 0.0000, 0.0000/
      DATA FIRST/.TRUE./
      SAVE FIRST,WTAB,ATAB,Y2
C
      IF(FIRST) THEN
          FIRST=.FALSE.
          CALL SPLINE(WTAB,ATAB,NPOINT,YWORK,Y2)
      ENDIF
      CALL SPLINT(WTAB,ATAB,Y2,NPOINT,AMAX1(WTAB(1),W),A)
      ALH2O=A
      RETURN
      END
C
*******************************************************************************
      REAL FUNCTION FFSAS(W)
*     REAL FUNCTION FFSAS(W)
*******************************************************************************
*
*     Relative partial molal free energy of H2SO4 
*     in sulfuric acid solution, (cal/mole), calculated by 
*     cubic spline fitting.
*
*     Source: Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.
*
      IMPLICIT NONE
      INTEGER NPOINT,I
      PARAMETER(NPOINT=109)
      REAL(KIND=4)  W,WTAB(NPOINT),FTAB(NPOINT),
     +              Y2(NPOINT),YWORK(NPOINT),F
      LOGICAL FIRST
      DATA (WTAB(I),I=1,NPOINT)/
     +0.08932,0.09819,0.10792,0.11980,0.13461,0.15360,0.16525,
     +0.17882,0.19482,0.21397,0.23728,0.26629,0.27999,0.29517,0.31209,
     +0.33107,0.35251,0.36430,0.37691,0.39043,0.40495,0.42059,0.43749,
     +0.44646,0.45580,0.46555,0.47572,0.48634,0.49745,0.50908,0.52126,
     +0.53405,0.54747,0.56159,0.57646,0.58263,0.58893,0.59537,0.60195,
     +0.60868,0.61557,0.62261,0.62981,0.63718,0.64472,0.65245,0.66037,
     +0.66847,0.67678,0.68530,0.69404,0.70300,0.71220,0.72164,0.73133,
     +0.73628,0.74129,0.74637,0.75152,0.75675,0.76204,0.76741,0.77286,
     +0.77839,0.78399,0.78968,0.79545,0.80130,0.80724,0.81327,0.81939,
     +0.82560,0.83191,0.83832,0.84482,0.85143,0.85814,0.86495,0.87188,
     +0.87892,0.88607,0.89334,0.90073,0.90824,0.91588,0.92365,0.93156,
     +0.93959,0.94777,0.95610,0.96457,0.97319,0.98196,0.99090,0.99270,
     +0.99452,0.99634,0.99725,0.99817,0.99835,0.99853,0.99872,0.99890,
     +0.99908,0.99927,0.99945,0.99963,0.99982, 1.0000/
      DATA (FTAB(I),I=1,NPOINT)/
     +15811.00,15654.00,15491.00,15299.00,15068.00,14786.00,14617.00,
     +14423.00,14191.00,13910.00,13564.00,13126.00,12915.00,12676.00,
     +12407.00,12097.00,11740.00,11537.00,11318.00,11079.00,10819.00,
     +10535.00,10225.00,10060.00, 9885.00, 9703.00, 9511.00, 9308.00,
     + 9093.00, 8867.00, 8627.00, 8373.00, 8105.00, 7821.00, 7520.00,
     + 7393.00, 7264.00, 7131.00, 6994.00, 6854.00, 6709.00, 6559.00,
     + 6404.00, 6244.00, 6079.00, 5908.00, 5732.00, 5551.00, 5364.00,
     + 5171.00, 4974.00, 4771.00, 4564.00, 4351.00, 4134.00, 4022.00,
     + 3910.00, 3796.00, 3681.00, 3564.00, 3446.00, 3325.00, 3203.00,
     + 3080.00, 2955.00, 2829.00, 2701.00, 2572.00, 2442.00, 2312.00,
     + 2183.00, 2054.00, 1926.00, 1800.00, 1677.00, 1555.00, 1437.00,
     + 1322.00, 1211.00, 1103.00, 1000.00,  902.00,  808.00,  719.00,
     +  635.50,  556.60,  481.50,  409.80,  341.60,  276.90,  215.70,
     +  157.80,  102.40,   49.49,   39.17,   28.89,   18.59,   14.33,
     +    8.18,    7.13,    6.08,    5.04,    4.03,    3.07,    2.19,
     +    1.44,    0.78,    0.25,    0.00/
      DATA FIRST/.TRUE./
      SAVE FIRST,WTAB,FTAB,Y2
C
      IF(FIRST) THEN
          FIRST=.FALSE.
          CALL SPLINE(WTAB,FTAB,NPOINT,YWORK,Y2)
      ENDIF
      CALL SPLINT(WTAB,FTAB,Y2,NPOINT,W,F)
      FFSAS=F
      RETURN
      END
*****************************************************************************
      REAL FUNCTION TCNAT(PPNA,PPWV)
*     REAL FUNCTION TCNAT(PPNA,PPWV)
*****************************************************************************
*
*     Function used to calculate the condensation temperature of
*     nitric acid trihydrate.
*
*     Input:
*     PPNA:   (Pa)        Partial pressure of HNO3
*     PPWV:   (Pa)        Partial pressure of H2O
*
*     Output:
*     TCNAT:  (K)         Condensation temperature of NAT
*
      IMPLICIT NONE
      REAL(KIND=4)  PPWV,PPNA
      REAL(KIND=4)  PNANAT
      INTEGER JMAX,J
      PARAMETER (JMAX=40)
      REAL(KIND=4)  TMAX,TMIN,T1,T2
      PARAMETER(TMAX=273.0E0,TMIN=160.0E0)
      REAL(KIND=4)  FMID,F,DX,XMID,TACC
      PARAMETER(TACC=0.01E0)
C
      REAL(KIND=4)  FUNC,T,PPN,PPW
      FUNC(T,PPN,PPW)=PPN-PNANAT(T,PPW)
C
      T1=TMIN
      T2=TMAX
C
      FMID=FUNC(T2,PPNA,PPWV)
      F=FUNC(T1,PPNA,PPWV)
CC    IF(F*FMID.GT.0.)
CC   +  PAUSE 'TCNAT - Root must be bracketed for bisection.'
      IF(F.LT.0.)THEN
        TCNAT=T1
        DX=T2-T1
      ELSE
        TCNAT=T2
        DX=T1-T2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5E0
        XMID=TCNAT+DX
        FMID=FUNC(XMID,PPNA,PPWV)
        IF(FMID.LE.0.) TCNAT=XMID
        IF(ABS(DX).LT.TACC.OR.FMID.EQ.0.) RETURN
11    CONTINUE
CC    PAUSE 'TCNAT -too many bisections'
      END
C
*****************************************************************************
      REAL FUNCTION TCNAD(PPNA,PPWV)
*     REAL FUNCTION TCNAD(PPNA,PPWV)
*****************************************************************************
*
*     Function used to calculate the condensation temperature of
*     nitric acid dihydrate.
*
*     Input:
*     PPNA:   (Pa)        Partial pressure of HNO3
*     PPWV:   (Pa)        Partial pressure of H2O
*
*     Output:
*     TCNAD:  (K)         Condensation temperature of NAD
*
      IMPLICIT NONE
      REAL(KIND=4)  PPWV,PPNA
      REAL(KIND=4)  PNANAD
      INTEGER JMAX,J
      PARAMETER (JMAX=40)
      REAL(KIND=4)  TMAX,TMIN,T1,T2
      PARAMETER(TMAX=273.0E0,TMIN=160.0E0)
      REAL(KIND=4)  FMID,F,DX,XMID,TACC
      PARAMETER(TACC=0.01E0)
C
      REAL(KIND=4)  FUNC,T,PPN,PPW
      FUNC(T,PPN,PPW)=PPN-PNANAD(T,PPW)
C
      T1=TMIN
      T2=TMAX
C
      FMID=FUNC(T2,PPNA,PPWV)
      F=FUNC(T1,PPNA,PPWV)
CC    IF(F*FMID.GT.0.)
CC   +  PAUSE 'TCNAD - Root must be bracketed for bisection.'
      IF(F.LT.0.)THEN
        TCNAD=T1
        DX=T2-T1
      ELSE
        TCNAD=T2
        DX=T1-T2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5E0
        XMID=TCNAD+DX
        FMID=FUNC(XMID,PPNA,PPWV)
        IF(FMID.LE.0.) TCNAD=XMID
        IF(ABS(DX).LT.TACC.OR.FMID.EQ.0.) RETURN
11    CONTINUE
CC    PAUSE 'TCNAD -too many bisections'
      END
C
*****************************************************************************
      REAL FUNCTION TCICE(PPWV)
*     REAL FUNCTION TCICE(PPWV)
*****************************************************************************
*
*     Function used to calculate the condensation temperature of ice.
*
*     Input:
*     PPWV:   (Pa)        Partial pressure of H2O
*
*     Output:
*     TCICE:  (K)         Condensation temperature of ice
*
      IMPLICIT NONE
      REAL(KIND=4)  PPWV
      REAL(KIND=4)  PWVICE
      INTEGER JMAX,J
      PARAMETER (JMAX=40)
      REAL(KIND=4)  TMAX,TMIN,T1,T2
      PARAMETER(TMAX=273.0E0,TMIN=150.0E0)
      REAL(KIND=4)  FMID,F,DX,XMID,TACC
      PARAMETER(TACC=0.01E0)
C
      REAL(KIND=4)  FUNC,T,PPW
      FUNC(T,PPW)=PPW-PWVICE(T)
C
      T1=TMIN
      T2=TMAX
C
      FMID=FUNC(T2,PPWV)
      F=FUNC(T1,PPWV)
CC    IF(F*FMID.GT.0.) 
CC   +  PAUSE 'TCICE - Root must be bracketed for bisection.'
      IF(F.LT.0.)THEN
        TCICE=T1
        DX=T2-T1
      ELSE
        TCICE=T2
        DX=T1-T2
      ENDIF
      DO J=1,JMAX
        DX=DX*.5E0
        XMID=TCICE+DX
        FMID=FUNC(XMID,PPWV)
        IF(FMID.LE.0.) TCICE=XMID
        IF(ABS(DX).LT.TACC.OR.FMID.EQ.0.) RETURN
      ENDDO
CC    PAUSE 'TCICE - too many bisections'
      END
C
C******************************************************************************
      SUBROUTINE SPLINE(X,Y,N,WORK,Y2)
C******************************************************************************
C     Routine to calculate 2.nd derivatives of tabulated function
C     Y(i)=Y(Xi), to be used for cubic spline calculation.
C
      IMPLICIT NONE
      INTEGER N,I
      REAL(KIND=4)  X(N),Y(N),WORK(N),Y2(N)
      REAL(KIND=4)  SIG,P,QN,UN,YP1,YPN
C
      YP1=(Y(2)-Y(1))/(X(2)-X(1))
      YPN=(Y(N)-Y(N-1))/(X(N)-X(N-1))
      IF(YP1.GT.99.0E+30) THEN
          Y2(1)=0.0
          WORK(1)=0.0
      ELSE
          Y2(1)=-0.5E0
          WORK(1)=(3.0E0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO I=2,N-1
          SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
          P=SIG*Y2(I-1)+2.0E0
          Y2(I)=(SIG-1.0E0)/P
          WORK(I)=(6.0E0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     +             /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*WORK(I-1))/P
      ENDDO
      IF(YPN.GT.99.0E+30) THEN
          QN=0.0
          UN=0.0
      ELSE
          QN=0.5E0
          UN=(3.0E0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*WORK(N-1))/(QN*Y2(N-1)+1.0E0)
      DO I=N-1,1,-1
          Y2(I)=Y2(I)*Y2(I+1)+WORK(I)
      ENDDO
C
      RETURN
      END
C
C******************************************************************************
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
C******************************************************************************
C     Cubic spline calculation
C
      IMPLICIT NONE
      INTEGER KLO,KHI,N,K
      REAL(KIND=4)  XA(N),YA(N),Y2A(N)
      REAL(KIND=4)  X,Y,H,A,B
C
      KLO=1
      KHI=N
 1    IF(KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XA(K).GT.X) THEN
              KHI=K
          ELSE
              KLO=K
          ENDIF
          GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     +        ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.0E0
C
      RETURN
      END
C******************************************************************************
      SUBROUTINE SPLIN2(XA,YA,Y2A,N,K,X,Y)
C******************************************************************************
C     Cubic spline calculation
C     Restriction
C     XA(K).LE. X .LT.XA(K+1)
C
      IMPLICIT NONE
      INTEGER KLO,KHI,N,K
      REAL(KIND=4)  XA(N),YA(N),Y2A(N)
      REAL(KIND=4)  X,Y,H,A,B
C
      KLO=MAX0(1,MIN0(K,N-1))
      KHI=KLO+1
C
      H=XA(KHI)-XA(KLO)
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     +        ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.0E0
C
      RETURN
      END
C******************************************************************************
CC    SUBROUTINE SPLVEC(N,XA,YA,Y2A,GROWTH,X,Y)
C******************************************************************************
C     Cubic spline calculation on array:
C     Restriction
C     IF GROWTH
C     XA(I).LE. X(I) .LT.XA(I+1)
C
C     IF .NOT.GROWTH:
C     XA(I-1).LE. X(I) .LT.XA(I)
C
C
CC    IMPLICIT NONE
CC    INTEGER N,I,J
CC    LOGICAL GROWTH
CC    REAL(KIND=4)  XA(N),YA(N),Y2A(N)
CC    REAL(KIND=4)  X(N),Y(N)
CC    DOUBLE PRECISION H,A,B
C
CC    IF(GROWTH) THEN
CC        DO 100 I=1,N-1
CC            J=I+1
CC            H=XA(J)-XA(I)
CC            A=(XA(J)-X(J))/H
CC            B=(X(J)-XA(I))/H
CC            Y(J)=A*YA(I)+B*YA(J)
CC   +                +((A**3-A)*Y2A(I)+(B**3-B)*Y2A(J))*(H**2)/6.0E0
CC100     CONTINUE
CC        Y(1)=YA(1)
CC          Y(1)=YA(1)-(Y(2)-YA(2))
CC    ELSE
CC        DO 200 I=1,N-1
CC            J=I+1
CC            H=XA(J)-XA(I)
CC            A=(XA(J)-X(I))/H
CC            B=(X(I)-XA(I))/H
CC            Y(I)=A*YA(I)+B*YA(J)
CC   +                +((A**3-A)*Y2A(I)+(B**3-B)*Y2A(J))*(H**2)/6.0E0
CC200     CONTINUE
CC          Y(N)=Y(N)-(Y(N-1)-YA(N-1))
CC        Y(N)=YA(N)
CC    ENDIF
C
CC    RETURN
CC    END
C******************************************************************************
      SUBROUTINE SPLVEC(N,XA,YA,Y2A,GROW,X,Y)
C******************************************************************************
C     Cubic spline calculation on array:
C     Restriction
C     IF GROWTH
C     XA(I).LE. X(I) .LT.XA(I+1)
C
C     IF .NOT.GROWTH:
C     XA(I-1).LE. X(I) .LT.XA(I)
C
C
      IMPLICIT NONE
      INTEGER N,I,J,K,L
      REAL(KIND=4)  GROW
      REAL(KIND=4)  XA(N),YA(N),Y2A(N)
      REAL(KIND=4)  X(N),Y(N)
      DOUBLE PRECISION H,A,B
C
      IF(GROW.GE.0.0) THEN
          L=INT(GROW)
          DO I=1,N-L-1
              J=I+1
              K=J+L
              H=XA(J)-XA(I)
              A=(XA(J)-X(K))/H
              B=(X(K)-XA(I))/H
              Y(K)=A*YA(I)+B*YA(J)
     +                +((A**3-A)*Y2A(I)+(B**3-B)*Y2A(J))*(H**2)/6.0E0
          ENDDO
          DO I=1,L+1
              Y(I)=YA(I)
          ENDDO
      ELSE
          L=IABS(INT(GROW))
          DO I=1+L,N-1
              J=I+1
              K=I-L
              H=XA(J)-XA(I)
              A=(XA(J)-X(K))/H
              B=(X(K)-XA(I))/H
              Y(K)=A*YA(I)+B*YA(J)
     +                +((A**3-A)*Y2A(I)+(B**3-B)*Y2A(J))*(H**2)/6.0E0
          ENDDO
          DO I=N-L,N
              Y(I)=YA(I)
          ENDDO
      ENDIF
C
      RETURN
      END
C******************************************************************************
      SUBROUTINE SPLINE_EQ(N,X,YA,Y,WORK)
C******************************************************************************
C     Cubic spline calculation on equidistant shifted array for PSC box:
C     Input:
C     N:        Array dimension
C     X:        Constant bin shift
C     YA:       Function values
C     Output:
C     Y:        Function values, interpolated at regular bins by cubic spline
C     WORK:     Work array
C

C
      IMPLICIT NONE
      INTEGER N,I,J,K
      REAL(KIND=4)  X,YA(N)
      REAL(KIND=4)  WORK(N,2),Y(N)
      REAL(KIND=4)  SIG,SIGM,P,QN,UN,YP1,YPN
      DOUBLE PRECISION A,B
      PARAMETER(SIG=0.5E0,SIGM=SIG-1.0E0)
      YP1=YA(2)-YA(1)
      YPN=YA(N)-YA(N-1)
      IF(YP1.GT.99.0E+30) THEN
          WORK(1,2)=0.0
          WORK(1,1)=0.0
      ELSE
          WORK(1,2)=-0.5E0
          WORK(1,1)=3.0E0*(YA(2)-YA(1)-YP1)
      ENDIF
      DO I=2,N-1
          J=I+1
          K=I-1
          P=SIG*WORK(K,2)+2.0E0
          WORK(I,2)=SIGM/P
          WORK(I,1)=(3.0E0*(YA(J)-YA(I)-YA(I)+YA(K))-SIG*WORK(K,1))/P
      ENDDO
      IF(YPN.GT.99.0E+30) THEN
          QN=0.0
          UN=0.0
      ELSE
          QN=0.5E0
          UN=3.0E0*(YPN-YA(N)+YA(N-1))
      ENDIF
      WORK(N,2)=(UN-QN*WORK(N-1,1))/(QN*WORK(N-1,2)+1.0E0)
      WORK(N,2)=AMAX1(-5.0E0,AMIN1(WORK(N,2),5.0E0))
      DO I=N-1,1,-1
          WORK(I,2)=WORK(I,2)*WORK(I+1,2)+WORK(I,1)
          WORK(I,2)=AMAX1(-5.0E0,AMIN1(WORK(I,2),5.0E0))
      ENDDO

CCC
      IF(X.GT.0.0) THEN
          DO I=1,N-1
              J=I+1
              A=X
              B=1.0E0-A
              Y(J)=A*YA(I)+B*YA(J)
     +         +((A*A-1.0E0)*A*WORK(I,2)+(B*B-1.0E0)*B*WORK(J,2))/6.0E0
          ENDDO
          Y(1)=YA(1)
      ELSE
          DO I=1,N-1
              J=I+1
              B=-X
              A=1.0E0-B
              Y(I)=A*YA(I)+B*YA(J)
     +         +((A*A-1.0E0)*A*WORK(I,2)+(B*B-1.0E0)*B*WORK(J,2))/6.0E0
          ENDDO
          Y(N)=YA(N)
      ENDIF
C
      RETURN
      END
C******************************************************************************
      SUBROUTINE SPLINE_NON_EQ(N,XA,YA,Y,WORK)
C******************************************************************************
C     Cubic spline calculation on non-equidistant shifted array for PSC box:
C     Input:
C     N:        Array dimension
C     XA:       Bin shift
C     YA:       Function values
C     Output:
C     Y:        Function values, interpolated at regular bins by cubic spline
C     WORK:     Work array
C
      IMPLICIT NONE
      INTEGER N,I,J,K
      REAL(KIND=4)  XA(N),YA(N)
      REAL(KIND=4)  WORK(N,2),Y(N)
      REAL(KIND=4)  SIG,P,QN,UN,YP1,YPN
      DOUBLE PRECISION A,B,H


      YP1=(YA(2)-YA(1))/(1.0E0+XA(2)-XA(1))
      YPN=(YA(N)-YA(N-1))/(1.0E0+XA(N)-XA(N-1))
      IF(YP1.GT.99.0E+30) THEN
          WORK(1,2)=0.0
          WORK(1,1)=0.0
      ELSE
          WORK(1,2)=-0.5E0
          WORK(1,1)=(3.0E0/(1.0E0+XA(2)-XA(1)))*
     +            ((YA(2)-YA(1))/(1.0E0+XA(2)-XA(1))-YP1)
      ENDIF
      DO I=2,N-1
          J=I+1
          K=I-1
          SIG=(1.0E0+XA(I)-XA(K))/(2.0E0+XA(J)-XA(K))
          P=SIG*WORK(K,2)+2.0E0
          WORK(I,2)=(SIG-1.0E0)/P
          WORK(I,1)=(6.0E0*((YA(J)-YA(I))/
     +            (1.0E0+XA(J)-XA(I))-(YA(I)-YA(K))
     +             /(1.0E0+XA(I)-XA(K)))/(2.0E0+XA(J)-XA(K))-
     +             SIG*WORK(K,1))/P
      ENDDO
      IF(YPN.GT.99.0E+30) THEN
          QN=0.0
          UN=0.0
      ELSE
          QN=0.5E0
          UN=(3.0E0/(1.0E0+XA(N)-XA(N-1)))*
     +       (YPN-(YA(N)-YA(N-1))/(1.0E0+XA(N)-XA(N-1)))
      ENDIF
      WORK(N,2)=(UN-QN*WORK(N-1,1))/(QN*WORK(N-1,2)+1.0E0)
      WORK(N,2)=AMAX1(-5.0E0,AMIN1(WORK(N,2),5.0E0))
      DO I=N-1,1,-1
          WORK(I,2)=WORK(I,2)*WORK(I+1,2)+WORK(I,1)
          WORK(I,2)=AMAX1(-5.0E0,AMIN1(WORK(I,2),5.0E0))
      ENDDO
C
      DO I=1,N
          IF(XA(I).GT.0.0) THEN
              IF(I.EQ.1) THEN
                  Y(I)=YA(I)
              ELSE
                  K=I-1
                  A=XA(I)
                  B=1.0E0-XA(K)
                  H=A+B
                  A=A/H
                  B=B/H
                  Y(I)=A*YA(K)+B*YA(I)
     +                 +((A*A-1.0E0)*A*WORK(K,2)+
     +                 (B*B-1.0E0)*B*WORK(I,2))*(H*H)/6.0E0
              ENDIF
          ELSE IF(XA(I).LT.0.0) THEN
              IF(I.EQ.N) THEN
                  Y(I)=YA(I)
              ELSE
                  J=I+1
                  B=-XA(I)
                  A=1.0E0+XA(J)
                  H=A+B
                  A=A/H
                  B=B/H
                  Y(I)=A*YA(I)+B*YA(J)
     +                 +((A*A-1.0E0)*A*WORK(I,2)+
     +                 (B*B-1.0E0)*B*WORK(J,2))*(H*H)/6.0E0
              ENDIF
          ELSE IF(XA(I).EQ.0.0) THEN
              Y(I)=YA(I)
          ENDIF
      END DO
C
      RETURN
      END
C******************************************************************************
      SUBROUTINE SPLINE_ARRAY(N,XA,YA,X,Y,WORK)
C******************************************************************************
C     Cubic spline calculation on non-equidistant shifted array for PSC box:
C     Input:
C     N:        Array dimension
C     XA:       Shifted x-values
C     YA:       Function values at XA
C     X:        X-values where interpolated y-values are to be calculated
C     Output:
C     Y:        Function values, interpolated by cubic spline at X-values
C     WORK:     Work array; NOTICE DIMENTION !
C
      IMPLICIT NONE
      INTEGER N,I,J,K,KLO,KHI
      REAL(KIND=4)  XA(N),YA(N)
      REAL(KIND=4)  WORK(N,2),X(N),Y(N)
      REAL(KIND=4)  SIG,P,QN,UN,YP1,YPN
      DOUBLE PRECISION A,B,H


      YP1=(YA(2)-YA(1))/(XA(2)-XA(1))
      YPN=(YA(N)-YA(N-1))/(XA(N)-XA(N-1))
      IF(YP1.GT.99.0E+30) THEN
          WORK(1,2)=0.0
          WORK(1,1)=0.0
      ELSE
          WORK(1,2)=-0.5E0
          WORK(1,1)=(3.0E0/(XA(2)-XA(1)))*
     +            ((YA(2)-YA(1))/(XA(2)-XA(1))-YP1)
      ENDIF
      DO I=2,N-1
          J=I+1
          K=I-1
          SIG=(XA(I)-XA(K))/(XA(J)-XA(K))
          P=SIG*WORK(K,2)+2.0E0
          WORK(I,2)=(SIG-1.0E0)/P
          WORK(I,1)=(6.0E0*((YA(J)-YA(I))/
     +            (XA(J)-XA(I))-(YA(I)-YA(K))
     +             /(XA(I)-XA(K)))/(XA(J)-XA(K))-
     +             SIG*WORK(K,1))/P
      ENDDO
      IF(YPN.GT.99.0E+30) THEN
          QN=0.0
          UN=0.0
      ELSE
          QN=0.5E0
          UN=(3.0E0/(XA(N)-XA(N-1)))*
     +       (YPN-(YA(N)-YA(N-1))/(XA(N)-XA(N-1)))
      ENDIF
      WORK(N,2)=(UN-QN*WORK(N-1,1))/(QN*WORK(N-1,2)+1.0E0)
      WORK(N,2)=AMAX1(-5.0E0,AMIN1(WORK(N,2),5.0E0))
      DO I=N-1,1,-1
          WORK(I,2)=WORK(I,2)*WORK(I+1,2)+WORK(I,1)
          WORK(I,2)=AMAX1(-5.0E0,AMIN1(WORK(I,2),5.0E0))
      ENDDO
C
      DO I=1,N
          KLO=1
          KHI=N
 1        IF(KHI-KLO.GT.1) THEN
              K=(KHI+KLO)/2
              IF(XA(K).GT.X(I)) THEN
                  KHI=K
              ELSE
                  KLO=K
              ENDIF
              GOTO 1
          ENDIF
          H=XA(KHI)-XA(KLO)
          A=(XA(KHI)-X(I))/H
          B=(X(I)-XA(KLO))/H
          Y(I)=A*YA(KLO)+B*YA(KHI)+
     +        ((A**3-A)*WORK(KLO,2)+(B**3-B)*WORK(KHI,2))*(H**2)/6.0E0
      ENDDO
      IF(XA(1).GT.X(1)) Y(1)=YA(1)
      IF(XA(N).LT.X(N)) Y(N)=YA(N)
C
      RETURN
      END
C******************************************************************************
      SUBROUTINE SPLINE_PART(X,Y,N,NLOW,NHIGH,WORK,Y2)
C******************************************************************************
C     Routine to calculate 2.nd derivatives of tabulated function
C     Y(i)=Y(Xi), to be used for cubic spline calculation.
C
      IMPLICIT NONE
      INTEGER N,NLOW,NHIGH,I,NH,NL
      REAL(KIND=4)  X(N),Y(N),WORK(N),Y2(N)
      REAL(KIND=4)  SIG,P,QN,UN,YP1,YPN
C
      NL=NLOW+1
      NH=NHIGH-1
      YP1=(Y(NL)-Y(NLOW))/(X(NL)-X(NLOW))
      YPN=(Y(NHIGH)-Y(NH))/(X(NHIGH)-X(NH))
      IF(YP1.GT.99.0E+30) THEN
          Y2(NLOW)=0.0
          WORK(NLOW)=0.0
      ELSE
          Y2(NLOW)=-0.5E0
          WORK(NLOW)=(3.0E0/(X(NL)-X(NLOW)))*
     +               ((Y(NL)-Y(NLOW))/(X(NL)-X(NLOW))-YP1)
      ENDIF
      DO I=NL,NH
          SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
          P=SIG*Y2(I-1)+2.0E0
          Y2(I)=(SIG-1.0E0)/P
          WORK(I)=(6.0E0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     +             /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*WORK(I-1))/P
      ENDDO
      IF(YPN.GT.99.0E+30) THEN
          QN=0.0
          UN=0.0
      ELSE
          QN=0.5E0
          UN=(3.0E0/(X(NHIGH)-X(NH)))*
     +       (YPN-(Y(NHIGH)-Y(NH))/(X(NHIGH)-X(NH)))
      ENDIF
      Y2(NHIGH)=(UN-QN*WORK(NH))/(QN*Y2(NH)+1.0E0)
      DO I=NH,NLOW,-1
          Y2(I)=Y2(I)*Y2(I+1)+WORK(I)
      ENDDO
C
      RETURN
      END
C******************************************************************************
      SUBROUTINE SPLINT_PART(XA,YA,Y2A,N,NLOW,NHIGH,X,Y)
C******************************************************************************
C     Cubic spline calculation
C
      IMPLICIT NONE
      INTEGER KLO,KHI,N,K,NLOW,NHIGH
      REAL(KIND=4)  XA(N),YA(N),Y2A(N)
      REAL(KIND=4)  X,Y,H,A,B
C
      KLO=NLOW
      KHI=NHIGH
 1    IF(KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XA(K).GT.X) THEN
              KHI=K
          ELSE
              KLO=K
          ENDIF
          GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     +        ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.0E0
C
      RETURN
      END
C******************************************************************************
      SUBROUTINE LIPVEC(N,XA,YA,GROWTH,X,Y)
C******************************************************************************
C     Linearly interpolation on array:
C     Restriction
C     IF GROWTH
C     XA(I).LE. X(I) .LT.XA(I+1)
C
C     IF .NOT.GROWTH:
C     XA(I-1).LE. X(I) .LT.XA(I)
C
C
      IMPLICIT NONE
      INTEGER N,I,J
      LOGICAL GROWTH
      REAL(KIND=4)  H,A,B,X(N),Y(N),XA(N),YA(N)
C
      IF(GROWTH) THEN
          DO I=1,N-1
              J=I+1
              H=XA(J)-XA(I)
              A=(XA(J)-X(J))/H
              B=(X(J)-XA(I))/H
              Y(J)=A*YA(I)+B*YA(J)
          ENDDO
CC          Y(1)=YA(1)-(Y(2)-YA(2))
CC          Y(1)=YA(1)
            Y(1)=YA(1)+(X(1)-XA(1))*(YA(2)-YA(1))/(XA(2)-XA(1))
      ELSE
          DO I=1,N-1
              J=I+1
              H=XA(J)-XA(I)
              A=(XA(J)-X(I))/H
              B=(X(I)-XA(I))/H
              Y(I)=A*YA(I)+B*YA(J)
          ENDDO
CC          Y(N)=Y(N)-(Y(N-1)-YA(N-1))
CC          Y(N)=YA(N)
            Y(N)=YA(N)+(X(N)-XA(N))*(YA(N)-YA(N-1))/(XA(N)-XA(N-1))  
      ENDIF
C
      RETURN
      END
C******************************************************************************
      SUBROUTINE LININT(N,XA,YA,X,Y)
C******************************************************************************
C     Linearly interpolation:
C
      IMPLICIT NONE
      INTEGER N,I
      REAL(KIND=4)  X,Y,XA(N),YA(N)
      LOGICAL GROWTH
C
      GROWTH=XA(N).GT.XA(1)
      IF(GROWTH) THEN          
          IF(X.LE.XA(1)) THEN
              Y=YA(1)
              RETURN
          ELSE IF(X.GE.XA(N)) THEN
              Y=YA(N)
              RETURN    
          ELSE
              DO I=2,N
                IF(X.LT.XA(I)) THEN
                  Y=YA(I-1)+(YA(I)-YA(I-1))*(X-XA(I-1))/(XA(I)-XA(I-1))
                    RETURN
                ENDIF
              ENDDO
          ENDIF
      ELSE
          IF(X.GE.XA(1)) THEN
              Y=YA(1)
              RETURN
          ELSE IF(X.LE.XA(N)) THEN
              Y=YA(N)
              RETURN    
          ELSE
              DO I=2,N
                IF(X.GT.XA(I)) THEN
                   Y=YA(I-1)+(YA(I)-YA(I-1))*(X-XA(I-1))/(XA(I)-XA(I-1))
                    RETURN
                ENDIF
              ENDDO
          ENDIF
      ENDIF  
      END
C******************************************************************************
      REAL FUNCTION RTSAFE(FUNCD,X1,X2,XACC)
C     Newton Ralphson zero point iteration:
C******************************************************************************
      IMPLICIT NONE
      REAL(KIND=4)  X1,X2,XACC,FL,FH,DF,XL,XH,SWAP,DXOLD,DX,TEMP,F
      INTEGER MAXIT,J
      PARAMETER (MAXIT=100)
      CALL FUNCD(X1,FL,DF)
      CALL FUNCD(X2,FH,DF)
C      IF(FL*FH.GE.0.) PAUSE 'root must be bracketed'
      IF(FL*FH.GE.0.) THEN
          RTSAFE=X2
          RETURN
      ENDIF
      IF(FL.LT.0.)THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
        SWAP=FL
        FL=FH
        FH=SWAP
      ENDIF
      RTSAFE=.5E0*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL FUNCD(RTSAFE,F,DF)
      DO J=1,MAXIT
        IF(((RTSAFE-XH)*DF-F)*((RTSAFE-XL)*DF-F).GE.0.E0
     +      .OR. ABS(2.E0*F).GT.ABS(DXOLD*DF)) THEN
          DXOLD=DX
          DX=0.5E0*(XH-XL)
          RTSAFE=XL+DX
          IF(XL.EQ.RTSAFE)RETURN
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=RTSAFE
          RTSAFE=RTSAFE-DX
          IF(TEMP.EQ.RTSAFE)RETURN
        ENDIF
        IF(ABS(DX).LT.XACC) RETURN
        CALL FUNCD(RTSAFE,F,DF)
        IF(F.LT.0.) THEN
          XL=RTSAFE
          FL=F
        ELSE
          XH=RTSAFE
          FH=F
        ENDIF
      ENDDO
C      PAUSE 'RTSAFE exceeding maximum iterations'
      RETURN
      END
*****************************************************************************
*      SUBROUTINE OPTICAL_PARAMETERS(TAIR,PAIR,
*     +                    NBINS,ND,LAMBDA,ANGLE,
*     +                    DIF_CROSS,EXT_CROSS,DEPOL,
*     +                    AEROSOL_BACKSCATTER_RATIO,
*     +                    EXTINCTION,
*     +                    DEPOLARIZATION)
*****************************************************************************
*     Subroutine for repeated calculation of aerosol backscatter ratio, extinction,
*     and depolarisation for a given differrential size distribution and
*     air temperature/pressure.
*     SUBROUTINE CROSS_SECTIONS is supposed to have been called once prior to
*     this subroutine for calculations of backscattering and extinction cross
*     sections at the same particle radii, wavelength, angle, and refractivive
*     indices. The cross sections must not be changed by the calling program
*     between calls to this subroutine.
*
*     Input:
*       TAIR:           (K)         Ambient air temperature
*       PAIR:           (Pa)        Air pressure
*       NBINS                       Number of size bins
*       ND:     (par/kg air)   Number density of aerosol particles
*       LAMBDA          (m)         Wavelength of incident light
*       ANGLE           (deg)       Backscatter angle
*       DIF_CROSS       (m**2 sr-1) Differential backscattering cross section
*       EXT_CROSS       (m**2)      Extinction cross section
*       DEPOL                       Linear depolarisation ratio
*
*     Output
*       AEROSOL_BACKSCATTER_RATIO   Particulate to molecular ratio
*       EXTINCTION      (1/m)       Aerosol extinction coefficient
*       DEPOLARIZATION              Linear depolarisation ratio
*
*
*****************************************************************************
      SUBROUTINE OPTICAL_PARAMETERS(TAIR,PAIR,
     +                    NBINS,ND,LAMBDA,ANGLE,
     +                    DIF_CROSS,EXT_CROSS,DEPOL,
     +                    AEROSOL_BACKSCATTER_RATIO,
     +                    EXTINCTION,
     +                    DEPOLARIZATION)
      IMPLICIT NONE
      INTEGER NBINS,I
      REAL TAIR,PAIR,ND(NBINS),LAMBDA,ANGLE,
     +     DIF_CROSS(NBINS),EXT_CROSS(NBINS),DEPOL(NBINS),
     +     AEROSOL_BACKSCATTER_RATIO,EXTINCTION,DEPOLARIZATION
C     Auxilary variables:
      REAL RHOAIR,RAY,RAYLGH
C     Physical constants:
      REAL MAIR,RGAS,DEPOL_M
      PARAMETER(
     +     DEPOL_M=0.014,
     +     MAIR=28.9644E-3,
     +     RGAS=8.31441)
C
      RHOAIR=MAIR*PAIR/(RGAS*TAIR)
      RAY=RAYLGH(LAMBDA,ANGLE,PAIR,TAIR)
C
      AEROSOL_BACKSCATTER_RATIO=0.0
      EXTINCTION=0.0
      DEPOLARIZATION=0.0
      DO I=1,NBINS
          AEROSOL_BACKSCATTER_RATIO=AEROSOL_BACKSCATTER_RATIO+
     +                              ND(I)*DIF_CROSS(I)
          EXTINCTION=EXTINCTION+ND(I)*EXT_CROSS(I)
          DEPOLARIZATION=DEPOLARIZATION+ND(I)*DEPOL(I)
      END DO

      DEPOLARIZATION=DEPOLARIZATION/AEROSOL_BACKSCATTER_RATIO
      AEROSOL_BACKSCATTER_RATIO=RHOAIR*
     +        AEROSOL_BACKSCATTER_RATIO/RAY
      DEPOLARIZATION=100.0*(DEPOLARIZATION*
     +               AEROSOL_BACKSCATTER_RATIO+DEPOL_M)/
     +          (AEROSOL_BACKSCATTER_RATIO+1.0)
      EXTINCTION=EXTINCTION*RHOAIR
C
      RETURN
      END SUBROUTINE
*****************************************************************************
*      SUBROUTINE CROSS_SECTIONS(
*     +                    T_MATRIX_DIR,NBINS,RADIUS,
*     +                    ASPECT,REF_IND,LAMBDA,ANGLE,
*     +                    DIF_CROSS,EXT_CROSS,DEPOL)
*****************************************************************************
*     Subroutine for Mie scattering (spherical partticles OR T-matrix calculations,
*     using database  of expansions of the elements of the scattering matrix in
*     generalised spherical functions (ALPHA1..ALPHA4, BETA1,BETA2 in Mishchenko's
*     notation).
*     The subroutine calculates the differential back scattering cross section
*     extinction cross section and depolarisation for a given wavelength, angel,
*     refractive index, and aspect ratio.
*     The recommendations of Mishchenko for input parameter settings have been
*     followed. The data base is calculated for volume-equivalent sizes
*     (RAT=1D0  in  Mishchenko's notation) and imaginary refractive index 1E-8.
*     See further description in original Mishchenko code:
*     Mishchenko and Travis, J.Quant.Spectrosc.Radiat.Transfer 60,no.3,309,1998.
*     The data base files are assumed to reside in a directory gives as input.
*
*
*     Input:
*       T_MATRIX_DIR                Character string holding the name of the
*                                   directory with the T-matrix database
*       NBINS                       Number of size bins
*       RADIUS          (m)         Equivalent sphere radius
*       ASPECT                      Ratio of horizontal to rotational axis
*                                   (0.5 - 2.0)
*       REF_IND                     Real part of refractive index
*                                   (1.29 - 1.61)
*       LAMBDA          (m)         Wavelength of incident light
*       ANGLE           (deg)       Backscatter angle
*
*     Output
*       DIF_CROSS       (m**2 sr-1) Differential backscattering cross section
*       EXT_CROSS       (m**2)      Extinction cross section
*       DEPOL                       Linear depolarisation ratio
*****************************************************************************
      SUBROUTINE CROSS_SECTIONS(
     +                    T_MATRIX_DIR,NBINS,RADIUS,
     +                    ASPECT,REF_IND,LAMBDA,ANGLE,
     +                    DIF_CROSS,EXT_CROSS,DEPOL)
      IMPLICIT NONE
      CHARACTER *(*) T_MATRIX_DIR
      INTEGER NBINS,I
      REAL ASPECT,REF_IND,LAMBDA,ANGLE,RADIUS(NBINS),
     +     DIF_CROSS(NBINS),EXT_CROSS(NBINS),DEPOL(NBINS)
C
      IF(ABS(ASPECT-1.0E0).LT.1.0E-5) THEN
C         Spherical particles; Mie scattering
          DO I=1,NBINS
             CALL MIEBCK(LAMBDA,RADIUS(I),
     +                   ANGLE,
     +                   REF_IND,1.0E-8,
     +                   DIF_CROSS(I))
             CALL MIEEXT(LAMBDA,RADIUS(I),
     +                   REF_IND,1.0E-8,
     +                   EXT_CROSS(I))
             DEPOL(I)=0.0
          END DO
      ELSE
C         Aspherical particles; T-matric calculations
          CALL T_MATRIX(T_MATRIX_DIR,NBINS,RADIUS,
     +                  ASPECT,
     +                  REF_IND,
     +                  LAMBDA,ANGLE,
     +                  DIF_CROSS,
     +                  EXT_CROSS,
     +                  DEPOL)
      ENDIF
      RETURN
      END SUBROUTINE
*****************************************************************************
*     SUBROUTINE T_MATRIX(T_MATRIX_DIR,NBINS,RADIUS,
*    +                    ASPECT,REF_IND,LAMBDA,ANGLE,
*    +                    DIF_CROSS,EXT_CROSS,DEPOL)
*
*****************************************************************************
*     Subroutine for T-matrix calculations, using database of expansions
*     of the elements of the scattering matrix in generalised spherical
*     functions (ALPHA1..ALPHA4, BETA1,BETA2 in Mishchenko's notation).
*     The recommendations of Mishchenko for input parameter settings have been
*     followed. The data base is calculated for volume-equivalent sizes
*     (RAT=1D0  in  Mishchenko's notation) and imaginary refractive index 1E-8.
*     See further description in original Mishchenko code:
*     Mishchenko and Travis, J.Quant.Spectrosc.Radiat.Transfer 60,no.3,309,1998.
*     The data base files are assumed to reside in a subdirectory of the current
*     directory .\t-matrix (this can be changed in the code below).
*
*     The subroutine also internally calculates the scattering matrix elements
*     (SCATTER), scattering cross section (QS) (m**2), single scattering albedo
*     (SA), and the asymmetry factor (AS) which can easily be output
*
*     Input (reals: double precision variables - REAL*8):
*       T_MATRIX_DIR                Character string holding the name of the
*                                   directory with the T-matrix database
*       NBINS                       Number of size bins
*       RADIUS          (m)         Equivalent sphere radius
*       ASPECT                      Ratio of horizontal to rotational axis
*                                   (0.5 - 2.0)
*       REF_IND                     Real part of refractive index
*                                   (1.29 - 1.61)
*       LAMBDA          (m)         Wavelength of incident light
*       ANGLE           (deg)       Backscatter angle
*
*     Output
*       DIF_CROSS       (m**2 sr-1) Differential backscattering cross section
*       EXT_CROSS       (m**2)      Extinction cross section
*       DEPOL                       Linear depolarisation ratio
*****************************************************************************
      SUBROUTINE T_MATRIX(T_MATRIX_DIR,NBINS,RADIUS,
     +                    ASPECT,REF_IND,LAMBDA,ANGLE,
     +                    DIF_CROSS,EXT_CROSS,DEPOL)
C----------------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER *(*) T_MATRIX_DIR
      INTEGER NBINS
      REAL ASPECT,REF_IND,LAMBDA,ANGLE,RADIUS(NBINS),
     +     DIF_CROSS(NBINS),EXT_CROSS(NBINS),DEPOL(NBINS)
CC    REAL*8 ASPECT,REF_IND,LAMBDA,ANGLE,RADIUS(NBINS),
CC   +     DIF_CROSS(NBINS),EXT_CROSS(NBINS),DEPOL(NBINS)
C
      INTEGER N_ASPECT,N_INDEX,IA,IR,MAXRAD,NPL,I,J,J_1,NRAD,K,LMAX
      PARAMETER (N_ASPECT=31,N_INDEX=17,MAXRAD=200,NPL=361)
      REAL*8 ASPECT_ARRAY(N_ASPECT), R_INDEX(N_INDEX),RAD(MAXRAD),
     +     TOT_EXT(MAXRAD),SCAT_EXT(MAXRAD),SCAT_ALBEDO(MAXRAD),
     +     ASYM_PAR(MAXRAD),ANG,SCATTER(6),DIFF,PI,QE,QS,SA,AS,
     +     RL,F
      REAL*8
     +     ALPH1(NPL,MAXRAD),ALPH2(NPL,MAXRAD),ALPH3(NPL,MAXRAD),
     +     ALPH4(NPL,MAXRAD),
     +     BET1(NPL,MAXRAD),BET2(NPL,MAXRAD),
     +     A1(NPL),A2(NPL),A3(NPL),A4(NPL),B1(NPL),B2(NPL)

      INTEGER NCOEF(MAXRAD)
      CHARACTER FILNAM*10, FRAME*100
      DATA FILNAM/'Aaaa.nrrrE'/
c
C     Aspect ratios:
      DATA ASPECT_ARRAY/
     + 0.50D0,0.55D0,
     + 0.60D0,0.65D0,
     + 0.70D0,0.75D0,
     + 0.80D0,0.85D0,
     + 0.90D0,0.95D0,
     + 1.0000000001D0,1.05D0,
     + 1.10D0,1.15D0,
     + 1.20D0,1.25D0,
     + 1.30D0,1.35D0,
     + 1.40D0,1.45D0,
     + 1.50D0,1.55D0,
     + 1.60D0,1.65D0,
     + 1.70D0,1.75D0,
     + 1.80D0,1.85D0,
     + 1.90D0,1.95D0,
     + 2.00D0/
C----------------------------------------------------------------------------
C     Real part of refractive indices:
      DATA R_INDEX/ 1.29,
     +1.31D0,1.33D0,1.35D0,
     +1.37D0,1.39D0,
     +1.41D0,1.43D0,1.45D0,
     +1.47D0,1.49d0,
     +1.51D0,1.53D0,1.55D0,
     +1.57D0,1.59D0,1.61D0/

      ANG=DBLE(ANGLE)
      DIFF=1.0E20
      K=1
      DO IA=1,N_ASPECT
          IF(ABS(ASPECT_ARRAY(IA)-DBLE(ASPECT)).LT.DIFF) THEN
               DIFF=DMIN1(ABS(ASPECT_ARRAY(IA)-DBLE(ASPECT)),DIFF)
               K=IA
          ENDIF
      ENDDO
      IA=K
C
      DIFF=1.0E20
      K=1
      DO IR=1,N_INDEX
          IF(ABS(R_INDEX(IR)-REF_IND).LT.DIFF) THEN
               DIFF=DMIN1(ABS(R_INDEX(IR)-DBLE(REF_IND)),DIFF)
               K=IR
          ENDIF
      ENDDO
      IR=K
C
      K=NINT(ASPECT_ARRAY(IA)*100.0)
      WRITE(FILNAM(2:4),'(I3.3)') K
      K=NINT(R_INDEX(IR)*100.0)
      WRITE(FILNAM(7:9),'(I3.3)') K
c      WRITE(*,*) T_MATRIX_DIR//FILNAM
c      PAUSE
      OPEN(1,FILE=T_MATRIX_DIR//FILNAM)
C
      PI=DACOS(-1.0D0)
      I=1
 10   CONTINUE
      READ(1,'(A)',ERR=200,END=100) FRAME
c      write(*,*) frame
      READ(FRAME(9:17),'(F9.4)') RAD(I)
      READ(1,*) TOT_EXT(I),SCAT_EXT(I),SCAT_ALBEDO(I),
     +          ASYM_PAR(I)
      TOT_EXT(I)=TOT_EXT(I)/(2.0D0*PI*RAD(I)*RAD(I))
      SCAT_EXT(I)=SCAT_EXT(I)/(2.0D0*PI*RAD(I)*RAD(I))
      READ(1,*) NCOEF(I)
      READ(1,*) (ALPH1(J,I),J=1,NCOEF(I))
      READ(1,*) (ALPH2(J,I),J=1,NCOEF(I))
      READ(1,*) (ALPH3(J,I),J=1,NCOEF(I))
      READ(1,*) (ALPH4(J,I),J=1,NCOEF(I))
      READ(1,*) (BET1(J,I),J=1,NCOEF(I))
      READ(1,*) (BET2(J,I),J=1,NCOEF(I))
c      WRITE(*,*) I,RAD(I),NCOEF(I)
      I=I+1
      GOTO 10
 100  CLOSE(1)
      NRAD=I-1
c      WRITE(*,*) NRAD,frame
c      PAUSE

      DO I=1,NBINS
          RL=DBLE(RADIUS(I)/LAMBDA)
          IF(RL.LE.RAD(1)) THEN
              QE=TOT_EXT(1)*(RL/RAD(1))**4
              QS=SCAT_EXT(1)*(RL/RAD(1))**4
              QE=QE*2.0D0*PI*RL*RL
              QS=QS*2.0D0*PI*RL*RL
              SA=SCAT_ALBEDO(1)
              AS=ASYM_PAR(1)
              DO K=1,NCOEF(1)
                  A1(K)=ALPH1(K,1)
                  A2(K)=ALPH2(K,1)
                  A3(K)=ALPH3(K,1)
                  A4(K)=ALPH4(K,1)
                  B1(K)=BET1(K,1)
                  B2(K)=BET2(K,1)
              END DO
              LMAX=NCOEF(1)
          ELSE IF(RL.GE.RAD(NRAD)) THEN
              QE=TOT_EXT(NRAD)
              QS=SCAT_EXT(NRAD)
              QE=QE*2.0D0*PI*RL*RL
              QS=QS*2.0D0*PI*RL*RL
              SA=SCAT_ALBEDO(NRAD)
              AS=ASYM_PAR(NRAD)
              DO K=1,NCOEF(NRAD)
                  A1(K)=ALPH1(K,NRAD)
                  A2(K)=ALPH2(K,NRAD)
                  A3(K)=ALPH3(K,NRAD)
                  A4(K)=ALPH4(K,NRAD)
                  B1(K)=BET1(K,NRAD)
                  B2(K)=BET2(K,NRAD)
              END DO
              LMAX=NCOEF(NRAD)
          ELSE
              J=1
              DO WHILE (RAD(J).LT.RL)
                  J=J+1
              END DO
              J_1=J-1
C              WRITE(*,*) J_1,J,RAD(J_1),RL,RAD(J)
              F=(RL-RAD(J_1))/(RAD(J)-RAD(J_1))
              QE=TOT_EXT(J_1)*(1.0D0-F)+TOT_EXT(J)*F
              QS=SCAT_EXT(J_1)*(1.0D0-F)+SCAT_EXT(J)*F
              QE=QE*2.0D0*PI*RL*RL
              QS=QS*2.0D0*PI*RL*RL
              SA=SCAT_ALBEDO(J_1)*(1.0D0-F)+SCAT_ALBEDO(J)*F
              AS=ASYM_PAR(J_1)*(1.0D0-F)+ASYM_PAR(J)*F
              DO K=1,MIN(NCOEF(J_1),NCOEF(J))
                  A1(K)=ALPH1(K,J_1)*(1.0D0-F)+ALPH1(K,J)*F
                  A2(K)=ALPH2(K,J_1)*(1.0D0-F)+ALPH2(K,J)*F
                  A3(K)=ALPH3(K,J_1)*(1.0D0-F)+ALPH3(K,J)*F
                  A4(K)=ALPH4(K,J_1)*(1.0D0-F)+ALPH4(K,J)*F
                  B1(K)=BET1(K,J_1)*(1.0D0-F)+BET1(K,J)*F
                  B2(K)=BET2(K,J_1)*(1.0D0-F)+BET2(K,J)*F
              END DO
              LMAX=MIN(NCOEF(J_1),NCOEF(J))
           ENDIF
           CALL MY_SCATTER(A1,A2,A3,A4,B1,B2,NPL,
     +                        LMAX-1,
     +                        ANG,SCATTER)
           DIF_CROSS(I)=REAL(QS*SCATTER(1))*LAMBDA*LAMBDA/
     +                      (4.0E0*REAL(PI))
           EXT_CROSS(I)=REAL(QE)*LAMBDA*LAMBDA
           DEPOL(I)=REAL(DIF_CROSS(I)*(SCATTER(1)-SCATTER(2))/
     +                   (SCATTER(1)+2.0D0*SCATTER(5)+SCATTER(2)))
      END DO
      RETURN
 200  CONTINUE
      WRITE(*,*) 'Error reading T-matrix file: ',T_MATRIX_DIR//FILNAM
      STOP
      END SUBROUTINE
C****************************************************************
 
C    CALCULATION OF THE SCATTERING MATRIX FOR GIVEN EXPANSION
C    COEFFICIENTS
 
C    A1,...,B2 - EXPANSION COEFFICIENTS
C    LMAX - NUMBER OF COEFFICIENTS MINUS 1
C    N - NUMBER OF SCATTERING ANGLES
C        THE CORRESPONDING SCATTERING ANGLES ARE GIVEN BY
C        180*(I-1)/(N-1) (DEGREES), WHERE I NUMBERS THE ANGLES
 
      SUBROUTINE MY_SCATTER(A1,A2,A3,A4,B1,B2,NPL,LMAX,
     +                      ANGLE,SCATTER)
C      INCLUDE 'tmq.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NPL,LMAX,L1MAX,L1,L
      REAL*8 A1(NPL),A2(NPL),A3(NPL),A4(NPL),B1(NPL),B2(NPL)
      REAL*8 ANGLE,SCATTER(6)
      L1MAX=LMAX+1
CC    DO 10 L1=1,L1MAX
CC       L=L1-1
CC 10 CONTINUE
      TAA=ANGLE*DACOS(-1D0)/180.0D0
      D6=DSQRT(6D0)*0.25D0
         U=DCOS(TAA)
         F11=0D0
         F2=0D0
         F3=0D0
         F44=0D0
         F12=0D0
         F34=0D0
         P1=0D0
         P2=0D0
         P3=0D0
         P4=0D0
         PP1=1D0
         PP2=0.25D0*(1D0+U)*(1D0+U)
         PP3=0.25D0*(1D0-U)*(1D0-U)
         PP4=D6*(U*U-1D0)
         DO 400 L1=1,L1MAX
            L=L1-1
            DL=DFLOAT(L)
            DL1=DFLOAT(L1)
            F11=F11+A1(L1)*PP1
            F44=F44+A4(L1)*PP1
            IF(L.EQ.LMAX) GO TO 350
            PL1=DFLOAT(2*L+1)
            P=(PL1*U*PP1-DL*P1)/DL1
            P1=PP1
            PP1=P
  350       IF(L.LT.2) GO TO 400
            F2=F2+(A2(L1)+A3(L1))*PP2
            F3=F3+(A2(L1)-A3(L1))*PP3
            F12=F12+B1(L1)*PP4
            F34=F34+B2(L1)*PP4
            IF(L.EQ.LMAX) GO TO 400
            PL2=DFLOAT(L*L1)*U
            PL3=DFLOAT(L1*(L*L-4))
            PL4=1D0/DFLOAT(L*(L1*L1-4))
            P=(PL1*(PL2-4D0)*PP2-PL3*P2)*PL4
            P2=PP2
            PP2=P
            P=(PL1*(PL2+4D0)*PP3-PL3*P3)*PL4
            P3=PP3
            PP3=P
            P=(PL1*U*PP4-DSQRT(DFLOAT(L*L-4))*P4)/DSQRT(DFLOAT(L1*L1-4))
            P4=PP4
            PP4=P
  400    CONTINUE
         F22=(F2+F3)*0.5D0
         F33=(F2-F3)*0.5D0
         SCATTER(1)=F11
         SCATTER(2)=F22
         SCATTER(3)=F33
         SCATTER(4)=F44
         SCATTER(5)=F12
         SCATTER(6)=F34
C        F22=F22/F11
C        F33=F33/F11
C        F44=F44/F11
C        F12=-F12/F11
C        F34=F34/F11
      RETURN
      END
*****************************************************************************
*     SUBROUTINE SAOPT(OPTICAL_PAR,NBINS,PND,PTSIZE,PAIR,TAIR,
*                      WAVEL,REF,OPT,RAY)
*****************************************************************************
*
*     CHARACTER*1 OPTICAL_PAR
*     INTEGER NBINS
*     REAL    PND(NBINS),PTSIZE(NBINS,3),
*             PAIR,TAIR,WAVEL,REFRE(NBINS),BKP,RAY
*
*     Input:
*         OPTICAL_PAR             Character flag, indicating calculation of
*                                 Molecular volume backscatter coefficient ('B') or
*                                 Aerosol extinction coefficient ('E')
*         NBINS:                  Number of particle radii bins
*         PND:     (par/kg air)   Number density of aerosol particles 
*         PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume, 
*                                 (not to be changed)
*         PAIR:    (Pa)           Air pressure
*         TAIR:    (K)            Ambient air temperature 
*         WAVEL:   (m)            Incident wavelength 
*         REF:                    Real part of refractive index
*
*     Output:
*         OPT:     (1/m sr)       Particle volume backscatter coefficient OR
*                  (1/m)          aerosol extinction coefficient
*         RAY:     (1/m sr)       Molecular volume backscatter coefficient
*
*
*     Subroutine used to calculate the particle volume backscatter coefficient
*     OR extinction coefficient of an ensemble of spherical particles, characterised
*     by a differential size distribution PND (calculated by SUBROUTINE PSCBOX),
*     at a specified incident wavelength (WAVEL). The particles are assumed to be
*     characterised by the refractive index (real part) (REF), whereas the imaginary
*     part is given as a parameter in the subroutine.
*     Depending on the value of the character flag OPTICAL_PAR the backscatter
*     coefficient or extinction coefficient is calculated.
*
*     Also the Rayleigh (molecular) volume backscatter coefficient (RAY) is
*     calculated at the same wavelength at the specified ambient atmospheric 
*     state, given by the air pressure (PAIR) and temperature (TAIR).
*
*     The subroutine is intended to be used together with SUBROUTINE PSCBOX
*     which calculates the size distribution of liquid stratospheric aerosols.
*     Subroutine SETBIN from the PSCBOX package must be called once prior to
*     any call to SUBROUTINE SAOPT.
*
*     An angel distribution of the backscattered light is assumed, 
*     corresponding to the sensitivity angel distribution of the University
*     of Wyoming backscatter sonde. Alternatively, the subroutine can easily
*     be changed by logical parameter SINGLE_ANGEL to calculate the backscatter
*     at a single specific angle.
*
*
      SUBROUTINE SAOPT(OPTICAL_PAR,NBINS,PND,PTSIZE,PAIR,TAIR,
     +                 WAVEL,REF,OPT,RAY)
C
      IMPLICIT NONE
      CHARACTER*1 OPTICAL_PAR
      INTEGER NBINS
      REAL PND(NBINS),PTSIZE(NBINS,3),PAIR,TAIR,WAVEL,REF(NBINS),OPT,RAY
C
C     Logical flag to specify if backscatter coefficient is to be calculated
C     at a specific ANGEL (.true.) or at the U. Wyoming backscatter sonde
C     angel distribution:
      LOGICAL SINGLE_ANGEL
C      PARAMETER (SINGLE_ANGEL=.TRUE.)
      PARAMETER (SINGLE_ANGEL=.FALSE.)
C
      REAL ANGEL
      PARAMETER(ANGEL=180.0)
c      PARAMETER(ANGEL=173.0)
C
C     Angel distribition of backscatter signal for the 
C     University of Wyoming backscatter sonde:
      INTEGER NANG
      PARAMETER(NANG=31)
      REAL ANGELS(NANG),ANGELW(NANG),W
      DATA ANGELW/
     +.022,.053,.058,.061,.064,.064,.064,.062,.059,.056,.052,.049,.044,
     +.041,.037,.033,.029,.025,.022,.019,.016,.014,.011,.009,.007,.006,
     +        .005,.004,.003,.002,.001/
      DATA ANGELS/
     +        180.,179.,178.,177.,176.,175.,174.,173.,172.,171.,
     +        170.,169.,168.,167.,166.,165.,164.,163.,162.,161.,
     +        160.,159.,158.,157.,156.,155.,154.,153.,152.,151.,150./
C
C     Imaginary part of refractive index:
      REAL REFIM
      PARAMETER(REFIM=1.0E-7)
C
      REAL MAIR,RGAS,RHOAIR
      PARAMETER(RGAS=8.31441,MAIR=28.9644E-3)
      INTEGER I,J
      REAL RAYLGH,M,MM,X
      SAVE ANGELW,ANGELS
      LOGICAL BACKSCAT
C
      BACKSCAT=INDEX('Bb',OPTICAL_PAR(1:1)).NE.0
      RHOAIR=MAIR*PAIR/(RGAS*TAIR)
C
      IF(SINGLE_ANGEL) THEN
C         Calculate Rayleigh backscattering:
          RAY=RAYLGH(WAVEL,ANGEL,PAIR,TAIR)
          X=0.0
C         Calculate Mie backscattering or extinction:
          DO I=1,NBINS
              IF(REF(I).GT.0.0) THEN
                  IF(BACKSCAT) THEN
                     CALL MIEBCK(WAVEL,PTSIZE(I,1),ANGEL,REF(I),REFIM,M)
                  ELSE
                     CALL MIEEXT(WAVEL,PTSIZE(I,1),REF(I),REFIM,M)
                  ENDIF
                  X=X+PND(I)*M
              ENDIF
          ENDDO
          OPT=X*RHOAIR
      ELSE
C         Calculate Rayleigh backscattering:
          X=0.0
          W=0.0
          DO J=1,NANG
              X=X+ANGELW(J)*RAYLGH(WAVEL,ANGELS(J),PAIR,TAIR)
              W=W+ANGELW(J)
          ENDDO
          RAY=X/W
C
          X=0.0
          DO I=1,NBINS
              IF(BACKSCAT) THEN
C                 Calculate Mie backscattering:
                  IF(REF(I).GT.0.0) THEN
                      MM=0.0
                      DO J=1,NANG
                          CALL MIEBCK(WAVEL,PTSIZE(I,1),ANGELS(J),
     +                                REF(I),REFIM,M)
                          MM=MM+M*ANGELW(J)
                      ENDDO
                      M=MM/W
                      X=X+PND(I)*M
                  ENDIF
              ELSE
C                 Calculate extinction:
                  IF(REF(I).GT.0.0) THEN
                      CALL MIEEXT(WAVEL,PTSIZE(I,1),REF(I),REFIM,M)
                      X=X+PND(I)*M
                  ENDIF
              ENDIF
          ENDDO
          OPT=X*RHOAIR
      ENDIF
C
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE MIEBCK(WAVEL,RADIUS,ANGEL,REFRE,REFIM,SBACK)
C-----------------------------------------------------------------------
C
C     This subroutine calculates the differential back scattering
C     cross section for a given wavelength, particle radius, angel,
C     and refractive index for a spherical particle.
C
C     The scattering matrix is calculated internally in the
C     subroutine and could be specified as output variable if desired.
C
C
C     Reference: Absorption & Scattering of light by small particles.
C                Craig F. Bohren & Donald R. Huffman, John Wiley, 1983.
C                Appendix A.
C
C     Input:
C         WAVEL:              Incident wavelength (same unit as RADIUS)
C         RADIUS:             Spherical particle radius (same unit as WAVEL)
C         ANGEL:  (deg)       Backscattering angel from forward direction
C         REFRE:              Real part of refractive index
C         REFIM:              Imaginary part of refractive index
C     Output:
C         SBACK:  (m**2/sr)   Differential backscattering cross section
C
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
      REAL AMU,THETA,PI,TAU,PI0,PI1,REFRE,REFIM,
     +        WAVEL,ANGEL,X,REFMED,
     +        RADIUS,PIMATH,
     +        XSTOP,YMOD,RN,CHI0,CHI1,APSI1,FN,APSI,CHI,SBACK
      INTEGER NSTOP,NMX,NN,N
      COMPLEX D(5000),Y,REFREL,XI,XI1,AN,BN,S1
CC      COMPLEX S2
      DOUBLE PRECISION PSI0,PSI1,PSI,DN,DX
C
C     Refractive index of medium:
      PARAMETER(REFMED=1.0)
C-----------------------------------------------------------------------
C
      PIMATH=4.0*ATAN(1.0)
      THETA=PIMATH*ANGEL/180.0
      REFREL=CMPLX(REFRE,REFIM)/REFMED
      X=2.0*PIMATH*RADIUS*REFMED/WAVEL
C
C-----------------------------------------------------------------------
C
      DX=X
      Y=X*REFREL
C
C-----------------------------------------------------------------------
C
C     Series terminated after NSTOP terms.
C
C-----------------------------------------------------------------------
C
      XSTOP=X+4.0*X**.3333+2.0
      NSTOP=XSTOP
      YMOD=CABS(Y)
      NMX=AMAX1(XSTOP,YMOD)+15
C
C-----------------------------------------------------------------------
C
      AMU=COS(THETA)
C
C-----------------------------------------------------------------------
C
C     Logarithmic derivative D(J) calculated by downward recurrence
C     beginning with initial value 0.0+I*0.0 at J = NMX.
C
C-----------------------------------------------------------------------
C
      D(NMX)=CMPLX(0.0,0.0)
      NN=NMX-1
C
C-----------------------------------------------------------------------
C
      DO 110 N=1,NN
          RN=NMX-N+1
          D(NMX-N)=(RN/Y)-(1.0/(D(NMX-N+1)+RN/Y))
110   CONTINUE
C
C-----------------------------------------------------------------------
C
      PI0=0.0
      PI1=1.0
C
C-----------------------------------------------------------------------
C
      S1=CMPLX(0.0,0.0)
CC      S2=CMPLX(0.0,0.0)
C
C-----------------------------------------------------------------------
C
C     Riccati-Bessel functions with real arguement X calculated by
C     upward recurrence.
C
C-----------------------------------------------------------------------
C
      PSI0=DCOS(DX)
      PSI1=DSIN(DX)
      CHI0=-SIN(X)
      CHI1=COS(X)
      APSI1=PSI1
      XI1=CMPLX(APSI1,-CHI1)
      N=1
C
C-----------------------------------------------------------------------
C
      DO 200 N=1,NSTOP
C
          DN=N
          RN=N
          FN=(2.0*RN+1.0)/(RN*(RN+1.0))
          PSI=(2.0*DN-1.0)*PSI1/DX-PSI0
          APSI=PSI
          CHI=(2.0*RN-1.)*CHI1/X-CHI0
          XI=CMPLX(APSI,-CHI)
          AN=(D(N)/REFREL+RN/X)*APSI-APSI1
          AN=AN/((D(N)/REFREL+RN/X)*XI-XI1)
          BN=(REFREL*D(N)+RN/X)*APSI-APSI1
          BN=BN/((REFREL*D(N)+RN/X)*XI-XI1)
C
C-----------------------------------------------------------------------
C
          PI=PI1
          TAU=RN*AMU*PI-(RN+1.)*PI0
          S1=S1+FN*(AN*PI+BN*TAU)
CC          S2=S2+FN*(AN*TAU+BN*PI)
C
C-----------------------------------------------------------------------
C
          PSI0=PSI1
          PSI1=PSI
          APSI1=PSI1
          CHI0=CHI1
          CHI1=CHI
          XI1=CMPLX(APSI1,-CHI1)
          RN=N+1
C
C-----------------------------------------------------------------------
C
          PI1=((2.0*RN-1.0)/(RN-1.0))*AMU*PI
          PI1=PI1-RN*PI0/(RN-1.0)
          PI0=PI
C
C-----------------------------------------------------------------------
C
 200  CONTINUE
C
C-----------------------------------------------------------------------
C
      SBACK=CABS(S1)*CABS(S1)*RADIUS*RADIUS/(X*X)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE MIEEXT(WAVEL,RADIUS,REFRE,REFIM,SEXT)
C-----------------------------------------------------------------------
C
C     This subroutine calculates the extinction cross section 
C     for a given wavelength, particle radius,
C     and refractive index for a spherical particle.
C
C     The scattering matrix is calculated internally in the
C     subroutine and could be specified as output variable if desired.
C
C
C     Reference: Absorption & Scattering of light by small particles.
C                Craig F. Bohren & Donald R. Huffman, John Wiley, 1983.
C                Appendix A.
C
C-----------------------------------------------------------------------
C     Input:
C         WAVEL:              Incident wavelength (same unit as RADIUS)
C         RADIUS:             Spherical particle radius (same unit as WAVEL)
C         REFRE:              Real part of refractive index
C         REFIM:              Imaginary part of refractive index
C     Output:
C         SEXT:   (m**2)      Extinction cross section
C
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
      REAL PI,TAU,PI0,PI1,REFRE,REFIM,
     +        WAVEL,X,REFMED,
     +        RADIUS,PIMATH,
     +        XSTOP,YMOD,RN,CHI0,CHI1,APSI1,FN,APSI,CHI,SEXT
      INTEGER NSTOP,NMX,NN,N
      COMPLEX D(5000),Y,REFREL,XI,XI1,AN,BN,S1
CC      COMPLEX S2
      DOUBLE PRECISION PSI0,PSI1,PSI,DN,DX
C
C     Refractive index of medium:
      PARAMETER(REFMED=1.0)
C-----------------------------------------------------------------------
C
      PIMATH=4.0*ATAN(1.0)
      REFREL=CMPLX(REFRE,REFIM)/REFMED
      X=2.0*PIMATH*RADIUS*REFMED/WAVEL
C
C-----------------------------------------------------------------------
C
      DX=X
      Y=X*REFREL
C
C-----------------------------------------------------------------------
C
C     Series terminated after NSTOP terms.
C
C-----------------------------------------------------------------------
C
      XSTOP=X+4.0*X**.3333+2.0
      NSTOP=XSTOP
      YMOD=CABS(Y)
      NMX=AMAX1(XSTOP,YMOD)+15
C
C-----------------------------------------------------------------------
C
C
C     Logarithmic derivative D(J) calculated by downward recurrence
C     beginning with initial value 0.0+I*0.0 at J = NMX.
C
C-----------------------------------------------------------------------
C
      D(NMX)=CMPLX(0.0,0.0)
      NN=NMX-1
C
C-----------------------------------------------------------------------
C
      DO 110 N=1,NN
          RN=NMX-N+1
          D(NMX-N)=(RN/Y)-(1.0/(D(NMX-N+1)+RN/Y))
110   CONTINUE
C
C-----------------------------------------------------------------------
C
      PI0=0.0
      PI1=1.0
C
C-----------------------------------------------------------------------
C
      S1=CMPLX(0.0,0.0)
CC      S2=CMPLX(0.0,0.0)
C
C-----------------------------------------------------------------------
C
C     Riccati-Bessel functions with real arguement X calculated by
C     upward recurrence.
C
C-----------------------------------------------------------------------
C
      PSI0=DCOS(DX)
      PSI1=DSIN(DX)
      CHI0=-SIN(X)
      CHI1=COS(X)
      APSI1=PSI1
      XI1=CMPLX(APSI1,-CHI1)
      N=1
C
C-----------------------------------------------------------------------
C
      DO 200 N=1,NSTOP
C
          DN=N
          RN=N
          FN=(2.0*RN+1.0)/(RN*(RN+1.0))
          PSI=(2.0*DN-1.0)*PSI1/DX-PSI0
          APSI=PSI
          CHI=(2.0*RN-1.)*CHI1/X-CHI0
          XI=CMPLX(APSI,-CHI)
          AN=(D(N)/REFREL+RN/X)*APSI-APSI1
          AN=AN/((D(N)/REFREL+RN/X)*XI-XI1)
          BN=(REFREL*D(N)+RN/X)*APSI-APSI1
          BN=BN/((REFREL*D(N)+RN/X)*XI-XI1)
C
C-----------------------------------------------------------------------
C
          PI=PI1
          TAU=RN*PI-(RN+1.)*PI0
          S1=S1+FN*(AN*PI+BN*TAU)
CC          S2=S2+FN*(AN*TAU+BN*PI)
C
C-----------------------------------------------------------------------
C
          PSI0=PSI1
          PSI1=PSI
          APSI1=PSI1
          CHI0=CHI1
          CHI1=CHI
          XI1=CMPLX(APSI1,-CHI1)
          RN=N+1
C
C-----------------------------------------------------------------------
C
          PI1=((2.0*RN-1.0)/(RN-1.0))*PI
          PI1=PI1-RN*PI0/(RN-1.0)
          PI0=PI
C
C-----------------------------------------------------------------------
C
 200  CONTINUE
C
C-----------------------------------------------------------------------
C
      SEXT=(4.0*PIMATH*RADIUS*RADIUS/(X*X))*REAL(S1)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MIEB_E(WAVEL,RADIUS,ANGEL,REFRE,REFIM,SBACK,SEXT)
C-----------------------------------------------------------------------
C
C     This subroutine calculates the differential back scattering cross section
C     and the extinction cross section
C     for a given wavelength, particle radius, angel,
C     and refractive index for a spherical particle.
C
C     The scattering matrix is calculated internally in the
C     subroutine and could be specified as output variable if desired.
C
C
C     Reference: Absorption & Scattering of light by small particles.
C                Craig F. Bohren & Donald R. Huffman, John Wiley, 1983.
C                Appendix A.
C
C     Input:
C         WAVEL:              Incident wavelength (same unit as RADIUS)
C         RADIUS:             Spherical particle radius (same unit as WAVEL)
C         ANGEL:  (deg)       Backscattering angel from forward direction
C         REFRE:              Real part of refractive index
C         REFIM:              Imaginary part of refractive index
C     Output:
C         SBACK:  (m**2/sr)   Differential backscattering cross section
C         SEXT:   (m**2)      Extinction cross section
C
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
      REAL AMU,THETA,
     +     BPI,BTAU,BPI0,BPI1,
     +     EPI,ETAU,EPI0,EPI1,
     +     REFRE,REFIM,
     +        WAVEL,ANGEL,X,REFMED,
     +        RADIUS,PIMATH,
     +        XSTOP,YMOD,RN,CHI0,CHI1,APSI1,FN,APSI,CHI,SBACK,SEXT
      INTEGER NSTOP,NMX,NN,N
      COMPLEX D(5000),Y,REFREL,XI,XI1,AN,BN,
     +        BS1,
     +        ES1
CC      COMPLEX S2
      DOUBLE PRECISION PSI0,PSI1,PSI,DN,DX
C
C     Refractive index of medium:
      PARAMETER(REFMED=1.0)
C-----------------------------------------------------------------------
C
      PIMATH=4.0*ATAN(1.0)
      THETA=PIMATH*ANGEL/180.0
      REFREL=CMPLX(REFRE,REFIM)/REFMED
      X=2.0*PIMATH*RADIUS*REFMED/WAVEL
C
C-----------------------------------------------------------------------
C
      DX=X
      Y=X*REFREL
C
C-----------------------------------------------------------------------
C
C     Series terminated after NSTOP terms.
C
C-----------------------------------------------------------------------
C
      XSTOP=X+4.0*X**.3333+2.0
      NSTOP=XSTOP
      YMOD=CABS(Y)
      NMX=AMAX1(XSTOP,YMOD)+15
C
C-----------------------------------------------------------------------
C
      AMU=COS(THETA)
C
C-----------------------------------------------------------------------
C
C     Logarithmic derivative D(J) calculated by downward recurrence
C     beginning with initial value 0.0+I*0.0 at J = NMX.
C
C-----------------------------------------------------------------------
C
      D(NMX)=CMPLX(0.0,0.0)
      NN=NMX-1
C
C-----------------------------------------------------------------------
C
      DO 110 N=1,NN
          RN=NMX-N+1
          D(NMX-N)=(RN/Y)-(1.0/(D(NMX-N+1)+RN/Y))
110   CONTINUE
C
C-----------------------------------------------------------------------
C
      BPI0=0.0
      BPI1=1.0

      EPI0=0.0
      EPI1=1.0
C
C-----------------------------------------------------------------------
C
      BS1=CMPLX(0.0,0.0)
C
      ES1=CMPLX(0.0,0.0)
CC      S2=CMPLX(0.0,0.0)
C
C-----------------------------------------------------------------------
C
C     Riccati-Bessel functions with real arguement X calculated by
C     upward recurrence.
C
C-----------------------------------------------------------------------
C
      PSI0=DCOS(DX)
      PSI1=DSIN(DX)
      CHI0=-SIN(X)
      CHI1=COS(X)
      APSI1=PSI1
      XI1=CMPLX(APSI1,-CHI1)
      N=1
C
C-----------------------------------------------------------------------
C
      DO 200 N=1,NSTOP
C
          DN=N
          RN=N
          FN=(2.0*RN+1.0)/(RN*(RN+1.0))
          PSI=(2.0*DN-1.0)*PSI1/DX-PSI0
          APSI=PSI
          CHI=(2.0*RN-1.)*CHI1/X-CHI0
          XI=CMPLX(APSI,-CHI)
          AN=(D(N)/REFREL+RN/X)*APSI-APSI1
          AN=AN/((D(N)/REFREL+RN/X)*XI-XI1)
          BN=(REFREL*D(N)+RN/X)*APSI-APSI1
          BN=BN/((REFREL*D(N)+RN/X)*XI-XI1)
C
C-----------------------------------------------------------------------
C
          BPI=BPI1
          BTAU=RN*AMU*BPI-(RN+1.)*BPI0
          BS1=BS1+FN*(AN*BPI+BN*BTAU)
C
          EPI=EPI1
          ETAU=RN*AMU*EPI-(RN+1.)*EPI0
          ES1=ES1+FN*(AN*EPI+BN*ETAU)
CC          S2=S2+FN*(AN*TAU+BN*PI)
C
C-----------------------------------------------------------------------
C
          PSI0=PSI1
          PSI1=PSI
          APSI1=PSI1
          CHI0=CHI1
          CHI1=CHI
          XI1=CMPLX(APSI1,-CHI1)
          RN=N+1
C
C-----------------------------------------------------------------------
C
          BPI1=((2.0*RN-1.0)/(RN-1.0))*AMU*BPI
          BPI1=BPI1-RN*BPI0/(RN-1.0)
          BPI0=BPI
C
          EPI1=((2.0*RN-1.0)/(RN-1.0))*AMU*EPI
          EPI1=EPI1-RN*EPI0/(RN-1.0)
          EPI0=EPI
C
C-----------------------------------------------------------------------
C
 200  CONTINUE
C
C-----------------------------------------------------------------------
C
      SBACK=CABS(BS1)*CABS(BS1)*RADIUS*RADIUS/(X*X)
      SEXT=(4.0*PIMATH*RADIUS*RADIUS/(X*X))*REAL(ES1)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
C
C-----------------------------------------------------------------------
      REAL FUNCTION RAYLGH(WAVEL,ANGEL,PAIR,TAIR)
C-----------------------------------------------------------------------
C
C     This function calculates the molecular Rayleigh backscattering
C     coefficient.
C
C     Source: Nicolet et al.: Planet Sp. Sci. 30,935,1982.
C
C     Input:
C         WAVEL:   (m)            Incident wavelength
C         ANGEL:   (deg)          Backscattering angel from 
C                                     forward direction
C         PAIR:    (Pa)           Ambient air pressure
C         TAIR     (K)            Temperature (K)
C
C     Output:
C         RAYLGH   (1/m sr)       Rayleigh backscattering coefficient
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL WAVEL,ANGEL,TAIR,PAIR
      REAL KBOLTZ,PI,C,A1,A2,A3,X,Y,W
      PARAMETER(KBOLTZ=1.380662E-23,
     +          C=2.346E-33,
     +          A1=3.916,A2=0.074,A3=0.05)
C-----------------------------------------------------------------------
CC      RAYLGH=(PAIR/(KBOLTZ*TAIR))*5.45E-32*(WAVEL/0.55E-6)**(-4)
C-----------------------------------------------------------------------
      W=WAVEL*1.0E6
      X=A1+A2*W+A3/W
      Y=W**X
      IF(ANGEL.EQ.180.0) THEN
          RAYLGH=(PAIR/(KBOLTZ*TAIR))*C*2.0/Y
      ELSE
          PI=4.0*ATAN(1.0)
          RAYLGH=(PAIR/(KBOLTZ*TAIR))*
     +            C*(1.0+(COS(ANGEL*PI/180.0))**2)/Y
      ENDIF
C
      RETURN
      END
C-----------------------------------------------------------------------
      REAL FUNCTION RAYLGH_EXT(WAVEL,PAIR,TAIR)
C-----------------------------------------------------------------------
C
C     This function calculates the molecular Rayleigh extinction
C     coefficient.
C
C     Source: Nicolet et al.: Planet Sp. Sci. 30,935,1982.
c
C
C     Input:
C         WAVEL:   (m)            Incident wavelength
C         PAIR:    (Pa)           Ambient air pressure
C         TAIR     (K)            Temperature (K)
C
C     Output:
C         RAYLGH   (1/m)       Rayleigh EXTINCTION coefficient
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL WAVEL,TAIR,PAIR
      REAL KBOLTZ,PI,C,A1,A2,A3,X,Y,W
      PARAMETER(KBOLTZ=1.380662E-23,
     +          C=2.346E-33,
     +          A1=3.916,A2=0.074,A3=0.05)
C-----------------------------------------------------------------------
CC      RAYLGH=(PAIR/(KBOLTZ*TAIR))*5.45E-32*(WAVEL/0.55E-6)**(-4)
C-----------------------------------------------------------------------
      W=WAVEL*1.0E6
      X=A1+A2*W+A3/W
      Y=W**X
      PI=4.0*ATAN(1.0)
      RAYLGH_EXT=16.0*PI*(PAIR/(KBOLTZ*TAIR))*C*2.0/(Y*3.0)
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SINDEX(WAVEL,W,RHOSAS,INDEX)
C-----------------------------------------------------------------------
C
C     This subroutine calculates real part refractive indices of sulfuric
C     acid solution of acid weight fraction W and bulk density RHOSAS
C     at wavelength WAVEL.
C
C     Source: Palmer and Williams, App. Opt. 14,208,1975
C
C     Input:
C
C     WAVEL:  (m)         Wavelength [360nm; 1000nm]
C     W:                  Sulfuric acid weight fraction [0.0;1.0]
C     RHOSAS: (kg/m**3)   Sulfuric acid density
C
C     Output:
C     INDEX:              Refractive index (real part)
C-----------------------------------------------------------------------
C
C
      IMPLICIT NONE
      REAL WAVEL,W,RHOSAS,INDEX
      REAL WL
      REAL A,ROSAS0,N0,N1,ROSAS,T0
      INTEGER NWGT,NWAVEL,I,J
      PARAMETER(NWGT=6,NWAVEL=31)
      REAL WTAB(NWGT),WAVTAB(NWAVEL),REFTAB(NWGT,NWAVEL)
      PARAMETER(T0=300.0)
      DATA (WTAB(I),I=1,NWGT)/0.25,0.38,0.50,0.75,0.845,0.956/
      DATA (WAVTAB(J),J=1,NWAVEL)/
     +        1.111,1.087,1.064,1.042,1.020,1.000,0.980,0.962,
     +        0.943,0.926,0.909,0.893,0.877,0.862,0.847,0.833,
     +        0.820,0.806,0.794,0.781,0.769,0.758,0.746,0.735,
     +        0.725,0.714,0.702,0.556,0.449,0.408,0.360/
C
      DATA ((REFTAB(I,J),I=1,NWGT),J=1,NWAVEL)/
     +357,372,386,418,423,425,358,373,387,419,424,425,
     +358,373,387,420,425,426,358,374,388,421,425,426,
     +358,375,389,421,426,427,359,375,389,422,427,427,
     +359,376,390,422,427,427,360,377,390,423,428,427,
     +360,377,391,423,428,427,360,377,391,424,429,428,
     +361,377,391,424,429,428,361,378,392,425,430,428,
     +361,378,392,425,430,429,361,379,392,425,431,429,
     +361,379,392,426,432,430,362,380,392,426,432,430,
     +362,380,392,427,432,430,362,380,392,427,432,430,
     +362,380,393,427,433,430,362,380,393,427,433,431,
     +362,381,393,427,434,431,362,381,393,427,434,431,
     +362,381,393,427,434,431,362,381,393,427,435,431,
     +363,381,394,427,435,432,363,381,394,427,435,432,
     +363,382,394,428,436,432,366,384,397,431,438,434,
     +369,387,402,432,442,438,373,392,408,438,448,443,
     +383,407,421,452,463,459/
C
C
C
      SAVE WTAB,WAVTAB,REFTAB
C
      DO 100 I=NWGT,1,-1
          IF(WTAB(I).LT.W) GOTO 200
 100  CONTINUE
      I=1
 200  CONTINUE
C
      WL=WAVEL*1.0E6
      WL=AMAX1(WAVTAB(NWAVEL),AMIN1(WL,WAVTAB(2)))
      DO 300 J=2,NWAVEL
          IF(WAVTAB(J).LE.WL) GOTO 400
 300  CONTINUE
      J=2
 400  CONTINUE
C
      IF(I.LT.NWGT) THEN
          N0=REFTAB(I,J)+
     +     (W-WTAB(I))*(REFTAB(I+1,J)-REFTAB(I,J))/(WTAB(I+1)-WTAB(I))
      ELSE
          N0=REFTAB(NWGT-1,J)+
     +      (W-WTAB(NWGT-1))*(REFTAB(NWGT,J)-REFTAB(NWGT-1,J))/
     +      (WTAB(NWGT)-WTAB(NWGT-1))
      ENDIF
      IF(I.LT.NWGT) THEN
          N1=REFTAB(I,J-1)+
     +            (W-WTAB(I))*(REFTAB(I+1,J-1)-REFTAB(I,J-1))/
     +            (WTAB(I+1)-WTAB(I))
      ELSE
          N1=REFTAB(NWGT-1,J-1)+
     +      (W-WTAB(NWGT-1))*(REFTAB(NWGT,J-1)-REFTAB(NWGT-1,J-1))/
     +      (WTAB(NWGT)-WTAB(NWGT-1))
      ENDIF
      N0=N0+(N1-N0)*(WL-WAVTAB(J))/(WAVTAB(J-1)-WAVTAB(J))
      N0=1.0+N0/1000.0
C
      ROSAS0=ROSAS(T0,W)
C
      A=(N0*N0-1.0)/((N0*N0+2.0)*ROSAS0)
      INDEX=SQRT((2.0*A*RHOSAS+1.0)/(1.0-A*RHOSAS))
c      write(*,*) 'N0= ',n0,'  Index= ',index
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE STS_INDEX(WAVEL,TAIR,WSA,WNA,INDEX)
C-----------------------------------------------------------------------
C
C     This subroutine calculates the real part of refractive indices of
C     liquid super cooled ternary solutions (H2SO4/HNO3/H2O).
C
C     Source:   Luo et al., GRL 23,3707,1996 (method).
C     Code by   U.K. Krieger, J. Moessinger, B. Luo, U. Weers and Th. Peter
C               Applied Optics 39, 3691, 2000.
C
C
C     Input:
C
C     WAVEL:  (m)         Wavelength [360nm; 1000nm]
C     TAIR:   (K)         Ambient air temperature
C     WSA                 Sulfuric acid weight fraction [0.0;1.0]
C     WNA                 Nitric acid weight fraction [0.0;1.0]
C
C     Output:
C     INDEX:              Refractive index (real part)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL WAVEL,TAIR,WSA,WNA,INDEX
      REAL*8 XL,T,WTS,WTN,REFRA,WT,RHO,RHO1,A,AH2O,ASUL,AHNO3
      XL=DBLE(WAVEL)*1.0D6
      T=DBLE(TAIR)
      WTS=DBLE(WSA)
      WTN=DBLE(WNA)
      IF (T.LT .185.OR.T.GT .370.OR.WTS+WTN.GT.
     * 0.7D0.OR.WTS+WTN.LT.0.05D0.OR.XL.LT.0.214
     * .OR.XL.GT.2.) then
            PRINT*, ' WARNING:  OUT OF
     * RANGE, THE RESULT MAY BE INCORRECT !'
       WRITE(*,*) t,wts+wtn,xl
      endif
      WT=WTS+WTN
      RHO=RHO1(WTS,WTN,T)
      A = (1.D0-WT)*RHO/18.016D0*AH2O(XL)+
     * WTS*RHO/98.08D0*ASUL(WT,XL)+WTN*
     * RHO/63.016D0*AHNO3(WT,XL)
      REFRA=DSQRT((1 .D0+2 .D0*A)/(1-A))
      INDEX=REAL(REFRA)
      RETURN
      END

C     MOLAR REFRACTIVITY OF HNO3
      FUNCTION AHNO3(WTN,XL)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(12)
      DATA X
     * /   9.9013,  -1.4456,   1.7326,  0.38526, -0.63329,  0.20787,
     *  3.61314E-3,  0.44435, -0.32184,  0.19673, -0.63184,  1.23058E-2/

      A0=X(1)+X(2)*WTN+X(3)*WTN**2+X(10)*WTN**3
      A1=X(4)+X(5)*WTN+X(6)*WTN**2+X(11)*WTN**3
      A2=X(7)+X(8)*WTN+X(9)*WTN**2+X(12)*WTN**3
      AHNO3= A0+A1/XL+A2/XL**2
      RETURN
      END

C     MOLAR REFRACTIVITY OF H2S04
      FUNCTION ASUL(WTS,XL)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XL, X(20)
      DATA X /
     *  12.4714,  0.45272,  0.208, -0.1258,  1.8034E-002,  -4.6424,
     *   5.6781,  -2.6354,  0.43043, -1.5631E-002 ,   10.149,  -7.6642,
     *   2.9535, -0.586,   6.0278E-002,  -22.475,    71.4748,
     *  -53.5743,   16.294,  -1.73956/

      XMU=1.D0/XL
      WW=(WTS-0.5D0)
      X1=X(1)+X(6)*WW+X(11)*WW**2+X(16)*WW**3
      X2=X(2)+X(7)*WW+X(12)*WW**2+X(17)*WW**3
      X3=X(3)+X(8)*WW+X(13)*WW**2+X(18)*WW**3
      X4=X(4)+X(9)*WW+X(14)*WW**2+X(19)*WW**3
      X5=X(5)+X(10)*WW+X(15)*WW**2+X(20)*WW**3
      ASUL=X1+X2*XMU+X3*XMU**2+X4*XMU**3+X5*XMU**4
      RETURN
      END

C     MOLAR REFRACTIVITY OF H20
      FUNCTION AH2O (XL)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 B1(5) ,B2(5)
      DATA B1 /1.4211,7.540441,-9.81512155,5.7218513,-1.237611/
      DATA B2 /3.5740152,5.820734E-2,1.35319E-2,-1.5603E-3,8.85744E-4/
      X=1./XL
      IF(X.LE.1.23) THEN
         AH2O=B1(1)+B1(2)*X+B1(3)*X**2+B1(4)*X**3+B1(5)*X**4
      ELSE
         AH2O=B2(1)+B2(2)*X+B2(3)*X**2+B2(4)*X**3+B2(5)*X**4
      ENDIF
      RETURN
      END

C     DENSITY OF TERNARY SOLUTION IN G/CM3
      FUNCTION RHO1(WTS,WTN,T)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL X(22)
      DATA X /2.393284E-02,-4.359335E-05,7.961181E-08,0.0,
     * -0.198716351,1 .39564574E-03,-2.020633E-06,
     * 0.51684706,-3.0539E-03,4.505475E-06,-0.30119511,
     * 1.840408E-03,-2 .7221253742E-06,-0.113316741169,
     * 8.47763E-04,-1 .22336185E-06,0.3455282,-2.2111E-03,
     * 3.503768245E-06,-0.2315332 ,1.60074E-03,-2.5827835E-06/
      W=WTS+WTN
      WTH=1.-W
      V1=X(1)+X(2)*T+X(3)*T**2+X(4)*T**3
      VS=X(5)+X(6)*T+X(7)*T**2+(X(8)+X(9)*T+X(10)*T**2)*W
     * +(X(11)+X(12)*T+X(13)*T**2)*W*W
      VN=X(14)+X(15)*T+X(16)*T**2+(X(17)+X(18)*T+X(19)*T**2) *W
     * +(X(20)+X(21)*T+X(22)*T**2)*W*W
      VMCAL=WTH/18.016D0*V1+VS*WTS/98.08D0+VN*WTN/63.016D0
      RHO1=1/VMCAL/1000 .D0
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LGNDIS(NDTOT,MEDIAN,GSTDEV,NBINS,RADIUS,ND)
C-----------------------------------------------------------------------
C
C     This subroutine calculates a lognormal size distribution,
C     characterised by three parameters NDTOT, MEDIAN, and GSTDEV,
C     discretisezed into a number NBINS on a lineary radius scale.
C
C
C     Input:
C         NDTOT:   (par/m**3)     Total number density of particles
C         MEDIAN:  (m)            Median particle radius
C         GSTDEV:                 Geometric standard deveation
C         NBINS:                  Number of radius bins on a lineary
C                                 radius scale
C         RADIUS:  (m)            Radius values for each bin
C     Output:
C         ND:      (par/m**3 m)   Differential number density of particles
C
C-----------------------------------------------------------------------
      INTEGER NBINS
      REAL NDTOT,MEDIAN,GSTDEV,
     +     RADIUS(NBINS),ND(NBINS)
C
      REAL LNGST,LNGST2,PI,STWOPI
      INTEGER I
C
      PI=4.0*ATAN(1.0)
      STWOPI=SQRT(2.0*PI)
      LNGST=ALOG(GSTDEV)
      LNGST2=LNGST*LNGST
      DO 100 I=1,NBINS
          ND(I)=(NDTOT/(STWOPI*RADIUS(I)*LNGST))*
     +           EXP(-((ALOG(RADIUS(I)/MEDIAN))**2)/(2.0*LNGST2))
 100  CONTINUE
C
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE BKSDST(NBINS,DR,ND,SBACK,BKS)
C----------------------------------------------------------------------------
C
C     This subroutine calculates the particle backscatter from a 
C     distribution of spherical particles.
C
C     Input:
C         NBINS:                  Number of radius bins on a radius scale
C         DR:      (m)            Radius increment
C         ND:      (par/m**3 m)   Differential number density of particles
C         SBACK:   (m**2/sr)      Differential backscattering 
C                                 cross sections
C
C     Output:
C         BKS:     (1/m sr)       Particle backscatter
C
C
C-----------------------------------------------------------------------
      INTEGER NBINS
      REAL DR,ND(NBINS),SBACK(NBINS),BKS
      REAL BPART
      INTEGER I
C
C-----------------------------------------------------------------------
C     Particle backscattering:
      BPART=0.0
      DO 100 I=1,NBINS
          BPART=BPART+SBACK(I)*ND(I)*DR
 100  CONTINUE
      BKS=BPART
      RETURN
      END
C
*******************************************************************************
*     SUBROUTINE SAGE_NUMBER_DENSITY(PAIR,LATITUDE,MEDIAN,GSTDEV,ND)
*******************************************************************************
*     Subroutine used to calculate aerosol number densities consistent
*     with SAGE II averaged extinction measurements, assuming constant
*     lognormal size distribution median radius and geometric standard
*     deviation.
*     Data source: Hitchman et al.,  JGR 99, 20689, 1994.
*
*     Input:
*               PAIR:    (Pa)        Air pressure
*               LATITUDE (deg)       Latitude [-90,+90]
*               MEDIAN   (m)         Lognormal median radius
*               GSTDEV               Lognormal geometric standar deviation
*     Output:
*               ND       (m**-3)     Aerosol number density
*
      SUBROUTINE SAGE_NUMBER_DENSITY(PAIR,LATITUDE,MEDIAN,GSTDEV,ND)
      REAL PAIR,LATITUDE,MEDIAN,GSTDEV,ND
      REAL FACTOR,AREA
      FACTOR=4.0*4.0*ATAN(1.0)*(MEDIAN**2)*EXP(2.0*(ALOG(GSTDEV))**2)
      CALL SAGE_AREA(PAIR,LATITUDE,AREA)
      ND=AREA/FACTOR
      RETURN
      END
*
*******************************************************************************
*     SUBROUTINE SAGE_AREA(PAIR,LATITUDE,AREA)
*******************************************************************************
*     Subroutine used to calculate aerosol surface area densities consistent
*     with SAGE II averaged extinction measurements.
*     Data source: Hitchman et al.,  JGR 99, 20689, 1994.
*
*     Input:
*               PAIR:    (Pa)        Air pressure
*               LATITUDE (deg)       Latitude [-90,+90]
*     Output:
*               AREA:    (m**2/m**3) Sulfate aerosol surface area density
*
      SUBROUTINE SAGE_AREA(PAIR,LATITUDE,AREA)
      REAL PAIR,LATITUDE,AREA
      INTEGER ILAT,IALT,A,L
      REAL H0,P0,P,H
      REAL S,Z,ZT,E0
      REAL HEIGHT
      PARAMETER(H0=6.5,P0=1000.0E2)
      PARAMETER(ILAT=19,IALT=14)
      REAL ALT(IALT),ZTROP(ILAT),LAT(ILAT),E(ILAT,IALT)
      DATA LAT/
     +-90.,-80.,-70.,-60.,-50.,-40.,-30.,-20.,-10.,
     + 0.,10.,20.,30.,40.,50.,60.,70.,80.,90./
      DATA ALT/8.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28.,30.,32.,34./
      DATA ZTROP/
     +8.,8.,8.,8.,8.,8.,10.,16.,16.,
     +16.,16.,14.,10.,8.,8.,8.,8.,8.,8./
      DATA E/
     +3.9,3.9,5.0,5.1,5.0,6.7,3.4,2.1,2.1,
     +2.7,2.5,1.7,4.0,6.0,5.8,6.6,5.9,6.1,6.1,
     +5.3,5.3,5.6,5.7,5.3,5.1,3.4,2.1,2.1,
     +2.7,2.5,1.7,4.0,7.1,7.7,8.6,7.9,7.8,7.8,
     +5.6,5.6,5.3,4.6,4.5,3.8,2.6,2.1,2.1,
     +2.7,2.5,1.7,3.0,5.0,6.0,5.7,7.1,6.6,6.6,
     +4.8,4.8,4.5,3.9,4.2,3.3,2.5,2.1,2.1,
     +2.7,2.5,1.7,2.8,4.2,5.3,4.6,5.8,5.4,5.4,
     +2.1,2.1,2.8,3.1,4.2,3.8,2.8,2.1,2.1,
     +2.7,2.5,2.3,3.0,4.4,4.8,4.1,4.4,4.1,4.1,
     +1.0,1.0,1.7,2.0,3.3,3.7,3.9,3.3,3.2,
     +2.9,2.8,3.2,4.1,4.1,3.6,3.1,2.6,2.4,2.4,
     +.62,.62,.91,1.1,2.0,2.5,3.1,4.0,5.5,
     +8.1,6.7,4.4,3.2,2.6,2.1,1.7,1.6,1.2,1.2,
     +.31,.31,.47,.56,1.0,1.3,1.7,2.6,4.3,
     +5.5,4.7,3.0,1.8,1.3,1.0,.76,.70,.54,.54,
     +.15,.15,.25,.28,.48,.66,.87,1.5,2.6,
     +3.1,2.5,1.5,.88,.60,.47,.33,.34,.28,.28,
     +.08,.08,.13,.14,.20,.29,.42,.74,1.3,
     +1.5,1.3,.80,.42,.26,.21,.13,.16,.16,.16,
     +.05,.05,.08,.07,.08,.12,.18,.34,.56,
     +.66,.57,.37,.18,.11,.09,.06,.08,.10,.10,
     +.05,.05,.05,.04,.04,.05,.08,.14,.22,
     +.24,.20,.14,.08,.05,.05,.03,.05,.06,.06,
     +.01,.01,.02,.02,.03,.03,.04,.06,.09,
     +.09,.07,.06,.04,.03,.03,.02,.02,.02,.02,
     +.01,.01,.01,.02,.02,.02,.02,.03,.04,
     +.04,.03,.03,.02,.02,.02,.02,.01,.01,.01/
c
      REAL AREA_LAT(ILAT,IALT),Y2AREA_LAT(ILAT,IALT),YLATWORK(ILAT)
      REAL AREA_ALT(IALT),Y2AREA_ALT(IALT),YALTWORK(IALT)
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      SAVE FIRST,AREA_LAT,Y2AREA_LAT,ALT,LAT
C
      H(P)=H0*ALOG(P0/P)
      S(Z,ZT,E0)=((0.8*(Z-ZT)/20.0)+0.7)*E0*1.0E-6
C
      IF(FIRST) THEN
          FIRST=.FALSE.
          DO A=1,IALT
              DO L=1,ILAT
                  AREA_LAT(L,A)=S(ALT(A),ZTROP(L),E(L,A))
              END DO
              CALL SPLINE(LAT,AREA_LAT(1,A),ILAT,
     +                    YLATWORK,Y2AREA_LAT(1,A))
          END DO
      ENDIF
C
      HEIGHT=AMAX1(ALT(1),AMIN1(H(PAIR),ALT(IALT)))
      DO A=1,IALT
          CALL SPLINT(LAT,AREA_LAT(1,A),Y2AREA_LAT(1,A),ILAT,
     +                LATITUDE,AREA_ALT(A))
      END DO
      CALL SPLINE(ALT,AREA_ALT,IALT,YALTWORK,Y2AREA_ALT)
      CALL SPLINT(ALT,AREA_ALT,Y2AREA_ALT,IALT,HEIGHT,AREA)
      RETURN
      END SUBROUTINE

