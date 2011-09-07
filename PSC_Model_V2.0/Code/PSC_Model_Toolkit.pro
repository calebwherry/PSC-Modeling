; Caleb Wherry
; calebwherry@gmail.com
;
; PSC_Model_Toolkit.pro
; IDL GUI Interface for PSC Model
; Created: October 8th, 2008
;
; Last Modified: October 27th, 2008
;


; Main procedure.
; All widgets intialized here.
PRO PSC_Model_Toolkit, GROUP_LEADER=wGroup

  ; Base Window
  baseWindow = Widget_Base( GROUP_LEADER=wGroup, UNAME='baseWindow'  $
      ,XOFFSET=250 ,YOFFSET=200 ,SCR_XSIZE=800 ,SCR_YSIZE=600 ,TITLE='PSC'+ $
      ' Model Toolkit' ,SPACE=3 ,XPAD=3 ,YPAD=3, TLB_FRAME_ATTR=1)

  ; Tab Widget
  tabs = Widget_Tab(baseWindow, UNAME='tabs' ,YOFFSET=4  $
      ,SCR_XSIZE=800 ,SCR_YSIZE=600)


  ; Tab 1: Trajectory Model Widget Intializations
  ;-------------------------------------------------------------------------------
  ;-------------------------------------------------------------------------------
  ; Named "tab2Base" because of intialization purposes.
  tab2Base = Widget_Base(tabs, UNAME='tab2Base' ,SCR_XSIZE=793  $
      ,SCR_YSIZE=539 ,TITLE='Trajectory Model' ,SPACE=3 ,XPAD=3  $
      ,YPAD=3)


  dateL = Widget_Label(tab2Base, UNAME='dateL'  $
      ,XOFFSET=205 ,YOFFSET=20 ,SCR_XSIZE=120 ,SCR_YSIZE=20  $
      ,/ALIGN_LEFT ,VALUE='-> Date (YYYYMMDD) :')


  dateTB = Widget_Text(tab2Base, UNAME='dateTB'  $
      ,XOFFSET=325 ,YOFFSET=17 ,SCR_XSIZE=120 ,SCR_YSIZE=20 ,/EDITABLE  $
      ,XSIZE=20 ,YSIZE=1)


  createImagesB = Widget_Button(tab2Base, UNAME='createImagesB', FRAME=1  $
      ,XOFFSET=460, YOFFSET=17, SCR_XSIZE=90, SCR_YSIZE=20  $
      ,/ALIGN_CENTER, VALUE='Create Images')


  previousImageB = Widget_Button(tab2Base, UNAME='previousImageB', FRAME=1  $
      ,XOFFSET=20 ,YOFFSET=125 ,VALUE='.\Images\left.bmp', /BITMAP)


  imageDrawbox = Widget_Draw(tab2Base, UNAME='imageDrawbox', XOFFSET=60  $
      ,YOFFSET=50, SCR_XSIZE=668, SCR_YSIZE=200, /BUTTON_EVENTS)


  clickXinvisL = Widget_Label(baseWindow, UNAME='clickXinvisL'  $
    ,XOFFSET=50 ,YOFFSET=50 ,/ALIGN_LEFT ,VALUE='', SCR_XSIZE=0, SCR_YSIZE=0)


  clickYinvisL = Widget_Label(baseWindow, UNAME='clickYinvisL'  $
    ,XOFFSET=50 ,YOFFSET=50 ,/ALIGN_LEFT ,VALUE='', SCR_XSIZE=0, SCR_YSIZE=0)


  loadingL = Widget_Label(tab2Base, UNAME='loadingL'  $
      ,XOFFSET=650 ,YOFFSET=30 ,SCR_XSIZE=100 ,SCR_YSIZE=20  $
      ,/ALIGN_LEFT ,VALUE='')


  nextImageB = Widget_Button(tab2Base, UNAME='nextImageB', FRAME=1  $
      ,XOFFSET=735, YOFFSET=125 ,VALUE='.\Images\right.bmp', /BITMAP)


  whichImageInvisL = Widget_Label(baseWindow, UNAME='whichImageInvisL'  $
    ,XOFFSET=50 ,YOFFSET=50 ,/ALIGN_LEFT ,VALUE='0', SCR_XSIZE=0, SCR_YSIZE=0)


  pointsL = Widget_Label(tab2Base, UNAME='pointsL'  $
      ,XOFFSET=310, YOFFSET=280, SCR_XSIZE=80, SCR_YSIZE=20  $
      , VALUE='Points Selected')


  clearPointsB = Widget_Button(tab2Base, UNAME='clearPointsB', FRAME=1  $
      ,XOFFSET=395, YOFFSET=278, SCR_XSIZE=80, SCR_YSIZE=17  $
      ,/ALIGN_CENTER, VALUE='Clear Points')


  pointsTB = Widget_Text(tab2Base, UNAME='pointsTB'  $
      ,XOFFSET=20 ,YOFFSET=300 ,SCR_XSIZE=745 ,SCR_YSIZE=170   $
      ,/SCROLL)


  trajTimeStepL = Widget_Label(tab2Base, UNAME='trajTimeStepL'  $
      ,XOFFSET=580 ,YOFFSET=486 ,SCR_XSIZE=100 ,SCR_YSIZE=20  $
      ,/ALIGN_LEFT ,VALUE='-> Time Step (min) :')


  trajTimeStepTB = Widget_Text(tab2Base, UNAME='trajTimeStepTB'  $
      ,XOFFSET=680 ,YOFFSET=483 ,SCR_XSIZE=40 ,SCR_YSIZE=20 ,/EDITABLE  $
      ,XSIZE=20 ,YSIZE=1, VALUE=['-15'])


  trajNumofDaysL = Widget_Label(tab2Base, UNAME='trajNumofDaysL'  $
      ,XOFFSET=542 ,YOFFSET=510 ,SCR_XSIZE=130 ,SCR_YSIZE=20  $
      ,/ALIGN_LEFT ,VALUE='-> Number of Days to Run :')


  trajNumofDaysTB = Widget_Text(tab2Base, UNAME='trajNumofDaysTB'  $
      ,XOFFSET=680 ,YOFFSET=507 ,SCR_XSIZE=40 ,SCR_YSIZE=20 ,/EDITABLE  $
      ,XSIZE=20 ,YSIZE=1, VALUE=['5'])


  ; Trajectory Direction base.
  trajDirectionBase = Widget_Base(tab2Base, UNAME='trajDirectionBase'  $
      ,XOFFSET=100 ,YOFFSET=480 ,SCR_XSIZE=100 ,SCR_YSIZE=50  $
      ,COLUMN=1 ,/EXCLUSIVE)


  directionRB1 = Widget_Button(trajDirectionBase, UNAME='directionRB1'  $
      ,/ALIGN_LEFT ,VALUE='Backwards (-)')


  ; Makes directionRB1 selected.
  WIDGET_CONTROL, directionRB1, /SET_BUTTON


  directionRB2 = Widget_Button(trajDirectionBase, UNAME='directionRB2'  $
      ,/ALIGN_LEFT ,VALUE='Forwards (+)')


  generate1B = Widget_Button(tab2Base, UNAME='generate1B', FRAME=1  $
      ,XOFFSET=290, YOFFSET=500, SCR_XSIZE=200, SCR_YSIZE=26  $
      ,/ALIGN_CENTER, VALUE='Generate *.config and *.inp Files')

  ; End of Trajectory Model Widget Intializations
  ;-------------------------------------------------------------------------------
  ;-------------------------------------------------------------------------------



  ; Tab 2: Microphysical Model Widget intializations
  ;-------------------------------------------------------------------------------
  ;-------------------------------------------------------------------------------
  ; Named "tab1Base" because of intialization purposes.
  tab1Base = Widget_Base(tabs, UNAME='tab1Base' ,SCR_XSIZE=793  $
      ,SCR_YSIZE=539 ,TITLE='Microphysical Model' ,SPACE=3 ,XPAD=3  $
      ,YPAD=3)


  generateB = Widget_Button(tab1Base, UNAME='generateB' ,FRAME=1  $
      ,XOFFSET=500 ,YOFFSET=465 ,SCR_XSIZE=100 ,SCR_YSIZE=26  $
      ,/ALIGN_CENTER ,VALUE='Generate *.inp File')


  openInputB = Widget_Button(tab1Base, UNAME='openInputB' ,FRAME=1  $
      ,XOFFSET=610 ,YOFFSET=465 ,SCR_XSIZE=125 ,SCR_YSIZE=26  $
      ,/ALIGN_CENTER ,VALUE='Open Existing *.inp File')


  runB = Widget_Button(tab1Base, UNAME='runB' ,FRAME=1 ,XOFFSET=550  $
      ,YOFFSET=500 ,SCR_XSIZE=134 ,SCR_YSIZE=26 ,/ALIGN_CENTER  $
      ,VALUE='Run Microphysical Model')


  inpFileNameL = Widget_Label(tab1Base, UNAME='inpFileNameL'  $
      ,XOFFSET=190 ,YOFFSET=12 ,SCR_XSIZE=91 ,SCR_YSIZE=18  $
      ,/ALIGN_LEFT ,VALUE='-> INP File Name :')


  inpFileNameTB = Widget_Text(tab1Base, UNAME='inpFileNameTB'  $
      ,XOFFSET=290 ,YOFFSET=9 ,SCR_XSIZE=120 ,SCR_YSIZE=20 ,/EDITABLE  $
      ,XSIZE=20 ,YSIZE=1, VALUE=['example.inp'])


  layersL = Widget_Label(tab1Base, UNAME='layersL' ,XOFFSET=226  $
      ,YOFFSET=36 ,/ALIGN_LEFT ,VALUE='-> Layers :')


  layersTB = Widget_Text(tab1Base, UNAME='layersTB' ,XOFFSET=290  $
      ,YOFFSET=35 ,SCR_XSIZE=40 ,SCR_YSIZE=20 ,/EDITABLE ,VALUE=[  $
      '-4' ] ,XSIZE=20 ,YSIZE=1)


  integrationStopTimeL = Widget_Label(tab1Base,  $
      UNAME='integrationStopTimeL' ,XOFFSET=155 ,YOFFSET=63  $
      ,/ALIGN_LEFT ,VALUE='-> Integration Stop Time :')


  integrationStopTimeTB = Widget_Text(tab1Base,  $
      UNAME='integrationStopTimeTB' ,XOFFSET=290 ,YOFFSET=60  $
      ,SCR_XSIZE=80 ,SCR_YSIZE=20 ,/EDITABLE ,VALUE=[ '150.00' ]  $
      ,XSIZE=20 ,YSIZE=1)


  timeIntervalL = Widget_Label(tab1Base, UNAME='timeIntervalL'  $
      ,XOFFSET=67 ,YOFFSET=92 ,/ALIGN_LEFT ,VALUE='-> Time Interval'+ $
      ' Between Plot/Print Output :')


  timeIntervalTB1 = Widget_Text(tab1Base, UNAME='timeIntervalTB1'  $
      ,XOFFSET=290 ,YOFFSET=91 ,SCR_XSIZE=40 ,SCR_YSIZE=20 ,/EDITABLE  $
      ,VALUE=[ '0.10' ] ,XSIZE=20 ,YSIZE=1)


  timeIntervalTB2 = Widget_Text(tab1Base, UNAME='timeIntervalTB2'  $
      ,XOFFSET=330 ,YOFFSET=91 ,SCR_XSIZE=40 ,SCR_YSIZE=20 ,/EDITABLE  $
      ,VALUE=[ '0.00' ] ,XSIZE=20 ,YSIZE=1)


  timeIntervalTB3 = Widget_Text(tab1Base, UNAME='timeIntervalTB3'  $
      ,XOFFSET=370 ,YOFFSET=91 ,SCR_XSIZE=40 ,SCR_YSIZE=20 ,/EDITABLE  $
      ,VALUE=[ '1.00' ] ,XSIZE=20 ,YSIZE=1)


  bottomLayerPotenTempL = Widget_Label(tab1Base,  $
      UNAME='bottomLayerPotenTempL' ,XOFFSET=70 ,YOFFSET=121  $
      ,SCR_XSIZE=208 ,SCR_YSIZE=15 ,/ALIGN_LEFT ,VALUE='-> Bottom'+ $
      ' Layer Potential Temperature (K) :')


  bottomLayerPotenTempTB = Widget_Text(tab1Base,  $
      UNAME='bottomLayerPotenTempTB' ,XOFFSET=290 ,YOFFSET=119  $
      ,SCR_XSIZE=80 ,SCR_YSIZE=20 ,/EDITABLE ,VALUE=[ '480.00' ]  $
      ,XSIZE=20 ,YSIZE=1)


  potenTempIncrL = Widget_Label(tab1Base, UNAME='potenTempIncrL'  $
      ,XOFFSET=85 ,YOFFSET=151 ,/ALIGN_LEFT ,VALUE='-> Potential'+ $
      ' Temperature Increment (K) :')


  potenTempIncrTB = Widget_Text(tab1Base, UNAME='potenTempIncrTB'  $
      ,XOFFSET=290 ,YOFFSET=148 ,SCR_XSIZE=80 ,SCR_YSIZE=20  $
      ,/EDITABLE ,/ALL_EVENTS ,VALUE=[ '40.00' ] ,XSIZE=20 ,YSIZE=1)


  tempCalcL = Widget_Label(tab1Base, UNAME='tempCalcL' ,XOFFSET=139  $
      ,YOFFSET=177 ,/ALIGN_LEFT ,VALUE='-> Temperature Calculation'+ $
      ' :')


  tempCalcBase = Widget_Base(tab1Base, UNAME='tempCalcBase'  $
      ,XOFFSET=285 ,YOFFSET=175 ,SCR_XSIZE=174 ,SCR_YSIZE=70  $
      ,TITLE='IDL' ,COLUMN=1 ,/EXCLUSIVE)


  tempCalcRB1 = Widget_Button(tempCalcBase, UNAME='tempCalcRB1'  $
      ,/ALIGN_LEFT ,VALUE='0: Ramp')

  ; Makes tempCalcRB1 selected.
  WIDGET_CONTROL, tempCalcRB1, /SET_BUTTON


  tempCalcRB2 = Widget_Button(tempCalcBase, UNAME='tempCalcRB2'  $
      ,/ALIGN_LEFT ,VALUE='1: Temperature Table')


  tempCalcRB3 = Widget_Button(tempCalcBase, UNAME='tempCalcRB3'  $
      ,/ALIGN_LEFT ,VALUE='2: Temperature/Pressure Table')

  ; Used to tell which radio button is pressed.
  tempCalcInvisL = Widget_Label(baseWindow, UNAME='tempCalcInvisL'  $
    ,XOFFSET=50 ,YOFFSET=50 ,/ALIGN_LEFT ,VALUE='0', SCR_XSIZE=0, SCR_YSIZE=0)


  tempPressHisFileL = Widget_Label(tab1Base,  $
      UNAME='tempPressHisFileL' ,XOFFSET=95 ,YOFFSET=261 ,/ALIGN_LEFT  $
      ,VALUE='-> Temperature/Pressure History File :')


  tempPressHisFileTB = Widget_Text(tab1Base,  $
      UNAME='tempPressHisFileTB' ,XOFFSET=290 ,YOFFSET=258  $
      ,SCR_XSIZE=120 ,SCR_YSIZE=20 ,XSIZE=20 ,YSIZE=1, VALUE='---')


  rampTime1L = Widget_Label(tab1Base, UNAME='rampTime1L' ,XOFFSET=159  $
      ,YOFFSET=291 ,/ALIGN_LEFT ,VALUE='-> Ramp Time/Temp 1 :')


  rampTime1TB1 = Widget_Text(tab1Base, UNAME='rampTime1TB1'  $
      ,XOFFSET=290 ,YOFFSET=288 ,SCR_XSIZE=60 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '0.00' ] ,XSIZE=20 ,YSIZE=1)


  rampTime1TB2 = Widget_Text(tab1Base, UNAME='rampTime1TB2'  $
      ,XOFFSET=350 ,YOFFSET=288 ,SCR_XSIZE=60 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '200.00' ] ,XSIZE=20 ,YSIZE=1)


  rampTime2L = Widget_Label(tab1Base, UNAME='rampTime2L' ,XOFFSET=159  $
      ,YOFFSET=319 ,/ALIGN_LEFT ,VALUE='-> Ramp Time/Temp 2 :')


  rampTime2TB1 = Widget_Text(tab1Base, UNAME='rampTime2TB1'  $
      ,XOFFSET=290 ,YOFFSET=315 ,SCR_XSIZE=60 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '50.00' ] ,XSIZE=20 ,YSIZE=1)


  rampTime2TB2 = Widget_Text(tab1Base, UNAME='rampTime2TB2'  $
      ,XOFFSET=350 ,YOFFSET=315 ,SCR_XSIZE=60 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '200.00' ] ,XSIZE=20 ,YSIZE=1)


  rampTime3L = Widget_Label(tab1Base, UNAME='rampTime3L' ,XOFFSET=159  $
      ,YOFFSET=346 ,/ALIGN_LEFT ,VALUE='-> Ramp Time/Temp 3 :')


  rampTime3TB1 = Widget_Text(tab1Base, UNAME='rampTime3TB1'  $
      ,XOFFSET=290 ,YOFFSET=343 ,SCR_XSIZE=60 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '100.00' ] ,XSIZE=20 ,YSIZE=1)


  rampTime3TB2 = Widget_Text(tab1Base, UNAME='rampTime3TB2'  $
      ,XOFFSET=350 ,YOFFSET=343 ,SCR_XSIZE=60 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '195.00' ] ,XSIZE=20 ,YSIZE=1)


  rampTime4L = Widget_Label(tab1Base, UNAME='rampTime4L' ,XOFFSET=159  $
      ,YOFFSET=373 ,/ALIGN_LEFT ,VALUE='-> Ramp Time/Temp 4 :')


  rampTime4TB1 = Widget_Text(tab1Base, UNAME='rampTime4TB1'  $
      ,XOFFSET=290 ,YOFFSET=369 ,SCR_XSIZE=60 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '150.00' ] ,XSIZE=20 ,YSIZE=1)


  rampTime4TB2 = Widget_Text(tab1Base, UNAME='rampTime4TB2'  $
      ,XOFFSET=350 ,YOFFSET=369 ,SCR_XSIZE=60 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '195.00' ] ,XSIZE=20 ,YSIZE=1)


  tempSineOsciL = Widget_Label(tab1Base, UNAME='tempSineOsciL'  $
      ,XOFFSET=68 ,YOFFSET=399 ,/ALIGN_LEFT ,VALUE='-> Temp Sine'+ $
      ' Oscillation Period/Amplitude :')


  tempSineOsciTB1 = Widget_Text(tab1Base, UNAME='tempSineOsciTB1'  $
      ,XOFFSET=290 ,YOFFSET=395 ,SCR_XSIZE=60 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '0.00' ] ,XSIZE=20 ,YSIZE=1)


  tempSineOsciTB2 = Widget_Text(tab1Base, UNAME='tempSineOsciTB2'  $
      ,XOFFSET=350 ,YOFFSET=395 ,SCR_XSIZE=60 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '0.00' ] ,XSIZE=20 ,YSIZE=1)


  mixH20RatioL = Widget_Label(tab1Base, UNAME='mixH20RatioL'  $
      ,XOFFSET=121 ,YOFFSET=425 ,/ALIGN_LEFT ,VALUE='-> Mixing Ratio'+ $
      ' of Water Vapor :')


  mixH20RatioTB1 = Widget_Text(tab1Base, UNAME='mixH20RatioTB1'  $
      ,XOFFSET=290 ,YOFFSET=421 ,SCR_XSIZE=30 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '5.00' ] ,XSIZE=20 ,YSIZE=1)


  mixH20RatioTB2 = Widget_Text(tab1Base, UNAME='mixH20RatioTB2'  $
      ,XOFFSET=320 ,YOFFSET=421 ,SCR_XSIZE=30 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '5.00' ] ,XSIZE=20 ,YSIZE=1)


  mixH20RatioTB3 = Widget_Text(tab1Base, UNAME='mixH20RatioTB3'  $
      ,XOFFSET=350 ,YOFFSET=421 ,SCR_XSIZE=30 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '5.00' ] ,XSIZE=20 ,YSIZE=1)


  mixH20RatioTB4 = Widget_Text(tab1Base, UNAME='mixH20RatioTB4'  $
      ,XOFFSET=380 ,YOFFSET=421 ,SCR_XSIZE=30 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '5.00' ] ,XSIZE=20 ,YSIZE=1)


  mixNARatioL = Widget_Label(tab1Base, UNAME='mixNARatioL'  $
      ,XOFFSET=132 ,YOFFSET=450 ,/ALIGN_LEFT ,VALUE='-> Mixing Ratio'+ $
      ' of Nitric Acid :')

  mixNARatioTB1 = Widget_Text(tab1Base, UNAME='mixNARatioTB1'  $
      ,XOFFSET=290 ,YOFFSET=448 ,SCR_XSIZE=30 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '16.00' ] ,XSIZE=20 ,YSIZE=1)

  mixNARatioTB2 = Widget_Text(tab1Base, UNAME='mixNARatioTB2'  $
      ,XOFFSET=320 ,YOFFSET=448 ,SCR_XSIZE=30 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '16.00' ] ,XSIZE=20 ,YSIZE=1)


  mixNARatioTB3 = Widget_Text(tab1Base, UNAME='mixNARatioTB3'  $
      ,XOFFSET=350 ,YOFFSET=448 ,SCR_XSIZE=30 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '16.00' ] ,XSIZE=20 ,YSIZE=1)


  mixNARatioTB4 = Widget_Text(tab1Base, UNAME='mixNARatioTB4'  $
      ,XOFFSET=380 ,YOFFSET=448 ,SCR_XSIZE=30 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '16.00' ] ,XSIZE=20 ,YSIZE=1)


  latitudeL = Widget_Label(tab1Base, UNAME='latitudeL' ,XOFFSET=215  $
      ,YOFFSET=477 ,SCR_XSIZE=57 ,SCR_YSIZE=18 ,/ALIGN_LEFT  $
      ,VALUE='-> Latitude :')


  latitudeTB = Widget_Text(tab1Base, UNAME='latitudeTB' ,XOFFSET=290  $
      ,YOFFSET=474 ,SCR_XSIZE=60 ,SCR_YSIZE=20 ,/EDITABLE ,VALUE=[  $
      '80.00' ] ,XSIZE=20 ,YSIZE=1)


  homogenNATNucFacL = Widget_Label(tab1Base,  $
      UNAME='homogenNATNucFacL' ,XOFFSET=20 ,YOFFSET=507 ,/ALIGN_LEFT  $
      ,VALUE='-> Homogeneous NAT Nucleation Correction Factor :')


  homogeNATFactoSB = Widget_Slider(tab1Base, UNAME='homogeNATFactoSB'  $
      ,XOFFSET=290 ,YOFFSET=505 ,SCR_XSIZE=120 ,SCR_YSIZE=18  $
      ,/SUPPRESS_VALUE ,MAXIMUM=10000 ,VALUE=100)


  homogeNATFactoL = Widget_Label(tab1Base, UNAME='homogeNATFactoL'  $
      ,XOFFSET=425 ,YOFFSET=508 ,SCR_XSIZE=28,SCR_YSIZE=20  $
      ,/ALIGN_LEFT ,VALUE='1.000')


  minRadiusL = Widget_Label(tab1Base, UNAME='minRadiusL' ,XOFFSET=560  $
      ,YOFFSET=13 ,/ALIGN_LEFT ,VALUE='-> Min Radius (Microns) :')


  minRadiusTB = Widget_Text(tab1Base, UNAME='minRadiusTB'  $
      ,XOFFSET=690 ,YOFFSET=10 ,SCR_XSIZE=80 ,SCR_YSIZE=21 ,/EDITABLE  $
      ,VALUE=[ '0.001' ] ,XSIZE=20 ,YSIZE=1)


  maxRadiusL = Widget_Label(tab1Base, UNAME='maxRadiusL' ,XOFFSET=557  $
      ,YOFFSET=39 ,/ALIGN_LEFT ,VALUE='-> Max Radius (Microns) :')


  maxRadiusTB = Widget_Text(tab1Base, UNAME='maxRadiusTB'  $
      ,XOFFSET=690 ,YOFFSET=36 ,SCR_XSIZE=80 ,SCR_YSIZE=21 ,/EDITABLE  $
      ,VALUE=[ '100.000' ] ,XSIZE=20 ,YSIZE=1)


  maxIntegTimeStepL = Widget_Label(tab1Base,  $
      UNAME='maxIntegTimeStepL' ,XOFFSET=484 ,YOFFSET=67 ,/ALIGN_LEFT  $
      ,VALUE='-> Max Integration Time Step (Seconds) :')


  maxIntegTimeStepTB = Widget_Text(tab1Base,  $
      UNAME='maxIntegTimeStepTB' ,XOFFSET=690 ,YOFFSET=64  $
      ,SCR_XSIZE=80 ,SCR_YSIZE=21 ,/EDITABLE ,VALUE=[ '360.00' ]  $
      ,XSIZE=20 ,YSIZE=1)


  timeUnitsL = Widget_Label(tab1Base, UNAME='timeUnitsL' ,XOFFSET=569  $
      ,YOFFSET=94 ,/ALIGN_LEFT ,VALUE='-> Time Units (h/d/m) :')


  timeUnitsTB = Widget_Text(tab1Base, UNAME='timeUnitsTB'  $
      ,XOFFSET=690 ,YOFFSET=91 ,SCR_XSIZE=80 ,SCR_YSIZE=21 ,/EDITABLE  $
      ,VALUE=[ 'h' ] ,XSIZE=20 ,YSIZE=1)


  integStartTimeL = Widget_Label(tab1Base, UNAME='integStartTimeL'  $
      ,XOFFSET=556 ,YOFFSET=122 ,/ALIGN_LEFT ,VALUE='-> Integration'+ $
      ' Start Time :')


  integStartTimeTB = Widget_Text(tab1Base, UNAME='integStartTimeTB'  $
      ,XOFFSET=690 ,YOFFSET=119 ,SCR_XSIZE=80 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '0.00' ] ,XSIZE=20 ,YSIZE=1)


  numDensSulfAeroL = Widget_Label(tab1Base, UNAME='numDensSulfAeroL'  $
      ,XOFFSET=505 ,YOFFSET=148 ,/ALIGN_LEFT ,VALUE='-> Number'+ $
      ' Density of Sulf Aerosols :')


  numDensSulfAeroTB1 = Widget_Text(tab1Base,  $
      UNAME='numDensSulfAeroTB1' ,XOFFSET=690 ,YOFFSET=144  $
      ,SCR_XSIZE=40 ,SCR_YSIZE=21 ,/EDITABLE ,VALUE=[ '10.00' ]  $
      ,XSIZE=20 ,YSIZE=1)


  numDensSulfAeroTB2 = Widget_Text(tab1Base,  $
      UNAME='numDensSulfAeroTB2' ,XOFFSET=730 ,YOFFSET=144  $
      ,SCR_XSIZE=40 ,SCR_YSIZE=21 ,/EDITABLE ,VALUE=[ '0.00' ]  $
      ,XSIZE=20 ,YSIZE=1)


  medianRadiusL = Widget_Label(tab1Base, UNAME='medianRadiusL'  $
      ,XOFFSET=539 ,YOFFSET=173 ,SCR_XSIZE=137 ,SCR_YSIZE=18  $
      ,/ALIGN_LEFT ,VALUE='-> Median Radius (Microns) :')


  medianRadiusTB1 = Widget_Text(tab1Base, UNAME='medianRadiusTB1'  $
      ,XOFFSET=690 ,YOFFSET=169 ,SCR_XSIZE=40 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '0.032' ] ,XSIZE=20 ,YSIZE=1)


  medianRadiusTB2 = Widget_Text(tab1Base, UNAME='medianRadiusTB2'  $
      ,XOFFSET=730 ,YOFFSET=169 ,SCR_XSIZE=40 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '0.150', '' ] ,XSIZE=20 ,YSIZE=1)


  geometricSTDDevL = Widget_Label(tab1Base, UNAME='geometricSTDDevL'  $
      ,XOFFSET=560 ,YOFFSET=198 ,/ALIGN_LEFT ,VALUE='-> Geometric'+ $
      ' Std. Dev. :')


  geometricSTDDevTB1 = Widget_Text(tab1Base,  $
      UNAME='geometricSTDDevTB1' ,XOFFSET=690 ,YOFFSET=194  $
      ,SCR_XSIZE=40 ,SCR_YSIZE=21 ,/EDITABLE ,VALUE=[ '1.60' ]  $
      ,XSIZE=20 ,YSIZE=1)


  geometricSTDDevTB2 = Widget_Text(tab1Base,  $
      UNAME='geometricSTDDevTB2' ,XOFFSET=730 ,YOFFSET=194  $
      ,SCR_XSIZE=40 ,SCR_YSIZE=21 ,/EDITABLE ,VALUE=[ '1.15' ]  $
      ,XSIZE=20 ,YSIZE=1)


  opticalCalcL = Widget_Label(tab1Base, UNAME='opticalCalcL'  $
      ,XOFFSET=536 ,YOFFSET=225 ,SCR_XSIZE=140 ,SCR_YSIZE=18  $
      ,/ALIGN_LEFT ,VALUE='-> Optical Calculations (T/F) :')


  opticalCalcTB = Widget_Text(tab1Base, UNAME='opticalCalcTB'  $
      ,XOFFSET=690 ,YOFFSET=221 ,SCR_XSIZE=40 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ 'T' ] ,XSIZE=20 ,YSIZE=1)


  particleSAorVL = Widget_Label(tab1Base, UNAME='particleSAorVL'  $
      ,XOFFSET=490 ,YOFFSET=251 ,SCR_XSIZE=186 ,SCR_YSIZE=18  $
      ,/ALIGN_LEFT ,VALUE='-> Particle Surface Area/Volume (2/3) :')


  particleSAorVTB = Widget_Text(tab1Base, UNAME='particleSAorVTB'  $
      ,XOFFSET=690 ,YOFFSET=249 ,SCR_XSIZE=40 ,SCR_YSIZE=21  $
      ,/EDITABLE ,VALUE=[ '3' ] ,XSIZE=20 ,YSIZE=1)


  commentsL = Widget_Label(tab1Base, UNAME='commentsL'  $
      ,XOFFSET=585, YOFFSET=285, SCR_XSIZE=100, SCR_YSIZE=20   $
      ,/ALIGN_LEFT, VALUE='Comments')


  commentsTB = Widget_Text(tab1Base, UNAME='commentsTB'  $
      ,XOFFSET=487, YOFFSET=305, SCR_XSIZE=260, SCR_YSIZE=20  $
      ,/EDITABLE, YSIZE=1, VALUE=['Enter Comments Here!'])


  graphingBase = Widget_Base(tab1Base, UNAME='graphingBase'  $
      ,XOFFSET=485 ,YOFFSET=400 ,SCR_XSIZE=275 ,SCR_YSIZE=25  $
      ,TITLE='IDL' ,COLUMN=2 ,/NONEXCLUSIVE)


  sizeDistroGraphsCB = Widget_Button(graphingBase,  $
      UNAME='sizeDistroGraphsCB' ,/ALIGN_LEFT ,VALUE='Size'+ $
      ' Distribution Graphs')
  ; Makes sizeDistroGraphsCB selected.
  WIDGET_CONTROL, sizeDistroGraphsCB, /SET_BUTTON
  ; Invisible label to tell if checkbox 1 is checked.
  sizeDistroGraphsInvisL = Widget_Label(baseWindow, UNAME='sizeDistroGraphsInvisL'  $
    ,XOFFSET=50 ,YOFFSET=50 ,/ALIGN_LEFT ,VALUE='1', SCR_XSIZE=0, SCR_YSIZE=0)


  historyTotalsGraphsCB = Widget_Button(graphingBase,  $
      UNAME='historyTotalsGraphsCB' ,/ALIGN_LEFT ,VALUE='History'+ $
      ' Totals Graphs')
  ; Makes historyTotalsGraphsCB selected.
  WIDGET_CONTROL, historyTotalsGraphsCB, /SET_BUTTON
  ; Invisible label to tell if checkbox 2 is checked.
  historyTotalsGraphsInvisL = Widget_Label(baseWindow, UNAME='historyTotalsGraphsInvisL'  $
    ,XOFFSET=50 ,YOFFSET=50 ,/ALIGN_LEFT ,VALUE='1', SCR_XSIZE=0, SCR_YSIZE=0)


  ; End of MicroPhysical Model Widget Intializations
  ;-------------------------------------------------------------------------------
  ;-------------------------------------------------------------------------------


  ; Realize all widgets.
  Widget_Control, baseWindow, /REALIZE

  ; Display toolkit on screen.
  XManager, 'PSC_Model_Toolkit', baseWindow, /NO_BLOCK

END