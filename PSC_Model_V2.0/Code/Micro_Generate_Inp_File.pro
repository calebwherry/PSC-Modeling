; Caleb Wherry
; Micro_Generate_Inp_File.pro
; Created: October 8th, 2008
;
; Last Modified October 27th, 2008
;


PRO Micro_Generate_Inp_File, ev

  ; Grab new *.inp file name
  id=WIDGET_INFO(ev.top, FIND_BY_UNAME='inpFileNameTB')
  WIDGET_CONTROL, id, GET_VALUE=file

  ; Create copy of input file name and manipulate it.
  file2=file
  file2 = STRSPLIT(file, '\.', /EXTRACT, /REGEX)
  file2 = STRCOMPRESS(STRING(file2[0]), /REMOVE_ALL)

  ; Test to see if a case built off of this filename already exists
  exists = FILE_TEST('.\Simulations\'+file2)

  IF exists EQ 1 THEN BEGIN
    yes_no = DIALOG_MESSAGE('Overwrite Existing File?', /QUESTION, /CENTER)
    IF yes_no EQ 'No' THEN GOTO, JUMP1 ; Line 254
  END

  ; Create folder for new *.inp file to reside in.
  FILE_MKDIR, FILEPATH(file2, ROOT_DIR=['.\'], SUBDIR=['Simulations'])

  ; Build file structure for simulation.
  FILE_MKDIR, FILEPATH('Size_Distro_Graphs', ROOT_DIR=['.\'], SUBDIR=['Simulations', file2])
  FILE_MKDIR, FILEPATH('History_Totals_Graphs', ROOT_DIR=['.\'], SUBDIR=['Simulations', file2])
  FILE_MKDIR, FILEPATH('Output', ROOT_DIR=['.\'], SUBDIR=['Simulations', file2])

  ; String Array to hold log values.
  logArr = STRARR(24)

  ; Name of case.
  logArr[0]=file2

  ; Open *.inp file.
  OPENW,1, FILEPATH(file, ROOT_DIR='.\', SUBDIR=['Simulations',file2])

    ; First line in *.inp file.
    PRINTF,1, 'Constant cooling - heating, PSC model E version, Homogeneous freezing'

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='layersTB')
      WIDGET_CONTROL, id, GET_VALUE=val
    ; Number of Layers
    PRINTF,1, val

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='minRadiusTB')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='maxRadiusTB')
      WIDGET_CONTROL, id, GET_VALUE=val2
    ; Min and Max radius
    PRINTF,1, val1+', '+val2

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='maxIntegTimeStepTB')
      WIDGET_CONTROL, id, GET_VALUE=val
    ; Max Integration Time Step
    PRINTF,1, val

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='timeUnitsTB')
      WIDGET_CONTROL, id, GET_VALUE=val
    ; Time Units
    PRINTF,1, val

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='integStartTimeTB')
      WIDGET_CONTROL, id, GET_VALUE=val
    ; Integration Start Time
    PRINTF,1, val

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='integrationStopTimeTB')
      WIDGET_CONTROL, id, GET_VALUE=val
    ; Integration Stop Time
    PRINTF,1, val

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='timeIntervalTB1')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='timeIntervalTB2')
      WIDGET_CONTROL, id, GET_VALUE=val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='timeIntervalTB3')
      WIDGET_CONTROL, id, GET_VALUE=val3
    ; Time Interval Between Plot/Print of Output
    PRINTF,1, val1+', '+val2+', '+val3

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='bottomLayerPotenTempTB')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='potenTempIncrTB')
      WIDGET_CONTROL, id, GET_VALUE=val2
    ; Bottom Layer Potential Temp and P.T. Increment
    PRINTF,1, val1+', '+val2
    logArr[23]=val1

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempCalcInvisL')
      WIDGET_CONTROL, id, GET_VALUE=val
    ; Temperature Calculation
    PRINTF,1, val

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempPressHisFileTB')
      WIDGET_CONTROL, id, GET_VALUE=val
    ; Temperature/Pressure History File
    PRINTF,1, val
    logArr[1]=val

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime1TB1')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime1TB2')
      WIDGET_CONTROL, id, GET_VALUE=val2
    ; Ramp Time/Temperature 1
    PRINTF,1, val1+', '+val2
    logArr[2]=val1
    logArr[3]=val2

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime2TB1')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime2TB2')
      WIDGET_CONTROL, id, GET_VALUE=val2
    ; Ramp Time/Temperature 2
    PRINTF,1, val1+', '+val2
    logArr[4]=val1
    logArr[5]=val2

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime3TB1')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime3TB2')
      WIDGET_CONTROL, id, GET_VALUE=val2
    ; Ramp Time/Temperature 3
    PRINTF,1, val1+', '+val2
    logArr[6]=val1
    logArr[7]=val2

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime4TB1')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime4TB2')
      WIDGET_CONTROL, id, GET_VALUE=val2
    ; Ramp Time/Temperature 4
    PRINTF,1, val1+', '+val2
    logArr[8]=val1
    logArr[9]=val2

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempSineOsciTB1')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempSineOsciTB2')
      WIDGET_CONTROL, id, GET_VALUE=val2
    ; Temp Sine Ocillation Period & Amplitude
    PRINTF,1, val1+', '+val2

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixH20RatioTB1')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixH20RatioTB2')
      WIDGET_CONTROL, id, GET_VALUE=val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixH20RatioTB3')
      WIDGET_CONTROL, id, GET_VALUE=val3
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixH20RatioTB4')
      WIDGET_CONTROL, id, GET_VALUE=val4
    ; Mixing Ratio of H20 In Each Layer
    PRINTF,1, val1+', '+val2+', '+val3+', '+val4
    logArr[10]=val1
    logArr[11]=val2
    logArr[12]=val3
    logArr[13]=val4

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixNARatioTB1')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixNARatioTB2')
      WIDGET_CONTROL, id, GET_VALUE=val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixNARatioTB3')
      WIDGET_CONTROL, id, GET_VALUE=val3
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixNARatioTB4')
      WIDGET_CONTROL, id, GET_VALUE=val4
    ; Mixing Ratio of NA In Each Layer
    PRINTF,1, val1+', '+val2+', '+val3+', '+val4
    logArr[14]=val1
    logArr[15]=val2
    logArr[16]=val3
    logArr[17]=val4

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='numDensSulfAeroTB1')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='numDensSulfAeroTB2')
      WIDGET_CONTROL, id, GET_VALUE=val2
    ; Number Density Sulf. Aerosols
    PRINTF,1, val1+', '+val2
    logArr[18]=val1

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='medianRadiusTB1')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='medianRadiusTB2')
      WIDGET_CONTROL, id, GET_VALUE=val2
    ; Median Radius
    PRINTF,1, val1+', '+val2
    logArr[19]=val1

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='geometricSTDDevTB1')
      WIDGET_CONTROL, id, GET_VALUE=val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='geometricSTDDevTB2')
      WIDGET_CONTROL, id, GET_VALUE=val2
    ; Geometric Std. Dev.
    PRINTF,1, val1+', '+val2
    logArr[20]=val1

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='latitudeTB')
      WIDGET_CONTROL, id, GET_VALUE=val
    ; Latitude
    PRINTF,1, val

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='opticalCalcTB')
      WIDGET_CONTROL, id, GET_VALUE=val
      val=STRCOMPRESS(STRING(val), /REMOVE_ALL)
    ; Optical Calculations
    IF (val EQ 'T') OR (val EQ 't') THEN PRINTF,1, '.t.'
    IF (val EQ 'F') OR (val EQ 'f') THEN PRINTF,1, '.f.'

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='particleSAorVTB')
      WIDGET_CONTROL, id, GET_VALUE=val
    ; Particle SA or Volume
    PRINTF,1, val

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='homogeNATFactoL')
      WIDGET_CONTROL, id, GET_VALUE=val
    ; Correction Factor for Homogeneous NAT Nucleation
    PRINTF,1, val
    logArr[21]=val

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='commentsTB')
      WIDGET_CONTROL, id, GET_VALUE=val
      val=STRCOMPRESS(STRING(val))
    IF val EQ '' THEN BEGIN
      PRINTF,1, 'Comments go here!'
      logArr[22]='No'
    ENDIF ELSE BEGIN
      ; Comments for this specific case.
      PRINTF,1, val
      logArr[22]='Yes'
    ENDELSE

  ; Close open file.
  CLOSE,1

  ; Create log file entry for new case.
  OPENU,1, 'Microphysical_Cases_Log.csv', /APPEND

    PRINTF,1, logArr[0]+', '+logArr[1]+', '+logArr[2]+', '+logArr[3]+', '+logArr[4]  $
              +', '+logArr[5]+', '+logArr[6]+', '+logArr[7]+', '+logArr[8]+', '+logArr[9]+', '+logArr[23] $
              +', '+logArr[10]+', '+logArr[11]+', '+logArr[12]+', '+logArr[13]+', '+logArr[14] $
              +', '+logArr[15]+', '+logArr[16]+', '+logArr[17]+', '+logArr[18]+', '+logArr[19] $
              +', '+logArr[20]+', '+logArr[21]+', '+logArr[22]

  CLOSE,1

  JUMP1:
END