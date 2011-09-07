; Caleb Wherry
; Micro_Open_Inp_File.pro
; Created: October 8th, 2008
;
; Last Modified: October 27th, 2008
;

PRO Micro_Open_Inp_File, ev

  file=''
  val=''

  file=DIALOG_PICKFILE(filter='*.inp', TITLE='Please Select An Input File:')

  ; Skip loop if no file is selected.
  IF file EQ '' THEN GOTO, JUMP2  ; Line 215

  ; Open file based on dialog box selection.
  OPENR,1, file

    ; Populate *.inp textbox with correct name.
    file=STRCOMPRESS(STRING(file), /REMOVE_ALL)
    file1=STRSPLIT(file, '\', /EXTRACT)
    num=SIZE(file1)
    id=WIDGET_INFO(ev.top, FIND_BY_UNAME='inpFileNameTB')
    WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(file1[num[1]-1]), /REMOVE_ALL)

    ; First line in *.inp file.
    READF,1, val

    ; Number of Layers
    READF,1, val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='layersTB')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(fix(val1)), /REMOVE_ALL)

    ; Min and Max radius
    READF,1, val1, val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='minRadiusTB')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.3)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='maxRadiusTB')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.3)'), /REMOVE_ALL)

    ; Max Integration Time Step
    READF,1, val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='maxIntegTimeStepTB')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Time Units
    READF,1, val
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='timeUnitsTB')
      val=STRSPLIT(val, /EXTRACT)
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val[0]), /REMOVE_ALL)
      val=''

    ; Integration Start Time
    READF,1, val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='integStartTimeTB')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Integration Stop Time
    READF,1, val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='integrationStopTimeTB')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Time Interval Between Plot/Print of Output
    READF,1, val1, val2, val3
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='timeIntervalTB1')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='timeIntervalTB2')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='timeIntervalTB3')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val3,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Bottom Layer Potential Temp and P.T. Increment
    READF,1, val1, val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='bottomLayerPotenTempTB')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='potenTempIncrTB')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Temperature Calculation
    READF,1, val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempCalcInvisL')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(fix(val1)), /REMOVE_ALL)
      IF val1 EQ 0 THEN BEGIN
        id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempCalcRB1')
        WIDGET_CONTROL, id, /SET_BUTTON
      END
      IF val1 EQ 1 THEN BEGIN
        id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempCalcRB2')
        WIDGET_CONTROL, id, /SET_BUTTON
      END
      IF val1 EQ 2 THEN BEGIN
        id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempCalcRB3')
        WIDGET_CONTROL, id, /SET_BUTTON
      END

    ; Temperature/Pressure History File
    READF,1, val
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempPressHisFileTB')
      val=STRSPLIT(val, /EXTRACT)
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val[0]), /REMOVE_ALL)
      val=''

    ; Ramp Time/Temperature 1
    READF,1, val1, val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime1TB1')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime1TB2')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Ramp Time/Temperature 2
    READF,1, val1, val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime2TB1')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime2TB2')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Ramp Time/Temperature 3
    READF,1, val1, val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime3TB1')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime3TB2')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Ramp Time/Temperature 4
    READF,1, val1, val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime4TB1')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='rampTime4TB2')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Temp Sine Ocillation Period & Amplitude
    READF,1, val1, val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempSineOsciTB1')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempSineOsciTB2')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Mixing Ratio of H20 In Each Layer
    READF,1, val1, val2, val3, val4
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixH20RatioTB1')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixH20RatioTB2')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixH20RatioTB3')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val3,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixH20RatioTB4')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val4,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Mixing Ratio of NA In Each Layer
    READF,1, val1, val2, val3, val4
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixNARatioTB1')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixNARatioTB2')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixNARatioTB3')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val3,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='mixNARatioTB4')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val4,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Number Density Sulf. Aerosols
    READF,1, val1, val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='numDensSulfAeroTB1')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='numDensSulfAeroTB2')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Median Radius
    READF,1, val1, val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='medianRadiusTB1')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.3)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='medianRadiusTB2')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.3)'), /REMOVE_ALL)

    ; Geometric Std. Dev.
    READF,1, val1, val2
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='geometricSTDDevTB1')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='geometricSTDDevTB2')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val2,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Latitude
    READF,1, val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='latitudeTB')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)

    ; Optical Calculations
    READF,1, val
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='opticalCalcTB')
      val=STRSPLIT(val, /EXTRACT)
      IF val[0] EQ '.t.' THEN WIDGET_CONTROL, id, SET_VALUE='T'
      IF val[0] EQ '.f.' THEN WIDGET_CONTROL, id, SET_VALUE='F'
      val=''

    ; Particle SA or Volume
    READF,1, val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='particleSAorVTB')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(fix(val1)), /REMOVE_ALL)

    ; Correction Factor for Homogeneous NAT Nucleation
    READF,1, val1
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='homogeNATFactoL')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val1,FORMAT='(F8.2)'), /REMOVE_ALL)
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='homogeNATFactoSB')
      WIDGET_CONTROL, id, SET_VALUE=val1*100

    ; Comments
    READF, 1, val
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='commentsTB')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val))

  CLOSE,1

  JUMP2:
END