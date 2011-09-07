; Caleb Wherry
; Micro_Size_Distros_Graphs.pro
; Microphysical Model Output Size Distribution Graphs
; Created: October 8th, 2008
;
; Last Modified: October 27th, 2008
;

; Function to retrieve the ending time of simulation.
; Used to set historyArr size.
FUNCTION TIME_STEP_CALC, inpFile, rootDir, subDir2
  file = inpFile+'.inp'

  file = FILEPATH(file, ROOT=[rootDir],SUBDIR=[subDir2])

  x1 = ''
  x5 = ''

  OPENR,1,file

  READF,1,x1
  READF,1,x2
  READF,1,x3
  READF,1,x4
  READF,1,x5
  READF,1,x6
  READF,1,x7

  value = FIX(x7)

  CLOSE,1

  RETURN, value
END
;
;

PRO Micro_Size_Distros_Graphs, val1

  ; Format inpFile file.
  val1 = STRSPLIT(val1, '\.', /EXTRACT, /REGEX)
  inpFile = val1[0]

  ; Set plotting to PostScript
  SET_PLOT, 'PS'

  ; Make IDL's plotting area hold 2 columns and 3 rows of plots:
  !P.MULTI = [0, 3, 2]

  ; Number of bin and time.
  binNum = 60
  time = 0.0

  ; Initialize and build *.dat file names for each layer from *.inp file name.
  layer1 = ''
  layer1H = ''
  layer1File = ''
  layer1FileH = ''

  layer2 = ''
  layer2H = ''
  layer2File = ''
  layer2FileH = ''

  layer3 = ''
  layer3H = ''
  layer3File = ''
  layer3FileH = ''

  layer4 = ''
  layer4H = ''
  layer4File = ''
  layer4FileH = ''

  ; D file names.
  layer1 = inpFile + '-D01.DAT'
  layer2 = inpFile + '-D02.DAT'
  layer3 = inpFile + '-D03.DAT'
  layer4 = inpFile + '-D04.DAT'

  ; H file names.
  layer1H = inpFile + '-H01.DAT'
  layer2H = inpFile + '-H02.DAT'
  layer3H = inpFile + '-H03.DAT'
  layer4H = inpFile + '-H04.DAT'

  ; Path where subdir resides.
  rootDir = '.\'
  ;rootDir = '.\Desktop\PSC_Model\'

  ; Subdirectory for *.dat files.
  subDir = '.\Simulations\'+inpFile+'\Output\'

  ; Subdirectory for *.inp files.
  subDir2 = '.\Simulations\'+inpFile

  ; File paths for D files.
  layer1File = FILEPATH(layer1, ROOT=[rootDir],SUBDIR=[subDir])
  layer2File = FILEPATH(layer2, ROOT=[rootDir],SUBDIR=[subDir])
  layer3File = FILEPATH(layer3, ROOT=[rootDir],SUBDIR=[subDir])
  layer4File = FILEPATH(layer4, ROOT=[rootDir],SUBDIR=[subDir])

  ; File paths for H files.
  layer1FileH = FILEPATH(layer1H, ROOT=[rootDir],SUBDIR=[subDir])
  layer2FileH = FILEPATH(layer2H, ROOT=[rootDir],SUBDIR=[subDir])
  layer3FileH = FILEPATH(layer3H, ROOT=[rootDir],SUBDIR=[subDir])
  layer4FileH = FILEPATH(layer4H, ROOT=[rootDir],SUBDIR=[subDir])

  ; 2D Array(binNum x 7) With Column Names Below:
  ; Radius, STS, SAT, NAT, ICE, Total DnDR, Anticumulative
  cloudArr =  FLTARR(binNum,7)

  value = TIME_STEP_CALC(inpFile, rootDir, SubDir2)

  ; 2D Array(value*10+1 x 69)
  historyARR = FLTARR(value*10+1, 69)

  ;-------------------------------------------------------------
  ; Calculations For Each Layer
  ;-------------------------------------------------------------

  ; FOR Loop for Each Layer: 1-4
  FOR num=1, 4 DO BEGIN

    ; IF block to open correct layer file: 1-4
    IF num EQ 1 THEN BEGIN
      OPENR,1,layer1File
      OPENR,2,layer1FileH
    END
    IF num EQ 2 THEN BEGIN
      OPENR,1,layer2File
      OPENR,2,layer2FileH
    END
    IF num EQ 3 THEN BEGIN
      OPENR,1,layer3File
      OPENR,2,layer3FileH
    END
    IF num EQ 4 THEN BEGIN
      OPENR,1,layer4File
      OPENR,2,layer4FileH
    END

    ; Covert num(integer) to layerNum(string) and take out spaces.
    layerNum = STRCOMPRESS(STRING(num), /REMOVE_ALL)

    ; Set plotting to PostScript file:
    DEVICE, FILENAME=rootDir+'Simulations\'+inpFile+'\Size_Distro_Graphs\'+inpFile+'-L'+layerNum+'.ps', /landscape

    ; WHILE loop to read in the history file.
    k=0
    WHILE NOT EOF(2) DO BEGIN
      READF,2, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20, $
               y21, y22, y23, y24, y25, y26, y27, y28, y29, y30, y31, y32, y33, y34, y35, y36, y37, y38, y39, y40, $
               y41, y42, y43, y44, y45, y46, y47, y48, y49, y50, y51, y52, y53, y54, y55, y56, y57, y58, y59, y60, $
               y61, y62, y63, y64, y65, y66, y67, y68, y69

      historyArr[k,0] = y1      ; 1: Time (Days, Hours, Minutes, Seconds; as specified in the input file)
      historyArr[k,1] = y2      ; 2: Temperature (K)
      historyArr[k,2] = y3      ; 3: Total number density of frozen sulfate particles (1/ccm)
      historyArr[k,3] = y4      ; 4: Total number density of PSC 1 particles (1/ccm)
      historyArr[k,4] = y5      ; 5: Total number density of PSC 2 particles (1/ccm)
      historyArr[k,5] = y6      ; 6: Saturation ratio of nitric acid over NAT
      historyArr[k,6] = y7      ; 7: Saturation ratio of water vapor over ICE
      historyArr[k,7] = y8      ; 8: Mixing raito of water vapor (ppmv)
      historyArr[k,8] = y9      ; 9: Mixing ratio of nitric acid vapor (ppbv)
      historyArr[k,9] = y10     ; 10: Total (gas phase and condensed phase) mix. ratio, water (ppmv)
      historyArr[k,10] = y11    ; 11: Total (gas phase and condensed phase) mix. ratio, nitric acid (ppbv)
      historyArr[k,11] = y12    ; 12: Ratio of total to initial mixing ratio of water (%)
      historyArr[k,12] = y13    ; 13: Ratio of total to initial mixing ratio of nitric acid (%)
      historyArr[k,13] = y14    ; 14: Total volume of frozen sulfate particles (micron**3/ccm)
      historyArr[k,14] = y15    ; 15: Total volume of PSC 1a particles (micron**3/ccm)
      historyArr[k,15] = y16    ; 16: Total volume of PSC 2 particles (micron**3/ccm)
      historyArr[k,16] = y17    ; 17: Total volume of all particles (micron**3/ccm)
      historyArr[k,17] = y18    ; 18: Total number density of PSC 1b (STS) particles (1/ccm)
      historyArr[k,18] = y19    ; 19: Total volume of PSC 1b (STS) particles (micron**3/ccm)
      historyArr[k,19] = y20    ; 20: Sulfuric acid weight fraction liquid PSC 1b aerosols [0;1]
      historyArr[k,20] = y21    ; 21: Mean radius liquid PSC 1b aerosols (micron)
      historyArr[k,21] = y22    ; 22: Median radius frozen sulfate aerosols (micron)
      historyArr[k,22] = y23    ; 23: Median radius PSC 1 (micron)
      historyArr[k,23] = y24    ; 24: Median radius PSC 2 (micron)
      historyArr[k,24] = y25    ; 25: Nitric acid weight fraction liquid PSC 1b aerosols [0;1]
      historyArr[k,25] = y26    ; 26: NAT condensation temperature (K)
      historyArr[k,26] = y27    ; 27: Ice frost point temperature (K)
      historyArr[k,27] = y28    ; 28: Air pressure (hPa)
      historyArr[k,28] = y29    ; 29: Total mixing ration of sulfuric acid (ppb)
      historyArr[k,29] = y30    ; 30: Optical parameter 1
      historyArr[k,30] = y31    ; 31: Optical parameter 2
      historyArr[k,31] = y32    ; 32: Optical parameter 3
      historyArr[k,32] = y33    ; 33: Moleratio
      historyArr[k,33] = y34    ; 34-46: Cumulated size distribution
      historyArr[k,34] = y35    ;
      historyArr[k,35] = y36    ;
      historyArr[k,36] = y37    ;
      historyArr[k,37] = y38    ;
      historyArr[k,38] = y39    ;
      historyArr[k,39] = y40    ;
      historyArr[k,40] = y41    ;
      historyArr[k,41] = y42    ;
      historyArr[k,42] = y43    ;
      historyArr[k,43] = y44    ;
      historyArr[k,44] = y45    ;
      historyArr[k,45] = y46    ;
                              ; 47-53: AEROSOL_BACKSCATTER_RATIO(J,L),J=1,N_WAVE)
      historyArr[k,46] = y47    ; 354 nm
      historyArr[k,47] = y48    ; 449 nm
      historyArr[k,48] = y49    ; 532 nm
      historyArr[k,49] = y50    ; 779 nm
      historyArr[k,50] = y51    ; 1545 nm
      historyArr[k,51] = y52    ; 1022 nm
      historyArr[k,52] = y53    ; 1064 nm
      historyArr[k,53] = y54    ; 53-60: EXTINCTION(J,L),J=1,N_WAVE)
      historyArr[k,54] = y55    ;
      historyArr[k,55] = y56    ;
      historyArr[k,56] = y57    ;
      historyArr[k,57] = y58    ;
      historyArr[k,58] = y59    ;
      historyArr[k,59] = y60    ;
      historyArr[k,60] = y61    ; 61-67: DEPOL(J,L),J=1,N_WAVE)
      historyArr[k,61] = y62    ;
      historyArr[k,62] = y63    ;
      historyArr[k,63] = y64    ;
      historyArr[k,64] = y65    ;
      historyArr[k,65] = y66    ;
      historyArr[k,66] = y67    ;
      historyArr[k,67] = y68    ; 68: Potential Temperature
      historyArr[k,68] = y69    ; 69: Altitude

      k=k+1

    ENDWHILE

    ; Close history file:
    CLOSE, 2

    ; WHILE loop to read in the data file.
    WHILE NOT EOF(1) DO BEGIN

      ; Read time variable and format it for while loop.
      ; FIX returns only integer part of argument.
      ; STRING converts argument to string
      ; STRCOMPRESS takes out all spaces when '/REMOVE_ALL' flag is included. Otherwise
      ;   it takes out all spaces except one preceeding space.
      READF,1,time
      time2 = STRCOMPRESS(STRING(FIX(time)), /REMOVE_ALL)

      ; Read rest of data.
      FOR index=0,binNum-1 DO BEGIN
        READF,1,x1,x2,x3,x4,x5,x6,x7

        cloudArr[index, 0] = x1
        cloudArr[index, 1] = x2
        cloudArr[index, 2] = x3
        cloudArr[index, 3] = x4
        cloudArr[index, 4] = x5
        cloudArr[index, 5] = x6
        cloudArr[index, 6] = x7
      ENDFOR

      !P.MULTI = [0, 3, 2]

      ; Plots 1-5
      PLOT_OI, cloudArr[*,0], cloudArr[*,1], XRANGE=[0.0001,100], TITLE='STS', XTITLE='Radius (microns)', YTITLE='dN/dr'
      PLOT_OI, cloudArr[*,0], cloudArr[*,2], XRANGE=[0.0001,100], TITLE='SAT', XTITLE='Radius (microns)', YTITLE='dN/dr'
      PLOT_OI, cloudArr[*,0], cloudArr[*,3], XRANGE=[0.0001,100], TITLE='NAT', XTITLE='Radius (microns)', YTITLE='dN/dr'
      PLOT_OI, cloudArr[*,0], cloudArr[*,4], XRANGE=[0.0001,100], TITLE='ICE', XTITLE='Radius (microns)', YTITLE='dN/dr'
      PLOT_OI, cloudArr[*,0], cloudArr[*,5], XRANGE=[0.0001,100], TITLE='Total dN/dr', XTITLE='Radius (microns)', YTITLE='dN/dr'

      ; Set text size for Plot 6:
      !P.CHARSIZE = .7

      ; 'Plot 6' Variables
      XYOUTS, /NORMAL,.69,.47, STRING(FORMAT='("Time:", A, " hours")', time2)
      XYOUTS, /NORMAL,.69,.44, STRING(FORMAT='("Temperature:", A, " K")', STRCOMPRESS(STRING(historyArr[time*10,1])))
      XYOUTS, /NORMAL,.69,.41, STRING(FORMAT='("Pressure:", A, " mb")', STRCOMPRESS(STRING(historyArr[time*10, 27])))
      XYOUTS, /NORMAL,.69,.38, STRING(FORMAT='("Potential Temperature:", A, " K")', STRCOMPRESS(STRING(historyArr[time*10, 67])))
      XYOUTS, /NORMAL,.69,.35, STRING(FORMAT='("Altitude:", A, " km")', STRCOMPRESS(STRING(historyArr[time*10,68])))
      XYOUTS, /NORMAL,.69,.32, STRING(FORMAT='("HNO!d3!n:", A, " ppbv")', STRCOMPRESS(STRING(historyArr[time*10, 8])))
      XYOUTS, /NORMAL,.69,.29, STRING(FORMAT='("H!d2!nO:", A, " ppmv")', STRCOMPRESS(STRING(historyArr[time*10, 7])))
      XYOUTS, /NORMAL,.69,.26, STRING(FORMAT='("V!dSTS!n:", A, " micron^3/cc")', STRCOMPRESS(STRING(historyArr[time*10, 18])))
      XYOUTS, /NORMAL,.69,.23, STRING(FORMAT='("V!dSAT!n:", A, " micron^3/cc")', STRCOMPRESS(STRING(historyArr[time*10, 13])))
      XYOUTS, /NORMAL,.69,.20, STRING(FORMAT='("V!dNAT!n:", A, " micron^3/cc")', STRCOMPRESS(STRING(historyArr[time*10, 14])))
      XYOUTS, /NORMAL,.69,.17, STRING(FORMAT='("V!dICE!n:", A, " micron^3/cc")', STRCOMPRESS(STRING(historyArr[time*10, 15])))
      XYOUTS, /NORMAL,.69,.14, STRING(FORMAT='("T!dNAT!n:", A, " K")', STRCOMPRESS(STRING(historyArr[time*10, 25])))
      XYOUTS, /NORMAL,.69,.11, STRING(FORMAT='("T!dICE!n:", A, " K")', STRCOMPRESS(STRING(historyArr[time*10, 26])))
      XYOUTS, /NORMAL,.69,.08, STRING(FORMAT='("R(532nm):", A, "")', STRCOMPRESS(STRING(historyArr[time*10, 48]+1)))
      XYOUTS, /NORMAL,.69,.05, STRING(FORMAT='("R(1064nm):", A, "")', STRCOMPRESS(STRING(historyArr[time*10, 52]+1)))
      XYOUTS, /NORMAL,.69,.02, STRING(FORMAT='("Aer. Depol.(532nm):", A, "")', STRCOMPRESS(STRING(historyArr[time*10, 62])))

      ; Change text size back to default.
      !P.CHARSIZE = 0

    ENDWHILE

    ; Close plotting device.
    DEVICE, /CLOSE

    ; Close open file.
    CLOSE,1

  ENDFOR
  ;--------------------------------------------------------------

; Closing message.
PRINT, 'Microphysical output size distribution graphs have been completed!'

END