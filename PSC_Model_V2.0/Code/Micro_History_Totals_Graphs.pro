; Caleb Wherry
; Micro_History_Totals_Graphs.pro
; Plots totals in history file as a function of time.
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

PRO Micro_History_Totals_Graphs, val2

  ; Format inpFile file.
  val2 = STRSPLIT(val2, '\.', /EXTRACT, /REGEX)
  inpFile = val2[0]

  ; Set plotting to PostScript
  SET_PLOT, 'PS'

  ; Make IDL's plotting area hold 2 columns and 2 rows of plots:
  !P.MULTI = [0, 2, 2]

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
  layer1FileH = FILEPATH(layer1H, ROOT=[rootDir],SUBDIR=[subDir])
  layer2FileH = FILEPATH(layer2H, ROOT=[rootDir],SUBDIR=[subDir])
  layer3FileH = FILEPATH(layer3H, ROOT=[rootDir],SUBDIR=[subDir])
  layer4FileH = FILEPATH(layer4H, ROOT=[rootDir],SUBDIR=[subDir])

  ; Retrieve value for historyArr length.
  value = TIME_STEP_CALC(inpFile, rootDir, SubDir2)

  ; 2D Array(value*10+1 x 69)
  historyARR = FLTARR(value*10+1, 69)

  FOR num=1, 4 DO BEGIN

    ; IF block to open correct layer file: 1-4
    IF num EQ 1 THEN BEGIN
      OPENR,1,layer1FileH
    END
    IF num EQ 2 THEN BEGIN
      OPENR,1,layer2FileH
    END
    IF num EQ 3 THEN BEGIN
      OPENR,1,layer3FileH
    END
    IF num EQ 4 THEN BEGIN
      OPENR,1,layer4FileH
    END

    ; Covert num(integer) to layerNum(string) and take out spaces.
    layerNum = STRCOMPRESS(STRING(num), /REMOVE_ALL)

    ; Set plotting to PostScript file:
    DEVICE, FILENAME=rootDir+'Simulations\'+inpFile+'\History_Totals_Graphs\'+inpFile+'-L'+layerNum+'-HisTot'+'.ps', /landscape

    ; WHILE loop to read in the history file.
    k=0
    WHILE NOT EOF(1) DO BEGIN
      READF,1, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20, $
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
      historyArr[k,46] = y47+1  ; 354 nm
      historyArr[k,47] = y48+1  ; 449 nm
      historyArr[k,48] = y49+1  ; 532 nm
      historyArr[k,49] = y50+1  ; 779 nm
      historyArr[k,50] = y51+1  ; 1545 nm
      historyArr[k,51] = y52+1  ; 1022 nm
      historyArr[k,52] = y53+1  ; 1064 nm
      historyArr[k,53] = y54    ; 54-60: EXTINCTION(J,L),J=1,N_WAVE)
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
    CLOSE, 1

    potentialTemp = STRCOMPRESS(STRING(historyArr[0,67]), /REMOVE_ALL)

    ; Page 1
    PLOT,historyArr[*,0],historyArr[*,1],/YNOZ,XTITLE='Hours',YTITLE='Temperature (K)', TITLE='Temperature'
    PLOT,historyArr[*,0],historyArr[*,27],/YNOZ,XTITLE='Hours',YTITLE='Pressure (mb)', TITLE='Pressure'
    PLOT,historyArr[*,0],historyArr[*,7],/YNOZ,XTITLE='Hours',YTITLE='Water Vapor Mixing Ratio (ppm)', TITLE='Mixing Ratio of Water Vapor'
    PLOT,historyArr[*,0],historyArr[*,8],/YNOZ,XTITLE='Hours',YTITLE='Nitric Acid Vapor Mixing Ratio (ppb)', TITLE='Mixing Ratio of Nitric Acid Vapor'
    XYOUTS, /NORMAL, .38, 1, 'Potential Temperature: ' + potentialTemp + ' K'

    ; Page 2
    PLOT,historyArr[*,0],historyArr[*,2],/YNOZ,XTITLE='Hours',YTITLE='Total Number Density (1/ccm', TITLE='Total Number Density of ICE'
    PLOT,historyArr[*,0],historyArr[*,3],/YNOZ,XTITLE='Hours',YTITLE='Total Number Density (1/ccm)', TITLE='Total Number Density of NAT'
    PLOT,historyArr[*,0],historyArr[*,4],/YNOZ,XTITLE='Hours',YTITLE='Total Number Density (1/ccm)', TITLE='Total Number Density of SAT'
    PLOT,historyArr[*,0],historyArr[*,17],/YNOZ,XTITLE='Hours',YTITLE='Total Number Density (1/ccm)', TITLE='Total Number Density of STS'
    XYOUTS, /NORMAL, .38, 1, 'Potential Temperature: ' + potentialTemp + ' K'

    ; Page 3
    PLOT,historyArr[*,0],historyArr[*,9],/YNOZ,XTITLE='Hours',YTITLE='Total Water Mixing Ratio (ppmv)', TITLE='Total Water Mixing Ratio'
    PLOT,historyArr[*,0],historyArr[*,10],/YNOZ,XTITLE='Hours',YTITLE='Total Nitric Acid Mixing Ratio (ppbv)', TITLE='Total Nitric Acid Mixing Ratio'
    PLOT,historyArr[*,0],historyArr[*,25],/YNOZ,XTITLE='Hours',YTITLE='NAT Condensation Temperature (K)', TITLE='NAT Condensation Temperature'
    PLOT,historyArr[*,0],historyArr[*,26],/YNOZ,XTITLE='Hours',YTITLE='ICE Condensation Temperature (K)', TITLE='ICE Condensation Temperature'
    XYOUTS, /NORMAL, .38, 1, 'Potential Temperature: ' + potentialTemp + ' K'

    ; Page 4
    PLOT,historyArr[*,0],historyArr[*,13],/YNOZ,XTITLE='Hours',YTITLE='V!dD!n (microns^3/ccm)', TITLE='Total Volume of SAT'
    PLOT,historyArr[*,0],historyArr[*,14],/YNOZ,XTITLE='Hours',YTITLE='V!dD!n (microns^3/ccm)', TITLE='Total Volume of NAT'
    PLOT,historyArr[*,0],historyArr[*,15],/YNOZ,XTITLE='Hours',YTITLE='V!dD!n (microns^3/ccm)', TITLE='Total Volume of ICE'
    PLOT,historyArr[*,0],historyArr[*,18],/YNOZ,XTITLE='Hours',YTITLE='V!dD!n (microns^3/ccm)', TITLE='Total Volume of STS'
    XYOUTS, /NORMAL, .38, 1, 'Potential Temperature: ' + potentialTemp + ' K'

    ; Page 5
    PLOT,historyArr[*,0],historyArr[*,48],/YNOZ,XTITLE='Hours',YTITLE='R (532 nm)', TITLE='BackScatter Ratio (532 nm)'
    PLOT,historyArr[*,0],historyArr[*,52],/YNOZ,XTITLE='Hours',YTITLE='R (1064 nm)', TITLE='BackScatter Ratio (1064 nm)'
    PLOT,historyArr[*,0],historyArr[*,62],/YNOZ,XTITLE='Hours',YTITLE='Aer. Depol. (532 nm) %', TITLE='Depolarization (532 nm)'
    PLOT,historyArr[*,0],historyArr[*,66],/YNOZ,XTITLE='Hours',YTITLE='Aer. Depol. (1064 nm) %', TITLE='Depolarization (1064 nm)'
    XYOUTS, /NORMAL, .38, 1, 'Potential Temperature: ' + potentialTemp + ' K'

    DEVICE, /CLOSE
  ENDFOR

; Closing message.
PRINT, 'Microphysical history totals graphs have been completed!'

END



