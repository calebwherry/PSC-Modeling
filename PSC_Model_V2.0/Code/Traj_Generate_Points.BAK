; Caleb Wherry
; Traj_Generate_Points.pro
; Created: October 8th, 2008
;
; Last Modified: October 30th, 2008
;

PRO Traj_Generate_Points, ev

  ; COMMON Block to share data between procedures
  ; ---------------------------------------------
  COMMON SHARE1, zalt,       $
  				 xlon1,      $
				 xlat1,      $
				 time1,      $
				 jday,       $
				 utcTime,    $
				 theta
  ; ----------------------------------------------

  ; COMMON Block to share data between Traj_Generate_Points.pro
  ;   and Traj_Generate_Dat_File.pro
  ; ---------------------------------------------
  COMMON SHARE2, altlev,  $
                 xlon2,   $
                 xlat2,   $
                 jday1,   $
                 utcTime1
  ; ---------------------------------------------

  X = ev.X
  Y = ev.Y
  PRESS = FIX(ev.PRESS)
  RELEASE = FIX(ev.RELEASE)

  ; Choose only the points in the wanted graph area.
  IF X GE 44 AND X LE 600 AND Y GE 44 AND Y LE 164 THEN BEGIN

	; When mouse "pressed", set labels to click value.
    IF PRESS EQ 1 THEN BEGIN
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='clickXinvisL')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(X))

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='clickYinvisL')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(Y))
    END

    ; When mouse "released", get previous X & Y and construct boxes.
    IF RELEASE EQ 1 THEN BEGIN
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='clickXinvisL')
      WIDGET_CONTROL, id, GET_VALUE=prevX

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='clickYinvisL')
      WIDGET_CONTROL, id, GET_VALUE=prevY

      prevX = FIX(prevX)
      prevY = FIX(prevY)

      num       = N_ELEMENTS(xlon1)
      num2      = N_ELEMENTS(zalt)

      altlev    = FLTARR(num2)
      xlon2     = FLTARR(num*num2)
      xlat2     = FLTARR(num*num2)
      jday1     = DBLARR(num*num2)
      utcTime1  = STRARR(num*num2)
      potenTemp = FLTARR(num*num2)

	  ; ID of pointsTB
	  id=WIDGET_INFO(ev.top, FIND_BY_UNAME='pointsTB')

      ; Event if only one point is selected
      IF prevX EQ X AND prevY EQ Y THEN BEGIN

		X -= 44
		Y -= 44

		; Relative index for orginal arrays from graphed arrays.
		relativeIndex = NINT((num/557.0)*FLOAT(X))

		; So the name isn't so long... :)
		ri = relativeIndex

		; Subset data arrays to selected points.
		altlev   = zalt[Y]
		xlon2    = xlon1[ri]
		xlat2    = xlat1[ri]
		jday1    = jday[X]
		utcTime1 = utcTime[X]

		; Construct fancy time format.
		time = STRCOMPRESS(STRING(utcTime1), /REMOVE_ALL)
		yy = STRMID(utcTime1,0,2)
		mm = STRMID(utcTime1,2,2)
		dd = STRMID(utcTime1,4,2)
		time = yy + ':' + mm + ':' + dd

        val = '            Lat: ' + STRCOMPRESS(STRING(xlat1[ri]), /REMOVE_ALL) +  $
              '            Lon: ' + STRCOMPRESS(STRING(xlon1[ri]), /REMOVE_ALL) +  $
              '            Alt: ' + STRCOMPRESS(STRING(zalt[Y]), /REMOVE_ALL) +  $
              '            UTC (HH:MM:SS): ' + STRCOMPRESS(STRING(time), /REMOVE_ALL) +  $
              '            Potential Temperature: ' + STRCOMPRESS(STRING(theta[Y,ri]), /REMOVE_ALL)

        WIDGET_CONTROL, id, SET_VALUE=val

      ; Event if more than one point is selected.
      ; **NOTE** Most comments are in first "box" because all 4 cases are identical
      ;    except for starting and ending points.
      END ELSE BEGIN

        ; Box: "Top Left Press -> Bottom Right Release"
        IF prevX LT X AND prevY GT Y THEN BEGIN

          ; Relative indexes from selected graph points to original arrays.
          relativeIndexBegin = NINT((num/557.0)*FLOAT(prevX-44))
		  relativeIndexEnd = NINT((num/557.0)*FLOAT(X-44))

		  ; So the names aren't so long... :)
	      riB = relativeIndexBegin
		  riE = relativeIndexEnd

		  ; Indexing numbers.
          k = 0
          m = 0

		  ; First For Loop to populate altitude array (Every third altitude).
          FOR i = Y, prevY, 3 DO BEGIN

			; "i-44" to make up for the black area on the screen.
			; Makes the bottom right of the graph (0,0) instead of (44,44)
            altlev[k] = zalt[i-44]

			; Second For Loop to populate lat, long, day, time, & poten temp arrays.
 		    FOR j = riB, riE DO BEGIN

			  ; Populate arrays with selected points.
			  xlon2[m] = xlon1[j]
			  xlat2[m] = xlat1[j]
			  jday1[m] = jday[j]
			  utcTime1[m] = utcTime[j]
			  potenTemp[m] = theta[i-44,j]

			  ; Construct fancy date format.
			  time = STRCOMPRESS(STRING(utcTime1[m]), /REMOVE_ALL)
			  yy = STRMID(utcTime1[m],0,2)
			  mm = STRMID(utcTime1[m],2,2)
			  dd = STRMID(utcTime1[m],4,2)
			  time = yy + ':' + mm + ':' + dd

			  ; Ouput to the points textbox beneath the drawing area in toolkit.
              val = '            Lat: ' + STRCOMPRESS(STRING(xlat1[j]), /REMOVE_ALL) +  $
                    '            Lon: ' + STRCOMPRESS(STRING(xlon1[j]), /REMOVE_ALL) +  $
                    '            Alt: ' + STRCOMPRESS(STRING(zalt[i-44]), /REMOVE_ALL) +  $
                    '            UTC (HH:MM:SS): ' + STRCOMPRESS(STRING(time), /REMOVE_ALL) +  $
                    '            Potential Temperature: ' + STRCOMPRESS(STRING(theta[i-44,j]), /REMOVE_ALL)

		      IF i EQ Y AND j EQ riB THEN BEGIN
			    WIDGET_CONTROL, id, SET_VALUE=val
			  ENDIF ELSE BEGIN
				WIDGET_CONTROL, id, SET_VALUE=val, /APPEND
			  END

			  m += 1
 			ENDFOR

            k += 1
 		  ENDFOR

		  ; Subset arrays to correct size
		  altlev    = altlev(0:k-1)
		  xlon2     = xlon2(0:m-1)
		  xlat2     = xlat2(0:m-1)
		  jday1     = jday1(0:m-1)
		  utcTime1  = utcTime1(0:m-1)
		  potenTemp = potenTemp(0:m-1)

		  ; Sort arrays and output the range of values:
		  altlevS    = SORT(altlev)
		  xlon2S     = SORT(xlon2)
		  xlat2S     = SORT(xlat2)
		  utcTime1S  = SORT(utcTime1)
		  potenTempS  = SORT(potenTemp)

		  val1 = '             ' + $
		         'Lat. Range                       '    + $
		         'Lon. Range                       ' + $
		         'Alt. Range                        '     + $
		         'UTC Range                 ' + $
		         'Poten. Temp. Range'

		  val2 = '      ' + $
		  		 '{' + STRCOMPRESS(STRING(xlat2[xlat2S[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(xlat2[xlat2S[m-1]]), /REMOVE_ALL) +'}' +  $
				 '        ' + $
				 '{' + STRCOMPRESS(STRING(xlon2[xlon2S[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(xlon2[xlon2S[m-1]]), /REMOVE_ALL) +'}' +  $
				 '           ' + $
				 '{' + STRCOMPRESS(STRING(altlev[altlevS[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(altlev[altlevS[k-1]]), /REMOVE_ALL) +'}' + $
                 '           ' + $
				 '{' + STRCOMPRESS(STRING(utcTime1[utcTime1S[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(utcTime1[utcTime1S[m-1]]), /REMOVE_ALL) +'}' + $
				 '           ' + $
				 '{' + STRCOMPRESS(STRING(potenTemp[potenTempS[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(potenTemp[potenTempS[m-1]]), /REMOVE_ALL) +'}'

		  WIDGET_CONTROL, id, SET_VALUE='', /APPEND
		  WIDGET_CONTROL, id, SET_VALUE=val1, /APPEND
		  WIDGET_CONTROL, id, SET_VALUE=val2, /APPEND
        END

        ; Box: "Bottom Left Press -> Top Right Release"
        IF prevX LT X AND prevY LT Y THEN BEGIN

          ; Relative indexes from selected graph points to original arrays.
          relativeIndexBegin = NINT((num/557.0)*FLOAT(prevX-44))
		  relativeIndexEnd = NINT((num/557.0)*FLOAT(X-44))

		  ; So the names aren't so long... :)
	      riB = relativeIndexBegin
		  riE = relativeIndexEnd

		  ; Indexing numbers.
          k = 0
          m = 0

          FOR i = prevY, Y, 3 DO BEGIN

            altlev[k] = zalt[i-44]

 		    FOR j = riB, riE DO BEGIN

			  xlon2[m] = xlon1[j]
			  xlat2[m] = xlat1[j]
			  jday1[m] = jday[j]
			  utcTime1[m] = utcTime[j]
			  potenTemp[m] = theta[i-44,j]

			  ; Construct fancy date format.
			  time = STRCOMPRESS(STRING(utcTime1[m]), /REMOVE_ALL)
			  yy = STRMID(utcTime1[m],0,2)
			  mm = STRMID(utcTime1[m],2,2)
			  dd = STRMID(utcTime1[m],4,2)
			  time = yy + ':' + mm + ':' + dd

              val = '            Lat: ' + STRCOMPRESS(STRING(xlat1[j]), /REMOVE_ALL) +  $
                    '            Lon: ' + STRCOMPRESS(STRING(xlon1[j]), /REMOVE_ALL) +  $
                    '            Alt: ' + STRCOMPRESS(STRING(zalt[i-44]), /REMOVE_ALL) +  $
                    '            UTC (HH:MM:SS): ' + STRCOMPRESS(STRING(time), /REMOVE_ALL) +  $
                    '            Potential Temperature: ' + STRCOMPRESS(STRING(theta[i-44,j]), /REMOVE_ALL)

		      IF i EQ Y AND j EQ riB THEN BEGIN
			    WIDGET_CONTROL, id, SET_VALUE=val
			  ENDIF ELSE BEGIN
				WIDGET_CONTROL, id, SET_VALUE=val, /APPEND
			  END

			  m += 1
 			ENDFOR

            k += 1
 		  ENDFOR

		  ; Subset arrays to correct size
		  altlev    = altlev(0:k-1)
		  xlon2     = xlon2(0:m-1)
		  xlat2     = xlat2(0:m-1)
		  jday1     = jday1(0:m-1)
		  utcTime1  = utcTime1(0:m-1)
		  potenTemp = potenTemp(0:m-1)

		  ; Sort arrays and output the range of values:
		  altlevS    = SORT(altlev)
		  xlon2S     = SORT(xlon2)
		  xlat2S     = SORT(xlat2)
		  utcTime1S  = SORT(utcTime1)
		  potenTempS  = SORT(potenTemp)

		  val1 = '             ' + $
		         'Lat. Range                       '    + $
		         'Lon. Range                       ' + $
		         'Alt. Range                        '     + $
		         'UTC Range                 ' + $
		         'Poten. Temp. Range'

		  val2 = '      ' + $
		  		 '{' + STRCOMPRESS(STRING(xlat2[xlat2S[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(xlat2[xlat2S[m-1]]), /REMOVE_ALL) +'}' +  $
				 '        ' + $
				 '{' + STRCOMPRESS(STRING(xlon2[xlon2S[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(xlon2[xlon2S[m-1]]), /REMOVE_ALL) +'}' +  $
				 '           ' + $
				 '{' + STRCOMPRESS(STRING(altlev[altlevS[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(altlev[altlevS[k-1]]), /REMOVE_ALL) +'}' + $
                 '           ' + $
				 '{' + STRCOMPRESS(STRING(utcTime1[utcTime1S[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(utcTime1[utcTime1S[m-1]]), /REMOVE_ALL) +'}' + $
				 '           ' + $
				 '{' + STRCOMPRESS(STRING(potenTemp[potenTempS[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(potenTemp[potenTempS[m-1]]), /REMOVE_ALL) +'}'

		  WIDGET_CONTROL, id, SET_VALUE='', /APPEND
		  WIDGET_CONTROL, id, SET_VALUE=val1, /APPEND
		  WIDGET_CONTROL, id, SET_VALUE=val2, /APPEND
        END

        ; Box: "Top Right Press -> Bottom Left Release"
        IF prevX GT X AND prevY GT Y THEN BEGIN

          ; Relative indexes from selected graph points to original arrays.
          relativeIndexBegin = NINT((num/557.0)*FLOAT(X-44))
		  relativeIndexEnd = NINT((num/557.0)*FLOAT(prevX-44))

		  ; So the names aren't so long... :)
	      riB = relativeIndexBegin
		  riE = relativeIndexEnd

		  ; Indexing numbers.
          k = 0
          m = 0

          FOR i = Y, prevY, 3 DO BEGIN

            altlev[k] = zalt[i-44]

 		    FOR j = riB, riE DO BEGIN

			  xlon2[m] = xlon1[j]
			  xlat2[m] = xlat1[j]
			  jday1[m] = jday[j]
			  utcTime1[m] = utcTime[j]
			  potenTemp[m] = theta[i-44,j]

			  ; Construct fancy date format.
			  time = STRCOMPRESS(STRING(utcTime1[m]), /REMOVE_ALL)
			  yy = STRMID(utcTime1[m],0,2)
			  mm = STRMID(utcTime1[m],2,2)
			  dd = STRMID(utcTime1[m],4,2)
			  time = yy + ':' + mm + ':' + dd

              val = '            Lat: ' + STRCOMPRESS(STRING(xlat1[j]), /REMOVE_ALL) +  $
                    '            Lon: ' + STRCOMPRESS(STRING(xlon1[j]), /REMOVE_ALL) +  $
                    '            Alt: ' + STRCOMPRESS(STRING(zalt[i-44]), /REMOVE_ALL) +  $
                    '            UTC (HH:MM:SS): ' + STRCOMPRESS(STRING(time), /REMOVE_ALL) +  $
                    '            Potential Temperature: ' + STRCOMPRESS(STRING(theta[i-44,j]), /REMOVE_ALL)

		      IF i EQ Y AND j EQ riB THEN BEGIN
			    WIDGET_CONTROL, id, SET_VALUE=val
			  ENDIF ELSE BEGIN
				WIDGET_CONTROL, id, SET_VALUE=val, /APPEND
			  END

			  m += 1
 			ENDFOR

            k += 1
 		  ENDFOR

		  ; Subset arrays to correct size
		  altlev    = altlev(0:k-1)
		  xlon2     = xlon2(0:m-1)
		  xlat2     = xlat2(0:m-1)
		  jday1     = jday1(0:m-1)
		  utcTime1  = utcTime1(0:m-1)
		  potenTemp = potenTemp(0:m-1)

		  ; Sort arrays and output the range of values:
		  altlevS    = SORT(altlev)
		  xlon2S     = SORT(xlon2)
		  xlat2S     = SORT(xlat2)
		  utcTime1S  = SORT(utcTime1)
		  potenTempS  = SORT(potenTemp)

		  val1 = '             ' + $
		         'Lat. Range                       '    + $
		         'Lon. Range                       ' + $
		         'Alt. Range                        '     + $
		         'UTC Range                 ' + $
		         'Poten. Temp. Range'

		  val2 = '      ' + $
		  		 '{' + STRCOMPRESS(STRING(xlat2[xlat2S[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(xlat2[xlat2S[m-1]]), /REMOVE_ALL) +'}' +  $
				 '        ' + $
				 '{' + STRCOMPRESS(STRING(xlon2[xlon2S[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(xlon2[xlon2S[m-1]]), /REMOVE_ALL) +'}' +  $
				 '           ' + $
				 '{' + STRCOMPRESS(STRING(altlev[altlevS[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(altlev[altlevS[k-1]]), /REMOVE_ALL) +'}' + $
                 '           ' + $
				 '{' + STRCOMPRESS(STRING(utcTime1[utcTime1S[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(utcTime1[utcTime1S[m-1]]), /REMOVE_ALL) +'}' + $
				 '           ' + $
				 '{' + STRCOMPRESS(STRING(potenTemp[potenTempS[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(potenTemp[potenTempS[m-1]]), /REMOVE_ALL) +'}'

		  WIDGET_CONTROL, id, SET_VALUE='', /APPEND
		  WIDGET_CONTROL, id, SET_VALUE=val1, /APPEND
		  WIDGET_CONTROL, id, SET_VALUE=val2, /APPEND
        END

        ; Box: "Bottom Right Press -> Top Left Release"
        IF prevX GT X AND prevY LT Y THEN BEGIN

          ; Relative indexes from selected graph points to original arrays.
          relativeIndexBegin = NINT((num/557.0)*FLOAT(X-44))
		  relativeIndexEnd = NINT((num/557.0)*FLOAT(prevX-44))

		  ; So the names aren't so long... :)
	      riB = relativeIndexBegin
		  riE = relativeIndexEnd

		  ; Indexing numbers.
          k = 0
          m = 0

          FOR i = prevY, Y, 3 DO BEGIN

            altlev[k] = zalt[i-44]

 		    FOR j = riB, riE DO BEGIN

			  xlon2[m] = xlon1[j]
			  xlat2[m] = xlat1[j]
			  jday1[m] = jday[j]
			  utcTime1[m] = utcTime[j]
			  potenTemp[m] = theta[i-44,j]

			  ; Construct fancy date format.
			  time = STRCOMPRESS(STRING(utcTime1[m]), /REMOVE_ALL)
			  yy = STRMID(utcTime1[m],0,2)
			  mm = STRMID(utcTime1[m],2,2)
			  dd = STRMID(utcTime1[m],4,2)
			  time = yy + ':' + mm + ':' + dd

              val = '            Lat: ' + STRCOMPRESS(STRING(xlat1[j]), /REMOVE_ALL) +  $
                    '            Lon: ' + STRCOMPRESS(STRING(xlon1[j]), /REMOVE_ALL) +  $
                    '            Alt: ' + STRCOMPRESS(STRING(zalt[i-44]), /REMOVE_ALL) +  $
                    '            UTC (HH:MM:SS): ' + STRCOMPRESS(STRING(time), /REMOVE_ALL) +  $
                    '            Potential Temperature: ' + STRCOMPRESS(STRING(theta[i-44,j]), /REMOVE_ALL)

		      IF i EQ Y AND j EQ riB THEN BEGIN
			    WIDGET_CONTROL, id, SET_VALUE=val
			  ENDIF ELSE BEGIN
				WIDGET_CONTROL, id, SET_VALUE=val, /APPEND
			  END

			  m += 1
 			ENDFOR

            k += 1
 		  ENDFOR

		  ; Subset arrays to correct size
		  altlev    = altlev(0:k-1)
		  xlon2     = xlon2(0:m-1)
		  xlat2     = xlat2(0:m-1)
		  jday1     = jday1(0:m-1)
		  utcTime1  = utcTime1(0:m-1)
		  potenTemp = potenTemp(0:m-1)

		  ; Sort arrays and output the range of values:
		  altlevS    = SORT(altlev)
		  xlon2S     = SORT(xlon2)
		  xlat2S     = SORT(xlat2)
		  utcTime1S  = SORT(utcTime1)
		  potenTempS  = SORT(potenTemp)

		  val1 = '             ' + $
		         'Lat. Range                       '    + $
		         'Lon. Range                       ' + $
		         'Alt. Range                        '     + $
		         'UTC Range                 ' + $
		         'Poten. Temp. Range'

		  val2 = '      ' + $
		  		 '{' + STRCOMPRESS(STRING(xlat2[xlat2S[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(xlat2[xlat2S[m-1]]), /REMOVE_ALL) +'}' +  $
				 '        ' + $
				 '{' + STRCOMPRESS(STRING(xlon2[xlon2S[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(xlon2[xlon2S[m-1]]), /REMOVE_ALL) +'}' +  $
				 '           ' + $
				 '{' + STRCOMPRESS(STRING(altlev[altlevS[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(altlev[altlevS[k-1]]), /REMOVE_ALL) +'}' + $
                 '           ' + $
				 '{' + STRCOMPRESS(STRING(utcTime1[utcTime1S[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(utcTime1[utcTime1S[m-1]]), /REMOVE_ALL) +'}' + $
				 '           ' + $
				 '{' + STRCOMPRESS(STRING(potenTemp[potenTempS[0]]), /REMOVE_ALL) + ', ' + STRCOMPRESS(STRING(potenTemp[potenTempS[m-1]]), /REMOVE_ALL) +'}'

		  WIDGET_CONTROL, id, SET_VALUE='', /APPEND
		  WIDGET_CONTROL, id, SET_VALUE=val1, /APPEND
		  WIDGET_CONTROL, id, SET_VALUE=val2, /APPEND
        END
      END
    END
  END

END