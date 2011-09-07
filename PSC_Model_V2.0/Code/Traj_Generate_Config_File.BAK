; Caleb Wherry
; Traj_Generate_Config_File.pro
; Created: October 8th, 2008
;
; Last Modified: October 27th, 2008
;

PRO Traj_Generate_Config_File, ev

  ; COMMON Block to share data between Traj_Generate_Points.pro
  ;   and Traj_Generate_Dat_File.pro
  ; ---------------------------------------------
  COMMON SHARE2, altlev,   $
                 xlon2,    $
                 xlat2,    $
                 jday1,    $
                 utcTime1
  ; ---------------------------------------------

  ; Get date.
    id=WIDGET_INFO(ev.top, FIND_BY_UNAME='dateTB')
  WIDGET_CONTROL, id, GET_VALUE=date

  ; Parse date.
  year = STRMID(date,2,2)
  month = STRMID(date,4,2)
  day = STRMID(date,6,2)

  ; How many days to run the trajectory model.
    id=WIDGET_INFO(ev.top, FIND_BY_UNAME='trajNumofDaysTB')
    WIDGET_CONTROL, id, GET_VALUE=num_days
  num_days = FIX(num_days)

  ; See if the run is a backwards or forwards trajectory run.
  ;---------------------------------------------------------
    id=WIDGET_INFO(ev.top, FIND_BY_UNAME='trajTimeStepTB')
    WIDGET_CONTROL, id, GET_VALUE=val
  val = FIX(val)

  IF val LT 0 THEN BEGIN
    utc = N_ELEMENTS(utcTime1)
    utc -= 1
  END
  IF val GT 0 THEN BEGIN
    utc = 0
  END
  ;---------------------------------------------------------


  ; Calculate Julian Day
  ;---------------------------------------------------------
  jday = 0
  mday = [31,28,31,30,31,30,31,31,30,31,30,31]

  IF month GT 1 THEN BEGIN
	  jday = TOTAL(mday(0:month-2)) + day
	ENDIF ELSE BEGIN
	  jday = day
	ENDELSE

	jday = FIX(jday)
  ;---------------------------------------------------------


  ; Calculations for metadata files: Very Long...
  ;--------------------------------------------------------------------------------------------
  ;--------------------------------------------------------------------------------------------
  ;--------------------------------------------------------------------------------------------
  hour = STRMID(utcTime1[utc], 0, 2)
  hour = FIX(hour)

  startTime1 = ''
  startTime2 = ''
  stopTime = ''

  ; Set start time to the correct metadata time.
  IF val LT 0 THEN BEGIN
    IF hour GE 0 AND hour Lt 6 THEN BEGIN
      startTime1 = date + '06'
      startTime2 = date + '_0600'

	    IF num_days LT FIX(day) THEN BEGIN
	      stopTime = STRING(FIX(date) - num_days)
	      stopTime = stopTime + '_0600'
	    END

      IF num_days GE FIX(day) THEN BEGIN
        month -= 1

        IF month EQ 0 THEN BEGIN
          year -= 1
          month = 12
          day = 31 - (num_days - FIX(day))
          stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0600', /REMOVE_ALL)
        END

				IF month EQ 1 OR month EQ 3 OR month EQ 5 OR month EQ 7 OR month EQ 9 OR month EQ 10 THEN BEGIN
				  day = 31 - (num_days - FIX(day))
				  stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0600', /REMOVE_ALL)
				END

				IF month EQ 4 OR month EQ 6 OR month EQ 8 OR month EQ 11 THEN BEGIN
				  day = 30 - (num_days - FIX(day))
				  stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0600', /REMOVE_ALL)
				END

				IF month EQ 2 THEN BEGIN
				  day = 28 - (num_days - FIX(day))
				  stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0600', /REMOVE_ALL)
				END
      END
    END

    IF hour GE 6 AND hour Lt 12 THEN BEGIN
      startTime1 = date + '12'
      startTime2 = date + '_1200'

	    IF num_days LT FIX(day) THEN BEGIN
	      stopTime = STRING(FIX(date) - num_days)
	      stopTime = stopTime + '_1200'
	    END

      IF num_days GE FIX(day) THEN BEGIN
        month -= 1

        IF month EQ 0 THEN BEGIN
          year -= 1
          month = 12
          day = 31 - (num_days - FIX(day))
          stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1200', /REMOVE_ALL)
        END

				IF month EQ 1 OR month EQ 3 OR month EQ 5 OR month EQ 7 OR month EQ 9 OR month EQ 10 THEN BEGIN
				  day = 31 - (num_days - FIX(day))
				  stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1200', /REMOVE_ALL)
				END

				IF month EQ 4 OR month EQ 6 OR month EQ 8 OR month EQ 11 THEN BEGIN
				  day = 30 - (num_days - FIX(day))
				  stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1200', /REMOVE_ALL)
				END

				IF month EQ 2 THEN BEGIN
				  day = 28 - (num_days - FIX(day))
				  stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1200', /REMOVE_ALL)
				END
      END
    END

    IF hour GE 12 AND hour Lt 18 THEN BEGIN
      startTime1 = date + '18'
      startTime2 = date + '_1800'

	    IF num_days LT FIX(day) THEN BEGIN
	      stopTime = STRING(FIX(date) - num_days)
	      stopTime = stopTime + '_1800'
	    END

      IF num_days GE FIX(day) THEN BEGIN
        month -= 1

        IF month EQ 0 THEN BEGIN
          year -= 1
          month = 12
          day = 31 - (num_days - FIX(day))
          stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1800', /REMOVE_ALL)
        END

				IF month EQ 1 OR month EQ 3 OR month EQ 5 OR month EQ 7 OR month EQ 9 OR month EQ 10 THEN BEGIN
				  day = 31 - (num_days - FIX(day))
				  stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1800', /REMOVE_ALL)
				END

				IF month EQ 4 OR month EQ 6 OR month EQ 8 OR month EQ 11 THEN BEGIN
				  day = 30 - (num_days - FIX(day))
				  stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1800', /REMOVE_ALL)
				END

				IF month EQ 2 THEN BEGIN
				  day = 28 - (num_days - FIX(day))
				  stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1800', /REMOVE_ALL)
				END
      END
    END

    IF hour GE 18 THEN BEGIN
      startTime1 = date + '00'
      startTime2 = date + '_0000'

	    IF num_days LT FIX(day) THEN BEGIN
	      stopTime = STRING(FIX(date) - num_days)
	      stopTime = stopTime + '_0000'
	    END

      IF num_days GE FIX(day) THEN BEGIN
        month -= 1

        IF month EQ 0 THEN BEGIN
          year -= 1
          month = 12
          day = 31 - (num_days - FIX(day))
          stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0000', /REMOVE_ALL)
        END

				IF month EQ 1 OR month EQ 3 OR month EQ 5 OR month EQ 7 OR month EQ 9 OR month EQ 10 THEN BEGIN
				  day = 31 - (num_days - FIX(day))
				  stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0000', /REMOVE_ALL)
				END

				IF month EQ 4 OR month EQ 6 OR month EQ 8 OR month EQ 11 THEN BEGIN
				  day = 30 - (num_days - FIX(day))
				  stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0000', /REMOVE_ALL)
				END

				IF month EQ 2 THEN BEGIN
				  day = 28 - (num_days - FIX(day))
				  stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0000', /REMOVE_ALL)
				END
      END
    END
  END

  IF val GT 0 THEN BEGIN
    IF hour GE 0 AND hour Lt 6 THEN BEGIN
      startTime1 = date + '00'
      startTime2 = date + '_0000'

		  IF month EQ 1 OR month EQ 3 OR month EQ 5 OR month EQ 7 OR month EQ 9 OR month EQ 10  OR month EQ 12 THEN BEGIN
		    IF (FIX(day)+num_days) LE 31 THEN BEGIN
		      day = FIX(day) + num_days
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0000', /REMOVE_ALL)
		    END

		    IF (FIX(day)+num_days) GT 31 THEN BEGIN
		      month += 1
		      day = (FIX(day) + num_days) - 31
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0000', /REMOVE_ALL)
		    END
			END

			IF month EQ 4 OR month EQ 6 OR month EQ 8 OR month EQ 11 THEN BEGIN
		    IF (FIX(day)+num_days) LE 30 THEN BEGIN
		      day = FIX(day) + num_days
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0000', /REMOVE_ALL)
		    END

		    IF (FIX(day)+num_days) GT 30 THEN BEGIN
		    	month += 1
		      day = (FIX(day) + num_days) - 30
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0000', /REMOVE_ALL)
		    END
			END

			IF month EQ 2 THEN BEGIN
		    IF (FIX(day)+num_days) LE 28 THEN BEGIN
		      day = FIX(day) + num_days
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0000', /REMOVE_ALL)
		    END

		    IF (FIX(day)+num_days) GT 28 THEN BEGIN
		    	month += 1
		      day = (FIX(day) + num_days) - 28
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0000', /REMOVE_ALL)
		    END
			END

			IF month EQ 13 THEN BEGIN
				year += 1
				month = 1
				day = (FIX(day) + num_days) - 31
				stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0000', /REMOVE_ALL)
			END
    END

    IF hour GE 6 AND hour Lt 12 THEN BEGIN
      startTime1 = date + '06'
      startTime2 = date + '_0600'

		  IF month EQ 1 OR month EQ 3 OR month EQ 5 OR month EQ 7 OR month EQ 9 OR month EQ 10  OR month EQ 12 THEN BEGIN
		    IF (FIX(day)+num_days) LE 31 THEN BEGIN
		      day = FIX(day) + num_days
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0600', /REMOVE_ALL)
		    END

		    IF (FIX(day)+num_days) GT 31 THEN BEGIN
		      month += 1
		      day = (FIX(day) + num_days) - 31
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0600', /REMOVE_ALL)
		    END
			END

			IF month EQ 4 OR month EQ 6 OR month EQ 8 OR month EQ 11 THEN BEGIN
		    IF (FIX(day)+num_days) LE 30 THEN BEGIN
		      day = FIX(day) + num_days
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0600', /REMOVE_ALL)
		    END

		    IF (FIX(day)+num_days) GT 30 THEN BEGIN
		    	month += 1
		      day = (FIX(day) + num_days) - 30
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0600', /REMOVE_ALL)
		    END
			END

			IF month EQ 2 THEN BEGIN
		    IF (FIX(day)+num_days) LE 28 THEN BEGIN
		      day = FIX(day) + num_days
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0600', /REMOVE_ALL)
		    END

		    IF (FIX(day)+num_days) GT 28 THEN BEGIN
		    	month += 1
		      day = (FIX(day) + num_days) - 28
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0600', /REMOVE_ALL)
		    END
			END

			IF month EQ 13 THEN BEGIN
				year += 1
				month = 1
				day = (FIX(day) + num_days) - 31
				stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_0600', /REMOVE_ALL)
			END
    END

    IF hour GE 12 AND hour Lt 18 THEN BEGIN
      startTime1 = date + '12'
      startTime2 = date + '_1200'

		  IF month EQ 1 OR month EQ 3 OR month EQ 5 OR month EQ 7 OR month EQ 9 OR month EQ 10  OR month EQ 12 THEN BEGIN
		    IF (FIX(day)+num_days) LE 31 THEN BEGIN
		      day = FIX(day) + num_days
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1200', /REMOVE_ALL)
		    END

		    IF (FIX(day)+num_days) GT 31 THEN BEGIN
		      month += 1
		      day = (FIX(day) + num_days) - 31
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1200', /REMOVE_ALL)
		    END
			END

			IF month EQ 4 OR month EQ 6 OR month EQ 8 OR month EQ 11 THEN BEGIN
		    IF (FIX(day)+num_days) LE 30 THEN BEGIN
		      day = FIX(day) + num_days
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1200', /REMOVE_ALL)
		    END

		    IF (FIX(day)+num_days) GT 30 THEN BEGIN
		    	month += 1
		      day = (FIX(day) + num_days) - 30
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1200', /REMOVE_ALL)
		    END
			END

			IF month EQ 2 THEN BEGIN
		    IF (FIX(day)+num_days) LE 28 THEN BEGIN
		      day = FIX(day) + num_days
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1200', /REMOVE_ALL)
		    END

		    IF (FIX(day)+num_days) GT 28 THEN BEGIN
		    	month += 1
		      day = (FIX(day) + num_days) - 28
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1200', /REMOVE_ALL)
		    END
			END

			IF month EQ 13 THEN BEGIN
				year += 1
				month = 1
				day = (FIX(day) + num_days) - 31
				stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1200', /REMOVE_ALL)
			END
    END

    IF hour GE 18 THEN BEGIN
      startTime1 = date + '18'
      startTime2 = date + '_1800'

		  IF month EQ 1 OR month EQ 3 OR month EQ 5 OR month EQ 7 OR month EQ 9 OR month EQ 10  OR month EQ 12 THEN BEGIN
		    IF (FIX(day)+num_days) LE 31 THEN BEGIN
		      day = FIX(day) + num_days
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1800', /REMOVE_ALL)
		    END

		    IF (FIX(day)+num_days) GT 31 THEN BEGIN
		      month += 1
		      day = (FIX(day) + num_days) - 31
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1800', /REMOVE_ALL)
		    END
			END

			IF month EQ 4 OR month EQ 6 OR month EQ 8 OR month EQ 11 THEN BEGIN
		    IF (FIX(day)+num_days) LE 30 THEN BEGIN
		      day = FIX(day) + num_days
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1800', /REMOVE_ALL)
		    END

		    IF (FIX(day)+num_days) GT 30 THEN BEGIN
		    	month += 1
		      day = (FIX(day) + num_days) - 30
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1800', /REMOVE_ALL)
		    END
			END

			IF month EQ 2 THEN BEGIN
		    IF (FIX(day)+num_days) LE 28 THEN BEGIN
		      day = FIX(day) + num_days
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1800', /REMOVE_ALL)
		    END

		    IF (FIX(day)+num_days) GT 28 THEN BEGIN
		    	month += 1
		      day = (FIX(day) + num_days) - 28
		      stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1800', /REMOVE_ALL)
		    END
			END

			IF month EQ 13 THEN BEGIN
				year += 1
				month = 1
				day = (FIX(day) + num_days) - 31
				stoptime = STRCOMPRESS(STRING(year)+STRING(month)+STRING(day)+'_1800', /REMOVE_ALL)
			END
    END
  END
  ;--------------------------------------------------------------------------------------------
  ;--------------------------------------------------------------------------------------------
  ;--------------------------------------------------------------------------------------------


  ; Fancy way of getting the directory above the PWD
  ;------------------------------------------------------------------
  filename = FILE_SEARCH('PSC_Model_ToolKit.prj',/FULLY_QUALIFY_PATH)
  filename = STRSPLIT(filename, '\', /EXTRACT)
  num = N_ELEMENTS(filename)
  file=''

  FOR c = 0, num-3 DO BEGIN
    file = file+filename[c]+'\'
  ENDFOR
  ;------------------------------------------------------------------

  path = file+'Trajectories\Input\'
  filename = 'GEOS5_'+STRING(year)+STRING(month)+STRCOMPRESS(STRING(day), /REMOVE_ALL)+'_'+utcTime1[utc]
  ext = '.config'
  file = path+filename+ext

  ; Open and write to file.
  OPENW,1, file

    PRINTF,1, "'Trajectories with GMAO4 analysis'"
    PRINTF,1, '4'
    PRINTF,1, '100,300,500,700'
    PRINTF,1, '.false.'
    PRINTF,1, '.true.'
    PRINTF,1, '.FALSE.'
    PRINTF,1, '.FALSE.'
    PRINTF,1, "'rfill.traj'"
    PRINTF,1, '2004022800'
    PRINTF,1, "'/private/var/automount/nfs/artisan/CALIPSO/GMAO/GEOS_510/D5OTVDYN/2007/'"
    PRINTF,1, '29'

    geos5FileStart = STRING(jday) + '/Processed/DAS.ops.asm.tavg3ddynv.GEOS510.' + STRCOMPRESS(startTime2, /REMOVE_ALL) + '.V01.nc3'
    ;geos5fileStop =
    PRINTF,1, "'" + STRCOMPRESS(STRING(geos5FileStart), /REMOVE_ALL) + "'"

    PRINTF,1, STRCOMPRESS(startTime1, /REMOVE_ALL)
    PRINTF,1, '6.'

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='trajTimeStepTB')
      WIDGET_CONTROL, id, GET_VALUE=val
      val = STRCOMPRESS(STRING(val), /REMOVE_ALL)
    PRINTF,1, val+'.'

    PRINTF,1, '.FALSE.'
    PRINTF,1, '.FALSE.'
    PRINTF,1, '21'
    PRINTF,1, '10.'
    PRINTF,1, '0.1'
    PRINTF,1, '0.1'
    PRINTF,1, '.FALSE.'
    PRINTF,1, '1'
    PRINTF,1, '.25'
    PRINTF,1, '-3.'
    PRINTF,1, '-1.'
    PRINTF,1, '.FALSE'
    PRINTF,1, '1'
    PRINTF,1, '.5'
    PRINTF,1, '195.'
    PRINTF,1, '-2.'
    PRINTF,1, '.FALSE.'
    PRINTF,1, '0.005'
    PRINTF,1, '0.005'
    PRINTF,1, '.FALSE'
    PRINTF,1, "'/usr38/users/fairlie/UW-Hybrid/Traject/'"
    PRINTF,1, "'henry_traj_20040228.dat'"
    PRINTF,1, '.false.'
    PRINTF,1, '65.'
    PRINTF,1, '212.'
    PRINTF,1, '2700.'
    PRINTF,1, '1'
    PRINTF,1, '1440'
    PRINTF,1, '.FALSE.'
    PRINTF,1, "'/usr38/users/fairlie/SOLVE/DIAL_data/Ozone_data/'"
    PRINTF,1, '1'
    PRINTF,1, "'dial_o3_20000315.diag'"
    PRINTF,1, '.FALSE.'
    PRINTF,1, "'/usr38/users/fairlie/SOLVE/DIAL_data/Ozone_data/'"
    PRINTF,1, '1'
    PRINTF,1, "'dial_o3_20000315.diag'"
    PRINTF,1, '.TRUE.'
    PRINTF,1, "'/private/var/automount/nfs/artisan/Users/mpitts/Trajectory/Orbits/'"
    PRINTF,1, '1'
    PRINTF,1, filename+'.dat'
    PRINTF,1, '.FALSE.'
    PRINTF,1, '1'
    PRINTF,1, '170.'
    PRINTF,1, '-45.'
    PRINTF,1, '135.'
    PRINTF,1, '-80.'
    PRINTF,1, '.FALSE.'
    PRINTF,1, "'/usr38/users/harvey/DIAL_data/Datfiles/'"
    PRINTF,1, '1'
    PRINTF,1, "'dial_o3_20000315.diag.dc8'"
    PRINTF,1, '1'
    PRINTF,1, "'/private/var/automount/nfs/artisan/Users/mpitts/Trajectory/Output/GEOS5_traject_20070526_06Z_bkwd_sts2.traj'"

  CLOSE,1

END