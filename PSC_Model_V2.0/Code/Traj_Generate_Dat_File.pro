; Caleb Wherry
; Traj_Generate_Dat_File.pro
; Procedure to create *.dat file.
; Created: October 8th, 2008
;
; Last Modified: October 27th, 2008
;

; Create Date File Procedure
; date is in the format of: yyyymmdd
PRO Traj_Generate_Dat_File, ev

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
  date = STRCOMPRESS(STRING(date), /REMOVE_ALL)

  ; Year with 4 digits: 2006
  year = STRMID(date,0,4)
  ; Year with 2 digits: 06
  year1 = STRMID(date,2,2)
  month = STRMID(date,4,2)
  day = STRMID(date,6,2)

  ; Size of arrays being passed in.
  nz2 = N_ELEMENTS(altlev)
  npt = N_ELEMENTS(xlon2)  ; Also size of xlat2.

  ; Turn everything into strings and take out spaces.
  nz2    = STRCOMPRESS(STRING(nz2), /REMOVE_ALL)
  altlev = STRCOMPRESS(STRING(altlev*1000), /REMOVE_ALL)
  npt    = STRCOMPRESS(STRING(npt), /REMOVE_ALL)
  xlon2  = STRCOMPRESS(STRING(xlon2), /REMOVE_ALL)
  xlat2  = STRCOMPRESS(STRING(xlat2), /REMOVE_ALL)
  jday1  = STRCOMPRESS(STRING(FLOAT(jday1)), /REMOVE_ALL)


  ; See if the run is a backwards or forwards trajectory.
  ;------------------------------------------------------
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
  ;------------------------------------------------------


  ; Finding & constructing path name from one file level up. Kinda tricky...
  ;-------------------------------------------------------------------------
  filename = FILE_SEARCH('PSC_Model_ToolKit.prj',/FULLY_QUALIFY_PATH)
  filename = STRSPLIT(filename, '\', /EXTRACT)
  num = N_ELEMENTS(filename)

  file=''
  FOR c = 0, num-3 DO BEGIN
   	file = file+filename[c]+'\'
  ENDFOR
  ;-------------------------------------------------------------------------


  ; Output date to unformatted file.
  OPENW,1, file+'Trajectories\Input\GEOS5_'+year1+month+day+'_'+utcTime1[utc]+'.dat',/F77_UNFORMATTED,/SWAP_ENDIAN
    WRITEU,1, date
    WRITEU,1, nz2
    WRITEU,1, altlev
    WRITEU,1, npt
    WRITEU,1, xlon2
    WRITEU,1, xlat2
    WRITEU,1, jday1
  CLOSE,1
END