; Mike Pitts
; Modified by: Caleb Wherry
; Traj_Graph_Features_Data.pro
; Procedure to plot *.dat file.
; Created: October 8th, 2008
;
; Last Modified: October 24th, 2008
;

; Rounds A to a long integer.
FUNCTION NINT, a

  ab = ABS(a)
  i  = WHERE(ab EQ 0,count)
  IF(count NE 0) THEN ab(i) = 1

  RETURN,LONG(a+0.5*(a/ab))
END


; Outputs graphs onto the window DRAW widget.
PRO Traj_Graph_Features_Data, ev

  ; COMMON Block to share data between procedures
  ; ---------------------------------------------
  COMMON SHARE1, zalt,       $
  				 xlon1,      $
				 xlat1,      $
				 time1,      $
				 jday,       $
				 utcTime,    $
				 theta
  ; ---------------------------------------------


  ; Load color table and color scale.
  ;--------------------------------------------
  DEVICE, DECOMPOSED = 0
  R = BYTARR(6)
  G = BYTARR(6)
  B = BYTARR(6)

  R = [0,255,255,0,0,255]
  G = [0,0,255,255, 0,255]
  B = [0,0,0,0,255,255]
  TVLCT, R, G, B
  LOADCT, 39
  ;--------------------------------------------


  ; Gat date.
    id=WIDGET_INFO(ev.top, FIND_BY_UNAME='dateTB')
    WIDGET_CONTROL, id, GET_VALUE=date
  date = STRCOMPRESS(STRING(date), /REMOVE_ALL)

  ; Error if no date is entered into the GUI.
  IF date EQ '' THEN BEGIN
    ans = DIALOG_MESSAGE('Please Enter A Date!', TITLE='Error: No Date!', /CENTER, /ERROR)
    GOTO, JUMP3  ; Line 260
  END

  ; Loading label...
  id=WIDGET_INFO(ev.top, FIND_BY_UNAME='loadingL')
  WIDGET_CONTROL, id, SET_VALUE='Loading...'

  ; Parse date.
  year = STRMID(date,0,4)
  month = STRMID(date,4,2)
  day = STRMID(date,6,2)

	;**********************************************
  ; Change To Represent where feature files are.
  ;**********************************************
  path = 'C:\PSCs\Data Files\Feature Files\Antarctic\'+year+' v3\'

  filename = 'PSC_features_'+month+day+'_v3.0_ffv3.1.dat'
  file = path+filename

  OPENR,1, file
    orbn1 = 'L1AtBkz01.10n060615-030333.hdf'
    num = 0L                     ; Number of Orbits
    nz = 0L                      ; Number of altitudes
    nh = 0L                      ; Number of X values(Lat/Lon)
    READU,1, num, nz, nh

    orbname = STRARR(num)
    FOR i=0,num-1 DO BEGIN
      orbname(i) = orbn1
    ENDFOR

    orb_time1 = DBLARR(num)      ; Orbit Start Times
    orb_time2 = DBLARR(num)      ; Orbit Stop Times
    rthresh   = FLTARR(4,5)
    pthresh   = FLTARR(4,5)
    rnpts     = FLTARR(4,5)
    pnpts     = FLTARR(4,5)
    nrat      = FLTARR(4,4)
    nperp     = FLTARR(4,4)
    zalt      = FLTARR(nz)		 ; Altitudes
    xlat      = FLTARR(nh)       ; Latitudes
    xlon      = FLTARR(nh) 		 ; Longitudes
    time      = DBLARR(nh)	     ; Times
    feature   = INTARR(nz,nh)    ; Features
    temp      = FLTARR(nz,nh)	 ; Temperature
    theta     = FLTARR(nz,nh)	 ; Potential Temperature
    lrat      = FLTARR(nz,nh)
    crat      = FLTARR(nz,nh)
    perp      = FLTARR(nz,nh)
    ldepol    = FLTARR(nz,nh)
    beta532   = FLTARR(nz,nh)
    aer532    = FLTARR(nz,nh)
    vorto     = FLTARR(nz,nh)
    vortc     = FLTARR(nz,nh)
    eqvl      = FLTARR(nz,nh)
    xtrop     = FLTARR(nh)
    xtropw    = FLTARR(nh)

    READU,1, orbname
    READU,1, orb_time1
    READU,1, orb_time2
    READU,1, rthresh
    READU,1, rnpts
    READU,1, pthresh
    READU,1, pnpts
    READU,1, nrat
    READU,1, nperp
    READU,1, zalt
    READU,1, xlat
    READU,1, xlon
    READU,1, time
    READU,1, feature
    READU,1, temp
    READU,1, theta
    READU,1, lrat
    READU,1, crat
    READU,1, perp
    READU,1, ldepol
    READU,1, beta532
    READU,1, aer532
    READU,1, vorto
    READU,1, vortc
    READU,1, eqvl
    READU,1, xtrop
    READU,1, xtropw
  CLOSE, 1

  ; Get the current orbit number that is displayed.
    id=WIDGET_INFO(ev.top, FIND_BY_UNAME='whichImageInvisL')
    WIDGET_CONTROL, id, GET_VALUE=orbit_num
  orbit_num = FIX(orbit_num)

  ; Error if the orbit number is less than 0.
  IF orbit_num LT 0 THEN BEGIN
    id=WIDGET_INFO(ev.top, FIND_BY_UNAME='whichImageInvisL')
    WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(0), /REMOVE_ALL)

    ans = DIALOG_MESSAGE('No More Images! Please Go Forward!', TITLE='Error: No More Images!', /CENTER, /ERROR)
    GOTO, JUMP3  ; Line 260
  END

  ; Error if the orbit number is greater then the given number of orbits
  ;   from the file read in.
  IF orbit_num GT num-1 THEN BEGIN
    id=WIDGET_INFO(ev.top, FIND_BY_UNAME='whichImageInvisL')
    WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(num-1), /REMOVE_ALL)

    ans = DIALOG_MESSAGE('No More Images! Please Go Backwards!', TITLE='Error: No More Images!', /CENTER, /ERROR)
    GOTO, JUMP3  ; Line 260
  END

  ;
  ; Plot the current graph on the screen according to the current orbit number.
  ;
  i = orbit_num

  ; Change longitude into a -180 to 180 degree E/W convention instead of full 360 degree convention.
  xx = WHERE(xlon GT 180,npts1)
  IF npts1 GT 0 THEN xlon(xx) = xlon(xx) - 360.

  ; Reverse altitude array to make beginning index lowest altitude.
  zalt    = REVERSE(zalt)

  ; Feature data is revered later in code? Unsure where. Do not reverese here.
  ;feature = REVERSE(feature,1)

  temp    = REVERSE(temp,1)
  theta   = REVERSE(theta,1)
  lrat    = REVERSE(lrat,1)
  crat    = REVERSE(crat,1)
  perp    = REVERSE(perp,1)
  ldepol  = REVERSE(ldepol,1)
  beta532 = REVERSE(beta532,1)
  aer532  = REVERSE(aer532,1)
  vorto   = REVERSE(vorto,1)
  vortc   = REVERSE(vortc,1)
  eqvl    = REVERSE(eqvl,1)

  ; Subset data into the time and altitudes we want.
  mz = WHERE(zalt GE 8 AND zalt LE 32,nz2)
  mm = WHERE(time GE orb_time1(i) AND time LE orb_time2(i) AND xlat NE -9999,npt)

  ; Subset xlon, xlat, and time to the above times, latitudes, and altitudes.
  IF npt GT 0 THEN BEGIN
    xlon1 = xlon(mm)
    xlat1 = xlat(mm)
    time1 = time(mm)
    theta = theta(*,mm)
  ENDIF

  ; Calculate UTC Time
  ;-------------------------------------------------------------
  ;-------------------------------------------------------------
  utcTime = STRARR(npt)
  jday = DBLARR(npt)
  month_day = FIX(month+day)

  FOR j = 0,npt-1 DO BEGIN
    j_hours  = DOUBLE((time1(j)-month_day)*24.)
    jday(j) = j_hours
    h1 = FIX(j_hours)
    j_minutes  = DOUBLE((j_hours-h1)*60.)
    m1 = FIX(j_minutes)
    j_seconds  = DOUBLE((j_minutes-m1)*60.)
    s1 = NINT(j_seconds)

    h1 = STRCOMPRESS(STRING(h1, FORMAT='(I02)'), /REMOVE_ALL)
    m1 = STRCOMPRESS(STRING(m1, FORMAT='(I02)'), /REMOVE_ALL)
    s1 =  STRCOMPRESS(STRING(s1, FORMAT='(I02)'), /REMOVE_ALL)
    utcTime(j) = h1+m1+s1
  ENDFOR
  ;--------------------------------------------------------------
  ;--------------------------------------------------------------


  ; Plotting Routine
  ;----------------------------------------------------------------
  ;----------------------------------------------------------------

  ; Populate feature2 2D array with correct number of elements.
  feature2 = INTARR(nz2,npt)
  FOR iz = 0,nz2-1 DO BEGIN
    FOR ih = 0,npt-1 DO BEGIN
      feature2(iz,ih) = feature(mz(iz),mm(ih))
    ENDFOR
  ENDFOR

  mlvls = [0,1,5,15]
  feature3 = TRANSPOSE(feature2)
  newa = BYTARR(npt,nz2)
  newa(*,*) = 0
  feet3 = feature3 MOD 100

  xx = WHERE(feet3 EQ 1 OR feet3 EQ 2 ,n1)
  IF n1 GT 0 THEN newa(xx) = 250

  xx = WHERE(feet3 EQ 3 OR feet3 EQ 4 ,n2)
  IF n2 GT 0 THEN newa(xx) = 180

  xx = WHERE(feet3 EQ 9 OR feet3 EQ 10,n3)
  IF n3 GT 0 THEN newa(xx) = 100

  xx = WHERE(feet3 EQ 27 OR feet3 EQ 28,n4)
  IF n4 GT 0 THEN newa(xx) = 40


  ; Set up IDL plot parameters.
  !ORDER      = 1
  !X.TICKLEN  = -0.02
  !X.TICKS    = 11

  ;Full size
  grid_position   = [0.065,.22,.9,.83]
  bar_position    = [.93,.22,.95,.83]
  xtitle          = '(Latitude; Longitude)'
  ytitle          = 'Altitude, km'
  title           = orbname(i)
  labsize         = 0.9
  !Y.RANGE    = [8,30]
  !X.RANGE    = [0,npt-1+100]
  !X.STYLE    = 1
  !Y.STYLE    = 1
  !P.POSITION = grid_position

  xlab1 = STRING(xlat(mm(0)),format='(F6.2)')+','+STRING(xlon(mm(0)),format='(F7.2)')
  xlab2 = STRING(xlat(mm(npt/2)),format='(F6.2)')+','+STRING(xlon(mm(npt/2)),format='(F7.2)')
  xlab3 = STRING(xlat(mm(npt-1)),format='(F6.2)')+','+STRING(xlon(mm(npt-1)),format='(F7.2)')

  PLOT,[0],xtitle=xtitle,ytitle=ytitle,charsize=labsize,title=title,xticks=[2],xtickname=[xlab1,xlab2,xlab3]

  ; Compute the device coordinates for registering the image on the grid.
  py0 = !P.CLIP(1)
  py0 = LONG(NINT(py0))+1
  py1 = !P.CLIP(1)+ nz*(!P.CLIP(3)-!P.CLIP(1))/nz
  py1 = LONG(NINT(py1))
  px0 = !P.CLIP(0)+1
  px0 = LONG(NINT(px0))
  px1 = !P.CLIP(0)+npt*(!P.CLIP(2)-!P.CLIP(0))/npt
  px1 = LONG(NINT(px1))

  EXPAND, newa, px1-px0, py1-py0, b

  ; Plot the graph on the screen.
  TVSCL,b,px0,py0,xsize=px1-px0,ysize=py1-py0

  PLOTS,[!P.CLIP(0),!P.CLIP(2)],[!P.CLIP(1),!P.CLIP(1)],/DEVICE,thick=1.5
  PLOTS,[!P.CLIP(2),!P.CLIP(2)],[!P.CLIP(1),!P.CLIP(3)],/DEVICE,thick=1.5
  PLOTS,[!P.CLIP(2),!P.CLIP(0)],[!P.CLIP(3),!P.CLIP(3)],/DEVICE,thick=1.5
  PLOTS,[!P.CLIP(0),!P.CLIP(0)],[!P.CLIP(1),!P.CLIP(3)],/DEVICE,thick=1.5

  ; Assign values that will appear in the legend's color scale.
  leg_min   = 0.
  leg_max   = 1.
  leg_data  = [.125,0.375,0.625,.875]

  ; Compute color bar constants.
  range = leg_max - leg_min

  !ORDER = 0
  !X.TICKS = 0
  !Y.TICKS = 0
  bar = 1 + FINDGEN(4)
  bar = REFORM(bar,[1,4])
  bar(0,0) = 250
  bar(0,1) = 170
  bar(0,2) = 100
  bar(0,3) = 40

  ; bar = 255-bar
  !Y.RANGE    = [leg_min,leg_max]
  !X.RANGE    = [0,1]
  !X.STYLE    = 5
  !Y.STYLE    = 5
  !P.POSITION = bar_position
  PLOT,[0],/nodata,/noerase

  ;Expand to fill the grid if non-scalable pixel devices.
  EXPAND,bar,!P.CLIP(2)-!P.CLIP(0),!P.CLIP(3)-!P.CLIP(1),bbar

  py0 = !P.CLIP(1)
  py0 = LONG(nint(py0))+1
  py1 = !P.CLIP(1)+ 3*(!P.CLIP(3)-!P.CLIP(1))/3
  py1 = LONG(nint(py1))

  px0 = !P.CLIP(0)+1
  px0 = LONG(nint(px0))
  px1 = !P.CLIP(0)+(!P.CLIP(2)-!P.CLIP(0))
  px1 = LONG(nint(px1))

  TVSCL, bbar, px0, py0, xsize=px1-px0, ysize=py1-py0

  ;Put a line border on the bar.
  PLOTS,[0,0],[leg_min,leg_max],thick=1
  PLOTS,[1,1],[leg_min,leg_max],thick=1
  PLOTS,[0,1],[leg_min,leg_min],thick=1
  PLOTS,[0,1],[leg_max,leg_max],thick=1

  ;Place the legend.
  leg_num=['5','15','45','135']
  XYOUTS,1.3,leg_data(0),leg_num(0),charsize=labsize
  PLOTS,[1,1.],[leg_data(0),leg_data(0)],linestyle=0
  XYOUTS,1.3,leg_data(1),leg_num(1),charsize=labsize
  XYOUTS,1.3,leg_data(2),leg_num(2),charsize=labsize
  PLOTS,[1,1.],[leg_data(2),leg_data(2)],linestyle=0
  XYOUTS,1.3,leg_data(3),leg_num(3),charsize=labsize
  PLOTS,[1,1.],[leg_data(3),leg_data(3)],linestyle=0

  XYOUTS,3,0.5,'Feature Averaging Scale (km)',charsize=labsize,orientation=90.,align=0.5

  ;--------------------------------------------------------------------------------------
  ;--------------------------------------------------------------------------------------
  ; End Plotting Routine

  ; Jump if no date is entered into the GUI.
  JUMP3:

  ; Set Loading label back to blank.
  id=WIDGET_INFO(ev.top, FIND_BY_UNAME='loadingL')
  WIDGET_CONTROL, id, SET_VALUE=''
END