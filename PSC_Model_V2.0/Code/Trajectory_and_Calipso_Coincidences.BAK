; Mike Pitts
; Trajectory_and_Calipso_Coincidences.pro
;
;

PRO Trajectory_and_Calipso_Coincidences

  ; Earth polar radius
	Re = 6357.

	; Read in trajectory data from here. Time, latitude, longitude, theta, etc.
	ifilet = 'GEOS5_traject_20070526_06Z_bkwd.traj'

	fdir = 'C:\Documents and Settings\Mike\Projects\CALIPSO Related\PSCs\Trajectory\Output\'
	indir = 'C:\Documents and Settings\Mike\Projects\CALIPSO Related\PSCs\Trajectory\Output\'
	outdir = indir
	outlf=''

	; In meters
	Alt_init=20655.

	OPENR, 1, fdir+ifilet, /F77_UNFORMATTED, /SWAP_ENDIAN
	  PRINT, ifilet

	  nmax=600000l
	  nrmax=181
	  ncmax=361
	  x=FLTARR(ncmax)
	  y=FLTARR(nrmax)
	  charexp='                                   '
	  nthp=0L
	  ukmo=' '
	  nmc=' '
	  ecmwf=' '
	  restart=' '
	  rfile='                                     '
	  nrday=0L
	  dir='                                       '
	  nfile=0L
	  wfile='                                               '
	  istime=0l
	  ictime=0l
	  dtflow=0.
	  dt=0.
	  igw=' '
	  stream=' '
	  nstrm=0L
	  ds1=0.
	  strm0=0.
	  dstrm=0.
	  pviso=' '
	  npviso=0L
	  ds2=0.
	  pv0=0.
	  dpv=0.
	  tiso=' '
	  ntiso=0L
	  ds3=0.
	  temp0=0.
	  dtemp=0.
	  space=' '
	  dxs=0.
	  dys=0.
	  henry=' '
	  dirh='                                     '
	  hfile='                                     '
	  range_fill=' '
	  stlat=0.
	  stlon=0.
	  er2_range=' '
	  nrings=0L
	  npts=0L
	  dial=' '
	  dird='                                     '
	  ndial=0l
	  dfile='                                     '
	  ; New option "curtain" - tdf 4/19/2K5
	  curtain=' '
	  dirc='                                     '
	  ncurt=0l
	  cfile='                                     '
	  ; New option "calipso" - tdf 8/22/2K5
	  calipso=' '
	  er2=' '
	  nleg=0L
	  xleg0=0.
	  yleg0=0.
	  xleg1=0.
	  yleg1=0.
	  dc8=' '
	  dirdc8='/usr73/SOLVEII/DIAL_data/Datfiles/'
	  ndc8=0L
	  dc8file='dial_o3_20000309.diag.dc8'
	  dtout=0.
	  ofile='                                          '
	  nr=0L
	  nc=0L
	  time=0.
	  ntraj=0l
	  char1='     '
	  char2='     '
	  char3='     '
	  char4='     '
	  char5='      '
	  char6='       '
	  char10='       '
	  char11='o3traj'
	  char12='      '
	  char13='qtraj'
	  char14='t0_traj'
	  char15='age_traj'
	  char16 = 'alt_traj'

	  ;
	  ; Read input parameters for this experiment
	  ;

	  READU,1,charexp
	  READU,1,nthp       ; these are now eta surfaces
	  thlevp=FLTARR(nthp)
	  READU,1,thlevp
	  READU,1,ukmo
	  READU,1,nmc
	  READU,1,ecmwf
	  READU,1,restart
	  READU,1,rfile
	  READU,1,nrday
	  READU,1,dir

	  READU,1,nfile
	  wfiles = STRARR(nfile)
	  FOR n=0, nfile-1 DO BEGIN
	    READU,1,wfile
	    wfiles(n) = STRTRIM(wfile)
	    PRINT,' wfile =', wfile
	  ENDFOR

	  READU,1,istime
	  READU,1,dtflow
	  READU,1,dt
	  READU,1,igw
	  READU,1,stream
	  READU,1,nstrm
	  READU,1,ds1
	  READU,1,strm0
	  READU,1,dstrm
	  READU,1,pviso
	  READU,1,npviso
	  READU,1,ds2
	  READU,1,pv0
	  READU,1,dpv
	  READU,1,tiso
	  READU,1,ntiso
	  READU,1,ds3
	  READU,1,temp0
	  READU,1,dtemp
	  READU,1,space
	  READU,1,dxs
	  READU,1,dys
	  ; Added new parameter  tdf
	  READU,1,henry
	  READU,1,dirh
	  READU,1,hfile
	  ; Added new parameter
	  READU,1,range_fill
	  READU,1,stlat
	  READU,1,stlon
	  READU,1,er2_range
	  READU,1,nrings
	  READU,1,npts
	  READU,1,dial
	  READU,1,dird

	  READU,1,ndial
	  FOR n=1,ndial DO BEGIN
	    READU,1,dfile
	  ENDFOR

	  ; New option "curtain" - tdf 4/19/2K5
	  READU,1,curtain
	  READU,1,dirc

	  READU,1,ncurt
	  FOR n=1,ncurt DO BEGIN
	    READU,1,cfile
	  ENDFOR

	  ; New option "calipso" - tdf 8/22/2K5
	  READU,1,calipso
	  READU,1,dirc

	  READU,1,ncurt
	  FOR n=1,ncurt DO BEGIN
	    READU,1,cfile
	  ENDFOR

	  READU,1,er2

	  READU,1,nleg
	  FOR n=1,nleg DO BEGIN
	    READU,1,xleg0
	    READU,1,yleg0
	    READU,1,xleg1
	    READU,1,yleg1
	  ENDFOR

	  READU,1,dc8
	  READU,1,dirdc8

	  READU,1,ndc8
	  FOR n=1,ndc8 DO BEGIN
	    READU,1,dc8file
	  ENDFOR

	  READU,1,dtout
	  READU,1,ofile

	  PRINT,'    Header read OK!'


	  ; Read trajectory output.
	  nfiles = 7L*24L+1L
	  ;nfiles=49l
	  ;***** Problem: Sometimes no trajectories are output.
	  ;               Solve this by keeping track of the number
	  ;               of outputs.
	  ;                                      tdf 5/10/2K5

	  times=FLTARR(nfiles)
	  traj_time = FLTARR(nfiles)
	  ctimes = LONARR(nfiles)
	  tlat = FLTARR(nfiles)
	  tlon = tlat
	  ttemp = tlat
	  talt = tlat
	  ttime = tlat
	  ttheta = tlat

	  ; Loop over files.
	  FOR n=0,nfiles-1 DO BEGIN

	    IF (EOF(1)) THEN GOTO, fend

	    READU,1,istime,ictime,time,ntraj
	    PRINT,n,istime,ictime,time,ntraj

	    stime=STRCOMPRESS(STRING(istime),/remove_all)
	    ctime=STRCOMPRESS(STRING(ictime),/remove_all)
	    mdate=STRMID(stime,0,8)

	    times(n) = time
	    ctimes(n) = ictime - 2007000000L
	    dfrac = FLOAT(ctimes(n))/100.
	    dfrac2 = (FLOAT(ctimes(n)) - FIX(dfrac)*100.)/24.
	    traj_time(n) = FIX(dfrac) + dfrac2

	    xn=FLTARR(ntraj)
	    yn=FLTARR(ntraj)
	    ; zn = eta level
	    zn=FLTARR(ntraj)
	    pvn=FLTARR(ntraj)
	    pn=FLTARR(ntraj)
	    ; thn = potential temperature        tdf 11/08/2K4
	    thn=FLTARR(ntraj)
	    qdfn=FLTARR(ntraj)
	    ; Montgomery S.F.                    tdf 11/08/2K4
	    o3n=FLTARR(ntraj)
	    ; shn = specific humidity            tdf 11/08/2K4
	    shn=FLTARR(ntraj)
	    t0n=FLTARR(ntraj)
	    qn=FLTARR(ntraj)
	    agen = FLTARR(ntraj)
	    altn=FLTARR(ntraj)

	    READU,1,char1
	    READU,1,xn
	    READU,1,char2
	    READU,1,yn
	    ; pn  = pressure
	    READU,1,char3
	    READU,1,pn
	    READU,1,char4
	    READU,1,pvn
	    ; thn = potential temperature        tdf 11/08/2K4
	    READU,1,char5
	    READU,1,thn
	    READU,1,char10
	    READU,1,qdfn
	    READU,1,char11
	    ; O3n = OZONE             tdf 11/08/2K4
	    READU,1,o3n
	    READU,1,char12
	    ; shn = specific humidity            tdf 11/08/2K4
	    READU,1,shn
	    READU,1,char13
	    READU,1,qn
	    READU,1,char14
	    READU,1,t0n
	    READU,1,char15
	    READU,1,agen
	    READU,1,char16
	    READU,1,altn

	    tmp = (1013./pn)^(-0.286)*thn

	    index=WHERE(xn GE 0. AND xn LE 360. AND yn GE -90. AND yn LE 90.)

	    aodn=0.*qdfn

	    IF n EQ 0 THEN BEGIN
	      tidx = INTARR(1)
	      ntraj0 = ntraj
	      ; Initial locations
	      x0=xn
	      y0=yn
	      z0=zn
	      p0=pn
	      pv0=pvn
	      th0=thn
	      alt0=altn
	      tmp0 = tmp
	      traj_alt1 = 23500.
	      traj_alt2 = 24500.
	      ;x1 = WHERE(y0 le -60. and x0 gt 70. and alt0 ge traj_alt1 and alt0 lt traj_alt2)
	      ;tidx(0) = max(x1)

	      ;x1 = WHERE(y0 le -65. and x0 gt 70. and alt0 ge traj_alt1 and alt0 lt traj_alt2)
	      ;tidx(1) = max(x1)

	      x1 = WHERE(y0 LE -80. AND x0 GT -90. AND alt0 GE traj_alt1 AND alt0 LT traj_alt2)
	      tidx(0) = MAX(x1)

	      ;x1 = WHERE(y0 le -65. and x0 gt 70. and alt0 ge traj_alt1 and alt0 lt traj_alt2)
	      ;tidx(0) = max(x1)

	      ; Initial points

	    ENDIF

	    IF ntraj EQ ntraj0 THEN BEGIN
		    tlat(n)= yn(tidx)
		    tlon(n) = xn(tidx)
		    talt(n) = altn(tidx)
		    ttemp(n) = tmp(tidx)
		    ttheta(n) = thn(tidx)
	    ENDIF ELSE BEGIN
		    tlat(n)= y0(tidx)
		    tlon(n) = x0(tidx)
		    talt(n) = alt0(tidx)
		    ttemp(n) = tmp0(tidx)
		   ttheta(n) = th0(tidx)
	    ENDELSE

	  ENDFOR

	CLOSE, 1

	fend:

	PRINT,'Time= ',time

	jumpall:

	ntpts = nfiles

	orbn_hno3 = 'MLS-Aura_L2GP-HNO3_v02-20-c01_2006d127.he5'
	orbn_h2o = 'MLS-Aura_L2GP-H2O_v02-20-c01_2006d127.he5'

	mday = [31,28,31,30,31,30,31,31,30,31,30,31]

	fresult = ' '
	fresult_mls  =''

	path_fname = 'G:\PSCs\Data Files\Feature Files\Antarctic\2007 v3\'

	filename = 'PSC_features_052*_v3.0_ffv3.1.dat'
	fresult = [fresult,FILE_SEARCH(path_fname,filename)]

	filename_mls = 'MLS_v2.2_2007052*_gas.dat'
	fresult_mls = [fresult_mls,FILE_SEARCH(path_fname,filename_mls)]

	nd = N_ELEMENTS(fresult)-1
	nd2 = N_ELEMENTS(fresult_mls)-1

	IF nd NE nd2 THEN STOP, 'Problem with files matching!'

	fresult = fresult(1:nd)
	fresult_mls = fresult_mls(1:nd)

	cal_lat = 0.
	cal_lon = 0.
	cal_time = 0L
	cal_lrat = 0.
	cal_ldepol = 0.
	cal_theta = 0.
	cal_temp = 0.
	cal_perp = 0.
	cal_feat = 0
	cal_hno3 = 0.
	cal_h2o = 0.
	cal_tnat = 0.

	FOR id = 0, nd-1 DO BEGIN
	  month= STRMID(fresult(id),64,2)
	  day = STRMID(fresult(id),66,2)

	  ; Calculate Julian day number
	  IF month GT 1 THEN BEGIN
	     jday = TOTAL(mday(0:month-2)) + day
	  ENDIF ELSE BEGIN
	     jday = day
	  ENDELSE

	  IF month EQ 5 AND day GE 20 AND day LE 26 THEN BEGIN

	    OPENR,2,fresult(id)   ; Open CALIPSO Feature File

	      orbn1 = 'L1AtBkz01.10n060615-030333.hdf'
		    num = 0L
		    nz = num
		    nh = num
		    READU,2,num,nz,nh
		    zalt = FLTARR(nz)
		    xlat = FLTARR(nh)
		    norbs = FLTARR(num)
		    xlon = xlat
		    xtrop = xlat
		    time = DBLARR(nh)
		    feature = INTARR(nz,nh)
		    temp = FLTARR(nz,nh)
		    lrat = temp
		    crat = temp
		    perp = temp
		    beta532 = temp
		    aer532 = temp
		    ldepol = temp
		    theta = temp
		    vorto = temp
		    vortc = temp
		    eqvl = temp

		    norbs(id) = num
	      IF num GT 0 THEN BEGIN  ; Can filter by number of orbits if desired
	        orbname = STRARR(num)
			    orbname(*) = orbn1
			    READU,2,orbname
			    orb_time1 = DBLARR(num)
			    orb_time2 = orb_time1
			    READU,2,orb_time1
			    READU,2,orb_time2
			    rthresh = FLTARR(4,5)
			    pthresh = rthresh
			    rnpts = rthresh
			    pnpts = rthresh
			    nrat = FLTARR(4,4)
			    nperp = nrat

			    READU,2,rthresh
			    READU,2,rnpts
		    	READU,2,pthresh
			    READU,2,pnpts
			    READU,2,nrat
	    		READU,2,nperp
	    		READU,2,zalt
	    		READU,2,xlat
			    READU,2,xlon
	    		READU,2,time
	    		READU,2,feature
	    		READU,2,temp
	    		READU,2,theta
		    	READU,2,lrat
		    	READU,2,crat
		    	READU,2,perp
		    	READU,2,ldepol
		    	READU,2,beta532
		    	READU,2,aer532
		    	READU,2,vorto
		    	READU,2,vortc
		      READU,2,eqvl

			    OPENR,2,fresult_mls(id)  ; Read in Aura MLS data interpolated to CALIPSO grid
	   		    num_mls = 0L
			      np1 = num_mls
			      np2 = num_mls
			      nh_mls = num_mls

			      READU,2,num_mls,nh_mls,np1,np2
			      xlat_mls = FLTARR(nh_mls)
			      xlon_mls = xlat_mls
		 	      time_mls = DBLARR(nh_mls)
			      p_hno3 = FLTARR(np1)
			      p_h2o = FLTARR(np2)
			      mls_hno3 = FLTARR(np1,nh_mls)
			      mls_h2o = FLTARR(np2,nh_mls)
			      mls_tnat = FLTARR(np1,nh_mls)
			      orbname_hno3 = STRARR(num_mls)
			      orbname_h2o = STRARR(num_mls)
	  		    orbname_hno3(*) = orbn_hno3
	  		    orbname_h2o(*) = orbn_h2o

			      READU,2,orbname_hno3
			      READU,2,orbname_h2o
			      READU,2,xlat_mls
			      READU,2,xlon_mls
			      READU,2,time_mls
			      READU,2,p_hno3
			      READU,2,p_h2o
			      READU,2,mls_hno3
			      READU,2,mls_h2o
			      READU,2,mls_tnat

			    CLOSE, 2

			    ; Interpolate MLS data to CALIPSO altitudes.
			    chno3 = FLTARR(nz,nh)
			    ch2o = chno3
			    ctnat = chno3
			    chno3(*,*) = -9999.
			    ch2o(*,*) = -9999.
			    ctnat(*,*) = -9999.

			    FOR ii = 0,nh-1 DO BEGIN
	   			  xobs = WHERE(FINITE(theta(*,ii)) NE 0 AND FINITE(temp(*,ii)) NE 0,nobs)

							; Should be at least 5 points in the profle.
	   			    IF nobs GT 5 THEN BEGIN
	   			      presc = 1013.25*(theta(xobs,ii)/temp(xobs,ii))^(-1./0.286)
	   		   	    xhno3 = WHERE(mls_hno3(*,ii) NE -9999.,nhno3)
	   			      xh2o = WHERE(mls_h2o(*,ii) NE -9999.,nh2o)
	   			      xtnat = WHERE(mls_tnat(*,ii) NE -9999.,ntnat)
	   			      IF nhno3 GT 5 THEN chno3(xobs,ii) = INTERPOL(mls_hno3(xhno3,ii),ALOG(p_hno3(xhno3)),ALOG(presc))
	   			      IF nh2o GT 5 THEN ch2o(xobs,ii) = INTERPOL(mls_h2o(xh2o,ii),ALOG(p_h2o(xh2o)),ALOG(presc))
	   			      IF ntnat GT 5 THEN ctnat(xobs,ii) = INTERPOL(mls_tnat(xtnat,ii),ALOG(p_hno3(xtnat)),ALOG(presc))
	   			    ENDIF
			    ENDFOR

	        ; Subset data for latitudes covered by grid.
	        x1 = WHERE(xlat LE -50. AND xlat GE -90. AND time GT 0.)
	        temp = temp(*,x1)
	        xlat = xlat(x1)
	        xlon = xlon(x1)
	        feat = feature(*,x1)
	        perp = perp(*,x1)
	        lrat = lrat(*,x1)
	        theta = theta(*,x1)
	        ldepol = ldepol(*,x1)
	        ctnat = ctnat(*,x1)
	        chno3 = chno3(*,x1)
	        ch2o = ch2o(*,x1)
	        aerback = aer532(*,x1)
	        time = time(x1)
	        nprofs = N_ELEMENTS(time)
	        temp3 = FLTARR(nprofs)
	        temp3(*)=-9999.
	        lrat3 = temp3
	        ldepol3 = temp3
	        theta3 = temp3
	        perp3 = temp3
	        feat3 = temp3
	        chno33 = temp3
	        ch2o3 = temp3
	        ctnat3 = temp3
	        x2 = WHERE(zalt LE 24.5 AND zalt GE 23.5,nn)
	        temp = temp(x2,*)
	        lrat = lrat(x2,*)
	        ldepol = ldepol(x2,*)
	        theta = theta(x2,*)
	        perp = perp(x2,*)
	        feat = feat(x2,*)
	        chno3 = chno3(x2,*)
	        ch2o = ch2o(x2,*)
	        ctnat = ctnat(x2,*)

	        FOR kk = 0,nprofs-1 DO BEGIN
	          x3 = WHERE(FINITE(temp(*,kk)) NE 0 AND FINITE(lrat(*,kk)) NE 0,nzg)

	          IF nzg GT 0 THEN BEGIN
		          temp3(kk) = TOTAL(temp(x3,kk))/FLOAT(nzg)
	    	      lrat3(kk) = TOTAL(lrat(x3,kk))/FLOAT(nzg)
	        	  ldepol3(kk) = TOTAL(ldepol(x3,kk))/FLOAT(nzg)
	   	    	  theta3(kk) = TOTAL(theta(x3,kk))/FLOAT(nzg)
	   	    	  perp3(kk) = TOTAL(perp(x3,kk))/FLOAT(nzg)
	        	  feat3(kk) = TOTAL(feat(x3,kk))/FLOAT(nzg)
	        	  chno33(kk) = TOTAL(chno3(x3,kk))/FLOAT(nzg)
	        	  ch2o3(kk) = TOTAL(ch2o(x3,kk))/FLOAT(nzg)
	        	  ctnat3(kk) = TOTAL(ctnat(x3,kk))/FLOAT(nzg)
	        	ENDIF
	        ENDFOR

	        cal_temp = [cal_temp,temp3]
	        cal_lrat = [cal_lrat,lrat3]
	        cal_ldepol = [cal_ldepol,ldepol3]
	        cal_theta = [cal_theta,theta3]
	        cal_lat = [cal_lat,xlat]
	        cal_lon = [cal_lon,xlon]
	        cal_time = [cal_time,time]
	        cal_perp = [cal_perp,perp3]
	        cal_feat = [cal_feat,feat3]
	        cal_hno3 = [cal_hno3,chno33]
	        cal_h2o = [cal_h2o,ch2o3]
	        cal_tnat = [cal_tnat,ctnat3]

	      ENDIF
	    CLOSE, 1
	  ENDIF
	ENDFOR

	nprofs = N_ELEMENTS(cal_temp)-1
	cal_temp = [cal_temp(1:nprofs)]
	cal_lrat = [cal_lrat(1:nprofs)]
	cal_ldepol = [cal_ldepol(1:nprofs)]
	cal_theta = [cal_theta(1:nprofs)]
	cal_lat = [cal_lat(1:nprofs)]
	cal_lon = [cal_lon(1:nprofs)]
	cal_time = [cal_time(1:nprofs)]
	cal_perp = [cal_perp(1:nprofs)]
	cal_feat = [cal_feat(1:nprofs)]
	cal_hno3 = [cal_hno3(1:nprofs)]
	cal_h2o = [cal_h2o(1:nprofs)]
	cal_tnat = [cal_tnat(1:nprofs)]
	xg = WHERE(cal_lrat GT -9999.)
	cal_temp = cal_temp(xg)
	cal_lrat = cal_lrat(xg)
	cal_ldepol = cal_ldepol(xg)
	cal_theta = cal_theta(xg)
	cal_lat = cal_lat(xg)
	cal_lon = cal_lon(xg)
	cal_time = cal_time(xg)
	cal_perp = cal_perp(xg)
	cal_feat = cal_feat(xg)
	cal_hno3 = cal_hno3(xg)
	cal_h2o = cal_h2o(xg)
	cal_tnat = cal_tnat(xg)

	xtraj = WHERE(cal_time GE MIN(traj_time) AND cal_time LE MAX(traj_time),npts)
	cal_lat2 = cal_lat(xtraj)
	cal_lon2 = cal_lon(xtraj)
	cal_temp2 = cal_temp(xtraj)
	cal_lrat2 = cal_lrat(xtraj)
	cal_ldepol2 = cal_ldepol(xtraj)
	cal_time2 = cal_time(xtraj)
	cal_feat2 = cal_feat(xtraj)
	cal_hno32 = cal_hno3(xtraj)
	cal_h2o2 = cal_h2o(xtraj)
	cal_tnat2 = cal_tnat(xtraj)
	cal_theta2 = cal_theta(xtraj)

	IF npts GT 0 THEN BEGIN

	  OPENW, 3, 'C:\Documents and Settings\Mike\Projects\CALIPSO Related\PSCs\Trajectory\Output\traj_20070526_23km_80s_nat.dat'

	    d = FLTARR(npts)
	    td = d
	    tempd = d
	    r_coin = FLTARR(ntpts)
	    r_coin(*) = -9999.
	    p_coin = r_coin
	    t_coin = r_coin
	    !x.style=1
	    !y.style=1
	    LOADCT,39
	    PLOT,1/r_coin,p_coin,xtitle='1/R',ytitle='Aer. Dep.',yrange=[-.2,1.2],xrange=[1.4,0],/nodata

	    FOR jj = 0,ntpts-1L DO BEGIN
	      FOR kk = 0L,npts-1L DO BEGIN
	        d(kk) = Re*ACOS((SIN(cal_lat2(kk)*!pi/180.)*SIN(tlat(jj)*!pi/180.))+$
	                   COS(cal_lat2(kk)*!pi/180.)*COS(tlat(jj)*!pi/180.)*$
	                   COS((cal_lon2(kk)-tlon(jj))*!pi/180.))
	        td(kk) = (cal_time2(kk) - traj_time(jj))*24.
	        tempd(kk) = cal_temp2(kk)-ttemp(jj)
	      ENDFOR

	      xxx = WHERE(d LT 250. AND ABS(td) LT 5. AND ABS(tempd) LT 5.,nnn)
	      PRINTF, 3,traj_time(jj),tlat(jj),tlon(jj),ttemp(jj),ttheta(jj),nnn

	      IF nnn GT 0 THEN BEGIN
	        rat_med = TOTAL(cal_lrat2(xxx))/FLOAT(nnn)
	        dep_med = TOTAL(cal_ldepol2(xxx))/FLOAT(nnn)
	        OPLOT,1./cal_lrat2(xxx),cal_ldepol2(xxx),psym=2

	        FOR jk = 0,nnn-1 DO BEGIN
	          PRINTF,3,zalt(x2(0)),cal_lat2(xxx(jk)),cal_lon2(xxx(jk)),cal_time2(xxx(jk)),$
	                    cal_lrat2(xxx(jk)),cal_ldepol2(xxx(jk)),cal_temp2(xxx(jk)),cal_theta2(xxx(jk)),$
	                    cal_feat2(xxx(jk)),$
	                    d(xxx(jk)),td(xxx(jk)),tempd(xxx(jk)),cal_hno32(xxx(jk)),cal_h2o2(xxx(jk)),cal_tnat2(xxx(jk)),$
	                    format = '(8e18.6,i8,6e18.6)'
	        ENDFOR
	      ENDIF
		  ENDFOR
	  CLOSE, 3
	ENDIF
	; Want to pick the CALIPSO points that are closest in both time and space
END