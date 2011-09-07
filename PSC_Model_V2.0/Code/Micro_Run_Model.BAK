; Caleb Wherry
; Micro_Run_Model.pro
; Created: October 8th, 2008
;
; Last Modified: October 27th, 2008
;

PRO Micro_Run_Model, ev

  ; Get name od input file
    id=WIDGET_INFO(ev.top, FIND_BY_UNAME='inpFileNameTB')
    WIDGET_CONTROL, id, GET_VALUE=file1
  file1=STRCOMPRESS(STRING(file1))

  ; Take off the extension of the file.
  file3=file1
  file3=STRSPLIT(file3, '\.', /EXTRACT, /REGEX)
  file3=STRCOMPRESS(STRING(file3[0]), /REMOVE_ALL)

  file2=''
  file2=file1

  ; Spawn command to run the microphysical model.
  cmd='.\pscmodel.exe '+file1+' '+file3
  SPAWN,cmd

    id=WIDGET_INFO(ev.top, FIND_BY_UNAME='sizeDistroGraphsInvisL')
    WIDGET_CONTROL, id, GET_VALUE=val
  IF val EQ 1 THEN Micro_Size_Distros_Graphs, file1

    id=WIDGET_INFO(ev.top, FIND_BY_UNAME='historyTotalsGraphsInvisL')
    WIDGET_CONTROL, id, GET_VALUE=val
  IF val EQ 1 THEN Micro_History_Totals_Graphs, file2
END