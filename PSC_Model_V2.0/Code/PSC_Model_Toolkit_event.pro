; Caleb Wherry
; PSC_Model_Toolkit_event.pro
; IDL GUI Event Handler
; Created: October 8th, 2008
;
; Last Modified: October 27th, 2008
;


; Events from main procedure.
PRO PSC_Model_Toolkit_event, ev

  wTarget = (WIDGET_INFO(ev.id,/NAME) eq 'TREE' ?  $
      WIDGET_INFO(ev.id, /tree_root) : ev.id)

  wWidget =  ev.top

  CASE wTarget OF

    ; Microphysical Model generate button event handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='generateB'): BEGIN
       Micro_Generate_Inp_File, ev
    END

    ; Microphysical Model run button event handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='runB'): BEGIN
       Micro_Run_Model, ev
    END

    ; Microphysical Model open input button even handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='openInputB'): BEGIN
       Micro_Open_Inp_File, ev
    END

    ; Microphysical Model slider event handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='homogeNATFactoSB'): BEGIN
       value = ev.value*.01
       id=WIDGET_INFO(ev.top, FIND_BY_UNAME='homogeNATFactoL')
       WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(value), /REMOVE_ALL)
    END

    ; Microphysical Model radio button 1 event handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='tempCalcRB1'): BEGIN
       id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempCalcInvisL')
       WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(0), /REMOVE_ALL)

       id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempPressHisFileTB')
       WIDGET_CONTROL, id, SET_VALUE='---'
    END

    ; Microphysical Model radio button 2 event handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='tempCalcRB2'): BEGIN
       id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempCalcInvisL')
       WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(1), /REMOVE_ALL)

       id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempPressHisFileTB')
       WIDGET_CONTROL, id, SET_VALUE=''
       WIDGET_CONTROL, id, /EDITABLE
    END

    ; Microphysical Model radio button 3 event handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='tempCalcRB3'): BEGIN
       id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempCalcInvisL')
       WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(2), /REMOVE_ALL)

       id=WIDGET_INFO(ev.top, FIND_BY_UNAME='tempPressHisFileTB')
       WIDGET_CONTROL, id, SET_VALUE=''
       WIDGET_CONTROL, id, /EDITABLE
    END

    ; Microphysical Model check box 1 event handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='sizeDistroGraphsCB'): BEGIN
       id=WIDGET_INFO(ev.top, FIND_BY_UNAME='sizeDistroGraphsInvisL')
       WIDGET_CONTROL, id, GET_VALUE=val
       IF val EQ 1 THEN WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(0), /REMOVE_ALL)
       IF val EQ 0 THEN WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(1), /REMOVE_ALL)
    END

    ; Microphysical Model check box 2 event handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='historyTotalsGraphsCB'): BEGIN
       id=WIDGET_INFO(ev.top, FIND_BY_UNAME='historyTotalsGraphsInvisL')
       WIDGET_CONTROL, id, GET_VALUE=val
       IF val EQ 1 THEN WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(0), /REMOVE_ALL)
       IF val EQ 0 THEN WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(1), /REMOVE_ALL)
    END

    ; Trajectory Model radio button 1 event handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='directionRB1'): BEGIN
       id=WIDGET_INFO(ev.top, FIND_BY_UNAME='trajTimeStepTB')
       WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(-15), /REMOVE_ALL)
    END

    ; Trajectory Model radio button 2 event handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='directionRB2'): BEGIN
       id=WIDGET_INFO(ev.top, FIND_BY_UNAME='trajTimeStepTB')
       WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(15), /REMOVE_ALL)
    END

    ; Trajectory Model create images button even handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='createImagesB'): BEGIN
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='whichImageInvisL')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(0), /REMOVE_ALL)

      Traj_Graph_Features_Data, ev
    END

    ; Trajectory Model previous image button even handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='previousImageB'): BEGIN
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='whichImageInvisL')
      WIDGET_CONTROL, id, GET_VALUE=val
      val = FIX(val)

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='whichImageInvisL')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val-1), /REMOVE_ALL)

      Traj_Graph_Features_Data, ev
    END

    ; Trajectory Model next image button even handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='nextImageB'): BEGIN
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='whichImageInvisL')
      WIDGET_CONTROL, id, GET_VALUE=val
      val = FIX(val)

      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='whichImageInvisL')
      WIDGET_CONTROL, id, SET_VALUE=STRCOMPRESS(STRING(val+1), /REMOVE_ALL)

      Traj_Graph_Features_Data, ev
    END


    ; Trajectory Model DrawBox events.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='imageDrawbox'): BEGIN
      Traj_Generate_Points, ev
    END

    ; Trajectory Model clearPointsB event handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='clearPointsB'): BEGIN
      id=WIDGET_INFO(ev.top, FIND_BY_UNAME='pointsTB')
      WIDGET_CONTROL, id, SET_VALUE=''
    END

    ; Trajectory Model generate button even handler.
    WIDGET_INFO(wWidget, FIND_BY_UNAME='generate1B'): BEGIN
      Traj_Generate_Dat_File, ev
      Traj_Generate_Config_File, ev
    END

    ELSE: ;PRINT, ' Whoops! Unhandled event...'
  ENDCASE

END