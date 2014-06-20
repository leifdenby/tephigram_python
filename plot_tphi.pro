;+
; NAME:
;	plot_tphi
;
; PURPOSE:
;	Tephigram atmospheric diagram and optionally
; 	The axes are temperature and theta, both in Celsius
;	(theta=potential temperature, actually it should be log(theta)).
; 	The axes are rotated by 45 degrees so the isobars are nearly horizontal.
;
; 	A temperature and dew point sounding may optionally be plotted
;	on the diagram.
;   
; CATEGORY:
;
; CALLING SEQUENCE:
;
;       plot_tphi,pres=pres,Temp=Temp,td=td,col_temp=col_temp,$
;	col_td=col_td,wspeed=wspeed,wangle=wangle,title=title,$
;	linestyle=linestyle,/over
;
; EXAMPLE:
;
; INPUTS:
;
; OPTIONAL INPUT PARAMETERS:
;  pres		Fltarr(nlev): pressures (hPa) of the sounding at nlev different levels
;  temp		Fltarr(nlev): sounding temperatures (deg C)
;  td	  	Fltarr(nlev): sounding dew point temperatures (deg C)
;  col_temp	Byte: color index for temperture line
;  col_td	Byte: color index for dew point temperature line
;  linestyle	Integer: linestyle index for sounding lines (default=0=solid line)
;  wspeed	Fltarr(nlev): wind speed at nlev levels (m/s)
;  wangle	Fltarr(nlev): wind angle at nlev levels (0 deg = wind
;			from the north, clockwise)
;  title	String: Plot title (default is 'Tephigram')
;
; KEYWORD INPUT PARAMETERS:
;  over		Set keyword to plot sounding over existing graph
;
; OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;	
; MODIFICATION HISTORY:
;
; Program in IDL      October 1997
; Frank Evans  University of Colorado    evans@nit.colorado.edu
; Dominik Brunner, KNMI, Feb 2000, converted to a subroutine.
;	The sounding can now be passed as paramters.
;-
PRO plot_tphi,pres=pres,Temp=Temp,td=td,col_temp=col_temp,$
	col_td=col_td,wspeed=wspeed,wangle=wangle,title=title,$
	linestyle=linestyle,over=over

; number of levels of vertical sounding
nlev=n_elements(pres)
; wind data supplied?
windok=(n_elements(wspeed) EQ nlev) AND (n_elements(wangle) EQ nlev) $
	AND (nlev GT 0)

; Define temperature range T1-T2
T1=-10  &  T2=40  &  theta1=-12

; overlay sounding over exising graph
IF keyword_set(over) THEN goto,soundingonly

psactive=!d.flags AND 1B	; check for device with scalable pixel
				; (postscript device)

IF n_elements(title) EQ 1 THEN plottitle=title ELSE plottitle='Tephigram'

; Set character sizes for title and labels
titlesize=1.0
labsize=0.8
legendsize=1.0

;  Set up plot area (reserved space for wind barbs if necessary
IF windok THEN xhi=0.89 ELSE xhi=0.95
xlo=0.15  &  ylo=0.14  &  yhi=0.95
IF (n_elements(wspeed) GT 0) AND (n_elements(wangle) EQ n_elements(wspeed)) $
THEN xhi=0.89 ELSE xhi=0.93
plot, [0], [0], position=[xlo,ylo,xhi,yhi], /nodata, ticklen=0.0, $
    xticks=1, xtickname=[' ',' '], yticks=1, ytickname=[' ',' '], $
    title='!5'+plottitle+'!3', charsize=titlesize, thick=1.5
!x.s = [0.,1.]	; same scale for device and normal coordinates
!y.s = [0.,1.]	; same scale for device and normal coordinates

;  Set up IDL scaling and rotation transformation matrix
scale=(xhi-xlo)*sqrt(0.5)/(T2-T1)
t3d, /reset, translate=-[T1,theta1,0.], scale=[scale,scale,1.]
t3d, rotate=[0,0,-45]
t3d, scale=[1.0,0.75,1.]
t3d, translate=[xlo,ylo,0.]

;  Initialize constants that control where lines drawn and labeled
thetalo=-20.     &  thetahi=110.
Tlo=-90.         &  Thi=40.
Tstart=-70.      &  Tend=+30.         &  Tstep=10.
Tstartlab=-60.   &  Tendlab=+30.      &  thetalabel=30.
thetastart=-30.  &  thetaend=+90.     &  thetastep=10.
thetastartlab=0. &  thetaendlab=+80.  &  Tlabel=-20.
Pstart=1000.     &  Pend=200.         &  Pstep=-100.

; Set character sizes for title and labels

;  Make the temperature lines and label them
for T = Tstart, Tend, Tstep do begin
    plots, [T,T], [thetalo,thetahi], /t3d,  noclip=0, thick=1.5
endfor
for T = Tstartlab, Tendlab, Tstep do begin
    xyouts, T-0.5,thetalabel, /t3d, alignment=-0.3, orient=90, noclip=0, $
    size=150*labsize, charthick=1.5, 'T='+string(T,format='(I3)'), /data
endfor

;  Make the theta (potential temperature) lines and label them
for th = thetastart, thetaend, thetastep do begin
    plots, [Tlo,Thi], [th,th], /t3d,  noclip=0, thick=1.5
endfor
for th = thetastartlab, thetaendlab, thetastep do begin
    xyouts, Tlabel, th+0.5, alignment=-0.1, orient=0, noclip=0,  $
    size=150*labsize, charthick=1.5, /t3d, '!7h!3='+string(th,format='(I2)')
endfor

; Make the curves of constant pressure
; Use theta = T*(1000/p)^0.286
Tarr=Tlo+(Thi-Tlo)*findgen(100)/100.0
for P = Pstart, Pend, Pstep do begin
    thetaarr = -273.15 + theta(P,Tarr)
    plots, Tarr, thetaarr, /t3d, noclip=0, linestyle=2
    x = (P/1000)^0.286
;    delT = (273+T2)*(1-x)/(1+x)
    delT = (273+T1)*(1-x)/(1+x)
;    T = T2 - delT  &    th = T2 + delT
    T = T1 - delT  &    th = T1 + delT
;    xyouts, T, th, alignment=0.1, orient=53.13,  $
;    size=150*labsize, charthick=1.5, /t3d, string(P,format='(I4)')

    xyouts, T-6., th-7.6, alignment=0.1, orient=50.,  $
    size=150*labsize, charthick=1.5, /t3d, string(P,format='(I4)')
endfor
;xyouts, xhi+0.06, ylo-0.01, alignment=0.5,  $
;    size=1.5*labsize, charthick=1.5, '!5p (hPa)!3', /normal
xyouts, xlo-0.1, (ylo+yhi)/2, alignment=0.5,orientation=90,  $
    size=1.5*labsize, charthick=1.5, '!5Pressure (hPa)!3', /normal

;  Make the saturated specific humidity curves
;    q_s = 0.622*e_s/p
qs=[30,20,15,10,7,5,3,2,1.5,1.0,0.7]
Tarr=-40+findgen(75)
for i = 0, n_elements(qs)-1 do begin
    Tk=Tarr+273.15
    Pqs = 0.622*esat(Tk)/(0.001*qs(i))
    thetaarr = -273.15 + theta(Pqs,Tarr)
    plots, Tarr, thetaarr, /t3d, noclip=0, linestyle=1
endfor
Pqs=[1070,1070,1070,1070,1070,1070,1070,1070,970,800,695]
Tarr=[33,26,21,15,9,4,-3,-8,-15,-22,-27]
thetaarr = -273.15 + theta(Pqs,Tarr)
xyouts, Tarr, thetaarr, alignment=0.0, orient=53,  $
    size=120*labsize, charthick=1.5, /t3d, string(qs,format='(F4.1)')
xyouts, (xlo+xhi)/2, ylo-0.05, alignment=1.0,  $
    size=1.5*labsize, charthick=1.5, '!5q!Ds!N (g/kg)!3', /normal

;  Make the saturated adiabat curves
;    Use dtheta = -L*theta/(Cp*T) dqs and integrate from T=-40 C
Tarr=fltarr(70)  &  thetaarr=fltarr(70)
for th = thetastartlab, thetaend, thetastep do begin
   T = -40
   P = 1000*( (T+273.15)/(th+273.15) )^3.5
   Tk=T+273.15
   qs0 = 0.622*esat(Tk)/P
   thetaarr(0) = th
   Tarr(0) = T
   for i = 1, 70-1 do begin
     T = T + 1.0
     P = 1000*( (T+273.15)/(thetaarr(i-1)+273.15) )^3.5
     Tk=T+273.15
     qs = 0.622*esat(Tk)/P
     thetaarr(i) = thetaarr(i-1) - $
          (2.5E6/1004)* (273.15+thetaarr(i-1))/(273.15+T) * (qs-qs0)
     qs0 = qs
     Tarr(i) = T
   endfor
   plots, Tarr, thetaarr, /t3d, noclip=0, linestyle=5
endfor

;!p.charthick=1.5
;!p.charsize=1.2
plots, [0.15,0.20], [0.045,0.045], /normal, linestyle=2
xyouts, 0.21, 0.04, alignment=0.0,  '!5isobars (pressure)',charsize=legendsize,$
	charthick=1.5,/normal
plots, [0.15,0.20], [0.025,0.025], /normal, linestyle=5
xyouts, 0.21, 0.02, alignment=0.0,  '!5saturated adiabats',charsize=legendsize,$
	charthick=1.5,/normal
plots, [0.45,0.50], [0.045,0.045], /normal, linestyle=1
xyouts, 0.51, 0.04, alignment=0.0,  '!5sat. specific humidity',charsize=legendsize,$
	charthick=1.5,/normal
plots, [0.45,0.50], [0.025,0.025], /normal, linestyle=0
xyouts, 0.51, 0.02, alignment=0.0,  '!5dry adiabats (!7h!5)',charsize=legendsize,$
	charthick=1.5,/normal

soundingonly:

; plot optional sounding
IF nlev GT 0 THEN BEGIN
  IF n_elements(linestyle) NE 1 THEN linestyle=0. ; solid line is default
  IF n_elements(temp) EQ nlev THEN BEGIN
     ;   Plot the temperature profile
     IF n_elements(col_temp) NE 0 THEN $
          plots, temp, -273.15+theta(pres,temp), /t3d, noclip=0,$
	  color=col_temp, thick=3, linestyle=linestyle $
     ELSE plots,temp,-273.15+theta(pres,temp),/t3d,noclip=0,$
	 	thick=3, linestyle=linestyle
  ENDIF
  IF n_elements(td) EQ nlev THEN BEGIN
     ;   Plot the dewpoint temperature profile
     IF n_elements(col_td) NE 0 THEN $
        plots, td, -273.15+theta(pres,td),/t3d,noclip=0, $
		color=col_td, thick=3 , linestyle=linestyle $
     ELSE plots,td,-273.15+theta(pres,td),/t3d,noclip=0,$
		thick=3, linestyle=linestyle
  ENDIF
ENDIF

; plot optional wind barbs
IF windok THEN BEGIN
   ; calculate x/y normal coordinates of the wind barbs
   x = (pres/1000)^0.286
   delT = (273+T2)*(1-x)/(1+x)
   T = T2 - delT -0.98  &    th = T2 + delT -0.98
   ; direct conversion to normal coordinates did not work
   ; therefore 2 steps: data->device->normal
   p=convert_coord(t,th,/to_device,/t3d)
   p=convert_coord(p,/device,/to_normal)
   index=WHERE(p[1,*] LE yhi)	; only wind barbs within plot range
   FOR i=0,n_elements(index)-1 DO wind_barb,wspeed[index[i]],size=1.7,$
	wangle[index[i]],p[0,index[i]],p[1,index[i]],/nocircle
ENDIF

END
