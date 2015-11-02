pro analyze

p = read_table("photometry.dat", /double)
phase = p[2,*]
t = sort(phase)
phase = phase[t]
mag = p[3,*]
mag = mag[t]
flt = p[0,*]
flt = flt[t]
fltuniq = flt[uniq(flt, sort(flt))]


m = read_table("model.dat", /double)
m_phase = m[0,*]
m_mod = m[4:*,*]

help, m_mod


!p.charsize=1.5
!p.thick=5
!x.thick=2
!y.thick=2
!p.charthick=1
!p.font=0


old_dev = !D.NAME
set_plot, "ps"
DEVICE, FILE='lc.eps',  xsize=25, ysize=35, /color, bits=8, encapsulated=1
DEVICE, /helvetica, /TT_FONT, font_index=3
DEVICE, /times, /italic, font_index = 4, /tt_font
DEVICE, /symbol, font_index = 5
DEVICE, /times,  font_index = 6, /tt_font
DEVICE, /helvetica, /bold, /TT_FONT, font_index=7

colorlist = [25, 50, 160, 210, 255, 120, 10, 40, 90, 130, 10, 40, 185, 25, 50, 160, 210,  230,  170,  0, 160, 0,0, 50,  25, 130,  90,  185, 25, 50, 160, 210,  160, 210, 255, 120]

fltlist =   ['!8U','B','V','R!D!6C!N','!8I!6!DC!N','!8J!6!DSAAO!N','!8H!6!DSAAO!N','!8K!6!DSAAO!N','!8R!6!DJ!N','!8I!6!DJ!N','!8J!6!DCTIO!N','!8H!6!DCTIO!6', $
  '!8K!6!DCTIO!N','C','M','T1','T2','Hp','BT','VT', '','',    '!8V!6!DW!N','!8B!6!DW!N','!8U!6!DW!N','!8W!6!DW!N','!8L!6!DW!N','F621M','','','','','C1','C2','C3','C4', '']
fltoffsets =[0.4,-0.1,+0.3,+0.3,-0.1,    0.3,   0.3,  -0.3, -0.05, -0.1, 0.2, 0.2, -0.1, 0, 0 ,0, 0, 0, 0,0,0,0,0.2,0.2,0.2,0.2,-0.3,0,0,0,0,0,0,0,0,0]
help, fltoffsets
plotsym, 0, /fill, 0.7


xi = 0.08
xf = 0.02
xbord = 0.07
yi = 0.1
yf = 0.01
ybord = 0.03
xwid = (1.0-xi-xf-xbord)/2.0
ywid = (1.0-yi-yf-2*ybord)/3.0



t = where( (flt ge 1) and (flt le 8))

loadct, 0, /silent
plot, [0],[0], /nodata, xrange=[0,1], yrange=[max(mag[t])+0.2,min(mag[t])-0.5], pos=[xi,1.0-yf-ywid, xi+xwid, 1.0-yf], ystyle=1, ytitle='mag'

loadct, 13, /silent

xlab = xi+0.02
ylab = 1.0-yf-0.02
wid = 0.0
for i=1,8 do begin
  t = where(flt eq i, count)
  if (count ge 1) then begin
    oplot, phase[t], mag[t], psym=8, color=colorlist[i-1]
    oplot, m_phase, m_mod[i,*], color=colorlist[i-1]
  endif
  xlab += 1.2*wid
  xyouts, xlab, ylab, fltlist[i-1], width=wid, color=colorlist[i-1], /normal, charsize=1.0
endfor

plotsym, 3, /fill
t = where(flt eq 28, count)
if (count ge 1) then begin
  oplot, phase[t], mag[t], psym=8, color=colorlist[28]
  oplot, m_phase, m_mod[28,*], color=colorlist[28]
endif
  xlab += 1.2*wid
  xyouts, xlab, ylab, fltlist[27], width=wid, color=colorlist[28], /normal, charsize=1.0

t = where( ((flt ge 9) and (flt le 13)) or ( (flt ge 33) and (flt le 36)))
plotsym, 0, /fill, 0.7
loadct, 0, /silent
plot, [0],[0], /nodata, xrange=[0,1], yrange=[max(mag[t])+0.2,min(mag[t])-0.2], pos=[xi+xwid+xbord,1.0-yf-ywid, 1.0-xf, 1.0-yf], ystyle=1, /noerase

loadct, 13, /silent

xlab = xi+xwid+xbord+0.02
ylab = 1.0-yf-0.02
wid = 0.0
for i=9,13 do begin
  t = where(flt eq i, count)
  if (count ge 1) then begin
    oplot, phase[t], mag[t], psym=8, color=colorlist[i-1]
    oplot, m_phase, m_mod[i,*], color=colorlist[i-1]
  endif
   xlab += 1.2*wid
  xyouts, xlab, ylab, fltlist[i-1], width=wid, color=colorlist[i-1], /normal, charsize=1.0
endfor

for i=33,36 do begin
  t = where(flt eq i, count)
  if (count ge 1) then begin
    oplot, phase[t], mag[t], psym=8, color=colorlist[i-1]
    oplot, m_phase, m_mod[i,*], color=colorlist[i-1]
  endif
  xlab += 1.2*wid
  xyouts, xlab, ylab, fltlist[i-1], width=wid, color=colorlist[i-1], /normal, charsize=1.0
endfor



t = where( ((flt ge 14) and (flt le 22)))
plotsym, 0, /fill, 0.7
loadct, 0, /silent
plot, [0],[0], /nodata, xrange=[0,1], yrange=[max(mag[t])+0.2,min(mag[t])-0.2], pos=[xi,yi+ywid+ybord, xi+xwid, yi+ywid+ybord+ywid], ystyle=1, /noerase, ytitle='mag'

loadct, 13, /silent


xlab = xi+0.02
ylab = 1.0-yf-ywid-ybord-0.02
wid = 0.0
for i=14,22 do begin
  t = where(flt eq i, count)
  if (count ge 1) then begin
    oplot, phase[t], mag[t], psym=8, color=colorlist[i-1]
    oplot, m_phase, m_mod[i,*], color=colorlist[i-1]
  endif
  xlab += 1.2*wid
  xyouts, xlab, ylab, fltlist[i-1], width=wid, color=colorlist[i-1], /normal, charsize=1.0
endfor


t = where( ((flt ge 23) and (flt le 27)) )
plotsym, 0, /fill, 0.7
loadct, 0, /silent
plot, [0],[0], /nodata, xrange=[0,1], yrange=[max(mag[t])+0.2,min(mag[t])-0.2], pos=[xi+xwid+xbord,yi+ywid+ybord, 1.0-xf, yi+ywid+ybord+ywid], ystyle=1, /noerase

loadct, 13, /silent
xlab = xi+xwid+xbord+0.02
ylab = 1.0-yf-ywid-ybord-0.02
wid = 0.0
for i=23,27 do begin
  t = where(flt eq i, count)
  if (count ge 1) then begin
    oplot, phase[t], mag[t], psym=8, color=colorlist[i-1]
    oplot, m_phase, m_mod[i,*], color=colorlist[i-1]
  endif
  xlab += 1.2*wid
  xyouts, xlab, ylab, fltlist[i-1], width=wid, color=colorlist[i-1], /normal, charsize=1.0
endfor

v = read_table('velocity.dat', /double)
if (n_elements(v) ge 1) then begin
  phase = v[1,*]
  vel = v[2,*]



  plotsym, 0, /fill, 0.7
  loadct, 0, /silent
  plot, phase,vel, psym=8, xrange=[0,1], pos=[xi, yi, 1.0-xf, yi+ywid], /noerase, ytitle='Velocity [km/s]', yminor=5, xtitle='Phase',/ynozero
  oplot, m_phase, m_mod[0,*]
endif


device, /close
set_plot, old_dev



end
