pro dat_plot

v=dblarr(4,25)

openr,1,"MF_FUNC_B.dat"
readf,1,v
close,1

set_plot,'ps'
device,filename= 'MF_FUNC_B_1.ps',/color
loadct,13
plot,v(0,*),v(1,*),color=43,xrange=[-3,3],yrange=[0,1]
;oplot,v(0,*),v(2,*),color=22
;oplot,v(0,*),v(3,*),color=40
device,/close

set_plot,'ps'
device,filename= 'MF_FUNC_B_2.ps',/color
loadct,13
plot,v(0,*),v(2,*),color=43,xrange=[-3,3],yrange=[0,50]
device,/close

set_plot,'ps'
device,filename= 'MF_FUNC_B_3.ps',/color
loadct,13
plot,v(0,*),v(3,*),color=43,xrange=[-3,3],yrange=[-5000,5000]
device,/close

end