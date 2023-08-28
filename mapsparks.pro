Pro mapsparks
device, retain=2, decomposed=0
loadct,0

maptype=TextBox(Title='Map Sparks', Label='<1>=all, <2>=repeat evts', Value='1', Cancel=cancelled)
maptype=fix(maptype)

  if maptype eq 2 then begin
    Print, "Select repeat event spark data file'
    file = DIALOG_PICKFILE(/read,/multiple, path='c:\data', filter = '*evtrpt.txt*')
    dataformat='(a,a,a,a,a)'
    readcol,file,rname,cx,cy,cx1,cy1,rptevt,format=dataformat
  endif
  
  if maptype eq 1 then begin
    print, 'Select all sparks data file'
    file = DIALOG_PICKFILE(/read,/multiple, path='c:\data', filter = '*mapallspks.txt*')
    dataformat='(a,a,a,a,a,a,a)'
    readcol,file,rname,cx,cy,bx,by,xxa,tta,format=dataformat
  endif

manselect=TextBox(Title='Manually select t-tubule ROI?', Label=' YES <1>, NO <2>', Value='1', Cancel=cancelled)
manselect=fix(manselect)

pixsz=TextBox(Title='Enter the Pixel Size', Label='Enter pixel size in microns/pixel',Value='0.06',Cancel=cancelled) 
  pixsz=float(pixsz)

print, 'Select t-tubule image file'
  tiffile=DIALOG_PICKFILE(/read,/multiple,path='c:\data')
  tiffarr=read_tiff(tiffile)
  tiffarr=float(tiffarr)
  imp=query_tiff(tiffile,t)
    xs=t.dimensions(0) & ys=t.dimensions(1)
    ix=(float(xs))
    iy=(float(ys))
    nframes=t.num_images
    
    arr=UINTARR(ix,iy,nframes)
    
   ;For I=0,nframes-1 do begin  ;this reads in a stack of images to create a z projection of the spark area.
    ;  arr(*,*,I)=READ_TIFF(tiffile,IMAGE_INDEX=I)
   ; endfor
    
    window,0,xs=ix,ys=iy
    
;    stop
   tvscl,tiffarr;(0,*,*)
STOP    
    color_lut,1
    k=n_elements(rname)
    for i= 1, k-1 do begin 
     ; stop    
        wset,0 & draw_box,cx(i)-2.5,cy(i)-2.5,5,5,200
          
    endfor
  ;stop
  
  
  
   iend=rstrpos(tiffile,'.')
   iname = strmid(tiffile, 0, iend)
   prevfile=file_search(iname+'*spkmap.tif') & evnum=n_elements(prevfile) ;& if prevfile eq '' then evnum=0
   ;stop 
    wset, 0 & im=tvrd(/true ) & write_tiff, iname+'_v'+string(evnum)+'_spkmap.tif',byte(im) 
    ;stop
    loadct,0
    
    if manselect eq 1 then goto, manslct
    
   ;get spark centers and output the  t-tubule image 10um x 10um centered at center of spark
    for i= 1, k-1 do begin
     
      ;stop
      ttbox=round(5/pixsz)
      ttbox=float(ttbox)    
      ttcxa=round(cx(i)-(ttbox/2)) & ttcxb=round(cx(i)+(ttbox/2)) & ttcyb=round(cy(i)+(ttbox/2)) & ttcya=round(cy(i)-(ttbox/2)) ;create image 10um square with center of image at the center of the spakr
      if ttcxa lt 0 then ttcxa=0 & if ttcya lt 0 then ttcya=0
      if ttcxb gt xs then ttcxb= xs-1 & if ttcyb gt ys then ttcyb=ys-1
    ; ttarr=tiffarr(0,ttcxa:ttcxb,ttcya:ttcyb) & txs=ttcxb-ttcxa & tys=ttcyb-ttcya i dont remember why the tiffarr was (0,*,*) so changed on 7/14/2023
     ttarr=tiffarr(ttcxa:ttcxb,ttcya:ttcyb) & txs=ttcxb-ttcxa & tys=ttcyb-ttcya
     color_lut,1
        wset,0 & spktt=tvrd(ttcxa,ttcya,ttbox,ttbox)
        window,3,xs=ttbox,ys=ttbox
        tvscl,spktt
        write_tiff, iname+'_v'+string(evnum)+'_spkmaptt.tif',byte(spktt)
      loadct,0   
      window,1,xs=ttcxb-ttcxa,ys=ttcyb-ttcya
      tvscl,ttarr
      write_tiff, iname+'_v'+string(i)+'_ttmap.tif',byte(ttarr);&tm=tvrd(/true) & write_tiff, iname+'_v'+string(i)+'_ttmap.tif',byte(tm)
      ;stop
      endfor
    
    manslct:
    ttbox=round(5/pixsz)
   
   
   
   
      device, cursor_standard=139
      ctr=0
      loop4:
      
      wset,0
      cursor, cx, cy, 3, /device
      
      case !err of  ; select different options (below) for mouse left (=1),
        ; center (=2) or right (=4)

        1: begin    ;  left mouse button was clicked (= select a spark)
          print, 'Select t-tubule ROI'
         
      ;stop   
      ttcxa=round(cx-(ttbox/2)) & ttcxb=round(cx+(ttbox/2)) & ttcya=round(cy-(ttbox/2)) & ttcyb=round(cy+(ttbox/2))
      if ttcxa lt 0 then ttcxa=0 & if ttcya lt 0 then ttcya=0
      if ttcxb gt xs then ttcxb= xs-1 & if ttcyb gt ys then ttcyb=ys-1
      draw_box,ttcxa,ttcya,ttbox,ttbox,200
       wset,0 & spktt=tvrd(ttcxa,ttcya,ttbox,ttbox,/true)
      ;ttarr=tiffarr(0,ttcxa:ttcxb,ttcya:ttcyb) & txs=ttcxb-ttcxa & tys=ttcyb-ttcya
      ttarr=tiffarr(ttcxa:ttcxb,ttcya:ttcyb) & txs=ttcxb-ttcxa & tys=ttcyb-ttcya
      
      window,2,xs=ttbox,ys=ttbox
      tvscl,spktt
      write_tiff, iname+'_v'+string(evnum)+'_'+string(ctr)+'_spkmaptt.tif',byte(spktt)
      window,1,xs=ttcxb-ttcxa,ys=ttcyb-ttcya
      loadct,0
      tvscl,ttarr
      ttim=TVRD(/true)
      ;stop
      write_tiff, iname+'_v'+string(ctr)+'_ttmap.tif',byte(ttim)
      wait,2
      ctr=ctr+1
    ;if manselect eq 1 then goto, loop4
    
    goto, loop4
  endcase

  2: begin      ;  middle mouse button was clicked

    print, ' No ROI selected !!'
    ;err_msg2='LongEvent'
    goto, done
  endcase
  ;endif

  ;if !mouse.button eq 4 then begin


  4: begin      ;  right mouse button was clicked
    device, cursor_standard=34
print, ' No ROI selected !!'

    ;STOP


    goto,   done
  endcase
endcase
    
    done:
if manselect eq 1 then begin
  wset,0
  maparr=TVRD(/true) & write_tiff, iname+'_v'+string(evnum)+'_maparr.tif',byte(maparr)
  
endif
close, /all, /file

wdelete,0
wdelete,1
wdelete,2
    end
    