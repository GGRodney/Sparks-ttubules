
PRO Spark_XY_111722
;automated spark detection program


;auto_spark_xy_d          8/22/00

;ILD program written by Chris Ward (with some suggestions from MFS) based on Peace Cheng's
;line scan (xt) event selection program

;This program is to identify and characterize local elevations in fluorescence in individual
;XY images in a series of images ( = "pic" file) of a given XY field.
;This is a modification of the automatic spark detection program (auto_spark) used for line
;scan images. Locating regions of increased fluorescence and the criteria for event
;acceptance uses Peace Cheng's program, now modified to work for XY images
;We have kept Peace's original code and comments, only making modifications to
;adapt from line scan to XY images
;Program flow involves two "big" loops over all pic files in a specified path
;   Second big loop has four "small" loops over all XY images in a given pic file
;   Before starting first big loop user defines name and path for input pic file
;   and name and path for output spark list file
;   A second (summary) output file has same name, with _summary appended
;Also specifies a line of text for comments/notes that will appear at top of spark list file

;Big loop #1:
;Creates an average image of each pic file, "movies" the file to allow its use or
;rejection and  lets the user define the background area and the fiber area
;Background area should NOT include fluorescent junk
;Fiber area should exclude fluorescent junk, and closely follow the fiber as carefully
;traced with the mouse

;For fiber AOI definition, drag mouse (left button) along

;Big loop #2:
;This loop automatically analyzes each pic file to identify and then characterize and
;save the parameters of all identified events within the selection criterion
;This big loop uses multiple small loop loops over all images in a given pic file to
;carry out the steps in the event identification and selection process
;First restore the average image and the fiber and off fiber masks for this pic file
;   All subsequent analyses are only applied to the predefined fiber area
;Small loop #1
;   Smoothes all images using Peace's smoothing routines
;   Defines potential event locations as anything greater or equal to more than 1.5 sd
;  above mean Note that original Peace routine used each pixel in average line
;(after smoothing)   as mean for that pixel location
;in our XY  adaptation we use average of all XY images (after smoothing) in each pic
;file to give the mean value of each pixel, and calc sd of each pixel relative to same pixel in mean image
;Small loop #2
;   Calculates mean excluding potential event locations
;Small loop #3
;   Calculates a dF/F0 image for each image in pic file using the mean excluding
; possible events  This dF/F0 array is 3x3 smoothed here or below before using it to
;extract any parameter values; however, spark regions are defined on Peace's smoothed images
;Small loop #4
;this is the guts of  Peace's event selection routine (still using his smoothed images)
;Routine systematically moves over an image identifying contiguous regions of pixels having fluorescence
;more than 2 sd above the mean (note in this version the mean is the smoothed mean image
;excluding potential event regions; potential
;event regions are defined as being greater than 1.5 sd above mean for that pixel, and are excluded
;from the calculation of the mean)Each identified region of contiguous pixels having values greater
;than 2 sd above "mean" is then checked to see if it has at least one pixel greater than the criterion
;For this program the criterion is 3 sd greater than the mean of a given pixel If an event fulfills the
;criterion, a green box is drawn around it in each of the 4 windows:
;upper L,  our smoothed dF/F0;  upper R, Peace's smoothed F/F0;
;lower L,  criterion ( = at least 3sd) greater than the mean of a given pixel If an event
;fulfills the criterion, a green box is drawn around it in each of the 4 windows:  upper L,
; our smoothed dF/F0;  upper R, Peace's smoothed F/F0 lower L,  criterion ( = at least 3sd);
;lower R, at least 2 sd

;Following parameters extracted for each event above the criterion during small loop #4
;   Name of pic file;  Overall # of spark in output list;  Image # in pic file;
;Spark # in this image;  Lower L corner x coordinate;  Lower L corner y
;coordinate;  x length of box (in pixels);  y length of box (in pixels);  peak dF/F;  area (um**2)
;@ 50% max dF/F;  mean dF/F @ 50% max;  area (um**2) @ 2 sd ;  mean dF/F @ 2 sd.
;These are put into successive columns of the output array after completion of big loop #2

;Second output file lists following for each pic file analyzed: file name of pic file,
;area of fiber AOI (in um**2), # of images in the pic file, total # of sparks in the pic file
; fixed an error in the area measurements on 3.28.01   pixels*psize changed to pixels * (psize*psize)

;11/16/2004 altered the program to allow user to analyze a subset of images within the stack.  The images should be selected using
;a typical numbering system, ie. 1 through end.  The program handles the conversion to IDL numbering.  Thus, the search for previous mask arrays
;was removed and the text file outname was corrected to search for previous file versions.

;04/04/2005 added the ability to only analyze a subset of the original stack size as well as write to the data file
;the fiber flourescence - mean background for each image within the stack,

;05/12/2005: included read_tiff for reading in Olympus FV 500 file structure GR

;02/01/2006 INCORPORATED Zeiss LSM file ability

;November 17, 2022, to map sparks on t-tuble image need to crop the spark image from 1864x1200 to 1840x1176. added option to do this.
;
; 
;***************************  this block added to read our images
;parameter input

;bell,/force
sspeed=2.0913              ;ms per line
ZOOM=1.63                 ;zoom set
psize=0.2934/zoom            ;
black=10                     ;used for cell edge detection. Cell region>black
human=0                         ;user interface 1=active, other=inactive
;crii=[3.5,3.8,4.2,4.8]        ;criteria for spark region folds of SD above mean
 crii=[2.0,2.5,3.0,3.8,4.2]   ;  this is my array************

iext = 'img'
iroot = ''
fname = ''
oname = ''
pname = ''
text = ''
tname = ''
rname = ''
sname = ''
lext = ''
str = ''
stg=''
flname=''
imgtype=''
svname=''
spknote=''
boxno=0
imgno=-1

tskc=0
acceptfile=0

device, retain=2, decomposed=0

;loadct,0
;color_lut, 1



start0:



slctscope:
;scope=textbox(Title='Select system used.', Label='BioRad Radiance<1> Olympus FV500<2> or Zeiss LSM<3>',Value='1',Cancel=cancelled);... default is <1>'
;   if scope eq '' or cancelled then begin
 ;      Bell,/allwrong
  ;     goto, slctscope
   ;endif
scope=2
scope=fix(scope)
;if scope eq 1 then scopeext='*_raw.pic'
if scope eq 2 then scopeext='*.tif'
;if scope eq 3 then scopeext='*.lsm'

entrpixsz:
if scope eq 2 then begin
    pixsz=TextBox(Title='Enter the Pixel Size', Label='Enter pixel size in microns/pixel',Value='0.21',Cancel=cancelled)
        if pixsz eq '' or cancelled then begin
        bell,/wrongbtn
        goto, entrpixsz
        endif
    pixsz=float(pixsz)
endif

if scope eq 2 or scope eq 3 then begin
entrsspeed:
    sspeed=TextBox(Title='Scan Speed',Label='Enter the line scan .',Value='1.26',Cancel=cancelled)
       if sspeed eq '' or cancelled then begin
         bell,/wrongbtn
         goto, entrsspeed
       endif
    sspeed=float(sspeed)
 endif

entrcri_SD:
cri_SD=TextBox(Title='Enter the SD value to use',Label='<1>=2.0, <2>=2.5, <3>=3.0,<4>=3.8, <5>=4.2',Value='2',Cancel=cancelled)
       if cri_SD eq '' or cancelled then begin
         bell,/allwrong
         goto, entrcri_SD
       endif
         cri_SD=fix(cri_SD)

drawbox=TextBox(Title=' Draw the box on the image?',Label='<1>=yes',Value='1',Cancel=cancelled)
    drawbox=fix(drawbox)

;oldmsk=TextBox(Title=' Evaluate pre-naalyzed fibers?',Label='<1>=no',Value='1')
;   oldmsk=fix(oldmsk)
 ;       if oldmsk ne 1 then goto, skipautoedge

saveimgs=TextBox(Title='Save these screen images???', Label=' <1> YES <2>=NO', Value='2', Cancel=cancelled)
    saveimgs=fix(saveimgs)

autoedge=TextBox(Title=' Use the autoedge routine?',Label=' <1>=YES, <2>=NO',Value='2',Cancel=cancelled)
    autoedge=fix(autoedge)

dorptevts:
;rptevts=TextBox(Title='Search for longevents and hotspots',Label='Enter the spatial criteria in microns',Value='0',Cancel=cancelled)
   ; If rptevts eq '' or cancelled then begin
    ; bell,/wrongbtn
    ; goto, dorptevts
   ; endif
    rptevts=FLOAT(0.5)
;stop
slcClrLut:
ClrLut=TextBox(Title='Display Images as GreyScale',Label='GreyScale <1> or RGB <2>',Value='1',Cancel=cancelled)
    ClrLut=fix(ClrLut)
    if ClrLut eq 1 then loadct,0 else color_lut, 1

cropimg:
croparr=TextBox(Title='Crop the spark image to match t-tubule image?', Label= 'YES <1> OR NO <2>', Value='2',Cancel=canceled)
  croparr=fix(croparr)
  

skipautoedge:
oldmsk = 1
if oldmsk ne 1 then scopeext = '*msk.tif'

print,'Select spark file'
file = DIALOG_PICKFILE(/read,/multiple, filter = scopeext,path='c:\data')

nimgs=n_elements(file)



; initialize some variables
;    Each of the following is an array for indicated spark info
  rnamea = ''   ;individual image name
  fiberfloura = 0.000  ;fiber flourescence - mean background for each image within the stack, added 04/04/2005 by GR
  numa = 0  ;spark # in list
  imgnuma = 0   ; # (starting from 0) of pic file in list of root*.pic
         ;pic input image file
  boxnoa = 0    ;# of the box in this *.pic image file
  cxa = 0   ;  this is the lower left x box coordinate
  cya = 0   ; this is the lower left y box coordinate
  bxa = 0   ; x size of box
  bya = 0   ;y size of box
  xxa = 0   ;x center of possible spark region
  tta = 0   ;y center of possible spark region
  DFF50area = 0    ; area at 50% dff
  DFFamp = 0     ;  DF/F amp in box
  meanDFFamp = 0  ; mean DF/F amp in box
  SD2area = 0   ; sd2 AREA
  SD2DFFamp = 0.    ; peak DF/F amplitude in 2SDarea in box

  naoi = 0    ; this not currently used
  iapp = 0
  nmax = 10000   ; max # of box locations available
  naois = 0




       ;if oldmsk ne 1 then begin
       ;  for k=0, nimgs-1 do begin
       ;     omsk=read_tiff(file(k))
        ;    sz = size(omsk) & ix = sz(1)  & iy = sz(2)
        ;    if nimgs eq 1 then iz=2 else iz=nimgs
       ;     if k eq 0 then xarr=fltarr(ix, iy, iz)
       ;     xarr(*, *, k)=omsk
       ;  endfor
       ;  acceptfile=nimgs-1

       ;  goto, Previous_Mask
       ;endif

for k=0, nimgs-1 do begin
 ;*** start of the first loop (Big loop #1)which loops over the number
                ;of pic files found in the path to view the
              ; files and determine the aoi's of each fiber.

ff=STRARR(n_elements(file))


          if k gt 0 then arr=0 &  xarraoi=0
          imgno=imgno+1



          if sname ne file(imgno) then begin  ; start of block for input file info

         ; if scope eq 1 then  begin
          ;  Rd_BioRad_xy_Pic, file(k), arr, nn, Header,Notes,NotesArray, RGBLUT, pixsz, pixdepth, sspeed, pixdwell, framerate,  ix, iy, nframes
          ;        arr=float(arr)

            ;       iok = rstrpos(file(k),'\')
            ;          iroot = strmid(file(k),0,iok+1)
            ;          nfiles = n_elements(file)
            ;          iend = rstrpos(file(k),'.')
            ;          iname = strmid(file(k), 0, iend)
            ;          iname = strmid(iname, iok+1)
               ;  read in images 1 to n into arr using
          ; endif

         if scope eq 2  then begin
         ;loadct,0
               tiffarr=read_tiff(file(k))
               tiffarr=float(tiffarr)
              imp = QUERY_TIFF(file(k),t)
              xs=t.dimensions(0) & ys=t.dimensions(1)
              ix=(float(xs))
              iy=(float(ys))
              nframes=t.num_images
              pixdepth=t.bits_per_sample
              framerate=float(00)
              pixdwell=float(00)
              iok = rstrpos(file(k),'\')
              iroot = strmid(file(k),0,iok+1)
              iend = rstrpos(file(k),'.')
               iname = strmid(file(k), 0, iend)
               iname = strmid(iname, iok+1)
               
               
           arr=UINTARR(ix,iy,nframes)

         For I=0,nframes-1 do begin
          arr(*,*,I)=READ_TIFF(file(k),IMAGE_INDEX=I)
         endfor
         
         if croparr eq 1 then begin
           print,'Select t-tubule image file'
           ttimage=DIALOG_PICKFILE(/read,/multiple, filter = scopeext,path='c:\data')
           ttarr=read_tiff(ttimage)
           ttarr=float(ttarr)
           ttimp=query_tiff(ttimage,tt)
           txs=tt.dimensions(0) & tys=tt.dimensions(1)
           tix=(float(txs))
           tiy=(float(tys))
           tnframes=tt.num_images
          ; stop
           dxs=(ix-tix)/2 & nxa=dxs & nxb=ix-dxs
           dys=(iy-tiy)/2 & nya=dys & nyb=iy-dys
           tarr=arr(nxa+1:nxb,nya:nyb-1,*)
           arr=tarr
          nix=ix-tix & niy=iy-tiy
          ix=ix-nix & iy=iy-niy
         endif
           endif
;stop
      ; if scope eq 3 then begin
      ;        LSM5_reader, file(k),a, ix, iy, nframes, pixsz_x, pixsz_y,scantime,channel,bitdepth
      ;        pixsz=pixsz_x*1000000
      ;           psize=pixsz
      ;           piysz=pixsz_y*1000000
      ;           ix=ix
      ;          arr=float(a)
      ;          pixdepth=bitdepth

      ;         nimgs=nframes

       ;       iok = rstrpos(file(k),'\')
       ;       iroot = strmid(file(k),0,iok+1)
       ;       iend = rstrpos(file(k),'.')
       ;        iname = strmid(file(k), 0, iend)
       ;        iname = strmid(iname, iok+1)

       ;  endif



                 if k eq 0 then aoiarr=fltarr(ix, iy, nimgs)
                    if k eq 0 then oldmskarr=fltarr(ix, iy, nimgs)
                    if k eq 0 then window,0, xs=ix, ys=iy

    ;OLD MASK SEARCH IS  to be COMMENTED OUT if you want to analyze SUBSET OF STACK, will need to fix this and allow for
    ; subset previous masks to be correctly found

         ;oldmsk = '' &  oldmsk=file_search(iroot+iname+'*_msk.tif')

              ;if oldmsk ne '' then begin
                  ;oldmskarr(*,*,k)=read_tiff(oldmsk)
                 ;acceptfile=acceptfile+1
                 ;slctimgs=2
                 ;firstimg=0
                 ;lastimg=fix(string(nframes))
                  ;goto, eol_rej
              ;endif

                 ; array aoi handler
                 ; -1 is the bkgd aoi
                 ; +1 is the fiber aoi
                 ;  0 is the rest


          window,0,xs=ix,ys=iy,title=file(k)

                   avgimg=rebin(arr, ix,iy, 1)
                   
                   if scope eq 1 then tv, avgimg else tvscl,avgimg

                   ;*******  this gains the image up to a user defined setting
                   gain=0
                   xgain=1
                   regain:
                   if scope eq 1 then begin
                        gain=TextBox(Title='Image Gain for viewer purposes',Label='Enter a number from <1-255> to "gain" the average image; type <0> to continue',Value='3')
                        gain=fix(gain)
                   endif
              ;window,0,xs=ix,ys=iy,title=file(k)
                   if scope eq 1 and gain ne 0 then begin


                     tv, avgimg*gain
                     xgain=gain

                     goto, regain
                   endif


                   ; *****  this movies the image stack so it can be accepted or rejected
                   removie:
                   xmovie=TextBox(Title='Movie through the stack of Images?', Label='enter <1> to preview the movie loop; <0> to exit',Value='0')
                   xmovie=fix(xmovie)
                   if xmovie eq 1 then begin
                 ;window,0,xs=ix,ys=iy,title=file(k)
                    wset,0

                    if scope eq 1 then movie, arr*xgain
                    if scope eq 2 or scope eq 3 then movie_mtiff,arr
                     ;goto, removie
                   endif

                 ;*** this will allow a user defined rejection of the pic file after viewing the movie
                 accpt=TextBox(Title='Use this file for Analysis?',Label='press  <1> to accept this file; <0> to reject',Value='1')
                   accpt=fix(accpt)

                   if accpt eq 0 then begin
                    ff(imgno)='rejected file'
                    goto, eol_rej
                   endif
                   acceptfile=acceptfile+1
;***This will allow for a reduction in the image stack to a specified number of images*********
          slctimgs=TextBox(Title='Do you want to analyze a subset of the image stack?',Label='<1>=YES   <2>=NO',Value='2')
          slctimgs=fix(slctimgs)
              if slctimgs eq 1 then begin
                 firstimg=TextBox(Title='What image do you want to start with?',Label='Enter first image, default is 0',Value='0')
                    if firstimg ne 0 then firstimg=firstimg-1 else firstimg=firstimg
                    firstimg=fix(firstimg)
                  lastimg=TextBox(Title='What image do you want to end with?',Label='Enter last image, defaults is'+string(nframes),Value=string(nframes))
                    ;lastimg=lastimg-1
                    lastimg=fix(lastimg)
              endif
              if slctimgs ne 1 then begin
                 firstimg=0
                 lastimg=fix(string(nframes))
              endif
;************************************************************************************************
         window,0,xs=ix,ys=iy,title=file(k)
         wset,0

              if scope eq 1 then tv, avgimg*xgain
           if scope eq 2 or scope eq 3 then tvscl,avgimg

                   xarraoi=fltarr(ix, iy) & xarraoi(*,*)=0.0

                   redobkgaoi:


                   print, 'select AOI for background'
                   bkgaoi=defroi(ix, iy)

                        device, cursor_standard=139
                        print,'*** Mouse buttons: <L> Accept this  <R> Re-do'
                        ; start of loop for mouse accepting of aoi's
                        cursor, ccx, ccy, 3, /device  ; wait for any mouse click and store
                           pnt=0
                           pnt=!mouse.button
                        ; cursor location in ccx,ccy
                          case !err of  ; select different options (below) for mouse left (=1),
                                 ; center (=2) or right (=4)
                   1: begin ;  left mouse button was clicked (= select a spark)
                      endcase

                   4: begin
                   bkgaoi=0
                   BELL,/bugs2
                   goto, redobkgaoi
                   endcase
                   endcase



                   redofibaoi:

                    wait,1
                    print, 'Select AOI for fiber'
                    fiberaoi=defroi(ix, iy)

                        device, cursor_standard=139
                        print,'*** Mouse buttons: <L> Accept this  <R> Re-do'
                        ; start of loop for mouse accepting of aoi's
                        cursor, ccx, ccy, 3, /device  ; wait for any mouse click and store
                           pnt=0
                           pnt=!mouse.button
                        ; cursor location in ccx,ccy
                            case !err of ; select different options (below) for mouse left (=1),
                                 ; center (=2) or right (=4)

                   1: begin
                    xarraoi(bkgaoi)=-1.  & xarraoi(fiberaoi)=1.
                    aoiarr(*,*,k)=xarraoi
                   endcase

                   4: begin
                    fiberaoi=0
                    BELL,/bugs2
                    goto, redofibaoi
                   endcase
                   endcase





          endif


          ;***************************



          eol_rej:

endfor ;*** end of the loop (Big loop #1)  over the number of pic files found in the path to view the
                   ;     files and determine the aoi's of each fiber.


         ;nimgs=imgno+1   ; this adjusts the pic file counter to the # of accepted images above
         ;imgno=-1      ; reset the imgno to -1


         arr=0 & xarraoi=1

         window,3,xs=ix,ys=iy,Title='cri*SD binary image'
         window,1,xs=ix,ys=iy,Title='2SD binary image'
         window,2,xs=ix,ys=iy,Title='original normalized image'
         window,0,xs=ix,ys=iy,Title='original normalized image overlayed with 2SD sparks'
         ;loadct,0
;stop
;*** Start of Big loop #2 which loops over the number of pic files
for k=0, acceptfile-1 do begin

                          ;set up some arrays for data storage, reset upon each new file.
                          rnamea = [rnamea, strarr(nmax)]      ; make room for nmax more entries
                    fiberfloura=[fiberfloura,fltarr(nmax)]
                           numa = [numa, intarr(nmax)]
                           imgnuma = [imgnuma, intarr(nmax)]
                           boxnoa = [boxnoa, intarr(nmax)]
                           cxa = [cxa, intarr(nmax)]
                           cya = [cya, intarr(nmax)]
                           bxa = [bxa, intarr(nmax)]
                           bya = [bya, intarr(nmax)]
                           xxa = [xxa, intarr(nmax)]
                           tta = [tta, intarr(nmax)]
                           DFF50area= [DFF50area, fltarr(nmax)]
                           DFFamp = [DFFamp, fltarr(nmax)]
                           meanDFFamp = [meanDFFamp, fltarr(nmax)]
                           SD2area = [SD2area, fltarr(nmax)]
                           SD2DFFamp = [SD2DFFamp, fltarr(nmax)]
                           naoi = [naoi, fltarr(nmax)]


                 print, 'File #  '+string(k)+' of   '+string(acceptfile)



              imgno=imgno+1


                   sname='dsdsdsds'  ; this in for testing............
                 if sname ne file(k) then begin  ; **************start of block for input file info

                               if ff(k) eq 'rejected file' then goto, eol
                               if scope eq 1 then begin
                                    Rd_BioRad_xy_Pic, file(k),arr, nn, Header,Notes,NotesArray, RGBLUT, pixsz, pixdepth, sspeed, pixdwell, framerate,  ix, iy, nframes
                                    arr=float(arr)
                               endif

                                if scope eq 2 then begin
                         ;loadct,0
                             tiffarr=read_tiff(file(k))
                             tiffarr=float(tiffarr)
                            imp = QUERY_TIFF(file(k),t)
                            xs=t.dimensions(0) & ys=t.dimensions(1)
                            ix=(float(xs))
                            iy=(float(ys))
                            nframes=t.num_images
                            iok = rstrpos(file(k),'\')
                            iroot = strmid(file(k),0,iok+1)
                            iend = rstrpos(file(k),'.')
                             iname = strmid(file(k), 0, iend)
                             iname = strmid(iname, iok+1)
                         arr=UINTARR(ix,iy,nframes)

                       For I=0,nframes-1 do begin
                        arr(*,*,I)=READ_TIFF(file(k),IMAGE_INDEX=I)
                       endfor
                       if croparr eq 1 then begin
                         dxs=(ix-tix)/2 & nxa=dxs & nxb=ix-dxs
                         dys=(iy-tiy)/2 & nya=dys & nyb=iy-dys
                         tarr=arr(nxa+1:nxb,nya:nyb-1,*)
                         arr=tarr
                         nix=ix-tix & niy=iy-tiy
                         ix=ix-nix & iy=iy-niy
                       endif

                         endif
;stop
                 if scope eq 3 then begin
              LSM5_reader, file(k),a, ix, iy, nframes, pixsz_x, pixsz_y,scantime,channel,bitdepth
              pixsz=pixsz_x*1000000
                 psize=pixsz
                 piysz=pixsz_y*1000000
                 ix=ix
                arr=float(a)
              pixdepth=bitdepth
               nimgs=nframes

              iok = rstrpos(file(k),'\')
              iroot = strmid(file(k),0,iok+1)
              iend = rstrpos(file(k),'.')
               iname = strmid(file(k), 0, iend)
               iname = strmid(iname, iok+1)

         endif
                   ;stop
                                if slctimgs eq 1 then n=lastimg-firstimg else n=nframes
                                 numframes=n ;this is the number of frames analyzed in the run, could be less than nframes if user
                                          ;selected to run a subset of the total number of frames in the pic file

                                incrarr=firstimg
                                narr=fltarr(ix,iy,n)
                              for kk=0, n-1 do begin
                                 narr(*,*,kk)=arr(*,*,incrarr)
                                 incrarr=incrarr+1
                              endfor

                             arr=narr

                                avgimg=rebin(arr, ix,iy, 1)

                                 iok = rstrpos(file(k),'\')
                                    iroot = strmid(file(k),0,iok+1)
                                    nfiles = n_elements(file)
                                    iend = rstrpos(file(k),'.')
                                    iname = strmid(file(k), 0, iend)
                                    iname = strmid(iname, iok+1)

                        prevmsk=0
                               if max(oldmskarr(*,*,k)) ge 2 then begin
                                    prevmsk=1  ; set a flag for a previous mask being used
                                    bkgaoi= where(oldmskarr(*,*,k) le 1.)
                                    fiberaoi=where(oldmskarr(*,*,k) gt 1.)
                                    maskimg=oldmskarr & maskimg(fiberaoi)=1.0 & maskimg(bkgaoi)=0.0
                                    aoi_matr=FLTARR(n-1) ; the matrix to store the avg. value of the AOI
                                    aoi_matr=TRANSPOSE(fiberaoi);  2d to 1d array for zoltans prog

                                 goto, skipaoi
                               endif


                                 bkgaoi= where(aoiarr(*,*,k) eq -1.)
                                 fiberaoi=where(aoiarr(*,*,k) eq 1.)

                                  xarraoi=fltarr(ix, iy) & xarraoi(*,*)=0.0

                                  xarraoi(bkgaoi)=-1.  & xarraoi(fiberaoi)=1.
                                  aoiarr(*,*,k)=xarraoi

                                 ; these are arrays set up for the edge detection, some are used in this program
                                 aoi_matr=FLTARR(n-1) ; the matrix to store the avg. value of the AOI
                                 aoi_matr=TRANSPOSE(fiberaoi);  2d to 1d array for zoltans prog
                                 img_aoi=bytarr(ix,iy) ; the matrix to store the current img.
                                 maskimg=bytarr(ix,iy) ; initialize the mask
                                 maskimg(aoi_matr)=1 ;keep pixels within fiber ROI only
                                 sub_img=bytarr(ix,iy) ;initialize sub-image
                                 newarr=arr
                                 mask=maskimg(*,*)

                                 if autoedge eq 1 then begin
;stop

                                           ; ***the following is edge detection routine
                                                 sub_img=reform(avgimg(*,*,*))*mask ;create sub-img based on rough mask,
                                                   ;to exlcude bright spots outside the cell before autoedge is applied
                                                 mask_auto=mask ;initialize auto edge mask
                                               threshold=2.0 ;the multiplying factor of MEAN:
                                                 mask_box=3 ;the box size used in calculating median in the following:
                                                 mask_auto(WHERE(MEDIAN(sub_img,mask_box) LT threshold*MEAN(sub_img)))=0.3 ;keep cell-area only,
                                                     ;based on intensity relative to MEAN, very sensitive to the factor
                                                     ;the factor of 0.7 seems to work well but it may have to be adjusted
                                                     ;for other intensity ranges
                                                 sub_img_auto=sub_img ;initialize new sub-img for autoedge result
                                                 sub_img_auto=sub_img*mask_auto ;modify masked img according to autoedge
                                                 spots=WHERE(((sub_img_auto EQ 0) AND (MEDIAN(sub_img_auto,mask_box*3) GT threshold*MEAN(sub_img))), count) ;use box twice the size of
                                                   ;in order to find spots that were blacked out, erraneously, within the cell
                                                 IF count NE 0 THEN mask_auto(spots)=1 ;restore mask only if there are spots
                                                 sub_img_auto=sub_img*mask_auto ;recalc sub img to correct for lost spots
                                                 xavgimg=sub_img_auto ;replace actual img w/ masked img


                                                xarraoi=fltarr(ix, iy) & xarraoi(*,*)=1.0
                                                bkgaoi=where(xavgimg le 0.)
                                                fiberaoi=where(xavgimg gt 0.)

                                 endif




                               ; these are two loops which are needed for the edge detection routine to suntract the Bkround F
                               ; and then mask the background to zero.

                               for i =0, n-1 do begin
                                     img_aoi=arr(*,*,i)
                                     aoi_matr(i)=total(img_aoi(bkgaoi))/n_elements(img_aoi(bkgaoi))
                               endfor

               skipaoi:
                                 xarr=arr
                               FOR i=0,n-1 DO BEGIN

                                        if oldmsk EQ '' THEN xarr(*,*,i)=(temporary(xarr(*,*,i))-aoi_matr(i))*((temporary(xarr(*,*,i))-aoi_matr(i))ge 0)
                                        if oldmsk ne '' THEN aoi_matr=TRANSPOSE(fiberaoi)
                                      xarr(*,*,i)=temporary(xarr(*,*,i))*maskimg

                               ENDFOR


                                      ;subroutine rd_pic.  Note image #
                                      window,0,xs=ix,ys=iy,Title='original normalized image overlayed with 2SD sparks for   '+file(k)          ; is printed from within subroutine
                               wset,0  & erase          ; is printed from within subroutine

                    if scope eq 1 then begin
                        if n eq 1 then tv, arr else tv, arr(*,*,0)  ; display first 8 bit image in
                                             ; input (*.pic) file
                    endif

                    if scope eq 2 then begin
                       if n eq 1 then tvscl, arr else tvscl, arr(*,*,0)
                    endif

                               arr = float(arr)  ;  convert input image to floating point

                               ;sz = size(arr)   ;  get input *.pic file  array size dimensions
                               if n eq 1 then begin
                                   ;  ix = sz(1)
                                    ; iy = sz(2)
                                   endif else begin
                                 ;  ix = sz(2)        ; assume arr(ix, iy, z)
                                   ;  iy = sz(3)
                                endelse
                   endif    ; end of the rd_pic end if above....





                             ;*****  these are arrays for the mask and images for each image in the 2nd looping
                             maskima=fltarr(ix,iy,n)
                             imbarr= fltarr(ix,iy,n)

                   ;****************************************  end of our addition for bio-rad images
                   ;*** begining of the auto-spark routine

                   ; *********** start of first small loop (#1) over the number of images in each pic file
                   pli=0  ; pic loop increment

                 for j = 0, n-1 do begin    ;added to increment the images in the pic file


                        ; this sets the criterion at the begining of the prog.
                 ccrit=fix(cri_SD-1)
                        cri=crii(ccrit)
;stop
                        ;input of image data

                        a=arr(*,*,j)
                        y=iy

                        ny=iy
                        nx=ix

                        ima=fltarr(ix,iy)
                        ima1=fltarr(ix, iy)


                        ima=arr(*,*, j)
                        ima=reform(ima)

                        ima1=arr(*,*, j)
                        ima1=reform(ima1)



                        ;setup LUT


                        tvlct,r,g,b,/get
                        ;****************  commeneted this out to look make the DFF color table stay intact
                        ;r(220)=255 & g(220)=0 & b(220)=0   ;red for spark region
                        ;r(221)=0 & g(221)=255 &b (221)=0
                        tvlct,r,g,b, /get


                        ;iterative computation of variance or SD

                        ;ima=float(median(ima,5))
                         ;ima=float(median(ima,5))$
                         ima(fiberaoi)=float(median(ima(fiberaoi),5))

                                            ;minimal median filter to remove data points at extremes
                        imm=ima                  ;keep a copy to be used later
                        ;ss=(0.8/pixsz)&st=(10./sspeed)&ct=0&imb=fltarr(nx,ny)
                                            ;0.8-um and 10-ms spatiotemproal smoothing filter
                        ss=(0.8/pixsz)&st=(0.8/pixsz)&ct=0&imb=fltarr(nx,ny) ;changed 6/16/05 to reflect a square pixel dimension for XY images

                        for ia=-fix(ss/2),fix(ss/2) do begin &imb=imb+shift(ima,ia,0)&ct=ct+1&endfor
                        for ib=-fix(st/2),fix(st/2) do begin &imb=imb+ shift(ima,0,ib)&ct=ct+1&endfor
                        ima=imb/ct
                        imb=ima
                        ;a=rebin(ima,nx,1)&a=where(a gt black)&pl=min(a)&pr=max(a)
                                            ;for linescans---cell edges detected and then stored as pl and pr
                        ;ima=ima/rebin(rebin(ima,nx,1),nx,ny)

                         ima=ima/avgimg

                        ;initial normalization, ima=ima/mean(ima) or ima=(ima/imagestack)



                            sd=stdev(ima(fiberaoi))       ;inital estimate SD of the image data

                   multsd=1.5

                            mask=bytarr(nx,ny)&mask(where(ima gt 1+multsd*sd))=1 & mask=median(mask,5)
                                                             ; not sure what this median does


                                         ;mask for potential spark regions >m+1.5SD
                        ima=imb*(1-mask)          ;excise potential spark regions
                        maskima(*,*, j)=mask          ; save mask for possible event locations in jth image
                        imbarr(*,*,j)= imb          ; save the jth smoothed image, both for use below

                        ;print, 'SD', sd
                   endfor       ; *********** end of first small loop (#1) over the numner of images in each pic file
                           ;** possible spark regions now identified in this set of images


                    ; **** Start of loops (#2-3) to exclude possible spark areas from average image
                    ;      and to calculate DF/F image stack


                    zarr=fltarr(ix, iy) ; arrays defined and set to 0
                    zmask=fltarr(ix, iy)
                    avgmsk=fltarr(ix, iy)& avgmsk(*,*)=1.

                    for j = 0, n-1 do begin
                          ;tv, (j, *,*)
                         zarr=zarr+arr(*,*,j)*(1-maskima(*,*,j))
                         zmask=zmask+(1-maskima(*,*,j))
                    endfor
                        ; this is the average of all images excluding potential spark areas
                        ; it is used below to normalize each imb image
                       avgmsk(fiberaoi)=zarr(fiberaoi)/zmask(fiberaoi)

                       dfarr=fltarr(ix,iy,n)
                       dfarr(*,*,*)=1.
                     for j = 0, n-1 do begin
                          ;make a DF/F array for use later
                               xxarr=xarr(*,*,j)
                               xxarr=(xxarr+aoi_matr(j)-(avgmsk))/(avgmsk)& xxarr=reform(xxarr)
                               fibxxarr=xxarr & fibxxarr(*,*) =0. & fibxxarr(fiberaoi)=xxarr(fiberaoi)
                               dfarr(*,*,j)=fibxxarr

                     endfor        ;**** end of loops (#2-3) to exclude possible spark areas from average image
                           ;     and to calculate DF/F image stack

                   ;******* Start of loop #4 to take care of spark areas over the number
                   ;       of images in each pic file

                   imb=fltarr(ix, iy)
                   mask=fltarr(ix, iy)

          curimg=firstimg-1

          FiberFlrA=fltarr(n);This array is to hold the background corrected fiber fluorescence with sparks
                FiberFimg=fltarr(n)       ;extracted.  This array is to be output to a separate text file  GR 05/04/2005
         imgBKG=fltarr(n)
         imz=fltarr(ix,iy,n)
                   for j = 0, n-1 do begin    ;added to increment the images in the pic file
                        skc=0                  ;initalize spark counter
                   curimg=curimg+1
                            ;get a file name to be able to write the files
                            ;if j lt 10 then outname=strtrim(iname+'_0'+strtrim(string(j), 1), 2)
                            ;if j ge 10 then outname=strtrim(iname+'_'+strtrim(string(j), 1), 2)
                            if curimg lt 10 then outname=strtrim(iname+'_0'+strtrim(string(curimg), 1), 2)
                            if curimg ge 10 then outname=strtrim(iname+'_'+strtrim(string(curimg), 1), 2)

                        ;stop


                        imb=imbarr(*,*, j)  ; select jth image in imb array
                        FiberFlra(j)=mean(imb(fiberaoi)) ;calculate the mean fiber fluorescence for image J

                        imgBKG(j)=mean(imb(bkgaoi))
                        FiberFimg(j)=FiberFlra(j)-imgBKG(j)
                        imb=reform(imb)
                        mask=maskima(*,*, j) ; select jth image in mask array
                        FiberFimg(j)=FiberFimg(j)-mean(mask(fiberaoi));calculated the meand fiber F minuse potential sparks for
                                                      ;image J
                        mask=reform(mask)

print,fiberfimg(0:j)


                            base = avgmsk
                            imb(fiberaoi)=imb(fiberaoi)/base(fiberaoi)  ; normalize imb by avgimg excluding possible events
                                                      ;  this is smooth version of F/Fo

                            tem=fltarr(ix, iy)& tem(*,*)=0.
                            tem(fiberaoi)=imb(fiberaoi)*(1.-mask(fiberaoi))


                            sd=stdev(tem(where(tem gt 0.)),mean)

                            corimb=imb & corimb(*,*)=0. & corimb(fiberaoi)=imb(fiberaoi) & imb=corimb
                            ; this corrects imb , setting pixels outside fiberarea to 1.

                        ;Spark detection

                        im=median(fix((imb-(1.+cri*sd)>0)*100000.<1),5)
                                               ;binary image of  cri*SD, spark sites
                                               ;median filter to reject subsize islets and
                                               ;single pixels in the raw cri*SD binary image
                        im(0:ss/2,*)=0&im(nx-ss/2-1:*,*)=0&im(*,0:st/2)=0&im(*,ny-st/2-1:*)=0
                                            ;take care the edge effect of the smoothing
                        ime=median(fix((imb-(1.+2.0*sd)>0)*100000.<1),5)
                                            ;2SD image for automated regional counting
                        ime(0:ss/2,*)=0&ime(nx-ss/2-1:*,*)=0&ime(*,0:st/2)=0&ime(*,ny-st/2-1:*)=0
                        imf=ime

                        ;**************  is the value for the numbers?????????





                        ;imb=imm/rebin(base,nx,ny)   ;this is from peace's origional prog


                        wset,0&tv,(imb-.5>0)*50<250    ;display the normalized image with contrast enhanced
                                            ; *****This image is F/F0 where F0 is avg line minus
                                            ;      potential sparks excluded
                                            ; for xy fo is average image minus potential sparks excluded
                        ;wset,2&tv,(imb-.5>0)*150<250  ; I commented this out***

                        rr=dfarr(*,*,j) &  rr=reform(rr) & rr=smooth(rr, 3) & rr=rr*maskimg

                        wset,2&tvscl, rr;*60+5>1  ;  tv the scaled DF/F image for viewing.........
                        ;stop
                        wset,3&tvscl,im           ;cri*SD image
                        wset,1&tvscl,ime                ;2-SD image


                        imz(*,*,j)=im
;stop

                        ;***** define the 2SD and the crit parameters
                        img=fltarr(ix, iy) & sd2=where(ime ge 0.1) & sd25=where(im ge 0.1) & img(sd2(*,*))=155.& img(sd25(*,*))=200.


                        ;stop

                        ;read, save, prompt='save these screen images??? <1> YES : <99> Stop'
                                 ;save=1
                                 ;save=fix(save)

                        tvlct, r, g, b, /get

                        ;if save eq 99 then stop
;stop
                        if saveimgs eq 1 then begin
                        ;wset, 2 & DFFimg=tvrd(/true) & write_tiff, strtrim(strcompress(iroot+outname+'_DFF.tif', /remove_all), 2), byte(DFFimg);, red=r, blue=b, green=g
                        wset, 2 & DFFimg=tvrd(/true) & write_tiff, iroot+outname+'_DFF.tif', byte(DFFimg)
                        ;write_tiff, strtrim(strcompress(iroot+outname+'_DFF.tif', /remove_all), 2), byte(rr*60+5>1), red=r, blue=b, green=g
                        ;wset, 1 & SD2img=tvrd(/true) & write_tiff, strtrim(strcompress(iroot+outname+'_SD2img.tif',/remove_all), 2), byte(SD2img);, red=r, blue=b, green=g
                        
                        wset, 1 & SD2img=tvrd(/true) & write_tiff, iroot+outname+'_SD2img.tif', byte(SD2img)
                        ;write_tiff, strtrim(strcompress(iroot+outname+'_SD2img.tif',/remove_all), 2), bytscl(ime), red=r, blue=b, green=g

                        ;wset, 3 & im=tvrd(/true ) & write_tiff, strtrim(strcompress(iroot+outname+'_ccrit.tif',/remove_all), 2), byte(SD3img);, red=r, blue=b, green=g
                        wset, 3 & im=tvrd(/true ) & write_tiff, iroot+outname+'_ccrit.tif', byte(SD3img)
                        
                        ;wset, 0 & Normimg=tvrd(/true) & write_tiff, strtrim(strcompress(iroot+outname+'_Norm_2SD_img.tif',/remove_all), 2), byte(Normimg);, red=r, blue=b, green=g
                        wset, 0 & Normimg=tvrd(/true) & write_tiff, iroot+outname+'_Norm_2SD_img.tif', byte(Normimg)
                        
                        write_tiff, iroot+outname+'_F.tif', byte(imb), red=r, blue=b, green=g
                        ;write_tiff, strtrim(strcompress(iroot+outname+'_F.tif',/remove_all), 2), byte(imb), red=r, blue=b, green=g
;stop
                        wset, 2 &  tv, img
                        ;write_tiff, strtrim(strcompress(fname+string(j)+'overlay.tif',/remove_all), 2), byte(img),red=r, blue=b, green=g
                        endif



                        ;stop




                        jump1:

                        ime=imf
                        ;***************** this is the begining of the loop to identify sparks in 1 image


                        While (total(im) ne 0) do begin
                        a=min(where(im eq 1))
                        tt=a/nx & xx=a mod nx

                            ;define search area
                        nnl=fix(min([4./pixsz, xx]))
                        nnr=fix(min([4./pixsz, nx-xx-1]))
                        mmb=fix(min([250./sspeed, tt]))
                        mme=fix(min([250./sspeed, ny-tt-1]))
                        ym=bytarr(nnl+nnr+1,mmb+mme+1)  ;array to hold growing points
                        ym(nnl,mmb)=1             ;initial seeding for growth
                        sk=ym                  ;array to hold  the spark as seen in the 2SD image

                        ;stop


                        ;*******************This loop finds contiguous points in each spark
                        for iii=0,500 do begin          ;surface growth generation count
                        yt=sk
                                 ;potential new surface points
                        ;yn=(ime(xx-nnl:xx+nnr,tt-mmb:tt+mme)) and fix(smooth(float(ym),3,edge=1)*100.<1)
                        yn=(ime(xx-nnl:xx+nnr,tt-mmb:tt+mme)) and fix(smooth(float(ym),3,/edge_truncate)*100.<1) ; changed to add the edge truncate command
                                            ;dilation of ym by 3*3 filter
                        if total(yn) eq 0 then goto, jump2
                                            ;no furth growth, stop and update spark count
                                            ;  **** jump out of iii loop

                        sk=sk>yn                 ;update sk
                        ime(xx-nnl:xx+nnr,tt-mmb:tt+mme)=ime(xx-nnl:xx+nnr,tt-mmb:tt+mme)-fix(sk)>0
                                            ;excise the points that already included in the cluster
                        im(xx-nnl:xx+nnr,tt-mmb:tt+mme)=im(xx-nnl:xx+nnr,tt-mmb:tt+mme)-fix(sk)>0
                        ym=yn-fix(yt)>0           ;true new surface growth point
                        if iii eq 500 then print,'WARNING: SPARK SEARCH AREA MAY BE TOO SMALL'


                        endfor



                        jump2:  ; *** jump here when no further growth of spark




                        a=rebin(float(sk),1,mmb+mme+1)&a=where(a ne 0)
                          ;at,bt: initiation and end time
                        at=max(a)&bt=min(a)
                        if at eq bt then at=at+1

                        a=rebin(float(sk),nnl+nnr+1,1)&a=where(a ne 0)
                          ;ax,bx: left and right edges
                        ax=max(a)&bx=min(a)
                        if ax eq bx then ax=ax+1


                        ;stop ;*******

                        wset,0&z=tvrd(xx-nnl,tt-mmb,nnl+nnr+1,mmb+mme+1)& if drawbox eq 1 then tv,z*(1-(sk))+(sk)*220,xx-nnl,tt-mmb
                        ytt=ym*0&ytt([bx,ax],bt:at)=1&ytt(bx:ax,[bt,at])=1
                        ;stop;******
                        wset,0&z=tvrd(xx-nnl,tt-mmb,nnl+nnr+1,mmb+mme+1)&if drawbox eq 1 then tv,z*(1-(ytt))+(ytt)*221,xx-nnl,tt-mmb

                        wset,1 & sd2z=tvrd(xx-nnl+bx,tt-1, ax-bx, at-bt)
                        sd2z=where(sd2z ge 1.) ; define area for the 2SD portion of the df/f image


                        if n_elements(sd2z) le 9  then goto,  skipspark  ; if the box is to small,  skip over the spark

                        skc=skc+1  ; **** increment the spark counter
                        if ccrit eq ccrit then tskc=tskc+1  ;*** increment the total spark counter




                        ;*************  extracting parameters  for a box just selected
                        ;if ccrit eq 2 then begin
                        if ccrit eq ccrit then begin
                               cxx=xx   ;** x center
                               cyy=tt         ;** y center
                               bxx=ax-bx     ;** x width
                               byy=at-bt     ;** y duration


                        ; *** the DF/F spark data in the box

                        if drawbox eq 1 then begin
                               wset, 1 & draw_box, xx-nnl+bx, tt-1, ax-bx, at-bt, !d.n_colors-255 ; *** this is the code for the draw box
                               wset, 2 & draw_box, xx-nnl+bx, tt-1, ax-bx, at-bt, !d.n_colors-255
                               wset, 3 & draw_box, xx-nnl+bx, tt-1, ax-bx, at-bt, !d.n_colors-255
                        endif

                        print, 'spark counter,   total spark counter', skc, tskc

                        ;stop

                        ;read, save, prompt='save these screen images??? <1> YES  <99> Stop '
                                 ;save=1
                                 ;save=fix(save)
                        tvlct, r, g, b, /get
                        if saveimgs eq 1 then begin
                          wset, 2 & DFFimg=tvrd(/true) & write_tiff, iroot+outname+'_DFF.tif', byte(DFFimg), red=r, blue=b, green=g
                          wset, 1 & SD2img=tvrd(/true) & write_tiff, iroot+outname+'_SD2img.tif', byte(SD2img), red=r, blue=b, green=g
                          wset, 3 & SD3img=tvrd(/true) & write_tiff, iroot+outname+'_ccrit.tif', byte(SD3img), red=r, blue=b, green=g
                          wset, 0 & Normimg=tvrd(/true) & write_tiff, iroot+outname+'_Norm_2SD_img.tif', byte(Normimg), red=r, blue=b, green=g
                          write_tiff, iroot+outname+'_F.tif', byte(imb);, red=r, blue=b, green=g
                        
                        
                        ;wset, 2 & DFFimg=tvrd(/true) & write_tiff, strtrim(strcompress(iroot+outname+'_DFF.tif', /remove_all), 2), byte(DFFimg), red=r, blue=b, green=g
                        ;wset, 1 & SD2img=tvrd(/true) & write_tiff, strtrim(strcompress(iroot+outname+'_SD2img.tif',/remove_all), 2), byte(SD2img), red=r, blue=b, green=g
                        ;wset, 3 & SD3img=tvrd(/true) & write_tiff, strtrim(strcompress(iroot+outname+'_ccrit.tif',/remove_all), 2), byte(SD3img), red=r, blue=b, green=g
                        ;wset, 0 & Normimg=tvrd(/true) & write_tiff, strtrim(strcompress(iroot+outname+'_Norm_2SD_img.tif',/remove_all), 2), byte(Normimg), red=r, blue=b, green=g
                        ;write_tiff, strtrim(strcompress(iroot+outname+'_F.tif',/remove_all), 2), byte(imb);, red=r, blue=b, green=g
                        ;wset, 2 &  tv, img
                        ;write_tiff, strtrim(strcompress(fname+string(j)+'overlay.tif',/remove_all), 2), byte(img),red=r, blue=b, green=g
                        endif



                        ;stop

                                  sarr=rr(xx-nnl+bx:xx-nnl+bx+(ax-bx), tt-1:tt-1+(at-bt))
                                  sd2_area=where(sk ge 1.0) ; extract pixels in SD*2 area

                                  sd2meanamp= total(sarr(sd2z))/n_elements(sd2z); mean amplitude of the 2sd potion of the DF/F box


                                   z=where(sarr ge 0.5*(max(sarr)))
;stop


                            ;********  make the image overlay for the 2SD, the area at 50%DFF and the crit

                                      ;  create a template image of the dff area and then make an array of the points of DFF area >50% max.
                                      imgdff=fltarr(ix, iy) & subimg=imgdff(xx-nnl+bx:xx-nnl+bx+(ax-bx), tt-1:tt-1+(at-bt)) & subimg(z(*,*))=200 & imgdff(xx-nnl+bx:xx-nnl+bx+(ax-bx), tt-1:tt-1+(at-bt))=subimg(*,*)
                                      dffarea=where(imgdff ge 1)

                                      ; SD2 and crit overlay
                                      wset, 3 & SD3img=tvrd(/true)
                                      img1=fltarr(ix, iy) & sd2=where(imf ge 0.1) & xsd2=where(imf le 0.1) & sd25=where(sd3img ge 0.1) & img1(sd2(*,*))=155.& img1(sd25(*,*))=100.

                                      ; SD2 and DFFamp
                                      img3=fltarr(ix, iy) & img3(sd2(*,*))=155.& img3(dffarea(*,*))=200 & img3(xsd2(*,*))=0.0

                                      ;  SD@ and DFFamp and crit
                                      img4=img3 & img4(sd25(*,*))=100. ; all three parameters overlayed

                                      tvlct, r, g, b, /get
                                           ;write_tiff, strtrim(strcompress(iroot+iname+'_2SD_crit_overlay.tif',/remove_all), 2), byte(img1),red=r, blue=b, green=g
                                           write_tiff, iroot+iname+'_2SD_crit_overlay.tif', byte(img1),red=r, blue=b, green=g
                                           ;write_tiff, strtrim(strcompress(iroot+iname+'_2SD_Amp_overlay.tif',/remove_all), 2), byte(img3),red=r, blue=b, green=g
                                           write_tiff, iroot+iname+'_2SD_Amp_overlay.tif', byte(img3),red=r, blue=b, green=g
                                           ;write_tiff, strtrim(strcompress(iroot+iname+'_3param_overlay.tif',/remove_all), 2), byte(img4),red=r, blue=b, green=g
                                           write_tiff, iroot+iname+'_3param_overlay.tif', byte(img4),red=r, blue=b, green=g


                                   mamp50=total(sarr(z))/n_elements(z)

                                   area=fix(n_elements(z)) ; area at 50% of max amp in DF/F box
                                      print, 'Area at 50% of peak in box is',      area*pixsz,'�m'


;stop
                                   sarramp=max(sarr)
                                      print, 'Peak DF/F in box is ',               sarramp,'DF/F', '     ', n_elements(sarr)


                                   print, 'Mean Amp at 50% of peaK in box is', mamp50 ,'DF/F'


                                     print, 'area at the 2*SD is ',                n_elements(sd2_area)*pixsz,'�m','   ', n_elements(sd2_area),'elements'

                                    print, 'Mean Amp of DF/F at 2*SD  ',         sd2meanamp,'DF/F'
;stop
                 sz=size(sarr)
                        window,5,xs=sz(1),ys=sz(2),xpos=500,ypos=0
                        device, retain=2, decomposed=0

                        color_lut, 1
                        ;stop
                        wset,5 & tv,sarr;*50+10>1
                        spkdff=sarr;*50+10>1
                        write_tiff, iroot+outname+'_spkDFF.tif', float(spkdff), red=rd, blue=bl, green=gr
                        if ClrLut eq 1 then loadct,0 else color_lut, 1
;stop

                                  ;wset, 0 & tvscl, sarr

                        ;*** end of DF/F sparks data

                        ; correct image  number (j) if less than 10
                            ;if j lt 10 then flname=ff(imgno)+'_0'+strtrim(string(j), 1)
                            ;if j ge 10 then flname=ff(imgno)+'_'+strtrim(string(j), 1)
                        ;stop
                            ;flname=strtrim(strmid(flname, rstrpos(flname, '\')+1),2)
                            rnamea(tskc-1) = fname
                            fiberflour=total(avgmsk(fiberaoi))/n_elements(avgmsk(fiberaoi))
                    FiberFlrA(j)=fiberflour

                        ;print, area*pixsz,sarramp, mamp50, sd2_area*pixsz, sd2meanamp

                        ;stop
                            ;rnamea(tskc-1)= strtrim(strcompress(iroot+outname),2)
                            rnamea(tskc-1)= iroot+outname
                            fiberfloura(tskc-1) = fiberflour
                            numa(tskc-1) = skc
                            imgnuma(tskc-1) = curimg
                            boxnoa(tskc-1) = skc
                            cxa(tskc-1) = xx-nnl+bx ;  this is the lower left x box coordinate
                            cya(tskc-1) = tt-1 ; this is the lower left y box coordinate
                            bxa(tskc-1) = ax-bx ; this is the x length
                            bya(tskc-1) = at-bt  ; this is the y length
                            xxa(tskc-1) = xx    ;this is the x center of the sarr array
                            tta(tskc-1) = tt    ;this is the y center of the sarr array
                            DFF50area(tskc-1)=  area*(pixsz*pixsz)
                            DFFamp(tskc-1) =  sarramp
                            meanDFFamp(tskc-1) = mamp50
                            SD2area(tskc-1) = sd2_area*(pixsz*pixsz)
                            SD2DFFamp(tskc-1) =sd2meanamp
                            naoi(tskc-1) = 0



                            skipspark:

                        endif

                        ;naoi=naoi+1
                        boxno=boxno+1  ; increment box number
                        endwhile  ; ********* end of loop to accept the sparks

                        ; *** this is where we will put the congruent box check

                        jump4: print, 'number of sparks detected at cri=', cri, 'is', skc



                   endfor    ;******* end of loop (#4) to take care of spark areas over the number of images in each pic file.


                    if k eq 0 then spkcnt=0

                    tvlct, r, g, b, /get


                    spkcnt=tskc

                   ;******************   out saving routine************************88

                    ;stop
                           save_text_file:     ; routine to terminate and save the spark info file

         ;print, 'writing file   '+strcompress(iroot+iname+'_FiberF.txt', /remove_all);text file of Fiber F for each img within the stack
       ;prevFilesav=file_search(strcompress(iroot+iname+'_*_FiberF.txt', /remove_all)) & Fvnum=n_elements(prevFilesav) &  if prevFilesav eq '' then Fvnum=0
       print, 'writing file   '+iroot+iname+'_FiberF.txt';text file of Fiber F for each img within the stack
       prevFilesav=file_search(iroot+iname+'_*_FiberF.txt') & Fvnum=n_elements(prevFilesav) &  if prevFilesav eq '' then Fvnum=0
;stop

       ;openw,u, strcompress(iroot+iname+'_v'+string(fvnum)+'_FiberF.txt', /remove_all), /get_lun
       openw,u, iroot+iname+'_v'+string(fvnum)+'_FiberF.txt', /get_lun
         printf,u,format='(a90,3a10,a20)','File','Frame #','Fiber F','BKG F','Corrected Fiber F'

          for g= 0, n-1 do begin
              framenum=firstimg+1+g
              printf,u,iroot+iname,framenum, FiberFlra(g),imgBKG(G),fiberfimg(g),FORMAT='(a90,i10,f10.2,f10.3,f20.3)
          endfor

         close,u
         free_lun,u

              print, 'writing file   '+iroot+iname+'_data.txt'

                    ;prevsav=file_search(strcompress(iroot+iname+'_*_data.txt', /remove_all)) & vnum=n_elements(prevsav) & if prevsav eq '' then vnum=0
                    prevsav=file_search(iroot+iname+'_*_data.txt') & vnum=n_elements(prevsav) & if prevsav eq '' then vnum=0
                      openw, v, iroot+iname+'_v'+string(vnum)+'_data.txt', /get_lun
                               ;openw, v, strcompress(iroot+iname+'_v'+string(vnum)+'_data.txt', /remove_all), /get_lun

if scope eq 3 then begin
    pixdepth='ND'
    framerate='ND'
    pixdwell='ND'
endif

                                                  printf,v, '# of frames =              ',numframes, firstimg+1, lastimg
                                                  printf,v, 'x size =                   ',ix
                                                  printf,v, 'y size =                   ',iy
                                                  printf,v, 'byte depth =               ',pixdepth
                                                  printf,v, 'pixel size =               ',pixsz
                                                  printf,v, 'Scan Speed =               ',sspeed
                                                  printf,v, 'FrameRate=                 ', framerate
                                                  printf,v, 'Pixel_dwell =              ',pixdwell
                                                  printf,v, 'Fiber area (�m2)=          ',(n_elements(fiberaoi))*(psize*psize)
                                                  printf,v,'Fiber Fluorescence =        ',total(avgmsk(fiberaoi))/n_elements(avgmsk(fiberaoi))
                                                  printf,v,'AVG background Flouresence= ', mean(avgimg(bkgaoi))
                                                  printf,v,'Fiber F-BKG F=             ',(total(avgmsk(fiberaoi))/n_elements(avgmsk(fiberaoi)))- (mean(avgimg(bkgaoi)))
                                                  printf,v, 'criterion value =     ',Cri
                                                  printf,v,'Mult SD value =       ',multsd
                                                  printf,v, ''

                          if tskc ge 1 then begin
                            numa = indgen(tskc)

                              tskc=tskc-1
                                   rnamea = rnamea(0:tskc)
                                   fiberfloura=fiberfloura(0:tskc)
                                   numa = numa(0:tskc)
                                   imgnuma = imgnuma(0:tskc)
                                   boxnoa = boxnoa(0:tskc)
                                   cxa = cxa(0:tskc)
                                   cya = cya(0:tskc)
                                   bxa = bxa(0:tskc)
                                   bya = bya(0:tskc)
                                   xxa = xxa(0:tskc)
                                   tta = tta(0:tskc)
                                   DFF50area = DFF50area(0:tskc)
                                   DFFamp = DFFamp(0:tskc)
                                   meanDFFamp = meanDFFamp(0:tskc)
                                   SD2area = SD2area(0:tskc)
                                   SD2DFFamp = SD2DFFamp(0:tskc)
                                   ;naoi=tskc

;***********************************ADDED 06/08/2005 BY GEORGE RODNEY & CHRIS WARD***************************************************
;This section determines whether an event occurs in sequential images (long event) and/or non-seqential images (hot-spot)
;by calculating the distance between each possible event.  The user inputs a spatial parameter (above at beginning of program) to
;exclude events that are not;close to each other spatially.   A large array is formed giving the x,y coordinates and image # of the
;first event with the x,y coord and image # of the second event along with the distance between these two events (in microns).
;***********************************************************************************************************************************
if rptevts ne 0 then begin
                   sxxa=xxa(sort(xxa)); sort the xxa array and list in ascending order
                   stta=tta(sort(xxa)); sort the tta aray in the same order as xxa to keep the centers of the events together
                   clustarr=fltarr(2,n_elements(xxa))
             spatial=rptevts/pixsz; spatial criterion in microns to indicate that the evts repeats
          if n_elements(xxa) gt 2 then begin
                   for g=0, n_elements(xxa)-1 do begin
                        clustarr(0,g)=xxa(g) & clustarr(1,g)=tta(g);create the sorted 2-D array of x and y locations

                   endfor


          endif

          imzarr=fltarr(ix,iy)

          seqevts=fltarr(3,n_elements(clustarr))

          prevevtrpt=file_search(iroot+iname+'_*_evtrpt.txt') & evnum=n_elements(prevevtrpt) & if prevevtrpt eq '' then evnum=0
          ;prevevtrpt=file_search(strcompress(iroot+iname+'_*_evtrpt.txt', /remove_all)) & evnum=n_elements(prevevtrpt) & if prevevtrpt eq '' then evnum=0

;stop
          rptfile=iroot+iname+'_v'+string(evnum)+'_evtrpt.txt'
          rptfile=strcompress(rptfile,/remove_all)
          
          openw,u, iroot+iname+'_v'+string(evnum)+'_evtrpt.txt', /get_lun
          printf,u,format='(a90,3a10,a40)','File','x-ctr','y-ctr','Reapt-#', 'img_nums'

                progressBar = Obj_New("PROGRESSBAR",Title='Be Patient  (:',/NOCANCEL)
                     progressBar->Start

           for ggg=0,tskc-1 do begin;loop through the array of identified sparks using the x and y centers to look for repeat events

                    evtctr=bytarr(n)
              spotvalsum=0
              evtctrimg=0
              spot=clustarr(0:1,ggg); get the x and y center for spark(ggg)
               spotx=spot(0) & spoty=spot(1); split the coordinates into separate x and y variables
               seqevts(0,ggg)=spotx & seqevts(1,ggg)=spoty
              if spotx ne 0 and spoty ne 0 then begin
                 for gg = 0, n-1 do begin ;loop through the number of images within the stack
                   imzarr(*,*)=imz(*,*,gg); get the gg image withing the stack
                   spotval=imzarr(spotx-spatial/2:spotx+spatial/2,spoty-spatial/2:spoty+spatial/2);this is an array of area determined by user input as
                                                                                         ;the area to look for repeat events. This area will be used for all n images to determine if an event is repeated
                   spotvalm=max(spotval);get the max value of the area of interest, should be 1 if an event occurs here
                      if spotvalm gt 0 then evtctr(gg)=1
                   spotvals=spotvalm
                   spotvalsum=spotvalsum+spotvals

                  endfor
               seqevts(2,ggg)=spotvalsum
               evtctrimg=where(evtctr eq 1)

               print,seqevts(0:2,ggg),evtctrimg,format='(2f10.2,i8,40i)'
               printf,u,rptfile,seqevts(0,ggg),seqevts(1,ggg),seqevts(2,ggg),evtctrimg,format='(a90,3i10,40i)'
                 counter=(ggg+1)*100/tskc
                   progressBar->Update, counter
               endif
            endfor

              progressBar->Destroy
              Obj_Destroy, progressBar

          close,u
          free_lun,u
          Obj_Destroy, progressBar
              allevtarr=fltarr(ix,iy)
                 for gg= 0, n-1 do begin
                   allevtarr(*,*)=allevtarr(*,*)+imz(*,*,gg)
                   ;draw_box,cxa(gg),cya(gg),bxa(gg),bya(gg)
                 endfor
              wset,3 & tvscl,allevtarr
              wset,3 & allevtstiff=tvrd() & write_tiff,iroot+iname+'_allevts.tif', byte(allevtstiff)

endif
;stop
                               printf, v,  'filename', 'Fiber_F', 'spk#', 'pic#', 'box#','cxa', 'cya','bxa','bya','x-cent', 'y-cent','DFFamp', 'DFF50area', 'mnDFFamp', 'sd2area', 'DFFampSD2', format='(a90,a12,7a5,7a10)'
                                   for i = 0, tskc do begin          ; read in the list of sparks

                                      rnamer = strcompress(rnamea(i),/remove_all)
                                         printf, v, rnamea(i), fiberfloura(i), numa(i),imgnuma(i),boxnoa(i),cxa(i),cya(i),bxa(i),bya(i),xxa(i),tta(i), $
                                             DFFamp(i), DFF50area(i), meanDFFamp(i), SD2area(i), SD2DFFamp(i), format='(a90,f12.4,7i5,7f10.4)'
                                   endfor
                               close, v & free_lun, v

                               print, 'Total number of sparks identified is', tskc
                               
                               openw,v, iroot+iname+'_v'+string(evnum)+'_mapallspks.txt', /get_lun
                               printf,v,format='(a90,7a5)','File','cxa', 'cya','bxa','bya','x-cent', 'y-cent'
                               ;stop
                                for i= 0, tskc do begin
                                    printf,v,strcompress(rnamea(i),/remove_all),cxa(i),cya(i),bxa(i),bya(i),xxa(i),tta(i),format='(a90,7i5)'
                                endfor
                                close, v & free_lun, v
                               
                               tskc=0
                          endif



                   eol:  ; escaped from the "rejected" file option



        write_tiff, iroot+iname+'_avgimg.tif', byte(reform(avgimg))
        write_tiff, iroot+iname+'_avgFimg_msk.tif', byte(avgmsk)

endfor ;*** end of big loop (K) which loops over the number of pic files found in the path











;stop
done:

close, /all, /file

wdelete,1
wdelete,2
wdelete,3
wdelete,0
wdelete,5

Print,'You have finished the selected files!'
;Bell,/speedy1
end
