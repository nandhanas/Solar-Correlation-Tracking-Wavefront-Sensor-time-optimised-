FUNCTION sushift, img,siftx,sifty

;NAME:
;		sushift 
;PURPOSE:
;		Shifting an image by fraction of a pixel
;SAMPLE CALLING SEQUENCE:
;		result = sushift(img,siftx,sifty) intentionally spelt wrongly (siftx,sifty)
;INPUTS:        img: 2-D image
;               siftx: desired shift in x direction
;               sifty: desired shift in y direction
;
;OUTPUTS:       img1 : shifted image
;                    
;
;RETURNS:       img1: 2-D image.
;
;
;HISTORY:       Written by R. Sridharan ; got the idea from Dr. Srikant, IIA.

!except=2
t = systime(1)

sz = size(img)
fimg = shift(fft(shift(img,sz(1)/2,sz(2)/2),-1,double=3),sz(1)/2,sz(2)/2)
fimg1 = dcomplexarr(sz(1),sz(2))
for i = 0,sz(1)-1 do begin
  for j = 0,sz(2)-1 do begin
     fimg1(i,j)  = fimg(i,j)*exp(complex(0,1)*(2*!pi*$
       ( ((i-(sz(1)/2.))*(-siftx))/sz(1)+((j-(sz(2)/2.))*(-sifty))/sz(2) ) ))
  endfor
endfor
img2 = float(shift(fft(shift(fimg1,sz(1)/2,sz(2)/2),1,double=3),sz(1)/2,sz(2)/2))
return, img2
end
