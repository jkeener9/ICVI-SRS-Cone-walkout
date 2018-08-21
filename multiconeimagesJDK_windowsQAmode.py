"""
Created on Oct 2017

@author: justin
"""
import glob
import numpy as np
import matplotlib.pyplot as plt
import os.path
import dicom
import pylab
from scipy.interpolate import UnivariateSpline
from matplotlib.backends.backend_pdf import PdfPages

d = "\\\\10.16.72.207\\public\\justin\\__JustinPythonCode__\\20180530_4mm\\"  #Enter Directory where dicom image files are saved
pp = PdfPages(d + "ConeAnalysisResults.pdf")
ConeDiamm=4  #Enter Cone Diameter in mm
BeamCenters=[] 

for k in glob.glob(d + '*.dcm'):
    dcmFiles = []
    dcmFiles=dicom.read_file(k)

    dcmFilespix=dcmFiles.pixel_array
    #plt.imshow(dcmFilespix)
    
    #invert due to acquistion in QA mode 
    orig_array=dcmFilespix
    #dcmFilespix=65535 - orig_array  #apparently does the same function as invert
    dcmFilespix=np.invert(orig_array)
    
    # Normalize the image to the max value in the image. 
    maxpix=np.max(dcmFilespix)
    dcmFilespixNorm=dcmFilespix/maxpix
    #plt.imshow(dcmFilespixNorm)
    
    Width=dcmFilespix.shape[1]
    Height=dcmFilespix.shape[0]
    halfWidth=int(Width/2)
    halfHeight=int(Height/2)
    pixelWidth = dcmFiles.ImagePlanePixelSpacing[1]
    pixelHeight = dcmFiles.ImagePlanePixelSpacing[0]

    SID=dcmFiles.RTImageSID
    ConeDiapix=ConeDiamm/pixelWidth*(SID/1000)   #for now assuming pixelWidth=PixelHeight

    CollAng = str(int(round(dcmFiles.BeamLimitingDeviceAngle)))

    left=int(halfWidth-round(ConeDiapix))
    right=int(halfWidth+round(ConeDiapix))
    top=int(halfHeight-round(ConeDiapix))
    bottom=int(halfHeight+round(ConeDiapix))

    dcmFilespixNormROI=dcmFilespixNorm[top:bottom, left:right]
    maxpixloc=np.where(dcmFilespixNormROI==1)
    maxpixloc=np.array(maxpixloc)
    maxpixloc=maxpixloc[:,0]

    WidthROI=dcmFilespixNormROI.shape[1]
    HeightROI=dcmFilespixNormROI.shape[0]
    y=np.linspace(left,left+WidthROI-1,WidthROI)
    x=np.linspace(top,top+HeightROI-1,HeightROI)

    maxfwhmy=0  #loop through columns
    for i in range(np.asscalar(maxpixloc[1])-round(0.4*ConeDiapix),np.asscalar(maxpixloc[1])+round(0.4*ConeDiapix)):
        splinefit=UnivariateSpline(y,dcmFilespixNormROI[:,i]-np.max(dcmFilespixNormROI[:,i])/2,s=0)  #see sbcrowe.net/python-exercise-field-size-calculator/
        if len(splinefit.roots())==2:  ## have had code crash sometimes due to root not being found correctly
            r1,r2=splinefit.roots()
            fwhm=abs(r1-r2)
            if fwhm > maxfwhmy:
                maxfwhmy=fwhm
                columnfwhm=i
                fwhmr1y=r1
                fwhmr2y=r2
                fwhmyloc=(r1+r2)/2  #location of cone determined as center of peak
   
    plt.figure()
    pylab.plot(y,dcmFilespixNormROI[:,columnfwhm], label="inplane")
    plt.axvspan(fwhmr1y,fwhmr2y,facecolor='b',alpha=0.4)

    maxfwhmx=0   #loop through rows
    for j in range(np.asscalar(maxpixloc[0])-round(0.4*ConeDiapix),np.asscalar(maxpixloc[0])+round(0.4*ConeDiapix)):
        splinefit=UnivariateSpline(x,dcmFilespixNormROI[j,:]-np.max(dcmFilespixNormROI[j,:])/2,s=0)  
        if len(splinefit.roots())==2:
            r1,r2=splinefit.roots()
            fwhm=abs(r1-r2)
            if fwhm > maxfwhmx:
                maxfwhmx=fwhm
                rowfwhm=j
                fwhmr1x=r1
                fwhmr2x=r2
                fwhmxloc=(r1+r2)/2

    pylab.plot(x,dcmFilespixNormROI[rowfwhm,:], label="crossplane") 
    plt.axvspan(fwhmr1x,fwhmr2x,facecolor='r',alpha=0.4)

    plt.title(str(ConeDiamm) + "mm Cone, COLL " + CollAng + "Deg Profiles")
    plt.xlabel("                                                                                     pixels")
    plt.legend()
    plt.grid()

    FWHMcross=(abs(fwhmr1x-fwhmr2x))*pixelWidth*1000/SID
    #print("FWHM crossplane projected to iso = " + str(FWHMcross) +"mm")
    FWHMin=(abs(fwhmr1y-fwhmr2y))*pixelHeight*1000/SID
    #print("FWHM inplane projected to iso = " + str(FWHMin) +"mm")
    captiontext=("FWHM crossplane projected to iso = " + str(FWHMcross) +"mm" + "\n" + "FWHM inplane projected to iso = " + str(FWHMin) +"mm")
    plt.figtext(0,0,captiontext)

    #plt.show()
    pp.savefig()
    plt.close()
    plt.figure()
    plt.imshow(dcmFilespixNormROI, cmap='Greys', extent=(left,right,bottom,top))
    plt.scatter(fwhmxloc,fwhmyloc, color='r', marker='+', s=500)
    plt.title(str(ConeDiamm) + "mm Cone, COLL " + CollAng + "Deg")
    #plt.show()
    pp.savefig()
    plt.close()

    
    BeamCenter=[fwhmxloc,fwhmyloc,int(round(dcmFiles.BeamLimitingDeviceAngle))]
    BeamCenters.append(BeamCenter)
        
npBeamCenters=np.array(BeamCenters)
#plt.figure()
#plt.scatter(npBeamCenters[:,0],npBeamCenters[:,1], color='r', marker='+', s=500)
#plt.xlabel("pixels")
#pp.savefig()
#plt.close()
averagex=np.mean(npBeamCenters[:,0])
averagey=np.mean(npBeamCenters[:,1])

xdist=np.vstack(npBeamCenters[:,0]-averagex)*pixelWidth*(1000/SID)
ydist=np.vstack(npBeamCenters[:,1]-averagey)*pixelWidth*(1000/SID)
CollAngles=np.vstack(npBeamCenters[:,2])
CenterDistfromAvgmm=np.hstack((xdist,ydist,CollAngles))   #still assuming pixel width=height
 
DistancefromAvg=np.sqrt(np.square(CenterDistfromAvgmm[:,0])+np.square(CenterDistfromAvgmm[:,1]))
DistancefromAvg=np.vstack(DistancefromAvg)
AverageDistancefromAvg=np.mean(DistancefromAvg)
MaxDistancefromAvg=np.max(DistancefromAvg)

plt.figure()
plt.scatter(CenterDistfromAvgmm[:,0],CenterDistfromAvgmm[:,1], color='r', marker='+', s=500)
plt.xlabel('                                                                            mm')
plt.ylabel('mm')
plt.title('Cone center coordinates relative to average, projected to isocenter')
captiontext=("Average Distance from Center = " + str(AverageDistancefromAvg) +"mm" + "\n" + "Max Distance from Center = " + str(MaxDistancefromAvg) +"mm")
plt.figtext(0,0,captiontext)
pp.savefig()
plt.close()
        
pp.close()

OutputTable=np.hstack((CollAngles,DistancefromAvg))
print("Collimator Angle & Distance from Center" + "\n" + str(OutputTable))
