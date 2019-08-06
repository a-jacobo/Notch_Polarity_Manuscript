#@File(label = 'Input directory',value = '/Volumes/Code/CodeForJDESH/Segmentation/Examples/Images',persist=false,style ='directory') inputPath
#@File(label = 'Output directory',value = '/Volumes/Code/CodeForJDESH/Segmentation/Examples/',persist=false,style ='directory') outputPath
#@Boolean(label = 'Modify segmentation parameters?', value = True) modify

from ij import IJ, ImagePlus, ImageStack
from ij.plugin import Converter, HyperStackConverter, ZProjector, ChannelSplitter
from ij.process import ImageStatistics, ImageProcessor
from ij.process import FloatProcessor, ByteProcessor, AutoThresholder
from ij.plugin.frame import RoiManager
from array import zeros
from ij.measure import Measurements, ResultsTable
import math
import functools
import os
import sys
from ij.gui import GenericDialog,WaitForUserDialog, Line, Roi, Plot, NonBlockingGenericDialog
import csv

def percentile(N, percent, key=lambda x:x):
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1

def getSegmentationParameters():
  gd = GenericDialog('Enter segmentation parameters')
  gd.addMessage('Slice selection thresholds (percentiles of average intensity): \n \n')
  gd.addNumericField('Lower threshold',0.25,2)
  gd.addNumericField('Upper threshold',0.65,2)
  gd.addMessage('Image processing parameters for cells: \n \n')
  gd.addNumericField('Width of Gaussian filter [px]',4,0)
  gd.addChoice('Auto threshold method', ['Huang','Huang2','Intermodes','IsoData','Li','MaxEntropy','Mean','MinError(I)','Minimum','Moments','Otsu','Percentile','RenyiEntropy','Shanbhag','Triangle','Yen']
,'Triangle')
  gd.addStringField('Particle size range [px^2]','10000-5000000',16)
  gd.addStringField('Particle circularity range','0.20-1.00')
  gd.addMessage('Image processing parameters for nuclei: \n \n')
  gd.addNumericField('Width of Gaussian filter [px]',4,0)
  gd.addChoice('Auto threshold method', ['Huang','Huang2','Intermodes','IsoData','Li','MaxEntropy','Mean','MinError(I)','Minimum','Moments','Otsu','Percentile','RenyiEntropy','Shanbhag','Triangle','Yen']
,'Otsu')
  gd.addStringField('Particle size range [px^2]','3500-10000',16)
  gd.addStringField('Particle circularity range','0.5-1.00')
  gd.addCheckbox('Run in testing mode?',0)
  gd.showDialog()
  if gd.wasCanceled():
    print 'User canceled dialog!'
# Read out inputs
  params = {}
  params['avgI_low']=gd.getNextNumber()
  params['avgI_high']=gd.getNextNumber()
  params['cellsigma'] = int(gd.getNextNumber())
  params['cellmethod'] = gd.getNextChoice()
  params['cellsize'] = gd.getNextString()
  params['cellcircularity'] = gd.getNextString()
  params['nucsigma'] = int(gd.getNextNumber())
  params['nucmethod'] = gd.getNextChoice()
  params['nucsize'] = gd.getNextString()
  params['nuccircularity'] = gd.getNextString()
  params['test'] = gd.getNextBoolean()
  return params



if modify:
    parameters = getSegmentationParameters()
else:
    parameters = {}
    parameters['avgI_low']=0.25
    parameters['avgI_high']=0.65
    parameters['cellsigma'] = 4
    parameters['cellmethod'] = 'Triangle'
    parameters['cellsize'] = '10000-5000000'
    parameters['cellcircularity'] = '0.20-1.00'
    parameters['nucsigma'] = 4
    parameters['nucmethod'] = 'Otsu'
    parameters['nucsize'] = '3500-10000'
    parameters['nuccircularity'] = '0.5-1.00'
    parameters['test']=0


#Class constructors
rm = RoiManager()
cs = ChannelSplitter()
rmi = rm.getInstance()

#Input and output directories
inputPath = str(inputPath)
outputPath = str(outputPath)

imagelist = [f for f in os.listdir(os.path.join(inputPath)) if f.lower()[0].isalnum()]

for sample in imagelist:
	rmi.reset()
	IJ.run('Close All')
	sample_label = sample.split('.')[0]
	outputFile = os.path.join(outputPath,sample_label+'.zip')
	if not os.path.exists(outputFile):	#Check that file does not yet exists
		print 'Segmenting ' + sample_label

		#Load images
		imp = IJ.openImage(os.path.join(inputPath,sample))

		#Split channels
		nucleistack = cs.getChannel(imp,1)
		cellstack = cs.getChannel(imp,2)

		#Select relevant slices by average intensity of cell marker
		avgI = []
		for sliceIndex in range(1,cellstack.size()+1): #loop over slices
			ip = cellstack.getProcessor(sliceIndex)
			avgI.append(ip.getStats().mean) #get average intensity
		avgth1=percentile(sorted(avgI),parameters['avgI_low'])
		avgth2=percentile(sorted(avgI),parameters['avgI_high'])
		sel = [avgI.index(i)+1 for i in avgI if (i >= avgth1 and i <= avgth2)]	#pick slices with intermediate average intensity
		print sel
		#Get ROIs
		rois = []
		for zslice in sel:	#Loop over slices
			#Get ROI from cell channel
			ip1 = cellstack.getProcessor(zslice)
			ip1.resetMinAndMax()
			ip1.blurGaussian(parameters['cellsigma'])
			bp1=ip1.convertToByte(1)
			impBin1 = ImagePlus('BinarizedCells '+str(zslice),bp1)
			impBin1.show()
			IJ.run('Auto Threshold', 'method='+parameters['cellmethod']+' white')
			IJ.run('Analyze Particles...', 'size='+parameters['cellsize']+' circularity='+parameters['cellcircularity']+' exclude add')
			if rmi.getCount()==0:
				print 'No cell ROI on slice '+str(zslice)
				impBin1.changes=False
				impBin1.close()
				continue


			#Get nuclei ROIs
			ip2 = nucleistack.getProcessor(zslice)
			ip2.resetMinAndMax()
			ip2.blurGaussian(parameters['nucsigma'])
			bp2=ip2.convertToByte(1)
			impBin2 = ImagePlus('BinarizedNuclei',bp2)
			impBin2.show()
			rmi.select(impBin2,rmi.getCount()-1)
			IJ.setBackgroundColor(0, 0, 0)
			IJ.run('Clear Outside')
			impBin2.deleteRoi()
			IJ.selectWindow('BinarizedNuclei')
			IJ.run('Auto Threshold', 'method='+parameters['nucmethod']+' ignore_black white')	#threshold the template region
			IJ.run('Watershed')
			#Remove template ROI
			rmi.reset()
			impBin2.killRoi()
			#Add nuclei to ROIs
			IJ.run('Analyze Particles...', 'size='+parameters['nucsize']+' circularity='+parameters['nuccircularity']+' exclude add')
			rois.append([zslice,rmi.getRoisAsArray()])
			print 'Added '+str(len(rmi.getRoisAsArray()))+' ROIs'
			#Close images
			impBin1.changes=False
			impBin1.close()
			impBin2.changes=False
			impBin2.close()

	#View and save ROIs
		imp.show()
		rmi.reset()
		print len(rois)
		for roisonslice in rois:
			imp.setPosition(1,roisonslice[0],1)
			for roi in roisonslice[1]:
				imp.setRoi(roi)
				rmi.runCommand('Add')
		rmi.runCommand('Save',outputFile)
		if parameters['test']:
			gd = NonBlockingGenericDialog('Continue?')
			gd.showDialog()
			if gd.wasCanceled():
				sys.exit('User cancelled script')
