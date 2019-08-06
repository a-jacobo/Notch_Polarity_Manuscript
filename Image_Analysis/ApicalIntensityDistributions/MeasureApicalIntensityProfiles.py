#@File(label = "Input directory",value = '/Volumes/Code/CodeForJDESH/ApicalIntensityDistributions/Examples/080818_TL_vangl2_spectrin_4dpf/',style ="directory",persist=false) inputpath
#@Integer(label = "Reference wavelength for detecting apical surface and cells", value = 488) ref_wavelength
#@Integer(label = "Noise tolerance threshold for detecting cells (5000 for TL and Notch mutant)", value = 200) noise_threshold
#@Integer(label = "Radius of apical surface (in pixels)", value = 13) radius

from ij import IJ, ImagePlus, ImageStack
from ij.plugin import Converter, HyperStackConverter, ZProjector, OverlayCommands
from ij.plugin import ZProjector, RGBStackMerge, RGBStackConverter, ChannelSplitter
from ij.process import ImageStatistics, ImageProcessor
from ij.process import FloatProcessor, ByteProcessor, AutoThresholder
from ij.plugin.frame import RoiManager
from array import zeros
from ij.measure import Measurements, ResultsTable
import math
import functools
import re
import csv
import os
import glob
import shutil
from ij.gui import GenericDialog,WaitForUserDialog, Line, Roi, Plot, NonBlockingGenericDialog
import math


def get_apical_channels(image_list,maxChannels):
	wave_lengths = []
	for c in range(1,maxChannels+1):
		w_re =re.compile('(?<=w'+str(c)+'iSIM)\w.*')
		for f in image_list:
			if (w_re.search(f) != None and int(w_re.search(f).group()[0:3]) != 405): #Leave out DAPI
				wave_lengths.append(int(w_re.search(f).group()[0:3]))
				break
	channel_no =  len(wave_lengths)
	print 'Detected '+str(channel_no) + ' channel(s) with wavelength(s) ' + str(wave_lengths)
	return wave_lengths,channel_no

oc = OverlayCommands()
rt = ResultsTable()
measured = rt.getResultsTable()
cs = ChannelSplitter()
cm = RGBStackMerge()
zp = ZProjector()
zp.setMethod(ZProjector.AVG_METHOD)


IJ.run("Close All")

#Find images and wavelengths
imagelist=[f for f in os.listdir(str(inputpath)) if (('.tif' in f.lower()) and ('thumb' not in f.lower()))]
wavelengths,nChannels=get_apical_channels(imagelist,6)
tm_wavelength = [x for x in wavelengths if not x==ref_wavelength][0]

#Sort image files
regexDetect = re.compile("^(?!.*(thumb)).*iSIM"+str(ref_wavelength)+'.*\.(?i)tif') 	#actin or spectrin stainings, will be used to detect apical surfaces and cells
regexMeasure = re.compile("^(?!.*(thumb)).*iSIM"+str(tm_wavelength)+'.*\.(?i)tif') 	#vangl2 or gnai stainings, will be used to measure distribution
detect_filelist = filter(regexDetect.search, os.listdir(str(inputpath)))
measure_filelist = filter(regexMeasure.search, os.listdir(str(inputpath)))

print('Number of neuromasts: '+str(len(measure_filelist)))

outputdir = os.path.join(str(inputpath),'Output')
if not os.path.exists(outputdir):
	os.mkdir(outputdir)

#Initialize empty stacks, these will contain ALL cells from ALL neuromasts
stacks = [ImageStack(60,60),ImageStack(60,60)]
#Initialize datafile
datafile = open(os.path.join(outputdir,'RadialIntensityProfile.csv'), 'w')
datawriter = csv.writer(datafile)
datawriter.writerow(['Angle','Intensity'])

for detect_file in detect_filelist:
	print detect_file
	try:
		sample_label = detect_file.split('w')[0]

		IJ.run("Close All")

		detect_imp=IJ.openImage(os.path.join(str(inputpath),detect_file))
		measure_imp = IJ.openImage(glob.glob(os.path.join(str(inputpath),sample_label+'*iSIM'+str(tm_wavelength)+'-***.TIF'))[0])

		#Find apical surface
		stack = detect_imp.getStack()
		dims = detect_imp.getDimensions()
		avgI = []
		for i in range(1,stack.size()+1): #loop over all images in the stack
			if i < 20:
				ip = stack.getProcessor(i)
				avgI.append(ip.getStats().mean) #get average intensity
			else:
				avgI.append(0)
			if len(avgI)==dims[3]:
				apicalZ = avgI.index(max(avgI)) #get position of apical surface
		if apicalZ==False:
			print('Selecting default apical slice...')
			apicalZ = 9

		#Get apical average projection of apical surface and find cell positions
		try:
			apical_detect_imp = zp.run(detect_imp,'avg',apicalZ-4, apicalZ+4)
			apical_measure_imp = zp.run(measure_imp,'avg',apicalZ-4, apicalZ+4)
		except:
			print('Projecting full apical part of stack...')
			apical_detect_imp = zp.run(detect_imp,'avg',5, stack.size())
			apical_measure_imp = zp.run(measure_imp,'avg',5, stack.size())
		apical_detect_imp.show()
		measured.reset()
		IJ.run("Find Maxima...", "noise="+str(noise_threshold)+" output=List exclude")
		measured = rt.getResultsTable()
		centroidX = measured.getColumn(measured.getColumnIndex('X'))
		centroidY = measured.getColumn(measured.getColumnIndex('Y'))
		print 'Detected '+str(len(centroidX))+' cells'

		apicalimps = [apical_measure_imp,apical_detect_imp]

		for cell in range(0,len(centroidX)):		#Loop over cells
			for ch in range(0,2):			#Loop over channels
				norm = 1/(apicalimps[ch].getProcessor().getStats().mean)*1000	#normalize to mean intensity of image
				apicalimps[ch].setRoi(int(centroidX[cell])-30,int(centroidY[cell])-30,60,60)
				toAdd = apicalimps[ch].crop().getProcessor()
				toAdd.multiply(norm)
				stacks[ch].addSlice(toAdd)
	except:
		print('There was a problem with sample '+sample_label)
		pass

#Merge channels and save image stacks

measureAll = ImagePlus('measure_all',stacks[0])
detectAll = ImagePlus('detect_all',stacks[1])

measureAvg = zp.run(measureAll,'avg')
detectAvg = zp.run(detectAll,'avg')

compositeAvg = cm.mergeChannels([measureAvg,detectAvg],False)  #this we'll save

compositeAll = cm.mergeChannels([measureAll,detectAll],False)

IJ.saveAs(compositeAll,"Tiff", os.path.join(outputdir,'Composite_All.tif'))
print('Saved ' + os.path.join(outputdir,'Composite_All.tif'))
IJ.saveAs(compositeAvg,"Tiff", os.path.join(outputdir,'Composite_AVG.tif'))
print('Saved ' + os.path.join(outputdir,'Composite_AVG.tif'))

#Measure radial profile of measure_avg
profile = []
phis = [x * 2*math.pi/100 for x in range(0, 100)]
for phi in phis:
	rprofile = []
	for r in range(radius - 5,radius + 5):
		rprofile.append(measureAvg.getProcessor().getPixel(30+int(r*math.cos(phi)),30+int(r*math.sin(phi))))
	datawriter.writerow([phi,sum(rprofile)])
	profile.append(sum(rprofile))
datafile.close()
print(profile)
