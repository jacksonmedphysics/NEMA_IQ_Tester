#!/usr/bin/env python




#import os,sys,numpy,dicom, math, datetime, argparse
import os,sys,numpy, math, datetime
try:
    import dicom
except:
    import pydicom as dicom
from scipy.misc import imsave
import matplotlib.pyplot as plt
import csv
from scipy import ndimage
import matplotlib

from scipy.interpolate import UnivariateSpline,splrep,sproot
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

#pyinstaller fixes
import FileDialog
import Tkinter

current_dir=os.path.dirname(os.path.realpath(__file__))

sys.path.append(current_dir)
import series_lister

#parser = argparse.ArgumentParser(description='Perform NEMA Image Quality Measurements on a Directory of PET images: Positive & Negative Contrast, Scatter Correction in Lung insert. Designed for either 4:1 or 8:1 ratio. Will output .pdf report and CSV files of raw ROI Data & contrast results. P Jackson 2015')
#parser.add_argument('-d', '--directory', type=str, help='Directory containing series of PET images, Use "Quotes" if given an error', required=True)
#parser.add_argument('-m', '--mirrored', help='Option to mirror BG ROIs, use any character to select', default=False)
#parser.add_argument('-min_slice', '--min_slice', help='Option to override minimum slice to scan for spheres', default=False)
#parser.add_argument('-max_slice', '--max_slice', help='Option to override maximum slice to scan for spheres', default=False)
#parser.add_argument('-mask_percent', '--mask_percent', help='Set fraction of voxels to use for masking spheres, adjust if detecting regions on margin of phantom. Default=0.0040', default=0.0040)
#parser.add_argument('-resample_size', '--resample_size', help='Set axial voxel width (mm) used when resampling image. Default = 1 mm', default=1)






input_dir=os.path.join(current_dir,'input')

unique_series_uids, series_names, unique_datetimes, series_paths_list, unique_modalities = series_lister.get_paths(input_dir)

for i in range(len(unique_modalities)):
    if unique_modalities[i].lower()=='pt':
        print 'PET Directory Found'
        filelist=series_paths_list[i]
        datetime=unique_datetimes[i]
        series_name=series_names[i]
        dcm=dicom.read_file(series_paths_list[i][0])
        station_name=dcm.StationName
        cal_time=dcm.RadiopharmaceuticalInformationSequence[0].RadiopharmaceuticalStartTime
        half_life=float(dcm.RadiopharmaceuticalInformationSequence[0].RadionuclideHalfLife)
        weight=float(dcm.PatientWeight)
        total_act=float(dcm.RadiopharmaceuticalInformationSequence[0].RadionuclideTotalDose)
        decay_time=(int(datetime.split(' ')[1][0:2])-int(cal_time[0:2]))*3600+(int(datetime.split(' ')[1][2:4])-int(cal_time[2:4]))*60+(int(datetime.split(' ')[1][4:6])-int(cal_time[4:6]))
        decay_factor=2**(-decay_time/half_life)



suv_factor=total_act*decay_factor/(weight*1000)        
                            
        

mirrored=False
min_scan=False
max_scan=False
good_percent=0.0040
resample_size=1.0
axial_voxel_width=1.0


#args = parser.parse_args()
#directory=args.directory

#mirrored=args.mirrored
#max_scan=args.max_slice
#min_scan=args.min_slice
#good_percent=float(args.mask_percent)
#axial_voxel_width=float(args.resample_size)
avw=axial_voxel_width
print mirrored, max_scan, min_scan, good_percent, axial_voxel_width

#directory='/Users/price/Documents/NEMA/20150717-Monash/20150717-NEMA-Monash/data/dicom/NEMA_IQ_PHANTOM_200150717/NEMA_IQ_8_1_20150717_114449_343000/NEMA_IQ_8_1_(AC)_0102'
#directory='/Users/price/Documents/NEMA/20150717-Monash/20150717-NEMA-Monash/data/dicom/NEMA_IQ_PHANTOM_200150717/NEMA_IQ_4_1_20150717_131740_296000/NEMA_IQ_4_1_(AC)_0102'

#output_csv=os.path.join(directory,'IQ_Data_'+str(datetime.datetime.now().day)+'-'+str(datetime.datetime.now().month)+'-'+str(datetime.datetime.now().year)+'.csv')
output_csv=os.path.join(current_dir,'IQ_Data_'+str(station_name)+'_'+str(datetime)+'.csv')
csv_out=open(output_csv,'wb')
wr=csv.writer(csv_out,quoting=csv.QUOTE_ALL)
wr.writerow(['ROI Name', 'Slice Number', 'X position', 'Y position', 'Mean 10mm', 'Mean 13mm', 'Mean 17mm', 'Mean 22mm', 'Mean 28mm', 'Mean 37mm'])
print 'Output written to: ', output_csv

#font = {'size'   : int(14/float(axial_voxel_width))}
font = {'size'   : 15}

matplotlib.rc('font', **font)

#axial_voxel_width=1 #Enter desired Voxel size in mm


#filelist=os.listdir(directory)
iterator=0

mirrortest=False
angle_correcting=False
#mirror_correcting=False
def measure_lung_ROI(x_center, y_center, voxel_width, image_slice):
    x,y=numpy.ogrid[:image_slice.shape[0],:image_slice.shape[1]]
    mask30=numpy.sqrt(((y-x_center)*voxel_width)**2+((x-y_center)*voxel_width)**2)<(30/2)
    mean30=ndimage.measurements.mean(image_slice,mask30)
    return mean30, mask30


def measure_background_ROIs(x_center, y_center, voxel_width, image_slice):
    x,y=numpy.ogrid[:image_slice.shape[0],:image_slice.shape[1]]
    mask10=numpy.sqrt(((y-x_center)*voxel_width)**2+((x-y_center)*voxel_width)**2)<(10/2)
    mean10=ndimage.measurements.mean(image_slice,mask10)

    mask13=numpy.sqrt(((y-x_center)*voxel_width)**2+((x-y_center)*voxel_width)**2)<(13/2)
    mean13=ndimage.measurements.mean(image_slice,mask13)
    
    mask17=numpy.sqrt(((y-x_center)*voxel_width)**2+((x-y_center)*voxel_width)**2)<(17/2)
    mean17=ndimage.measurements.mean(image_slice,mask17)

    mask22=numpy.sqrt(((y-x_center)*voxel_width)**2+((x-y_center)*voxel_width)**2)<(22/2)
    mean22=ndimage.measurements.mean(image_slice,mask22)

    mask28=numpy.sqrt(((y-x_center)*voxel_width)**2+((x-y_center)*voxel_width)**2)<(28/2)
    mean28=ndimage.measurements.mean(image_slice,mask28)

    mask37=numpy.sqrt(((y-x_center)*voxel_width)**2+((x-y_center)*voxel_width)**2)<(37/2)
    mean37=ndimage.measurements.mean(image_slice,mask37)

    
    #plt.imshow(mask10)
    return mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37
    

slice_locations=numpy.array([])
indices=numpy.array([])

instance_list=[]
newlist=[]
for filename in filelist:
	if filename.split('.')[-1]=='IMA' or filename.split('.')[-1]=='ima' or filename.split('.')[-1]=='dcm' or filename.split('.')[-1]=='DCM':
		newlist.append(filename)
		#instance_list.append(int(dicom.read_file(os.path.join(directory,filename)).InstanceNumber))
		instance_list.append(int(dicom.read_file(filename).InstanceNumber))
		print instance_list
zipped=sorted(zip(instance_list,newlist))
print zipped
sorted_list=[x for (y,x) in zipped]
print sorted_list


for filename in sorted_list:
    if filename.split('.')[-1]=='IMA' or filename.split('.')[-1]=='ima' or filename.split('.')[-1]=='dcm' or filename.split('.')[-1]=='DCM':
        
            
        dcm=dicom.read_file(filename)
        image_data=dcm.pixel_array*dcm.RescaleSlope+dcm.RescaleIntercept
        if iterator==0:
            image_volume=image_data
            xdim=float(dcm.PixelSpacing[0])
            ydim=float(dcm.PixelSpacing[1])
            zdim=float(dcm.SliceThickness)
        else:
            image_volume=numpy.dstack((image_volume,image_data))
        iterator+=1
        slice_locations=numpy.append(slice_locations,float(dcm.SliceLocation))
        indices=numpy.append(indices,int(dcm.ImageIndex))



print slice_locations
print indices
key=indices.argsort()
print key

slice_locations2=slice_locations
filelist2=filelist
indices2=indices
for i in range(len(key)):
	slice_locations2[i]=slice_locations[key[i]]
	filelist2[i]=filelist[key[i]]
	indices2[i]=indices[key[i]]
slice_locations=slice_locations2
filelist=filelist2

print slice_locations
print indices
#print filelist


z_widths=numpy.array([])




for i in range(len(slice_locations)-1):
	z_widths=numpy.append(z_widths, abs(slice_locations[i]-slice_locations[i+1]))
true_slice_thickness=z_widths.mean()

plt.figure(figsize=(20,40))
if mirrortest==True:
    image_volume=numpy.fliplr(image_volume)
sums=numpy.array([])       
for i in range(image_volume.shape[2]):
    sums=numpy.append(sums,numpy.sum(image_volume[:,:,i]))

maxes=numpy.array([])       
for i in range(image_volume.shape[2]):
    maxes=numpy.append(maxes,numpy.max(image_volume[:,:,i]))

plt.subplot(421)
maxes_x=numpy.arange(1,image_volume.shape[2]+1,1)
plt.plot(maxes_x,maxes)
plt.xlabel('Slice #')
plt.ylabel('Max Voxel Activity (Bq/ml)')

nslices=float(len(slice_locations))


max_slice=maxes.argmax()
if min_scan==False:
    low_thresh=int(0.15*nslices)  # To account for scanners with very noisy images near edge of FOV
else:
    low_thresh=int(min_scan)
if max_scan==False:  
    high_thresh=int(0.85*nslices)
else:
    high_thresh=int(max_scan)
print "sensitive slice range: ", low_thresh, high_thresh
if max_slice<low_thresh or max_slice>high_thresh:
    print "Max detected outside of range, calculating inner max..."
    max_slice=maxes[low_thresh:high_thresh].argmax()+low_thresh

spline_min=max_slice-10
spline_max=max_slice+10

search=7

spline=UnivariateSpline(maxes_x[max_slice-search:max_slice+search],maxes[max_slice-search:max_slice+search],k=3, s=1)
spline_x=numpy.linspace(maxes_x[max_slice-search],maxes_x[max_slice+search],10000)
spline_y=spline(spline_x)
#plt.plot(spline_x,spline_y)

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

popt,pconv=curve_fit(gaus,maxes_x[max_slice-search:max_slice+search],maxes[max_slice-search:max_slice+search],p0=[maxes.max(),max_slice,5])
plt.plot(spline_x,gaus(spline_x,*popt))
max_gaus=gaus(spline_x,*popt).argmax()
max_slice=int(spline_x[max_gaus])


    
plt.title('Central Slice #: ' + str(max_slice)+', Peak Voxel Activity: '+str(round(maxes.max()/1000,2))+' (kBq/ml)')

print 'Peak Slice: ', str(max_slice)
#plt.plot(maxes)
#plt.show()

#need to find the threshold value to get 
#plt.figure(figsize=(30,10))

#plt.subplot(131)
plt.subplot(423)
plt.title('Axial View of Central Slice')
best_slice=image_volume[:,:,max_slice]
resampled_slice=ndimage.zoom(best_slice,(xdim/axial_voxel_width,ydim/axial_voxel_width),order=0) #resamples image voxels to desired square voxels in axial plane



plt.imshow(resampled_slice, vmin=0, vmax=maxes.max())

found_threshold=False
fraction=1.0
maximum_value=numpy.max(resampled_slice)
#good_percent=0.0025
#good_percent=0.0040


while found_threshold==False:
    tested_value=maximum_value*fraction
    mask_lo=best_slice>tested_value
    current_percent=float(numpy.sum(mask_lo))/float((mask_lo.size))
    if current_percent>good_percent:
        found_threshold=True
        mask_val=tested_value
    fraction=fraction-0.0001


mask=resampled_slice>tested_value
    
print 'Max Voxel Value: ', str(maximum_value)
print 'Threshold: ',str(tested_value)
print 'Percent: ',str(current_percent)
print 'Mask Value: ', str(mask_val)
print 'Number of Thresholded Voxels: ', str(numpy.sum(mask))

#plt.subplot(132)
#plt.subplot(422)
#plt.imshow(mask,cmap=plt.cm.gray)


label_im,nb_labels=ndimage.label(mask)
label_im_lo,nb_labels_lo=ndimage.label(mask_lo)

print 'Number of label regions found: ', nb_labels

#plt.subplot(133)
plt.subplot(422)
plt.imshow(label_im, cmap=plt.cm.spectral)
plt.title('Masks and Centroids for Localisation, Threshold Value: '+str(round(tested_value,1))+' (Bq/ml)')
#plt.subplots_adjust(wspace=0.02, hspace=0.02, top=0.9, bottom=0.1, left=0.1, right=0.8)
#plt.subplots_adjust(wspace=0.02, hspace=0.02)

label_list=[]
label_sizes=[]
for i in range(nb_labels):
    label_list.append(i+1)
    print 'Measuring Size of label: ',i ,
    label_sizes.append(ndimage.measurements.sum(mask_lo,labels=label_im_lo,index=[i+1])[0])
    print label_sizes[-1]

label_list_sorted=list(reversed([x for (y,x) in sorted(zip(label_sizes,label_list))]))
label_sizes_sorted=list(reversed([y for (y,x) in sorted(zip(label_sizes,label_list))]))

print 'Label Sizes: ',label_sizes_sorted

centers_of_mass=[]
counter=1
for i in label_list_sorted:
    if counter < 5:
        print 'Computing center of mass for label # ', i-1,
        centers_of_mass.append(ndimage.measurements.center_of_mass(mask,labels=label_im,index=[i])[0])
        print centers_of_mass[-1]
    counter+=1
#centers_of_mass=ndimage.measurements.center_of_mass(mask,labels=label_im,index=[1,2,3,4])
print "Centers of Mass: ", centers_of_mass
xpoint=numpy.array([])
ypoint=numpy.array([])
for point in centers_of_mass:
    xpoint=numpy.append(xpoint,point[1])
    ypoint=numpy.append(ypoint,point[0])
plt.scatter(xpoint,ypoint)


x_center=(centers_of_mass[0][1]+centers_of_mass[3][1])/2
y_center=(centers_of_mass[0][0]+centers_of_mass[3][0])/2
print 'Central Point: ', x_center, y_center
plt.scatter(x_center,y_center,color='r')

neg_sphere1=((2*x_center-xpoint[1]),2*y_center-ypoint[1])
neg_sphere2=((2*x_center-xpoint[2]),2*y_center-ypoint[2])
print neg_sphere1
plt.scatter(*neg_sphere1,color='g')
plt.scatter(*neg_sphere2,color='g')

#plt.subplot(132)
plt.subplot(422)
plt.scatter(xpoint,ypoint)
plt.scatter(x_center,y_center,color='r')
plt.scatter(*neg_sphere1,color='g')
plt.scatter(*neg_sphere2,color='g')


#plt.subplot(131)
#plt.subplot(421)
#plt.scatter(xpoint,ypoint)
#plt.scatter(x_center,y_center,color='r')
#plt.scatter(*neg_sphere1,color='g')
#plt.scatter(*neg_sphere2,color='g')

mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37=measure_background_ROIs(x_center, y_center, axial_voxel_width, resampled_slice)
print 'Means for circles at image center value: ', mean10, mean13, mean17, mean22, mean28

bg_mask=resampled_slice<0
bg_mask_offcenter=resampled_slice<0

#### Getting ROI data for Hot & Cold Spheres on Central Slice ####
mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37=measure_background_ROIs(neg_sphere1[0], neg_sphere1[1], axial_voxel_width, resampled_slice)
wr.writerow(['Negative Sphere 1', str(max_slice), neg_sphere1[0], neg_sphere1[1], mean10, mean13, mean17, mean22, mean28, mean37])
neg_sphere_1=mean37
bg_mask=numpy.ma.mask_or(bg_mask,mask37)

mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37=measure_background_ROIs(neg_sphere2[0], neg_sphere2[1], axial_voxel_width, resampled_slice)
wr.writerow(['Negative Sphere 2', str(max_slice), neg_sphere2[0], neg_sphere2[1], mean10, mean13, mean17, mean22, mean28, mean37])
neg_sphere_2=mean28
bg_mask=numpy.ma.mask_or(bg_mask,mask28)



mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37=measure_background_ROIs(xpoint[0], ypoint[0], axial_voxel_width, resampled_slice)
wr.writerow(['Positive Sphere 1', str(max_slice), xpoint[0], ypoint[0], mean10, mean13, mean17, mean22, mean28, mean37])
pos_sphere_1=mean22
bg_mask=numpy.ma.mask_or(bg_mask,mask22)

mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37=measure_background_ROIs(xpoint[1], ypoint[1], axial_voxel_width, resampled_slice)
wr.writerow(['Positive Sphere 2', str(max_slice), xpoint[1], ypoint[1], mean10, mean13, mean17, mean22, mean28, mean37])
pos_sphere_2=mean17
bg_mask=numpy.ma.mask_or(bg_mask,mask17)

mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37=measure_background_ROIs(xpoint[2], ypoint[2], axial_voxel_width, resampled_slice)
wr.writerow(['Positive Sphere 3', str(max_slice), xpoint[2], ypoint[2], mean10, mean13, mean17, mean22, mean28, mean37])
pos_sphere_3=mean13
bg_mask=numpy.ma.mask_or(bg_mask,mask13)

mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37=measure_background_ROIs(xpoint[3], ypoint[3], axial_voxel_width, resampled_slice)
wr.writerow(['Positive Sphere 4', str(max_slice), xpoint[3], ypoint[3], mean10, mean13, mean17, mean22, mean28, mean37])
pos_sphere_4=mean10
bg_mask=numpy.ma.mask_or(bg_mask,mask10)



#mean10, mean13, mean17, mean22, mean28, mask10, mask13, mask17, mask22, mask28 = measure_background_ROIs(343, 434, axial_voxel_width, resampled_slice)
#plt.imshow(mask28, alpha=0.4)
#plt.imshow(mask22, alpha=0.4)
#plt.imshow(mask17, alpha=0.4)
#plt.imshow(mask13, alpha=0.4)
#plt.imshow(mask10, alpha=0.4)
#print mean10, mean13, mean17, mean22, mean28
phantom_angle=math.degrees(-math.atan2(neg_sphere1[1]-ypoint[1],neg_sphere1[0]-xpoint[1]))
angle_sphere1to4=math.degrees(-math.atan2(ypoint[0]-ypoint[3],xpoint[0]-xpoint[3]))

#mirrored=False
#if xpoint[0]>xpoint[3]:
#    mirrored=True
#
#if mirrored==True:
#    phantom_angle=phantom_angle-180
#    angle_sphere1to4=angle_sphere1to4 + 60
#
#if mirrored==False:
#    angle_sphere1to4=angle_sphere1to4+120

#phantom_angle=(phantom_angle+360) % 360
#angle_sphere1to4=(angle_sphere1to4) % 360
#print 'Angle from filled sphere 2 to centre: ', phantom_angle
#print 'Angle from filled spheres 1 & 4: ', angle_sphere1to4
#print 'Mirroring BG points? ', mirrored



distances=[83.6115559967, 98.7666931599, 110.6302981583, 115.2079259482, 103.8301091001, 92.7612801398, 88.8528003793, 81.1923877424, 89.4044451475, 112.4353322199, 117.0900211562, 107.8459556288]
angles=[266.82,230.74,211.15,192.41,173.48,152.06,128.29,91.14,38.82,348.26,329.79,311.47]

print 'Mirrored: ',str(mirrored)
if mirrored!=False:
    print 'Mirroring BGs...'
    for i in range(len(angles)):
        angles[i]=(180-angles[i]+360)%360

print angles

#if angle_correcting==True:
#    for i in range(len(angles)):
#        angles[i]=(angles[i]+angle_sphere1to4+360)%360


#print angles




def find_xy(x_center,y_center,axial_voxel_width,distance,angle):
    shiftx=distance*math.cos(math.radians(angle))/axial_voxel_width
    shifty=-distance*math.sin(math.radians(angle))/axial_voxel_width
    x=x_center+shiftx
    y=y_center+shifty
    return x,y

bg_counter=1
means10=numpy.array([])
means13=numpy.array([])
means17=numpy.array([])
means22=numpy.array([])
means28=numpy.array([])
means37=numpy.array([])

for i in range(len(distances)):
    x,y=find_xy(x_center,y_center,axial_voxel_width,distances[i],angles[i])
    mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37=measure_background_ROIs(x, y, axial_voxel_width, resampled_slice)
    bg_mask=numpy.ma.mask_or(bg_mask,mask37)
    bg_mask_offcenter=numpy.ma.mask_or(bg_mask_offcenter,mask37)
    wr.writerow(['Background '+ str(bg_counter), str(max_slice), x, y, mean10, mean13, mean17, mean22, mean28, mean37])
    means10=numpy.append(means10,mean10)
    means13=numpy.append(means13,mean13)
    means17=numpy.append(means17,mean17)
    means22=numpy.append(means22,mean22)
    means28=numpy.append(means28,mean28)
    means37=numpy.append(means37,mean37)
    bg_counter+=1

mean30_center, mask30=measure_lung_ROI(x_center, y_center, axial_voxel_width, resampled_slice)
bg_mask=numpy.ma.mask_or(bg_mask,mask30)
bg_mask_offcenter=numpy.ma.mask_or(bg_mask_offcenter,mask30)
wr.writerow(['Lung Center', str(max_slice), x_center, y_center, '','','','', mean30_center,''])

### Repeating +2cm
slice_plus_2=int(round(max_slice+20/zdim))
resampled_slice=ndimage.zoom(image_volume[:,:,slice_plus_2],(xdim/axial_voxel_width,ydim/axial_voxel_width),order=0)
for i in range(len(distances)):
    x,y=find_xy(x_center,y_center,axial_voxel_width,distances[i],angles[i])
    mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37=measure_background_ROIs(x, y, axial_voxel_width, resampled_slice)
    #bg_mask=numpy.ma.mask_or(bg_mask,mask37)
    wr.writerow(['Background '+ str(bg_counter), str(slice_plus_2), x, y, mean10, mean13, mean17, mean22, mean28, mean37])
    means10=numpy.append(means10,mean10)
    means13=numpy.append(means13,mean13)
    means17=numpy.append(means17,mean17)
    means22=numpy.append(means22,mean22)
    means28=numpy.append(means28,mean28)
    means37=numpy.append(means37,mean37)
    bg_counter+=1

ave_background=means37.mean()
ave_suv=ave_background/suv_factor


plt.subplot(424)
plt.imshow(resampled_slice, vmin=0, vmax=maxes.max())
plt.title('Slice +2cm')
#plt.set_clim(vmin=0, vmax=maximum_value)

mean30_plus2, mask30=measure_lung_ROI(x_center, y_center, axial_voxel_width, resampled_slice)
wr.writerow(['Lung Center +2cm', str(slice_plus_2), x_center, y_center, '','','','', mean30_plus2,''])
#bg_mask=numpy.ma.mask_or(bg_mask,mask30_center)

### Repeating +1cm
slice_plus_1=int(round(max_slice+10/zdim))
resampled_slice=ndimage.zoom(image_volume[:,:,slice_plus_1],(xdim/axial_voxel_width,ydim/axial_voxel_width),order=0)
for i in range(len(distances)):
    x,y=find_xy(x_center,y_center,axial_voxel_width,distances[i],angles[i])
    mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37=measure_background_ROIs(x, y, axial_voxel_width, resampled_slice)
    #bg_mask=numpy.ma.mask_or(bg_mask,mask37)
    wr.writerow(['Background '+ str(bg_counter), str(slice_plus_1), x, y, mean10, mean13, mean17, mean22, mean28, mean37])
    means10=numpy.append(means10,mean10)
    means13=numpy.append(means13,mean13)
    means17=numpy.append(means17,mean17)
    means22=numpy.append(means22,mean22)
    means28=numpy.append(means28,mean28)
    means37=numpy.append(means37,mean37)

    bg_counter+=1
    
plt.subplot(425)
plt.imshow(resampled_slice, vmin=0, vmax=maxes.max())
plt.title('Slice +1cm')

mean30_plus1, mask30=measure_lung_ROI(x_center, y_center, axial_voxel_width, resampled_slice)
wr.writerow(['Lung Center +1cm', str(slice_plus_1), x_center, y_center, '','','','', mean30_plus1,''])

### Repeating -1cm
slice_minus_1=int(round(max_slice-10/zdim))
resampled_slice=ndimage.zoom(image_volume[:,:,slice_minus_1],(xdim/axial_voxel_width,ydim/axial_voxel_width),order=0)
for i in range(len(distances)):
    x,y=find_xy(x_center,y_center,axial_voxel_width,distances[i],angles[i])
    mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37=measure_background_ROIs(x, y, axial_voxel_width, resampled_slice)
    #bg_mask=numpy.ma.mask_or(bg_mask,mask37)
    wr.writerow(['Background '+ str(bg_counter), str(slice_minus_1), x, y, mean10, mean13, mean17, mean22, mean28, mean37])
    means10=numpy.append(means10,mean10)
    means13=numpy.append(means13,mean13)
    means17=numpy.append(means17,mean17)
    means22=numpy.append(means22,mean22)
    means28=numpy.append(means28,mean28)
    means37=numpy.append(means37,mean37)
    bg_counter+=1

plt.subplot(426)
plt.imshow(resampled_slice, vmin=0, vmax=maxes.max())
plt.title('Slice -1cm')

mean30_minus1, mask30=measure_lung_ROI(x_center, y_center, axial_voxel_width, resampled_slice)
wr.writerow(['Lung Center -1cm', str(slice_minus_1), x_center, y_center, '','','','', mean30_minus1,''])

### Repeating -2cm
slice_minus_2=int(round(max_slice-20/zdim))
resampled_slice=ndimage.zoom(image_volume[:,:,slice_minus_2],(xdim/axial_voxel_width,ydim/axial_voxel_width),order=0)
for i in range(len(distances)):
    x,y=find_xy(x_center,y_center,axial_voxel_width,distances[i],angles[i])
    mean10, mean13, mean17, mean22, mean28, mean37, mask10, mask13, mask17, mask22, mask28, mask37=measure_background_ROIs(x, y, axial_voxel_width, resampled_slice)
    #bg_mask=numpy.ma.mask_or(bg_mask,mask37)
    wr.writerow(['Background '+ str(bg_counter), str(slice_minus_2), x, y, mean10, mean13, mean17, mean22, mean28, mean37])
    means10=numpy.append(means10,mean10)
    means13=numpy.append(means13,mean13)
    means17=numpy.append(means17,mean17)
    means22=numpy.append(means22,mean22)
    means28=numpy.append(means28,mean28)
    means37=numpy.append(means37,mean37)
    bg_counter+=1

plt.subplot(427)
plt.imshow(resampled_slice, vmin=0, vmax=maxes.max())
plt.title('Slice -2cm')

mean30_minus2, mask30=measure_lung_ROI(x_center, y_center, axial_voxel_width, resampled_slice)
wr.writerow(['Lung Center -2cm', str(slice_minus_2), x_center, y_center, '','','','', mean30_minus2,''])
plt.subplot(423)
plt.imshow(bg_mask,cmap=plt.cm.gray, alpha=0.3)

plt.subplot(424)
plt.imshow(bg_mask_offcenter,cmap=plt.cm.gray, alpha=0.3)
plt.subplot(425)
plt.imshow(bg_mask_offcenter,cmap=plt.cm.gray, alpha=0.3)
plt.subplot(426)
plt.imshow(bg_mask_offcenter,cmap=plt.cm.gray, alpha=0.3)
plt.subplot(427)
plt.imshow(bg_mask_offcenter,cmap=plt.cm.gray, alpha=0.3)
csv_out.close()


#output_csv=os.path.join(directory,'IQ_Results_'+str(datetime.datetime.now().day)+'-'+str(datetime.datetime.now().month)+'-'+str(datetime.datetime.now().year)+'.csv')
output_csv=os.path.join(current_dir,'IQ_Results_'+str(station_name)+'_'+str(datetime)+'.csv')
csv_out=open(output_csv,'wb')
wr=csv.writer(csv_out,quoting=csv.QUOTE_ALL)
wr.writerow(['ROI Name', 'Diameter mm', 'Average Counts for Sphere', 'Average Background Counts', 'Percent Contrast (4:1)', 'Percent Contrast (8:1)', 'Percent Contrast Cold Sphere', 'Percent Background Variability', 'Relative Error Lung Insert %'])

wr.writerow(['Negative Sphere 1', str(37), neg_sphere_1, means37.mean(), '', '', (1-(neg_sphere_1/means37.mean()))*100, ''])
wr.writerow(['Negative Sphere 2', str(28), neg_sphere_2, means28.mean(), '', '', (1-(neg_sphere_2/means28.mean()))*100, ''])
wr.writerow(['Positive Sphere 1', str(22), pos_sphere_1, means22.mean(), ((pos_sphere_1/means22.mean())-1)/((4)-1)*100, ((pos_sphere_1/means22.mean())-1)/((8)-1)*100, '', '',''])
wr.writerow(['Positive Sphere 2', str(17), pos_sphere_2, means17.mean(), ((pos_sphere_2/means17.mean())-1)/((4)-1)*100, ((pos_sphere_2/means17.mean())-1)/((8)-1)*100, '', '',''])
wr.writerow(['Positive Sphere 3', str(13), pos_sphere_3, means13.mean(), ((pos_sphere_3/means13.mean())-1)/((4)-1)*100, ((pos_sphere_3/means13.mean())-1)/((8)-1)*100, '', '',''])
wr.writerow(['Positive Sphere 4', str(10), pos_sphere_4, means10.mean(), ((pos_sphere_4/means10.mean())-1)/((4)-1)*100, ((pos_sphere_4/means10.mean())-1)/((8)-1)*100, '', '',''])
wr.writerow(['Background Variability 1', str(37), '', '', '', '', '', means37.std()/means37.mean()*100])
wr.writerow(['Background Variability 2', str(28), '', '', '', '', '', means28.std()/means28.mean()*100])
wr.writerow(['Background Variability 3', str(22), '', '', '', '', '', means22.std()/means22.mean()*100])
wr.writerow(['Background Variability 4', str(17), '', '', '', '', '', means17.std()/means17.mean()*100])
wr.writerow(['Background Variability 5', str(13), '', '', '', '', '', means13.std()/means13.mean()*100])
wr.writerow(['Background Variability 6', str(10), '', '', '', '', '', means10.std()/means10.mean()*100])
wr.writerow(['Lung Error Slice +2cm', str(30), mean30_plus2, means37.mean(), '', '', '', '', (mean30_plus2/means37.mean())*100])
wr.writerow(['Lung Error Slice +1cm', str(30), mean30_plus1, means37.mean(), '', '', '', '', (mean30_plus1/means37.mean())*100])
wr.writerow(['Lung Error Slice Centre', str(30), mean30_center, means37.mean(), '', '', '', '', (mean30_center/means37.mean())*100])
wr.writerow(['Lung Error Slice -1cm', str(30), mean30_minus1, means37.mean(), '', '', '', '', (mean30_minus1/means37.mean())*100])
wr.writerow(['Lung Error Slice -2cm', str(30), mean30_minus2, means37.mean(), '', '', '', '', (mean30_minus2/means37.mean())*100])
csv_out.close()
ratio8=False
if ((pos_sphere_1/means22.mean())-1)/((4-1)-1)*100>120:
    print '8:1 Ratio Detected'
    ratio8=True
else:
    print '4:1 Ratio Detected'

#plt.subplot(427)
#plt.imshow(resampled_slice, vmin=0, vmax=maxes.max())

## Generating Coronal View Centered on 17mm Sphere
## Position in resampled slice is xpoint[1],ypoint[1]
## Position in original slice is xpoint[1]*(axial_voxel_width/xdim), ypoint[1]*(axial_voxel_width/ydim)
coronal_y_pos=int(round(ypoint[1]*(axial_voxel_width/ydim)))

coronal_slice=image_volume[coronal_y_pos,:,:]
coronal_resampled=ndimage.zoom(coronal_slice,(xdim/axial_voxel_width,true_slice_thickness/axial_voxel_width),order=0) #resamples image voxels to desired square voxels in axial plane
plt.subplot(428)
plt.imshow(coronal_resampled.T, vmin=0, vmax=maxes.max())
plt.title('Coronal View Centred on 17mm Sphere')

plt.subplot(423)
plt.text(10,int(20/avw),'Contrast Values:',color='w')
plt.text(10,int(40/avw),'37mm Negative Sphere: '+str(round((1-(neg_sphere_1/means37.mean()))*100,1))+' %',color='w')
plt.text(10,int(60/avw),'28mm Negative Sphere: '+str(round((1-(neg_sphere_2/means28.mean()))*100,1))+' %',color='w')
if ratio8==False:
    plt.text(10,int(80/avw),'22mm Positive Sphere: '+str(round(((pos_sphere_1/means22.mean())-1)/((4)-1)*100,1))+' %',color='w')
    plt.text(10,int(100/avw),'17mm Positive Sphere: '+str(round(((pos_sphere_2/means17.mean())-1)/((4)-1)*100,1))+' %',color='w')
    plt.text(10,int(120/avw),'13mm Positive Sphere: '+str(round(((pos_sphere_3/means13.mean())-1)/((4)-1)*100,1))+' %',color='w')
    plt.text(10,int(140/avw),'10mm Positive Sphere: '+str(round(((pos_sphere_4/means10.mean())-1)/((4)-1)*100,1))+' %',color='w')

if ratio8==True:
    plt.text(10,int(80/avw),'22mm Positive Sphere: '+str(round(((pos_sphere_1/means22.mean())-1)/((8)-1)*100,1))+' %',color='w')
    plt.text(10,int(100/avw),'17mm Positive Sphere: '+str(round(((pos_sphere_2/means17.mean())-1)/((8)-1)*100,1))+' %',color='w')
    plt.text(10,int(120/avw),'13mm Positive Sphere: '+str(round(((pos_sphere_3/means13.mean())-1)/((8)-1)*100,1))+' %',color='w')
    plt.text(10,int(140/avw),'10mm Positive Sphere: '+str(round(((pos_sphere_4/means10.mean())-1)/((8)-1)*100,1))+' %',color='w')

plt.text(10,int(520/avw),'Mean Background SUV: '+str(round(ave_suv,4)),color='w')

plt.subplot(424)
plt.text(10,int(20/avw),'37mm Background ROIs SD: '+str(round(means37.std()/means37.mean()*100,2))+' %', color='w')
plt.text(10,int(40/avw),'28mm Background ROIs SD: '+str(round(means28.std()/means28.mean()*100,2))+' %', color='w')
plt.text(10,int(60/avw),'22mm Background ROIs SD: '+str(round(means22.std()/means22.mean()*100,2))+' %', color='w')
plt.text(10,int(80/avw),'17mm Background ROIs SD: '+str(round(means17.std()/means17.mean()*100,2))+' %', color='w')
plt.text(10,int(100/avw),'13mm Background ROIs SD: '+str(round(means13.std()/means13.mean()*100,2))+' %', color='w')
plt.text(10,int(120/avw),'10mm Background ROIs SD: '+str(round(means10.std()/means10.mean()*100,2))+' %', color='w')

try:
    injected_dose=float(dcm.RadiopharmaceuticalInformationSequence[0].RadionuclideTotalDose)/1000000
except:
    injected_dose=0.0
plt.subplot(425)
#plt.text(10,20,'Injected Dose: '+str(round(injected_dose,1))+' MBq',color='w')

#start_time=dcm.RadiopharmaceuticalInformationSequence[0].RadiopharmaceuticalStartTime
#start_hour=float(start_time[0:2])
#start_minute=float(start_time[2:4])
#plt.text(10,40,'Calibration Time: '+str(int(start_hour))+':'+str(int(start_minute)),color='w')
#acquisition_hour=float(dcm.AcquisitionTime[0:2])
#acquisition_minute=float(dcm.AcquisitionTime[2:4])
#plt.text(10,60,'Calibration Time: '+str(int(acquisition_hour))+':'+str(int(acquisition_minute)),color='w')
#time_pi=60*(acquisition_hour-start_hour)+(acquisition_minute-start_minute)
#plt.text(10,80,'Decay Time (min): '+str(int(time_pi)),color='w')
#halflife=float(dcm.RadiopharmaceuticalInformationSequence[0].RadionuclideHalfLife)/60
#decay_factor=2**(-time_pi*halflife)

bg_concentration=means37.mean()
plt.text(10,int(20/avw),'Mean Background Concentration: '+str(round(means37.mean()/1000,2))+' kBq/ml', color='w')
scan_length=abs(slice_locations[0]-slice_locations[-1])
plt.text(10,int(40/avw),'Transaxial Scan Length: '+str(round(scan_length,0))+' mm',color='w')
plt.text(10,int(60/avw),'Axial Step Size: '+str(round(true_slice_thickness,1))+' mm',color='w')
try:
    frame_duration=(float(dcm.ActualFrameDuration)/1000)/60
except:
    frame_duration=0.0
plt.text(10,int(80/avw),'Frame Duration: '+str(round(frame_duration,1))+' minutes',color='w')

#print start_hour, ':', start_minute
try:
    plt.text(10,int(100/avw),'Convolution Kernel: '+dcm.ConvolutionKernel,color='w')
except:
    plt.text(10,int(100/avw),'Convolution Kernel: ',color='w')
try:
    plt.text(10,int(120/avw),'Corrections Applied: '+str(dcm.CorrectedImage),color='w')
except:
    plt.text(10,int(120/avw),'Corrections Applied: ',color='w')    

plt.subplot(426)
try:
    plt.text(10,int(20/avw),'Pixel Spacing (mm): '+str(dcm.PixelSpacing),color='w')
except:
    plt.text(10,int(20/avw),'Pixel Spacing (mm): ',color='w')
plt.text(10,int(40/avw),'Slice Thickness : '+str(round(zdim,1))+' mm',color='w')
plt.text(10,int(60/avw),'Image Matrix Size: '+str(image_volume.shape),color='w')
try:
    plt.text(10,int(80/avw),'Reconstruction Method: '+dcm.ReconstructionMethod,color='w')
except:
    plt.text(10,int(80/avw),'Reconstruction Method: ',color='w')
try:
    plt.text(10,int(100/avw),'Scatter Correction Method: '+dcm.ScatterCorrectionMethod,color='w')
except:
    plt.text(10,int(100/avw),'Scatter Correction Method: ',color='w')

plt.subplot(427)
plt.text(10,int(20/avw),'Lung Insert Relative Error +2cm: '+str(round((mean30_plus2/means37.mean())*100,1))+' %',color='w')
plt.text(10,int(40/avw),'Lung Insert Relative Error +1cm: '+str(round((mean30_plus1/means37.mean())*100,1))+' %',color='w')
plt.text(10,int(60/avw),'Lung Insert Relative Error +0cm: '+str(round((mean30_center/means37.mean())*100,1))+' %',color='w')
plt.text(10,int(80/avw),'Lung Insert Relative Error -1cm: '+str(round((mean30_minus1/means37.mean())*100,1))+' %',color='w')
plt.text(10,int(100/avw),'Lung Insert Relative Error -2cm: '+str(round((mean30_minus2/means37.mean())*100,1))+' %',color='w')



#wr.writerow(['Lung Error Slice +1cm', str(30), mean30_plus1, means37.mean(), '', '', '', '', (mean30_plus1/means37.mean())*100])
#wr.writerow(['Lung Error Slice Centre', str(30), mean30_center, means37.mean(), '', '', '', '', (mean30_center/means37.mean())*100])
#wr.writerow(['Lung Error Slice -1cm', str(30), mean30_minus1, means37.mean(), '', '', '', '', (mean30_minus1/means37.mean())*100])
#wr.writerow(['Lung Error Slice -2cm', str(30), mean30_minus2, means37.mean(), '', '', '', '', (mean30_minus2/means37.mean())*100])
try:
    AcquisitionDate=dcm.AcquisitionDate
except:
    AcquisitionDate=''
plt.tight_layout()

if ratio8==False:
    plt.savefig(os.path.join(current_dir,'IQ_segmented_4_1_'+str(station_name)+'_'+str(datetime)+'.pdf'))
else:
    plt.savefig(os.path.join(current_dir,'IQ_segmented_8_1_'+str(station_name)+'_'+str(datetime)+'.pdf'))


#plt.show()
