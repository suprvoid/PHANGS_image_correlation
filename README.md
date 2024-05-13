 # Instructions #

- "maptype" is '7.7', 'CO', 'Ha', or '21'. more can be added but fits files are unpacked differently
- "imgtype" is 'original', 'normalized', 'dog', or 'mask'
- "projection" is 'sky' or 'polar'

 to run the analysis for a new galaxy:
 1. create an instance of Galaxy()
 2. fetch the predicted scale data from mega tables using Galaxy.getscales(filepath)
 3. input observation data from fits file using Galaxy.readdata(filepath, maptype)
 4. fetch the galaxy geometry from sample tables using Galaxy.setgeometry(id, filepath)
 5. run the image processing with Galaxy.picture_prop(Rmax, ...)
    - you will have to try different Rmax values for each galaxy, it is the maximum radius to which the image is deprojected
    - try and see whether tlog = True or False preserves more detail. check the results with Galaxy.show(maptype, 'original', 'sky') and Galaxy.show(maptype, 'normalized', 'sky')
    - tlog = True works better for galaxies with bright bulges while False works better for more uniform brightness.
 6. run the measurement with Galaxy.auto_scales(maptype)

 different maps can be cross-correlated. just input more observation data, run the processing and then Galaxy.auto_scales(maptype1, maptype2)
 correlation can be 'linear' (an elementwise product of the image with the offset image)
 or 'spearman'. the latter is much slower, produces a smoother signal and usually gives less self-consistent results. default is 'linear'.
 images can be offset radially when correltating and can produce 2D signals. this is not used in the default measurement method.

 the find_peaks parameters can also be changed. change peak_height to set the minimum height. (the range is 0-1)

 to verify results visually, use Galaxy.plot_peaks_subplots(maptype) to see every correlation signal.
 in my experience when measurements are inconsistent (with predictions and/or themselves) the peaks in the correspoonding signals are usually very small.

 the image processing is saved at different stages:
 1. a log normalization is applied which can be toggled with tlog=False. saved as 'normalized'.
 2. a difference of gaussians algorithm is applied. saved as 'dog'
    the values of the two parameters, dog_sigma1 and dog_sigma2 can be changed to pick out different size features.
 3. a threshold is applied using a hyperbolic tangent function. saved as 'mask'
    the threshold value and slope, mask_thres and mask_grad can be changed but generally don't need to be.
 each of these is deprojected into polar coordinates. 'original' is not deprojected.
 the measurement can be ran with any of these images by changing the imgtype parameter. default is 'mask'



note to self: difference of gaussians acts as a band pass filter for different scale lengths. it might not be useful at all since it will always introduce a bias towards the scales that get through the band pass filter.

running the analysis with just the 'normalized' images should pick out the most dominant scale?
I've disabled all image processing including tlog (deprojection is still needed, obviously)
