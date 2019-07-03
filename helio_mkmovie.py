#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: William Wainwright

This is a small library of functions used make "quick and dirty" movies from helioviewer 
data in the specified aia filters + hmi data, over a specified time range. Requires 
FFMPEG as well as the python packages as follows:
    Sunpy==0.9.9
    datetime
    dateutil
    astropy
    numpy
    PIL
    matplotlib
    os
    subprocess
"""




def get_png (
start='2015/01/17 07:00:00'
,end='2015/01/17 08:30:00'
,wavelength=[94,131,171,193,211,304,335,1600,'hmi']
,basedir='/data/reu/wwainwri/Data/solar/'
,xcenter=941.39
,ycenter=746.31
,width=650
,height=400
,skipdays=0
,skiphours=0):
    
    """Sunpy-Helioviewer client based method to get png files cropped to a 
    region of interest over a given timeframe for any AIA/HMI wavelengths.
    All arguments are by keyword. Script functions in sunpy 0.9.9, it has
    been tested in 1.0.0+, but is potentially slower, buggier, and spams 
    the console. It should still work though. Skipdays and skiphours are 
    intended for use when a large download has ceased in the middle of a 
    time range for a single wavelength and you do not wish to start from
    the beginning. The time skip accounts for the rotation of the sun
    to make sure an object at xcenter and ycenter at start time will be
    in view at the skipped to time.

    Args:
        start (str): Start time string for the download, expected format YYYY/mm/dd HH:MM:SS
        
        end (str): End time string for the download, expected format YYYY/mm/dd HH:MM:SS
        
        wavelength (list): List object of wavelengths, for example [171,304,hmi] <---hmi is str object, wavelengths may be string or int objects
        
        basedir (str): Directory where the data will be downloaded to, subdirectory will be made based on date, expects full /path/to/dir/ format
        
        xcenter (float): Float or int for the center x coordinate of the image, expects pixel coordinates as determined by the test_coords() function
        
        ycenter (float): Float or int for the center y coordinate of the image, expects pixel coordinates as determined by the test_coords() function
        
        width (int): Integer value for the width of the image in pixels
        
        height (int): Integer value for the height of the image in pixels
        
        skipdays (float): Number of days to skip from the start date, intended for use in large downloads where a crash has occurred before the end
        
        skiphours (float): Number of hours to skip from the start date, intended for use in large downloads where a crash has occurred before the end
    """
    
    #Imports
    from sunpy.net import helioviewer
    from datetime import datetime,timedelta
    from dateutil.rrule import rrule, SECONDLY
    import os
    import numpy as np
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from sunpy.coordinates import frames
    
    #Convert xcenter and ycenter pixel coordinates to arcsec offput from solar center coordinates (sunpy makes this unnecessarily complicated)
    arcperpx=np.float128(0.600714)
    xcenter=(xcenter-600)*arcperpx*4096/1200
    ycenter=(ycenter-600)*arcperpx*4096/1200
    
    xcenterinit=float(xcenter)
    
    #Datetime objects
    startdt=datetime.strptime(start,'%Y/%m/%d %H:%M:%S')
    enddt=datetime.strptime(end,'%Y/%m/%d %H:%M:%S')
    endstr=datetime.strftime(enddt,'%Y_%m_%d_%H_%M_%S_')
    
    #Check for bad directory format
    if not basedir[-1] == '/':
        basedir = basedir + '/'  
    
    
    #Sub directory path from date input
    subdir=start.split(' ')[0]+'/'
    
    #Open helioviewer client instance
    hv = helioviewer.HelioviewerClient()


    #akm=np.arctan(1/3600)*1.496e8 <-- The way it SHOULD be done, but it works better the other way
    akm=1.496e8/1230.262271999999939
    c1 = SkyCoord(akm*xcenter*u.km,akm*ycenter*u.km,6.9e5*u.km,obstime=startdt,frame=frames.Heliocentric)
    lat=c1.transform_to('heliographic_stonyhurst').lat.value
    
    #Differential rotation in deg per day based on solar latitude
    A=14.713
    B=-2.396
    C=-1.787
    xoffsetdegperday=A+B*np.sin(lat/60)**2+C*np.sin(lat/60)**4
    
    rotdays=180/xoffsetdegperday
    xoffset=2460/rotdays/2/3600 #Arcsec per 12 sec

    
    for wv in wavelength:
        #Create directory for each wavelength if not already there
        path=basedir+subdir+str(wv)+'/'
        if not os.path.isdir(path):
            os.makedirs(path)
        
        xcenter = float(xcenterinit)
        
        
        #Skips wavelengths that have already been downloaded to save time
        if os.path.isfile(f"{basedir}{subdir}{wv}/{endstr}AIA_{wv}.png") or os.path.isfile(f"{basedir}{subdir}{wv}/{endstr}HMI_Mag.png"):
            continue
        
        #If the program crashed in the middle of downloading, you can skip the the nearest hour and day that needs to be downloaded
        #This math applies the rotation offset needed to keep your object in frame where it was before the program failed
        startdttemp = startdt + timedelta(days=skipdays,hours=skiphours)
        xcenter = xcenter + (skipdays*24+skiphours)*300*xoffset
        
        skipdays=0
        skiphours=0
        
        #Download png files     
        for t in rrule(SECONDLY, dtstart=startdttemp, until=enddt):
            if t.second%12 == 0:
                xcenter = xcenter + xoffset
                if wv == 'hmi':
                    hv.download_png(t, 0.6, "[SDO,HMI,HMI,magnetogram,1,100]", x0=xcenter, y0=ycenter, width=width, height=height,directory=path,watermark=False,overwrite=True)
                else:
                    hv.download_png(t, 0.6, "[SDO,AIA,AIA,"+str(wv)+",1,100]", x0=xcenter, y0=ycenter, width=width, height=height,directory=path,watermark=False,overwrite=True)
    

def timestamp(
start='2015/01/17 07:00:00'
,end='2015/01/17 08:30:00'
,basedir='/data/reu/wwainwri/Data/solar/'
,wavestamp='171'
,font="/Library/Fonts/arial.ttf"
,fontsize=30):
   
    """Uses Python Image Library (PIL) to annotate the images with a timestamp over the given time range.
    The script timestamps every picture in the folder given by the start day. The output is sent to a folder
    in the same directory as the input wavelength folder, with a '-timestamped' appended to the directory name.
    Expects input files named in a YYYY_mm_dd_HH_MM_SS_AIA_wavelength.png format.
    
    Args:
        start (str): Start time string for the download, expected format YYYY/mm/dd HH:MM:SS
        
        end (str): End time string for the download, expected format YYYY/mm/dd HH:MM:SS
        
        basedir (str): Directory where the data will be downloaded to, subdirectory will be made based on date, expects full /path/to/dir/ format
        
        wavestamp (str): String for the specific wavelength to timestamp, default is the bottom left corner in a 3x3 grid (AIA 171)
        
        font (str): Path to and filename of the font you want to use, default "/Library/Fonts/arial.ttf"
        
        fontsize (int): Font size you want to use, default 30
    """
    
    #Imports
    from PIL import Image,ImageDraw,ImageFont
    import os
    
    #Type Error
    wavestamp=str(wavestamp)
    
    
    #Datetime objects
    
    #Check for bad directory format
    if not basedir[-1] == '/':
        basedir = basedir + '/'  
    
    
    #Sub directory path from date input
    subdir=start.split(' ')[0]+'/'
    
      
    #Make the directory
    path = basedir+subdir+wavestamp+'/'
    if not os.path.isdir(basedir+subdir+wavestamp+'-timestamped/'):
        os.makedirs(basedir+subdir+wavestamp+'-timestamped/')
    
    #Timestamp the png files for the specified wavelength
    for fn in os.listdir(path):
        im=Image.open(path+fn)
        draw = ImageDraw.Draw(im)
        fonts = ImageFont.truetype(str(font), fontsize)
        ts = fn.split('.')[0].split('_')
        #Format the timestamp from filename, expects YYYY_mm_dd_HH_MM_SS_AIA_171.png format
        timestamp = ts[0]+'/'+ts[1]+'/'+ts[2]+'   '+ts[3]+':'+ts[4]+':'+ts[5]
        draw.text((20, im.size[1]-50), timestamp, font=fonts)
        im.save(basedir+subdir+wavestamp+'-timestamped/'+fn)





def mk_movie (
start='2015/01/17 07:00:00'
,end='2015/01/17 08:30:00'
,wavelength=[94,131,171,193,211,304,335,1600,'hmi']
,basedir='/data/reu/wwainwri/Data/solar/'
,outdir='/home/wwainwri/Programs/wainwright/aiamovie/'
,timestamp=True
,wavestamp='171'
,mode='3x3'):
    
    """FFMpeg based movie maker. Relies on the specific format of png file downloaded
    by the get_png() function. For the best results, please use the same parameters
    as used in the get_png() function. All movies are compiled on a per-wavelength 
    basis before being assembled into the specified grid size.
    
    Args:
        start (str): Start time string for the movie to pull data from, expected format YYYY/mm/dd HH:MM:SS
        
        end (str): End time string for the movie to pull data from, expected format YYYY/mm/dd HH:MM:SS
        
        wavelength (list): List object of wavelengths, for example [171,304,'hmi'], wavelengths may be string or int objects
        
        basedir (str): Directory where the data will be pulled from, *USE THE SAME DIRECTORY FROM get_png()*, expects full /path/to/dir/ format
        
        outdir (str): Directory where the movies will be created, subdirectory will be made based on date, expects full /path/to/dir/ format
        
        timestamp (bool): Boolean of whether or not to use the timestamped image set, default True, do not use True if you called get_png() with False
        
        wavestamp (str): String for the specific wavelength to timestamp, default is the bottom left corner in a 3x3 grid (AIA 171)
        
        mode (str): Video mode string, allowed: '3x3', '2x2', '1x1' or 'singular' - anything not matching a pattern will default to single only
    
    
    """
    #Imports
    from datetime import datetime
    import os
    import subprocess
    
    #Datetime objects
    startdt=datetime.strptime(start,'%Y/%m/%d %H:%M:%S')
    enddt=datetime.strptime(end,'%Y/%m/%d %H:%M:%S')
    
    #Formats a string for the start and end time from the input times
    dtname=datetime.strftime(startdt,'%Y_%m_%d_%H%M-')
    dtname = dtname + datetime.strftime(enddt,'%Y_%m_%d_%H%M')
    
    #Check for bad directory format
    if not basedir[-1] == '/':
        basedir = basedir + '/'
    if not outdir[-1] == '/':
        outdir = outdir + '/'    
    
    #Check for wavelength list length
    if mode == '3x3' and not len(wavelength) == 9:
        raise Exception('Wavelength array must be exactly 9 long')
     #Check for wavelength list length
    if mode == '2x2' and not len(wavelength) == 4:
        raise Exception('Wavelength array must be exactly 4 long')
    
    #Extend output directory with date
    outdir = outdir + dtname + '/'
    
    #Sub directory path from date input
    subdir=start.split(' ')[0]+'/'
    
    #Make sure the folder exists
    if not os.path.isdir(outdir+'individual/'):
        os.makedirs(outdir+'individual/')
    
    #Empty String
    readin=''
    
    #Type Error
    wavestamp = str(wavestamp)
    
    for wv in wavelength:    
        if str(wv) == wavestamp and timestamp:
            #Checks the timestamp folder instead of the regular one
            command = f"ffmpeg -pattern_type glob -i '{basedir}{subdir}{wv}-timestamped/*.png' -pix_fmt yuv420p {outdir}individual/aia{wv}movie_timestamp.mp4"
            subprocess.call(['/bin/tcsh','-c',command])
            #Indexes the location and name of the file ffmpeg just made for the purposes of tiling in the next step
            readin = readin + f"-i {outdir}individual/aia{wv}movie_timestamp.mp4 "
        else:
            #The general command for creating the rest of the mp4s
            command = f"ffmpeg -pattern_type glob -i '{basedir}{subdir}{wv}/*.png' -pix_fmt yuv420p {outdir}individual/aia{wv}movie.mp4"
            subprocess.call(['/bin/tcsh','-c',command])
            #Indexes the location and name of the file ffmpeg just made for the purposes of tiling in the next steo
            readin = readin + f"-i {outdir}individual/aia{wv}movie.mp4 "
    
    if mode == '3x3':
        #3x3 mosaic of videos - Calls ffmpeg to tile the movies it just created / already exist
        command = f"ffmpeg {readin} -filter_complex \"[0:v][1:v][2:v]vstack=inputs=3[col0];[3:v][4:v][5:v]vstack=inputs=3[col1];[6:v][7:v][8:v]vstack=inputs=3[col2];[col0][col1][col2]hstack=inputs=3[v]\" -map \"[v]\" {outdir}aiapanel_{mode}_movie.mp4"
        subprocess.call(['/bin/tcsh','-c',command])
        
    if mode == '2x2':
        command = f"ffmpeg {readin} -filter_complex \"[0:v][1:v]vstack=inputs=2[col0];[2:v][3:v]vstack=inputs=2[col1];[col0][col1]hstack=inputs=2[v]\" -map \"[v]\" {outdir}aiapanel_{mode}_movie.mp4"
        subprocess.call(['/bin/tcsh','-c',command])



def test_coords (
time='2015/01/17 07:30:00'
,wavelength='171'
,outdir='~/Desktop/'
,delete=True):
    
    """Test coordinate function for use prior to the get_png() function.
    Use by specifying the start time of your desired observation time window
    and any wavelength you can identify your target in. Before running the script,
    it may be necessary to type '%matplotlib qt' without quotes to get the window
    to pop up. A matplotlib plot should show your coordinates where your cursor is.
    Use these coordinates when calling the get_png() function as your xcenter and ycenter
    coordinates at the same starting time. If you change your starting time you will
    need to adjust the coordiates.
    
    Args:
        time (str): Start time string for test image to be pulled from, expected format YYYY/mm/dd HH:MM:SS
        
        wavelength (str): String or int wavelength that you intend to observe, preferably one where you can see your target
        
        outdir (str): Directory of where to (temporarily) download the file to, probably safe to leave this default
        
        delete (bool): Unless you specifically want to keep the test file, then leave this to True, will remove the test image after viewing
    """
    
    #A note to anyone reading this in the future - If the image the get_png method downloads is not in the same spot
    #as what you selected from this test function, then it is probable that sunpy or helioviewer changed its methods and
    #the math that is required to convert from px to arcsec offset.
    
    #Imports
    from sunpy.net import helioviewer
    from datetime import datetime
    import os
    import matplotlib.pyplot as plt
    
    #Datetime object
    timedt = datetime.strptime(time,'%Y/%m/%d %H:%M:%S')
    
    #Check for bad directory format
    if not outdir[-1] == '/':
        outdir = outdir + '/'    
    
    #Create the directory if it is not there
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    #Confirm wavelength variable type
    wavelength = str(wavelength)
    
    #Creates an instance of the helioviewer client
    hv = helioviewer.HelioviewerClient()
    
    #The parameters required to show a full disk image of the sun
    height=4096
    width=4096
    x0=0
    y0=0
    
    
    #Downloads and opens test image
    filename= hv.download_png(timedt, 0.6, "[SDO,AIA,AIA,"+wavelength+",1,100]", x0=x0, y0=y0, width=width, height=height, directory=outdir, watermark=False,overwrite=True)
    im = plt.imread(filename)
    fig,ax=plt.subplots()
    ax.imshow(im)
    #Deletes the test image unless otherwise specified
    if(delete):
        os.remove(filename)
        
    







