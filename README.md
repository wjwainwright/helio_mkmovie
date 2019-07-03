# helio_mkmovie
Uses the sunpy to download cropped pngs from helioviewer and calls ffmpeg to assemble a movie grid of different AIA wavelengths


Check the docstrings and comments for more info

Example use:

import helio_mkmovie as helio

helio.test_coords(time='2017/09/17 00:00:00') 

This time should be the time you intend to use as the start time for the get_png function. Make sure to use %matplotlib qt. Pops up a window where you can zoom in and get the center coordinates for your object, write down the x and y value in the matplotlib coordinates.

helio.get_png (start='2017/07/09 00:00:00', end='2017/07/16 00:00:00', wavelength=[94,131,171,193,211,304,335,1600,'hmi'], basedir='~/data/helio/', xcenter=307, ycenter=675, width=800, height=600, skipdays=0, skiphours=0)

The purpose of changing skipdays and skiphours instead of the start time should a download fail is so that the tracking catches up to where the object is in the new "start time". This method also ensures that the new files end up in the same folder as the the others, which is named for the start date.

helio.timestamp(start='2017/07/09 00:00:00', end='2017/07/16 00:00:00', basedir='~/data/helio/', wavestamp='171', font="/Library/Fonts/arial.ttf", fontsize=30)

This timestamps the wavelength of choice. Ignore the font options unless you care to change them, all variables have set defaults in the script.

helio.mk_movie (start='2017/07/09 00:00:00', end='2017/07/16 00:00:00', wavelength=[94,131,171,193,211,304,335,1600,'hmi'], basedir='~/data/helio/', outdir='~/data/helio/movies/', timestamp=True, wavestamp='171', mode='3x3')
