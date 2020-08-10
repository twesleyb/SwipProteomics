# SWIP Predicted Protein Structure
A structure of the SWIP protein was generated using [Phyre]().
The SWIP protein has no recognizable domains, but is highly helical. 
SWIP interacts directly with Strumpellin, probably through its C-terminus (Jia
et al., Figure S5).

<p align='center'>
    <img src='../figs/github/SWIP_P1019R.png' height = '250' />
    </p>

In pymol, generate a gif of the SWIP protein can be generated using the following commands:
```
# Generate a gif.
mset 1,180 
util.mroll 1, 180, 1
set ray_trace_frames, 1 
set cache_frames, 0

# Save as a video: File > save as video > PNG images
# The image stack was made into a gif with the img2gif shell script.

```

The full use of Pymol without a liscense is restricted.
A better alternative to pymol is chimera. Convert `protein.pse` (pymol specific format) 
to the univeral .pdb format and create a video in chimera.

```
# Convert the video to a stack  of images.
ffmpeg -i infile.avi -f image2 image-%03d.jpg

# Convert the images to a gif.
img2gif JPG/ 1 -1 wt.gif
```
