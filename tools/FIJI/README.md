# FIJI User Tutorial

FIJI ([https://imagej.net/Fiji](https://imagej.net/Fiji)) is an image processing package—a "batteries-included" distribution of ImageJ, bundling a lot of plugins which facilitate scientific image analysis.

FIJI can be used to post-process raw X-ray CT images of particles to obtain "cleaned" binary (black/white) stack images. Then the voxel coordinates of the stack images can be exported for constructing the bonded-sphere DEM particle geometry templates as input for LIGGGHTS-INL.

A concise user tutorial is provided.


## Import image sequence

File -> Import -> Image Sequence

<img src="figs/fig_fiji_step01.png">

## Crop image sequence (optional)

Select region. Then go to Image -> Crop

<img src="figs/fig_fiji_step02.png">

## Decide thrshold

### Auto local threshold

<img src="figs/fig_fiji_step03.png">

### Manual threshold

Image -> Adjust -> Threshold

<img src="figs/fig_fiji_step04.png">

## Binarize images

Process -> Binary -> Make Binary

<img src="figs/fig_fiji_step05.png">

## Denoise (optional)

### Salt and pepper

Process -> Noise -> Despeckle

<img src="figs/fig_fiji_step06.png">

### Manual denoise

<img src="figs/fig_fiji_step07.png">

## Subtract background

Process -> Subtract Background

<img src="figs/fig_fiji_step08.png">

## Check the 3D view (optional)

Plug-ins -> 3D Viewer

<img src="figs/fig_fiji_step09.png">

## Save image stack and image sequence

File -> Save As -> Image Sequence

<img src="figs/fig_fiji_step10.png">

## Export *.obj file and visualize in ParaView

File -> Save As -> Waterfront OBJ

<img src="figs/fig_fiji_step11.png">



