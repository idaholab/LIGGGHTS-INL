# FIJI User Tutorial

FIJI ([https://imagej.net/Fiji](https://imagej.net/Fiji)) is an image processing package—a "batteries-included" distribution of ImageJ, bundling a lot of plugins which facilitate scientific image analysis.

FIJI can be used to post-process raw X-ray CT images of particles to obtain "cleaned" binary (black/white) stack images. Then the voxel coordinates of the stack images can be exported for constructing the bonded-sphere DEM particle geometry templates as input for LIGGGHTS-INL.

A concise user tutorial is provided.


## Import image sequence

`File -> Import -> Image Sequence`


<img src="figs/fig_fiji_step01.png">

## Crop image sequence (optional)

Select region. Then do
`Image -> Crop`

We recommend saving the cropped image sequence in TIFF and load it again for the following operations. This saves computer memory if the memory limit is a concern.

<img src="figs/fig_fiji_step02.png">

## Decide threshold

### Auto local threshold

Only good with high-contrast images. Low-contrast images like those of biomass should not use any automatic thresholding.

<img src="figs/fig_fiji_step03.png">

### Manual threshold

`Image -> Adjust -> Threshold`

Low-contrast images like those of biomass should require manual thresholding.

<img src="figs/fig_fiji_step04.png">

## Binarize images

`Process -> Binary -> Make Binary`

Equivalent to auto threshholding. Should not apply to low-contrast images like those of biomass.

<img src="figs/fig_fiji_step05.png">

## Denoise (optional)

### Despeckle

`Process -> Noise -> Despeckle`

<img src="figs/fig_fiji_step06.png">

## Subtract background (optional)

`Process -> Subtract Background`

Do not subtract background if already binarized.

<img src="figs/fig_fiji_step08.png">

## Check the 3D view (optional)

`Plug-ins -> 3D Viewer`

<img src="figs/fig_fiji_step09.png">

## Save image stack and image sequence

`File -> Save As -> Image Sequence`

<img src="figs/fig_fiji_step10.png">

## Export XYZ coordinates and visualize in ParaView

`Analyze -> Tools -> Save XY Coordinates (save all slices)`


