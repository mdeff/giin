# GIIN: Graph-based Image Inpainting

## Installation

### Requirements

* [GSPbox](https://lts2research.epfl.ch/gsp/): graph signal processing toolbox
* [UNLocBoX](http://unlocbox.sourceforge.net/): convex optimization toolbox

### Optional speed up

The FLANN library implements fast algorithms to approximate nearest neighbors
search. It can be used to speed up the graph construction. It is optional and
if you do not want to use it you can comment `param.nnparam.use_flann = 1;` in
`lib/giin_patch_graph.m`.  Otherwise compile it this way:

```sh
cd gspbox/3rdparty/sources/flann-1.8.4-src/
mkdir build && cd build
cmake ..
make
```

Make sure that the MATLAB bin folder is in your path to compile the MEX. If it
did work you should have a file called `nearest_neighbors.mexa64` in
`build/src/matlab`.

## Usage

1. Place a file `image.png` in the `data` sub-folder. The masked area should be
   bright green, i.e. RGB [0,255,0]. You may want to generate some synthetic
   images with e.g. `giin_image('vertical');`.

2. Launch the inpainting process from MATLAB:

   ```m
   inpaint('vertical');
   ```

   Or from the shell:

   ```sh
   ./launch.sh inpaint vertical
   ```

   You may then want to adjust the paths to the toolboxes in `launch.sh`.
