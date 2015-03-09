#!/bin/bash
export MATLABPATH=$MATLABPATH:'/media/share/toolboxes/unlocbox/'
export MATLABPATH=$MATLABPATH:'/media/share/toolboxes/gspbox/'
CMD="$@; exit"
echo "MATLAB command: $CMD"
matlab -nodisplay -nojvm -r "$CMD" -logfile results/$2.log