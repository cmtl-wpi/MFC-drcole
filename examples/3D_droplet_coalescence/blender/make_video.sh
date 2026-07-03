#!/bin/bash
set -u
cd /tmp/claude-1000/-home-daveygravy-repos-MFC/d2645757-1127-4573-b24b-1c9c631af71f/scratchpad
B="$PWD/blender/blender"
NF=101                      # timesteps 0..101
rm -f ply_anim/*.ply frames_anim/*.png
mkdir -p ply_anim frames_anim

echo "=== [1/4] extracting $((NF+1)) isosurfaces (parallel) ==="
seq 0 $NF | xargs -P 16 -I{} bash -c \
  'n=$(printf "%03d" "$1"); python3 extract_iso.py "$1" ply_anim/f$n.ply alpha1 2.0' _ {} > extract.log 2>&1
echo "extracted: $(ls ply_anim/*.ply 2>/dev/null | wc -l) plys"

echo "=== [2/4] fixed camera scale from max radius ==="
RMAX=$(grep -oE 'rmax=[0-9.]+' extract.log | cut -d= -f2 | sort -gr | head -1)
FIXED=$(python3 -c "print(f'{$RMAX*1.08:.3f}')")
echo "max rmax=$RMAX -> MFC_FIXED_R=$FIXED"
export MFC_FIXED_R=$FIXED

echo "=== [3/4] rendering (parallel, 8x30 threads) ==="
seq 0 $NF | xargs -P 8 -I{} bash -c \
  'n=$(printf "%03d" "$1"); "'"$B"'" -b -t 30 -P render_blender.py -- ply_anim/f$n.ply frames_anim/f$n.png 128 3q opaque >/dev/null 2>&1; echo -n .' _ {}
echo ""
echo "rendered: $(ls frames_anim/*.png 2>/dev/null | wc -l) frames"

echo "=== [4/4] encoding mp4 ==="
ffmpeg -y -framerate 20 -i frames_anim/f%03d.png \
  -c:v libx264 -pix_fmt yuv420p -crf 18 -movflags +faststart \
  coalescence.mp4 >ffmpeg.log 2>&1
ls -la coalescence.mp4 && echo "DONE"
