# --- VMD 4K RENDER SCRIPT (v16 - Manual Frame) ---
#
# This script renders a 4K-resolution image based
# on your *current* manually-set view.
# It uses high-quality settings for publication.
# -----------------------------------------------------------------

puts "--- Sourcing 4K Render Script (v16) ---"
puts "--- Mode: Manual Framing (4K Quality) ---"

# --- 1. Tune Lighting ---
puts "Setting high-quality AO..."
display aoambient 0.8
display aodirect 0.3

# --- 2. Get Current View Composition ---
# Read the user's current zoom level
set user_zoom [display get height]

# Read the user's current window size to get the aspect ratio
lassign [display get size] win_x win_y
if {$win_y == 0} { error "Window height is 0." }
set aspect_ratio [expr {double($win_x) / $win_y}]

puts [format "DEBUG: Read user zoom (height): %.2f" $user_zoom]
puts [format "DEBUG: Read window aspect ratio: %.2f" $aspect_ratio]

# --- 3. 4K Render (Preserving Composition) ---
puts "--- Render Stage ---"

# Set a fixed 4K+ height for the high-res render
set res_y 2160
# Calculate the width based on the user's aspect ratio
set res_x [expr {round($res_y * $aspect_ratio)}]
if {$res_x < 10} { set res_x 10 }

# --- 4. Apply Settings in Correct Order ---
# 1. Resize the render buffer FIRST
display resize $res_x $res_y
puts "DEBUG: Set resolution to ${res_x}x${res_y}"

# 2. Re-apply the user's zoom level SECOND
display height $user_zoom
display update
puts "DEBUG: Re-applied user zoom."

# 3. Set render quality
render aasamples TachyonInternal 16
render aosamples TachyonInternal 64
puts "DEBUG: Set 4K quality (AA=16, AO=64)"

set output_filename "My_4K_Render_v16_manual.tga"

# 4. Render the final image
puts "--- Starting 4K render ($res_x x $res_y) ---"
puts "This may take several minutes."
puts "Output file: $output_filename"
render TachyonInternal $output_filename
puts "--- 4K Render complete! ---"
puts "Image saved to $output_filename"
