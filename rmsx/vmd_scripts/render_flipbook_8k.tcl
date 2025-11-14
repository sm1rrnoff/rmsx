# --- VMD 8K RENDER SCRIPT (v17 - Manual Frame) ---
#
# This is an 8K version of the manual framing script.
# It will be very slow and produce a huge file.
# -----------------------------------------------------------------

puts "--- Sourcing 8K Render Script (v17) ---"
puts "--- Mode: Manual Framing (8K Quality) ---"

# --- 1. Tune Lighting ---
puts "Setting high-quality AO..."
display aoambient 0.8
display aodirect 0.3

# --- 2. Get Current View Composition ---
set user_zoom [display get height]
lassign [display get size] win_x win_y
if {$win_y == 0} { error "Window height is 0." }
set aspect_ratio [expr {double($win_x) / $win_y}]

puts [format "DEBUG: Read user zoom (height): %.2f" $user_zoom]
puts [format "DEBUG: Read window aspect ratio: %.2f" $aspect_ratio]

# --- 3. 8K Render (Preserving Composition) ---
puts "--- Render Stage ---"

# !!! CHANGED: Set a fixed 8K height
set res_y 4320
set res_x [expr {round($res_y * $aspect_ratio)}]
if {$res_x < 10} { set res_x 10 }

# --- 4. Apply Settings in Correct Order ---
display resize $res_x $res_y
puts "DEBUG: Set resolution to ${res_x}x${res_y}"

display height $user_zoom
display update
puts "DEBUG: Re-applied user zoom."

render aasamples TachyonInternal 16
render aosamples TachyonInternal 64
puts "DEBUG: Set 8K quality (AA=16, AO=64)"

# !!! CHANGED: New filename
set output_filename "My_8K_Render_v17_manual.tga"

# 4. Render the final image
puts "--- Starting 8K render ($res_x x $res_y) ---"
puts "!!! WARNING: This will take a very long time. !!!"
puts "Output file: $output_filename"
render TachyonInternal $output_filename
puts "--- 8K Render complete! ---"
puts "Image saved to $output_filename"
