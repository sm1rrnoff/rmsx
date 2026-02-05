# -----------------------------------------------------------
#                VMD IDLE LOADER SCRIPT
# -----------------------------------------------------------
#
# This script waits until VMD is fully loaded (after .vmdrc)
# and then runs your main script, passing all arguments to it.

# --- CONFIGURATION ---
# ▼▼▼ EDIT THIS LINE ▼▼▼
# Set the path to your *main* script.
# IMPORTANT: Resolve relative to *this* loader script so it works from any CWD
# (e.g., when launched from Jupyter notebooks).
# You can override with environment variable RMSX_VMD_MAIN.
set script_path [file normalize [info script]]
set script_dir [file dirname $script_path]
if {[info exists env(RMSX_VMD_MAIN)] && $env(RMSX_VMD_MAIN) ne ""} {
    set main_script_path $env(RMSX_VMD_MAIN)
} else {
    set main_script_path [file join $script_dir "grid_color_scale_centered_xaxis_hotkeys.tcl"]
}
# --- END CONFIGURATION ---


# -----------------------------------------------------------
#                VMD IDLE LOADER SCRIPT (with Delay)
# -----------------------------------------------------------
#
# This script waits until VMD is fully loaded (after .vmdrc)
# *and* a specified time has passed, then runs your main script.


# ▼▼▼ SET YOUR DELAY HERE (in milliseconds) ▼▼▼
# 1000 = 1 second
set additional_delay_ms 1500
# --- END CONFIGURATION ---


# --- FIX ---
# Set the modulation environment variable IMMEDIATELY.
# This must be done *before* VMD becomes idle, or it will be ignored.
# RMSX's main script uses:
#  - 'user'  for color
#  - 'user2' for thickness modulation
set env(VMDMODULATERIBBON) user2
set env(VMDMODULATENEWTUBE) user2
set env(VMDMODULATENEWCARTOON) user2
puts "--- Set env(VMDMODULATE*) to 'user2' (during init) ---"
# --- END FIX ---

# This procedure will run after the delay
proc run_main_script_after_delay {} {
    # Import the global command-line arguments
    global argc argv

    # Import the main script path we set above
    global main_script_path

    # Check if the file exists before trying to run it
    if {![file exists $main_script_path]} {
        puts "ERROR: Main script not found at '$main_script_path'"
        puts "Please edit wait_to_load.tcl and fix the 'main_script_path' variable."
        return
    }

    puts "--- VMD is idle and delay is over. Calling main script: $main_script_path ---"

    # Run your main script.
    source $main_script_path

    puts "--- Main script has finished ---"
}

# Arm the loader:
# 1) wait until VMD becomes idle (startup + .vmdrc complete)
# 2) wait an additional delay for slower systems / lots of files
puts "--- Loader armed. Main script will run after VMD is idle + $additional_delay_ms ms. ---"
after idle [list after $additional_delay_ms run_main_script_after_delay]
