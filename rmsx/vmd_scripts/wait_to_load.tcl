# -----------------------------------------------------------
#                VMD IDLE LOADER SCRIPT
# -----------------------------------------------------------
#
# This script waits until VMD is fully loaded (after .vmdrc)
# and then runs your main script, passing all arguments to it.

# --- CONFIGURATION ---
# ▼▼▼ EDIT THIS LINE ▼▼▼
# Set the full path to your *main* script
set main_script_path "./rmsx/vmd_scripts/grid_color_scale_centered_xaxis_hotkeys.tcl"
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
set env(VMDMODULATENEWTUBE) user
puts "--- Set env(VMDMODULATENEWTUBE) to 'user' (during init) ---"
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
        puts "Please edit loader.tcl and fix the 'main_script_path' variable."
        return
    }

    puts "--- VMD is idle and delay is over. Calling main script: $main_script_path ---"

    # Run your main script.
    source $main_script_path

    puts "--- Main script has finished ---"
}

# --- MODIFIED COMMAND ---
# Instead of 'after idle', we use our millisecond delay.
# This will wait for VMD to be idle AND for $additional_delay_ms to pass.
puts "--- Loader script armed. Main script will run after $additional_delay_ms ms. ---"
after $additional_delay_ms run_main_script_after_delay


