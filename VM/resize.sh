#!/bin/bash
RESOLUTIONS=("1280x720" "1600x900" "1920x1080" "2560x1440" "3840x2160")

function set_standard_mode() {
    # Generate Modeline using cvt and strip the frequency from the label
    # Extract the actual mode configuration (after the label) to keep only resolution details
    for RESOLUTION in "${RESOLUTIONS[@]}"; do 

        MODECONFIG=$(cvt $(echo $RESOLUTION | tr 'x' ' ') | grep Modeline | sed -E 's/Modeline "[^"]+" //')
        # Apply xrandr commands with the clean label
        cmd="xrandr --newmode \"$RESOLUTION\" $MODECONFIG"
        eval "$cmd"
        cmd="xrandr --addmode screen $RESOLUTION"
        eval "$cmd"
    done
}

# Get all unique modes (resolutions) from the first column of xrandr output
modes=$(xrandr | awk '{ print $1 }' | grep '^[0-9]*x[0-9]*$' | sort -u)

# Check if we found any modes
if [ -z "$modes" ]; then
    zenity --error --text="No available screen resolutions found."
    exit 1
fi
# check if there is 1 mode
if [ $(echo "$modes" | wc -l) -eq 1 ]; then
    zenity --info --text="Only one screen resolution available: $modes - will set standard modes"
    set_standard_mode
    modes=$(xrandr | awk '{ print $1 }' | grep '^[0-9]*x[0-9]*$' | sort -u)
fi


# Display the unique list of available modes in zenity
selected_mode=$(echo "$modes" | zenity --list --title="Select Screen Resolution" \
  --text="Choose a screen resolution" \
  --column="Resolution" --width=300 --height=400)

# Check if the user pressed cancel or didn't select a mode
if [ -z "$selected_mode" ]; then
    echo "No resolution selected. Exiting."
    exit 1
fi

# Prompt the user to select an output if multiple displays are connected
outputs=$(xrandr | grep " connected" | awk '{ print $1 }')
output_count=$(echo "$outputs" | wc -l)

if [ "$output_count" -gt 1 ]; then
    # Let the user choose an output if more than one is connected
    selected_output=$(echo "$outputs" | zenity --list --title="Select Display Output" \
      --text="Choose a display output to apply the resolution" \
      --column="Output" --width=300 --height=200)

    # Check if the user pressed cancel or didn't select an output
    if [ -z "$selected_output" ]; then
        echo "No display output selected. Exiting."
        exit 1
    fi
else
    # Use the single output if only one display is connected
    selected_output="$outputs"
fi

# Apply the selected resolution to the chosen output using xrandr
if xrandr --output "$selected_output" --mode "$selected_mode" 2>/dev/null; then
    zenity --info --text="Resolution set to $selected_mode on $selected_output successfully."
else
    zenity --error --text="Failed to set resolution. $selected_mode may not be supported on $selected_output."
    exit 1
fi

# Check if the background image path is stored in ~/.fehbg
if [ -f ~/.fehbg ]; then
    # Source ~/.fehbg to reapply the same background image
    source ~/.fehbg
else
    zenity --warning --text="No existing background image found to reapply."
fi
