#!/bin/bash

# Find the Xorg process for the current user
XORG_PID=$(ps -u "$USER" -o pid,cmd | grep -E 'Xorg|X ' | grep -v grep | awk '{print $1}')

# Check if Xorg process was found
if [ -z "$XORG_PID" ]; then
    zenity --info --text="No Xorg process found for user $USER. Exiting."
    exit 1
fi

# Prompt the user for confirmation using zenity
if zenity --question --text="Are you sure you want to terminate your X11 session?"; then
    # User clicked "Yes"
    zenity --info --text="Terminating Xorg process (PID: $XORG_PID)..."
    kill "$XORG_PID"

    # Verify if the process is still running, and if so, forcefully kill it
    if ps -p "$XORG_PID" > /dev/null; then
        zenity --warning --text="Process did not terminate gracefully. Forcing termination..."
        kill -9 "$XORG_PID"
    fi
    zenity --info --text="Xorg process terminated."
    stty echo
else
    # User clicked "No"
    zenity --info --text="Operation canceled by the user."
fi
