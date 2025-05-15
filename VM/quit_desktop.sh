#!/bin/bash

terminate_process() {
    local PID=$1
    local PROCESS_NAME=$2

    if kill "$PID"; then
        echo "$PROCESS_NAME (PID: $PID) terminated."
    else
        echo "$PROCESS_NAME (PID: $PID) did not terminate on first attempt. Retrying..."
    fi

    # Retry in a loop until the process is confirmed terminated
    for attempt in {1..5}; do
        if ps -p "$PID" > /dev/null; then
            sleep 1
            kill "$PID"
        else
            echo "$PROCESS_NAME (PID: $PID) terminated."
            return
        fi
    done

    # Forcefully kill if still running
    if ps -p "$PID" > /dev/null; then
        echo "$PROCESS_NAME (PID: $PID) did not terminate gracefully. Forcing termination..."
        kill -9 "$PID"
    fi
}

# Prompt the user for confirmation to terminate X11 session
if ! zenity --question --text="Are you sure you want to terminate your X11 session?"; then
    zenity --info --text="Operation canceled by the user."
    exit 0
fi

# Check for environment variables for x11vnc, Xvfb, and fluxbox PIDs
if [[ -n "$FLUXBOX_XVFB_PID" ]]; then
    # Terminate x11vnc and Xvfb first, then fluxbox last
    zenity --info --text="Terminating framebuffer"
    terminate_process "$FLUXBOX_XVFB_PID" "Xvfb"
    exit 0
else
    FLUXBOX_XVFB_PID=$(ps -u "$USER" -o pid,cmd | grep -E 'Xvfb ' | grep -v grep | awk '{print $1}')
    if [[ -n "$FLUXBOX_XVFB_PID" ]]; then
        pid_array=($FLUXBOX_XVFB_PID)
        for pid in "${pid_array[@]}"; do
            zenity --info --text="Terminating Xvfb process (PID: $pid)..."
            terminate_process "$pid" "Xvfb"
        done
    	exit 0
    fi
fi

# If  Xvfb is not running, check for an Xorg process
XORG_PID=$(ps -u "$USER" -o pid,cmd | grep -E 'Xorg|X ' | grep -v grep | awk '{print $1}')

# Check if Xorg process was found
if [ -z "$XORG_PID" ]; then
    zenity --info --text="No Xorg process found for user $USER. Exiting."
    exit 1
fi

# Terminate the Xorg process if found
zenity --info --text="Terminating Xorg process (PID: $XORG_PID)..."
terminate_process "$XORG_PID" "Xorg"
zenity --info --text="Xorg process terminated."
stty echo
