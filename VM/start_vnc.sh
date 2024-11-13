#!/bin/bash

# Define display, default resolution, and port
DISPLAY_NUM=1               # Virtual display number
DEFAULT_RESOLUTION="1920x1080x24" # Default screen resolution
VNC_PORT=5091               # Port for VNC server
PID_FILE="/tmp/vnc_pids.txt" # File to store PIDs

# Function to gracefully stop all services
cleanup() {
    echo "Stopping x11vnc, fluxbox, and Xvfb..."

    if [[ -f "$PID_FILE" ]]; then
        # Read PIDs from the file
        X11VNC_PID=$(grep "x11vnc" "$PID_FILE" | awk '{print $2}')
        FLUXBOX_PID=$(grep "fluxbox" "$PID_FILE" | awk '{print $2}')
        XVFB_PID=$(grep "Xvfb" "$PID_FILE" | awk '{print $2}')

        # Terminate each process if it's running
        kill -SIGTERM "$X11VNC_PID" "$FLUXBOX_PID" "$XVFB_PID" 2>/dev/null

        # Remove PID file after cleanup
        rm -f "$PID_FILE"
    fi
    echo "All services stopped."
    exit 0
}

# Trap SIGINT and SIGTERM signals to call cleanup on script termination
trap cleanup SIGINT SIGTERM

# Start Xvfb virtual display with the default resolution and save the PID
echo "Starting virtual display :$DISPLAY_NUM with resolution $DEFAULT_RESOLUTION..."
Xvfb :$DISPLAY_NUM -screen 0 $DEFAULT_RESOLUTION +extension RANDR &
XVFB_PID=$!
echo "Xvfb $XVFB_PID" > "$PID_FILE"

# Wait for Xvfb to start by checking for the X display socket
for i in {1..10}; do
    if xdpyinfo -display :$DISPLAY_NUM > /dev/null 2>&1; then
        echo "Xvfb started successfully on display :$DISPLAY_NUM"
        break
    else
        echo "Waiting for Xvfb to start..."
        sleep 1
    fi
done

# Exit if Xvfb failed to start
if ! xdpyinfo -display :$DISPLAY_NUM > /dev/null 2>&1; then
    echo "Failed to start Xvfb on display :$DISPLAY_NUM. Exiting."
    cleanup
fi

# Start fluxbox and pass it the Xvfb PID explicitly
echo "Starting fluxbox on display :$DISPLAY_NUM with Xvfb PID $XVFB_PID..."
DISPLAY=:$DISPLAY_NUM FLUXBOX_XVFB_PID=$XVFB_PID fluxbox &
FLUXBOX_PID=$!  # Save the PID of fluxbox

# Append fluxbox PID to the file
echo "fluxbox $FLUXBOX_PID" >> "$PID_FILE"

# Wait for fluxbox to start by checking for the process
for i in {1..10}; do
    if ps -p "$FLUXBOX_PID" > /dev/null; then
        echo "Fluxbox started successfully on display :$DISPLAY_NUM"
        break
    else
        echo "Waiting for fluxbox to start..."
        sleep 1
    fi
done

# Exit if fluxbox failed to start
if ! ps -p "$FLUXBOX_PID" > /dev/null; then
    echo "Failed to start fluxbox on display :$DISPLAY_NUM. Exiting."
    cleanup
fi

# Start x11vnc on the virtual display and set it to listen on the specified port
echo "Starting x11vnc server on port $VNC_PORT..."
x11vnc -display :$DISPLAY_NUM -rfbport $VNC_PORT -forever -auth guess -shared &
X11VNC_PID=$!  # Save the PID of x11vnc

# Append x11vnc PID to the file
echo "x11vnc $X11VNC_PID" >> "$PID_FILE"

echo "x11vnc server started on port $VNC_PORT. Use SSH to forward this port and connect via VNC client."

# Monitoring loop to ensure all processes are running
while true; do
    if ! ps -p "$X11VNC_PID" > /dev/null || ! ps -p "$XVFB_PID" > /dev/null || ! ps -p "$FLUXBOX_PID" > /dev/null; then
        echo "One of the essential processes (x11vnc, Xvfb, or fluxbox) has terminated. Exiting..."
        cleanup  # Call the cleanup function to terminate all processes
    fi
    sleep 2  # Check every 2 seconds
done

