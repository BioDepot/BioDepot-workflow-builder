sudo Xvfb :1 -screen 0 "1920x1080"x24 &
sudo x11vnc -display :1 -rfbport 5900 -listen 0.0.0.0 -noshm -noipv6 -noxdamage &> /dev/null &
sudo fluxbox -display :1 &
