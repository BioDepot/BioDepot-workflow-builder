sudo Xvfb :1 -screen 0 "1920x1080"x24 &
sudo x11vnc -display :1 -forever -xkb -noxrecord -shared &
sudo fluxbox -display :1 &
