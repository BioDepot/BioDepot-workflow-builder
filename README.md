# GUIdock-lite-orange-qt5

Based on alpine fluxbox and GUIdock-novnc. PyQt5 and Qt5-svg libraries were compiled from source. 

To use start the container after it has been built 
    sudo docker  run --rm -it -p 6080:6080 -v ${PWD}:/local GUIdock-lite-orange /bin/bash
    
 Then start the server by typing in the container terminal
    ../startup.sh &
    DISPLAY=:1 xterm -e orange-canvas

Finally - use a browser and go to localhost:6080 (for windows docker-machine you will need to know the ip of the container)


 
