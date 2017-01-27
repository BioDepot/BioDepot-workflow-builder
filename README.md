# BioDepot-workflow-builder

Based on alpine fluxbox and GUIdock-novnc. PyQt5 and Qt5-svg libraries were compiled from source. 

To use start the container after it has been built 
```docker run --rm -p 6080:6080 -v /var/run/docker.sock:/var/run/docker.sock -ti -v ~/Desktop/:/data biodepot/bwb```

Finally - use a browser and go to localhost:6080 (for windows docker-machine you will need to know the ip of the container)

Please refer to the User Mannual for more information 

Demo - [Link](https://drive.google.com/file/d/0B6xuS_tbRDJ0RzN6NlJ0T1U4VUU/view?usp=sharing)

 
