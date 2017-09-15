![BwB](https://github.com/BioDepot/BioDepot-workflow-builder/raw/master/Media/logo.png)

[![](https://images.microbadger.com/badges/image/biodepot/bwb.svg)](https://microbadger.com/images/biodepot/bwb "Get your own image badge on microbadger.com")  [![](https://images.microbadger.com/badges/version/biodepot/bwb.svg)](https://microbadger.com/images/biodepot/bwb "Get your own version badge on microbadger.com")


# BioDepot Workflow Builder (BwB)

BioDepot is a self-contained tool with graphical user interface for bioinformatic workflows. The package is based on Orange 3 by Biolab and NoVnc. Widgets are mainly written in Python (Qt5, Docker-Py, PyQt5).

[![Gif](https://j.gifs.com/58jYKR.gif)](https://youtu.be/VY1peA4ITog)

## Requirements  

- Docker 1.13.0 or above
- Internet Browser  

## Running BwB
Currently, BwB uses docker sock binding, to run BwB:

1. Install Docker   
2. On Docker-enabled machines run:  
``` 
docker run --rm -p 6080:6080 -v ${PWD}:/data -v /var/run/docker.sock:/var/run/docker.sock biodepot/bwb
```
3. BwB can be accessed from browser: http://localhost:6080 (for windows machines, you'll need to know the address of your docker-machine)  

## Sample Workflows
Some sample workflows are available: 
1. [Adaptation of Michael Love's et al. _Aligment at Gene Level_ Alignment](Sample_Workflows/Airway). Alignment (STAR), Count matrix computation, differential gene analyses (DESeq).
2. [DToxS RNA-Seq: Burroughs-Wheeler Aligner for alignment and EdgeR for differential gene expression](Sample_Workflow/DToxS_RNASeq)
3. [More sample workflows](Sample_Workflow)

## Developing Widget
To add a widget to BioDepot:

1. Download BioDepot source codes:  `https://github.com/BioDepot/BioDepot-workflow-builder/tree/master/biodepot`    
2. Locate to folder orangebiodepot (eg: `~/Desktop/biodepot/orangebiodepot`)   
3. Write the widget (or copy any of existing widgets and modify as necessary)   
4. On BwB, run xterm and locate to the shared widget location (e.g: /data/biodepot)   
5. Rebuild BioDepot and rerun Orange:  

```
    pip3 uninstall BioDepot
    pip3 install -e .
    orange-canvas
```

![Screenshot](https://github.com/BioDepot/BioDepot-workflow-builder/raw/master/Media/Screenshot.png)

## More resources:
- Demo: [Link](https://youtu.be/VY1peA4ITog)
- Manual: [Manual](https://github.com/kristiyanto/BioDepot-workflow-builder/blob/master/simplified_manual.pdf)

