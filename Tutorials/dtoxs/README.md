
# Tutorial: Dtoxs workflow

## Preparation 

## Run [setup] workflow

The setup workflow widget downloads the tutorial files (for example the human reference data) needed to run Dtoxs workflow. This avoids carrying large data files inside the container. The widget is a custom container widget that calls a bash script inside the Bwb container. 

 The setup widget is stored in a mini-workflow file. Start up the Bwb application by following the instructions in the [manual](https://github.com/BioDepot/BioDepot-workflow-builder)

  
 After Bwb launches, open `/root/tutorials/dtoxs/Setup.ows` workflow, then run this workflow which will prepare the data for dtoxs tutorial. You should see a single widget.
 
 ![setup1](media/dtoxs_tutorial_setup.png)
 
 After double-clicking on the icon you should see the custom condtainer form below.


![setup2](media/dtoxs_tutorial_setup2.png)

The fields are already filled in. The Docker image for the container for the widget is the container that is running Bwb. The Run command points to the bash script inside this container that will by run by the widget inside the container. So when the widget is started it will run the bash script inside the same container that is running Bwb.

Click on the start button to run the script and download all the files needed to run the tutorial.

Once the data has been downloaded by the widget you should see a 'Finished!' flag under the widget. We are ready to run the DetoxS pipeline.





## Run DetoxS pipeline

Open `/data/tutorials/dtoxs/dtoxs.ows` workflow, the workflow will be started automatically (Be sure the setup workflow was run before the Dtoxs pipeline, otherwise there will be no results since the input data is not ready). Once the alignment was finished, start the analysis step by double click the [Dtoxs Analysis] icon and [Run].

![run](media/dtoxs_tutorial_runpipeline.png)
