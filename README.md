- [Introduction](#introduction)
- [Terminology](#terminology)
- [Overview of R-RANSAC](#overview-of-r-ransac)
  * [Data Management](#data-management)
  * [Model Initialization](#model-initialization)
  * [Model Management](#model-management)
- [Setup](#setup)
- [Documentation](#documentation)
- [Developer's Guide](#developer-s-guide)
- [References](#references)



## Introduction
R-RANSAC stands for recursive random sample consensus. It is a modular multiple target tracking (MTT) paradigm that autonomously initializes tracks without knowing the number of tracks a priori. 

The purpose of this documentation is to give a high level overview of R-RANSAC, a brief overview of the code structure and documentation, a setup guide, and a developers guide.

## Terminology 
Before giving an overview of R-RANSAC, we will present some common terminology that is frequently encountered. Some of this terminology has been taken from the RADAR community. 

* **Phenomenon**: Something that produces an observable signal. In the case of target tracking, the phenomenon is referred to as a **target**, which is an object that can move physically.

* **Measurement Source**: A sensor equipped with an algorithm that captures information from the environment and produces meaningful measurements used to observe the target.

* **Surveillance Region**: The portion of the environment that is observable by the sensors. There is a local surveillance region for each sensor and a global surveillance region that is a union of all the local surveillance regions. For example, the surveillance region of a camera is everything in it's view, and the global surveillance region would be the area viewed by multiple cameras. Note that there can be many local surveillance regions but only one global surveillance region. 

* **Frame of reference**: Consists of an abstract coordinate system and the set of physical reference points that uniquely fix (locate and orient) the coordinate system and standardize measurements within that frame. We will often refer to a frame of reference as just **frame**.

* **Local Frame**: The frame that coincides with a local surveillance region.

* **Global Frame**: The frame that coincides with the global surveillance region. It is possible that the global frame is the same as a local frame. 

* **Sensor Scan**: When a sensor observes its surveillance region and extracts meaningful data. For example, the sensor scan of a camera occurs when the camera produces a new image of its surveillance region. 

* **False Measurement**: A measurement extracted from a sensor scan that does not correspond to a phenomenon of interest. For example, motion in a camera can generate false measurements due to parallax depending on the algorithm. Another example is just noisy sensors. 

* **True Measurement**: A measurement extracted from a sensor scan that does correspond to a phenomenon of interest. 

* **Model**: This is simply a model of the phenomenon. In regards to target tracking, a model is referred to as a **track**.

* **Model Hypothesis**: This is a hypothetical model of the phenomenon (i.e. a possible model) created by the RANSAC algorithm. A model hypothesis that meets certain criteria becomes a model. In regards to target tracking, a model hypothesis is referred to as a **track hypothesis**. We will often abbreviate the term and mention only hypothesis.

* **Good Model**: A model that is deemed very likely to correctly describe a phenomenon, based on predefined criteria, becomes a good model. In regards to target tracking, a good model is referred to as a **good track**. 

* **Poor Model**: A model that is not a good model. In regards to target tracking, a poor model is referred to as a **poor track**. 

* **Time Window**: An interval of time extending into the past from the current time. 

* **Expired Measurement**: A measurement that was observed in the past outside the time window. 

* **Measurement Source**: An algorithm that takes sensor data and produces a measurement. We will often refer to a measurement source as just a source. 




## Overview of R-RANSAC

The figure below gives a high level depiction of the data dn process flow of R-RANSAC. The left side of the image is not part of R-RANSAC but depicts the flow of information from the sensors to R-RANSAC. The cylinders on the left represent different sensors with a total of "n" sensors. Each sensor has its own local surveillance region and local frame. When a sensor scan occurs, the sensor captures data (measurements) and transforms the data to the global frame. Once measurements are represented properly in the global frame, they are passed into R-RANSAC along with a transformation denoted as "T". The transformation "T" is used when the global frame moves, and contains the data necessary to transform all previous measurements and models in R-RANSAC to the new global frame. This way, all of the measurements and models, old and new, are represented in the current global frame and not in different global frames. 

![RRANSAC Overview: A depiction of the data and process flow](/images/Overview.png )

### Data Management

All of the new measurements and transformation are passed into the data manager. The data manager removes old measurements (measurements whose time stamps are beyond the time window) from the data tree, clusters and consensus sets. The data tree is a data structure based on the R*-tree data structure and stores all measurements that do not pertain to clusters or a consensus sets. There is only one data tree. A cluster is a data structure that contains a set of neighboring measurements. There can be multiple clusters, one for each set of neighboring measurements. The consensus set is a data structure that contains measurements affiliated with a model. Each model has one consensus set. 

Once all of the old measurements are removed, all of the measurements in the data tree, clusters, and consensus sets, except for the newly received measurements, and models are transformed to the new global frame provided that a transformation was given. If no transformation was given, it is assumed that the global frame did not move. 

Once all of the measurements and models are transformed to the current global frame, the new measurements are associated with existing models. Typically, each measurement is checked to see if it falls within a validation region of a model using a data association method. Associated measurements are used to update the model they are associated with, and then are added to the model's consensus set.

The measurements not associated with any model are checked to see if they pertain to a cluster. If they do, they are added to the cluster. The remaining measurements are used to seed a potential new cluster. Taking a new measurements, we find all of the neighboring measurements and add them to the cluster. We then find all of the neighboring measurements to the measurements just added to the cluster and add them to the cluster. We repeat the process recursively until there are no other neighboring measurements to add to the cluster. If the cluster is large enough, then it is kept and the associated measurements are removed from the tree; otherwise, the cluster is discarded. 

All remaining new measurements are added to the data tree. This completes the data management process. 

### Model Initialization 

The model initializer is based on the RANSAC algorithm. For each cluster, measurements are randomly sampled to form a model hypothesis. Each model hypothesis is validated using other measurements from the same cluster. The best model hypothesis that reaches certain criteria, from each cluster, is used to create a new model, and all of the measurements that supported the model hypothesis are added to the model's consensus set and removed from the cluster they came from. 

Any cluster that had measurements removed from it, is checked to see if it is still a valid cluster. If not, it is dissolved and the measurements are put back onto the data tree. 

### Model Management

With the generation of new models or additional information gathered from the new measurements, existing models can begin to coalesce. The coalescing models are merged together. Models that are no longer probable models (usually the models that haven't received new measurements for a while), are pruned from the list of models. The model's probability is then checked to see if it should be escalated to a good model or demoted to a regular model. Any model that is promoted to the status of a good model receives a unique ID number to identify it. 


## Setup

## Documentation

The code is documented using Doxygen. To generate the documentation, navigate to the root directory of the project and run the following command in a terminal
```
doxygen Doxyfile
```

The documentation folder should now be populated with a sub-folder named html. There are many ways to view the documentation. One method uses phython 3. In the same terminal, enter the command

```
cd documentation/html
python3 -m http.server
```
Then open up a web browser and type in 
```
localhost:8000
```
into the web browser.

Another method is to point a HTML browser to the index.html file. To do this, open up a search engine and type the following in the web browser
```
file://localhost/<path_to_rransac>/rransac/documentation/html/index.html
```
The main page is the README.md file seen on github. The tabs Classes and Files can be used to navigate through the code. 


## Developer's Guide

Google style guide.

## References 

Various researchers have contributed to the development of R-RANSAC. Below is a list of the more relevant published papers. We recommend reading them in the order listed.

1. Niedfeldt, P. C., & Beard, R. W. (2014). Multiple target tracking using recursive RANSAC. Proceedings of the American Control Conference, 3393–3398. https://doi.org/10.1109/ACC.2014.6859273
2. Niedfeldt, P. C., Ingersoll, K., & Beard, R. W. (2017). Comparison and Analysis of Recursive-RANSAC for Multiple Target Tracking. IEEE Transactions on Aerospace and Electronic Systems, 53(1), 461–476. https://doi.org/10.1109/TAES.2017.2650818
3. Yang, F., Tang, W., & Lan, H. (2017). A density-based recursive RANSAC algorithm for unmanned aerial vehicle multi-target tracking in dense clutter. IEEE International Conference on Control and Automation, ICCA, (k 1), 23–27. https://doi.org/10.1109/ICCA.2017.8003029
4. Defranco, P. C., Beard, R. W., Warnick, K. F., & Mclain, T. W. (2015). Detecting and Tracking Moving Objects from a Small Unmanned Air Vehicle. All Theses and Dissertations, (March).



