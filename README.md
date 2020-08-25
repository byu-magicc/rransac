# R-RANSAC

- [Introduction](#introduction)
- [Terminology](#terminology)  
- [Overview of R-RANSAC](#overview-of-r-ransac)

## Introduction
R-RANSAC stands for recursive random sample consensus. It is a modular multiple target tracking (MTT) paradigm that autonomously initializes tracks without knowing the number of tracks a priori. 

The purpose of this documentation is to give a high level overview of R-RANSAC, a brief overview of the code structure and documentation, a setup guide, and a developers guide.

## Terminology 
Before giving an overview of R-RANSAC, we will present some common terminology that is frequently encountered. Some of this terminology has been taken from the RADAR community. 

* **Phenomenon**: Something that produces an observable signal. In the case of target tracking, the phenomenon is referred to as a **target**, which is an object that can move physically.

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

* **Time Window**: An interval of time extending into the past from the current time. 




## Overview of R-RANSAC

The figure below gives a high level depiction of the data dn process flow of R-RANSAC. The left side of the image is not part of R-RANSAC but depicts the flow of information from the sensors to R-RANSAC. The cylinders on the left represent different sensors with a total of "n" sensors. Each sensor has its own local surveillance region and local frame. When a sensor scan occurs, the sensor captures data (measurements) and transforms the data to the global frame. Once measurements are represented properly in the global frame, they are passed into R-RANSAC along with a transformation denoted as "T". The transformation "T" is used when the global frame moves, and contains the data necessary to transform all previous measurements and models in R-RANSAC to the new global frame. This way, all of the measurements and models, old and new, are represented in the current global frame and not in different global frames. 

![RRANSAC Overview: A depiction of the data and process flow](/images/Overview.png )



## References 

Various researchers have contributed to the development of R-RANSAC. Below is a list of the more relevant published papers. We recommend reading them in the order listed.

1. Niedfeldt, P. C., & Beard, R. W. (2014). Multiple target tracking using recursive RANSAC. Proceedings of the American Control Conference, 3393–3398. https://doi.org/10.1109/ACC.2014.6859273
2. Niedfeldt, P. C., Ingersoll, K., & Beard, R. W. (2017). Comparison and Analysis of Recursive-RANSAC for Multiple Target Tracking. IEEE Transactions on Aerospace and Electronic Systems, 53(1), 461–476. https://doi.org/10.1109/TAES.2017.2650818
3. Yang, F., Tang, W., & Lan, H. (2017). A density-based recursive RANSAC algorithm for unmanned aerial vehicle multi-target tracking in dense clutter. IEEE International Conference on Control and Automation, ICCA, (k 1), 23–27. https://doi.org/10.1109/ICCA.2017.8003029
4. Defranco, P. C., Beard, R. W., Warnick, K. F., & Mclain, T. W. (2015). Detecting and Tracking Moving Objects from a Small Unmanned Air Vehicle. All Theses and Dissertations, (March).



