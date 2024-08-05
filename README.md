# DynAMoS - Dynamic Adaptive Movement Segmentation

```DynAMoS``` is an algorithm developed for segmenting upper-limb movements during functional task execution (i.e., drinking task) from a wrist-worn single Inertial Measurement Unit (IMU). ```DynAMoS``` is based on adaptive thresholding with statistics-based post-processing aimed to reduce erroneous segmentation by applying statistical consideration to the movement duration histograms. The algorithm was developed on the analysis of the angular velocity norm derived from the gyroscope recording from the IMU.

## DynAMoS description
First, movement onset and offset are identified by applying an adaptive threshold to the angular velocity norm $Ω$. In particular, the adaptive threshold was defined as $Th_{DynAMoS}= 0.11 Ω_{max}$, where $Ω_{max}$ represents the maximum value of the angular velocity norm recorded during the task. Then, the duration of each segmented movement was calculated as the difference between the offset and onset time instants. Finally, the statistics-based post-processing is iteratively applied to the movement durations according to the following steps:
1. Definition of the movement duration histogram and computation of the median ($M$) value;
2. Starting from the median value $M$, the lower threshold $αM$ (with $0<α<1$) and the upper threshold $βM$ (with $1<β<2$) are obtained. The two thresholds will be used to identify data distribution outliers (i.e., movements characterized by “atypical” durations). The algorithm analyzes movements with "atypical" durations. If a movement has a duration lower than the threshold $αM$, the algorithm tries to merge it with the preceding or following movement. Meanwhile, if a movement has a duration longer than the upper threshold $βM$, the algorithm tries to split it into two movements. In the case of a duration lower than the threshold $αM$, the algorithm attempts merging the movement under analysis with the preceding or the following, separately. The first attempt is performed with the merging candidate closest in time. If it fails, the other movement is considered. Merging fails if the new movement duration is lower than $αM$ or higher than the $βM$ thresholds. If none of the attempts satisfies the thresholds, the movements are not merged. In case of merging success, the extremities of the “parent” movements are used as starting and ending points. In the case of a movement with a duration longer than the upper threshold $βM$, the algorithm tries to split it into two movements. To this purpose, the algorithms use local minima points in the signal as possible splitting points. Starting from the minimum point with the lowest value, it splits the movement in two, with the new ending and new starting at the identified minimum. If these movements have a length longer than αM and shorter than βM thresholds, the split is accepted and the two new movements are created. Otherwise, the algorithm moves to another minimum point, if it exists;
3. After each splitting or merging event, $M$, $αM$, and $βM$ values are updated considering the new movements;
4. The algorithm runs iteratively until all the outliers in movement duration are processed.

## DynAMoS Functionalities
The current version of ```DynAMoS``` consists of a library of functions that allow the execution of the algorithm from the preprocessing (filtering and angular velocity norm) to the extraction of the onset and offset and the application of the iterative algorithm.
The functions are the following:
1. ```[fnorm] = preprocess(signals,fs)```: As input the signals from the 3-axis of the gyroscope and the sampling frequency fs.
2. ```[int] = intervals(act)```: This function gets as input (act) a binary vector that represents the movement (1) and non-movement (0) states identified from the thresholding and returns a $nx2$ matrix. The first column represents the onset sample and the second is the offset sample.
3. ```[ints] = ints_correction(ints,signal)```: This function executes the iterative process of the algorithm as described above. Inputs are the onset and offset obtained from the ```intervals``` function and the norm of the angular velocity (signal).
4. ```[onoff] = dynamos(fnorm)```: This function, with input the norm of the angular velocity, performs ```DynAMoS``` on the norm of the angular velocity. It returns a $nx2$ matrix where the first column represents the onset sample and the second is the offset sample.

## Hands-on Example
In the data.mat file sample data are available to test the algorithm. Data sampling frequency $fs = 100Hz$.
