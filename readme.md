\## ECSE6560 - Spotfi based localization



This project is about implementing the Spotfi (https://web.stanford.edu/~skatti/pubs/sigcomm15-spotfi.pdf) algorithm on a publicly available wildv2 dataset (https://www.kaggle.com/competitions/wildv2/overview). To run these codes, downlad the wildv2 dataset, then paste the dataset as folder 'wildv2' in the same directory as the codes.



\# MATLAB files description

* cdf\_fft.m

Run this code to plot the CDF for a range of user locations using FFT algorithm. The number of user locations and number of wifi APs to use can be selected.

* cdf\_music.m

This code implements the music algorithm and shows a CDF plot, just like in cdf\_fft.m

* compare\_aoa\_fft\_music.m

This code outputs a figure that has location of all wifi APs, the location of the particular user (which can be varied), the ground truth direction of the user from each AP, the estimated direction from each AP using FFT and MUSIC algorithm. The mean error of each method is also printed.

* compare\_cdf.m

this is a standard plot to load the saved cdf from 'cdf\_fft.m' or 'cdf\_music.m' and plots multiple CDFs in a single figure.

* load\_h5\_structure.m

This is a helper function to load the wildv2 dataset which is in h5 format. It can be quite confusing to navigate through the dataset (in MATLAB, the size of the arrays are opposite to what the website provies. It is just how MATLAB handles h5 files)

* plot\_localization\_scenario.m

This is also a helper function to plot the number of wifi APs and the user location.

