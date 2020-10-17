# Video compression projects
 The video filesare saved in as images and can be loaded by the programm.   
 The fold `final project` is the full implementation version. and the alogrithm has a difference with the previous one.  
 e.g. the final project doesn't use APBT CALVC and doesn't apply 4x4 block algorithms  
 ## how to use
 **Run optimization.m** 
keep the foreman images in the fold and do not change inter_intra_parameters.mat  
**Three arguments should be given from the user.**   
**intra_mode** :1 to enable the intra-prediction and 0 to disable. default as 0  
**scales**: change the bitrate and the D R curve. default as [0.07,0.2,0.4,0.8,1.0,1.5,2,3,4,4.5]  
**EOB**: set a large interger as the EOB default as 4000    
The programm automatically load 'inter_intra_parameters.mat' the mat file saves the bitrate and PSNR of the baseline in the chapter 5 and chapter 4.  inter_still_image.m can calcuate the data.  
**output arguments**:  
**PSNR_mean**: save the PSNR of the optimization shape: (len(qScale),1)  
**bitPerPixel**: save the bitrate of the optimization shape: (len(qScale),1)
**figures**  
**titlel 'D R plot with restricted 0-4 bits' plot of the optimization and baseline encoder 0-4 bits**
**title' D R plot' plot of the optimization and baseline encoder len(qScale)**
 
