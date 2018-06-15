# SWELL-Data-ECG-Algorithm
Analysing heartbeat dynamics using point process adaptive filtering on the SWELL Dataset.

## Refer:
1. Analysis of Heartbeat Dynamics by Point Process Adaptive Filtering pdf file
2. Analysis of work-related stress using SWELL dataset pdf file

## Execution:
1. Open main.m
2. Set address to path of file (.S00)
3. Execute main.m

## Description of files:
### myDoReadData()
Reads data from a Portilab Poly5 document (.S00)
Reads all data (portiHRdata) and seperately stores the data relevant to further processing (HeartPre, SkinRaw)

### peakfinder()
Identifying peaks and location of spikes
min_dist should be multiple of 1/2

### init()
Initializing required vectors

### del_parts()
Dividing [0,T] into J intervals

### mean_rate()
Calculating mean heart rate

### f()
Calculating pdf for History dependent Inverse Gaussian Distribution

### integ_f()
Calculating integral of f()

### cif()
Calculating lambda using the conditional intensity function

### df()
Computing (df/dtheta) for all thetas

### d_lambda()
Calculating 1st derivative of lambda with respect to theta 

### d_loglambda()
Calculating 1st derivative of log(lambda) with respect to theta 

### dsquaref()
Calculating 2nd derivative of f() with respect to theta 

### dsquare_lambda()
Calculating 2nd derivative of lambda with respect to theta 

### dsquare_loglambda()
Calculating 2nd derivative of log(lambda) with respect to theta 
