import os
import glob
import h5py
from gwpy.timeseries import TimeSeries
import numpy as np
import scipy
from scipy.stats import chi2
import matplotlib.pyplot as plt

def retrieve_file_paths(base_path, station, day, month, gps_start=None, gps_end=None):
    # Get all the paths for a given station in a base_path
    # To do: convert gps segment to UTC and format the path appropriately.
    path_glob = os.path.join(base_path, \"/GNOMEDrive/gnome/serverdata/%s/2017/%s/%s_2017%s%s_*.hdf5"\% (station, month, day, station, month, day))
    return glob.glob(path_glob)

def retrieve_data_timeseries(hfile, setname):
    dset = hfile[setname]
    sample_rate = dset.attrs["SamplingRate(Hz)"]
    gps_epoch = contruct_utc_from_metadata(dset.attrs["Date"], dset.attrs["t0"])
    data = retrieve_channel_data(hfile, setname)
    ts_data = TimeSeries(data, sample_rate=sample_rate, epoch=gps_epoch)
    return ts_data

def retrieve_channel_data(hfile, setname):
    return hfile[setname][:]

def construct_utc_from_metadata(datestr, t0str):
    instr = "%d-%d-%02dT" % tuple(map(int, datestr.split('/')))
    instr += t0str
    from astropy import time
    t = time.Time(instr, format='isot',scale='utc')
    return t.gps

def joinedTimeSeries(data_list, setname="MagneticFields"):
    hfile = h5py.File(data_list[0],"r")
    full_data = retrieve_data_timeseries(hfile, "MagneticFields").copy(order = 'C')
        for seg in data_list[1:]:
            hfile = h5py.File(seg, "r")
            nextMinute = retrieve_data_timeseries(hfile, "MagneticFields").copy(order='C')
            full_data.append(nextMinute, gap = "ignore")
            hfile.close()
    return full_data

# This function calculates the excess power matrix for each 2*deltaT-second interval, with noise calculated as the sum o fthe spectrogram matrix elements in the preceeding and proceeding 1*deltaT seconds

def getLocalEP(timeseries, fsamp = 512, deltaT = 1, deltaF = 1, showplot = False, logplot = True, figtitle = ''):
     fNyq = int(famp/2) # This is the nyquist Frequency
     tLim = np.size(timeseries)/(fsamp*60) # duration of signal in minutes
     dt = famp*deltaT;
     dk = int(deltaF*deltaT); # width of frequency block in data point num
     numf = int(fNyq/deltaF); # number of frequency blocks
     numt = int(np.size(timeseries)/dt)-4; # numt is the numer of time bloacks -4 because the first and last 2*dt s interval can't be counted
     f,t,Sxx = scipy.signal.spectrogram(timeseries, fs=fsamp, nperseg = fsamp*deltaT, scaling = 'spectrum', mode = 'psd', noverlap = 10);
     Sxx = Sxx[1:];
     Sxx = np.transpose(Sxx);
     TF = np.zeros((numpt, numf))
    for i in range(0,numt):
        for k in range(0, numf): #summing in each frequency band
            ep = 2.0*np.sum((np.abs(Sxx[i+1][k*dk:(k+1)*dk]+Sxx[i+2][k*dk:(k+1)*dk]))/(Sxx[i][k*dk:(k+1)*dk]+Sxx[i+3][k*dk:(k+1)*dk]));
                Tf[i][k] = ep;
    Tfrep = np.transpose(TF);
    if showplot == True:
        f,(a0, a1) = plt.subplots(2,1,gridspec_kw = {'height_ratios':[3,1]}, figsize = [12,6], sharex =False)
        tme = np.linspace(0,tLim, np.size(timeseries))
        a1.plot(tme, timeseries)
        if logplot == True:
            a0.imshow(20*np.log(Tfrep),extent = [0,tLim, 0, fNyq], aspect = 'auto', origin ='lower');
        else:
            a0.imshow(TFrep, extent = [0, tLim, 0, fNyq], aspect = 'auto', origin = 'lower');
        a1.set_xlabel('time (min)');
        a1.set_ylabel('amplitude');
        a1.set_xlim([0,tLim])
        a0.set_ylabel('freq (Hz)')
        a0.set_title('Excess Power: Delta T = '+str(deltaT)+' min. Delta f = '+str(deltaF)+ 'Hz.')
        plt.suptitle(figtitle)
        plt.show()
    return TFrep;
# Compute and plot a histogram of the elements of the Excess Power Mat.
# At the same time, an analytic form of the $\chi^2$ distribution with dof degrees is plotted.
# In the absence of a signal, the elements of Mat should follow a chi-square distribution.
# dof = 2*dt*df, where dt and df are the sizes of the tiles in seconds and Hz respectively.
# f_k(\epsilon) = chi2.pdf(\epsilon, k, loc = 0) is the chi-squared distribution for k degrees of freedom for random variable \epsilon,the excess power.
# nbins can take a nomber or an arracy indicating the location of the bins:
# eg. nbins = 20 evenly divides the x axis to 20 equally spaced segments.
# nbins = np.linspace(0, 500, 200) generates 200 bins equally spaced between 0 and 500.
def plotExPowHist2(Mat, dof, log = False, chiSquare = False, nbins = 20):
    Mvals = np.reshape(Mat, -1);
    x = np.linspace(0, 15*dof, 500);
    y = chi2.pdf(x, dof, loc = 0)
    plt.figure(figsize = [12, 4])
    if log == False:
        plt.hist(Mvals, density = True, bins = nbins, rwidth = 0.9);
        if chiSquare == True:
            plt.plot(x, y, 'r-', label = 'Chi-Square dis (dof, = '+str(dof)+')');
            plt.legend()
        plt.xlabel('Excess Power, $\\epsilon$')
        plt.ylabel('$f_k(\\epsilon)$')
    elif log == True:
        plt.hist(Mvals, density = True, bins = nbins, rwidth = 0.9, log = True)
        if chiSquare == True:
            plt.semilogy(x, y, 'r-', label = 'Chi-Square dis (dof, k = '+str(dof)+')')
            plt.legend()
        plt.xlabel('Excess Power, $\\epsilon$')
        plt.ylabel('$f_k(\\epsilon)$')
    plt.show()
