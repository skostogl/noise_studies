import cl2pd
from cl2pd import dotdict
from cl2pd import importData
from cl2pd import plotFunctions
from cl2pd import variablesDF
import pandas as pd, seaborn as sns, matplotlib.pyplot as plt
from utils import *
from dateutil import tz
import os
import datetime
import calendar
import h5py as h5

ct = dotdict.dotdict
pd = importData.pd
np = importData.np
cals = importData.cals

plt.rcdefaults()
params = {'legend.fontsize': 16, 'axes.labelsize': 16, 
   'axes.titlesize': 16, 
   'xtick.labelsize': 14, 
   'ytick.labelsize': 14, 
   'image.cmap': 'jet', 
   'lines.linewidth': 2, 
   'lines.markersize': 5}
plt.rcParams.update(params)

def get_fill_info(fill_number):

    """
    Places information about the fill into a dictionary
    Input: fill_number
    Output: dictionary with start-end time of each interval
    """

    info = importData.LHCFillsByNumber(fill_number)
    info.set_index('mode', inplace=True)
    info.index = info.index.where(~info.index.duplicated(), info.index + '2')
    return info.to_dict(orient='index')


def get_HS_tunes(fill_number, mode=None):

    """
    Returns tunes from HS BBQ
    """

    tunes = [
     'LHC.BQBBQ.CONTINUOUS_HS.B%:TUNE_H', 'LHC.BQBBQ.CONTINUOUS_HS.B%:TUNE_V']
    if mode is None:
        tune_data = importData.LHCCals2pd(tunes, fill_number)
    else:
        tune_data = importData.LHCCals2pd(tunes, fill_number, mode)
    return tune_data


def plot_tune_evolution(tune_data, info, vertical_lines=True):

    """
    Plot of the tune evolution furing the cycle
    Input: tune data from "get_HS_tunes", info dictionary of the fill
    """

    vars = [
     'B1:TUNE_H', 'B1:TUNE_V', 'B2:TUNE_H', 'B2:TUNE_V']
    for var in vars:
        data = tune_data['LHC.BQBBQ.CONTINUOUS_HS.%s' % var].dropna()
        fig1, ax1 = plt.subplots(figsize=(12, 6))
        ax1.plot_date(data.index.to_pydatetime(), np.array(data), ms=0.5, c='k')
        jet = plt.get_cmap('jet')
        modes = ['INJPHYS', 'PRERAMP', 'RAMP', 'FLATTOP', 'SQUEEZE', 'ADJUST', 'STABLE', 'ADJUST2', 'BEAMDUMP']
        colors = iter(jet(np.linspace(0, 1, len(modes) + 1)))
        if vertical_lines:
            for mode in modes:
                try:
                    plt.axvline(info[mode]['startTime'], c=next(colors), label=mode)
                except:
                    continue

            ax1.legend(loc=4)
        plt.ylabel(var)

def plot_vlines(info):
    
    """
    Plots vertical lines at start of each period
    """

    jet = plt.get_cmap('jet')
    modes = ['INJPHYS', 'PRERAMP', 'RAMP', 'FLATTOP', 'SQUEEZE', 'ADJUST', 'STABLE', 'ADJUST2', 'BEAMDUMP']
    colors = iter(jet(np.linspace(0, 1, len(modes) + 1)))
    for mode in modes:
        try:
            plt.axvline(info[mode]['startTime'], c=next(colors), label=mode)
        except:
            continue

    plt.legend(loc=4)


def flattenoverlap(v,timestamps, frf,test=100,start=0):

  """
  Remove overlap in BBQ data, from pytimber https://github.com/rdemaria/pytimber
  """

  out=[v[0]]
  out2=[timestamps[0] for i in v[0]]
  out3=[frf[0] for i in v[0]]
  stat=[]
  print("Flatten: ...")
  for j in range(1,len(v)):
    v1=v[j-1]
    v2=v[j]
    newi=0
    for i in range(start,len(v2)-test):
      s=sum(v1[-test:]-v2[i:i+test])
      if s==0:
        newi=i+test
        break
    if newi==0:
      print("Warning: no overlap for chunk %d,%d"%((j-1,j)))
    out.append(v2[newi:])
    out2.append([timestamps[j] for k in v2[newi:]])
    out3.append([frf[j] for k in v2[newi:]])
    stat.append(newi)
  print("average overlap %.2f samples"%np.average(stat))
  return np.hstack(out), np.hstack(out2), np.hstack(out3)

def get_data(modes, time, rename_duplicates=False, remove_overlap=False, n=8000):

    """
    TbT data for the modes and time specified
    If the same mode is specified more than once, it will be renamed
    If remove_overlap is True, it will remove the overlap and combine the data with a sliding window of 8000 turns
    Input: list of modes & dictionary time with start and end time for each key
    Output: df-> beam ('B1, B2') -> plane ('H', 'V') - > 'tbt' (Tbt data & frev interpolated) 
    """

    df     = dotdict.dotdict()
    beams  = ['B1', 'B2']
    planes = ['H', 'V']
    for beam in beams:
        df[beam] = dotdict.dotdict()
        for plane in planes:
            df[beam][plane] = dotdict.dotdict()
            df[beam][plane]['tbt'] = dotdict.dotdict()
            var = ['LHC.BQBBQ.CONTINUOUS_HS.%s:ACQ_DATA_%s' % (beam, plane), 'ALB.SR4.%s:FGC_FREQ' % beam]
            if rename_duplicates:
              counter_dup = 0
            for mode in modes:
                if mode in df[beam][plane]['tbt'] and rename_duplicates:
                  print "Renaming key..."
                  counter_dup +=1
                  new_mode = '%s%s' %(mode, str(counter_dup))
                  df[beam][plane]['tbt'][new_mode] = dotdict.dotdict()
                else:
                  df[beam][plane]['tbt'][mode] = dotdict.dotdict()
                  new_mode = mode
                if time[new_mode][0] == 'all':
                    raw_data = importData.LHCCals2pd(var, time[mode][1], beamModeList=mode)
                else:
                    t1 = time[new_mode][0]
                    t2 = time[new_mode][1]
                    raw_data = importData.cals2pd(var, t1, t2)
                raw_data['status'] = new_mode
                raw_data['ALB.SR4.%s:FGC_FREQ' % beam] = raw_data['ALB.SR4.%s:FGC_FREQ' % beam].interpolate(limit_direction='both')
                if not remove_overlap:
                  df[beam][plane]['tbt'][new_mode] = raw_data
                else:
                  raw_data = raw_data.dropna(subset=[var[0]]) 
                  m = []
                  for i in raw_data[var[0]]:
                    m.append(i)
                  m = np.array(m)
                  test2 = tuple([np.array(raw_data.index), m, np.array(raw_data[var[1]])]) 
                  test={var[0]:test2}
                  flatten={}
                  for name,(timestamps,values, values2) in test.items():
                      flatten[name], timestamps2, frf2=flattenoverlap(values, timestamps, values2)
                  step=1
                  #n = 8000
                  turns = np.arange(0, len(flatten[var[0]]))
                  chunk_t = [turns[x:x+n] for x in xrange(0, len(turns)-n, step)]
                  chunk_var = [flatten[var[0]][x:x+n] for x in xrange(0, len(flatten[var[0]])-n, step)]
                  chunk_time = [timestamps2[x:x+n] for x in xrange(0, len(timestamps2)-n, step)]
                  chunk_frf = [frf2[x:x+n] for x in xrange(0, len(frf2)-n, step)]
                  raw_data2 = pd.DataFrame({ var[0]:chunk_var, 'turns':chunk_t,'timestamps': chunk_time , var[1]: chunk_frf, 'status': new_mode } )
                  df[beam][plane]['tbt'][new_mode] = raw_data2
                   
    return df


def plot_50(df, qxmin=0.275, qxmax=0.31, qymin=0.307, qymax=0.317):
  
  """
  Stem plots for the 50Hz lines - computed with the interpolated frev - and the theoretical limits of the tune between injection and stable.
  Input: main df, (optional) theoretical tunes from injection to stable
  Ouput: df -> beam -> plane -> 'crossing' 
  """

  h = 35640
  for beam in df.keys():
    fig1, ax1 = plt.subplots(figsize=(12, 6))
    fig2, ax2 = plt.subplots(figsize=(12, 6))
    for plane in df[beam].keys():
      df[beam][plane]['crossing'] = pd.DataFrame()
      c   = ['b', 'r']
      amp = [-1, 1]
      counter = 0
      for status in df[beam][plane]['tbt'].keys():
        data = (df[beam][plane]['tbt'][status]).iloc[0]
        frf  = data['ALB.SR4.%s:FGC_FREQ' % beam]
        harm = 50 / (frf / h)
        c1   = c[counter]
        indexes = np.array([ harm * k for k in range(1, 113) ])
        if plane == 'H':
          df2 = pd.DataFrame({'beam': beam, 'status': status, 'horizontal': indexes[np.where((indexes >= qxmin) & (indexes <= qxmax))], 'horizontal #': np.where((indexes >= qxmin) & (indexes <= qxmax))[0]})
          ax1.stem(indexes[np.where((indexes <= qxmin) | (indexes >= qxmax))], [ amp[counter] for i in indexes[np.where((indexes <= qxmin) | (indexes >= qxmax))] ], c1, markerfmt='%so' % c1, label=status)
          ax1.stem(indexes[np.where((indexes >= qxmin) & (indexes <= qxmax))], [ amp[counter] for i in indexes[np.where((indexes >= qxmin) & (indexes <= qxmax))] ], c='g', markerfmt='go')
        else:
          df2 = pd.DataFrame({'beam': beam, 'status': status, 'vertical': indexes[np.where((indexes >= qymin) & (indexes <= qymax))], 'vertical #': np.where((indexes >= qymin) & (indexes <= qymax))[0]})
          ax2.stem(indexes[np.where((indexes <= qymin) | (indexes >= qymax))], [ amp[counter] for i in indexes[np.where((indexes <= qymin) | (indexes >= qymax))] ], c1, markerfmt='%so' % c1, label=status)
          ax2.stem(indexes[np.where((indexes >= qymin) & (indexes <= qymax))], [ amp[counter] for i in indexes[np.where((indexes >= qymin) & (indexes <= qymax))] ], c='g', markerfmt='go')
        df[beam][plane]['crossing'] = pd.concat([df[beam][plane]['crossing'],df2])
        counter += 1

    ax1.vlines(qxmax, ymin=-1.2, ymax=1.2, color='k', linestyle='--', label='$\\rm Q_H$ from injection to stable', linewidth=3.5)
    ax1.vlines(qxmin, ymin=-1.2, ymax=1.2, color='k', linestyle='--', linewidth=3.5)
    ax1.hlines(1.2, xmin=qxmin, xmax=qxmax, color='k', linestyle='--', linewidth=3.5)
    ax1.hlines(-1.2, xmin=qxmin, xmax=qxmax, color='k', linestyle='--', linewidth=3.5)
    ax2.vlines(qymax, ymin=-1.2, ymax=1.2, color='k', linestyle='--', label='$\\rm Q_V$ from injection to stable', linewidth=3.5)
    ax2.vlines(qymin, ymin=-1.2, ymax=1.2, color='k', linestyle='--', linewidth=3.5)
    ax2.hlines(1.2, xmin=qymin, xmax=qymax, color='k', linestyle='--', linewidth=3.5)
    ax2.hlines(-1.2, xmin=qymin, xmax=qymax, color='k', linestyle='--', linewidth=3.5)
    ax1.set_xlim([0.0, 0.5])
    ax2.set_xlim([0.0, 0.5])
    ax1.set_title('%s H' % beam)
    ax2.set_title('%s V' % beam)
    ax1.set_xlabel('Tune')
    ax2.set_xlabel('Tune')
    ax1.set_ylabel('50 Hz lines')
    ax2.set_ylabel('50 Hz lines')
    ax1.legend()
    ax2.legend()
    fig1.tight_layout()
    fig2.tight_layout()



def get_fft(df, fr=50, df_fft = None):

  """
  Computes FFT and organizes by harmonics of 50Hz if fr=50 (h1, h2..) or by bins if fr is a list (f1, f2 ..)
  If df_fft is provided new frequency branches will be appended to the existing df_fft
  If the index turns is detected(slidng window) the first turn of the sliding window will be saved in a new column of the dataframe
  """
  
  h = 35640
  if df_fft is None:
    df_fft = dotdict.dotdict()
  if type(fr) is not list:
    lim = 112 ## int(frev/50/2)
  else:
    lim = len(np.arange(fr[0], fr[-1]))
  lists = [[] for j in range(lim)]
  flag_t = False
  for beam in df.keys():
    for plane in df[beam].keys():
      for status in df[beam][plane]['tbt'].keys():
        group = df[beam][plane]['tbt'][status].dropna()
        print beam, plane, status
        for i in range(len(group)):
            data = group.iloc[i]
            if 'turns' in data.index:
              date = data.timestamps[0]
              turns = data.turns[0]
            else:
              date = data.name
            fourier= np.fft.fft(data['LHC.BQBBQ.CONTINUOUS_HS.%s:ACQ_DATA_%s' %(beam, plane)])
            freq = np.fft.fftfreq(len(data['LHC.BQBBQ.CONTINUOUS_HS.%s:ACQ_DATA_%s' %(beam, plane)]))
            if type(fr) is not list:
                if 'turns' in data.index:
                  frf = data['ALB.SR4.%s:FGC_FREQ' %beam][0]
                else:
                  frf = data['ALB.SR4.%s:FGC_FREQ' %beam]
                harm = fr / (frf / h)
                indexes = [ int(k * harm * len(fourier)) for k in range(1, lim+1) ]
            else:
                indexes = np.arange(fr[0], fr[-1])
            for j in range(len(indexes)):
              if 'turns' in data.index:
                flag_t=True
                lists[j].append([status,j, indexes[j],freq[indexes[j]], fourier[indexes[j]], date, turns, beam, plane])
              else:
                lists[j].append([status,j, indexes[j],freq[indexes[j]], fourier[indexes[j]], date, beam, plane])
  
  if type(fr) is not list:
    key = 'h'
  else:
    key='f'
  for j in range(lim):
    if flag_t:
      df_fft['%s%s' %(key,str(j))] = pd.DataFrame(data = lists[j], columns = ['status', 'h', 'bin','f','fourier', 'timestamps', 'turns', 'beam', 'plane'])
    else:
      df_fft['%s%s' %(key,str(j))] = pd.DataFrame(data = lists[j], columns = ['status', 'h', 'bin','f','fourier', 'timestamps','beam', 'plane'])
      df_fft['%s%s' %(key,str(j))].set_index(['timestamps'],inplace=True)
  return df_fft


def clean_fft(df_fft, search='f'):

  """
  Cleans the fft data for h%: harmonics of 50Hz, f%: FFT bins, 'all' will clean the whole fft tree
  """
  if search == 'all':
    clean_l  = [f for f in df_fft.keys() ]
  else:
    clean_l  = [f for f in df_fft.keys() if f.startswith(search)]
  for i in clean_l:
    df_fft.pop(i, None)
  

def get_abs_phase(df, search='h'):

  """
  Computes the absolute values of the FFT and the angle for h: harmonics of 50Hz, f: FFT bins
  """ 

  for beam in df.keys():
    print beam
    for plane in df[beam].keys():
      print plane
      for status in df[beam][plane]['fft'].keys():
        for harm in [k for k in df[beam][plane]['fft'][status].keys() if k.startswith(search)]:
          try:
            df[beam][plane]['fft'][status][harm].set_index(['timestamps'],inplace=True)
          except:
            continue
          df[beam][plane]['fft'][status][harm]['abs'] = df[beam][plane]['fft'][status][harm]['fourier'].abs()
          df[beam][plane]['fft'][status][harm]['angle'] = np.unwrap(np.angle(df[beam][plane]['fft'][status][harm]['fourier']))
          df[beam][plane]['fft'][status][harm]['un_angle'] = (np.angle(df[beam][plane]['fft'][status][harm]['fourier']))
 
def get_redundant_pairs(df):

  """
  Get diagonal and lower triangular pairs of correlation matrix
  """

  pairs_to_drop = set()
  cols = df.columns
  for i in range(0, df.shape[1]):
    for j in range(0, i + 1):
      pairs_to_drop.add((cols[i], cols[j]))
  return pairs_to_drop



def get_top_abs_correlations(df, n1=0, n2=5, threshold=None):

  """
  Returns strongest n2-n1 absolute correlations if threshold is None, otherwise it will return correlations above threshold 
  """

  au_corr = df.corr().abs().unstack()
  labels_to_drop = get_redundant_pairs(df)
  au_corr = au_corr.drop(labels=labels_to_drop).sort_values(ascending=False)
  if threshold is not None:
    print "Threshold set to ", threshold
    return au_corr[au_corr>threshold]
  print n2-n1, " strongest absolute correlations"
  return au_corr[n1:n2]

def calculate_pvalues(df):

  """
  Returns of p-values of the pearson correlation
  """
  from scipy.stats import pearsonr
  df = df.dropna()._get_numeric_data()
  dfcols = pd.DataFrame(columns=df.columns)
  pvalues = dfcols.transpose().join(dfcols, how='outer')
  for r in df.columns:
    for c in df.columns:
      pvalues[r][c] = round(pearsonr(df[r], df[c])[1], 4)
  return pvalues


def heatmaps(df, status = 'all',beam = 'all', plane='all', mode='amplitude', threshold=None, pval_threshold=None, search = 'h', flag_set_lim=True, harms='all', ax = None):
  
  """
  Heatmaps with correlations
  Input: -df
         -status eg. 'all', ['RAMP'], ['FLATTOP', 'RAMP']
         -beam eg. ['B1', 'B2'], 'all'
         -plane eg. ['H'], 'all'
         -mode eg. 'amplitude', 'angle'
         -threshold: sets a threshold to the correlation values
         -pval_threshold: sets a threshold to the p-values
         -search: 'h' means harmonics of 50Hz, 'f' means bins
         -flag_set_lim: set limit to the colorbar of the heatmap
         -harms: 'all', ['h1', 'h2']
         -ax: pass specific subplot
  Output: -correlation matrix
          -strongest absolute correlations, 300 first, or if threshold is defined the correlations > threshold
          -bins of harmonics
          -matrix with p-values
  """

  corr_tot = dotdict.dotdict()
  pval_tot = dotdict.dotdict()
  strongest = dotdict.dotdict()

  if ax is not None:
    counter_beam = 0
    counter_plane = 0

  if beam == 'all':
    beam = df[df.keys()[0]]['beam'].unique()
  for beams in beam:
    if plane == 'all':
      plane = df[df.keys()[0]]['plane'].unique()
    for planes in plane:
      if planes == plane[0] and ax is not None:
        counter_plane=0 
      if status == 'all':
        stat = df[df.keys()[0]]['status'].unique()
      for stat in status:
        print beams, planes, stat
        dfnew = pd.DataFrame()
        if harms=='all':
          max_lim = len([k for k in df.keys() if k.startswith(search)])
          min_lim = 0
        elif type(harms) is list:
          min_lim = harms[0]
          max_lim = harms[1]
        bins_tot = []
        for h in range(min_lim,max_lim):
          group = df['%s%s' %(search,h)]
          group = group[ (group['beam'] == beams) & (group['plane'] == planes) & (group['status'] == stat)]
          bins_tot.append(group['bin'].iloc[0])
          if mode == 'amplitude':
            dfnew = pd.concat([dfnew, group['fourier'].abs()], axis=1, ignore_index=True)
          elif mode=='angle':
            signal1 = np.unwrap(np.angle(group['fourier']))
            turns = group['turns']
            x0 = np.angle(group['fourier'])[0]
            f = group['f'].iloc[0]
            signal2 = x0 + 2.0*np.pi*f*turns
            dfnew = pd.concat([dfnew, signal1-signal2], axis=1, ignore_index=True)
        corr = dfnew.corr()
        pval = calculate_pvalues(dfnew)
        if ax is None:
          fig1, ax1 = plt.subplots(figsize=(9, 6))
        else:
          try:
            plt.sca(ax[counter_plane, counter_beam])
          except:
            plt.sca(ax[counter_beam])
        plt.title('%s , %s%s , %s' % (mode, beams, planes , stat), fontsize=16)
        plt.xlabel('h')
        plt.ylabel('h')
        
        for i in range(corr.shape[0]):
          corr.iloc[i, i] = 0.0
        if threshold is not None:
          if pval_threshold is None:
            corr[corr.abs()<threshold] =0.
          else:
            mask = (pval<pval_threshold) & (corr.abs()>threshold) 
            corr = corr*mask
        elif pval_threshold is not None:
          mask = (pval<pval_threshold) 
          corr = corr*mask
        if flag_set_lim:
          axs = sns.heatmap(corr, xticklabels=15, yticklabels=15, vmin=-1, vmax=1, cbar_kws={'label': 'Correlation coefficients'}, cmap='seismic')
        else:
          axs = sns.heatmap(corr, xticklabels=15, yticklabels=15, cbar_kws={'label': 'Correlation coefficients'}, cmap='seismic', ax = ax[counter_beam, counter_plane])
        plt.tight_layout()
        #plt.show()
        a = get_top_abs_correlations(dfnew, n1=0, n2=300, threshold = threshold)
        print 'Top Absolute Correlations'
        print a
        corr_tot['%s_%s_%s_%s' %(beams, planes, stat, mode)] = corr
        pval_tot['%s_%s_%s_%s' %(beams, planes, stat, mode)] = pval
        strongest['%s_%s_%s_%s' %(beams, planes, stat, mode)] = a
      if ax is not None:
        counter_plane+=1
    if ax is not None:
      counter_beam+=1
  return corr_tot, strongest, bins_tot, pval_tot


def plot_harm(hs, df, status = 'all',beam = 'all', plane='all', mode='amplitude', normalise=True, search='h', fit_slope=False, remove_slope=False):

  """
  Plots the time evolution of the amplitude/angle of the harmonics specified
  Input: -hs eg. [23, 42]
         -df
         -normalise: normalise to max
         -search: 'h' for harmonics of 50Hz, 'f' for bins  
  """


  if beam == 'all':
    beam = df[df.keys()[0]]['beam'].unique()
  if plane == 'all':
    plane = df[df.keys()[0]]['plane'].unique()
  if status == 'all':
    status = df[df.keys()[0]]['status'].unique()
  for beams in beam:
    for planes in plane:
      fig1, ax1 = plt.subplots(figsize=(12,6))
      plt.title('%s, %s' %(beams,planes))
      for stat in status:
        print beams, planes, stat
        jet= plt.get_cmap('jet')
        colors = iter(jet(np.linspace(0,1,len(hs) + 1)))  
        for h in hs:
          c1 = next(colors)
          label = 'h%s' %h
          group = df['%s%s' %(search,h)]
          group = group[ (group['beam'] == beams) & (group['plane'] == planes) & (group['status'] == stat)]
          if mode == 'amplitude':
            ax1.grid()
            if normalise:
           
              a = group['fourier'].abs().max()
            else: 
              a = 1.0
            if stat==status[0] :
              (group['fourier'].abs()/a).plot(ax=ax1, c=c1, label=label, legend=True)
              #ax1.scatter(group['fourier'].index, group['fourier'].abs()/a, c=c1, label=label)
              ax1.legend()  
              ax1.grid()
            else:
              (group['fourier'].abs()/a).plot(ax=ax1, c=c1, legend=False)
              #ax1.scatter(group['fourier'].index, group['fourier'].abs()/a, c=c1)
            ax1.set_ylabel('Amplitude')  
          elif mode == 'angle':
            ax1.grid()
            ax1.set_ylabel('Angle')  
            if not remove_slope:
              if stat==status[0] :
                ax1.plot(group['turns'], np.unwrap(np.angle(group['fourier'])), c=c1, label = label)
                ax1.legend()
              else:
                ax1.plot(group['turns'], np.unwrap(np.angle(group['fourier'])), c=c1)
              if fit_slope:
                turns = np.array(group['turns'])
                x0 = np.angle(group['fourier'])[0]
                f = group['f'].iloc[0]
                signal = x0 + 2.0*np.pi*f*turns
                ax1.plot(turns, signal, linestyle = '--', linewidth=2, c='r')  
            else:
              turns = np.array(group['turns'])
              x0 = np.angle(group['fourier'])[0]
              f = group['f'].iloc[0]
              signal1 = np.unwrap(np.angle(group['fourier']))
              signal2 = x0 + 2.0*np.pi*f*turns
              ax1.plot(turns, signal1-signal2, linestyle = '--', linewidth=2, c=c1, label = label)
            ax1.legend() 

def ob_tbt(filename, fill_number):
  df = importData.LHCFillsByNumber(fill_number)
  df = df[df['mode']!='FILL']
  df = df.reset_index(drop=True) 

  fi = h5.File(filename, 'r')
  date = (filename.split('_')[-2])
  time = (filename.split('_')[-1]).split('.')[0]
  beam = (filename.split('_')[-4])[0:2]
  plane = (filename.split('_')[-4])[2:3]
  if plane == 'H':
    plane = 'horizontal'
  else:
    plane = 'vertical'
  date2 = ((datetime.datetime.strptime(date+time, "%Y%m%d%Hh%Mm%Ss")))
  date_current = ((datetime.datetime.strptime(date+time, "%Y%m%d%Hh%Mm%Ss")).replace(tzinfo=tz.gettz('CET'))).astimezone(tz.gettz('UTC'))
  print "File at: ",  date2, ' (CET) ', date_current, ' (UTC) ' , ' Status: ' , df[(df['startTime']<=date_current)  & (df['endTime']>=date_current)]['mode'].values
  alldat = fi[beam][plane]
  return alldat
