import pandas as pd, numpy as np


#mode options are
# df:  return the slice dataframes (default)
# len: return the numbers of entries in each slice
class BinIterator:
    def __init__(self,df,xvar,min,max,bins,mode='df'):
        self._df = df
        self._xvar = xvar
        self._min = min
        self._max = max
        self._bins = bins
        self._i = 0
        self._mode = mode
    def __iter__(self):
        self._i = 0
        return self
    def __next__(self):
        if(self._i >= self._bins):
            raise StopIteration
        dx = (self._max-self._min)/self._bins
        #case:  single dataframe
        if type(self._df)==pd.DataFrame:
            ret = self._df.query("%s >= %s and %s < %s"\
                       %(self._xvar,
                         self._min + self._i*dx,
                         self._xvar,
                         self._min +(self._i+1)*dx))
            if self._mode == 'len':
                ret = len(ret)
        #case:  list of dataframes
        else: 
            ret = [df.query("%s >= %s and %s < %s"\
                       %(self._xvar,
                         self._min + self._i*dx,
                         self._xvar,
                         self._min +(self._i+1)*dx)) for df in self._df]
            if self._mode == 'len':
                ret = [len(r) for r in ret]
        self._i+=1
        return self._min+(self._i-1/2)*dx, ret


def query_or_all(df,q):
    if q != "":
        return df.query(q)
    else:
        return df

from scipy.optimize import curve_fit 
def getmeanstd(df, query,nsigma=1.5):

    a = df.eval(query)
    x0, sigma = np.mean(a), np.std(a)
    #print(x0,sigma)
    nbins = 50
    y,x = np.histogram(a, bins=nbins,range=(x0-nsigma*sigma,x0+nsigma*sigma))
    del a
    x = np.add(x[1:],x[:-1])/2
    
    def gaus(x,a,x0,sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
    
    popt,pcov = curve_fit(gaus,x,y,p0=[1,x0,sigma])

    return popt[1],abs(popt[2]), np.sqrt(pcov[1][1]), np.sqrt(pcov[2][2])
    
    
