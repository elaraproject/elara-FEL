import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

class hgroup:
    """
    Class for harmonic groups that perhaps overlap.
    """
    
    css = []
    x = np.array([])
    y = np.array([])
    
    def __init__ (self, cs=None):
        self.css = []
        self.x = np.array([])
        self.y = np.array([])
        if cs is not None:
            self.addhgcs(cs)
        return

    def add (self, x, y):
        return self.addhgcs(CubicSpline(x, y))
        
    def addhgcs (self, hgcs):
        """
        add a cs only if it overlaps or is the first one.
        
        return: None
        """
        if isinstance(hgcs, hgroup):
            if len(self.css) != 0 and not self.overlaps(hgcs):
                raise ValueError('this cs does not overlap and cannot be added')
                
            # This will cover correct 'step' plotting
            if self.overlapstop(hgcs):
                self.x = np.append(self.x, self.x[0]-1e-6)
                self.x = np.append(self.x, hgcs.x[-1]+1e-6)
            if self.overlapsbot(hgcs):
                self.x = np.append(self.x, self.x[1]+1e-6)
                self.x = np.append(self.x, hgcs.x[0]-1e-6)

            self.x = np.concatenate([self.x, hgcs.x[:]])
            self.x.sort()
            self.x = np.unique(self.x)
            for cs in hgcs.css:

                self.css.append(cs)
            self.y = [np.nanmax([csa(x, extrapolate=False) for csa in self.css], ) for x in self.x]
                
            return
        if len(self.css) != 0 and not self.overlaps(hgcs):
            raise ValueError('this cs does not overlap and cannot be added')

        self.x = np.concatenate([self.x, hgcs.x[:]])
        self.x.sort()
        self.x = np.unique(self.x)
        self.css.append(hgcs)
        self.y = [np.nanmax([csa(x, extrapolate=False) for csa in self.css], ) for x in self.x]
        return
    
    def overlaps(self, hgcs):
        """
        Check if this cs or hgroup overlaps with the existing ones
        
        return: True if there is any overlap, False otherwise
        """
        
        a = hgcs.x >= self.x[0]
        b = hgcs.x <= self.x[-1]
        
        return sum(a & b) > 0

    def overlapstop(self, hgcs):
        """
        Check if the high end of this cs or hgroup overlaps with the existing ones
        
        return: True if there is overlap on the high end, False otherwise
        """
        
        a = hgcs.x >= self.x[0]
        
        return sum(a) > 0
    
    def overlapsbot(self, hgcs):
        """
        Check if the low end of this cs or hgroup overlaps with the existing ones
        
        return: True if there is overlap on the low end, False otherwise
        """
        
        b = hgcs.x <= self.x[-1]
        
        return sum(b) > 0
    
    def plot (self, **args):
        return plt.plot(self.x, self.y, **args)
    
    def print(self):
        return print(len(self.x), len(self.y))


class henvelope:
    hgs = []
    name = None
    
    def __init__ (self, hgcs=None, name=''):
        self.hgs = []
        self.name = name
        if hgcs is not None:
            self.addhgcs(hgcs)
        return

    def add (self, x, y):
        return self.addhgcs(CubicSpline(x, y))

    def addhgcs (self, hgcs):
        if isinstance(hgcs, CubicSpline):
            olaps = [ind for ind, hg in enumerate(self.hgs) if hg.overlaps(hgcs)]
            if len(olaps) == 0:
                self.hgs.append(hgroup(hgcs))
            else:
                thishg = hgroup(hgcs)
                for iolap in olaps:
                    thishg.addhgcs(self.hgs[iolap])
                self.hgs = [i for j, i in enumerate(self.hgs) if j not in olaps]
                self.hgs.append(thishg)

        elif isinstance(hgcs, hgroup):
            olaps = [ind for ind, hg in enumerate(self.hgs) if hg.overlaps(hgcs)]
            if len(olaps) == 0:
                self.hgs.append(hgcs)
            else:
                for iolap in olaps:
                    hgcs.addhgcs(self.hgs[iolap])
                self.hgs = [i for j, i in enumerate(self.hgs) if j not in olaps]
                self.hgs.append(hgcs)

        return
    
    def plot (self, **args):
        if len(self.hgs) > 0:
            for hg in self.hgs:
                ret = hg.plot(**args)
            return ret
        return

