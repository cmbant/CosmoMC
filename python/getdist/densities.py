import numpy as np
from scipy.interpolate import splrep, splev, bisplrep, bisplev

class DensitiesError(Exception):
    pass

defaultContours = [0.68, 0.95]

def getContourLevels(inbins, contours=defaultContours, missing_norm=0, half_edge=True):
    """
     Get countour levels encolosing "countours" of the probability.
     inbins is the density. If half_edge, then edge bins are only half integrated over in each direction.
     missing_norm accounts of any points not included in inbins (e.g. points in far tails that are not in inbins)
     Works for any dimension bins array
    """
    contour_levels = np.zeros(len(contours))
    if half_edge:
        abins = inbins.copy()
        lastindices = [-1] + [slice(None, None, None) for d in abins.shape[1:]]
        firstindices = [0] + [slice(None, None, None) for d in abins.shape[1:]]
        for _ in abins.shape:
            abins[tuple(lastindices)] /= 2
            abins[tuple(firstindices)] /= 2
            lastindices = np.roll(lastindices, 1)
            firstindices = np.roll(firstindices, 1)
    else:
        abins = inbins
    norm = np.sum(abins)
    targets = (1 - np.array(contours)) * norm - missing_norm
    if True:
        bins = abins.reshape(-1)
        indexes = inbins.reshape(-1).argsort()
        sortgrid = bins[indexes]
        cumsum = np.cumsum(sortgrid)
        ixs = np.searchsorted(cumsum, targets)
        for i, ix in enumerate(ixs):
            if ix == 0:
                raise DensitiesError("Contour level outside plotted ranges")
            h = cumsum[ix] - cumsum[ix - 1]
            d = (cumsum[ix] - targets[i]) / h
            contour_levels[i] = sortgrid[ix] * (1 - d) + d * sortgrid[ix - 1]
    else:
        # old method, twice or so slower
        try_t = np.max(inbins)
        lastcontour = 0
        for i, (contour, target) in enumerate(zip(contours, targets)):
            if contour < lastcontour: raise DensitiesError('contour levels must be decreasing')
            lastcontour = contour
            try_b = 0
            lasttry = -1
            while True:
                try_sum = np.sum(abins[inbins < (try_b + try_t) / 2])
                if try_sum > target:
                    try_t = (try_b + try_t) / 2
                else:
                    try_b = (try_b + try_t) / 2
                if try_sum == lasttry: break
                lasttry = try_sum
            contour_levels[i] = (try_b + try_t) / 2
    return contour_levels


class GridDensity(object):
    def normalize(self, by='integral', in_place=False):
        if by == 'integral':
            norm = self.integral()
        elif by == 'max':
            norm = np.max(self.P)
            if norm == 0:
                raise DensitiesError('no samples in bin')
        else:
            raise DensitiesError("Density: unknown normalization")
        if in_place: self.P /= norm
        else: self.setP(self.P / norm)
        self.spl = None
        self.missing_norm /= norm
        return self

    def setP(self, P=None):
        if P is not None:
            for size, ax in zip(P.shape, self.axes):
                if size <> ax.size:
                    raise DensitiesError("Array size mismatch in Density arrays: P %s, axis %s" % (size, ax.size))
            self.P = P
        else:
            self.P = np.zeros([ax.size for ax in self.axes])
        self.spl = None

    def bounds(self):
        # give bounds in order x, y, z..
        b = [(ax[0], ax[-1]) for ax in self.axes]
        b.reverse()
        return b

    def getContourLevels(self, contours=defaultContours):
        return getContourLevels(self.P, contours, missing_norm=self.missing_norm)

class Density1D(GridDensity):

    def __init__(self, x, P=None, missing_norm=0):
        self.n = x.size
        self.axes = [x]
        self.x = x
        self.spacing = x[1] - x[0]
        self.missing_norm = missing_norm
        self.setP(P)

    def bounds(self):
        return self.x[0], self.x[-1]

    def _initSpline(self):
        self.spl = splrep(self.x, self.P, s=0)

    def Prob(self, x, derivative=0):
        if self.spl is None: self._initSpline()
        if isinstance(x, np.ndarray):
            return splev(x, self.spl, derivative)
        else:
            return splev([x], self.spl, derivative)

    def integral(self):
        return  ((self.P[0] + self.P[-1]) / 2 + np.sum(self.P[1:-1])) * self.spacing

    def initLimitGrids(self, factor=None):
        class InterpGrid(object): pass

        if self.spl is None: self._initSpline()
        g = InterpGrid()
        if factor is None:
            g.factor = max(2, 20000 // self.n)
        else:
            g.factor = factor
        g.bign = (self.n - 1) * g.factor + 1
        vecx = self.x[0] + np.arange(g.bign) * self.spacing / g.factor
        g.grid = splev(vecx, self.spl)

        norm = np.sum(g.grid)
        g.norm = norm - (0.5 * self.P[-1]) - (0.5 * self.P[0])

        g.sortgrid = np.sort(g.grid)
        g.cumsum = np.cumsum(g.sortgrid)
        return g

    def getLimits(self, p, interpGrid=None, accuracy_factor=None):
        g = interpGrid or self.initLimitGrids(accuracy_factor)
        parr = np.atleast_1d(p)
        targets = (1 - parr) * g.norm - self.missing_norm * g.factor
        ixs = np.searchsorted(g.cumsum, targets)
        results = []
        for ix, target in zip(ixs, targets):
            trial = g.sortgrid[ix]
            if ix > 0:
                d = g.cumsum[ix] - g.cumsum[ix - 1]
                frac = (g.cumsum[ix] - target) / d
                trial = (1 - frac) * trial + frac * g.sortgrid[ix + 1]

            finespace = self.spacing / g.factor
            lim_bot = (g.grid[0] >= trial)
            if lim_bot:
                mn = self.x[0]
            else:
                i = np.argmax(g.grid > trial)
                d = (g.grid[i] - trial) / (g.grid[i] - g.grid[i - 1])
                mn = self.x[0] + (i - d) * finespace

            lim_top = (g.grid[-1] >= trial)
            if lim_top:
                mx = self.x[-1]
            else:
                i = g.bign - np.argmax(g.grid[::-1] > trial) - 1
                d = (g.grid[i] - trial) / (g.grid[i] - g.grid[i + 1])
                mx = self.x[0] + (i + d) * finespace
            if parr is not p: return mn, mx, lim_bot, lim_top
            results.append((mn, mx, lim_bot, lim_top))
        return results


class Density2D(GridDensity):

    def __init__(self, x, y, P=None, missing_norm=0):
        self.x = x
        self.y = y
        self.axes = [y, x]
        self.missing_norm = missing_norm
        self.setP(P)

    def integral(self):
        norm = np.sum(self.P[1:-1, 1:-1]) + (self.P[0, 0] + self.P[0, -1] + self.P[-1, 0] + self.P[-1, -1]) / 4.0 \
            + (np.sum(self.P[1:-1, 0]) + np.sum(self.P[0, 1:-1]) + np.sum(self.P[1:-1, -1]) + np.sum(self.P[-1, 1:-1])) / 2.0
        norm *= (self.x[1] - self.x[0]) * (self.y[1] - self.y[0])
        return norm

    def _initSpline(self):
        self.spl = bisplrep(self.x, self.y, self.P, s=0)

    def Prob(self, x, y, dx=0, dy=0):
        if self.spl is None: self._initSpline()
        return bisplev(x, y, self.spl, dx, dy)
