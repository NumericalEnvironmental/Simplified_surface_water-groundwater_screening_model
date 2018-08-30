###########################################################################
#
# RiverRun.py - a script to help quantify how groundwater banking responds
# to an engineered discharge event into a dry wash as artificial recharge
#
###########################################################################


from numpy import *
from pandas import *
from scipy.integrate import quad
from scipy.special import *
import matplotlib.pyplot as plt


### classes ###


class WellSet:              # well locations (for monitoring groundwater mound time series)

    def __init__(self, aquifer):
        self.name = []
        self.x = []
        self.y = []
        self.tSeries = []       # containers for time-series head perturbations
        self.hSeries = []
        lineInput = []
        inputFile = open('wells.txt','r')
        for line in inputFile: lineInput.append(line.split())
        inputFile.close()
        for i in range(1, len(lineInput)):
            self.name.append(lineInput[i][0])
            self.x.append(float(lineInput[i][1]))
            self.y.append(float(lineInput[i][2]))
            t0 = float(lineInput[i][3])
            tf = float(lineInput[i][4])
            numInt = int(lineInput[i][5])
            self.tSeries.append(linspace(t0, tf, numInt, endpoint=True))
            self.hSeries.append(zeros(numInt, float) + aquifer.h0)
        self.x = array(self.x)
        self.y = array(self.y)
        print('Read well attributes.')  
        self.name

        
    def Rotate(self, reachItem):
        # offset and clockwise rotation of well locations through angle theta relative to reach orientation
        xp = (self.x-reachItem.x)*cos(pi-reachItem.orient) - (self.y-reachItem.y)*sin(pi-reachItem.orient)
        yp = (self.x-reachItem.x)*sin(pi-reachItem.orient) + (self.y-reachItem.y)*cos(pi-reachItem.orient)
        return xp, yp


    def TimeSeries(self, aquifer, reachSubset):
        # create time series response to recharge through each reach
        print('Calculating resultant groundwater mounding.')
        for rch in reachSubset:
            v = (rch.qGWavg/aquifer.Sy)
            travel = aquifer.bz/v
            xp, yp = self.Rotate(rch)                           # rotated positions of wells with respect to reach
            for i in range(len(xp)):                            # for each well location
                tStart = self.tSeries[i] - travel - rch.tGW0
                tEnd = self.tSeries[i] - travel - rch.tGWf
                for j, tOn in enumerate(tStart):
                    if tOn > 0.:                    
                        self.hSeries[i][j] += aquifer.dh(xp[i], yp[i], tOn, rch)
                for j, tOff in enumerate(tEnd):
                    if tOff > 0.:
                        self.hSeries[i][j] -= aquifer.dh(xp[i], yp[i], tOff, rch)        


    def WriteTimeSeries(self):
        # write the well mounding time histories to output files
        for i, well in enumerate(self.name):
            fileName = well + '.csv'
            SeriesFile = open(fileName,'w')
            SeriesFile.writelines(['time', ',', 'GW_elev'])
            SeriesFile.writelines(['\n'])
            for j in range(len(self.tSeries[i])):
                SeriesFile.writelines([str(self.tSeries[i][j]), ',', str(self.hSeries[i][j])])
                SeriesFile.writelines(['\n'])
            SeriesFile.close()    
        

    def PlotTimeSeries(self):
        # well mounding time histories plots
        colors = ['red', 'blue', 'green', 'magenta', 'cyan', 'olive']
        for i, well in enumerate(self.name):
            plt.plot(self.tSeries[i], self.hSeries[i], color = colors[i], label = well)
        plt.xlabel('Time')
        plt.ylabel('GW Elevation')
        plt.legend(loc=4)
        plt.show()        


class Aquifer:              # aquifer properties for groundwater mounding calculations

    def __init__(self):
        lineInput = []        
        inputFile = open('aquifer.txt','r')
        for line in inputFile: lineInput.append(line.split())
        inputFile.close()
        self.K = float(lineInput[0][1])        # hydraulic conductivity
        self.Sy = float(lineInput[1][1])       # specific yield
        b = float(lineInput[2][1])             # aquifer thickness
        self.v = self.K*b/self.Sy              # aquifer diffusivity
        self.h0 = float(lineInput[3][1])       # background saturated thickness
        self.bz = float(lineInput[4][1])       # unsaturated zone thickness
        print('Read aquifer properties.')


    def dh(self, x, y, t, basin):       # Hantush (1967) model
        z = basin.qGWavg/(2.*self.K) * self.v*t \
            * (self.S((basin.l + x)/sqrt(4*self.v*t), (basin.a + y)/sqrt(4*self.v*t))
            + self.S((basin.l + x)/sqrt(4*self.v*t), (basin.a - y)/sqrt(4*self.v*t))
            + self.S((basin.l - x)/sqrt(4*self.v*t), (basin.a + y)/sqrt(4*self.v*t))
            + self.S((basin.l - x)/sqrt(4*self.v*t), (basin.a - y)/sqrt(4*self.v*t)))
        return sqrt(z + self.h0**2.) - self.h0
    
    
    def S(self, alpha, beta):
        return quad(self.Integrand, 0., 1., args=(alpha, beta))[0]
    
    
    def Integrand(self, tau, alpha, beta):
        return erf(alpha/sqrt(tau)) * erf(beta/sqrt(tau))

    
class Tally:        # data structure for tracking cumulative fluxes through time

    def __init__(self):
        self.t = [0.]
        self.inflow = [0.]
        self.outflow = [0.]
        self.GWflow = [0.]
        self.streamVol = [0.]
        
    def Balance(self, reach, connection, source, srcList, t, dt):
        # water volume distributions
        inflow = 0.
        GWflow = 0.
        streamVol = 0.
        for i, rch in enumerate(reach):
            streamVol += rch.h * rch.floorArea
            GWflow += rch.qGW * rch.floorArea * dt
        for i in range(len(srcList)):         # external source at time t, if applicable
            inflow += source[i].Flux(t-dt) * dt
        self.t.append(t)
        self.inflow.append(inflow + self.inflow[-1])
        self.GWflow.append(GWflow + self.GWflow[-1])
        self.outflow.append(connection[0].Q*dt + self.outflow[-1])
        self.streamVol.append(streamVol)        

    def WriteBalances(self):
        # write balance summaries to tally file
        tallyFile = open('tally.csv','w')
        tallyFile.writelines(['time', ',', 'inflow', ',', 'outflow', ',', 'GW', ',', 'stream'])
        tallyFile.writelines(['\n'])
        for i in range(len(self.t)):
            tallyFile.writelines([str(self.t[i]), ',', str(self.inflow[i]), ',', str(self.outflow[i]), ',', str(self.GWflow[i]), ',', str(self.streamVol[i])])
            tallyFile.writelines(['\n'])
        tallyFile.close()       
        

class Params:       # general simulated scenario attributes
    
    def __init__(self):
        lineInput = []        
        inputFile = open('params.txt','r')
        for line in inputFile: lineInput.append(line.split())
        inputFile.close()
        self.hMin = float(lineInput[0][1])           # minimum water level in reach required water for flux out of reach
        self.dt0 = float(lineInput[1][1])            # starting time step
        self.dhMax = float(lineInput[2][1])          # max water level change in any reach per time step
        self.dtMax = float(lineInput[3][1])          # max time step, regardless of head change
        self.dtf = float(lineInput[4][1])            # time step adjustment factor
        self.tEnd = float(lineInput[5][1])           # end time of simulation
        self.qGWMin = float(lineInput[6][1])         # minimum groundwater discharge (L/T) to qualify 
        print('Read model parameters.')
        

class Source:       # external source term

    def __init__(self, t, Qext):
        self.tEnd = array(t)                     # respective end time for each flux (array)    
        self.Q = Qext                       # volumetric fluxes (array)
 
    def Flux(self, t):
        tIndex = sum(t > self.tEnd)
        return self.Q[tIndex]


class Connection:   # reach connection

    def __init__(self, r1, r2, reach):
        self.r1 = r1                                                # connected reaches
        self.r2 = r2
        self.n = 0.5*reach[r1].n + 0.5*reach[r2].n      # average Manning coefficient       
        self.w = 0.5*reach[r1].w + 0.5*reach[r2].w      # average stream width
        self.Q = 0.                                                 # flux from up-reach to down (placeholder)
        
    def Manning(self, reach, dt, params):
        # instantaneous volumetric flux from r1 into r2, based on h at start of time step
        S = ((reach[self.r1].h+reach[self.r1].z) - (reach[self.r2].h+reach[self.r2].z)) \
            / (0.5*reach[self.r1].L + 0.5*reach[self.r2].L)       # water surface slope
        dirS = sign(S)
        hMean = 0.5*reach[self.r1].h + 0.5*reach[self.r2].h       # water depth at interface
        Qmag = (1./self.n) * (hMean*self.w) * (hMean*self.w/(2.*hMean + self.w))**0.6667 * abs(S)**0.5
        if S > 0.:
            Q = dirS * Qmag * (reach[self.r1].h > params.hMin)       # flow from r1 to r2
        else:
            Q = dirS * Qmag * (reach[self.r2].h > params.hMin)       # flow from r2 to r1
        return Q

    
class Reach:        # rectangular subsection along stream

    def __init__(self, x, y, z, L, w, orient, h, n, K, fixed, connect, dirx):
        self.x = x                  # location of centerpoint
        self.y = y
        self.z = z                  # mean stream bed elevation       
        self.L = L                  # length
        self.w = w                  # width
        self.floorArea = L * w      # area of reach footprint
        self.orient = orient        # directional orientation, radians E of due N
        self.h = h                  # water depth
        self.n = n                  # Manning roughness coefficient
        self.K = K                  # hydraulic conductivity of stream bed  #debug
        self.fixed = fixed          # true = fixed-condition boundary / false = flow-thru reach
        self.connect = connect      # list of connection indices associated with this reach
        self.dirx = dirx            # -1 if reach is r1; +1 if reach is r2
        self.qGW = 0.               # instantaneous discharge rate to groundwater (L/T)
        self.qGWCum = 0.            # time-weighted cumulative discharge to groundwater (L)
        self.tGW0 = 0.              # starting and ending times for discharge period to GW
        self.tGWf = 0.
        self.qGWavg = 0.            # average discharge to groundwater (L/T) between tGW0 and tGWf
        self.l = 0.5 * L            # source area half-lengths (to be passed to aquifer mounding model)
        self.a = 0.5 * w
        
    def Drain(self, dt, t, params):
        # groundwater discharge bookkeeping
        hDrain = min(self.K*dt, 0.999*self.h)           # height of column of water discharged to groundwater over dt
        qGW = hDrain/dt                                 # implied groundwater recharge rate
        if hDrain > params.qGWMin:
            if (self.tGW0==0.): self.tGW0 = t-dt
            self.tGWf = t
            self.h -= hDrain                            # adjust water column height in response to discharge
            self.qGW = qGW                        
            self.qGWCum += hDrain
        
        
### utility functions ###

       
def VolBal(rch, indexRch, connection, source, srcList, t, dt):     
    # implied water level change via explicit solution of volume balance equation
    Qtot = 0.
    for i, cn in enumerate(rch.connect):          # stream fluxes
        Qtot += rch.dirx[i] * connection[cn].Q
    if indexRch in srcList:         # external source at time t, if applicable
        i = srcList.index(indexRch)
        Qtot += source[i].Flux(t)
    dh = dt * Qtot/rch.floorArea
    return dh
 
 
def Updates(reach, dhNew, dt, t, params):
    # update water levels
    for i, rch in enumerate(reach):
        rch.h += dhNew[i]                   # update water levels (1st pass)
        rch.Drain(dt, t, params)            # address any discharge to groundwater        

        
def ReadReaches():
    # read in table of reach properties
    reach = []
    lineInput = []        
    inputFile = open('reaches.csv','r')
    for line in inputFile: lineInput.append(line.split(','))
    inputFile.close()    
    for i in range(1, len(lineInput)):
        x = float(lineInput[i][0])
        y = float(lineInput[i][1])
        z = float(lineInput[i][2])
        L = float(lineInput[i][3])
        w = float(lineInput[i][4])
        orient = float(lineInput[i][5])
        h = float(lineInput[i][6])
        n = float(lineInput[i][7])
        K = float(lineInput[i][8])
        fixed = (lineInput[i][9].strip()=='TRUE')
        reach.append(Reach(x, y, z, L, w, orient, h, n, K, fixed, [], []))
    print('Read reach table.')
    return reach


def ProcessConnects(reach):
    # read connections table; updated reach connects lists
    connection = []
    lineInput = []        
    inputFile = open('connections.csv','r')
    for line in inputFile: lineInput.append(line.split(','))
    inputFile.close()
    for i in range(1, len(lineInput)):
        r1 = int(lineInput[i][0])
        r2 = int(lineInput[i][1])
        connection.append(Connection(r1, r2, reach))
        reach[r1].connect.append(i-1)
        reach[r1].dirx.append(-1)        
        reach[r2].connect.append(i-1)
        reach[r2].dirx.append(1)        
    print('Read connection table.')    
    return connection, reach


def ReadSources():
    # source term flux tables, by affected reach
    source = []
    sources_df = read_csv('sources.txt', sep='\t')
    srcList = list(set(list(sources_df['reach'])))
    for src in srcList: 
        tEnd = list(sources_df[sources_df['reach']==src]['t'])
        Qext = list(sources_df[sources_df['reach']==src]['Q'])
        source.append(Source(tEnd, Qext))
    print('Read source term information.') 
    return srcList, source


def ReadPrints():
    # read in output time list
    printTimes = []
    lineInput = []        
    inputFile = open('print_times.txt','r')
    for line in inputFile: lineInput.append(line.split())
    inputFile.close()
    for i in range(len(lineInput)): printTimes.append(float(lineInput[i][0]))   
    return printTimes

    
def WriteReaches(reach):
    # summarize reach properties
    leakyReach = []
    reachFile = open('reaches_summary.csv','w')
    reachFile.writelines(['reach', ',', 'x', ',', 'y', ',', 'L', ',', 'w', ',', 
        'orient', ',', 'q', ',', 'tStart', ',', 'tEnd'])
    reachFile.writelines(['\n'])
    for i, rch in enumerate(reach):
        tSpan = rch.tGWf - rch.tGW0
        if tSpan > 0.:
            rch.qGWavg = rch.qGWCum/tSpan
            leakyReach.append(rch)
            reachFile.writelines([str(i), ',', str(rch.x), ',', str(rch.y), ',', str(rch.L), ',', str(rch.w), ',',
                str(rch.orient), ',', str(rch.qGWavg), ',', str(rch.tGW0), ',', str(rch.tGWf)])
            reachFile.writelines(['\n'])
    reachFile.close() 
    return leakyReach
    

def CheckConvg(reach, dhNew, params):
    h0 = zeros(len(reach), float)
    for i, rch in enumerate(reach): h0[i] = rch.h
    c1 = (h0 + dhNew < 0.)
    c2 = (abs(dhNew) >= params.dhMax)
    if sum(c1) + sum(c2) > 0:
        convFlag = False
    else:
        convFlag = True
    return convFlag
   
    
def Run(mode):

    print('***** Stream flow model. *****')
    params = Params()                               # read model parameters
    reach = ReadReaches()                           # read reach features
    connection, reach = ProcessConnects(reach)      # read reach connection table; populate reach connection lists
    srcList, source = ReadSources()                 # read external source terms
    printTimes = ReadPrints()                       # list of output times
    
    # initialize
    printTimes.append(params.tEnd)    
    printTimes = array(printTimes)    
    dhNew = zeros(len(reach), float)
    qGW = zeros(len(reach), float)
    t = 0.
    dt = params.dt0
    tally = Tally()
    
    # open reach output file
    outputFile = open('reach_output.csv','w')
    outputFile.writelines(['time', ',', 'reach', ',', 'x', ',', 'y', ',', 'z', ',', 'h', ',', 'qGW'])
    outputFile.writelines(['\n'])
    
    # explicit time stepping through simulation
    while t < params.tEnd:
        currentForcedTime = printTimes[sum(t >= printTimes)]
        convFlag = False
        while convFlag == False:
            for cn in connection: cn.Q = cn.Manning(reach, dt, params)      # stream flux estimates
            for indexRch, rch in enumerate(reach):
                if not rch.fixed:
                    dhNew[indexRch] = VolBal(rch, indexRch, connection, source, srcList, t, dt)         # stream flux balance
            convFlag = CheckConvg(reach, dhNew, params)          # check convergence; return time step reduction factor, if applicable
            if convFlag == False:
                dt /= params.dtf
                assert(dt > params.dt0)                                           # kill simulation for obvious iteration failure
            else:
                t += dt
                Updates(reach, dhNew, dt, t, params)           # update reach object array variables                
                tally.Balance(reach, connection, source, srcList, t, dt)         # update volume balances
                dtPrior = dt                
                dt = min(params.dtf*dt, currentForcedTime-t, params.dtMax)        # select new time step
        if t in printTimes:
            print('t =', str(t))
            for i, rch in enumerate(reach):                             # append results to output file
                outputFile.writelines([str(t), ',', str(i), ',', str(rch.x), ',', str(rch.y), ',', str(rch.z), ',', str(rch.h), ',', str(rch.qGW)])
                outputFile.writelines(['\n'])      
            dt = dtPrior
    
    outputFile.close()
    tally.WriteBalances()
    leakyReach = WriteReaches(reach)
    
    if mode:
        print('***** Groundwater mounding model. *****')
        aquifer = Aquifer()                                               # read aquifer parameters    
        wellSet = WellSet(aquifer)                                        # read well locations and history delineation
        wellSet.TimeSeries(aquifer, leakyReach)
        wellSet.WriteTimeSeries()
        wellSet.PlotTimeSeries()
    
    print('Finished.')

    
### run script ###

Run(1)  # run mode: 1 = run groundwater model