from mpl_toolkits.mplot3d import Axes3D
import scipy as scp
import numpy as np
import time
import datetime
#from datetime import timedelta
import matplotlib.pyplot as plt

print("libraries imported")

def getTime(mode):
    global start_time0, start_time, end_time0, end_time
    match mode:
        case 'integrateStart':
            start_time0 = time.time()
            start_time = datetime.datetime.fromtimestamp(start_time0)
            return start_time
        case 'integrateEnd':
            end_time0 = time.time()
            end_time = datetime.datetime.fromtimestamp(end_time0)
            return end_time
        case 'current':
            current_time0 = time.time()
            current_time = datetime.datetime.fromtimestamp(current_time0)
            return current_time0, current_time
def printAlligned(message, variable, unit, msgspaces=20, varspaces=30):
    output = print("{:{msgspace}s}= {:{varspace}s} {}".format(message, str(variable), unit, msgspace=msgspaces, varspace=varspaces))
    return output
def debug(level='partial'):
    match level:
        case 'partial':
            print('   DEBUG LEVEL PARTIAL')
            print("Iteration {} at current time = {} seconds".format(stepN, t))
            printAlligned('S position',iniSpos, 'm')
            printAlligned('S velocity', iniSvel, 'm s^-1')
        case 'full':
            print('   DEBUG LEVEL FULL')
            print("List of globals:")
            print("Iteration {} at current time = {} seconds".format(stepN, t))
            printAlligned('S position',iniSpos, 'm')
            printAlligned('S velocity', iniSvel, 'm s^-1')
        case 'constants':
            print('   DEBUG LEVEL CONSTANTS')
            localT = getTime('current')[1]
            printAlligned('Current time', localT, '')
        case 'start':
            print('   DEBUG LEVEL START')
            localT = getTime('integrateStart')
            printAlligned('Start time', localT, 'seconds', msgspaces=15)
            printAlligned('S position',iniSpos, 'm')
            printAlligned('S velocity', iniSvel, 'm s^-1')
        case 'end':
            print('   DEBUG LEVEL END')
            end_time = time.time()
            printAlligned('end_time',end_time, 'seconds', msgspaces=15)
            printAlligned('Wall runtime',end_time-start_time, msgspaces=15)
        case 'none':
            return None

def initialize(mode,inSpos=[0,0],inSvel=[0,0],startdebugmode='start'): 
    # Constants
    match mode:
        case 'constants':
            global ConstantsList, GravConst, dEM, mE, mM, mT, rE, rM, PeriodEM, l2Original, l2VelOriginal, drOriginal
            GravConst = 6.6726e-11
            dEM = 3.844e8                                     #Earth-Moon separation
            mE = 5.9742e24                                  #Earth Mass
            mM = 7.35e22                                    #Moon mass
            mT = mE + mM                                    #Barycentre mass
            rE = dEM*mM/mT                                    #Barycentre-Earth separation
            rM = dEM*mE/mT                                    #Barycentre-Moon separation
            PeriodEM = np.sqrt(4 * np.pi**2 * dEM**3 / GravConst / mT)    #Orbital period
            l2Original = rM + dEM*(mM/3/mE)**(1/3)
            l2VelOriginal = 2*np.pi * l2Original / PeriodEM
            drOriginal = dEM*(mM/3/mE)**(1/3)
            ConstantsList = [GravConst, dEM, mE, mM, mT, rE, rM, PeriodEM, l2Original, l2VelOriginal, drOriginal]
            return None
        case 'inputs':
            global iniSpos, iniSvel
            iniSpos = np.array(inSpos)
            iniSvel = np.array(inSvel)
            debug(startdebugmode)
            return None
    return 'FAILED TO INITIALIZE'
##initial velocity ???
initialize('constants')


'''unused debug
print(start_time)
printAlligned('S position',Spos)
print("{:20s}= {}".format('S position',Spos))
print("S velocity     = {}".format(Svel))
print("S acceleration = {}".format(Sacc))
print("P = {}".format(P))
print("HR P = {}".format(timedelta(seconds = P)))
'''

def calcEpos(time, dim=2, constantsList=[]):
    # GravConst, dEM, mE, mM, mT, rE, rM, PeriodEM, l2Original, l2VelOriginal, drOriginal
    loc_rE = constantsList[5]
    loc_PeriodEM = constantsList[7]
    match dim:
        case 2:
            output = np.array([-loc_rE*np.cos(2*np.pi * time/loc_PeriodEM),-loc_rE*np.sin(2*np.pi * time/loc_PeriodEM)])
            return output
        case 3:
            output = np.array([-loc_rE*np.cos(2*np.pi * time/loc_PeriodEM),-loc_rE*np.sin(2*np.pi * time/loc_PeriodEM),0])
            return output
        case _:
            print('ERROR: UNKNOWN DIMENSION')
            return None
def calcMpos(time, dim=2, constantsList=[]):
    #GravConst, dEM, mE, mM, mT, rE, rM, PeriodEM, l2Original, l2VelOriginal, drOriginal
    loc_rM = constantsList[6]
    loc_PeriodEM = constantsList[7]
    match dim:
        case 2:
            output = np.array([loc_rM*np.cos(2*np.pi * time/loc_PeriodEM),loc_rM*np.sin(2*np.pi * time/loc_PeriodEM)])
            return output
        case 3:
            output = np.array([loc_rM*np.cos(2*np.pi * time/loc_PeriodEM),loc_rM*np.sin(2*np.pi * time/loc_PeriodEM),0])
            return output
        case _:
            print('ERROR: UNKNOWN DIMENSION')
            return None
def d2t_r(time, inPos, constantsList=[]):
    # GravConst, dEM, mE, mM, mT, rE, rM, PeriodEM, l2Original, l2VelOriginal, drOriginal
    dim = len(inPos)
    loc_GravConst= constantsList[0]
    loc_mE = constantsList[2]
    loc_mM = constantsList[3]
    Mpos = calcMpos(time, dim, constantsList)
    Epos = calcEpos(time, dim, constantsList)
    SumE = 0
    SumM = 0
    for i in range(dim):
        SumE = SumE + (inPos[i] - Epos[i])**2
        SumM = SumM + (inPos[i] - Mpos[i])**2
    dE = np.sqrt(SumE)
    dM = np.sqrt(SumM)
    Sacc = []
    for i in range(dim):
        Sacc.append(-loc_GravConst*loc_mE*(inPos[i] - Epos[i])/dE**3 - loc_GravConst*loc_mM*(inPos[i] - Mpos[i])/dM**3)
    Sacc = np.array(Sacc)
    return Sacc

def TaylorODEsolve3BP(mode, t_0, t_f, dt, limits=[0,0], progress=True):
    LocalConstantsList = ConstantsList
    match mode:
        case 'solve':
            Spos = iniSpos
            Svel = iniSvel
            steps = int((t_f - t_0) / dt)
            tList = np.arange(t_0, t_f, dt)
            SposList = np.zeros((steps + 1, 2))
            EposList = np.zeros((steps + 1, 2))
            MposList = np.zeros((steps + 1, 2))
            
            t = t_0
            SposList[0] = Spos
            EposList[0] = calcEpos(t_0, 2, LocalConstantsList)
            MposList[0] = calcMpos(t_0, 2, LocalConstantsList)
            for i in range(steps):
                match progress:
                    case True:
                        if i% (int(steps/10)) == 0:
                            n = round(i / steps, 2) * 100
                            printAlligned('{}%'.format(n), Spos,'m')
                            printAlligned('{}%'.format(n), Svel,'m s^-1')
                    case False:
                        pass
                Sacc = d2t_r(t, Spos, LocalConstantsList)
                Spos = Spos + dt * Svel + Sacc * dt**2 / 2
                Svel = Svel + dt * Sacc
                t = tList[i+1]
                SposList[i+1] = Spos
                EposList[i+1] = calcEpos(t, 2, LocalConstantsList)
                MposList[i+1] = calcMpos(t, 2, LocalConstantsList)

            SposList = np.transpose(SposList)
            EposList = np.transpose(EposList)
            MposList = np.transpose(MposList)
            return SposList, MposList, EposList, tList
        case 'limits':
            Spos = iniSpos
            Svel = iniSvel
            steps = int((t_f - t_0) / dt)
            tList = np.arange(t_0, t_f, dt)
            SposList = np.zeros((steps + 1, 2))
            EposList = np.zeros((steps + 1, 2))
            MposList = np.zeros((steps + 1, 2))
            
            t = t_0
            SposList[0] = Spos
            EposList[0] = calcEpos(t_0, 2, LocalConstantsList)
            MposList[0] = calcMpos(t_0, 2, LocalConstantsList)
            for i in range(steps):
                '''#Progress bar, not used
                if i% (int(steps/10)) == 0:
                    n = round(i / steps, 2) * 100
                    printAlligned('{}%'.format(n), Spos,'m')
                    printAlligned('{}%'.format(n), Svel,'m s^-1')
                '''
                
                Sacc = d2t_r(t, Spos, LocalConstantsList)
                Spos = Spos + dt * Svel + Sacc * dt**2 / 2
                Svel = Svel + dt * Sacc
                t = tList[i+1]
                SposList[i+1] = Spos
                EposList[i+1] = calcEpos(t, 2, LocalConstantsList)
                MposList[i+1] = calcMpos(t, 2, LocalConstantsList)

                Sdist = np.sqrt(Spos[0]**2 + Spos[1]**2)
                if Sdist < limits[0]:
                    move = 1
                    SposList = np.transpose(SposList)
                    EposList = np.transpose(EposList)
                    MposList = np.transpose(MposList)
                    return SposList, MposList, EposList, tList, move, i
                if Sdist > limits[1]:
                    move = -1
                    SposList = np.transpose(SposList)
                    EposList = np.transpose(EposList)
                    MposList = np.transpose(MposList)
                    return SposList, MposList, EposList, tList, move, i

            SposList = np.transpose(SposList)
            EposList = np.transpose(EposList)
            MposList = np.transpose(MposList)
            move = 0
            return SposList, MposList, EposList, tList, move, steps

def RK4ODEsolve3BP(mode, t_0, t_f, dt, limits=[0,0], progress=True): #t_0: time initial [seconds], t_f: time final [seconds], dt: time step
    LocalConstantsList = ConstantsList
    match mode:
        case 'solve':
            # GravConst, dEM, mE, mM, mT, rE, rM, PeriodEM, l2Original, l2VelOriginal, drOriginal
            Spos = iniSpos
            Svel = iniSvel
            steps = int((t_f - t_0) / dt)
            tList = np.arange(t_0, t_f, dt)
            SposList = np.zeros((steps + 1, 2))
            EposList = np.zeros((steps + 1, 2))
            MposList = np.zeros((steps + 1, 2))
            
            t = t_0
            SposList[0] = Spos
            EposList[0] = calcEpos(t_0, 2, LocalConstantsList)
            MposList[0] = calcMpos(t_0, 2, LocalConstantsList)
            for i in range(steps):
                w0 = d2t_r(t, Spos, LocalConstantsList)
                u1 = Spos + dt * Svel / 2 #u, v, w denote position, velocity, and acceleration (respectively) of test points
                v1 = Svel + dt * w0 / 2
                w1 = d2t_r(t + dt/2, u1, LocalConstantsList)
                u2 = Spos + dt * v1 / 2
                v2 = Svel + dt * w1 / 2
                w2 = d2t_r(t + dt/2, u2, LocalConstantsList)
                u3 = Spos + dt * v2
                v3 = Svel + dt * w2
                w3 = d2t_r(t + dt, u3, LocalConstantsList)
                
                match progress:
                    case True:
                        if i% (int(steps/10)) == 0:
                            n = round(i / steps, 2) * 100
                            printAlligned('{}%'.format(n), Spos,'m')
                            printAlligned('{}%'.format(n), Svel,'m s^-1')
                    case False:
                        pass
                
                Spos = Spos + dt * (Svel + 2*v1 + 2*v2 + v3) / 6
                Svel = Svel + dt * (w0 + 2*w1 + 2*w2 + w3) / 6
                t = tList[i+1]
                SposList[i+1] = Spos
                EposList[i+1] = calcEpos(t, 2, LocalConstantsList)
                MposList[i+1] = calcMpos(t, 2, LocalConstantsList)

            SposList = np.transpose(SposList)
            EposList = np.transpose(EposList)
            MposList = np.transpose(MposList)
            return SposList, MposList, EposList, tList
        case 'limits':
            Spos = iniSpos
            Svel = iniSvel
            steps = int((t_f - t_0) / dt)
            tList = np.arange(t_0, t_f, dt)
            SposList = np.zeros((steps + 1, 2))
            EposList = np.zeros((steps + 1, 2))
            MposList = np.zeros((steps + 1, 2))
            
            t = t_0
            SposList[0] = Spos
            EposList[0] = calcEpos(t_0, 2, LocalConstantsList)
            MposList[0] = calcMpos(t_0, 2, LocalConstantsList)
            for i in range(steps):
                w0 = d2t_r(t, Spos, LocalConstantsList)
                u1 = Spos + dt * Svel / 2 #u, v, w denote position, velocity, and acceleration (respectively) of test points
                v1 = Svel + dt * w0 / 2
                w1 = d2t_r(t + dt/2, u1, LocalConstantsList)
                u2 = Spos + dt * v1 / 2
                v2 = Svel + dt * w1 / 2
                w2 = d2t_r(t + dt/2, u2, LocalConstantsList)
                u3 = Spos + dt * v2
                v3 = Svel + dt * w2
                w3 = d2t_r(t + dt, u3, LocalConstantsList)
                
                '''#progress bar, not in use
                if i% (int(steps/10)) == 0:
                    n = round(i / steps, 2) * 100
                    printAlligned('{}%'.format(n), Spos,'m')
                    printAlligned('{}%'.format(n), Svel,'m s^-1')
                '''
                
                Spos = Spos + dt * (Svel + 2*v1 + 2*v2 + v3) / 6
                Svel = Svel + dt * (w0 + 2*w1 + 2*w2 + w3) / 6
                Sdist = np.sqrt(Spos[0]**2 + Spos[1]**2)
                
                t = tList[i+1]
                SposList[i+1] = Spos
                EposList[i+1] = calcEpos(t, 2, LocalConstantsList)
                MposList[i+1] = calcMpos(t, 2, LocalConstantsList)
                if Sdist < limits[0]:
                    move = 1
                    SposList = np.transpose(SposList)
                    EposList = np.transpose(EposList)
                    MposList = np.transpose(MposList)
                    return SposList, MposList, EposList, tList, move, i
                if Sdist > limits[1]:
                    move = -1
                    SposList = np.transpose(SposList)
                    EposList = np.transpose(EposList)
                    MposList = np.transpose(MposList)
                    return SposList, MposList, EposList, tList, move, i

            SposList = np.transpose(SposList)
            EposList = np.transpose(EposList)
            MposList = np.transpose(MposList)
            move = 0
            return SposList, MposList, EposList, tList, move, steps
        case '3dSolve':

            Spos = iniSpos
            Svel = iniSvel
            steps = int((t_f - t_0) / dt)
            tList = np.arange(t_0, t_f, dt)
            SposList = np.zeros((steps + 1, 3))
            EposList = np.zeros((steps + 1, 3))
            MposList = np.zeros((steps + 1, 3))
            
            t = t_0
            SposList[0] = Spos
            EposList[0] = calcEpos(t_0, 3, LocalConstantsList)
            MposList[0] = calcMpos(t_0, 3, LocalConstantsList)
            for i in range(steps):
                w0 = d2t_r(t, Spos, LocalConstantsList)
                u1 = Spos + dt * Svel / 2 #u, v, w denote position, velocity, and acceleration (respectively) of test points
                v1 = Svel + dt * w0 / 2
                w1 = d2t_r(t + dt/2, u1, LocalConstantsList)
                u2 = Spos + dt * v1 / 2
                v2 = Svel + dt * w1 / 2
                w2 = d2t_r(t + dt/2, u2, LocalConstantsList)
                u3 = Spos + dt * v2
                v3 = Svel + dt * w2
                w3 = d2t_r(t + dt, u3, LocalConstantsList)
                
                match progress:
                    case True:
                        if i% (int(steps/10)) == 0:
                            n = round(i / steps, 2) * 100
                            printAlligned('{}%'.format(n), Spos,'m')
                            printAlligned('{}%'.format(n), Svel,'m s^-1')
                    case False:
                        pass
                
                Spos = Spos + dt * (Svel + 2*v1 + 2*v2 + v3) / 6
                Svel = Svel + dt * (w0 + 2*w1 + 2*w2 + w3) / 6
                t = tList[i+1]
                SposList[i+1] = Spos
                EposList[i+1] = calcEpos(t, 3, LocalConstantsList)
                MposList[i+1] = calcMpos(t, 3, LocalConstantsList)

            SposList = np.transpose(SposList)
            EposList = np.transpose(EposList)
            MposList = np.transpose(MposList)
            return SposList, MposList, EposList, tList

def test(mode, t_0, t_f,dt, progress=True):
    match mode:
        case '3dSolve':
            Spos = iniSpos
            Svel = iniSvel
            steps = int((t_f - t_0) / dt)
            tList = np.arange(t_0, t_f, dt)
            SposList = np.zeros((steps + 1, 3))
            EposList = np.zeros((steps + 1, 3))
            MposList = np.zeros((steps + 1, 3))
            
            t = t_0
            SposList[0] = Spos
            EposList[0] = calcEpos(t_0)
            MposList[0] = calcMpos(t_0)
            for i in range(steps):
                w0 = d2t_r(t, Spos)
                u1 = Spos + dt * Svel / 2 #u, v, w denote position, velocity, and acceleration (respectively) of test points
                v1 = Svel + dt * w0 / 2
                w1 = d2t_r(t + dt/2, u1)
                u2 = Spos + dt * v1 / 2
                v2 = Svel + dt * w1 / 2
                w2 = d2t_r(t + dt/2, u2)
                u3 = Spos + dt * v2
                v3 = Svel + dt * w2
                w3 = d2t_r(t + dt, u3)
                
                match progress:
                    case True:
                        if i% (int(steps/10)) == 0:
                            n = round(i / steps, 2) * 100
                            printAlligned('{}%'.format(n), Spos,'m')
                            printAlligned('{}%'.format(n), Svel,'m s^-1')
                    case False:
                        pass
                
                Spos = Spos + dt * (Svel + 2*v1 + 2*v2 + v3) / 6
                Svel = Svel + dt * (w0 + 2*w1 + 2*w2 + w3) / 6
                t = tList[i+1]
                SposList[i+1] = Spos
                EposList[i+1] = calcEpos(t)
                MposList[i+1] = calcMpos(t)

            SposList = np.transpose(SposList)
            EposList = np.transpose(EposList)
            MposList = np.transpose(MposList)
            return SposList, MposList, EposList, tList

def plot(data=None,label='',mode='plot', plotLims=[-5e8,5e8,-5e8,5e8]):
    global ax
    match mode:
        case 'plot':
            ax.plot(*data, label=label)
        case 'axes':
            ax = plt.figure().add_subplot()
            plotLims = ax.axis(plotLims)
            ax.set_aspect('equal',adjustable='box')
            print('plot limits:{}'.format(plotLims))
            plt.axhline(0, color='black')
            plt.axvline(0, color='black')
        case 'show':
            ax.legend()
            plt.show()
def set_axes_equal(ax):
    """
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    """

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def findFFactor(mode, fGuess, dPlaces, tolerance = 1e6, t_f = PeriodEM):
    fGuessStr = str(fGuess)
    print('guess: {}'.format('.' in fGuessStr))
    if ('.' in fGuessStr):
        dPlace = len(fGuessStr.split('.')[1])
    else:
        dPlace = 0
    f = fGuess
    fo = f
    #find length of decimal part
    digit = 0
    n = 1 #factor finding iteration number
    flist = []
    flist.append(f)
    pointslist = []
    match mode:
        case 'RK4':
            while dPlace <= dPlaces:
                dr = drOriginal *f
                l2 = rM + dr
                l2Vel = 2*np.pi * l2 / PeriodEM
                limits = [l2 - tolerance, l2 + tolerance]
                initialize('inputs',inSpos=[l2,0],inSvel=[0, l2Vel],startdebugmode='none')
                out = RK4ODEsolve3BP('limits', 0, t_f, 100, limits= limits) #solve in limit mode
                
                move = out[4]
                steps = out[5]
                if move == 1:
                    digit+=1
                    f = fo + digit * 10** (-dPlace)
                    print('{:10s}{:>2.2s}  {:10s}{:10}'.format('move:', str(move), 'steps:',steps))
                elif move == -1:
                    digit-=1
                    f = fo + digit * 10** (-dPlace)
                    print('{:10s}{:>2.2s}  {:10s}{:10}  {:2s}{:18.12f}  {:2s}{:18}'.format('move:', str(move), 'steps:', steps, 'f:', f,'n:',n))
                    dPlace+=1
                    digit = 0
                    fo = f
                elif move ==0:
                    print('{:10s}{:>2.2s}  {:10s}{:10}  {:2s}{:18.12f}  {:2s}{:18}'.format('move:', str(move), 'steps:', steps, 'f:', f,'n:',n))
                    break
                else:
                    print('ERROR')
                    break
                flist.append(f)
                n+=1
        case 'Taylor':
            while dPlace <= dPlaces:
                dr = drOriginal *f
                l2 = rM + dr
                l2Vel = 2*np.pi * l2 / PeriodEM
                limits = [l2 - tolerance, l2 + tolerance]
                initialize('inputs',inSpos=[l2,0],inSvel=[0, l2Vel],startdebugmode='none')
                out = TaylorODEsolve3BP('limits', 0, t_f, 10, limits= limits) #solve in limit mode
                pointslist = [out[i] for i in range(3)]
                move = out[4]
                steps = out[5]
                if move == 1:
                    digit+=1
                    f = fo + digit * 10** (-dPlace)
                    print('{:10s}{:>2.2s}  {:10s}{:10}'.format('move:', str(move), 'steps:',steps))
                elif move == -1:
                    digit-=1
                    f = fo + digit * 10** (-dPlace)
                    print('{:10s}{:>2.2s}  {:10s}{:10}  {:2s}{:18.12f}  {:2s}{:18}'.format('move:', str(move), 'steps:', steps, 'f:', f,'n:',n))
                    dPlace+=1
                    digit = 0
                    fo = f
                elif move ==0:
                    print('{:10s}{:>2.2s}  {:10s}{:10}  {:2s}{:18.12f}  {:2s}{:18}'.format('move:', str(move), 'steps:', steps, 'f:', f,'n:',n))
                    break
                else:
                    print('ERROR')
                    break
                flist.append(f)
                n+=1
    dr = drOriginal *f
    l2 = rM + dr
    l2Vel = 2*np.pi * l2 / PeriodEM
    limits = [l2 - tolerance, l2 + tolerance]
    initialize('inputs',inSpos=[l2,0],inSvel=[0, l2Vel],startdebugmode='none')
    out = RK4ODEsolve3BP('limits', 0, t_f, 100, limits= limits)
    pointslist = [out[i] for i in range(4)]
    print('f={} in {} steps'.format(f,n))
    return f, n, flist, pointslist

fRK4    = 1.048615605975
dr = drOriginal *fRK4
l2 = rM + dr
l2vel = 2*np.pi * l2 / PeriodEM
initialize('inputs', [l2+2.5848769097e1,-2e2,0], [0,l2vel,4e-3])
SposList, MposList, EposList, tList = RK4ODEsolve3BP('3dSolve', 0, 2*PeriodEM, 100, progress=True)

L2PosList = np.zeros(shape=(len(tList),3)) #Stopped to sleep here
L2PosList = np.array([[l2*np.cos(2*np.pi*tList[i]/PeriodEM),l2*np.sin(2*np.pi*tList[i]/PeriodEM),0] for i in range(len(tList))])
print(np.shape(L2PosList))
L2PosList = np.transpose(L2PosList)
print(np.shape(L2PosList))

L2FrameSposList = SposList - L2PosList
L2RotFrameSposList = np.array([np.cos(2*np.pi*tList/PeriodEM)*L2FrameSposList[0]+np.sin(2*np.pi*tList/PeriodEM)*L2FrameSposList[1],-np.sin(2*np.pi*tList/PeriodEM)*L2FrameSposList[0]+np.cos(2*np.pi*tList/PeriodEM)*L2FrameSposList[1],L2FrameSposList[2]])
print(np.shape(L2RotFrameSposList))
#output = np.array([-loc_rE*np.cos(2*np.pi * time/loc_PeriodEM),-loc_rE*np.sin(2*np.pi * time/loc_PeriodEM),0])
fig = plt.figure(1)
ax = fig.add_subplot(projection = '3d', proj_type = 'persp')
ax.plot(*L2RotFrameSposList, label='Satellite Trajectory')
ax.plot(0, 0, marker="o", markersize=1, markeredgecolor="black")
ax.xaxis.get_offset_text().set_visible(False)
ax.yaxis.get_offset_text().set_visible(False)
ax.zaxis.get_offset_text().set_visible(False)
ax.set_xlabel('$Radius$')
ax.set_ylabel('$X$')
ax.set_zlabel('$Y$')

expList = set_axes_equal(ax)
print(expList)
ax.legend()
plt.show()

fig = plt.figure(2)
ax = fig.add_subplot(projection = '3d', proj_type = 'persp')
ax.plot(*L2FrameSposList, label='Satellite Trajectory')
ax.plot(0, 0, marker="o", markersize=1, markeredgecolor="black")
ax.xaxis.get_offset_text().set_visible(False)
ax.yaxis.get_offset_text().set_visible(False)
ax.zaxis.get_offset_text().set_visible(False)
ax.set_xlabel('$X$')
ax.set_ylabel('$Y$')
ax.set_zlabel('$Z$')

expList = set_axes_equal(ax)
print(expList)
ax.legend()
plt.show()

ax = plt.figure(3).add_subplot(projection='3d',proj_type='persp')
ax.plot(*SposList, label='Satellite Trajectory')
ax.plot(*MposList, label='Moon Trajectory')
ax.plot(*EposList, label='Earth Trajectory')
aspect = ax.get_box_aspect()
print(aspect)
ax.legend()
set_axes_equal(ax)
aspect = ax.get_box_aspect()
print(aspect)

plt.show()