import sensors
reload(sensors)
from soar.graphing import *
from soar.widgets import *
import pickle
import math

def distToWall(sensorA, sensorB):
   num = sensorA**2*sensorB**2*(math.sin(math.pi/7))**2
   denom = sensorB**2 + sensorA**2 - 2*sensorA*sensorB*math.cos(math.pi/7)
   return math.sqrt(num/denom)


def getLR():
##    (sonarDistances, pose)= sensors
   left = distToWall(sonarDistances()[0], sonarDistances()[1])
   right = distToWall(sonarDistances()[6], sonarDistances()[7])
   #if one sensor doesnt return a value, then use 5.0 and the other reading value.
   if sonarDistances()[6] == 5.0 and sonarDistances()[7] == 5.0 and sonarDistances()[0] == 5.0 and sonarDistances()[1] == 5.0:
       minDist = min(sonarDistances())

       i = sonarDistances().index(minDist)
       if i <= 3:
           left = minDist
           right =6.0
       else:
           left = 6.0
           right = minDist

       if minDist == 5:
           left = 6.0 #max distance estimated is set to arbitrary number of 6
           right = 6.0

   elif (sonarDistances()[6] == 5.0 and sonarDistances()[7] == 5.0) or (sonarDistances()[0] == 5.0 and sonarDistances()[1]):
       if sonarDistances()[6] == 5.0:
           right = 6.0
       else:
           left = 6.0

   return [left, right]

#Question 4
def getLRT():
   ratioLeft = getLR()[0]/sonarDistances()[0]
   ratioRight = getLR()[1] / sonarDistances()[7]
   minRatio = min([ratioLeft, ratioRight])
   theta = math.acos(minRatio)
   z = getLR()
   z.append(theta)
   return z



runLength = 100
makeGraph = True
showSensors = True
d_desired = 0

proportional = 1
K = 1

proportionalPlusDelay = 2
K1 = 36
K2 = -34

proportionalPlusAngle = 3
K3 = 0
K4 = 0

# select a controller
controller = proportionalPlusDelay

def setup():
    if showSensors:
        robot.w = DrawingWindow(300, 300, -1, 1, -1, 1, "draw")
    else:
        robot.w = False
    robot.e = 0
    robot.ee = 0
    robot.step = 0
    robot.D = []

def step():
#    l,r = sensors.getLR(robot.w)
    print getLR()
    l,r,t = getLRT()
    # l = distance to left wall (meters)
    # r = distance to right wall (meters)
    # t = angle of robot (radians)
    d = 0.5*(r-l)

    if controller==proportional:
        angularVelocity = K*(d_desired-d)

    if controller==proportionalPlusDelay:
        robot.ee = robot.e
        robot.e = d_desired-d
        angularVelocity = K1*robot.e+K2*robot.ee

    if controller==proportionalPlusAngle:
        pass

    motorOutput(.1,angularVelocity)
    if makeGraph:
        # Gather data and make a graph at the end
        robot.D.append(l)
        robot.step += 1
        if robot.step == runLength:
            (minVal, maxVal) = (min(robot.D), max(robot.D))
            rangeIncr = (maxVal - minVal)*0.05
            g = GraphingWindow(500, 300, 0, runLength,
                               minVal-rangeIncr, maxVal+rangeIncr, "data")
            # Graph the data values as a function of their indices
            g.graphDiscrete(lambda x: robot.D[x])

            # Save the data into a file
            f = open("leftValues.py", "w")
            pickle.dump(robot.D, f)
            f.close()
