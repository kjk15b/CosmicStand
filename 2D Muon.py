# This program should simulate and graph muons interacting with the detector

import math
import random as ran
from math import atan, degrees
import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


def muonEvent(maxCount):
    class GrowingList(list):
        def __setitem__(self, index, value):
            if index >= len(self):
                self.extend([None] * (index + 1 - len(self)))
            list.__setitem__(self, index, value)
    count = 0
    scintillators = [0, 0, 0, 0, 0, 0, 0, 0]
    topFlag = 0
    botFlag = 0
    topEvent = 0
    botEvent = 0
    coincidence = 0
    topXpos = GrowingList()
    topYpos = GrowingList()
    botXpos = GrowingList()
    botYpos = GrowingList()
    slope = GrowingList()
    gamma = GrowingList()
    i = 0
    while i <= maxCount:
        slope[i] = 0
        gamma[i] = 0
        topXpos[i] = 0
        topYpos[i] = 0
        botXpos[i] = 0
        botYpos[i] = 0
        i += 1
    while count < maxCount:
        topMuonX = ran.randint(0, 200)
        topMuonY = ran.randint(0, 200)
        botMuonX = ran.randint(0, 200)
        botMuonY = ran.randint(0, 200)
        muDir = ran.randint(0, 100)
        topXpos[count] = topMuonX
        topYpos[count] = topMuonY
        botXpos[count] = botMuonX
        botYpos[count] = botMuonY
        if (10 <= topMuonX <= 81) & (10 <= topMuonY <= 34):
            topFlag = 1
        if topFlag == 1:
            topEvent += 1
            if 10 <= topMuonY <= 16:
                scintillators[7] += 1
            elif 16 < topMuonY <= 22:
                scintillators[6] += 1
            elif 22 < topMuonY <= 28:
                scintillators[5] += 1
            elif 28 < topMuonY <= 34:
                scintillators[4] += 1
        if (10 <= botMuonX <= 81) & (10 <= botMuonY <= 34):
            botFlag = 1
        if botFlag == 1:
            botEvent += 1
            if 10 <= botMuonY <= 16:
                scintillators[3] += 1
            elif 16 < botMuonY <= 22:
                scintillators[2] += 1
            elif 22 < botMuonY <= 28:
                scintillators[1] += 1
            elif 28 < botMuonY <= 34:
                scintillators[0] += 1
        if (botFlag & topFlag) == 1:
            coincidence += 1
        if muDir > 15:
            if (botMuonX - topMuonX) == 0:
                slope[count] = "undefined"
            elif (botMuonX - topMuonX) != 0:
                slope[count] = (botMuonY - topMuonY) / (botMuonX - topMuonX)
            gamma[count] = degrees(math.atan((botMuonY - topMuonY) / 48))
        if 0 <= muDir <= 15:
            if (topMuonX - botMuonX) == 0:
                slope[count] = "undefined"
            elif (topMuonX - botMuonX) != 0:
                slope[count] = (topMuonY - botMuonY) / (topMuonX - botMuonX)
            gamma[count] = degrees(math.atan((topMuonY - botMuonY) / 48))
        count += 1
    print("Top Events Occured:", "\n", topEvent, "\n", "Bottom Events Occured:", "\n", botEvent, "\n", "Coincidences",
          "\n", coincidence, "\n", "Hits on Each Scintillator:", "\n", scintillators, "\n", "Slopes of Muons:", "\n",
          slope, "\n", "Incident Angles:", "\n", gamma, "\n")
    today = datetime.datetime.now()
    txt = open("Monte Carlo Sim.txt",'w')
    txt.write("\n")
    txt.write("Current Run Time:")
    txt.write("\n")
    txt.write(str(today))
    txt.write("\n")
    txt.write("Top Events Triggered:")
    txt.write("\n")
    txt.write(str(topEvent))
    txt.write("\n")
    txt.write("Bottom Events Triggered:")
    txt.write("\n")
    txt.write(str(botEvent))
    txt.write("\n")
    txt.write("Number of Coincidences:")
    txt.write("\n")
    txt.write(str(coincidence))
    txt.write("\n")
    txt.write("Hits on Top Plane Scintillators:")
    txt.write("\n")
    i = 4
    while i <= 7:
        txt.write(str(scintillators[i]))
        txt.write("\n")
        i += 1
    txt.write("Hits on Bottom Plane Scintillators:")
    txt.write("\n")
    i = 0
    while i <= 3:
        txt.write(str(scintillators[i]))
        txt.write("\n")
        i += 1
    txt.write("Slopes of Individual Cosmic Rays:")
    txt.write("\n")
    i = 0
    while i <= maxCount:
        txt.write(str(slope[i]))
        txt.write("\n")
        i += 1
    txt.write("Incident Angles of Cosmic Rays:")
    txt.write("\n")
    i = 0
    while i <= maxCount:
        txt.write(str(gamma[i]))
        txt.write("\n")
        i += 1
    txt.close()

    x = np.arange(8)
    plt.figure(1)
    plt.bar(x, scintillators)
    plt.xticks(x,('1', '2', '3', '4', '5', '6', '7', '8'))
    plt.xlabel('Detector')
    plt.ylabel('Events')
    plt.title('Combined')
    plt.figure(2)
    plt.hist2d(topXpos, topYpos, bins=maxCount, norm=LogNorm())
    plt.colorbar()
    plt.show()


muonEvent(10000)
