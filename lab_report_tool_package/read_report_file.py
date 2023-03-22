##########################################################################################################
# Copyright: QIN Qijun
# 
# This document is the property of Qin Qijun, please cite the source if you need to quote it.
# Workshop: https://github.com/QuintinUmi
# Contact:  qqj030212@gmail.com
##########################################################################################################

import os
import sys

##
# Author QIN Qijun
#
def read_file_split_data(filePath):

    if(not os.path.isfile(filePath)):
        return -1

    f = open(filePath, "r")

    dataLine = []
    while(True):

        dataStr = f.readline()
        if(dataStr == ""):
            break

        dataSplit = dataStr.split()
        # dataSplit = dataSplit[0: len(dataSplit) - 1]
        isdata = True
        # print(dataSplit)

        if(len(dataSplit) == 0):
            isdata = False
        for dataCheck in dataSplit:
            if(not isnumber(dataCheck)):
                isdata = False
            else:
                dataSplit[dataSplit.index(dataCheck)] = float(dataCheck)

        if(isdata):
            dataLine.append(dataSplit)
            # print("digit: ", dataSplit)

    return dataLine

def get_colume_data(dataLine, colume):

    dataColume = []
    for datum in dataLine:
        dataColume.append(datum[colume])

    return dataColume

def isnumber(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
 
    return False


if __name__ == "__main__":
    data = read_file_split_data("C:\\Users\\qqj03\\Desktop\\Lab Result\\G04_Acrylic.txt")
    print(data)