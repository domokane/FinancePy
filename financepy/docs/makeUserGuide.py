# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:55:21 2017

@author: Dominic
"""

VERSION=0.15

fileName = "FinPyManualV_" + str(VERSION)
userGuideFileName = "./" + fileName + ".tex"
headFile = "./head.tex"
tailFile = "./tail.tex"
introFile = "./intro.tex"
import glob

################################################################################

def buildHead():
    ''' Start latex file with a header that sets all of the document 
    properties. '''

    f = open(headFile,'r')
    lines = f.readlines()
    f.close()

    f = open(userGuideFileName,'w')
    f.writelines(lines)
    f.close()

################################################################################

def buildTail():
    ''' Add on end latex to close document. '''
    f = open(tailFile,'r')
    lines = f.readlines()
    f.close()

    f = open(userGuideFileName,'a')
    f.writelines(lines)
    f.close()

################################################################################

def buildIntro():
    ''' Add on end latex to close document. '''
    f = open(introFile,'r')
    lines = f.readlines()
    f.close()

    f = open(userGuideFileName,'a')
    f.writelines(lines)
    f.close()

################################################################################

def buildChapter(foldername):
    ''' Parse a folder by loading up all of the modules in that folder that 
    start with Fin. '''

    readMeFile = foldername + "//" + "README.txt"
    f = open(readMeFile,'r')
    lines = f.readlines()
    f.close()

    f = open(userGuideFileName,'a')
    f.write("\n")
    f.writelines("\\chapter{" + foldername + "}")
    f.write("\n")
    f.writelines("\\section{Introduction}")
    f.write("\n")
    f.writelines(lines)
    f.write("\n")
    f.close()

    modules = glob.glob(foldername + "//Fin*.py")

    for module in modules:
        moduleName = module.split("\\")[-1]
        f = open(userGuideFileName,'a')
        f.write("\\newpage\n")
        f.write("\\section{" + moduleName[0:-3] + "}\n")
        f.write("\n")
        f.close()
        parseModule(module)

    modules = glob.glob(foldername + "//TestFin*.py")

    for module in modules:
        moduleName = module.split("\\")[-1]
        f = open(userGuideFileName,'a')
        f.write("\n")
        f.write("\\newpage\n")
        f.write("\\section{" + moduleName[0:-3] + "}\n")
        f.write("\n")
        f.close()
        parseModule(module)

################################################################################

def parseModule(moduleName):

    ''' Parse a module looking for classes, functions and classes for 
    enumerated types. Functions inside classes are parsed inside the class. '''

    f = open(moduleName,'r')
    lines = f.readlines()
    f.close()

    numEnums = 0
    numClasses = 0
    numFunctions = 0
    
    startEnumLines = []
    startClassLines = []
    startFunctionLines = []

    # Module level classes and functions
    numRows = len(lines)

    for rowNum in range(0,numRows):

        line = lines[rowNum]

        if line.find("class",0,5) != -1 and line.find("Enum") != -1:
            numEnums += 1
            startEnumLines.append(rowNum)

        if line.find("class",0,5) != -1 and line.find("Enum") == -1:
            numClasses += 1
            startClassLines.append(rowNum)

        if line.find("def",0,4) != -1:
            numFunctions += 1
            startFunctionLines.append(rowNum)

    startEnumLines.append(numRows)
    startClassLines.append(numRows)
    startFunctionLines.append(numRows)

    print("startClassLines",startClassLines)

    f = open(userGuideFileName,'a')

    for c in range(0,numEnums):
        newLines = parseEnum(lines,startEnumLines[c],startEnumLines[c+1])

        for newLine in newLines:
            f.writelines(newLine)

    for c in range(0,numClasses):
        newLines = parseClass(lines,startClassLines[c],startClassLines[c+1])

        for newLine in newLines:
            f.writelines(newLine)

    for c in range(0,numFunctions):
        newLines = parseFunction(lines,startFunctionLines[c],startFunctionLines[c+1],False)

        for newLine in newLines:
            f.writelines(newLine)

        f.write("\n")

    f.close()

################################################################################

def parseClass(lines,startLine,endLine):
    ''' Parse a Python class consisting of data members and function methods. '''

    n1 = lines[startLine].find("class")

    newLines = []
       
    className = lines[startLine].split(" ")[1]
    className = className.replace(":","")
    className = className.replace("\n","")

    print(className,startLine,endLine)

    newLines.append("\\subsection{Class: " + className + "}")
    newLines.append("\n")

    ##################################################
    # COMMENT AFTER CLASS BUT BEFORE FIRST FUNCTION
    ##################################################

    startCommentRow = startLine
    endCommentRow = endLine
    commentEndLine = endLine

    for rowNum in range(startLine,endLine):
        line = lines[rowNum]
        if line.find(" def ") > 0:
            commentEndLine = rowNum
            print("End Comment Row",rowNum)
            break

    classComment = ""
    startComment = False
    endComment = False

    for rowNum in range(startLine,commentEndLine):
        line = lines[rowNum]
        if line.find("'''") > 0:
            startCommentRow = rowNum
            startComment = True
            startLine = rowNum+1
            break

    for rowNum in range(startLine,commentEndLine):
        line = lines[rowNum]
        if line.find("'''") > 0:
            endCommentRow = rowNum
            endComment = True
            break

    if startComment == True and endComment == False:
        # assume it's a one-line comment
        endCommentRow = startCommentRow
        endComment = True

    if endComment == True:
        print(startCommentRow,endCommentRow)
        for rowNum in range(startCommentRow,endCommentRow+1):
            line = lines[rowNum]
            line = line.replace("'","")
            line = line.replace("_","\_")
            line = line.replace("\n","")
            line = line.replace("#","\#")
            line = line.lstrip()
            classComment += line

    newLines.append(classComment)
    newLines.append("\n")
    newLines.append("\n")

    ##################################################
    # Now get the data members

    newLines.append("\\subsubsection{Class Data Members}\n")

    dataMembers = set()

    for rowNum in range(startLine,endLine):
        row = lines[rowNum]
        row = row.replace(" ","")
        row = row.replace("_","\_")
        row = row.replace("!","")
        row = row.replace("<","")
        row = row.replace(">","")

        n1 = row.find("self.")
        n2 = row.find("=")

        if n1 != -1 and n2 > n1:
            dataMember = row[n1:n2]
            dataMember = dataMember.strip(" ")
            dataMember = dataMember.strip(")")
            dataMembers.add(dataMember)

    if len(dataMembers) > 0:
        newLines.append("\\begin{itemize}\n")

        for dataMember in dataMembers:
            newLines.append("\\item{" + dataMember + "}\n")
        newLines.append("\\end{itemize}")
        newLines.append("\n")
        newLines.append("\n")
    else:
        newLines.append("No data members found.")
        newLines.append("\n")
        newLines.append("\n")

    newLines.append("\\subsubsection{Class Functions}\n")
    newLines.append("\n")

    # Now get the functions
    numClassFunctions = 0
    startClassFunctionLines = []

    print(startLine,endLine)
    for rowNum in range(startLine,endLine):

        line = lines[rowNum]

        if line.find(" def ") != -1:
            numClassFunctions +=1
            startClassFunctionLines.append(rowNum)

    startClassFunctionLines.append(endLine)

    for c in range(0,numClassFunctions):
        newLines += parseFunction(lines,startClassFunctionLines[c],startClassFunctionLines[c+1],True)
        newLines += "\n"

    return newLines

################################################################################

def parseFunction(lines,startLine,endLine,classFlag):
    ''' Given a set of lines and a start line I extract the function definiton 
    and any comment that goes below. 
    TODO: Add parsing of function arguments and any comments.'''
    
    functionLine = lines[startLine]
    leftCol = functionLine.find("def ")
    n2 = functionLine.find("(")
    functionName = functionLine[leftCol+4:n2]
    functionName = functionName.replace("_","\_")

    # Ensure function stops before any class
    for rowNum in range(startLine+1,endLine):
        line = lines[rowNum]
        if line.find("class ") >= 0:
            endLine = rowNum #update start line to after function signature
            break

    # Ensure function stops before any other function
    for rowNum in range(startLine+1,endLine):
        line = lines[rowNum]
        if line.find("def ") >= 0:
            endLine = rowNum #update start line to after function signature
            break

    functionSignature = ""
    for rowNum in range(startLine,endLine):
        line = lines[rowNum]
        functionSignature += str(line)
        if line.find(":") >= 0:
            startLine = rowNum #update start line to after function signature
            break

    functionComment = ""
    startCommentRow = startLine
    endCommentRow = startLine
    startComment = False
    endComment = False

    for rowNum in range(startLine,endLine):
        line = lines[rowNum]
        if line.find("'''") > 0:
            startCommentRow = rowNum
            startComment = True
            startLine = rowNum+1
            break

    for rowNum in range(startLine,endLine):
        line = lines[rowNum]
        if line.find("'''") > 0:
            endCommentRow = rowNum
            endComment = True
            break

    if startComment == True and endComment == False:
        # assume it's a one-line comment
        endCommentRow = startCommentRow
        endComment = True

    if endComment == True:
        print(startCommentRow,endCommentRow)
        for rowNum in range(startCommentRow,endCommentRow+1):
            line = lines[rowNum]
            line = line.replace("_","\_")
            line = line.replace("'","")
            line = line.replace("\n","")
            line = line.replace("#","\#")
            line = line.lstrip()
            functionComment += line

    # LATEX FORMATTING
    if classFlag == True:
        functionDescription = "\\subsection{Class Method {\it " + functionName + "}}\n"
    else:
        functionDescription = "\\subsection{Function {\it " + functionName + "}}\n"

    functionDescription += functionComment + "\n"
    functionDescription += "\n"
    functionDescription += "\\begin{lstlisting}\n"
    functionDescription += functionSignature
    functionDescription += "\\end{lstlisting}\n"

    print(functionDescription)
    return functionDescription

################################################################################

def parseEnum(lines,startLine,endLine):
    ''' Parse a Class that implements an Enumerated type. '''
    enumDescription = [] 

    enumLine = lines[startLine]
    n1 = enumLine.find("class")
    n2 = enumLine.find("(")
    enumName = enumLine[n1+5:n2]

    enumTypes = []
    for rowNum in range(startLine+1,endLine):
        line = lines[rowNum]
        line = line.replace(" ","")
        n = line.find("=")
        if n != -1:
            enumType = line[0:n]
            enumType = enumType.replace("_","\_")
            enumTypes.append(enumType)
        else:
            break

    enumDescription.append("\\subsubsection{Enumerated Type:" + enumName + "}")
    enumDescription.append("\n")
    enumDescription.append("\\begin{itemize}")
    enumDescription.append("\n")
    for enumType in enumTypes:
        enumDescription.append("\\item{" + enumType + "}")
        enumDescription.append("\n")
    enumDescription.append("\\end{itemize}")
    enumDescription.append("\n")
    enumDescription.append("\n")

    return enumDescription

################################################################################

buildHead()
buildIntro()
buildChapter("..//finutils")
buildChapter("..//products//equities")
buildChapter("..//products//credit")
buildChapter("..//products//bonds")
buildChapter("..//products//libor")
buildChapter("..//products//fx")
buildChapter("..//models")
buildChapter("..//portfolio")
buildChapter("..//risk")
buildChapter("..//market//curves")
buildChapter("..//tests")
buildChapter("..//docs")
buildTail()

if 1==1:
    import os
    # Do it twice for the TOC
    os.system("pdflatex " + userGuideFileName)
    os.system("pdflatex " + userGuideFileName)
    
    os.startfile(fileName+".pdf")