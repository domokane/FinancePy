# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:55:21 2017

@author: Dominic
"""

import glob
import os
import sys
import subprocess

VERSION = 0.15

fileName = "FinancePyManualV_" + str(VERSION)
userGuideFileName = "./" + fileName + ".tex"
headFile = "./head.tex"
tailFile = "./tail.tex"

verbose = False

##########################################################################

def parseMarkdown(lines):

    if verbose == True:
        print("MARKDOWN IN")
        print(lines)

    parsedLines = ["\n"]

    bulletListActive = False
    numberedListActive = False

    for line in lines:

        lineFound = False

        if line[0] == "*":
            lineFound= True
            if bulletListActive is False:
                parsedLines.append("\\begin{itemize}")
                parsedLines.append("\n")
                bulletListActive = True
                parsedLines.append("\\item{" + line[1:] + "}")
                parsedLines.append("\n")
            else:
                parsedLines.append("\\item{" + line[1:] + "}")
                parsedLines.append("\n")
        elif bulletListActive is True:
            bulletListActive = False
            parsedLines.append("\\end{itemize}")
            parsedLines.append("\n")

        if line[0].isdigit() is True:
            lineFound= True
            col = line.find(".") + 1
            if numberedListActive is False:
                parsedLines.append("\\begin{enumerate}")
                parsedLines.append("\n")
                numberedListActive = True
                parsedLines.append("\\item{" + line[col:] + "}")
                parsedLines.append("\n")
            else:
                parsedLines.append("\\item{" + line[col:] + "}")
                parsedLines.append("\n")
        elif numberedListActive is True:
            numberedListActive = False
            parsedLines.append("\\end{enumerate}")
            parsedLines.append("\n")

        if line.find("###") > -1:
            lineFound= True
            content = line.replace("#","")
            content = content.replace("\n","")
            content = content.strip()
            parsedLines.append("\n")
            parsedLines.append("\\subsubsection*{" + content + "}")
            parsedLines.append("\n")
        elif line.find("##") > -1:
            lineFound= True
            content = line.replace("#","")
            content = content.replace("\n","")
            content = content.strip()
            parsedLines.append("\n")
            parsedLines.append("\\subsection*{" + content + "}")
            parsedLines.append("\n")
        elif line.find("#") > -1:
            lineFound= True
            content = line.replace("#","")
            content = content.replace("\n","")
            content = content.strip()
            parsedLines.append("\n")
            parsedLines.append("\\section*{" + content + "}")
            parsedLines.append("\n")

        if lineFound is False:
            parsedLines.append(line)

    if bulletListActive is True:
        bulletListActive = False
        parsedLines.append("\\end{itemize}")
        parsedLines.append("\n")

    if numberedListActive is True:
        numberedListActive = False
        parsedLines.append("\\end{enumerate}")
        parsedLines.append("\n")

    if verbose:
        print("MARKDOWN OUT")
        print(parsedLines)

    return parsedLines

##########################################################################

def addToList(listName, item):

    for x in listName:
        if item == x:
            return listName

    listName.append(item)
    return listName

##########################################################################

def open_file(filename):
    if sys.platform == "win32":
        os.startfile(filename)
    else:
        opener = "open" if sys.platform == "darwin" else "xdg-open"
        subprocess.call([opener, filename])

##########################################################################


def buildHead():
    ''' Start latex file with a header that sets all of the document
    properties. '''

    f = open(headFile, 'r')
    lines = f.readlines()
    f.close()

    f = open(userGuideFileName, 'w')
    f.writelines(lines)
    f.close()

##########################################################################


def buildTail():
    ''' Add on end latex to close document. '''
    f = open(tailFile, 'r')
    lines = f.readlines()
    f.close()

    f = open(userGuideFileName, 'a')
    f.writelines(lines)
    f.close()

##########################################################################


def buildIntro(introfile):
    ''' Add on end latex to close document. '''
    print("Building introduction from file:", introfile)

    f = open(introfile, 'r')
    lines = f.readlines()
    f.close()

    parsedLines = parseMarkdown(lines)

    f = open(userGuideFileName, 'a')

    f.write("\chapter{Introduction to FinancePy}")
    f.writelines(parsedLines)
    f.close()

##########################################################################


def buildChapter(folderName):
    ''' Parse a folder by loading up all of the modules in that folder that
    start with the threee letters - Fin. '''

    print("Building chapter in folder:", folderName)

    readMeFile = folderName + "//" + "README.md"
    f = open(readMeFile, 'r')
    readMeLines = f.readlines()
    f.close()

    chapterName = folderName.replace("//",".")
    chapterName = chapterName.replace("...","")

    newLines = []
    newLines.append("\n")
    newLines.append("\\chapter{" + chapterName + "}")
    newLines.append("\n")
    newLines.append("\\section{Introduction}")
    newLines.append("\n")

    f = open(userGuideFileName, 'a')
    f.writelines(newLines)
    f.close()

    # validIntro = False
    # for line in readMeLines:

    #     if line[0:3] == "Fin":
    #         validIntro = True

    # if validIntro is True:
    #     for line in readMeLines[0:1]:
    #         line = line.replace("#","")
    #         line = line.replace("$","\\$")
    #         newLines.append(line)

    #     newLines.append("\\begin{itemize}\n")

    #     for line in readMeLines[1:]:
    #         line = line.replace("#","")
    #         line = line.replace("$","\\$")
    #         if line[0:3] == "Fin":
    #             newLines.append("\\item{" + line + "}\n")

    #     newLines.append("\\end{itemize}")

    #     newLines.append("\n")

    readMeLines = parseMarkdown(readMeLines)

    f = open(userGuideFileName, 'a')
    f.writelines(readMeLines)
    f.close()

    modules = glob.glob(folderName + "//Fin*.py")

    for module in modules:
        moduleName = module.split("\\")[-1]
        f = open(userGuideFileName, 'a')
        f.write("\\newpage\n")
        f.write("\\section{" + moduleName[0:-3] + "}\n")
        f.write("\n")
        f.close()
        parseModule(module)

    modules = glob.glob(folderName + "//TestFin*.py")

    for module in modules:
        moduleName = module.split("\\")[-1]
        f = open(userGuideFileName, 'a')
        f.write("\n")
        f.write("\\newpage\n")
        f.write("\\section{" + moduleName[0:-3] + "}\n")
        f.write("\n")
        f.close()
        parseModule(module)

##########################################################################


def parseModule(moduleName):
    ''' Parse a module looking for classes, functions and classes for
    enumerated types. Functions inside classes are parsed inside the class. '''
    print(moduleName)
    f = open(moduleName, 'r')
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

    for rowNum in range(0, numRows):

        line = lines[rowNum]

        if line.find("class", 0, 5) != -1 and line.find("Enum") != -1:
            numEnums += 1
            startEnumLines.append(rowNum)

        if line.find("class", 0, 5) != -1 and line.find("Enum") == -1:
            numClasses += 1
            startClassLines.append(rowNum)

        if line.find("def", 0, 4) != -1:
            numFunctions += 1
            startFunctionLines.append(rowNum)

    startEnumLines.append(numRows)
    startClassLines.append(numRows)
    startFunctionLines.append(numRows)

    # print("startClassLines", startClassLines)

    f = open(userGuideFileName, 'a')

    for c in range(0, numEnums):
        newLines = parseEnum(lines, startEnumLines[c], startEnumLines[c + 1])

        for newLine in newLines:
            f.writelines(newLine)

    for c in range(0, numClasses):
        newLines = parseClass(lines,
                              startClassLines[c],
                              startClassLines[c + 1])

        for newLine in newLines:
            f.writelines(newLine)

    for c in range(0, numFunctions):
        newLines = parseFunction(
            lines, startFunctionLines[c], startFunctionLines[c + 1], False)

        for newLine in newLines:
            f.writelines(newLine)

        f.write("\n")

    f.close()

##########################################################################


def parseClass(lines, startLine, endLine):
    ''' Parse a Python class consisting of data members and functions. '''

    n1 = lines[startLine].find("class")

    newLines = []

    className = lines[startLine].split(" ")[1]
    className = className.replace(":", "")
    className = className.replace("\n", "")

#    print(className, startLine, endLine)

    newLines.append("\\subsection*{Class: " + className + "}")
    newLines.append("\n")

    ##################################################
    # COMMENT AFTER CLASS BUT BEFORE FIRST FUNCTION
    ##################################################

    startCommentRow = startLine
    endCommentRow = endLine
    commentEndLine = endLine

    for rowNum in range(startLine, endLine):
        line = lines[rowNum]
        if line.find(" def ") > 0:
            commentEndLine = rowNum
            break

    classComment = ""
    startComment = False
    endComment = False

    for rowNum in range(startLine, commentEndLine):
        line = lines[rowNum]
        if line.find("'''") > 0:
            startCommentRow = rowNum
            startComment = True
            startLine = rowNum + 1
            break

    for rowNum in range(startLine, commentEndLine):
        line = lines[rowNum]
        if line.find("'''") > 0:
            endCommentRow = rowNum
            endComment = True
            break

    if startComment is False and endComment is False:
        # assume it's a one-line comment
        endCommentRow = startCommentRow
        endComment = True

    if endComment:
#        print(startCommentRow, endCommentRow)
        for rowNum in range(startCommentRow, endCommentRow + 1):
            line = lines[rowNum]
            line = line.replace("'", "")
            line = line.replace("_", r"\_")
            line = line.replace("\n", "")
            line = line.replace("#", r"\#")
            line = line.lstrip()
            classComment += line + " "

    newLines.append(classComment)
    newLines.append("\n")
    newLines.append("\n")

    ##################################################
    # Now get the data members

    newLines.append("\\subsubsection*{Data Members}\n")
#    newLines.append("The class data members are:")

    dataMembers = []

    for rowNum in range(startLine, endLine):
        row = lines[rowNum]
        row = row.replace(" ", "")
        row = row.replace("_", r"\_")
        row = row.replace("!", "")
        row = row.replace("<", "")
        row = row.replace(">", "")

        n1 = row.find("self.")
        n2 = row.find("=")
        n3 = row.find("[")

        if n1 != -1 and n2 > n1 and n3 == -1:
            dataMember = row[n1:n2]
            dataMember = dataMember.strip(" ")
            dataMember = dataMember.strip(")")
            dataMember = dataMember.replace("self.","")
            dataMembers = addToList(dataMembers,dataMember)

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

    newLines.append("\\subsection*{Functions}\n")
    newLines.append("\n")

    # Now get the functions
    numClassFunctions = 0
    startClassFunctionLines = []

#    print(startLine, endLine)
    for rowNum in range(startLine, endLine):

        line = lines[rowNum]

        if line.find(" def ") != -1:
            numClassFunctions += 1
            startClassFunctionLines.append(rowNum)

    startClassFunctionLines.append(endLine)

    for c in range(0, numClassFunctions):
        newLines += parseFunction(lines,
                                  startClassFunctionLines[c],
                                  startClassFunctionLines[c + 1],
                                  True)
        newLines += "\n"

    return newLines

##########################################################################


def parseFunction(lines, startLine, endLine, classFlag):
    ''' Given a set of lines and a start line I extract the function definiton
    and any comment that goes below.
    TODO: Add parsing of function arguments and any comments.'''

    functionLine = lines[startLine]
    leftCol = functionLine.find("def ")
    n2 = functionLine.find("(")
    functionName = functionLine[leftCol + 4:n2]
    functionName = functionName.replace("_", r"\_")

    # Ensure function stops before any class
    for rowNum in range(startLine + 1, endLine):
        line = lines[rowNum]
        if line.find("class ") >= 0:
            endLine = rowNum  # update start line to after function signature
            break

    # Ensure function stops before any other function
    for rowNum in range(startLine + 1, endLine):
        line = lines[rowNum]
        if line.find("def ") >= 0:
            endLine = rowNum  # update start line to after function signature
            break

    functionSignature = ""
    for rowNum in range(startLine, endLine):
        line = lines[rowNum]
        functionSignature += str(line)
        if line.find(":") >= 0:
            startLine = rowNum  # update start line to after function signature
            break

    functionComment = ""
    startCommentRow = startLine+1
    endCommentRow = startLine
    startComment = False
    endComment = False

    for rowNum in range(startLine, endLine):
        line = lines[rowNum]

        if line.count("'''") == 1:
            startCommentRow = rowNum
            startComment = True
            startLine = rowNum + 1
            break

        if line.count("'''") == 2:
            startCommentRow = rowNum
            startComment = True
            startLine = rowNum
            break

    for rowNum in range(startLine, endLine):
        line = lines[rowNum]
        if line.find("'''") > 0:
            endCommentRow = rowNum
            endComment = True
            break

    if startComment is False and endComment is False:
        # assume it's a one-line comment
        endCommentRow = startCommentRow
        endComment = True

    if endComment:
#        print(startCommentRow, endCommentRow)
        for rowNum in range(startCommentRow, endCommentRow + 1):
            line = lines[rowNum]
            line = line.replace("_", r"\_")
            line = line.replace("'", "")
            line = line.replace("\n", "")
            line = line.replace("#", r"\#")
            line = line.lstrip()
            # This is because we remove trailing whitespace
            functionComment += line + " "


    if functionComment == " ":
        functionComment = "PLEASE ADD A FUNCTION DESCRIPTION"

    # LATEX FORMATTING
    if classFlag:
        functionDescription = r"\subsubsection*{{\bf " + \
            functionName + "}}\n"
    else:
        functionDescription = r"\subsubsection*{{\bf " + \
            functionName + "}}\n"

    functionDescription += functionComment + "\n"
    functionDescription += "\n"
    functionDescription += "\\begin{lstlisting}\n"
    functionDescription += functionSignature
    functionDescription += "\\end{lstlisting}\n"

    return functionDescription

##########################################################################


def parseEnum(lines, startLine, endLine):
    ''' Parse a Class that implements an Enumerated type. '''
    enumDescription = []

    enumLine = lines[startLine]
    n1 = enumLine.find("class")
    n2 = enumLine.find("(")
    enumName = enumLine[n1 + 5:n2]

    enumTypes = []
    for rowNum in range(startLine + 1, endLine):
        line = lines[rowNum]
        line = line.replace(" ", "")
        n = line.find("=")
        if n != -1:
            enumType = line[0:n]
            enumType = enumType.replace("_", r"\_")
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

##########################################################################


buildHead()
buildIntro("..//README.md")

if 1==1:
    buildChapter("..//financepy//finutils")
    buildChapter("..//financepy//products//equity")
    buildChapter("..//financepy//products//credit")
    buildChapter("..//financepy//products//bonds")
    buildChapter("..//financepy//products//libor")
    buildChapter("..//financepy//products//fx")
    buildChapter("..//financepy//models")
    buildChapter("..//financepy//portfolio")
    buildChapter("..//financepy//risk")
    buildChapter("..//financepy//market//curves")
#    buildChapter(".//financepy//tests")
#    buildChapter(".//financepy//docs")

buildTail()

if 1 == 1:
    # Do it twice for the TOC
    os.system("pdflatex " + userGuideFileName)
    open_file(fileName + ".pdf")
