###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


import glob
import os
import sys
import subprocess
import shutil
import fileinput
from re import sub

with open("../version.py", "r") as fh:
    version_number = fh.read()
    start = version_number.find("\"")
    end = version_number[start+1:].find("\"")
    VERSION = str(version_number[start+1:start+end+1])
    VERSION = VERSION.replace('\n', '')

# fileName = "FinancePyManualV_" + str(VERSION)
fileName = "FinancePyManual"
userGuideFileName = "./" + fileName + ".tex"
headFile = "./head.tex"
tailFile = "./tail.tex"
newHeadFile = "./head_" + str(VERSION) + ".tex"

shutil.copyfile(headFile, newHeadFile)

with fileinput.FileInput(newHeadFile, inplace=True, backup='.bak') as file:
    for line in file:
        print(line.replace("VERSION_NUMBER_TO_BE_REPLACED", VERSION), end='')

verbose = False
parseDataMembers = False

###############################################################################


def parse_markdown(lines):

    if verbose is True:
        print("MARKDOWN IN")
        print(lines)

    parsedLines = ["\n"]

    bulletListActive = False
    numberedListActive = False

    for line in lines:

        lineFound = False

        line = sub("_", "\\_", line)

        if line[0] == "*":
            lineFound = True
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
            lineFound = True
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
            lineFound = True
            content = line.replace("#", "")
            content = content.replace("\n", "")
            content = content.strip()
            parsedLines.append("\n")
            parsedLines.append("\\subsubsection*{" + content + "}")
            parsedLines.append("\n")
        elif line.find("##") > -1:
            lineFound = True
            content = line.replace("#", "")
            content = content.replace("\n", "")
            content = content.strip()
            parsedLines.append("\n")
            parsedLines.append("\\subsection*{" + content + "}")
            parsedLines.append("\n")
        elif line.find("#") > -1:
            lineFound = True
            content = line.replace("#", "")
            content = content.replace("\n", "")
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


def add_to_list(listName, item):

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


def build_head():
    """ Start latex file with a header that sets all of the document
    properties. """

    f = open(newHeadFile, 'r')
    lines = f.readlines()
    f.close()

    f = open(userGuideFileName, 'w', encoding="utf-8")
    f.writelines(lines)
    f.close()

##########################################################################


def build_tail():
    """ Add on end latex to close document. """
    f = open(tailFile, 'r')
    lines = f.readlines()
    f.close()

    f = open(userGuideFileName, 'a', encoding="utf-8")
    f.writelines(lines)
    f.close()

##########################################################################


def build_intro(introfile):
    """ Add on end latex to close document. """
    print("Building introduction from file:", introfile)

    f = open(introfile, 'r')
    lines = f.readlines()
    f.close()

    parsedLines = parse_markdown(lines)

    f = open(userGuideFileName, 'a', encoding="utf-8")

    f.write("\\chapter{Introduction to FinancePy}")
    f.writelines(parsedLines)
    f.close()

##########################################################################


def build_chapter(folderName):
    """ Parse a folder by loading up all of the modules in that folder that
    start with the three letters - Fin. """

    print("Building chapter in folder:", folderName)

    readMeFile = folderName + "//" + "README.md"
    f = open(readMeFile, 'r', encoding="utf-8")
    readMeLines = f.readlines()
    f.close()

    chapterName = folderName.replace("//", ".")
    chapterName = chapterName.replace("...", "")

    newLines = []
    newLines.append("\n")
    newLines.append("\\chapter{" + chapterName + "}")
    newLines.append("\n")
#    newLines.append("\\section{Introduction}")
#    newLines.append("\n")

    f = open(userGuideFileName, 'a', encoding="utf-8")
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

    readMeLines = parse_markdown(readMeLines)

    f = open(userGuideFileName, 'a', encoding="utf-8")
    f.writelines(readMeLines)
    f.close()

    modules = glob.glob(folderName + "//*.py")

    for module in modules:
        moduleName = module.split("\\")[-1]
        escapedModuleName = sub("_", "\\_", moduleName[0:-3])
        f = open(userGuideFileName, 'a', encoding="utf-8")
        f.write("\\newpage\n")
        f.write("\\section{" + escapedModuleName + "}\n")
        f.write("\n")
        f.close()
        parse_module(module)

    modules = glob.glob(folderName + "//TestFin*.py")

    for module in modules:
        moduleName = module.split("\\")[-1]
        escapedModuleName = sub("_", "\\_", moduleName[0:-3])
        f = open(userGuideFileName, 'a', encoding="utf-8")
        f.write("\n")
        f.write("\\newpage\n")
        f.write("\\section{" + escapedModuleName + "}\n")
        f.write("\n")
        f.close()
        parse_module(module)

##########################################################################


def parse_module(moduleName):
    """ Parse a module looking for classes, functions and classes for
    enumerated types. Functions inside classes are parsed inside the class. """
    print(moduleName)
    f = open(moduleName, 'r', encoding="utf-8")
    lines = f.readlines()
    f.close()

    lines = [sub(r"\\", r"\\\\", line) for line in lines]
    lines = [sub("_", "\\_", line) for line in lines]

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

    f = open(userGuideFileName, 'a', encoding="utf-8")

    for c in range(0, numEnums):
        newLines = parse_enum(lines, startEnumLines[c], startEnumLines[c + 1])

        for newLine in newLines:
            f.writelines(newLine)

    for c in range(0, numClasses):
        newLines = parse_class(lines,
                              startClassLines[c],
                              startClassLines[c + 1])

        for newLine in newLines:
            f.writelines(newLine)

    for c in range(0, numFunctions):
        newLines = parse_function(
            lines, startFunctionLines[c], startFunctionLines[c + 1])

        for newLine in newLines:
            f.writelines(newLine)

        f.write("\n")

    f.close()

##########################################################################


def parse_class(lines, startLine, endLine):
    """ Parse a Python class consisting of data members and functions. """

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
        if line.find('"""') > 0:
            startCommentRow = rowNum
            startComment = True
            startLine = rowNum + 1
            break

    for rowNum in range(startLine, commentEndLine):
        line = lines[rowNum]
        if line.find('"""') > 0:
            endCommentRow = rowNum
            endComment = True
            break

    if startComment is False and endComment is False:
        # assume it's a one-line comment
        endCommentRow = startCommentRow
        endComment = True

    if endComment:
        #  print(startCommentRow, endCommentRow)
        for rowNum in range(startCommentRow, endCommentRow + 1):
            line = lines[rowNum]

            line = line.replace('"""', "")
            line = line.replace("'", "")
            line = line.replace("\n", "\n")
            line = line.replace("#", r"\#")
            line = line.lstrip()

            if len(line) == 0:
                line = "\n"

            classComment += line + " "

    newLines.append(classComment)
    newLines.append("\n")
    newLines.append("\n")

    ##################################################
    # Now get the data members

    if parseDataMembers:
        newLines.append("\\subsubsection*{Data Members}\n")

        dataMembers = []

        for rowNum in range(startLine, endLine):
            row = lines[rowNum]
            row = row.replace(" ", "")
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
                dataMember = dataMember.replace("self.", "")
                dataMembers = addToList(dataMembers, dataMember)

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

    # Remove inheritance name from className
    endClassName = className.find("(")
    if (endClassName != -1):
        className = className[:endClassName]

    for c in range(0, numClassFunctions):
        newLines += parse_function(lines,
                                  startClassFunctionLines[c],
                                  startClassFunctionLines[c + 1],
                                  className)
        newLines += "\n"

    return newLines

##########################################################################


def parse_function(lines, startLine, endLine, className=""):
    """ Given a set of lines and a start line I extract the function definition
    and any comment that goes below.
    TODO: Add parsing of function arguments and any comments."""

    functionLine = lines[startLine]
    leftCol = functionLine.find("def ")
    indent = leftCol + 4

    # Do not include a commented out function
    hashCol = functionLine.find("#")
    if hashCol < leftCol and hashCol != -1:
        return ""

    n2 = functionLine.find("(")
    functionName = functionLine[leftCol + 4:n2]

    # If the function name starts with a _ and is not init then ignore it
    if functionName[0] == "_" and functionName != "__init__":
        return ""

    # Functions beginning with underscores ('_') are not to be parsed
    isPrivate = (functionLine.find("def _") != -1)
    if isPrivate and functionName != r"\_\_init\_\_":
        return ""

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
        line = lines[rowNum][indent:]
        functionSignature += str(line)
        if line.find("):") >= 0:
            startLine = rowNum  # update start line to after function signature
            break

    # Replace `__init__` with className and remove `self` from signatures
    if className != "":
        # Replace '__init__' with the function's class name
        if functionName == "\\_\\_init\\_\\_":
            functionName = className
            functionSignature = functionSignature.replace("\\_\\_init\\_\\_",
                                                          className)

            functionSignature = functionSignature.replace("def ", "")

            missingSpaces = len(className) - len("__init__")
            if (missingSpaces >= 0):
                functionSignature = functionSignature.replace(
                    "\n ", "\n " + " " * (missingSpaces))
            else:
                functionSignature = functionSignature.replace(
                    "\n" + " " * (-missingSpaces), "\n")

        # Remove 'self' and any whitespace following it
        functionSignature = functionSignature.replace("self", "")

        unchanged = ""
        while unchanged != functionSignature:
            unchanged = functionSignature
            functionSignature = functionSignature.replace("(,", "(")
            functionSignature = functionSignature.replace("( ", "(")
            functionSignature = functionSignature.replace("(\n", "(")

    functionComment = ""
    startCommentRow = startLine+1
    endCommentRow = startLine
    endComment = False

    for rowNum in range(startLine+1, endLine):
        line = lines[rowNum]

        if line.count("'''") == 1 or line.count('"""') == 1:
            if line.count("'''") == 1:
                commentInit = "'''"
            else:
                commentInit = '"""'

            startCommentRow = rowNum
            for rowNum in range(rowNum+1, endLine):
                line = lines[rowNum]
                if line.find(commentInit) > 0:
                    endCommentRow = rowNum
                    endComment = True
                    break
            break

        if line.count("'''") == 2 or line.count('"""') == 2:
            if line.count("'''") == 2:
                commentInit = "'''"
            else:
                commentInit = '"""'

            startCommentRow = rowNum
            endCommentRow = rowNum
            endComment = True
            break

    if endComment:
        #  print(startCommentRow, endCommentRow)
        for rowNum in range(startCommentRow, endCommentRow + 1):
            line = lines[rowNum]
            line = line.replace("$", "\\$")
            line = line.replace(commentInit, "")
            line = line.replace("\n", "\n")
            line = line.replace("#", r"\#")
            line = line.lstrip()
            # This is because we remove trailing whitespace
            functionComment += line + " "

    if functionComment == "":
        functionComment = "PLEASE ADD A FUNCTION DESCRIPTION"

    paramDescription = extract_params(functionSignature)

    # Inside lstlisting, backslashes used for escaping are interpreted as backslashes
    # However, must be after `extract_params` where escaping is required
    functionSignature = functionSignature.replace("\\_", "_")

    # LATEX FORMATTING
    if className != "":
        functionDescription = r"\subsubsection*{{\bf " + \
            functionName + "}}\n"
    else:
        functionDescription = r"\subsubsection*{{\bf " + \
            functionName + "}}\n"

    functionDescription += "{\\it "
    functionDescription += functionComment
    functionDescription += "}"
    functionDescription += "\n"
    functionDescription += "\\vspace{0.25cm}\n"
    functionDescription += "\\begin{lstlisting}\n"
    functionDescription += functionSignature
    functionDescription += "\\end{lstlisting}\n"
    functionDescription += "\\vspace{0.25cm}\n"
    functionDescription += "\\noindent \n"
    functionDescription += "The function arguments are described in the following table.\n"
    functionDescription += "\\vspace{0.25cm}\n"
    functionDescription += paramDescription
    return functionDescription

##########################################################################


def parse_enum(lines, startLine, endLine):
    """ Parse a Class that implements an Enumerated type. """
    enumDescription = []

    enumLine = lines[startLine]
    n1 = enumLine.find("class")
    n2 = enumLine.find("(")
    # len("class ") == 6
    enumName = enumLine[n1 + 6:n2]

    enumTypes = []
    for rowNum in range(startLine + 1, endLine):
        line = lines[rowNum]
        line = line.replace(" ", "")
        n = line.find("=")
        if n != -1:
            enumType = line[0:n]
            enumTypes.append(enumType)
        else:
            break

    enumDescription.append("\\subsubsection*{Enumerated Type: " + enumName+"}")
    enumDescription.append("\n")
    enumDescription.append("This enumerated type has the following values:\n")
    enumDescription.append("\\begin{itemize}[nosep]")
    enumDescription.append("\n")
    for enumType in enumTypes:
        enumDescription.append("\\item{" + enumType + "}")
        enumDescription.append("\n")
    enumDescription.append("\\end{itemize}")
    enumDescription.append("\n")
    enumDescription.append("\n")

    return enumDescription

##########################################################################


def extract_params(functionSignature):
    """ Parse a function signature into a table containing each function
    argument's name, type, description and default value"""
    # A good example to look at for testing is `BondConvertible`

    functionSignature = functionSignature.replace("%", "\%")

    # Remove information that isn't to do with the parameters
    stripedSignature = functionSignature.split(
        "(", 1)[1].replace("):", "").strip()
    if stripedSignature == "":
        # The function has no parameters
        return ""

    paramDescription = "\\begin{center}\n"
    paramDescription += "\\begin{tabular}{ c  c  c  c }\n"
    paramDescription += "\\hline\n"
    paramDescription += "{ \\bf Argument Name} & { \\bf Type} & {\\bf Description} & {\\bf Default Value} \\\\\n"
    paramDescription += "\\hline\n"

    lines = stripedSignature.split("\n")
    for line in lines:
        # Find comment
        # If multiple arguments are on the same line as a comment,
        # the comment will be used for each argument on that line.
        commentLocation = line.find("#")
        pComment = "-"
        if commentLocation != -1:
            pComment = line[commentLocation+1:].strip()
            line = line[:commentLocation]

        line = line.strip()
        # Split by comma while leaving commas that are in square brackets '[]'.
        # This allows us to parse 'Union[Date, str]' for maturity_date_or_tenor.
        if line.find("[") != -1 or line.find("(") != -1:
            # https://stackoverflow.com/questions/26808913/split-string-at-commas-except-when-in-bracket-environment
            params = []
            p = []
            bracketLevel = 0
            for c in line + ",":
                if c == "," and bracketLevel == 0:
                    params.append("".join(p))
                    p = []
                else:
                    if c == "[":
                        bracketLevel += 1
                    elif c == "]":
                        bracketLevel -= 1

                    if c == "(":
                        bracketLevel += 1
                    elif c == ")":
                        bracketLevel -= 1

                    p.append(c)
        else:
            params = line.split(",")

        for param in params:
            param = param.strip()
            if param == "":
                continue

            # Find default value
            pDefault = "-"
            defaultLocation = param.find("=")
            if defaultLocation != -1:
                pDefault = param[defaultLocation+1:].strip()

                # Rip the type name out if it's an enumerated type
                if pDefault[0:3] == "Fin":
                    dotCol = pDefault.find(".")
                    pDefault = pDefault[dotCol+1:]

                param = param[:defaultLocation]

            # Find type
            pType = "-"
            typeLocation = param.find(':')
            if typeLocation != -1:
                pType = param[typeLocation+1:].strip()
                pType = parse_type(pType)
                param = param[:typeLocation].strip()

            # Everything remaining must be the name
            pName = param

            paramDescription += f"{pName} & {pType} & {pComment} & {pDefault} \\\\\n"
            paramDescription += "\\hline\n"

    paramDescription += "\\end{tabular}"
    paramDescription += "\\end{center}\n"

    return paramDescription

###############################################################################


def parse_type(pType):
    pType = pType.replace(" ", "")
    u = pType.find("Union")
    b = pType.find("(")
    if u != -1 and b == -1:
        lb = pType.find("[")
        rb = pType.find("]")
        cm = pType.find(",")
        s = pType[lb+1:cm] + " or " + pType[cm+1:rb]
    elif u == -1 and b != -1:
        # Problem as list has a comma in it and this has already been used to
        # split the line of arguments above
        lb = pType.find("(")
        rb = pType.find(")")
        cm = pType.find(",")
        s = pType[lb+1:cm] + " or " + pType[cm+1:rb]
    else:
        s = pType

    return s

###############################################################################


build_head()
build_intro("..//README.md")

if 1 == 1:
    build_chapter("..//financepy//utils")
    build_chapter("..//financepy//market//curves")
    build_chapter("..//financepy//market//volatility")
    build_chapter("..//financepy//products//equity")
    build_chapter("..//financepy//products//credit")
    build_chapter("..//financepy//products//bonds")
    build_chapter("..//financepy//products//rates")
    build_chapter("..//financepy//products//fx")
    build_chapter("..//financepy//models")

#    build_chapter("..//financepy//portfolio")
#    build_chapter("..//financepy//risk")
#    build_chapter(".//financepy//tests")
#    build_chapter(".//financepy//docs")
    pass

build_tail()

print("Latex filename:", userGuideFileName)

if 1 == 1:
    # Do it twice for the TOC
    print("pdflatex " + userGuideFileName)
    os.system("pdflatex " + userGuideFileName)
    print("Doing it again for TOC")
    os.system("pdflatex " + userGuideFileName)

    pdfFileName1 = fileName + ".pdf"
    pdfFileName2 = '..\\' + pdfFileName1

    # TODO: Only works if you have financepy-examples-git
    # Maybe add `financepy-examples-git` as a submodule?

    print("Removing unneeded files.")
    os.remove(fileName + ".out")
#    os.remove(fileName + ".tex")
    os.remove(fileName + ".toc")
    os.remove(fileName + ".aux")
    os.remove(fileName + ".log")
    os.remove(newHeadFile)
    os.remove(newHeadFile + ".bak")

#    print("Moving ", pdfFileName1, " to ", pdfFileName2)
#    shutil.move(pdfFileName1, pdfFileName2)
#    print(pdfFileName2)
    open_file(pdfFileName1)
