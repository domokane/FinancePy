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

# file_name = "FinancePyManualV_" + str(VERSION)
file_name = "FinancePyManual"
user_guide_file_name = "./" + file_name + ".tex"
head_file = "./head.tex"
tail_file = "./tail.tex"
new_head_file = "./head_" + str(VERSION) + ".tex"

shutil.copyfile(head_file, new_head_file)

with fileinput.FileInput(new_head_file, inplace=True, backup='.bak') as file:
    for line in file:
        print(line.replace("VERSION_NUMBER_TO_BE_REPLACED", VERSION), end='')

VERBOSE = False
parseDataMembers = False

###############################################################################


def parse_markdown(lines):

    if VERBOSE is True:
        print("MARKDOWN IN")
        print(lines)

    parsed_lines = ["\n"]

    bullet_list_active = False
    numbered_list_active = False

    for line in lines:

        line_found = False

        line = sub("_", "\\_", line)

        if line[0] == "*":
            line_found = True
            if bullet_list_active is False:
                parsed_lines.append("\\begin{itemize}")
                parsed_lines.append("\n")
                bullet_list_active = True
                parsed_lines.append("\\item{" + line[1:] + "}")
                parsed_lines.append("\n")
            else:
                parsed_lines.append("\\item{" + line[1:] + "}")
                parsed_lines.append("\n")
        elif bullet_list_active is True:
            bullet_list_active = False
            parsed_lines.append("\\end{itemize}")
            parsed_lines.append("\n")

        if line[0].isdigit() is True:
            line_found = True
            col = line.find(".") + 1
            if numbered_list_active is False:
                parsed_lines.append("\\begin{enumerate}")
                parsed_lines.append("\n")
                numbered_list_active = True
                parsed_lines.append("\\item{" + line[col:] + "}")
                parsed_lines.append("\n")
            else:
                parsed_lines.append("\\item{" + line[col:] + "}")
                parsed_lines.append("\n")
        elif numbered_list_active is True:
            numbered_list_active = False
            parsed_lines.append("\\end{enumerate}")
            parsed_lines.append("\n")

        if line.find("###") > -1:
            line_found = True
            content = line.replace("#", "")
            content = content.replace("\n", "")
            content = content.strip()
            parsed_lines.append("\n")
            parsed_lines.append("\\subsubsection*{" + content + "}")
            parsed_lines.append("\n")
        elif line.find("##") > -1:
            line_found = True
            content = line.replace("#", "")
            content = content.replace("\n", "")
            content = content.strip()
            parsed_lines.append("\n")
            parsed_lines.append("\\subsection*{" + content + "}")
            parsed_lines.append("\n")
        elif line.find("#") > -1:
            line_found = True
            content = line.replace("#", "")
            content = content.replace("\n", "")
            content = content.strip()
            parsed_lines.append("\n")
            parsed_lines.append("\\section*{" + content + "}")
            parsed_lines.append("\n")

        if line_found is False:
            parsed_lines.append(line)

    if bullet_list_active is True:
        bullet_list_active = False
        parsed_lines.append("\\end{itemize}")
        parsed_lines.append("\n")

    if numbered_list_active is True:
        numbered_list_active = False
        parsed_lines.append("\\end{enumerate}")
        parsed_lines.append("\n")

    if VERBOSE:
        print("MARKDOWN OUT")
        print(parsed_lines)

    return parsed_lines

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

    f = open(new_head_file, 'r')
    lines = f.readlines()
    f.close()

    f = open(user_guide_file_name, 'w', encoding="utf-8")
    f.writelines(lines)
    f.close()

##########################################################################


def build_tail():
    """ Add on end latex to close document. """
    f = open(tail_file, 'r')
    lines = f.readlines()
    f.close()

    f = open(user_guide_file_name, 'a', encoding="utf-8")
    f.writelines(lines)
    f.close()

##########################################################################


def build_intro(intro_file):
    """ Add on end latex to close document. """
    print("Building introduction from file:", intro_file)

    f = open(intro_file, 'r')
    lines = f.readlines()
    f.close()

    parsed_lines = parse_markdown(lines)

    f = open(user_guide_file_name, 'a', encoding="utf-8")

    f.write("\\chapter{Introduction to FinancePy}")
    f.writelines(parsed_lines)
    f.close()

##########################################################################


def build_chapter(folder_name):
    """ Parse a folder by loading up all of the modules in that folder that
    start with the three letters - Fin. """

    print("Building chapter in folder:", folder_name)

    read_me_file = folder_name + "//" + "README.md"
    f = open(read_me_file, 'r', encoding="utf-8")
    read_me_lines = f.readlines()
    f.close()

    chapter_name = folder_name.replace("//", ".")
    chapter_name = chapter_name.replace("...", "")

    new_lines = []
    new_lines.append("\n")
    new_lines.append("\\chapter{" + chapter_name + "}")
    new_lines.append("\n")
#    new_lines.append("\\section{Introduction}")
#    new_lines.append("\n")

    f = open(user_guide_file_name, 'a', encoding="utf-8")
    f.writelines(new_lines)
    f.close()

    # validIntro = False
    # for line in read_me_lines:

    #     if line[0:3] == "Fin":
    #         validIntro = True

    # if validIntro is True:
    #     for line in read_me_lines[0:1]:
    #         line = line.replace("#","")
    #         line = line.replace("$","\\$")
    #         new_lines.append(line)

    #     new_lines.append("\\begin{itemize}\n")

    #     for line in read_me_lines[1:]:
    #         line = line.replace("#","")
    #         line = line.replace("$","\\$")
    #         if line[0:3] == "Fin":
    #             new_lines.append("\\item{" + line + "}\n")

    #     new_lines.append("\\end{itemize}")

    #     new_lines.append("\n")

    read_me_lines = parse_markdown(read_me_lines)

    f = open(user_guide_file_name, 'a', encoding="utf-8")
    f.writelines(read_me_lines)
    f.close()

    modules = glob.glob(folder_name + "//*.py")

    for module in modules:
        module_name = module.split("\\")[-1]
        escaped_module_name = sub("_", "\\_", module_name[0:-3])
        f = open(user_guide_file_name, 'a', encoding="utf-8")
        f.write("\\newpage\n")
        f.write("\\section{" + escaped_module_name + "}\n")
        f.write("\n")
        f.close()
        parse_module(module)

    modules = glob.glob(folder_name + "//TestFin*.py")

    for module in modules:
        module_name = module.split("\\")[-1]
        escaped_module_name = sub("_", "\\_", module_name[0:-3])
        f = open(user_guide_file_name, 'a', encoding="utf-8")
        f.write("\n")
        f.write("\\newpage\n")
        f.write("\\section{" + escaped_module_name + "}\n")
        f.write("\n")
        f.close()
        parse_module(module)

##########################################################################


def parse_module(module_name):
    """ Parse a module looking for classes, functions and classes for
    enumerated types. Functions inside classes are parsed inside the class. """
    print(module_name)
    f = open(module_name, 'r', encoding="utf-8")
    lines = f.readlines()
    f.close()

    lines = [sub(r"\\", r"\\\\", line) for line in lines]
    lines = [sub("_", "\\_", line) for line in lines]

    num_enums = 0
    num_classes = 0
    num_functions = 0

    start_enum_lines = []
    start_class_lines = []
    start_function_lines = []

    # Module level classes and functions
    numRows = len(lines)

    for row_num in range(0, numRows):

        line = lines[row_num]

        if line.find("class", 0, 5) != -1 and line.find("Enum") != -1:
            num_enums += 1
            start_enum_lines.append(row_num)

        if line.find("class", 0, 5) != -1 and line.find("Enum") == -1:
            num_classes += 1
            start_class_lines.append(row_num)

        if line.find("def", 0, 4) != -1:
            num_functions += 1
            start_function_lines.append(row_num)

    start_enum_lines.append(numRows)
    start_class_lines.append(numRows)
    start_function_lines.append(numRows)

    # print("start_class_lines", start_class_lines)

    f = open(user_guide_file_name, 'a', encoding="utf-8")

    for c in range(0, num_enums):
        new_lines = parse_enum(lines, start_enum_lines[c],
                               start_enum_lines[c + 1])

        for newLine in new_lines:
            f.writelines(newLine)

    for c in range(0, num_classes):
        new_lines = parse_class(lines,
                                start_class_lines[c],
                                start_class_lines[c + 1])

        for newLine in new_lines:
            f.writelines(newLine)

    for c in range(0, num_functions):
        new_lines = parse_function(
            lines, start_function_lines[c], start_function_lines[c + 1])

        for newLine in new_lines:
            f.writelines(newLine)

        f.write("\n")

    f.close()

##########################################################################


def parse_class(lines, start_line, end_line):
    """ Parse a Python class consisting of data members and functions. """

    n1 = lines[start_line].find("class")

    new_lines = []

    class_name = lines[start_line].split(" ")[1]
    class_name = class_name.replace(":", "")
    class_name = class_name.replace("\n", "")

#    print(class_name, start_line, end_line)

    new_lines.append("\\subsection*{Class: " + class_name + "}")
    new_lines.append("\n")

    ##################################################
    # COMMENT AFTER CLASS BUT BEFORE FIRST FUNCTION
    ##################################################

    start_comment_row = start_line
    end_comment_row = end_line
    comment_end_line = end_line

    for row_num in range(start_line, end_line):
        line = lines[row_num]
        if line.find(" def ") > 0:
            comment_end_line = row_num
            break

    class_comment = ""
    start_comment = False
    end_comment = False

    for row_num in range(start_line, comment_end_line):
        line = lines[row_num]
        if line.find('"""') > 0:
            start_comment_row = row_num
            start_comment = True
            start_line = row_num + 1
            break

    for row_num in range(start_line, comment_end_line):
        line = lines[row_num]
        if line.find('"""') > 0:
            end_comment_row = row_num
            end_comment = True
            break

    if start_comment is False and end_comment is False:
        # assume it's a one-line comment
        end_comment_row = start_comment_row
        end_comment = True

    if end_comment:
        #  print(start_comment_row, end_comment_row)
        for row_num in range(start_comment_row, end_comment_row + 1):
            line = lines[row_num]

            line = line.replace('"""', "")
            line = line.replace("'", "")
            line = line.replace("\n", "\n")
            line = line.replace("#", r"\#")
            line = line.lstrip()

            if len(line) == 0:
                line = "\n"

            class_comment += line + " "

    new_lines.append(class_comment)
    new_lines.append("\n")
    new_lines.append("\n")

    ##################################################
    # Now get the data members

    if parseDataMembers:
        new_lines.append("\\subsubsection*{Data Members}\n")

        data_members = []

        for row_num in range(start_line, end_line):
            row = lines[row_num]
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
                data_members = add_to_list(data_members, dataMember)

        if len(data_members) > 0:
            new_lines.append("\\begin{itemize}\n")

            for dataMember in data_members:
                new_lines.append("\\item{" + dataMember + "}\n")
            new_lines.append("\\end{itemize}")
            new_lines.append("\n")
            new_lines.append("\n")
        else:
            new_lines.append("No data members found.")
            new_lines.append("\n")
            new_lines.append("\n")

        new_lines.append("\\subsection*{Functions}\n")
        new_lines.append("\n")

    # Now get the functions
    num_class_functions = 0
    start_class_function_lines = []

#    print(start_line, end_line)
    for row_num in range(start_line, end_line):

        line = lines[row_num]

        if line.find(" def ") != -1:
            num_class_functions += 1
            start_class_function_lines.append(row_num)

    start_class_function_lines.append(end_line)

    # Remove inheritance name from class_name
    end_class_name = class_name.find("(")
    if end_class_name != -1:
        class_name = class_name[:end_class_name]

    for c in range(0, num_class_functions):
        new_lines += parse_function(lines,
                                    start_class_function_lines[c],
                                    start_class_function_lines[c + 1],
                                    class_name)
        new_lines += "\n"

    return new_lines

##########################################################################


def parse_function(lines, start_line, end_line, class_name=""):
    """ Given a set of lines and a start line I extract the function definition
    and any comment that goes below.
    TODO: Add parsing of function arguments and any comments."""

    function_line = lines[start_line]
    leftCol = function_line.find("def ")
    indent = leftCol + 4

    # Do not include a commented out function
    hashCol = function_line.find("#")
    if hashCol < leftCol and hashCol != -1:
        return ""

    n2 = function_line.find("(")
    function_name = function_line[leftCol + 4:n2]

    # If the function name starts with a _ and is not init then ignore it
    if function_name[0] == "_" and function_name != "__init__":
        return ""

    # Functions beginning with underscores ('_') are not to be parsed
    isPrivate = function_line.find("def _") != -1
    if isPrivate and function_name != r"\_\_init\_\_":
        return ""

    # Ensure function stops before any class
    for row_num in range(start_line + 1, end_line):
        line = lines[row_num]
        if line.find("class ") >= 0:
            end_line = row_num  # update start line to after function signature
            break

    # Ensure function stops before any other function
    for row_num in range(start_line + 1, end_line):
        line = lines[row_num]
        if line.find("def ") >= 0:
            end_line = row_num  # update start line to after function signature
            break

    function_signature = ""
    for row_num in range(start_line, end_line):
        line = lines[row_num][indent:]
        function_signature += str(line)
        if line.find("):") >= 0:
            start_line = row_num  # update start line to after function signature
            break

    # Replace `__init__` with class_name and remove `self` from signatures
    if class_name != "":
        # Replace '__init__' with the function's class name
        if function_name == "\\_\\_init\\_\\_":
            function_name = class_name
            function_signature = function_signature.replace("\\_\\_init\\_\\_",
                                                          class_name)

            function_signature = function_signature.replace("def ", "")

            missingSpaces = len(class_name) - len("__init__")
            if missingSpaces >= 0:
                function_signature = function_signature.replace(
                    "\n ", "\n " + " " * (missingSpaces))
            else:
                function_signature = function_signature.replace(
                    "\n" + " " * (-missingSpaces), "\n")

        # Remove 'self' and any whitespace following it
        function_signature = function_signature.replace("self", "")

        unchanged = ""
        while unchanged != function_signature:
            unchanged = function_signature
            function_signature = function_signature.replace("(,", "(")
            function_signature = function_signature.replace("( ", "(")
            function_signature = function_signature.replace("(\n", "(")

    function_comment = ""
    start_comment_row = start_line+1
    end_comment_row = start_line
    end_comment = False

    for row_num in range(start_line+1, end_line):
        line = lines[row_num]

        if line.count("'''") == 1 or line.count('"""') == 1:
            if line.count("'''") == 1:
                comment_init = "'''"
            else:
                comment_init = '"""'

            start_comment_row = row_num
            for row_num in range(row_num+1, end_line):
                line = lines[row_num]
                if line.find(comment_init) > 0:
                    end_comment_row = row_num
                    end_comment = True
                    break
            break

        if line.count("'''") == 2 or line.count('"""') == 2:
            if line.count("'''") == 2:
                comment_init = "'''"
            else:
                comment_init = '"""'

            start_comment_row = row_num
            end_comment_row = row_num
            end_comment = True
            break

    if end_comment:
        #  print(start_comment_row, end_comment_row)
        for row_num in range(start_comment_row, end_comment_row + 1):
            line = lines[row_num]
            line = line.replace("$", "\\$")
            line = line.replace(comment_init, "")
            line = line.replace("\n", "\n")
            line = line.replace("#", r"\#")
            line = line.lstrip()
            # This is because we remove trailing whitespace
            function_comment += line + " "

    if function_comment == "":
        function_comment = "PLEASE ADD A FUNCTION DESCRIPTION"

    param_description = extract_params(function_signature)

    # Inside lstlisting, backslashes used for escaping are interpreted as
    # backslashes
    # However, must be after `extract_params` where escaping is required
    function_signature = function_signature.replace("\\_", "_")

    # LATEX FORMATTING
    if class_name != "":
        func_description = r"\subsubsection*{{\bf " + \
            function_name + "}}\n"
    else:
        func_description = r"\subsubsection*{{\bf " + \
            function_name + "}}\n"

    func_description += "{\\it "
    func_description += function_comment
    func_description += "}"
    func_description += "\n"
    func_description += "\\vspace{0.25cm}\n"
    func_description += "\\begin{lstlisting}\n"
    func_description += function_signature
    func_description += "\\end{lstlisting}\n"
    func_description += "\\vspace{0.25cm}\n"
    func_description += "\\noindent \n"
    func_description += "The function arguments are described in the following table.\n"
    func_description += "\\vspace{0.25cm}\n"
    func_description += param_description
    return func_description

##########################################################################


def parse_enum(lines, start_line, end_line):
    """ Parse a Class that implements an Enumerated type. """
    enum_description = []

    enumLine = lines[start_line]
    n1 = enumLine.find("class")
    n2 = enumLine.find("(")
    # len("class ") == 6
    enum_name = enumLine[n1 + 6:n2]

    enum_types = []
    for row_num in range(start_line + 1, end_line):
        line = lines[row_num]
        line = line.replace(" ", "")
        n = line.find("=")
        if n != -1:
            enum_type = line[0:n]
            enum_types.append(enum_type)
        else:
            break

    enum_description.append("\\subsubsection*{Enumerated Type: " + enum_name+"}")
    enum_description.append("\n")
    enum_description.append("This enumerated type has the following values:\n")
    enum_description.append("\\begin{itemize}[nosep]")
    enum_description.append("\n")
    for enum_type in enum_types:
        enum_description.append("\\item{" + enum_type + "}")
        enum_description.append("\n")
    enum_description.append("\\end{itemize}")
    enum_description.append("\n")
    enum_description.append("\n")

    return enum_description

##########################################################################


def extract_params(function_signature):
    """ Parse a function signature into a table containing each function
    argument's name, type, description and default value"""
    # A good example to look at for testing is `BondConvertible`

    function_signature = function_signature.replace("%", "\%")

    # Remove information that isn't to do with the parameters
    stripedSignature = function_signature.split(
        "(", 1)[1].replace("):", "").strip()
    if stripedSignature == "":
        # The function has no parameters
        return ""

    param_description = "\\begin{center}\n"
    param_description += "\\begin{tabular}{ c  c  c  c }\n"
    param_description += "\\hline\n"
    param_description += "{ \\bf Argument Name} & { \\bf Type} & {\\bf Description} & {\\bf Default Value} \\\\\n"
    param_description += "\\hline\n"

    lines = stripedSignature.split("\n")
    for line in lines:
        # Find comment
        # If multiple arguments are on the same line as a comment,
        # the comment will be used for each argument on that line.
        comment_location = line.find("#")
        pComment = "-"
        if comment_location != -1:
            pComment = line[comment_location+1:].strip()
            line = line[:comment_location]

        line = line.strip()
        # Split by comma while leaving commas that are in square brackets '[]'.
        # This allows us to parse 'Union[Date, str]' for maturity_date_or_tenor
        if line.find("[") != -1 or line.find("(") != -1:
            # https://stackoverflow.com/questions/26808913/split-string-at-commas-except-when-in-bracket-environment
            params = []
            p = []
            bracket_level = 0
            for c in line + ",":
                if c == "," and bracket_level == 0:
                    params.append("".join(p))
                    p = []
                else:
                    if c == "[":
                        bracket_level += 1
                    elif c == "]":
                        bracket_level -= 1

                    if c == "(":
                        bracket_level += 1
                    elif c == ")":
                        bracket_level -= 1

                    p.append(c)
        else:
            params = line.split(",")

        for param in params:
            param = param.strip()
            if param == "":
                continue

            # Find default value
            p_default = "-"
            default_location = param.find("=")
            if default_location != -1:
                p_default = param[default_location+1:].strip()

                # Rip the type name out if it's an enumerated type
                if p_default[0:3] == "Fin":
                    dotCol = p_default.find(".")
                    p_default = p_default[dotCol+1:]

                param = param[:default_location]

            # Find type
            p_type = "-"
            typeLocation = param.find(':')
            if typeLocation != -1:
                p_type = param[typeLocation+1:].strip()
                p_type = parse_type(p_type)
                param = param[:typeLocation].strip()

            # Everything remaining must be the name
            p_name = param

            param_description += f"{p_name} & {p_type} & {pComment} & {p_default} \\\\\n"
            param_description += "\\hline\n"

    param_description += "\\end{tabular}"
    param_description += "\\end{center}\n"

    return param_description

###############################################################################


def parse_type(p_type):
    p_type = p_type.replace(" ", "")
    u = p_type.find("Union")
    b = p_type.find("(")
    if u != -1 and b == -1:
        lb = p_type.find("[")
        rb = p_type.find("]")
        cm = p_type.find(",")
        s = p_type[lb+1:cm] + " or " + p_type[cm+1:rb]
    elif u == -1 and b != -1:
        # Problem as list has a comma in it and this has already been used to
        # split the line of arguments above
        lb = p_type.find("(")
        rb = p_type.find(")")
        cm = p_type.find(",")
        s = p_type[lb+1:cm] + " or " + p_type[cm+1:rb]
    else:
        s = p_type

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

print("Latex filename:", user_guide_file_name)

if 1 == 1:
    # Do it twice for the TOC
    print("pdflatex " + user_guide_file_name)
    os.system("pdflatex " + user_guide_file_name)
    print("Doing it again for TOC")
    os.system("pdflatex " + user_guide_file_name)

    pdf_filename_1 = file_name + ".pdf"
    pdf_filename_2 = '..\\' + pdf_filename_1

    # TODO: Only works if you have financepy-examples-git
    # Maybe add `financepy-examples-git` as a submodule?

    print("Removing unneeded files.")
    os.remove(file_name + ".out")
#    os.remove(file_name + ".tex")
    os.remove(file_name + ".toc")
    os.remove(file_name + ".aux")
#    os.remove(file_name + ".fls")
    os.remove(file_name + ".fdb_latexmk")
    os.remove(file_name + ".log")
    os.remove(new_head_file)
    os.remove(new_head_file + ".bak")

    print("Moving ", pdf_filename_1, " to ", pdf_filename_2)
    shutil.move(pdf_filename_1, pdf_filename_2)
    print(pdf_filename_2)
    open_file(pdf_filename_1)
