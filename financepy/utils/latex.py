###############################################################################
# Some code to convert output to latex tables
###############################################################################

def convertToLatexTable(txt, sep=" ", header_list=[]):

    num_header_cols = len(header_list)

    txt = txt.replace("_", "\_")
    txt = txt.replace("\t", "")
    txt = txt.replace("      ", " ")
    txt = txt.replace("     ", " ")
    txt = txt.replace("    ", " ")
    txt = txt.replace("   ", " ")
    txt = txt.replace("  ", " ")

    col_str = "{"
    for i in range(0, num_header_cols):
        col_str += "r"
    col_str += "}"
    
    table_str = "\\begin{table}[htbp]\n"    
    table_str += "\\begin{center}\n"       
    table_str += "\\begin{tabular} " + col_str + "\n"    

    if header_list != []:
        header_str = header_list[0]
        for i in range(1, num_header_cols):
            header_str += " & " + header_list[i]
        header_str += "\\\ \n"    
        table_str += header_str

    rows = txt.split("\n")

    num_rows= len(rows)

    for i in range(0, num_rows):
        
        row_str = rows[i]
        cols = row_str.split(sep)
        
        num_cols = len(cols)
        
        if num_header_cols > 0:
            if num_cols != num_header_cols:
                print("Num row cols " + str(num_cols) + \
                      " is not the same as the number of header cols " \
                          + str(num_header_cols))
                return ""

        col_str = cols[0]
        for i in range(1, num_cols):
            col_str += " & " + cols[i]
        
        table_str += col_str + "\\\ \n"

    table_str += "\end{tabular}\n"
    table_str += "\end{center}\n"    
    table_str += "\end{table}\n"    

    return table_str

###############################################################################

txt = "NOV-15-2017      11875.00\n\
MAY-15-2018      11875.00\n\
NOV-15-2018      11875.00\n\
MAY-15-2019      11875.00\n\
NOV-15-2019      11875.00\n\
MAY-15-2020      11875.00\n\
NOV-15-2020      11875.00\n\
MAY-15-2021      11875.00\n\
NOV-15-2021      11875.00\n\
MAY-15-2022      11875.00\n\
NOV-15-2022      11875.00\n\
MAY-15-2023      11875.00\n\
NOV-15-2023      11875.00\n\
MAY-15-2024      11875.00\n\
NOV-15-2024      11875.00\n\
MAY-15-2025      11875.00\n\
NOV-15-2025      11875.00\n\
MAY-15-2026      11875.00\n\
NOV-15-2026      11875.00\n\
MAY-15-2027    1011875.00"

latex_str = convertToLatexTable(txt, " ", ["Dates","Flows"])
print(latex_str)

 
txt="OBJECT TYPE: Bond\n\
ISSUE DATE: MAY-15-2010\n\
MATURITY DATE: MAY-15-2027\n\
COUPON: 0.02375\n\
FREQUENCY: FrequencyTypes.SEMI_ANNUAL\n\
ACCRUAL TYPE: DayCountTypes.ACT_ACT_ICMA\n\
FACE AMOUNT: 1000000"

latex_str = convertToLatexTable(txt, ":", ["FIELD", "VALUE"])
print(latex_str)
