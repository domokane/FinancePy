import re


def print_table_start(columns_str):
    text = (
        "\\begin{center}\n"
        "\\begin{tabular}{" + columns_str + "}\n"
        "\\toprule\n"
    )
    print(text, end="")


PRECISION = 6


def print_table_header(*args):

    values = []

    for arg in args:
        text = "\\textbf{" + str(arg) + "}"
        values.append(text)

    print(" & ".join(values) + r" \\")
    print(r"\midrule")


def print_table_row(*args):

    values = []

    for arg in args:

        if isinstance(arg, float):
            text = f"{arg:.{PRECISION}f}"
        elif isinstance(arg, int):
            text = str(arg)
        else:
            text = str(arg)

        values.append(text)

    print(" & ".join(values) + r" \\")


def print_table_end():
    text = (
        "\\bottomrule\n"
        "\\end{tabular}\n"
        "\\end{center}\n"
        "\\noindent\n"
    )
    print(text, end="")


def print_latex_text(obj):
    text = str(obj)

    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }

    text = re.sub(
        r"[\\&%$#_{}~^]",
        lambda m: replacements[m.group()],
        text
    )

    print(r"\par\noindent")

    for line in text.splitlines():
        print(line + r"\par")

    print(r"\medskip")


def print_object_to_latex(obj, caption=None, label=None):
    """Convert a FinancePy object's string representation to a LaTeX table."""

    name = type(obj).__name__

    def escape_latex(s):
        chars = {
            "\\": r"\textbackslash{}",
            "&": r"\&",
            "%": r"\%",
            "$": r"\$",
            "#": r"\#",
            "_": r"\_",
            "{": r"\{",
            "}": r"\}",
            "~": r"\textasciitilde{}",
            "^": r"\textasciicircum{}",
        }
        return "".join(chars.get(c, c) for c in str(s))

    rows = []

    for line in str(obj).splitlines():
        if ":" not in line:
            continue

        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()

        if key:
            rows.append((
                escape_latex(key),
                escape_latex(format_value(value))
            ))

    # Exclude OBJECT TYPE
    rows = [(k, v) for k, v in rows if k != "OBJECT TYPE"]

    mid = (len(rows) + 1) // 2
    left_rows = rows[:mid]
    right_rows = rows[mid:]

    name_line = rf"\begin{{termsheet}}{{{name}}}"

    lines = [
        name_line,
        r"\centering",
        r"\renewcommand{\arraystretch}{1.05}",
        r"\begin{tabular}{@{}ll@{\hspace{1.5em}}ll@{}}"
    ]

    for i in range(mid):
        key1, value1 = left_rows[i]

        if i < len(right_rows):
            key2, value2 = right_rows[i]
        else:
            key2, value2 = "", ""

        lines.append(
            f"{key1} & \\textbf{{{value1}}} & "
            f"{key2} & \\textbf{{{value2}}} " + r"\\"
        )

    lines.append(r"\end{tabular}")

    if caption is not None:
        lines.append(r"\caption{" + caption + "}")

    if label is not None:
        lines.append(r"\label{" + label + "}")

    lines.append(r"\end{termsheet}")

    out = "\n".join(lines)
    print(out)


def format_value_old(value):
    value = value.strip()

    match = re.fullmatch(
        r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)(\s*[a-zA-Z%]*)",
        value,
    )

    if not match:
        return value

    number, unit = match.groups()
    x = float(number)

    if x.is_integer():
        return f"{int(x):,}{unit}"

    return f"{x:,.5f}{unit}"


DAY_COUNT_NAMES = {
    "ACT_360": "ACT/360",
    "ACT_365F": "ACT/365F",
    "ACT_ACT_ISDA": "ACT/ACT (ISDA)",
    "ACT_ACT_ICMA": "ACT/ACT (ICMA)",
    "ACT_ACT_AFB": "ACT/ACT (AFB)",
    "THIRTY_360_BOND": "30/360 (Bond)",
    "THIRTY_E_360": "30E/360",
    "THIRTY_E_360_ISDA": "30E/360 (ISDA)",
    "THIRTY_E_PLUS_360": "30E+/360",
}


def format_value(value):
    value = value.strip()

    if value == "True":
        return "Yes"
    if value == "False":
        return "No"

    # Remove FinancePy enum type prefix
    # e.g. DayCountTypes.ACT_360 -> ACT_360
    value = re.sub(r"\b\w+Types\.", "", value)

    value = re.sub(r"\b\w+Types\.", "", value)

    if value in DAY_COUNT_NAMES:
        return DAY_COUNT_NAMES[value]

    match = re.fullmatch(
        r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)(\s*[a-zA-Z%]*)",
        value,
    )

    if not match:
        return value

    number, unit = match.groups()
    x = float(number)

    if x.is_integer():
        return f"{int(x):,}{unit}"

    return f"{x:,.5f}{unit}"
