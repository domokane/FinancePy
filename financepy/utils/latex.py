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

    name_line = rf"\begin{{termsheet}}{{{name}}}"

    lines = [
        name_line,
        r"\centering",
        r"\begin{tabular}{ll}",
    ]

    # We omit Object Type as it is in header of term sheet
    for key, value in rows[1:]:
        lines.append(f"{key} & {value} " + r"\\")

    lines.extend([
        r"\end{tabular}",
    ])

    if caption is not None:
        lines.append(r"\caption{" + caption + "}")

    if label is not None:
        lines.append(r"\label{" + label + "}")

    lines.append(r"\end{termsheet}")

    out = "\n".join(lines)
    print(out)


def format_value(value):
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
