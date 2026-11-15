##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################
"""Optional, context-local valuation reporting for FinancePy.

Example::

    with ValuationLog(verbosity=3) as log:
        pv = product.value(...)

    log.print()
    print(log.to_markdown())
    print(log.to_latex())

Use log_value(), log_table(), log_row(), and log_section() within pricing
methods. All functions are no-ops if no valuation log is active.
"""

from contextlib import contextmanager
from numbers import Real
from contextvars import ContextVar
import sys
import re


@contextmanager
def suppress_logging():
    token = _active_log.set(None)
    try:
        yield
    finally:
        _active_log.reset(token)


_active_log = ContextVar("financepy_valuation_log", default=None)


def _caller_info():
    """Find the class and method calling a public log helper."""
    frame = sys._getframe(2)
    try:
        obj = frame.f_locals.get("self")
        return {
            "class": type(obj).__name__ if obj is not None else None,
            "function": frame.f_code.co_name,
        }
    finally:
        del frame


def _check_level(level):
    if not isinstance(level, int) or isinstance(level, bool) or not 1 <= level <= 5:
        raise ValueError("level must be an integer between 1 and 5")


def _escape_markdown(value):
    return str(value).replace("\\", "\\\\").replace("|", "\\|").replace("\n", " ")


def _escape_latex(value):
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
    return re.sub(r"[\\&%$#_{}~^]", lambda match: replacements[match.group()], str(value))


class ValuationLog:
    """Record valuation results and render them as text, Markdown or LaTeX.

    verbosity: 1 (summary) through 5 (maximum detail).
    Only entries at or below the selected verbosity are retained.
    """

    def __init__(self, decimals=6, min_width=12, verbosity=3):
        _check_level(verbosity)
        if not isinstance(decimals, int) or decimals < 0:
            raise ValueError("decimals must be a non-negative integer")
        if not isinstance(min_width, int) or min_width < 1:
            raise ValueError("min_width must be a positive integer")
        self.entries = []
        self.decimals = decimals
        self.min_width = min_width
        self.verbosity = verbosity
        self._tokens = []
        self._table_open = False
        self._table_enabled = False
        self._table_columns = 0

    def __enter__(self):
        self._tokens.append(_active_log.set(self))
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        _active_log.reset(self._tokens.pop())
        return False

    def clear(self):
        """Discard collected results without leaving the active context."""
        self.entries.clear()
        self._table_open = False
        self._table_enabled = False
        self._table_columns = 0

    def _format(self, value):
        if isinstance(value, Real) and not isinstance(value, bool):
            return f"{float(value):.{self.decimals}f}"
        return str(value)

    def _blocks(self, show_methods):
        """Yield (kind, data) for renderers, preserving calculation order."""
        previous_method = None
        i = 0
        while i < len(self.entries):
            entry = self.entries[i]
            kind = entry["type"]

            if kind == "section":
                yield "section", entry["title"]
                previous_method = None
                i += 1
                continue

            if kind == "row":
                raise ValueError("Row encountered without a table header")

            method = f"{entry.get('class') or '<module>'}.{entry.get('function') or '<unknown>'}"
            if show_methods and method != previous_method:
                yield "method", method
            previous_method = method

            if kind == "value":
                pairs = list(entry["values"])
                i += 1
                # Consecutive entries from this method share a single table.
                while i < len(self.entries):
                    next_entry = self.entries[i]
                    if next_entry["type"] != "value":
                        break
                    next_method = (
                        f"{next_entry.get('class') or '<module>'}."
                        f"{next_entry.get('function') or '<unknown>'}"
                    )
                    if show_methods and next_method != method:
                        break
                    pairs.extend(next_entry["values"])
                    i += 1
                yield "values", pairs
            elif kind == "table":
                rows = []
                i += 1
                while i < len(self.entries) and self.entries[i]["type"] == "row":
                    rows.append(self.entries[i]["values"])
                    i += 1
                yield "table", (entry["headers"], rows)
            else:
                raise ValueError(f"Unknown valuation log entry type: {kind}")

    def _table_lines(self, headers, rows):
        widths = [
            max(
                self.min_width,
                len(str(header)),
                *(len(self._format(row[i])) for row in rows),
            )
            for i, header in enumerate(headers)
        ]

        def format_cells(values, is_header=False):
            parts = []
            for i, value in enumerate(values):
                display = str(value) if is_header else self._format(value)
                is_numeric = isinstance(value, Real) and not isinstance(value, bool)
                if is_header or is_numeric:
                    parts.append(f"{display:>{widths[i]}}")
                else:
                    parts.append(f"{display:<{widths[i]}}")
            return "  ".join(parts)

        yield format_cells(headers, is_header=True)
        yield "-" * (sum(widths) + 2 * (len(widths) - 1))
        for row in rows:
            yield format_cells(row)

    def to_string(self, show_methods=True):
        """Return a readable plain-text valuation report."""
        lines = []
        for kind, data in self._blocks(show_methods):
            if kind in ("section", "method"):
                if lines:
                    lines.append("")
                lines.append(data)
            elif kind == "values":
                for name, value in data:
                    lines.append(f"{name}: {self._format(value)}")
            elif kind == "table":
                headers, rows = data
                lines.extend(self._table_lines(headers, rows))
        return "\n".join(lines)

    def print(self, show_methods=True):
        print(self.to_string(show_methods=show_methods))

    def to_markdown(self, show_methods=False):
        """Return Markdown tables suitable for a book manuscript."""
        lines = []
        for kind, data in self._blocks(show_methods):
            if kind in ("section", "method"):
                if lines:
                    lines.append("")
                lines.append("## " + _escape_markdown(data))
                lines.append("")
            elif kind == "values":
                lines.extend(("| Parameter | Value |", "|:--|--:|"))
                for name, value in data:
                    lines.append(
                        f"| {_escape_markdown(name)} | {_escape_markdown(self._format(value))} |"
                    )
                lines.append("")
            elif kind == "table":
                headers, rows = data
                lines.append("| " + " | ".join(_escape_markdown(x) for x in headers) + " |")
                lines.append("| " + " | ".join("---:" for _ in headers) + " |")
                for row in rows:
                    lines.append("| " + " | ".join(
                        _escape_markdown(self._format(x)) for x in row
                    ) + " |")
                lines.append("")
        return "\n".join(lines).strip()

    def to_latex(self, show_methods=False):
        """Return standalone LaTeX tabular fragments (no formula generation)."""
        lines = []

        def table(headers, rows, show_header=True):
            alignment = "".join(
                "r" if rows and all(
                    isinstance(row[i], Real) and not isinstance(row[i], bool)
                    for row in rows
                ) else "l"
                for i in range(len(headers))
            )

            out = [
                r"\begin{tabular}{" + alignment + "}",
                r"\hline",
            ]

            if show_header:
                out.append(" & ".join(_escape_latex(x) for x in headers) + r" \\")
                out.append(r"\hline")

            for row in rows:
                out.append(
                    " & ".join(_escape_latex(self._format(x)) for x in row) + r" \\"
                )

            out.extend((r"\hline", r"\end{tabular}"))
            return out

        for kind, data in self._blocks(show_methods):
            if kind == "section":
                lines.extend((r"\textbf{" + _escape_latex(data) + "}", ""))
                lines.append(r"\par\medskip")
            elif kind == "method":
                lines.extend((r"\textbf{" + _escape_latex(data) + "}", ""))
            elif kind == "values":
                lines.extend(table(("Parameter", "Value"), data, show_header=False))
                lines.append(r"\par\medskip")
            elif kind == "table":
                headers, rows = data
                lines.extend(table(headers, rows))
                lines.append(r"\par\medskip")

        return "\n".join(lines).strip()


def log_section(title, level=1):
    """Start a semantic report section, such as 'Forward CDS Valuation'."""

    logger = _active_log.get()
    if logger is None:
        return
    _check_level(level)
    logger._table_open = False
    if level <= logger.verbosity:
        logger.entries.append({"type": "section", "title": str(title), "level": level})


def log_value(*args, level=1):
    """Record one or more name/value pairs at the requested verbosity level."""

    logger = _active_log.get()
    if logger is None:
        return
    _check_level(level)
    logger._table_open = False
    if level > logger.verbosity:
        return
    if not args or len(args) % 2:
        raise ValueError("log_value requires one or more name/value pairs")
    logger.entries.append({
        **_caller_info(),
        "type": "value",
        "level": level,
        "values": [(args[i], args[i + 1]) for i in range(0, len(args), 2)],
    })


def log_table(*headers, level=4):
    """Start a table. Consecutive log_row() calls append rows to it."""

    logger = _active_log.get()
    if logger is None:
        return
    _check_level(level)
    if not headers:
        raise ValueError("log_table requires at least one header")
    logger._table_open = True
    logger._table_enabled = level <= logger.verbosity
    logger._table_columns = len(headers)
    if logger._table_enabled:
        logger.entries.append({
            **_caller_info(),
            "type": "table",
            "level": level,
            "headers": tuple(headers),
        })


def log_row(*values):
    """Append a row to the last table, including when that table is hidden."""

    logger = _active_log.get()
    if logger is None:
        return
    if not logger._table_open:
        raise ValueError("log_row() requires a preceding log_table()")
    if not logger._table_enabled:
        return
    if len(values) != logger._table_columns:
        raise ValueError(
            f"Table row has {len(values)} values; expected {logger._table_columns}"
        )
    logger.entries.append({"type": "row", "values": tuple(values)})


def log_enter(level=1):

    logger = _active_log.get()

    if logger is None or level > logger.verbosity:
        return

    logger.entries.append({
        **_caller_info(),
        "type": "enter",
        "level": level,
    })


def log_amount(name, value, decimals=2, level=1):

    log_value(name, f"{value:,.{decimals}f}", level=level)


def log_percentage(name, value, decimals=2, level=1):

    log_value(name, f"{100 * value:.{decimals}f}%", level=level)


def log_bool(name, value, level=1):

    log_value(name, bool(value), level=level)


def log_float(name, value, decimals=6, level=1):

    log_value(name, f"{value:.{decimals}f}", level=level)


def log_int(name, value, level=1):

    log_value(name, f"{value:.0f}", level=level)


def log_date(name, value, level=1):
    """Log a date using its string representation."""

    log_value(name, str(value), level=level)


def log_bps(name, value, decimals=2, level=1):
    """Log a spread or rate expressed in basis points."""

    log_value(name, f"{value*10000:.{decimals}f} bps", level=level)
