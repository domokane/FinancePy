from hmac import new
from operator import is_
import os
import re
import ast
import difflib
from typing import Dict, Set

# Configuration constants
HASH_LINE_LENGTH = 88  # Standard line width for PEP 8
OUTPUT_DIR = "./pep8_output"
SHOW_DIFF = False  # Toggle to show diff output


def to_snake_case(name: str) -> str:
    """Convert a string to snake_case."""
    s1 = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    s2 = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s1)
    return s2.lower()


def rename_functions_and_calls(source: str) -> str:
    """Rename function definitions and calls to snake_case, excluding classes."""
    try:
        tree = ast.parse(source)
    except SyntaxError as e:
        print(f"Syntax error in source code: {e}")
        return source

    function_names: Set[str] = set()
    class_names: Set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            if not (node.name.startswith("__") and node.name.endswith("__")):
                function_names.add(node.name)
        elif isinstance(node, ast.ClassDef):
            class_names.add(node.name)
    function_names -= class_names

    rename_map: Dict[str, str] = {
        name: (
            "_" + to_snake_case(name[1:])
            if name.startswith("_")
            else to_snake_case(name)
        )
        for name in function_names
    }

    # Replace function definitions
    source = re.sub(
        r"\bdef\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(",
        lambda m: f"def {rename_map.get(m.group(1), m.group(1))}(",
        source,
    )

    # Replace function calls
    def replace_call(match):
        name = match.group(1)
        if name in rename_map and name not in class_names:
            return f"{rename_map[name]}("
        return match.group(0)

    source = re.sub(
        r"(?<!\.)(\b[A-Za-z_][A-Za-z0-9_]*)\s*\(",
        replace_call,
        source,
    )
    return source


def rename_variables(source: str) -> str:
    """Rename variable names in assignments to snake_case, excluding classes and functions."""
    try:
        tree = ast.parse(source)
    except SyntaxError as e:
        print(f"Syntax error in source code: {e}")
        return source

    function_names: Set[str] = set()
    class_names: Set[str] = set()
    variable_names: Set[str] = set()

    # Collect function, class, and variable names
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            if not (node.name.startswith("__") and node.name.endswith("__")):
                function_names.add(node.name)
        elif isinstance(node, ast.ClassDef):
            class_names.add(node.name)
        elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            if not (node.id.startswith("__") and node.id.endswith("__")):
                variable_names.add(node.id)

    variable_names -= class_names | function_names  # Exclude class and function names

    rename_map: Dict[str, str] = {
        name: (
            "_" + to_snake_case(name[1:])
            if name.startswith("_")
            else to_snake_case(name)
        )
        for name in variable_names
    }

    # Replace variable names (standalone identifiers, not in attribute chains or function calls)
    def replace_variable(match):
        name = match.group(1)
        if name in rename_map and name not in class_names | function_names:
            return rename_map[name]
        return match.group(0)

    source = re.sub(
        r"(?<!\.)(\b[A-Za-z_][A-Za-z0-9_]*)\b(?!\s*\()",
        replace_variable,
        source,
    )
    return source


def ensure_blank_line_after_defs_and_classes(source: str) -> str:
    """Ensure a blank line after each def or class line."""
    lines = source.splitlines()
    new_lines = []
    n = len(lines)
    i = 0
    while i < n:
        line = lines[i]
        new_lines.append(line)
        stripped = line.lstrip()
        if stripped.startswith(("def ", "class ")):
            if i + 1 < n and lines[i + 1].strip() != "":
                new_lines.append("")
        i += 1
    return "\n".join(new_lines) + "\n"


def strip_hash_lines(source: str) -> str:
    """Remove lines consisting only of # characters."""
    lines = source.splitlines()
    cleaned_lines = [
        line for line in lines if not (line.strip() and set(line.strip()) == {"#"})
    ]
    return "\n".join(cleaned_lines) + "\n"


def remove_blank_lines_before_defs(source: str) -> str:
    """Remove blank lines before function or class definitions."""
    lines = source.splitlines()
    new_lines = []
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith(("def ", "class ")):
            while new_lines and new_lines[-1].strip() == "":
                new_lines.pop()
        new_lines.append(line)
    return "\n".join(new_lines) + "\n"


def normalize_hashes_and_functions(source: str) -> str:
    """Normalize code with hash lines before functions/classes."""
    lines = source.splitlines()
    new_lines = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        stripped = line.lstrip()
        if stripped.startswith("#") and set(stripped.strip()) == {"#"}:
            i += 1
            continue

        is_vectorize = stripped.startswith("@vectorize")

        if is_vectorize:
            while new_lines and new_lines[-1].strip() == "":
                new_lines.pop()

            new_lines.extend(["", "#" * HASH_LINE_LENGTH, "", ""])

            while lines[i].find("def ") == -1 and i < n:
                print("VECTORIZE LINE:", lines[i])
                new_lines.append(lines[i])
                i += 1

                continue

            while lines[i].strip() != "" and i < n:
                new_lines.append(lines[i])
                i += 1
                continue

            line = lines[i]
            stripped = line.lstrip()

        is_njit = stripped.startswith("@njit")

        if is_njit:
            while new_lines and new_lines[-1].strip() == "":
                new_lines.pop()
            new_lines.extend(["", "#" * HASH_LINE_LENGTH, "", ""])
            new_lines.append(lines[i])
            i += 1
            new_lines.append(lines[i])
            i += 1
            continue

        is_def = stripped.startswith("def ")

        if is_def:
            while new_lines and new_lines[-1].strip() == "":
                new_lines.pop()
            new_lines.extend(["", "#" * HASH_LINE_LENGTH, "", ""])
            new_lines.append(line)
            i += 1
            continue

        is_class = stripped.startswith("class ")

        if is_class:
            while new_lines and new_lines[-1].strip() == "":
                new_lines.pop()
            new_lines.extend(["", "#" * HASH_LINE_LENGTH, "", ""])
            new_lines.append(line)
            i += 1
            continue

        new_lines.append(line)
        i += 1
    if not new_lines or new_lines[-1].strip() != "":
        new_lines.append("")
    return "\n".join(new_lines)


def add_hash_before_if_main(source: str) -> str:
    """Insert a hash line before 'if __name__ == "__main__":'."""
    lines = source.splitlines()
    new_lines = []
    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("if __name__") and "__main__" in stripped:
            indent = len(line) - len(stripped)
            new_lines.extend([" " * indent + "#" * HASH_LINE_LENGTH, ""])
        new_lines.append(line)
    new_lines.append("")
    return "\n".join(new_lines)


def add_hash_on_dedent_to_0(source: str) -> str:
    """Insert a hash line when indentation drops from >=4 to 0, only if no def/class follows."""
    lines = source.splitlines()
    new_lines = []
    prev_indent = 0
    n = len(lines)
    for i, line in enumerate(lines):  # Fixed: Use enumerate correctly
        stripped = line.lstrip()
        if stripped == "":
            new_lines.append(line)
            continue
        curr_indent = len(line) - len(stripped)
        # Check if there are any def or class statements in the remaining lines
        has_future_defs = False
        for j in range(i + 1, n):
            future_stripped = lines[j].lstrip()
            if future_stripped.startswith(("def ", "class ")):
                has_future_defs = True
                break
        # Insert hash line only if indentation drops to 0, no def/class follows, and not a def/class line
        if (
            prev_indent >= 4
            and curr_indent == 0
            and not has_future_defs
            and not stripped.startswith(("def ", "class "))
        ):
            new_lines.extend(["#" * HASH_LINE_LENGTH, ""])
        new_lines.append(line)
        prev_indent = curr_indent
    return "\n".join(new_lines) + "\n"


def transform(source: str) -> str:
    """Apply all transformations to the source code."""
    source_new = strip_hash_lines(source)
    source_new = rename_functions_and_calls(source_new)
    source_new = rename_variables(source_new)
    source_new = ensure_blank_line_after_defs_and_classes(source_new)
    source_new = remove_blank_lines_before_defs(source_new)
    source_new = normalize_hashes_and_functions(source_new)
    source_new = add_hash_before_if_main(source_new)
    source_new = add_hash_on_dedent_to_0(source_new)
    return source_new


def convert_file(filepath: str) -> None:
    """Convert a single file, overwriting any existing _2.py file."""
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            source = f.read()
        print(f"Converting file: {filepath}")
        new_source = transform(source)

        # Show diff if enabled
        if SHOW_DIFF and source != new_source:
            diff = difflib.unified_diff(
                source.splitlines(),
                new_source.splitlines(),
                fromfile="original",
                tofile="pep8",
                lineterm="",
            )
            print(f"Changes in {filepath}:\n{'\n'.join(diff)}\n")

        os.makedirs(OUTPUT_DIR, exist_ok=True)
        base, ext = os.path.splitext(os.path.basename(filepath))

        bkp_filepath = os.path.join(OUTPUT_DIR, "../backup/", f"{base}.py")
        with open(bkp_filepath, "w", encoding="utf-8") as f:
            f.write(source)
        print(f"Backup Written to {bkp_filepath}")

        new_filepath = os.path.join(OUTPUT_DIR, "../converted/", f"{base}.py")
        with open(new_filepath, "w", encoding="utf-8") as f:
            f.write(new_source)
        print(f"Written to {new_filepath}")

    except (IOError, OSError) as e:
        print(f"Error processing {filepath}: {e}")


def convert_folder(folder_path: str) -> None:
    """Process all .py files in a folder."""
    try:
        for filename in os.listdir(folder_path):
            if filename.endswith(".py"):
                convert_file(os.path.join(folder_path, filename))
    except (IOError, OSError) as e:
        print(f"Error accessing folder {folder_path}: {e}")


if __name__ == "__main__":
    TARGET_FOLDER_RELATIVE_PATH = "../golden_tests"
    TARGET_FOLDER_RELATIVE_PATH = "../unit_tests"
    TARGET_FOLDER_RELATIVE_PATH = "../financepy/market/curves"
    TARGET_FOLDER_RELATIVE_PATH = "../financepy/market/volatility"
    TARGET_FOLDER_RELATIVE_PATH = "../financepy/models"
    convert_folder(TARGET_FOLDER_RELATIVE_PATH)
    print("PEP8 conversion complete.")
