from hmac import new
import os
import re
import ast
import difflib


def to_snake_case(name: str) -> str:
    s1 = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    s2 = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s1)
    return s2.lower()


def rename_functions_and_calls(source: str) -> str:
    """
    Rename functions (definitions and calls) to snake_case.
    Class names and calls to classes are NOT renamed.
    Preserves original line breaks and parentheses.
    """
    tree = ast.parse(source)

    function_names = set()
    class_names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            if not (node.name.startswith("__") and node.name.endswith("__")):
                function_names.add(node.name)
        elif isinstance(node, ast.ClassDef):
            class_names.add(node.name)

    function_names -= class_names  # exclude class names

    rename_map = {}
    for name in function_names:
        if name.startswith("_"):
            rename_map[name] = "_" + to_snake_case(name[1:])
        else:
            rename_map[name] = to_snake_case(name)

    # Replace function definitions
    source = re.sub(
        r"\bdef\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(",
        lambda m: f"def {rename_map.get(m.group(1), m.group(1))}(",
        source,
    )

    # Replace function calls (bare names, not attributes)
    source = re.sub(
        r"\b([A-Za-z_][A-Za-z0-9_]*)\s*\(",
        lambda m: f"{rename_map.get(m.group(1), m.group(1))}(",
        source,
    )

    return source


def ensure_blank_line_after_defs_and_classes(source: str) -> str:
    """
    Ensure there is a blank line after each 'def' or 'class' line.
    If there is already a blank line, leave it as-is.
    """
    lines = source.splitlines()
    new_lines = []
    n = len(lines)

    i = 0
    while i < n:
        line = lines[i]
        new_lines.append(line)

        stripped = line.lstrip()
        if stripped.startswith(("def ", "class ")):
            # Check if next line exists and is not blank
            if i + 1 < n and lines[i + 1].strip() != "":
                new_lines.append("")  # insert blank line
        i += 1

    return "\n".join(new_lines) + "\n"


def strip_hash_lines(source: str) -> str:
    """
    Remove all full-line hashes (lines consisting only of # characters),
    but leave normal comments intact.
    """
    lines = source.splitlines()
    cleaned_lines = [
        line for line in lines if not (line.strip() and set(line.strip()) == {"#"})
    ]
    return "\n".join(cleaned_lines) + "\n"


def remove_blank_lines_before_defs(source: str) -> str:
    """
    Remove blank lines immediately before a function or class definition.
    """
    lines = source.splitlines()
    new_lines = []
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith(("def ", "class ")):
            # Remove preceding blank lines
            while new_lines and new_lines[-1].strip() == "":
                new_lines.pop()
        new_lines.append(line)
    return "\n".join(new_lines) + "\n"


def normalize_hashes_and_functions(source: str) -> str:
    """
    Normalize Python code:
    - Remove existing hash lines
    - Remove blank lines before function/class definitions
    - Add a single blank line + 88 '#' + 2 blank lines before each function/class
    - Ensure file ends with a single newline
    """
    lines = source.splitlines()
    new_lines = []
    i = 0
    n = len(lines)

    while i < n:
        line = lines[i]
        stripped = line.lstrip()

        # Skip existing full-line hashes
        if stripped.startswith("#") and set(stripped.strip()) == {"#"}:
            i += 1
            continue

        is_def_or_class = stripped.startswith(("def ", "class "))

        if is_def_or_class:
            # Remove preceding blank lines
            while new_lines and new_lines[-1].strip() == "":
                new_lines.pop()

            # Add one blank line, 88 '#'s, 2 blank lines
            new_lines.append("")
            new_lines.append("#" * 88)
            new_lines.append("")
            new_lines.append("")

            # Add the def/class line
            new_lines.append(line)
            i += 1
            continue

        # Normal line
        new_lines.append(line)
        i += 1

    # Ensure file ends with single newline
    if not new_lines or new_lines[-1].strip() != "":
        new_lines.append("")

    return "\n".join(new_lines)


def add_hash_before_if_main(source: str) -> str:
    """
    Insert a full line of 88 '#' two spaces before 'if __name__ == "__main__":'.
    """
    lines = source.splitlines()
    new_lines = []

    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("if __name__") and "__main__" in stripped:
            indent = max(0, len(line) - len(stripped) - 2)
            new_lines.append(" " * indent + "#" * 88)
            new_lines.append("")
        new_lines.append(line)

    new_lines.append("")  # Ensure file ends with a newline
    return "\n".join(new_lines) + "\n"


def add_hash_on_dedent_to_0(source: str) -> str:
    """
    Insert a line of 88 '#' whenever indentation drops from >=4 spaces to 0.
    """
    lines = source.splitlines()
    new_lines = []

    prev_indent = 0
    for line in lines:
        stripped = line.lstrip()
        curr_indent = len(line) - len(stripped)

        # Skip empty lines
        if stripped == "":
            new_lines.append(line)
            continue

        # Insert hash line only if dedent is from >=4 to 0
        if prev_indent >= 4 and curr_indent == 0:
            new_lines.append("#" * 88)
            new_lines.append("")

        new_lines.append(line)
        prev_indent = curr_indent

    return "\n".join(new_lines) + "\n"


def transform(source: str) -> str:
    source_new = strip_hash_lines(source)
    source_new = rename_functions_and_calls(source_new)
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

        if source != new_source and False:
            diff = difflib.unified_diff(
                source.splitlines(),
                new_source.splitlines(),
                fromfile="original",
                tofile="pep8",
                lineterm="",
            )
            print(f"Changes in {filepath}:\n{'\n'.join(diff)}\n")

        folder = os.path.dirname(filepath)
        filename = os.path.basename(filepath)
        base, ext = os.path.splitext(filename)
        target_folder = os.path.join(folder, "../pep8_output")
        os.makedirs(target_folder, exist_ok=True)
        new_filepath = os.path.join(target_folder, base + "_2.py")

        with open(new_filepath, "w", encoding="utf-8") as f:
            f.write(new_source)

        print(f"Written to {new_filepath} (overwritten if existed)")

    except (IOError, OSError) as e:
        print(f"Error processing {filepath}: {e}")


def convert_folder(folder_path: str) -> None:
    """Process all .py files in a folder, skipping already converted _2.py files."""
    try:
        for filename in os.listdir(folder_path):
            if filename.endswith(".py") and not filename.endswith("_2.py"):
                convert_file(os.path.join(folder_path, filename))
    except (IOError, OSError) as e:
        print(f"Error accessing folder {folder_path}: {e}")


if __name__ == "__main__":
    FOLDER_PATH = "../pep8_tests"
    convert_folder(FOLDER_PATH)
    print("PEP8 conversion complete.")
