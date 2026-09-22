# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

"""Generate HTML API documentation for FinancePy.

The generator scans the FinancePy source package and creates a browsable
HTML API reference from public modules, classes, functions, signatures,
and docstrings.

The generated documentation mirrors the FinancePy package hierarchy.

For example:

    financepy.products
        -> html/products/index.html

    financepy.products.bonds
        -> html/products/bonds/index.html

    financepy.products.bonds.bond
        -> html/products/bonds/bond.html

Run from the repository root using:

    python scripts/generate_api_docs.py

The generated HTML files should not be edited manually.
"""

from __future__ import annotations

import html
import importlib
import inspect
import os
import pkgutil
import shutil
import sys
from pathlib import Path


###############################################################################
# DIRECTORIES
###############################################################################

SCRIPT_FILE = Path(__file__).resolve()

# financepy-git/scripts/
SCRIPT_DIR = SCRIPT_FILE.parent

# financepy-git/
PROJECT_ROOT = SCRIPT_DIR.parent

# financepy-git/financepy/
FINANCEPY_DIR = PROJECT_ROOT / "financepy"

# financepy-git/html/
HTML_DIR = PROJECT_ROOT / "html"

PACKAGE_NAME = "financepy"


###############################################################################
# CONFIGURATION
###############################################################################

EXCLUDED_DIRECTORIES = {
    "__pycache__",
    "html",
    "scripts",
}

EXCLUDED_MODULES = {
    "__init__",
}


###############################################################################
# HTML STYLE
###############################################################################

STYLE = """
<style>

body {
    font-family: Arial, Helvetica, sans-serif;
    margin: 0;
    padding: 0;
    color: #222;
    background: #ffffff;
}

header {
    background: #20232a;
    color: white;
    padding: 20px 40px;
}

header h1 {
    margin: 0;
    font-size: 28px;
}

header p {
    margin: 8px 0 0 0;
    color: #cccccc;
}

main {
    max-width: 1200px;
    margin: 0 auto;
    padding: 30px 40px 80px 40px;
}

nav {
    margin-bottom: 30px;
    font-size: 95%;
}

nav a {
    margin-right: 6px;
}

a {
    color: #0057a8;
    text-decoration: none;
}

a:hover {
    text-decoration: underline;
}

h1 {
    margin-top: 0;
}

h2 {
    margin-top: 40px;
    padding-bottom: 6px;
    border-bottom: 1px solid #dddddd;
}

h3 {
    margin-top: 30px;
}

h4 {
    margin-top: 24px;
    margin-bottom: 8px;
}

code {
    font-family: Consolas, "Courier New", monospace;
}

pre {
    background: #f5f5f5;
    border: 1px solid #dddddd;
    border-radius: 4px;
    padding: 14px;
    overflow-x: auto;
    line-height: 1.4;
}

.signature {
    background: #f5f5f5;
    border-left: 4px solid #777777;
    padding: 12px;
    margin: 10px 0 18px 0;
    overflow-x: auto;
}

.docstring {
    white-space: pre-wrap;
    line-height: 1.5;
}

.package-list,
.module-list {
    line-height: 1.9;
}

.member {
    margin-bottom: 35px;
}

.inheritance {
    color: #555555;
    margin: 8px 0 14px 0;
}

.source-name {
    color: #666666;
    font-family: Consolas, "Courier New", monospace;
    font-size: 90%;
}

.card-grid {
    display: grid;
    grid-template-columns: repeat(
        auto-fit,
        minmax(220px, 1fr)
    );
    gap: 14px;
    margin-top: 20px;
}

.card {
    border: 1px solid #dddddd;
    border-radius: 5px;
    padding: 16px;
    background: #fafafa;
}

.card h3 {
    margin-top: 0;
    margin-bottom: 8px;
}

.card p {
    margin-bottom: 0;
    color: #555555;
}

.generated {
    margin-top: 60px;
    padding-top: 15px;
    border-top: 1px solid #dddddd;
    color: #777777;
    font-size: 90%;
}

</style>
"""


###############################################################################
# BASIC UTILITIES
###############################################################################


def escape(value) -> str:
    """HTML-escape a value."""

    return html.escape(str(value))


###############################################################################


def display_name(name: str) -> str:
    """Convert an internal Python name into a readable display name."""

    return name.replace("_", " ").title()


###############################################################################


def is_public_name(name: str) -> bool:
    """Return True if a Python name should appear in the documentation."""

    return not name.startswith("_")


###############################################################################


def get_docstring(obj) -> str:
    """Return a cleaned docstring."""

    doc = inspect.getdoc(obj)

    if doc is None:
        return ""

    return doc


###############################################################################


def get_signature(obj) -> str:
    """Return an object's signature if available."""

    try:
        return str(inspect.signature(obj))
    except (TypeError, ValueError):
        return ""


###############################################################################


def relative_link(from_file: Path, to_file: Path) -> str:
    """Return a relative HTML link between generated files."""

    return Path(
        os.path.relpath(
            to_file,
            start=from_file.parent,
        )
    ).as_posix()


###############################################################################
# MODULE AND PACKAGE PATHS
###############################################################################


def module_parts(module_name: str) -> list[str]:
    """Return module components excluding the financepy prefix."""

    parts = module_name.split(".")

    if parts and parts[0] == PACKAGE_NAME:
        parts = parts[1:]

    return parts


###############################################################################


def module_to_output_path(module_name: str) -> Path:
    """Return the generated HTML path for a Python module.

    Example
    -------
    financepy.products.bonds.bond

    becomes

    html/products/bonds/bond.html
    """

    parts = module_parts(module_name)

    return HTML_DIR.joinpath(*parts).with_suffix(".html")


###############################################################################


def package_to_output_path(package_name: str) -> Path:
    """Return the generated index path for a Python package.

    Example
    -------
    financepy.products.bonds

    becomes

    html/products/bonds/index.html
    """

    parts = module_parts(package_name)

    if not parts:
        return HTML_DIR / "index.html"

    return HTML_DIR.joinpath(*parts) / "index.html"


###############################################################################
# COMMON HTML
###############################################################################


def page_header(title: str, subtitle: str = "") -> str:
    """Return the common HTML page header."""

    subtitle_html = ""

    if subtitle:
        subtitle_html = f"<p>{escape(subtitle)}</p>"

    return f"""<!DOCTYPE html>
<html lang="en">

<head>

<meta charset="utf-8">

<meta
    name="viewport"
    content="width=device-width, initial-scale=1"
>

<title>{escape(title)} - FinancePy</title>

{STYLE}

</head>

<body>

<header>

<h1>FinancePy</h1>

{subtitle_html}

</header>

<main>
"""


###############################################################################


def page_footer() -> str:
    """Return the common HTML page footer."""

    return """
<div class="generated">

Generated automatically from the FinancePy source code.
Do not edit this file manually.

</div>

</main>

</body>

</html>
"""


###############################################################################
# DISCOVERY
###############################################################################


def discover_modules() -> list[str]:
    """Discover public FinancePy modules and packages."""

    modules = []

    for module_info in pkgutil.walk_packages(
        [str(FINANCEPY_DIR)],
        prefix=f"{PACKAGE_NAME}.",
    ):

        module_name = module_info.name
        parts = module_name.split(".")

        if any(
            part in EXCLUDED_DIRECTORIES
            for part in parts
        ):
            continue

        if parts[-1] in EXCLUDED_MODULES:
            continue

        modules.append(module_name)

    return sorted(modules)


###############################################################################


def discover_packages(
    modules: list[str],
) -> set[str]:
    """Return all packages implied by the discovered module hierarchy."""

    packages = {PACKAGE_NAME}

    for module_name in modules:

        parts = module_name.split(".")

        # Every prefix except the final component is a package.
        for i in range(2, len(parts)):
            packages.add(
                ".".join(parts[:i])
            )

    return packages


###############################################################################
# OBJECT DISCOVERY
###############################################################################


def get_module_classes(module):
    """Return public classes defined by a module."""

    classes = []

    for name, obj in inspect.getmembers(
        module,
        inspect.isclass,
    ):

        if not is_public_name(name):
            continue

        if obj.__module__ != module.__name__:
            continue

        classes.append((name, obj))

    return classes


###############################################################################


def get_module_functions(module):
    """Return public functions defined by a module."""

    functions = []

    for name, obj in inspect.getmembers(
        module,
        inspect.isfunction,
    ):

        if not is_public_name(name):
            continue

        if obj.__module__ != module.__name__:
            continue

        functions.append((name, obj))

    return functions


###############################################################################


def get_class_methods(cls):
    """Return public methods defined directly by a class."""

    methods = []

    for name, obj in cls.__dict__.items():

        if name == "__init__":
            continue

        if not is_public_name(name):
            continue

        if isinstance(obj, staticmethod):
            obj = obj.__func__

        elif isinstance(obj, classmethod):
            obj = obj.__func__

        if inspect.isfunction(obj):
            methods.append((name, obj))

        elif inspect.ismethoddescriptor(obj):
            methods.append((name, obj))

        elif isinstance(obj, property):
            methods.append((name, obj))

    return methods


###############################################################################
# DOCSTRING FORMATTING
###############################################################################


def format_docstring(docstring: str) -> str:
    """Convert a plain docstring into safe HTML."""

    if not docstring:

        return (
            "<p><em>"
            "No description available."
            "</em></p>"
        )

    return (
        '<div class="docstring">'
        + escape(docstring)
        + "</div>"
    )


###############################################################################
# FUNCTIONS
###############################################################################


def format_function(name: str, func) -> str:
    """Generate HTML documentation for a function."""

    signature = get_signature(func)
    docstring = get_docstring(func)

    return f"""
<div class="member">

<h3>{escape(name)}</h3>

<div class="signature">

<code>
{escape(name)}{escape(signature)}
</code>

</div>

{format_docstring(docstring)}

</div>
"""


###############################################################################
# METHODS
###############################################################################


def format_method(name: str, method) -> str:
    """Generate HTML documentation for a method or property."""

    if isinstance(method, property):

        signature = ""
        docstring = get_docstring(method)

    else:

        signature = get_signature(method)
        docstring = get_docstring(method)

    signature_html = ""

    if signature:

        signature_html = f"""
<div class="signature">

<code>
{escape(name)}{escape(signature)}
</code>

</div>
"""

    return f"""
<div class="member">

<h4>{escape(name)}</h4>

{signature_html}

{format_docstring(docstring)}

</div>
"""


###############################################################################
# CLASSES
###############################################################################


def format_class(name: str, cls) -> str:
    """Generate HTML documentation for a class."""

    signature = get_signature(cls)
    docstring = get_docstring(cls)

    bases = [
        base.__name__
        for base in cls.__bases__
        if base is not object
    ]

    inheritance_html = ""

    if bases:

        inheritance_html = (
            '<div class="inheritance">'
            "Inherits from: "
            + ", ".join(
                escape(base)
                for base in bases
            )
            + "</div>"
        )

    methods = get_class_methods(cls)

    methods_html = ""

    if methods:

        methods_html += "<h3>Methods</h3>"

        for method_name, method in methods:

            methods_html += format_method(
                method_name,
                method,
            )

    return f"""
<section>

<h2>{escape(name)}</h2>

<div class="signature">

<code>
{escape(name)}{escape(signature)}
</code>

</div>

{inheritance_html}

{format_docstring(docstring)}

{methods_html}

</section>
"""


###############################################################################
# BREADCRUMBS
###############################################################################


def make_breadcrumbs(
    current_file: Path,
    module_name: str,
) -> str:
    """Create breadcrumb navigation from a module/package name."""

    parts = module_name.split(".")

    links = []

    root_file = HTML_DIR / "index.html"

    links.append(
        f'<a href="{relative_link(current_file, root_file)}">'
        "FinancePy"
        "</a>"
    )

    if parts and parts[0] == PACKAGE_NAME:
        parts = parts[1:]

    package_parts = []

    for part in parts[:-1]:

        package_parts.append(part)

        package_name = (
            PACKAGE_NAME
            + "."
            + ".".join(package_parts)
        )

        package_file = package_to_output_path(
            package_name
        )

        links.append(
            f'<a href="'
            f'{relative_link(current_file, package_file)}'
            f'">{escape(display_name(part))}</a>'
        )

    if parts:

        links.append(
            escape(display_name(parts[-1]))
        )

    return (
        "<nav>"
        + " &rsaquo; ".join(links)
        + "</nav>"
    )


###############################################################################
# MODULE PAGES
###############################################################################


def generate_module_page(
    module_name: str,
) -> Path | None:
    """Generate documentation for one Python module."""

    try:

        module = importlib.import_module(
            module_name
        )

    except Exception as exc:

        print(
            f"WARNING: Could not import "
            f"{module_name}: "
            f"{type(exc).__name__}: {exc}"
        )

        return None

    classes = get_module_classes(module)
    functions = get_module_functions(module)

    # Package modules themselves usually have no useful
    # classes/functions. Their index pages are generated separately.
    if not classes and not functions:
        return None

    output_file = module_to_output_path(
        module_name
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    module_doc = get_docstring(module)

    title = display_name(
        module_name.split(".")[-1]
    )

    body = page_header(
        title,
        "FinancePy API Reference",
    )

    body += make_breadcrumbs(
        output_file,
        module_name,
    )

    body += f"""
<h1>{escape(title)}</h1>

<p class="source-name">
{escape(module_name)}
</p>
"""

    if module_doc:

        body += format_docstring(
            module_doc
        )

    if classes:

        body += "<h1>Classes</h1>"

        for class_name, cls in classes:

            body += format_class(
                class_name,
                cls,
            )

    if functions:

        body += "<h1>Functions</h1>"

        for function_name, func in functions:

            body += format_function(
                function_name,
                func,
            )

    body += page_footer()

    output_file.write_text(
        body,
        encoding="utf-8",
    )

    return output_file


###############################################################################
# PACKAGE RELATIONSHIPS
###############################################################################


def direct_child_packages(
    package_name: str,
    packages: set[str],
) -> list[str]:
    """Return immediate child packages of a package."""

    prefix = package_name + "."

    depth = package_name.count(".") + 1

    children = []

    for candidate in packages:

        if not candidate.startswith(prefix):
            continue

        if candidate.count(".") != depth:
            continue

        children.append(candidate)

    return sorted(children)


###############################################################################


def direct_child_modules(
    package_name: str,
    generated_files: dict[str, Path],
) -> list[str]:
    """Return immediate child modules of a package."""

    prefix = package_name + "."

    depth = package_name.count(".") + 1

    children = []

    for module_name in generated_files:

        if not module_name.startswith(prefix):
            continue

        if module_name.count(".") != depth:
            continue

        children.append(module_name)

    return sorted(children)


###############################################################################
# PACKAGE PAGES
###############################################################################


def generate_package_page(
    package_name: str,
    packages: set[str],
    generated_files: dict[str, Path],
):
    """Generate an index page for a package."""

    output_file = package_to_output_path(
        package_name
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    package_short_name = (
        package_name.split(".")[-1]
    )

    title = display_name(
        package_short_name
    )

    body = page_header(
        title,
        "FinancePy API Reference",
    )

    if package_name != PACKAGE_NAME:

        body += make_breadcrumbs(
            output_file,
            package_name,
        )

    body += f"""
<h1>{escape(title)}</h1>

<p class="source-name">
{escape(package_name)}
</p>
"""

    child_packages = direct_child_packages(
        package_name,
        packages,
    )

    child_modules = direct_child_modules(
        package_name,
        generated_files,
    )

    if child_packages:

        body += """
<h2>Packages</h2>

<div class="card-grid">
"""

        for child_package in child_packages:

            child_name = (
                child_package.split(".")[-1]
            )

            child_file = (
                package_to_output_path(
                    child_package
                )
            )

            link = relative_link(
                output_file,
                child_file,
            )

            body += f"""
<div class="card">

<h3>
<a href="{link}">
{escape(display_name(child_name))}
</a>
</h3>

<p>
{escape(child_package)}
</p>

</div>
"""

        body += "</div>"

    if child_modules:

        body += """
<h2>Modules</h2>

<div class="module-list">

<ul>
"""

        for module_name in child_modules:

            module_short_name = (
                module_name.split(".")[-1]
            )

            module_file = (
                generated_files[module_name]
            )

            link = relative_link(
                output_file,
                module_file,
            )

            body += f"""
<li>
<a href="{link}">
{escape(display_name(module_short_name))}
</a>
</li>
"""

        body += """
</ul>

</div>
"""

    if not child_packages and not child_modules:

        body += """
<p>
No public modules were found in this package.
</p>
"""

    body += page_footer()

    output_file.write_text(
        body,
        encoding="utf-8",
    )


###############################################################################
# ROOT INDEX
###############################################################################


def generate_root_index(
    packages: set[str],
    generated_files: dict[str, Path],
):
    """Generate the main FinancePy API index."""

    output_file = HTML_DIR / "index.html"

    body = page_header(
        "FinancePy API Reference",
        "Automatically generated from the FinancePy source code.",
    )

    body += """
<h1>FinancePy API Reference</h1>

<p>
This reference is generated automatically from the FinancePy source
code, function signatures, type annotations and docstrings.
</p>

<p>
For worked examples showing how to construct and value FinancePy
products, see the <code>examples</code> directory in the FinancePy
repository.
</p>
"""

    child_packages = direct_child_packages(
        PACKAGE_NAME,
        packages,
    )

    if child_packages:

        body += """
<h2>Contents</h2>

<div class="card-grid">
"""

        for child_package in child_packages:

            short_name = (
                child_package.split(".")[-1]
            )

            child_file = (
                package_to_output_path(
                    child_package
                )
            )

            link = relative_link(
                output_file,
                child_file,
            )

            body += f"""
<div class="card">

<h3>
<a href="{link}">
{escape(display_name(short_name))}
</a>
</h3>

<p>
{escape(child_package)}
</p>

</div>
"""

        body += "</div>"

    body += page_footer()

    output_file.write_text(
        body,
        encoding="utf-8",
    )


###############################################################################
# OUTPUT
###############################################################################


def clean_output(force: bool = False):
    """Prepare the documentation output directory."""

    if force and HTML_DIR.exists():

        try:
            shutil.rmtree(HTML_DIR)

        except PermissionError as exc:
            print()
            print("WARNING: Could not completely clean HTML directory.")
            print("A file may be open in a browser, Explorer, or Dropbox.")
            print()
            print(f"    {exc}")
            print()
            print("Continuing without a full clean.")
            print()

    HTML_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )


###############################################################################
# BUILD
###############################################################################


def build():
    """Generate the complete FinancePy HTML API documentation."""

    project_root = str(PROJECT_ROOT)

    if project_root not in sys.path:

        sys.path.insert(
            0,
            project_root,
        )

    clean_output()

    print(
        "Discovering FinancePy modules..."
    )

    modules = discover_modules()

    print(
        f"Found {len(modules)} modules."
    )

    packages = discover_packages(
        modules
    )

    generated_files: dict[str, Path] = {}

    for module_name in modules:

        print(
            f"Documenting {module_name}"
        )

        output_file = generate_module_page(
            module_name
        )

        if output_file is not None:

            generated_files[module_name] = (
                output_file
            )

    print()
    print("Generating package indexes...")

    # Generate deepest packages first.
    sorted_packages = sorted(
        packages,
        key=lambda x: (
            -x.count("."),
            x,
        ),
    )

    for package_name in sorted_packages:

        if package_name == PACKAGE_NAME:
            continue

        print(
            f"Indexing {package_name}"
        )

        generate_package_page(
            package_name,
            packages,
            generated_files,
        )

    generate_root_index(
        packages,
        generated_files,
    )

    print()
    print("=" * 72)
    print(
        "FinancePy API documentation generated"
    )
    print("=" * 72)
    print()
    print("Output directory:")
    print(f"    {HTML_DIR}")
    print()
    print("Open:")
    print(
        f"    {HTML_DIR / 'index.html'}"
    )
    print()


###############################################################################


if __name__ == "__main__":
    build()