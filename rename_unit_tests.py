
import re
import subprocess
import sys
import uuid
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent / "tests" / "unit"


def to_snake_case(name):
    # Remove Fin prefix, but only when followed by an uppercase letter.
    name = re.sub(r"^Fin(?=[A-Z])", "", name)

    # Separate acronyms from following words.
    name = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", name)

    # Separate lowercase/digit from uppercase.
    name = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", name)

    return name.lower()


def new_filename(filename):
    if not filename.startswith("test_") or not filename.endswith(".py"):
        return filename

    body = filename[5:-3]
    return "test_" + to_snake_case(body) + ".py"


def main():
    apply = "--apply" in sys.argv

    if not TEST_DIR.is_dir():
        raise SystemExit(f"Directory not found: {TEST_DIR}")

    changes = []
    for path in sorted(TEST_DIR.glob("test_*.py")):
        target = path.with_name(new_filename(path.name))
        if path.name != target.name:
            changes.append((path, target))

    # Check case-insensitive filename conflicts (important on Windows).
    final_names = {}
    for path in TEST_DIR.iterdir():
        target = next(
            (dst for src, dst in changes if src == path),
            path,
        )
        key = target.name.casefold()
        if key in final_names:
            raise SystemExit(
                f"Filename collision: {final_names[key]} and {target.name}"
            )
        final_names[key] = target.name

    for source, target in changes:
        print(f"{source.name} -> {target.name}")

    print(f"\n{len(changes)} files to rename.")

    if not apply:
        print("Preview only. Run with --apply to rename.")
        return

    # Two-stage rename supports case-only changes on Windows.
    staged = []
    for source, target in changes:
        temporary = source.with_name(
            f"__rename_{uuid.uuid4().hex}.py"
        )
        subprocess.run(
            ["git", "mv", str(source), str(temporary)],
            check=True,
        )
        staged.append((temporary, target))

    for temporary, target in staged:
        subprocess.run(
            ["git", "mv", str(temporary), str(target)],
            check=True,
        )

    print("\nRenaming complete.")


if __name__ == "__main__":
    main()
