import ast
from pathlib import Path


ROOT = Path("financepy")


def is_public_method(node: ast.FunctionDef) -> bool:
    if node.name == "__init__":
        return True

    return not node.name.startswith("_")


for filename in ROOT.rglob("*.py"):

    source = filename.read_text(encoding="utf-8")
    tree = ast.parse(source)

    for node in ast.walk(tree):

        if not isinstance(node, ast.ClassDef):
            continue

        for item in node.body:

            if not isinstance(
                item,
                (ast.FunctionDef, ast.AsyncFunctionDef),
            ):
                continue

            if not is_public_method(item):
                continue

            missing = []

            arguments = (
                item.args.posonlyargs
                + item.args.args
                + item.args.kwonlyargs
            )

            for arg in arguments:

                if arg.arg in ("self", "cls"):
                    continue

                if arg.annotation is None:
                    missing.append(arg.arg)

            if item.args.vararg is not None:
                if item.args.vararg.annotation is None:
                    missing.append(
                        "*" + item.args.vararg.arg
                    )

            if item.args.kwarg is not None:
                if item.args.kwarg.annotation is None:
                    missing.append(
                        "**" + item.args.kwarg.arg
                    )

            if item.returns is None:
                missing.append("return")

            if missing:
                print(
                    f"{filename}: "
                    f"{node.name}.{item.name}: "
                    f"{', '.join(missing)}"
                )
