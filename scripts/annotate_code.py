##############################################################################
# FinancePy safe annotation utility
##############################################################################
#
# Adds high-confidence type annotations using FinancePy naming conventions.
#
# Rules, in precedence order:
#
#     *_dt                         -> Date
#     *_curves                     -> list[DiscountCurve]
#     *_curve                      -> DiscountCurve
#     model                        -> Model
#
#     corr_matrix                  -> np.ndarray
#     *_times                      -> np.ndarray
#     *vector                      -> np.ndarray
#     contains "volatilities"      -> np.ndarray
#     contains "correlations"      -> np.ndarray
#
#     contains "num"               -> int
#     contains "node_index"        -> int
#
#     contains "rate"              -> float
#     contains "price"             -> float
#     contains "bump"              -> float
#     contains "recovery"          -> float
#     contains "sigma"             -> float
#     contains "volatility"        -> float
#     contains "cpn"               -> float
#
#     __init__ return              -> None
#
# Existing annotations are NEVER changed.
# Imports are NEVER changed.
# Function bodies are NEVER deliberately changed.
# Function calls are NEVER deliberately changed.
# Original source files are NEVER overwritten.
#
# CRLF/LF line endings are preserved.
#
# Modified files are written into an "annotated" subfolder.
#
# Usage:
#
#     python annotate_code.py . --dry-run
#     python annotate_code.py .
#
##############################################################################

import argparse
import ast
from pathlib import Path

import libcst as cst


# ============================================================================
# TYPE ANNOTATION TRANSFORMER
# ============================================================================


class TypeAnnotationTransformer(cst.CSTTransformer):
    """Add high-confidence FinancePy type annotations."""

    def __init__(self) -> None:

        self.annotations_added = 0

        self.date_annotations_added = 0
        self.curve_annotations_added = 0
        self.curve_list_annotations_added = 0
        self.model_annotations_added = 0

        self.int_annotations_added = 0
        self.float_annotations_added = 0
        self.ndarray_annotations_added = 0

        self.return_annotations_added = 0

    # ------------------------------------------------------------------------
    # PARAMETER ANNOTATION
    # ------------------------------------------------------------------------

    def _annotate_param(
        self,
        param: cst.Param,
    ) -> cst.Param:
        """Annotate a parameter using FinancePy naming conventions."""

        # Never overwrite an existing annotation.
        if param.annotation is not None:
            return param

        name = param.name.value

        # --------------------------------------------------------------------
        # Date
        # --------------------------------------------------------------------

        if name.endswith("_dt"):

            annotation = "Date"
            self.date_annotations_added += 1

        # --------------------------------------------------------------------
        # List of discount curves
        #
        # Must come before the singular curve rule.
        # --------------------------------------------------------------------

        elif name.endswith("_curves"):

            annotation = "list[DiscountCurve]"
            self.curve_list_annotations_added += 1

        # --------------------------------------------------------------------
        # Single discount curve
        # --------------------------------------------------------------------

        elif name.endswith("_curve"):

            annotation = "DiscountCurve"
            self.curve_annotations_added += 1

        # --------------------------------------------------------------------
        # Model
        # --------------------------------------------------------------------

        elif name == "model":

            annotation = "Model"
            self.model_annotations_added += 1

        # --------------------------------------------------------------------
        # NumPy arrays
        #
        # These rules must come before the general float rules.
        #
        # In particular:
        #
        #     volatility   -> float
        #     volatilities -> np.ndarray
        #
        # --------------------------------------------------------------------

        elif (
            name == "corr_matrix"
            or name.endswith("_times")
            or name.endswith("vector")
            or "volatilities" in name
            or "correlations" in name
        ):

            annotation = "np.ndarray"
            self.ndarray_annotations_added += 1

        # --------------------------------------------------------------------
        # Integer parameters
        # --------------------------------------------------------------------

        elif (
            "num" in name
            or "node_index" in name
        ):

            annotation = "int"
            self.int_annotations_added += 1

        # --------------------------------------------------------------------
        # Float parameters
        # --------------------------------------------------------------------

        elif (
            "rate" in name
            or "price" in name
            or "bump" in name
            or "recovery" in name
            or "sigma" in name
            or "volatility" in name
            or "cpn" in name
        ):

            annotation = "float"
            self.float_annotations_added += 1

        # --------------------------------------------------------------------
        # No known convention
        # --------------------------------------------------------------------

        else:

            return param

        # --------------------------------------------------------------------
        # Add annotation
        # --------------------------------------------------------------------

        self.annotations_added += 1

        return param.with_changes(
            annotation=cst.Annotation(
                cst.parse_expression(annotation)
            )
        )

    # ------------------------------------------------------------------------
    # FUNCTION DEFINITIONS
    # ------------------------------------------------------------------------

    def leave_FunctionDef(
        self,
        original_node: cst.FunctionDef,
        updated_node: cst.FunctionDef,
    ) -> cst.FunctionDef:
        """Modify only parameter and return annotations."""

        params = updated_node.params

        # --------------------------------------------------------------------
        # Annotate positional-only parameters
        # --------------------------------------------------------------------

        new_posonly_params = [
            self._annotate_param(param)
            for param in params.posonly_params
        ]

        # --------------------------------------------------------------------
        # Annotate normal parameters
        # --------------------------------------------------------------------

        new_params_list = [
            self._annotate_param(param)
            for param in params.params
        ]

        # --------------------------------------------------------------------
        # Annotate keyword-only parameters
        # --------------------------------------------------------------------

        new_kwonly_params = [
            self._annotate_param(param)
            for param in params.kwonly_params
        ]

        # --------------------------------------------------------------------
        # Annotate *args if its name happens to match a rule
        # --------------------------------------------------------------------

        new_star_arg = params.star_arg

        if isinstance(
            new_star_arg,
            cst.Param,
        ):

            new_star_arg = self._annotate_param(
                new_star_arg
            )

        # --------------------------------------------------------------------
        # Annotate **kwargs if its name happens to match a rule
        # --------------------------------------------------------------------

        new_star_kwarg = params.star_kwarg

        if new_star_kwarg is not None:

            new_star_kwarg = self._annotate_param(
                new_star_kwarg
            )

        # --------------------------------------------------------------------
        # Rebuild parameter list
        # --------------------------------------------------------------------

        new_params = params.with_changes(
            posonly_params=new_posonly_params,
            params=new_params_list,
            kwonly_params=new_kwonly_params,
            star_arg=new_star_arg,
            star_kwarg=new_star_kwarg,
        )

        # --------------------------------------------------------------------
        # Return annotation
        #
        # Only __init__ is inferred.
        # --------------------------------------------------------------------

        new_return = updated_node.returns

        if (
            updated_node.name.value == "__init__"
            and new_return is None
        ):

            new_return = cst.Annotation(
                cst.Name("None")
            )

            self.return_annotations_added += 1
            self.annotations_added += 1

        # --------------------------------------------------------------------
        # IMPORTANT
        #
        # Only params and returns are changed.
        #
        # The function body is untouched.
        # --------------------------------------------------------------------

        return updated_node.with_changes(
            params=new_params,
            returns=new_return,
        )


# ============================================================================
# AST SAFETY CHECK
# ============================================================================


class AnnotationStripper(ast.NodeTransformer):
    """Remove function annotations from an AST.

    The original and transformed ASTs are compared after annotations have
    been removed.

    If anything other than annotations changed structurally, the output file
    is rejected.
    """

    def visit_FunctionDef(
        self,
        node: ast.FunctionDef,
    ):

        self.generic_visit(node)

        self._strip_arguments(
            node.args
        )

        node.returns = None
        node.type_comment = None

        return node

    def visit_AsyncFunctionDef(
        self,
        node: ast.AsyncFunctionDef,
    ):

        self.generic_visit(node)

        self._strip_arguments(
            node.args
        )

        node.returns = None
        node.type_comment = None

        return node

    @staticmethod
    def _strip_arguments(
        args: ast.arguments,
    ) -> None:
        """Remove annotations from function arguments."""

        all_args = (
            list(args.posonlyargs)
            + list(args.args)
            + list(args.kwonlyargs)
        )

        if args.vararg is not None:

            all_args.append(
                args.vararg
            )

        if args.kwarg is not None:

            all_args.append(
                args.kwarg
            )

        for arg in all_args:

            arg.annotation = None
            arg.type_comment = None


# ============================================================================
# AST HELPERS
# ============================================================================


def stripped_ast(
    source: str,
) -> str:
    """Return AST representation with function annotations removed."""

    tree = ast.parse(
        source
    )

    tree = AnnotationStripper().visit(
        tree
    )

    ast.fix_missing_locations(
        tree
    )

    return ast.dump(
        tree,
        include_attributes=False,
    )


def safety_check(
    original_source: str,
    modified_source: str,
) -> bool:
    """Verify that only function annotations changed structurally."""

    try:

        original_ast = stripped_ast(
            original_source
        )

        modified_ast = stripped_ast(
            modified_source
        )

    except SyntaxError:

        return False

    return original_ast == modified_ast


# ============================================================================
# NEWLINE HANDLING
# ============================================================================


def detect_newline(
    raw_source: bytes,
) -> str:
    """Detect the newline convention used by the source file."""

    if b"\r\n" in raw_source:

        return "\r\n"

    return "\n"


def restore_newlines(
    source: str,
    newline: str,
) -> str:
    """Restore the original newline convention."""

    # ------------------------------------------------------------------------
    # Normalize to LF first
    # ------------------------------------------------------------------------

    source = source.replace(
        "\r\n",
        "\n",
    )

    # ------------------------------------------------------------------------
    # Restore CRLF if required
    # ------------------------------------------------------------------------

    if newline == "\r\n":

        source = source.replace(
            "\n",
            "\r\n",
        )

    return source


# ============================================================================
# PROCESS ONE FILE
# ============================================================================


def process_file(
    source_file: Path,
    output_file: Path,
    dry_run: bool,
) -> tuple[int, int, int, int, int, int, int, int, int]:
    """Process one Python file while preserving line endings."""

    # ------------------------------------------------------------------------
    # Read exact bytes from disk
    # ------------------------------------------------------------------------

    try:

        raw_source = source_file.read_bytes()

        source = raw_source.decode(
            "utf-8"
        )

    except Exception as exc:

        print(
            f"ERROR reading {source_file}: {exc}"
        )

        return (0,) * 9

    # ------------------------------------------------------------------------
    # Detect original newline convention
    # ------------------------------------------------------------------------

    newline = detect_newline(
        raw_source
    )

    # ------------------------------------------------------------------------
    # Parse with LibCST
    # ------------------------------------------------------------------------

    try:

        module = cst.parse_module(
            source
        )

    except Exception as exc:

        print(
            f"ERROR parsing {source_file}: {exc}"
        )

        return (0,) * 9

    # ------------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------------

    transformer = TypeAnnotationTransformer()

    updated_module = module.visit(
        transformer
    )

    # ------------------------------------------------------------------------
    # Collect statistics
    # ------------------------------------------------------------------------

    total = transformer.annotations_added

    dates = transformer.date_annotations_added

    curves = transformer.curve_annotations_added

    curve_lists = (
        transformer.curve_list_annotations_added
    )

    models = transformer.model_annotations_added

    ints = transformer.int_annotations_added

    floats = transformer.float_annotations_added

    arrays = transformer.ndarray_annotations_added

    returns = transformer.return_annotations_added

    # ------------------------------------------------------------------------
    # Nothing changed
    # ------------------------------------------------------------------------

    if total == 0:

        return (0,) * 9

    # ------------------------------------------------------------------------
    # Generate modified source
    # ------------------------------------------------------------------------

    modified_source = updated_module.code

    # ------------------------------------------------------------------------
    # Restore original newline convention
    # ------------------------------------------------------------------------

    modified_source = restore_newlines(
        modified_source,
        newline,
    )

    # ------------------------------------------------------------------------
    # SAFETY CHECK
    #
    # Once annotations are removed, the original and transformed ASTs must
    # be identical.
    # ------------------------------------------------------------------------

    if not safety_check(
        source,
        modified_source,
    ):

        print()
        print("=" * 78)

        print(
            f"SAFETY CHECK FAILED: {source_file}"
        )

        print(
            "File was NOT written."
        )

        print(
            "Something other than function annotations "
            "appears to have changed."
        )

        print("=" * 78)
        print()

        return (0,) * 9

    # ------------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------------

    print(
        f"{source_file}"
        f" -> {output_file}"
        f"  "
        f"({dates} Date, "
        f"{curves} DiscountCurve, "
        f"{curve_lists} list[DiscountCurve], "
        f"{models} Model, "
        f"{ints} int, "
        f"{floats} float, "
        f"{arrays} np.ndarray, "
        f"{returns} return)"
    )

    # ------------------------------------------------------------------------
    # Write output
    #
    # write_bytes() is deliberate.
    #
    # This avoids Python text-mode newline conversion and helps ensure that
    # Git sees only the annotation lines as changed.
    # ------------------------------------------------------------------------

    if not dry_run:

        output_file.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        output_file.write_bytes(
            modified_source.encode(
                "utf-8"
            )
        )

    return (
        total,
        dates,
        curves,
        curve_lists,
        models,
        ints,
        floats,
        arrays,
        returns,
    )


# ============================================================================
# MAIN
# ============================================================================


def main() -> None:

    # ------------------------------------------------------------------------
    # Command-line arguments
    # ------------------------------------------------------------------------

    parser = argparse.ArgumentParser(
        description=(
            "Safely add FinancePy type annotations "
            "while preserving source line endings."
        )
    )

    parser.add_argument(
        "folder",
        type=Path,
        help=(
            "Source folder to process recursively."
        ),
    )

    parser.add_argument(
        "--output",
        default="annotated",
        help=(
            "Output subfolder name "
            "(default: annotated)."
        ),
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Show proposed changes without writing files."
        ),
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------------
    # Source folder
    # ------------------------------------------------------------------------

    source_folder = args.folder.resolve()

    if not source_folder.exists():

        raise SystemExit(
            f"Folder does not exist: "
            f"{source_folder}"
        )

    if not source_folder.is_dir():

        raise SystemExit(
            f"Not a folder: "
            f"{source_folder}"
        )

    # ------------------------------------------------------------------------
    # Output folder
    # ------------------------------------------------------------------------

    output_folder = (
        source_folder
        / args.output
    ).resolve()

    if output_folder == source_folder:

        raise SystemExit(
            "Output folder cannot be the same "
            "as the source folder."
        )

    # ------------------------------------------------------------------------
    # Collect Python files
    # ------------------------------------------------------------------------

    source_files = sorted(
        source_folder.rglob("*.py")
    )

    # ------------------------------------------------------------------------
    # Totals
    # ------------------------------------------------------------------------

    files_changed = 0

    total_annotations = 0

    total_dates = 0
    total_curves = 0
    total_curve_lists = 0
    total_models = 0

    total_ints = 0
    total_floats = 0
    total_arrays = 0

    total_returns = 0

    # ------------------------------------------------------------------------
    # Process files
    # ------------------------------------------------------------------------

    for source_file in source_files:

        # --------------------------------------------------------------------
        # Never process generated output
        # --------------------------------------------------------------------

        if output_folder in source_file.parents:

            continue

        # --------------------------------------------------------------------
        # Never process this script
        # --------------------------------------------------------------------

        if (
            source_file.resolve()
            == Path(__file__).resolve()
        ):

            continue

        # --------------------------------------------------------------------
        # Preserve directory structure
        # --------------------------------------------------------------------

        relative_path = (
            source_file.relative_to(
                source_folder
            )
        )

        output_file = (
            output_folder
            / relative_path
        )

        # --------------------------------------------------------------------
        # Process file
        # --------------------------------------------------------------------

        (
            total,
            dates,
            curves,
            curve_lists,
            models,
            ints,
            floats,
            arrays,
            returns,
        ) = process_file(
            source_file,
            output_file,
            args.dry_run,
        )

        # --------------------------------------------------------------------
        # Accumulate totals
        # --------------------------------------------------------------------

        if total > 0:

            files_changed += 1

            total_annotations += total

            total_dates += dates

            total_curves += curves

            total_curve_lists += curve_lists

            total_models += models

            total_ints += ints

            total_floats += floats

            total_arrays += arrays

            total_returns += returns

    # ------------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)

    if args.dry_run:

        print(
            "DRY RUN - no files were written"
        )

    else:

        print(
            f"Output folder               : "
            f"{output_folder}"
        )

    print(
        f"Files generated             : "
        f"{files_changed}"
    )

    print(
        f"Date annotations            : "
        f"{total_dates}"
    )

    print(
        f"DiscountCurve annotations   : "
        f"{total_curves}"
    )

    print(
        f"Curve-list annotations      : "
        f"{total_curve_lists}"
    )

    print(
        f"Model annotations           : "
        f"{total_models}"
    )

    print(
        f"Integer annotations         : "
        f"{total_ints}"
    )

    print(
        f"Float annotations           : "
        f"{total_floats}"
    )

    print(
        f"np.ndarray annotations      : "
        f"{total_arrays}"
    )

    print(
        f"Return annotations          : "
        f"{total_returns}"
    )

    print(
        f"Total annotations           : "
        f"{total_annotations}"
    )

    print("=" * 78)


# ============================================================================
# ENTRY POINT
# ============================================================================


if __name__ == "__main__":

    main()
