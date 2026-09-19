"""Keep the console open when a FinancePy example exits.

Every example imports and installs this handler near startup.  The pause is
registered with :mod:`atexit`, so it also runs after most unhandled Python
exceptions and lets Windows users read the traceback before closing the
console window.
"""
from __future__ import annotations

import atexit


def install_double_click_pause() -> None:
    """Pause for Enter when the Python interpreter is about to exit."""

    def _pause() -> None:
        try:
            input("\nPress Enter to close this window...")
        except (EOFError, KeyboardInterrupt):
            pass

    atexit.register(_pause)
