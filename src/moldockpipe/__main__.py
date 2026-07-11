from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from .project import ProjectRepository


def configure_qt_runtime() -> None:
    """Prefer PyQt's bundled Qt files when Conda also exposes Qt DLLs."""
    try:
        import PyQt6
    except ImportError:
        return
    qt_root = Path(PyQt6.__file__).resolve().parent / "Qt6"
    plugins = qt_root / "plugins"
    binaries = qt_root / "bin"
    if plugins.is_dir():
        os.environ["QT_PLUGIN_PATH"] = str(plugins)
    if binaries.is_dir():
        if hasattr(os, "add_dll_directory"):
            os.add_dll_directory(str(binaries))
        os.environ["PATH"] = str(binaries) + os.pathsep + os.environ.get("PATH", "")


def main() -> int:
    parser = argparse.ArgumentParser(prog="moldockpipe")
    subparsers = parser.add_subparsers(dest="command")
    create = subparsers.add_parser("create", help="Create a portable project")
    create.add_argument("path", type=Path)
    subparsers.add_parser("ui", help="Open the desktop UI")
    args = parser.parse_args()
    if args.command == "create":
        repo = ProjectRepository.create(args.path)
        print(repo.root)
        return 0
    try:
        configure_qt_runtime()
        from PyQt6.QtWidgets import QApplication
        from .ui.main_window import MainWindow
    except ImportError as exc:
        parser.error(f"PyQt6 is required for the desktop application: {exc}")
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
