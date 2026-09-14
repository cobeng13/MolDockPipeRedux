from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from .project import ProjectRepository


_DLL_DIRECTORY_HANDLES: list[object] = []


def configure_qt_runtime() -> None:
    """Prefer PyQt's bundled Qt files when Conda also exposes Qt DLLs."""
    try:
        import PyQt6
    except ImportError:
        return
    pyqt_root = Path(PyQt6.__file__).resolve().parent / "Qt6"
    conda_root = Path(sys.prefix) / "Library"
    qt_root = next((root for root in (pyqt_root, conda_root)
                    if (root / "plugins").is_dir() or (root / "bin").is_dir()), pyqt_root)
    plugins = qt_root / "plugins"
    binaries = qt_root / "bin"
    if plugins.is_dir():
        os.environ["QT_PLUGIN_PATH"] = str(plugins)
    if binaries.is_dir():
        if hasattr(os, "add_dll_directory") and not _DLL_DIRECTORY_HANDLES:
            # The returned object owns the Windows DLL-directory registration.
            # Retain it for the process lifetime; dropping it immediately can
            # remove the path again and cause native Qt load/unload failures.
            _DLL_DIRECTORY_HANDLES.append(os.add_dll_directory(str(binaries)))
        os.environ["PATH"] = str(binaries) + os.pathsep + os.environ.get("PATH", "")


def main() -> int:
    parser = argparse.ArgumentParser(prog="moldockpipe")
    subparsers = parser.add_subparsers(dest="command")
    create = subparsers.add_parser("create", help="Create a portable project")
    create.add_argument("path", type=Path)
    run = subparsers.add_parser("run", help="Run pipeline stages without the desktop UI")
    run.add_argument("project", type=Path)
    run.add_argument(
        "--stages", nargs="+", choices=("screening", "molscrub", "meeko", "vina", "postdock"),
        default=("screening", "molscrub", "meeko", "vina", "postdock"),
        help="Stages to run in order (default: the complete pipeline)",
    )
    subparsers.add_parser("ui", help="Open the desktop UI")
    args = parser.parse_args()
    if args.command == "create":
        repo = ProjectRepository.create(args.path)
        print(repo.root)
        return 0
    if args.command == "run":
        if not (args.project / "project.sqlite").is_file():
            parser.error(f"Not a MolDockPipe project: {args.project}")
        from .pipeline import PipelineRunner

        repo = ProjectRepository(args.project)
        repo.recover_interrupted_runs()
        runner = PipelineRunner(repo)
        methods = {
            "screening": runner.run_screening,
            "molscrub": runner.run_molscrub,
            "meeko": runner.run_meeko,
            "vina": runner.run_vina,
            "postdock": runner.run_postdock,
        }
        summary: dict[str, object] = {}
        failed = 0
        for stage in args.stages:
            result = methods[stage]()
            summary[stage] = result
            if isinstance(result, tuple) and len(result) > 1:
                failed += int(result[1])
        print(json.dumps(summary, sort_keys=True))
        return 1 if failed else 0
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
