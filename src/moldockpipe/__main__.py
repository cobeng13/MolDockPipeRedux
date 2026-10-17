from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from .project import ProjectRepository
from .execution_lock import project_execution_lock


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


def _run_project(project: Path, stages: list[str]) -> int:
    from .pipeline import PipelineRunner

    # Protect import and crash recovery too: neither may run concurrently with
    # a GUI worker that is still writing this project's state and artifacts.
    with project_execution_lock(project):
        repo = ProjectRepository(project)
        repo.migrate()
        runner = PipelineRunner(repo)
        if "vina" in stages:
            runner.validate_docking_configuration()
        repo.import_project_input()
        repo.recover_interrupted_runs()
        summary = {}
        failed = 0
        for stage in stages:
            result = getattr(runner, f"run_{stage}")()
            summary[stage] = result
            if isinstance(result, tuple) and len(result) > 1:
                failed += int(result[1])
        print(json.dumps(summary, sort_keys=True))
        return 1 if failed else 0


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
    doctor = subparsers.add_parser("doctor", help="Test chemistry packages and docking engines")
    doctor.add_argument("--project", type=Path, default=Path.cwd(), help="Project root for executable discovery")
    doctor.add_argument("--json", action="store_true", help="Print the report as JSON")
    doctor.add_argument("--output", type=Path, help="Also save the report to this file")
    args = parser.parse_args()
    if args.command == "doctor":
        from .diagnostics import format_report, run_diagnostics

        report = run_diagnostics(args.project)
        rendered = json.dumps(report, indent=2) if args.json else format_report(report)
        if args.output:
            args.output.write_text(rendered, encoding="utf-8")
        print(rendered)
        return int(any(check["status"] == "FAIL" for check in report["checks"]))
    if args.command == "create":
        repo = ProjectRepository.create(args.path)
        print(repo.root)
        return 0
    if args.command == "run":
        if not (args.project / "project.sqlite").is_file():
            parser.error(f"Not a MolDockPipe project: {args.project}")
        return _run_project(args.project, args.stages)
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
