# Large ligand campaign freeze investigation — 2026-10-24

The investigation below records the baseline at `3403bce`; its source line
references refer to that commit. The fixes described at the end are now applied
to the working tree.

The strongest evidence points to dashboard work blocking the Qt GUI thread. A
synthetic campaign at 20,000 parents and 32 states per parent stalls Qt for about
2.9 seconds per live dashboard refresh, even with the current indexes. This
matches the architecture behind the report that Vina kept running: docking uses
a background QThread and an external subprocess; the dashboard slot runs on the
GUI thread. This does not establish the cause on either user's machine without
their version, database, and resource measurements.

## Measured behavior

Benchmarked current source at commit `3403bce`, Windows, using the
`moldockpipe-clean` Python environment, offscreen Qt, and temporary synthetic
SQLite projects. One enabled receptor, one conformer per state, nine poses per
completed docking run; all screening and preparation succeeded. These are
measurements of dashboard refreshes, not a full chemistry/Vina workload.

| Parent ligands | States per parent | Pose records | Median refresh | Maximum Qt heartbeat gap | Private memory before → after |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 1 | 9,000 | 22 ms | 37 ms | 31.62 → 32.59 MiB |
| 5,000 | 1 | 45,000 | 48 ms | 62 ms | 32.29 → 32.13 MiB |
| 20,000 | 1 | 180,000 | 150 ms | 161 ms | 32.62 → 33.39 MiB |
| 20,000 | 8 | 1,440,000 | 712 ms | 748 ms | 34.32 → 33.85 MiB |
| 20,000 | 32 | 5,760,000 | 2,872 ms | 2,921 ms | 32.11 → 33.18 MiB |

First four cases: 15 repeated refreshes; last case: five. The timer is a 10 ms Qt
timer, with events serviced between synchronous refreshes. Its gaps show the
direct blockage from a refresh, not end-to-end responsiveness during a concurrent
Vina run. The last case ran in a separate process. Small memory changes do not
prove absence of a long-running leak. Temporary databases were removed by the
benchmark; existing user projects were not migrated or changed.

Raw results are saved beside this report in `dashboard-scale-results.json`.
Reproduce using `tools/benchmark_dashboard.py` with the application's Python,
`PYTHONPATH=src`, and `QT_QPA_PLATFORM=offscreen`:

```text
python tools/benchmark_dashboard.py --output results.json --repeats 15
python tools/benchmark_dashboard.py --output results-32.json --repeats 5 --case 20000 32
```

## Findings and mechanisms

1. **Live aggregation blocks the GUI; refresh throttling can collapse.**
   `ui/main_window.py:748` handles queued progress on the GUI thread and invokes
   `_refresh_live_dashboard` (`:600`), which executes project-wide counts and
   calls `_update_dashboard` (`:482`) for further queries and failure lists.
   At `:802`, `_last_live_refresh_at` is set *before* executing the refresh.
   If queries take longer than the one-second interval, the next queued event
   can immediately perform another full refresh. Routine progress can arrive
   up to ten times per second. This can make the window appear continuously
   frozen and delay delivery of Stop clicks while the pipeline keeps working.
   The 32-state result demonstrates the condition, rather than proving that the
   affected users generated that many states. Slower hardware, additional
   receptors, CPU contention, or storage latency can lower the threshold.

2. **An old database without the history index has a severe scaling hazard.**
   `ProjectRepository.molscrub_failures` (`project.py:115`) performs a correlated
   latest-history lookup for every parent. On isolated benchmark databases after
   dropping `idx_stage_runs_entity`, that lookup took 61 ms at 1k parents and
   1.49 seconds at 5k; at 20k it was still running when interrupted at three
   seconds. Current migration already creates this index (`project.py:261`).
   Therefore this is a version/migration check, not an unimplemented fix in this
   checkout. Several local archived projects lack it, but none of those has a
   large ligand library. The largest local project has 5,043 parents and 2,047
   states and has the index. There is no local 20k end-user database to inspect.

3. **Failure-related progress can bypass the routine rate limit.**
   `ui/workers.py:43` always delivers `item_failed`, and skipped-item events with
   a nonzero `failed` value. That field is cumulative. In a controlled probe with
   monotonic time held constant, 20,000 skipped events with `failed=1` emitted
   20,000 signals; the same events with `failed=0` emitted two (first and final).
   Known-failure resumes intentionally populate the cumulative count. Large
   failure batches can flood the Qt queue; UI failure logging adds work even
   though the log document itself is capped at 1,000 blocks. This probe measures
   signal selection, not actual queue memory consumption.

4. **Memory pressure remains plausible during preparation, but was not observed
   as a dashboard leak.** `parallel.py:16` defaults preparation workers to all
   logical CPUs and has no available-memory cap. MolScrub enumeration starts a
   separate Python/native chemistry process per active job; returned native
   molecules and Meeko work also consume memory. In-flight futures are bounded
   to the worker count, so the code does not launch 20k chemistry jobs at once.
   Pipeline candidate rows are still fetched into whole-project lists, and
   state expansion increases their size. Vina CPU docking runs sequentially per
   state/receptor, with one external Vina process at a time and CPU count
   defaulting to all logical CPUs. No complete process-tree memory benchmark or
   live measurements from the affected users were available. The small process
   sizes in the table must not be treated as full-campaign RAM requirements.

5. **Loading and Results have additional synchronous work.** CSV preview/import,
   project activation, idle checkpoint refreshes, and export completeness scans
   occur on the GUI thread. `ui/results_dialog.py:162` ranks all poses and builds
   every parent row, eagerly for every receptor tab. The compound selector also
   creates all rows and checks export paths synchronously. These are additional
   potential freezes if users open those views during docking; they were not
   timed in this investigation. Main dashboard data tables are hidden and
   skipped, live export scans are already disabled, normal progress is already
   throttled, and history/state/conformer/pose indexes already exist in current
   source and `build/lib`. Those existing protections are not new findings to fix.

## Recommended implementation order

1. Move dashboard database queries and filesystem checks into a dedicated
   background snapshot worker. Keep at most one request in flight, use a
   separate SQLite connection there, coalesce pending updates, and apply only
   the latest snapshot on the GUI thread. Record completion time for scheduling.
   Use progress counters for cheap live labels rather than querying per event.
2. Preserve every failure in durable pipeline history while batching/coalescing
   GUI notifications, distinguishing new failures from cumulative totals.
3. Page or lazily load Results and compound tables; move CSV parsing and export
   completeness checks off the GUI thread.
4. Add a memory-aware preparation concurrency cap and optional CPU headroom.
5. Add opt-in diagnostic capture: version and executable path, database index
   presence, parent/state/conformer/run/pose counts, dashboard duration, GUI
   heartbeat lag, app private bytes, child-process memory, system available RAM
   and commit usage, and CPU utilization. This can distinguish a blocked GUI
   from system-wide paging on the affected machines.

## Validation and scope

Existing targeted tests: **39 passed** across `test_resume_ui_handoff.py`,
`test_parallel.py`, and `test_bulk_planning.py`. The first attempt lacked the
Conda `Library/bin` DLL path and produced RDKit import errors; rerunning with
that path configured passed all 39. No production behavior was changed.

Next evidence needed from affected users: installed build/date and Python path,
RAM, generated-state count, number of receptors, exact phase/view at freeze,
whether output files continue advancing, and ideally an anonymized database
and process-tree resource capture. No hard 20,000-ligand limit was found; work
scales with states, receptor pairs, poses, history, and storage speed.

## Applied fixes and subsequent measurements

The dashboard now collects a cancellable snapshot in a background QThread,
keeps one request in flight, rejects results for another project or running
state, and waits one second after completion before collecting again. Idle
export checks stream records rather than materializing every export row.
Project creation/opening/migration, CSV import, and input refresh also run in
background tasks. Closing waits for outstanding tasks without blocking Qt.

Routine failure/skip notifications share the item rate limit; detailed failures
remain persisted by the pipeline and final stage totals still arrive. Results
load only the visible receptor in 200-row pages with cancellable SQL and safe
thread shutdown. Compound selection uses a QTableView model with virtual cells;
leaderboard aggregation and export-path checks run in the background.

Automatic preparation uses at most eight workers, leaves one logical CPU free,
and caps concurrency against available physical/commit memory. All requested
worker counts share a conservative 512 MiB/job budget, reserving at least 2 GiB
or 25% of available memory. Unknown memory availability falls back to two.
This estimate reduces risk; it cannot guarantee a bound on individual native
chemistry processes. Existing docking CPU choices are preserved.

**Tools → Save performance report** writes an opt-in local JSON snapshot with
version/source hashes, counts, indexes, dashboard duration, GUI heartbeat gaps,
available memory, and Windows app/child process memory and cumulative CPU times.
It contains local paths and does not include ligand structures.

The post-fix benchmark used 20,000 parents, 32 states per parent, 5,760,000 poses,
and three snapshots. Maximum GUI request time was **3.1 ms** and maximum Qt
heartbeat gap was **18.4 ms**, compared with **2,921 ms** at baseline. Snapshot
collection itself still took 5.04 seconds median in the background; the fix
removes that work from the GUI rather than promising faster SQL. Private memory
was 33.22 → 33.21 MiB over the measured snapshots. Other tests were running on
the host during this follow-up benchmark, and no Vina/chemistry was launched
by the benchmark. Raw follow-up data is in `dashboard-scale-after-fix.json`.

Final validation: **288 passed, 3 skipped** in the full suite; the final focused
GUI run also passed **35 tests, 2 skipped**, including Stop during a slow
snapshot, coalescing 20,000 requests, stale-result rejection, closing during
table SQL, bounded result pages, and a 20k-row virtual compound selector.
The skips cover optional native WebEngine checks and the existing interactive
window smoke test. Initial combined runs exposed a Qt test-lifecycle issue:
selector tests recreated QApplication and old dialogs could outlive their
workers. Tests now share a session application and dispose dialogs on the GUI
thread after workers finish; the full suite passed with that correction.
The installed `moldockpipe-clean` environment was verified to import this
checkout directly. Restart the application to load these source changes.
