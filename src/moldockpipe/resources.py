"""Small OS resource probes without a mandatory monitoring dependency."""
import os


def process_tree_snapshot():
    """Best-effort Windows process tree sample; tolerate children exiting."""
    if os.name != "nt":
        return []
    import ctypes
    from ctypes import wintypes

    class Entry(ctypes.Structure):
        _fields_ = [("size", wintypes.DWORD), ("usage", wintypes.DWORD),
            ("pid", wintypes.DWORD), ("heap", ctypes.c_size_t), ("module", wintypes.DWORD),
            ("threads", wintypes.DWORD), ("parent", wintypes.DWORD),
            ("priority", wintypes.LONG), ("flags", wintypes.DWORD), ("name", wintypes.WCHAR * 260)]

    class Memory(ctypes.Structure):
        _fields_ = [("size", wintypes.DWORD), ("faults", wintypes.DWORD)] + [
            (name, ctypes.c_size_t) for name in ("peak_ws", "ws", "peak_paged", "paged",
                "peak_nonpaged", "nonpaged", "pagefile", "peak_pagefile", "private")]

    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.CreateToolhelp32Snapshot.argtypes = [wintypes.DWORD, wintypes.DWORD]
    kernel.CreateToolhelp32Snapshot.restype = wintypes.HANDLE
    kernel.Process32FirstW.argtypes = kernel.Process32NextW.argtypes = [wintypes.HANDLE, ctypes.POINTER(Entry)]
    kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel.OpenProcess.restype = wintypes.HANDLE
    kernel.CloseHandle.argtypes = [wintypes.HANDLE]
    kernel.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.POINTER(ctypes.c_uint64)] * 4
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Memory), wintypes.DWORD]
    snapshot = kernel.CreateToolhelp32Snapshot(2, 0)
    if snapshot == ctypes.c_void_p(-1).value:
        return []
    entries = {}
    try:
        entry = Entry()
        entry.size = ctypes.sizeof(entry)
        more = kernel.Process32FirstW(snapshot, ctypes.byref(entry))
        while more:
            entries[entry.pid] = (entry.parent, entry.name)
            more = kernel.Process32NextW(snapshot, ctypes.byref(entry))
    finally:
        kernel.CloseHandle(snapshot)
    family = {os.getpid()}
    while True:
        children = {pid for pid, (parent, _) in entries.items() if parent in family}
        if children <= family:
            break
        family |= children
    result = []
    for pid in sorted(family):
        parent, name = entries.get(pid, (None, ""))
        item = {"pid": pid, "parent_pid": parent, "name": name}
        handle = kernel.OpenProcess(0x1000 | 0x10, False, pid)
        if handle:
            try:
                memory = Memory()
                memory.size = ctypes.sizeof(memory)
                if psapi.GetProcessMemoryInfo(handle, ctypes.byref(memory), memory.size):
                    item.update(working_set_bytes=memory.ws, private_bytes=memory.private,
                                peak_working_set_bytes=memory.peak_ws)
                times = [ctypes.c_uint64() for _ in range(4)]
                if kernel.GetProcessTimes(handle, *(ctypes.byref(value) for value in times)):
                    item["cpu_seconds"] = (times[2].value + times[3].value) / 10**7
            finally:
                kernel.CloseHandle(handle)
        result.append(item)
    return result


def available_memory_bytes():
    if os.name == "nt":
        import ctypes

        class MemoryStatus(ctypes.Structure):
            _fields_ = [("length", ctypes.c_uint32), ("load", ctypes.c_uint32)] + [
                (name, ctypes.c_uint64) for name in (
                    "total_physical", "available_physical", "total_pagefile",
                    "available_pagefile", "total_virtual", "available_virtual", "extended")]

        status = MemoryStatus()
        status.length = ctypes.sizeof(status)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
            return min(status.available_physical, status.available_pagefile)
    else:
        try:
            # Linux MemAvailable includes reclaimable cache.
            from pathlib import Path
            for line in Path("/proc/meminfo").read_text().splitlines():
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
        except (OSError, ValueError):
            pass
    return None
