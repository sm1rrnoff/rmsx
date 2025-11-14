import os
import platform
import shutil
import re
from pathlib import Path


# Find the best available VMD executable on this system.
# Checks an explicit path first, then $VMD_EXECUTABLE, then PATH,
# then common install locations (macOS, Linux, Windows).
# Returns the full path to the VMD executable if found, otherwise raises FileNotFoundError.

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def parse_version_from_str(s):
    """
    Extract version numbers from a string like:
      "VMD 2.0.0.app"            -> (2, 0, 0)
      "VMD 2.0.0a7-pre2.app"     -> (2, 0, 0)
      "vmd-1.9.4a57"             -> (1, 9, 4)
    If no version is found, return (0,).
    """
    m = re.search(r"(\d+(?:\.\d+)+)", s)
    if not m:
        return (0,)
    return tuple(int(x) for x in m.group(1).split("."))


def candidate(version, path, source_priority):
    """
    Create a sortable candidate tuple: (source_priority, version_tuple, Path)
    Higher source_priority wins. Then higher version wins.
    """
    return (source_priority, version, Path(path))


# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------

def find_vmd_executable(explicit_path=None):
    """
    Locate the best available VMD executable on this system.

    Priority order:
      1. Explicit path provided by user.
      2. Environment variable VMD_EXECUTABLE.
      3. Real versioned installs in typical locations (macOS .app bundles, /opt/vmd-x.y.z, Windows “Program Files University of Illinois VMD”, etc.).
      4. Executable found on PATH (“vmd” command).
      5. Generic fallback known install paths.

    If none is found, raises FileNotFoundError with guidance to install VMD or adjust PATH.
    """
    candidates = []

    # 1. Explicit user-given path (highest priority)
    if explicit_path:
        candidates.append(candidate((9999,), explicit_path, source_priority=100))

    # 2. Environment variable override
    env_path = os.environ.get("VMD_EXECUTABLE")
    if env_path:
        candidates.append(candidate((9998,), env_path, source_priority=90))

    # 3. OS-specific known locations
    system = platform.system()

    if system == "Darwin":
        # macOS: look for VMD*.app bundles in standard Applications roots
        for base in [Path("/Applications"), Path.home() / "Applications"]:
            if not base.exists():
                continue
            for app in base.glob("VMD*.app"):
                version = parse_version_from_str(app.name)

                # Known plausible internal layouts of VMD app bundles:
                candidate_paths = [
                    app / "Contents" / "Resources" / "VMD.app" / "Contents" / "MacOS" / "VMD",
                    app / "Contents" / "MacOS" / "VMD",
                    app / "Contents" / "vmd2" / "lib" / "vmd_MACOSXARM64",
                    app / "Contents" / "vmd2" / "vmd_MACOSXARM64",
                ]

                for exe in candidate_paths:
                    candidates.append(candidate(version, exe, source_priority=50))

    elif system == "Linux":
        # Versioned installs under /opt
        opt = Path("/opt")
        if opt.exists():
            for d in opt.glob("vmd*"):
                version = parse_version_from_str(d.name)
                exe = d / "bin" / "vmd"
                candidates.append(candidate(version, exe, source_priority=50))
        # Typical system paths
        for loc in ["/usr/local/bin/vmd", "/usr/bin/vmd"]:
            candidates.append(candidate((0,), loc, source_priority=10))
        # Additional library/data directory (for version detection even if binary elsewhere) per guide:
        # e.g. /usr/local/lib/vmd is often used as INSTALLLIBDIR on Unix systems. :contentReference[oaicite:0]{index=0}

    elif system == "Windows":
        # Program Files directories
        pf       = Path(os.environ.get("ProgramFiles",      r"C:\Program Files"))
        pf_x86   = Path(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)"))

        # Look under both roots for “VMD*” folders (including versioned ones) and the known default directory
        for root in [pf, pf_x86]:
            if not root.exists():
                continue
            # Known default library/data install path: C:\Program Files \ Uiniversity of Illinois \ VMD :contentReference[oaicite:1]{index=1}
            default_dir = root / "University of Illinois" / "VMD"
            candidates.append(candidate(parse_version_from_str(default_dir.name), default_dir / "vmd.exe", source_priority=60))

            # Search versioned VMD folders
            for d in root.glob("VMD*"):
                version = parse_version_from_str(d.name)
                exe = d / "vmd.exe"
                candidates.append(candidate(version, exe, source_priority=50))

        # Generic fallback
        candidates.append(candidate((0,), pf       / "VMD" / "vmd.exe", source_priority=10))
        candidates.append(candidate((0,), pf_x86   / "VMD" / "vmd.exe", source_priority=10))

    # 4. PATH‐based discovery (lowest real priority for binary in PATH)
    path_exe = shutil.which("vmd")
    if path_exe:
        candidates.append(candidate((0,), path_exe, source_priority=5))

    # 5. Sort candidates and pick the best valid one
    # Sort by (source_priority desc, version desc)
    candidates.sort(key=lambda c: (c[0], c[1]), reverse=True)

    for source_priority, version, path in candidates:
        try:
            if path.is_file() and os.access(path, os.X_OK):
                return str(path)
        except Exception:
            # In case Path.is_file or os.access throws; ignore and continue
            pass

    # No valid executable found
    raise FileNotFoundError(
        "Could not find a VMD executable in standard locations.\n"
        "If you already have VMD installed at a custom location, one of the following will help:\n"
        "  • Set the environment variable VMD_EXECUTABLE to the full binary path\n"
        "  • Add the directory containing the VMD binary to your system PATH\n"
        "And then re-run your script."
    )
