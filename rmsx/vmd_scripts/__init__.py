"""
VMD helper scripts shipped with RMSX.

This subpackage contains:
- Python helpers for locating VMD on the host system
- Tcl scripts that VMD sources to build FlipBook-style layouts
"""

from .vmd_finder import find_vmd_executable

__all__ = ["find_vmd_executable"]

