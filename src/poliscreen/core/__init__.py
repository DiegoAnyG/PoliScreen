"""PoliScreen core: chemistry and engines, with no user interface.

Each engine is exposed behind a stable interface so the implementation can be
changed (or isolated in its own container) without touching its callers.
"""

from .design import AdmelabBridge, AdmelabError, DesignResult

__all__ = ["AdmelabBridge", "AdmelabError", "DesignResult"]
