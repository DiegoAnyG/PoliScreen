"""Nucleo de PoliScreen: química y motores, sin interfaz de usuario.

Cada motor se expone detras de una interfaz estable para poder cambiar la
implementación (o aislarla en su propio contenedor) sin tocar a quien la usa.
"""

from .design import AdmelabBridge, AdmelabError, DesignResult

__all__ = ["AdmelabBridge", "AdmelabError", "DesignResult"]
