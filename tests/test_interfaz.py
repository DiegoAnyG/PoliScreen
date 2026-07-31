"""Pruebas de la interfaz con AppTest: recorre cada etapa y cada modo buscando excepciones.

Un fallo en una etapa solo se manifiesta al dibujarla, y la aplicación solo dibuja la activa. Sin
esto, un cambio en Resultados puede romper Ligandos y no notarse hasta que alguien lo abre.
"""
import os
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

APP = str(Path(__file__).resolve().parent.parent / "src" / "poliscreen" / "ui" / "streamlit_app.py")
ETAPAS = ["Receptores", "Ligandos", "Ejecutar", "Resultados"]
MODOS = ["Construir por reacción", "Generar péptidos", "Subir estructuras"]


def _app(tmp_path, **estado):
    at = AppTest.from_file(APP, default_timeout=600)
    at.run()
    at.session_state["proj_dir"] = str(tmp_path)
    for k, v in estado.items():
        at.session_state[k] = v
    at.run()
    return at


def _sin_excepcion(at, que):
    if at.exception:
        detalle = "; ".join(str(getattr(e, "value", e))[:200] for e in at.exception)
        pytest.fail(f"{que}: {detalle}")


@pytest.mark.parametrize("etapa", ETAPAS)
def test_cada_etapa_se_dibuja(tmp_path, etapa):
    _sin_excepcion(_app(tmp_path, etapa=etapa), f"etapa {etapa}")


@pytest.mark.parametrize("modo", MODOS)
def test_cada_modo_de_ligandos_se_dibuja(tmp_path, modo):
    _sin_excepcion(_app(tmp_path, etapa="Ligandos", modo_ligandos=modo), f"modo {modo}")


def test_la_carpeta_vacia_no_rompe_ninguna_etapa(tmp_path):
    """Un proyecto recien creado no tiene tablas ni poses: todas las vistas deben tolerarlo."""
    for etapa in ETAPAS:
        _sin_excepcion(_app(tmp_path, etapa=etapa), f"proyecto vacio en {etapa}")


def test_las_dos_vistas_del_visualizador(tmp_path):
    for vista in ("Resumen", "Complejo 3D"):
        _sin_excepcion(_app(tmp_path, etapa="Resultados", vis_res_vista=vista), f"vista {vista}")


def test_el_estado_sobrevive_al_cambiar_de_etapa(tmp_path):
    """Streamlit descarta los widgets no dibujados; sin reasignarlos se perdian los parámetros."""
    at = _app(tmp_path, etapa="Ligandos", modo_ligandos="Generar péptidos", pep_len=7)
    at.session_state["etapa"] = "Resultados"
    at.run()
    at.session_state["etapa"] = "Ligandos"
    at.run()
    assert at.session_state["modo_ligandos"] == "Generar péptidos"
    assert at.session_state["pep_len"] == 7
    _sin_excepcion(at, "ida y vuelta entre etapas")


def test_una_ruta_de_windows_no_crea_carpetas_con_barras(tmp_path):
    at = AppTest.from_file(APP, default_timeout=600)
    at.run()
    at.session_state["proj_dir"] = str(tmp_path).replace("/", "\\")
    at.run()
    _sin_excepcion(at, "ruta con barras invertidas")
    assert not [p for p in Path.cwd().iterdir() if "\\" in p.name]


def test_el_dialogo_de_descargas_se_abre(tmp_path):
    (tmp_path / "ranking.csv").write_text("compound,best_dock\na,-8.0\n")
    at = _app(tmp_path, etapa="Resultados", _abrir_descargas=True)
    _sin_excepcion(at, "dialogo de descargas")
    assert [c for c in at.checkbox if c.key and c.key.startswith("dl_")]
