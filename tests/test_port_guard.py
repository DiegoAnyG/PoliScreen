"""Starting the interface on a port someone else already holds.

Streamlit just exits, and in the Windows installer the console closes with it, so the reason is
gone before it can be read. What answered the browser then was another PoliScreen entirely — a
container Docker Desktop had restarted on its own — serving a different filesystem, which is why
the project paths pointed somewhere impossible.
"""
import socket

import pytest

from poliscreen import cli


@pytest.fixture
def occupied_port():
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("127.0.0.1", 0))
    s.listen(1)
    yield s.getsockname()[1]
    s.close()


@pytest.fixture
def free_port():
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


def test_a_busy_port_is_seen(occupied_port):
    assert cli.port_in_use(occupied_port)


def test_a_free_port_is_not(free_port):
    assert not cli.port_in_use(free_port)


def test_the_interface_refuses_rather_than_starting_a_second_one(occupied_port, capsys, monkeypatch):
    import subprocess
    started = []
    monkeypatch.setattr(subprocess, "call", lambda *a, **k: started.append(a) or 0)

    args = type("A", (), {"port": occupied_port, "expose": False})()
    assert cli.cmd_ui(args) == 1, "it returned success while refusing to serve"
    assert not started, "it launched Streamlit anyway, on a port already taken"

    said = capsys.readouterr().err
    assert "docker ps" in said, "the usual cause is not even mentioned"
    assert str(occupied_port + 1) in said, "no way out is offered"
    assert "NOT the copy you just launched" in said, "the dangerous part is left unsaid"


def test_the_browser_waits_for_the_page_instead_of_a_guessed_delay(free_port, monkeypatch):
    """Double-clicking the launcher has to end in a browser.

    The address printed on the console is clickable only with Ctrl, which nobody outside a terminal
    knows about; and a first start that has to unpack Streamlit takes longer than any delay worth
    hard-coding, so the page is opened when the port answers, not when a timer says so.
    """
    import webbrowser
    opened = []
    monkeypatch.setattr(webbrowser, "open", opened.append)
    answers = iter([False, False])
    monkeypatch.setattr(cli, "port_in_use", lambda p: next(answers, True))

    url = f"http://localhost:{free_port}"
    cli._open_when_serving(url, free_port).join(timeout=10)

    assert opened == [url], "the interface was left as a link nobody knows how to click"
