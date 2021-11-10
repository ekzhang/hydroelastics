"""Hydroelastic contact force simulations in Taichi."""

import typer

app = typer.Typer(add_completion=False, help=__doc__)


@app.command()
def julia_set():
    """Runs the official demo animated render of a Julia Set."""
    from .julia_set import run

    run()


@app.command()
def test():
    """Runs the tests"""
    import subprocess

    subprocess.run("python -m pytest", shell=True)


@app.command()
def cloth():
    """Runs a simple cloth simulation based on linked springs."""
    from .cloth import run

    run()


@app.command()
def spring_pendulum():
    """Runs a simple spring-pendulum physics simulation."""
    from .spring_pendulum import run

    run()


if __name__ == "__main__":
    app()
