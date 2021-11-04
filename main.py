import typer

app = typer.Typer(add_completion=False)


@app.command()
def hello(name: str):
    """Prints a hello message."""
    typer.echo(f"Hello {name}")


@app.command()
def spring_pendulum():
    """Runs a simple spring-pendulum physics simulation."""
    raise NotImplementedError()


if __name__ == "__main__":
    app()
