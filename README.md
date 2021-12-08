# Hydroelastics.jl

This repository contains code for a contact force simulator for soft objects, using a hydroelastic pressure field intersection model.

<p align="center">
  <img src="https://i.imgur.com/g54zdzJ.gif" width="720">
</p>

The simulator is implemented using the [Julia](https://julialang.org/) programming language.

## Setup

The following requirements must be present in your environment.

- [Julia 1.7+](https://julialang.org/)
- [Just](https://github.com/casey/just): a command runner

Then, run the following command to instantiate the Julia package environment.

```bash
just install
```

To run the interactive [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) or open a [Pluto](https://github.com/fonsp/Pluto.jl) notebook server, use the following commands:

```bash
just repl      # open Julia REPL (can add packages here)
just notebook  # start server (try opening examples/cube.jl)
```

To unit tests specified in the `test/` folder, use the command:

```bash
just test
```

## Development

We use [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) to automatically format Julia code. To format all of the files in the project, you can run the command:

```bash
just format
```

Note that there is a CI job that automatically checks formatting and test cases before pull requests can be merged.

## Acknowledgements

Made by [Vincent Huang](https://twitter.com/vvhuang_), [Franklyn Wang](https://twitter.com/franklyn_wang), and [Eric Zhang](https://twitter.com/ekzhang1). Licensed under the [MIT license](LICENSE).
