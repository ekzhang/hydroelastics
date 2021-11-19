# Hydroelastics.jl

This repository contains code for a contact force simulator for soft objects, using a hydroelastic intersection model.

The simulator is differentiable and implemented using the [Julia](https://julialang.org/) programming language.

## Setup

The following requirements must be present in your environment.

- [Julia 1.6+](https://julialang.org/) and [Pkg](https://docs.julialang.org/en/v1/stdlib/Pkg/)
- [Just](https://github.com/casey/just): a command runner
- OpenGL-enabled graphics card used by [GLMakie](https://makie.juliaplots.org/dev/documentation/backends_and_output/)

Then, run the following command to instantiate the Julia package environment.

```bash
just install
```

To run an interactive Julia REPL or run a file, use the following command:

```bash
just run                       # repl
just run examples/pendulum.jl  # run an example
```

To unit tests specified in the `test/` folder, use the command:

```bash
just test
```

## Development

We use [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) to automatically format Python code. To format all of the files in the project, you can run the command:

```bash
just format
```

Note that there is a CI job that automatically checks formatting and test cases before pull requests can be merged.
