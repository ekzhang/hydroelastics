# Hydroelastics

This repository contains code for a contact force simulator for soft objects, using a hydroelastic intersection model.

The simulator is differentiable and implemented using [Taichi](https://taichi.graphics/).

## Setup

The following requirements must be present in your environment.

- [Python 3.9](https://python.org/) and [Pip](https://pypi.org/project/pip/)
- [Just](https://github.com/casey/just): a command runner
- GPU drivers for CUDA, OpenGL, Metal, or Vulkan

Then, run the following command to install Python dependencies, including the Taichi runtime.

```bash
just install
```

To run the code in `main.py`, use the command:

```bash
just run
```

We will provide a command-line interface with options in the future.

## Development

We use [Black](https://black.readthedocs.io/en/stable/) to automatically format Python code. To format all of the files in the project, you can run the command:

```bash
just format
```

Note that there is a CI job that automatically checks formatting (and hopefully in the future, test cases) before pull requests can be merged.
