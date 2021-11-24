julia := "julia --project=."

install:
    {{julia}} -e "import Pkg; Pkg.instantiate()"

repl:
    {{julia}}

notebook:
    {{julia}} -e "import Pluto; Pluto.run()"

format:
    {{julia}} -e "using JuliaFormatter; format(\".\", verbose=true)"

check:
    {{julia}} -e "using JuliaFormatter; exit(!format(\".\", verbose=true))"

test:
    {{julia}} -e "import Pkg; Pkg.test()"
