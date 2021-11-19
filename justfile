julia := "julia --project=."

install:
    {{julia}} -e "import Pkg; Pkg.instantiate()"

run *ARGS:
    {{julia}} {{ARGS}}

format:
    {{julia}} -e "using JuliaFormatter; format(\".\", verbose=true)"

check:
    {{julia}} -e "using JuliaFormatter; exit(!format(\".\", verbose=true))"

test:
    {{julia}} -e "import Pkg; Pkg.test()"
