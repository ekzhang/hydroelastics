python := "python3"

install:
    {{python}} -m pip install -r requirements.txt

run *ARGS:
    {{python}} -m hydroelastics.main {{ARGS}}

format:
    {{python}} -m black .

check:
    {{python}} -m black --check .

test:
    {{python}} -m pytest
