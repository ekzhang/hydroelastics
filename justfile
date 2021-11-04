python := "python3"

install:
    {{python}} -m pip install -r requirements.txt

run *ARGS:
    {{python}} main.py {{ARGS}}

format:
    {{python}} -m black .

check:
    {{python}} -m black --check .
