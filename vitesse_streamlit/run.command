#!/bin/zsh
set -eu
APP_DIR="$(cd -- "$(dirname -- "$0")" && pwd)"
cd "$APP_DIR"
APP_BUNDLED_PY="$HOME/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3"
APP_BUNDLED_DEPS="$APP_DIR/../../work/app_dependencies"
if [[ -x "$APP_BUNDLED_PY" && -d "$APP_BUNDLED_DEPS/streamlit" ]]; then
  export PYTHONPATH="$APP_BUNDLED_DEPS${PYTHONPATH:+:$PYTHONPATH}"
  exec "$APP_BUNDLED_PY" -m streamlit run app.py --server.address 127.0.0.1 --server.port 8501
fi
if [[ ! -x "$APP_DIR/.venv/bin/python" ]]; then
  APP_PYTHON="${VITESSE_PYTHON:-python3}"
  "$APP_PYTHON" -c 'import sys; sys.exit("Python 3.12+ required. Set VITESSE_PYTHON to a newer interpreter.") if sys.version_info < (3,12) else None'
  "$APP_PYTHON" -m venv "$APP_DIR/.venv"
  "$APP_DIR/.venv/bin/python" -m pip install -r requirements.txt
fi
exec "$APP_DIR/.venv/bin/python" -m streamlit run app.py --server.address 127.0.0.1 --server.port 8501
