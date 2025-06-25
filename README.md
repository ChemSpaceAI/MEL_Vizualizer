# MEL_Vizualizer

MEL_Vizualizer is a Python Flask-based web application designed for visualizing and processing molecular data. The app generates images dynamically and manages the storage of generated images by cleaning up old files automatically.

## Features

- Flask web application with modular structure
- Dynamic image generation and storage
- Automatic cleanup of generated images when more than 30 PNGs accumulate in the `static/generated_images` folder
- Uses virtual environment for dependency management
- Supports deployment with WSGI servers (e.g. Gunicorn)
- Configuration and utility functions organized in `app` folder

## Project Structure

```
MEL_Visualizer/
│
├── .env                   # Environment variables (if any)
├── run_app.sh             # Shell script to run the app (optional)
├── test.ipynb             # Testing notebooks
│
├── app/                   # Main Flask application package
│   ├── __init__.py        # Flask app factory initialization
│   ├── config.py          # App configuration settings
│   │
│   ├── routes/            # Flask route/view functions
│   │   ├── __init__.py
│   │   └── main.py        # Main route handlers (index, views, etc.)
│   │
│   ├── functions/         # Core business logic and utilities
│   │   ├── __init__.py
│   │   ├── Chemistry_logic.py
│   │   ├── Enumeration_MEL_main.py
│   │   ├── Enumeration_MEL.py
│   │   ├── Enumerators.py
│   │   ├── Generation_of_MEL.py
│   │   ├── MinimalCaps_and_Synthons.py
│   │   ├── Molecules.py
│   │   └── utills.py
│   │
│   ├── static/            # Static files (images, css, js)
│   │   ├── Chemspace_logo.png
│   │   ├── generated_images/  # Generated image outputs
│   │   └── styles.css
│   │
│   └── templates/         # Jinja2 HTML templates
│       ├── base.html
│       ├── enumeration_view.html
│       ├── index.html
│       ├── macros.html
│       └── mel_type.html
│
├── wsgi.py                # WSGI entry point for production server
├── requirements.txt       # Python dependencies
└── structure.txt          # Optional: textual project structure overview

```

## Setup & Installation

1. Clone the repository:

```bash
git clone git@github.com:ChemSpaceAI/MEL_Vizualizer.git
cd MEL_Vizualizer
```

2. Create and activate a Python virtual environment:

```bash
python3 -m venv MelViz_Virtual_Env
source MelViz_Virtual_Env/bin/activate
```

3. Upgrade pip and install dependencies:

```bash
python -m ensurepip --upgrade
pip install --upgrade pip setuptools wheel
pip install -r requirements.txt
```

4. Run the application:

```bash
python wsgi.py
```

Or use the run script:

```bash
bash run_app.sh
```

## Usage

- The app automatically manages the generated images folder.
- If the number of PNG images in `/app/static/generated_images` exceeds 30, older images will be deleted to free space.

## Contributing

Feel free to open issues or submit pull requests to improve the project.

## License

Specify your project license here.

---

If you want me to include anything else like API docs, usage examples, or detailed config instructions, just say!
