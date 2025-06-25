from flask import Flask
from dotenv import load_dotenv
import os

def create_app():
    load_dotenv()  # Load from .env file

    app = Flask(__name__)
    app.config.from_object("app.config.Config")

    # Register blueprints
    from app.routes.main import main_bp
    app.register_blueprint(main_bp)

    return app
