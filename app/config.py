import os

class Config:
    # This ensures the path is absolute and environment-configurable
    BASE_DIR = os.path.abspath(os.path.dirname(__file__))

    # Static image directory path
    IMG_DIR = os.getenv(
        "IMG_DIR", 
        os.path.join(BASE_DIR, "static", "generated_images")
    )

    SECRET_KEY = os.getenv("SECRET_KEY", "default-secret")
    DEBUG = os.getenv("DEBUG", "false").lower() == "true"

    BASE_DATA_PATH = os.getenv("MEL_DATA_PATH", "/chemai_narnia/MEL_Package/MEL_project_app/data/")    
    IDS_FILE_PATH = os.getenv("IDS_FILE", "/chemai_narnia/MEL_Package/MEL_project_app/data/REACTION_file_for_MEL.tsv")
    SYNTHON_DATA_PATH = os.getenv("SYNTHON_DATA_PATH", "/chemai_narnia/MEL_Package/MEL_project_app/data/SYNTHONS_by_reaction_id_modified_for_MEL")
    MIN_CAP_DATA_PATH = os.getenv("MIN_CAP_DATA_PATH", "/chemai_narnia/MEL_Package/MEL_project_app/data/Minimal_Caps_by_reaction_id_for_generation_MEL")
    MEL_PRODUCT_DATA_PATH = os.getenv("MEL_PRODUCT_DATA_PATH", "/chemai_narnia/MEL_Package/MEL_project_app/data/Generated_MEL/enumerated_MEL")
    GENERATED_IMAGES_FOLDER = os.getenv("GENERATED_IMAGES_FOLDER", "/chemai_narnia/MEL_Vizualizer/app/static/generated_images" )