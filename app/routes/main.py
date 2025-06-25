from flask import Blueprint, render_template, request, redirect, url_for

from app.functions.utills import get_reaction_ids, cleanup_generated_images
from app.functions.MinimalCaps_and_Synthons import get_random_synthons_for_reaction_id, get_minimal_caps
from app.functions.Generation_of_MEL import get_generation
from app.functions.Enumeration_MEL_main import get_enumeration
from app.config import Config
main_bp = Blueprint("main", __name__)

@main_bp.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        scope = request.form.get("reaction_scope", "")
        mel_type = request.form.get("mel_type", "")
        if mel_type and scope:
            return redirect(url_for("main.view_by_mel_type", mel_type=mel_type, reaction_scope=scope))
    return render_template("index.html")

@main_bp.route("/reaction_id", methods=["GET", "POST"])
def view_by_mel_type():

    cleanup_generated_images(Config.GENERATED_IMAGES_FOLDER, max_files=30)  # Clean old images
    mel_type = request.values.get("mel_type")
    reaction_id = request.values.get("reaction_id")
    action = request.values.get("action", "")
    reaction_scope = request.values.get("reaction_scope", "")

    if not mel_type:
        # Optionally handle missing mel_type, e.g. redirect or error
        return redirect(url_for("main.index"))

    reaction_ids = get_reaction_ids(mel_type) if mel_type else []

    result = {}
    synthons_info = {}
    caps_info = {}
    generation_info = {}
    enum_caps_info = {}
    enumeration_info = {}
    enum_tables_info = {}
    if reaction_id:

        if action == "show_REAL_reaction":
            synthons_info, _ = get_random_synthons_for_reaction_id(reaction_id)
        elif action == "show_generation_mincaps":
            caps_info, _ = get_minimal_caps(reaction_id, mode="generation")
        elif action == "show_enumeration_mincaps":
            enum_caps_info, _ = get_minimal_caps(reaction_id, mode="enumeration")
        elif action == "show_generation_example":
            generation_info, _ = get_generation(reaction_id)
        elif action == "show_enumeration_example":
            enumeration_info,  _,  enum_tables_info = get_enumeration(reaction_id)

    return render_template(
        "mel_type.html",
        result=result,
        reaction_id=reaction_id,
        mel_type=mel_type,
        reaction_scope=reaction_scope,
        reaction_ids=reaction_ids,
        
        action=action,

        synthons_info=synthons_info,
        caps_info=caps_info,
        enum_caps_info=enum_caps_info,
        generation_info=generation_info,
        enumeration_info=enumeration_info,
        enum_tables_info=enum_tables_info,

    )

    