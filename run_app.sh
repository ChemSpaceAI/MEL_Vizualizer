#!/bin/bash

export FLASK_APP=wsgi.py
export FLASK_ENV=development
nohup gunicorn -w 1 --threads 5 -b 10.244.32.19:3937 wsgi:app > gunicorn.log 2>&1 &