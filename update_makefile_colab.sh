#!/bin/bash
git add .
git rm update_makefile_colab
git commit -m "Update Makefile"
git pull
git push origin main
