#!/bin/bash
## Convenience script to add, commit and push updates to gh-pages branch on GitHub 
## Author: Thomas Girke
## Last update: Apr 9, 2017

## (1) Makes sure you are in gh-pages branch
git checkout gh-pages  

## (2) Commit edits made in gh-pages branch 
git add -A :/
git commit -am "some edits"
git push -u origin gh-pages
echo "Committed/pushed changes to gh-pages branch on GitHub"


